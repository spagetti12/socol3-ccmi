 !*************************************************************************
 MODULE messy_o3orig
 !*************************************************************************

   ! MODULE FOR CALCULATIN DIAGNOSTIC OZONE ORIGIN TRACER
   !
   ! MESSy-SMCL
   !
   ! Description: Grewe, V., The origin of ozone, ACP 6, 1495-1511, 2006. 
   !              http://www.atmos-chem-phys.net/6/1495/2006/acp-6-1495-2006.pdf
   !  Briefly:   d/dt O3(i)=PO3(i) - DO3/O3 * O3(i)
   !             O3(i) Diagnostic ozone tracer for region i  [molec/cm3]
   !             PO3(i) Ozone production in region i         [molec/cm3/s]
   !             DO3   Ozone loss                            [molec/cm3/s]
   !             
   !   11.04.2006: Status
   !       Input from Mecca    O3prod Brutto Ozone Production [molec/cm3/s]
   !                           o3m1 pre mecca ozone           [molec/cm3]
   !                           o3   post mecca ozone          [molec/cm3]
   !                           cO3 Chemical ozone change      [molec/cm3/timestep]
   !       Input from file     o3_orig_regions                [number]
   !       Input from Tropop   pressure of tropopause         [Pa]
   !       Input from ECHAM    press  pressure of levels      [Pa]    
   !       Output              o3orig1-5 bzw 9                [molec/cm3]
   !       Output              pos/neg Error due to transport [molec/cm3]
   !       Output              o3loss                         [molec/cm3/s]
   !       Output              status
   !
   !       Integration method: Euler implicit backward      
   !               Ozone loss (here pos. definit) is calculated from chemical ozone change and 
   !                          ozone production so that ozone change calculated
   !                          by MECCA gives the sameresults as an integration with 
   !                          the integration scheme used for O3ORIG.
   !               O3Prod is adapted whenever o3loss is negative.
   !       Initialisation is performed when sum of o3orig = 0.: 
   !                          Equally distributed. (Could be done better of course)
   !       Scaling: a) Prior to each integration step sum of all o3origin tracers is 
   !                          scaled to chemical ozone (o3m1), to prevet divergence between diagnostic
   !                          o3orig fields and 'chemical' ozone. 
   !                          The calculated deficit gives a source for negative error estimate
   !                          The calculated excess gives a source for positive error estimate
   !                b) After each integration step: sum of all o3origin tracers is scaled to 
   !                          chemical ozone (o3). This shouldn't be necessary since ozone loss 
   !                          is derived implicitely. However, it doesn't harm.
   ! Authors:
   !    Volker Grewe, DLR, March 2006 and October 2010
   !    Patrick Joeckel, MPICH, March 2007
   !    Stefanie Meul and Sophie Oberlaender, FUB, October 2010
   !    Hella Garny, DLR, November 2012 (adjust regions; no use of TP )

   ! ----------- >

   USE messy_main_constants_mem, ONLY: DP, SP, I4, I8

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: DP, SP, I4, I8

   ! ----------- <

   ! GLOBAL PARAMETERS
   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'o3orig'
   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.2'

 ! REAL(DP),PUBLIC :: M=0.5_dp
   INTEGER,PUBLIC,PARAMETER     ::  n_trac_orig=16
   INTEGER,PUBLIC  :: i_integrate=1, i_trac_orig=16, i_trac_err   !sm_01092010
   LOGICAL,PUBLIC  :: l_err=.true.                     ! Switch on error diagnostics ?
   CHARACTER(LEN=20), DIMENSION(n_trac_orig),PUBLIC   :: sn_o3orig         ! Short Names of the variables
   ! SUBROUTINES
   PUBLIC :: o3orig_integrate                      ! integration of ozone origin tracers

   PUBLIC :: o3orig_read_nml

 CONTAINS

 ! ----------------------------------------------------------------------
   SUBROUTINE o3orig_integrate(kproma,nlev,           & ! Input  (Grid)
   dtime,                                             & ! Input  (Timestep) 
   o3_orig_regions,tp,                                & ! Input  (Regions)
                 press,                                             & ! Input  (Regions)
   o3loss,o3prod,o3,o3m1,                             & ! Input  (Chem)
   po3orig,                                            & ! In- Output (O3orig)
   status                                             & ! Output
   )
 ! 2nd try: Take only Prod and O3, O3m1 as known values.
 !        Calculate from that o3loss
     IMPLICIT NONE

     ! I/O
     INTEGER,                   INTENT(IN) :: kproma,nlev           ! Input
     REAL(DP),                  INTENT(IN) :: dtime                 ! Input [s]
     REAL(DP), DIMENSION(:),    INTENT(IN) :: o3_orig_regions,tp    ! Input [s]
     REAL(DP), DIMENSION(:,:),  INTENT(IN) :: press                 ! Input [Pa]
     REAL(DP), DIMENSION(:,:),  INTENT(INOUT) :: o3prod             ! Input [molec/cm3/s]
     REAL(DP), DIMENSION(:,:),  INTENT(OUT):: o3loss                ! Output [molec/cm3/s]
     REAL(DP), DIMENSION(:,:),  INTENT(IN) :: o3,o3m1               ! Input [molec/cm3]
     REAL(DP), DIMENSION(:,:,:),INTENT(INOUT):: po3orig             ! In/Output [molec/cm3]
     INTEGER,                   INTENT(OUT):: status

     ! LOCAL
     REAL (DP), DIMENSION(kproma,nlev)     :: rescaleO3     ! rescaling parameter 
     REAL (DP), DIMENSION(kproma,nlev)     :: totO3_orig    ! total ozone for calculation with the tropopause as boundary 
     
     REAL (DP), DIMENSION(kproma,nlev)     :: prod_err, loss_err  ! prod and loss for error analysis [molec/cm3/timestep]
     REAL (DP), DIMENSION(0:14)            :: prod          ! HARD WIRED !!!!  !sm_01092010
     REAL (DP)                             :: nenner
     INTEGER                               :: i,j,k,icount1,icount2, i_trac

 ! INIT
     status=0 ! No Error
     icount1=0 ! No initialization
     icount2=0 ! No error in initialization
     
 ! Preset variables
 ! Sum of ozone from different origin after transport, should be equal to o3m1 
 !                   in case of perfect transport scheme
     totO3_orig(:,:)=0.
     do i_trac=1,i_trac_orig-i_trac_err
      totO3_orig(:,:)=totO3_orig(:,:)+po3orig(:,:,i_trac)
     enddo
  
  
 ! Difference in ozone between 'chemical' ozone and sum of ozone origin tracers gives 
 !   production term for positive/negative error estimate.
 !   First: preset values, later redefine

     prod_err(:,:)=totO3_orig(:,:)-o3m1(:,:)
     loss_err(:,:)=0._dp

        
     do k=1,nlev
     do i=1,kproma
!             redefine production term for errors according to sign
       if (prod_err(i,k).lt.0._dp) then
          loss_err(i,k)=-prod_err(i,k)
          prod_err(i,k)=0._dp
       endif

!             Initialise ozone tracers at first time step, i.e. tot_o3==0.
       if (totO3_orig(i,k).eq.0._dp .and.o3m1(i,k).gt.0._dp) then
          icount1=icount1+1

          if (int(o3_orig_regions(i)).eq.1) then
             if (press(i,k).ge.20000._dp) then !NHTS
                     po3orig(i,k,1)=o3m1(i,k)                
             elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NPLS   
                     po3orig(i,k,4)=o3m1(i,k)
             elseif (press(i,k).lt.3000._dp) then  !NPUS
                     po3orig(i,k,10)=o3m1(i,k)
             endif
          elseif (int(o3_orig_regions(i)).eq.2) then
             if (press(i,k).ge.20000._dp) then !NHTS
                     po3orig(i,k,1)=o3m1(i,k)              
             elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                     po3orig(i,k,5)=o3m1(i,k)
             elseif (press(i,k).lt.3000._dp) then  !NMUS
                     po3orig(i,k,11)=o3m1(i,k)
             endif
          elseif (int(o3_orig_regions(i)).eq.3) then
             if (press(i,k).ge.10000._dp) then !TRTS
                     po3orig(i,k,2)=o3m1(i,k)                
             elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                     po3orig(i,k,5)=o3m1(i,k)
             elseif (press(i,k).lt.3000._dp) then  !NMUS
                     po3orig(i,k,11)=o3m1(i,k)
             endif
          elseif (int(o3_orig_regions(i)).eq.4) then
             if (press(i,k).ge.10000._dp) then !TRTS
                     po3orig(i,k,2)=o3m1(i,k)                
             elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !TRLS   
                     po3orig(i,k,6)=o3m1(i,k)
             elseif (press(i,k).lt.3000._dp.and.press(i,k).ge.500._dp) then  !TRMS
                     po3orig(i,k,7)=o3m1(i,k)
             elseif (press(i,k).lt.500._dp) then   !TRUS
                     po3orig(i,k,12)=o3m1(i,k)
             endif
          elseif (int(o3_orig_regions(i)).eq.5) then
             if (press(i,k).ge.10000._dp) then !TRTS
                     po3orig(i,k,2)=o3m1(i,k)                
             elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                     po3orig(i,k,8)=o3m1(i,k)
             elseif (press(i,k).lt.3000._dp) then  !SMUS
                     po3orig(i,k,13)=o3m1(i,k)
             endif
          elseif (int(o3_orig_regions(i)).eq.6) then
             if (press(i,k).ge.20000._dp) then !SHTS
                     po3orig(i,k,3)=o3m1(i,k)                
             elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                     po3orig(i,k,8)=o3m1(i,k)
             elseif (press(i,k).lt.3000._dp) then  !SMUS
                     po3orig(i,k,13)=o3m1(i,k)
             endif
          elseif (int(o3_orig_regions(i)).eq.7) then
             if (press(i,k).ge.20000._dp) then !SHTS
                     po3orig(i,k,3)=o3m1(i,k)                
             elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SPLS   
                     po3orig(i,k,9)=o3m1(i,k)
             elseif (press(i,k).lt.3000._dp) then  !SPUS
                     po3orig(i,k,14)=o3m1(i,k)
             endif

          endif 
!sm_so_17092010-
       elseif (totO3_orig(i,k).eq.0._dp .and.o3m1(i,k).eq.0._dp) then
          icount2=icount2+1
       else
                 ! rescale ozone origin tracers to avoid differences to 'chemical' ozone
          IF (totO3_orig(i,k) > 0._dp) THEN
             rescaleO3(i,k)=o3m1(i,k)/totO3_orig(i,k)
          ELSE
             rescaleO3(i,k)= 0._dp
          ENDIF
                po3orig(i,k,1:i_trac_orig-i_trac_err)=po3orig(i,k,1:i_trac_orig-i_trac_err)*rescaleO3(i,k)
       endif

     enddo
     enddo

! Output in case of initialistion
!    if (icount1.gt.0) write(*,*) modstr,' Grid points initialized: ',icount1
!    if (icount2.gt.0) write(*,*) modstr,  & 
!              ' Grid points uninitializable: ',icount2
!    if (icount1.gt.0) write(*,*) 'Regions',(o3_orig_regions(i),i=1,kproma)

! Use only ozone production calculated by mecca and derive ozone loss consistent 
!      with time integration scheme, so that chemical ozone changes (o3-o3m1) 
!      would be exactly reproduced 
! Determine ozone loss rate     positive definit ?!
! Euler implicit backward
!um_ak_20110719+
!     o3loss(:,:)=( o3prod(:,:) - (o3(:,:)-o3m1(:,:)) /dtime ) * o3m1(:,:) / o3(:,:)
     WHERE (o3(:,:) > 0._dp)
        o3loss(:,:)=( o3prod(:,:) - &
             (o3(:,:)-o3m1(:,:))/dtime ) * o3m1(:,:) &
             / o3(:,:)
     ELSEWHERE
        o3loss(:,:)= 0._dp
     END WHERE
!um_ak_20110719-
! euler  not yet implemented
!    o3loss(:,:)=( o3prod(:,:) - (o3(:,:)-o3m1(:,:)) /dtime )   

     do k=1,nlev
     do i=1,kproma
        
! Adjust Production in case of underestimation. Min values better than 0?
        if (o3loss(i,k).le.0._dp) then
                  o3prod(i,k)=o3prod(i,k)-o3loss(i,k)
                  o3loss(i,k)=0.
        endif

! Determine Production region
!               Prod(:)=0.
!               if (press(i,k).ge.tp(i)) then ! Troposphere
!                 Prod(:)=0.
!          Prod(int(o3_orig_regions(i)))=o3prod(i,k)
!               else
!                 Prod(5)=o3prod(i,k)
!               endif
!               Prod(4)=Prod(0)          ! Sollte besser gemacht werden Rest =4 nicht 0

   !sm_01092010+
                        
      Prod(:)=0._dp
      if (int(o3_orig_regions(i)).eq.1) then
         if (press(i,k).ge.20000._dp) then !NHTS
                     Prod(1)=o3prod(i,k)                
         elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NPLS   
                    Prod(4)=o3prod(i,k)
         elseif (press(i,k).lt.3000._dp) then  !NPUS
                    Prod(10)=o3prod(i,k)
         endif
      elseif (int(o3_orig_regions(i)).eq.2) then
         if (press(i,k).ge.20000._dp) then !NHTS
                    Prod(1)=o3prod(i,k)                
         elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                    Prod(5)=o3prod(i,k)
         elseif (press(i,k).lt.3000._dp) then  !NMUS
                    Prod(11)=o3prod(i,k)
         endif
      elseif (int(o3_orig_regions(i)).eq.3) then
         if (press(i,k).ge.10000._dp) then !TRTS
                     Prod(2)=o3prod(i,k)                
         elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                     Prod(5)=o3prod(i,k)
         elseif (press(i,k).lt.3000._dp) then  !NMUS
                     Prod(11)=o3prod(i,k)
         endif
      elseif (int(o3_orig_regions(i)).eq.4) then
         if (press(i,k).ge.10000._dp) then !TRTS
                     Prod(2)=o3prod(i,k)                
         elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !TRLS   
                     Prod(6)=o3prod(i,k)
         elseif (press(i,k).lt.3000._dp.and.press(i,k).ge.500._dp) then  !TRMS
                     Prod(7)=o3prod(i,k)
         elseif (press(i,k).lt.500._dp) then   !TRUS
                     Prod(12)=o3prod(i,k)
         endif
      elseif (int(o3_orig_regions(i)).eq.5) then
         if (press(i,k).ge.10000._dp) then !TRTS
                     Prod(2)=o3prod(i,k)                
         elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                     Prod(8)=o3prod(i,k)
         elseif (press(i,k).lt.3000._dp) then  !SMUS
                     Prod(13)=o3prod(i,k)
         endif
      elseif (int(o3_orig_regions(i)).eq.6) then
         if (press(i,k).ge.20000._dp) then !SHTS
                     Prod(3)=o3prod(i,k)                
         elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                     Prod(8)=o3prod(i,k)
         elseif (press(i,k).lt.3000._dp) then  !SMUS
                    Prod(13)=o3prod(i,k)
         endif
      elseif (int(o3_orig_regions(i)).eq.7) then
         if (press(i,k).ge.20000._dp) then !SHTS
                     Prod(3)=o3prod(i,k)                
         elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SPLS   
                     Prod(9)=o3prod(i,k)
         elseif (press(i,k).lt.3000._dp) then  !SPUS
                     Prod(14)=o3prod(i,k)
         endif

      endif 

       
    !sm_01092010-  
                  
                          

! Integrate: only Euler backward implemented so far

! Ozone origin
      IF (o3m1(i,k)> 0._dp) THEN
         nenner = 1._dp+o3loss(i,k)/o3m1(i,k)*dtime
      ELSE
         nenner = 1._dp
      ENDIF
      po3orig(i,k,1:i_trac_orig-i_trac_err)=     &
           (po3orig(i,k,1:i_trac_orig-i_trac_err)+Prod(1:i_trac_orig-i_trac_err)*dtime)/nenner    

      if (l_err) then
! Error analysis
         IF (o3m1(i,k)> 0._dp) THEN
            nenner=1._dp+(o3loss(i,k)+loss_err(i,k)/o3m1(i,k)/dtime)/o3m1(i,k)*dtime
         ELSE
            nenner = 1._dp
         ENDIF
        po3orig(i,k,i_trac_orig-1)=(po3orig(i,k,i_trac_orig-1)+prod_err(i,k))/nenner ! positive error 
        po3orig(i,k,i_trac_orig)  =(po3orig(i,k,i_trac_orig)  +loss_err(i,k))/nenner ! negative error 
      endif

! Euler forward explicit
!      o3orig_1(i,k)=o3orig_1(i,k)+Prod(1)*dtime-o3orig_1(i,k)/o3m1(i,k)*o3loss(i,k)*dtime

    enddo ! kproma
    enddo ! levs

! make sure that ozone and tagged ozone experience same changes, shouldn't be necessary!
    totO3_orig(:,:)=0._dp
    do i_trac=1,i_trac_orig-i_trac_err
      totO3_orig(:,:)=totO3_orig(:,:)+po3orig(:,:,i_trac)
    enddo
   
    WHERE (totO3_orig(:,:) > 0.0_dp)
       rescaleO3(:,:)=o3(:,:)/totO3_orig(:,:)
    ELSEWHERE
       rescaleO3(:,:)=0.0_dp
    END WHERE

    do i_trac=1,i_trac_orig-i_trac_err
      po3orig(:,:,i_trac)=po3orig(:,:,i_trac)*rescaleO3(:,:)
    enddo
    
    if (icount1.gt.0 .and. l_err) then
      po3orig(:,:,i_trac_orig-1)=0._dp  !sm_02092010
      po3orig(:,:,i_trac_orig)=0._dp    !sm_02092010
    endif
     
    status = icount2

  END SUBROUTINE o3orig_integrate
! ----------------------------------------------------------------------

! ************************************************************************
! O3ORIG
! ************************************************************************
! ----------------------------------------------------------------------
  SUBROUTINE o3orig_read_nml(status, iou)

    ! Radon MODULE ROUTINE (CORE)
    !
    ! READ o3orig NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/  i_integrate,  &       ! 1: Euler backward   2: euler forward
                                           !  besser �ber CPL tp_4_o3_orig,   &     ! 1: WMO; 2: PV
                     i_trac_orig, &        ! Number of tracers 
                     sn_o3orig,   &        ! Short name of the tracers
                     l_err                 ! Error diagnostic on?
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'o3orig_read_nml'
    LOGICAL                 :: lex          ! file exists ?
    INTEGER                 :: fstat        ! file status
    INTEGER                 :: i_trac

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    ! CHECK TRACER INITIALIZATION
      SELECT CASE(i_integrate)
       CASE(1)     ! Euler backward
         WRITE(*,*) '... Euler backward integration scheme for ozone origin tracers'
       CASE(2)     ! Euler forward
         WRITE(*,*) '... Euler forward integration scheme for ozone origin tracers'
         WRITE(*,*) '... Not yet implemented'
         RETURN
       CASE DEFAULT
         WRITE(*,*) '... Incorrect value for i_integrate  ',i_integrate
         RETURN
       END SELECT
  
       IF (L_ERR) then
         WRITE(*,*) 'Error diagnostics switched on (=L_ERR)' 
         i_trac_err=2
       ELSE
         WRITE(*,*) 'Error diagnostics switched off (=L_ERR)'
         i_trac_err=0
       ENDIF
       write(*,*) '...',i_trac_orig-i_trac_err,'regions and ',i_trac_err,' error-tracers in total:',i_trac_orig,'(=i_trac_orig)'
    !
    ! CHECK CONSTAN/ONLINE/OFFLINE FLUX



    do i_trac=1,i_trac_orig
       write(*,*) 'O3ORIG-Tracer',i_trac,': ',sn_o3orig(i_trac)
    enddo

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    END SUBROUTINE o3orig_read_nml

!*************************************************************************
END MODULE messy_o3orig
!*************************************************************************

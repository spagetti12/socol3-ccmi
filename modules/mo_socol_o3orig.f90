         !*************************************************************************
         MODULE mo_socol_o3orig
         !*************************************************************************

           ! MODULE FOR CALCULATING DIAGNOSTIC OZONE ORIGIN TRACER
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

           USE mo_kind, ONLY: DP
           USE mo_memory_base,   ONLY: new_stream, add_stream_element,  &
                                      default_stream_setting, add_stream_reference, &
                                      delete_stream, t_stream
           USE mo_exception,     ONLY: message, message_text, finish

           IMPLICIT NONE

           PRIVATE

           INTEGER,PUBLIC,PARAMETER ::  n_trac_orig=21
           INTEGER,PUBLIC  :: i_integrate=1, i_trac_orig=21, i_trac_err   !sm_01092010
           LOGICAL,PUBLIC  :: l_err=.true.                     ! Switch on error diagnostics ?
           CHARACTER(LEN=20), DIMENSION(n_trac_orig),PUBLIC   :: sn_o3orig         ! Short Names of the variables

          ! streams

          REAL(DP), DIMENSION(:,:,:), PUBLIC, POINTER :: o3=>NULL()                ! kproma*nlev
          !!ASREAL(DP), DIMENSION(:,:),   POINTER :: tp=>NULL()                ! kproma*ngl    Tropopause in Pa
          REAL(DP), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: o3orig_regions   ! kproma*ngl
          REAL(DP), DIMENSION(:,:),   PUBLIC, POINTER :: o3orig_regions2=>NULL()    ! kproma*ngl

          REAL(DP), DIMENSION(:,:,:), PUBLIC, POINTER :: o3prod=>NULL()
          REAL(DP), DIMENSION(:,:,:), PUBLIC, POINTER :: o3losst=>NULL(),o3lossd=>NULL()  ! kproma*nlev t=tracer d=diagnosed
          REAL(DP), DIMENSION(:,:,:), PUBLIC, POINTER :: o3m1=>NULL()              ! kproma*nlev
          REAL(DP), DIMENSION(:,:,:), PUBLIC, POINTER :: co3=>NULL()              ! kproma*nlev

          REAL(DP), DIMENSION(:,:,:), PUBLIC, POINTER :: o3rescale=>NULL() 

          REAL(DP), DIMENSION(:,:,:), PUBLIC, POINTER :: o3_orig=>NULL(),o3_origm1=>NULL()          !nproma*nlev*n_trac_orig

        ! tracer_numbers
          INTEGER, PUBLIC, DIMENSION(n_trac_orig) :: itrac_o3orig

        ! PTR TO OZONE TRACER/CHANNEL OBJECT

          INTEGER               ::  idx_o3losst, idx_o3prod   !   Tracer identifier for tracer_gp channel

          TYPE (t_stream), PUBLIC, POINTER :: o3orig ! the o3orig stream

           ! SUBROUTINES
           PUBLIC :: o3orig_integrate             ! integration of ozone origin tracers
           PUBLIC :: o3orig_initialize
           PUBLIC :: o3orig_new_tracer
           PUBLIC :: o3orig_init_memory           ! allocate memory
           PUBLIC :: o3orig_physc                 ! integrate one time step
           PUBLIC :: o3orig_free_memory
           PUBLIC :: o3orig_read_regions

         CONTAINS

         ! ----------------------------------------------------------------------
           SUBROUTINE o3orig_integrate(kproma,nlev,jrow,           & ! Input  (Grid)
           dtime,                                             & ! Input  (Timestep) 
           !o3_orig_regions,tp,                  & ! Input  (Regions)
           o3_orig_regions,                      & ! Input  (Regions)
                      press,                                             & ! Input  (Regions)
           o3loss,o3prod,o3,o3m1,                             & ! Input  (Chem)
           po3orig,                                            & ! In- Output (O3orig)
           status                                             & ! Output
           )
         ! 2nd try: Take only Prod and O3, O3m1 as known values.
         !        Calculate from that o3loss
             IMPLICIT NONE

             ! I/O
             INTEGER,                   INTENT(IN) :: kproma,nlev,jrow           ! Input
             REAL(DP),                  INTENT(IN) :: dtime                 ! Input [s]
             REAL(DP), DIMENSION(:),    INTENT(IN) :: o3_orig_regions !,tp    ! Input
             REAL(DP), DIMENSION(:,:),  INTENT(IN) :: press                 ! Input [Pa]
             REAL(DP), DIMENSION(:,:),  INTENT(INOUT) :: o3prod           ! Input [molec/cm3/s]
             REAL(DP), DIMENSION(:,:),  INTENT(OUT):: o3loss                 ! Output [molec/cm3/s]
             REAL(DP), DIMENSION(:,:),  INTENT(IN) :: o3, o3m1               ! Input [molec/cm3]
             REAL(DP), DIMENSION(:,:,:),INTENT(INOUT):: po3orig             ! In/Output [molec/cm3]
             INTEGER,                   INTENT(OUT):: status

             ! LOCAL
             REAL (DP), DIMENSION(kproma,nlev)     :: rescaleO3     ! rescaling parameter 
             REAL (DP), DIMENSION(kproma,nlev)     :: totO3_orig    ! total ozone for calculation with the tropopause as boundary 
             
             REAL (DP), DIMENSION(kproma,nlev)     :: prod_err, loss_err  ! prod and loss for error analysis [molec/cm3/timestep]
             REAL (DP), DIMENSION(0:21)            :: prod          ! HARD WIRED !!!!  !sm_01092010
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
                     if (press(i,k).ge.85000._dp) then !NPBL
                             po3orig(i,k,1)=o3m1(i,k)
                     elseif (press(i,k) .lt. 85000._dp .and. press(i,k).ge.20000._dp) then ! NPFT
                             po3orig(i,k,6)=o3m1(i,k)
                     elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NPLS   
                             po3orig(i,k,11)=o3m1(i,k)
                     elseif (press(i,k).lt.3000._dp) then  !NPUS
                             po3orig(i,k,17)=o3m1(i,k)
                     endif
                  elseif (int(o3_orig_regions(i)).eq.2) then
                     if (press(i,k).ge.85000._dp) then !NMBL
                             po3orig(i,k,2)=o3m1(i,k)
                     elseif (press(i,k).lt.85000._dp .and. press(i,k).ge.20000._dp) then !NMFT
                             po3orig(i,k,7)=o3m1(i,k)
                     elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                             po3orig(i,k,12)=o3m1(i,k)
                     elseif (press(i,k).lt.3000._dp) then  !NMUS
                             po3orig(i,k,18)=o3m1(i,k)
                     endif
                  elseif (int(o3_orig_regions(i)).eq.3) then
                     if (press(i,k).ge.85000._dp) then !TRBL
                             po3orig(i,k,3)=o3m1(i,k)
                     elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.10000._dp) then !TRFT
                            po3orig(i,k,8)=o3m1(i,k)
                     elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                             po3orig(i,k,12)=o3m1(i,k)
                     elseif (press(i,k).lt.3000._dp) then  !NMUS
                             po3orig(i,k,18)=o3m1(i,k)
                     endif
                  elseif (int(o3_orig_regions(i)).eq.4) then
                     if (press(i,k).ge.85000._dp) then !TRBL
                             po3orig(i,k,3)=o3m1(i,k)                
                     elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.10000._dp) then !TRFT   
                             po3orig(i,k,8)=o3m1(i,k)
                     elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !TRLS
                             po3orig(i,k,13)=o3m1(i,k)
                     elseif (press(i,k).lt.3000._dp.and.press(i,k).ge.500._dp) then  !TRMS
                             po3orig(i,k,14)=o3m1(i,k)
                     elseif (press(i,k).lt.500._dp) then   !TRUS
                             po3orig(i,k,19)=o3m1(i,k)
                     endif
                  elseif (int(o3_orig_regions(i)).eq.5) then
                     if (press(i,k).ge.85000._dp) then !TRBL
                             po3orig(i,k,3)=o3m1(i,k)
                     elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.10000._dp) then !TRFT
                            po3orig(i,k,8)=o3m1(i,k)
                     elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                             po3orig(i,k,15)=o3m1(i,k)
                     elseif (press(i,k).lt.3000._dp) then  !SMUS
                             po3orig(i,k,20)=o3m1(i,k)
                     endif
                  elseif (int(o3_orig_regions(i)).eq.6) then
                     if (press(i,k).ge.85000._dp) then !SMBL
                             po3orig(i,k,4)=o3m1(i,k)
                     elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.20000._dp) then !SMFT
                             po3orig(i,k,9)=o3m1(i,k)
                     elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                             po3orig(i,k,15)=o3m1(i,k)
                     elseif (press(i,k).lt.3000._dp) then  !SMUS
                             po3orig(i,k,20)=o3m1(i,k)
                     endif
                  elseif (int(o3_orig_regions(i)).eq.7) then
                     if (press(i,k).ge.85000._dp) then !SPBL
                             po3orig(i,k,5)=o3m1(i,k)
                     elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.20000._dp) then !SPFT
                             po3orig(i,k,10)=o3m1(i,k)
                     elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SPLS   
                             po3orig(i,k,16)=o3m1(i,k)
                     elseif (press(i,k).lt.3000._dp) then  !SPUS
                             po3orig(i,k,21)=o3m1(i,k)
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
                        o3rescale(i,k,jrow) = rescaleO3(i,k)
               endif

             enddo
             enddo

        ! Output in case of initialistion
!!$             if (icount1.gt.0) write(message_text,*) 'Grid points initialized: ',icount1
!!$             CALL message('O3orig: ',TRIM(message_text))
!!$             if (icount2.gt.0) write(message_text,*)  'O3orig: Grid points uninitializable: ',icount2
!!$             CALL message('O3orig: ',TRIM(message_text))

             ! if (icount1.gt.0) write(*,*) 'Regions',(o3_orig_regions(i),i=1,kproma)

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
                 if (press(i,k).ge.85000._dp) then !NPBL
                    Prod(1)=o3prod(i,k)
                 elseif (press(i,k).lt.85000._dp .and. press(i,k).ge.20000._dp) then !NPFT
                    Prod(6)=o3prod(i,k)                
                 elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NPLS   
                    Prod(11)=o3prod(i,k)
                 elseif (press(i,k).lt.3000._dp) then  !NPUS
                    Prod(17)=o3prod(i,k)
                 endif
              elseif (int(o3_orig_regions(i)).eq.2) then
                 if (press(i,k).ge.85000._dp) then !NMBL
                    Prod(2)=o3prod(i,k)
                 elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.20000._dp) then !NMFT
                    Prod(7)=o3prod(i,k)
                 elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                    Prod(12)=o3prod(i,k)
                 elseif (press(i,k).lt.3000._dp) then  !NMUS
                    Prod(18)=o3prod(i,k)
                 endif
              elseif (int(o3_orig_regions(i)).eq.3) then
                 if (press(i,k).ge.85000._dp) then !TRBL
                    Prod(3)=o3prod(i,k)
                 elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.10000._dp) then !TRFT
                    Prod(8)=o3prod(i,k)
                 elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !NMLS   
                    Prod(12)=o3prod(i,k)
                 elseif (press(i,k).lt.3000._dp) then  !NMUS
                    Prod(18)=o3prod(i,k)
                 endif
              elseif (int(o3_orig_regions(i)).eq.4) then
                 if (press(i,k).ge.85000._dp) then !TRBL
                    Prod(3)=o3prod(i,k)
                 elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.10000._dp) then !TRFT
                    Prod(8)=o3prod(i,k)
                 elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !TRLS   
                    Prod(13)=o3prod(i,k)
                 elseif (press(i,k).lt.3000._dp.and.press(i,k).ge.500._dp) then  !TRMS
                    Prod(14)=o3prod(i,k)
                 elseif (press(i,k).lt.500._dp) then   !TRUS
                    Prod(19)=o3prod(i,k)
                 endif
              elseif (int(o3_orig_regions(i)).eq.5) then
                 if (press(i,k).ge.85000._dp) then !TRBL
                    Prod(3)=o3prod(i,k)
                 elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.10000._dp) then !TRFT
                    Prod(8)=o3prod(i,k)
                 elseif (press(i,k).lt.10000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                    Prod(15)=o3prod(i,k)
                 elseif (press(i,k).lt.3000._dp) then  !SMUS
                    Prod(20)=o3prod(i,k)
                 endif
              elseif (int(o3_orig_regions(i)).eq.6) then
                 if (press(i,k).ge.85000._dp) then !SMBL
                    Prod(4)=o3prod(i,k)
                 elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.20000._dp) then !SMFT
                    Prod(9)=o3prod(i,k)
                 elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SMLS   
                    Prod(15)=o3prod(i,k)
                 elseif (press(i,k).lt.3000._dp) then  !SMUS
                    Prod(20)=o3prod(i,k)
                 endif
              elseif (int(o3_orig_regions(i)).eq.7) then
                 if (press(i,k).ge.85000._dp) then !SPBL
                    Prod(5)=o3prod(i,k)
                 elseif (press(i,k).lt.85000._dp.and.press(i,k).ge.20000._dp) then !SPFT
                    Prod(10)=o3prod(i,k)
                 elseif (press(i,k).lt.20000._dp.and.press(i,k).ge.3000._dp) then !SPLS   
                    Prod(16)=o3prod(i,k)
                 elseif (press(i,k).lt.3000._dp) then  !SPUS
                    Prod(21)=o3prod(i,k)
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

        ! ----------------------------------------------------------------------
          SUBROUTINE o3orig_initialize

            ! Reads namelist *o3origctl*
            ! *o3orig_initialize* is called from *call_init_submodels*, src/call_submodels.f90

           USE mo_control,            ONLY: nlev, nproma
            USE mo_doctor,             ONLY: nout, nerr
            USE mo_exception,       ONLY: finish
            USE mo_mpi,                 ONLY: p_io, p_parallel, p_parallel_io, p_bcast
            USE mo_namelist,         ONLY: position_nml, nnml, POSITIONED, MISSING, &
                                             LENGTH_ERROR, READ_ERROR

            IMPLICIT NONE

            ! Local
            INTEGER     :: ierr  ! error return value from position_nml
            INTEGER     :: i_trac

            NAMELIST /o3origctl/  i_integrate,  &       ! 1: Euler backward   2: euler forward
                                  i_trac_orig, &        ! Number of tracers
                                  sn_o3orig,   &        ! Short name of the tracers
                                  l_err                 ! Error diagnostic on?


            IF (p_parallel_io) THEN
              CALL position_nml ('O3ORIGCTL', status=ierr)
              SELECT CASE (ierr)
              CASE (POSITIONED)
                READ (nnml, o3origctl)
              CASE(MISSING)
                 CALL finish ('o3orig_initialize', &
                      'namelist o3origctl is missing in namelist.echam')
              CASE(LENGTH_ERROR)
                 CALL finish ('o3orig_initialize', &
                      'namelist o3origctl is longer than in the definition')
              CASE(READ_ERROR)
                 CALL finish ('o3orig_initialize', 'general read error when reading namelist')
              END SELECT
            ENDIF

            IF (p_parallel) THEN

              CALL p_bcast(i_integrate, p_io)
              CALL p_bcast(i_trac_orig, p_io)
              do i_trac=1,i_trac_orig
                CALL p_bcast(sn_o3orig(i_trac), p_io)
              enddo
              CALL p_bcast(l_err, p_io)

              if (i_trac_orig.gt.n_trac_orig) then
                CALL finish('o3orig_initialize','Number of regions exceeds limit')
              endif

            ENDIF

            IF (p_parallel_io) THEN
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

              do i_trac=1,i_trac_orig
                 write(*,*) 'O3ORIG-Tracer',i_trac,': ',sn_o3orig(i_trac)
              enddo

            END IF

          END SUBROUTINE o3orig_initialize
        ! ------------------------------------------------------------------------

        ! ------------------------------------------------------------------------
          SUBROUTINE o3orig_new_tracer
            ! Defines additional tracers for O3 orig diagnostics 

            ! *o3orig_new_tracer* is called from *call_request_tracer*,
            ! src/call_submodels.f90.

            USE mo_advection,       ONLY: iadvec   ! selected advection scheme
            USE mo_tracdef,         ONLY: trlist, t_trlist, t_trinfo
            USE mo_tracer,          ONLY: ntrac, & ! number of tracers
                                          new_tracer, & ! request resources for a tracer
                                          GAS, & ! phase indicators
                                          OFF, ON, &
                                          trlist, RESTART, INITIAL, CONSTANT

            IMPLICIT NONE

            ! LOCAL
            INTEGER :: i_trac, idt
            INTEGER :: trac_tran  ! switch for transport in ECHAM5
            INTEGER :: trac_vdiff ! switch for vertical diffusion in ECHAM5
            INTEGER :: trac_conv  ! switch for convection in ECHAM5

            ! Executable statements:

            trac_tran=iadvec
            trac_vdiff=ON
            trac_conv=ON

        !      write(*,*) 'start to define',i_trac_orig,'tracers for O3ORIG'
               do i_trac=1,i_trac_orig
         
        !        write(*,*) "i_trac=",i_trac,"Name:",TRIM(sn_o3orig(i_trac)),"#",status,sum_status
                 CALL new_tracer(TRIM(sn_o3orig(i_trac)),'MEZON'    &
                       , units='mol/mol', idx =  itrac_o3orig(i_trac)                 &
                       , longname='O3orig tracer - '//TRIM(sn_o3orig(i_trac))         &
                       , nwrite=OFF &
                       , ninit=RESTART+CONSTANT, vini = 0.0_dp, nrerun=ON &
                       , code = 100+i_trac, table = 131  &
                       , ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv &
                       , nphase=GAS  )

               end do

          END SUBROUTINE o3orig_new_tracer
        ! ------------------------------------------------------------------------

        ! ------------------------------------------------------------------------
          SUBROUTINE o3orig_init_memory

            ! Allocates streams for o3orig

            ! *o3orig_init_memory* is called from *call_init_submodel_memory*,
            ! src/call_submodels.f90.

          USE mo_control,      ONLY: nlev, nproma
          USE mo_doctor,        ONLY: nout
          USE mo_kind,          ONLY: dp
          USE mo_linked_list,   ONLY: HYBRID, NETCDF
          USE mo_mpi,           ONLY: p_parallel_io

          !--- Create new stream:

            CALL new_stream (o3orig ,'o3orig')

          ! add entries for geopotential and log surface pressure (used by the
          ! afterburner):
            CALL add_stream_reference (o3orig, 'geosp'   ,'g3b'   ,lpost=.FALSE.)
            CALL add_stream_reference (o3orig, 'lsp'     ,'sp'    ,lpost=.FALSE.)

            CALL default_stream_setting (o3orig, units ='', &
                                               lpost = .FALSE., &
                                               lrerun    = .FALSE. , &
                                               leveltype = HYBRID , &
                                               table     = 199,     &
                                               laccu     = .FALSE.)

            CALL add_stream_element(o3orig      &
                 , 'O3Losst'                        &
                 , o3losst                          &
                 , longname='Ozone loss from tracer' &
                 , units = 'mol/mol/s'              &
                 , code=221 &
                 )

            IF (p_parallel_io) WRITE(*,*) ' ... O3Loss from tracer added to stream o3orig'

            CALL add_stream_element(o3orig      &
                 , 'O3Lossd'                        &
                 , o3lossd                          &
                 , longname='Ozone loss diagnosed  ' &
                 , units = 'mol/mol/s'              &
                 , code=222 &
                 )

            IF (p_parallel_io) WRITE(*,*) ' ... O3Loss diagnosed added to stream o3orig'

            CALL add_stream_element(o3orig      &
                 , 'O3Prod'                         &
                 , o3prod                           &
                 , longname='Ozone production     ' &
                 , units = 'mol/mol/s'              &
                 , code=224 &
                 )

            IF (p_parallel_io) WRITE(*,*) ' ... O3Prod added to stream o3orig'


            CALL add_stream_element(o3orig      &
                 , 'O3M1'                           &
                 , o3m1                             &
                 , longname='Pre-Chem ozone       ' &
                 , units = 'mol/mol'                &
                 , code=225 &
                 )

            IF (p_parallel_io) WRITE(*,*) ' ... O3M1   added to stream o3orig'


            CALL add_stream_element(o3orig      &
                 , 'O3'                             &
                 , o3                               &
                 , longname='Post-Chem ozone      ' &
                 , units = 'mol/mol'                &
                 , code=226 &
                 )

            IF (p_parallel_io) WRITE(*,*) ' ... O3     added to stream o3orig'



            CALL add_stream_element(o3orig      &
                 , 'cO3'                            &
                 , co3                              &
                 , longname='Chem ozone tendency  ' &
                 , units = 'mol/mol/s'              &
                 , code=227 &
                 )

            IF (p_parallel_io) WRITE(*,*) ' ... cO3    added to stream o3orig'


            CALL add_stream_element(o3orig      &
                 , 'Reg'                            &
                 , o3orig_regions2                  &
                 , longname='Regions              ' &
                 , units = 'none   '                &
                 , code=228 &
                 , lpost = .FALSE. &
                 )
            
            IF (p_parallel_io) WRITE(*,*) ' ... Reg    added to stream o3orig'

             CALL add_stream_element(o3orig      &
                 , 'Rescale'                            &
                 , o3rescale                 &
                 , longname='Rescaling factor O3              ' &
                 , units = '-'                &
                 , code=231 &
                 )

            ALLOCATE(o3_orig(nproma,nlev,n_trac_orig))
            ALLOCATE(o3_origm1(nproma,nlev,n_trac_orig))

  END SUBROUTINE o3orig_init_memory
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE o3orig_physc(i_o3orig_ctrl,kproma,kbdim,jrow,ntrac,pxtm1,pxtte,press)

    ! ECHAM5
    USE mo_control,            ONLY: ngl,nlev,nlon
    USE mo_decomposition,      ONLY: dcl => local_decomposition
    USE mo_exception,          ONLY: finish
    USE mo_time_control,       ONLY: time_step_len
    USE mo_socol_tracers,      ONLY: idt_o3

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)            :: i_o3orig_ctrl, &             ! 1: pre mecca; 2: post Mecca
                                      jrow, kbdim, kproma, ntrac
    INTEGER     :: status

    REAL(DP) :: pxtte(kbdim,nlev,ntrac), pxtm1(kbdim,nlev,ntrac)

    REAL(DP) :: press(kbdim,nlev)

    REAL(DP), DIMENSION (kbdim,nlev)          ::  o3l, o3p

!   REAL(DP), DIMENSION(kproma,nlev,n_trac_orig)     :: o3origM1  !sm_18082010 sm_02092010

    INTEGER                               :: jt
    INTEGER ::  idt ! um_ak_20110719

    status=0

    SELECT CASE(i_o3orig_ctrl)
    CASE(1)                    ! Pre-MECCA
                               ! Set ozone and tendencies

       o3m1(1:kproma,:,jrow) = pxtm1(1:kproma,:,idt_o3) & 
                 + pxtte(1:kproma,:,idt_o3) * time_step_len

        o3orig_regions2(1:kproma,jrow) = o3orig_regions(1:kproma,jrow)

    CASE(2)                    ! Post-MECCA
     
        o3(1:kproma,:,jrow)   = pxtm1(1:kproma,:,idt_o3) & 
                   + pxtte(1:kproma,:,idt_o3) * time_step_len
        co3(1:kproma,:,jrow)=(o3(1:kproma,:,jrow)-o3m1(1:kproma,:,jrow))/time_step_len  ! per second!
        do jt=1,i_trac_orig
           o3_orig(1:kproma,:,jt)  = pxtm1(1:kproma,:,itrac_o3orig(jt))  &
              + pxtte(1:kproma,:,itrac_o3orig(jt)) * time_step_len
        enddo
       
        o3_origM1(:,:,:)=o3_orig(:,:,:)

        !!o3losst(1:kproma,:,jrow)   = pxtte(1:kproma,:,idx_O3losst) 
        !!o3prod(1:kproma,:,jrow)    = pxtte(1:kproma,:,idx_O3prod) 

        !!o3l(1:kproma,:)    = pxtte(1:kproma,:,idx_O3losst)
        !!o3p(1:kproma,:)    = pxtte(1:kproma,:,idx_O3prod)

        o3l(1:kproma,:)    = o3losst(1:kproma,:,jrow)
        o3p(1:kproma,:)   = o3prod(1:kproma,:,jrow)

        CALL o3orig_integrate(kproma,nlev,jrow,               & ! Input  (Grid)
             time_step_len,                              & ! Input  (Timestep)
             !!o3orig_regions(1:kproma,jrow),tp(1:kproma,jrow), & ! Input  (Regions; tropopause pressure)
             o3orig_regions(1:kproma,jrow), & ! Input  (Regions; tropopause pressure)
             press(1:kproma,:),                          & ! Input  (Pressure)
             !!press_3d(1:kproma,_RI_YC_),                 & ! Input  (Pressure)
             o3l(1:kproma,:),o3p(1:kproma,:),            & ! Input  (Chem)
             o3(:,:,jrow),o3m1(:,:,jrow),            & ! Input  (Chem)
             !!o3(1:kproma,_RI_YC_),o3m1(1:kproma,_RI_YC_),& ! Input  (Chem)
             o3_orig(1:kproma,:,:),                      & ! In/Output (O3Orig)
             status)
        
        o3lossd(1:kproma,:,jrow)    = o3l(1:kproma,:)
        o3prod(1:kproma,:,jrow)     = o3p(1:kproma,:)

         do jt=1,i_trac_orig
           pxtte(1:kproma,:,itrac_o3orig(jt))=pxtte(1:kproma,:,itrac_o3orig(jt))  &
             + (o3_orig(1:kproma,:,jt)-o3_origM1(1:kproma,:,jt)) / time_step_len
         enddo
 
    CASE DEFAULT
         status=1

    IF (status/=0) THEN
       CALL finish(&
         'o3orig_physc', 'incorrect value for i_o3orig_ctrl')
    ENDIF
    END SELECT

  END SUBROUTINE o3orig_physc
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE o3orig_free_memory

  implicit none

  DEALLOCATE(o3_orig)
  DEALLOCATE(o3_origm1)

  END SUBROUTINE o3orig_free_memory
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE o3orig_read_regions

    ! reads in information about horizontal dependence of o3orig regions
    ! *o3orig_read_regions* is called from stepon

    USE mo_control,       ONLY: no3orig
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_decomposition,      ONLY: lc => local_decomposition, global_decomposition

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:)
    REAL(dp), POINTER :: gl_orig(:,:)

    LOGICAL       :: lex
    INTEGER       :: nvarid

    IF (p_parallel_io) THEN

       CALL message('','O3 orig: read in latitude bands.')
       INQUIRE (no3orig, exist=lex)
       WRITE(message_text,*) 'lex: ', lex
       CALL message('o3orig_read_regions',message_text)
       IF (lex) THEN
          o3nc1%format = NETCDF
          CALL IO_open_unit (no3orig, o3nc1, IO_READ)
       ELSE
          CALL finish ('o3orig_read_regions', 'Could not open input file')
       ENDIF

    ENDIF

    !     Allocate memory for o3orig_regions per PE

    IF (.NOT. ALLOCATED(o3orig_regions)) ALLOCATE (o3orig_regions(lc%nproma, lc%ngpblks))

    !     Read file
    IF (p_parallel_io) THEN

       !     Allocate memory

       ALLOCATE (zin(lc%nlon,lc%nlat))

       CALL IO_INQ_VARID (o3nc1%file_id, 'orig', nvarid)
       CALL IO_GET_VAR_DOUBLE (o3nc1%file_id, nvarid, zin(:,:))

    END IF

    NULLIFY (gl_orig)
    IF (p_pe == p_io) gl_orig => zin(:,:)
    CALL scatter_gp (gl_orig, o3orig_regions(:,:), global_decomposition)

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       CALL IO_close(o3nc1)

    END IF

  END SUBROUTINE o3orig_read_regions

! **********************************************************************
END MODULE mo_socol_o3orig
! **********************************************************************

MODULE mo_socol_scav

  ! *mo_socol_scav* contains subroutines to calculate and diagnose
  ! wet deposition of several chemical species
  !
  ! Andrea Stenke, ETH Zurich, April 2010

  USE mo_kind,          ONLY: dp      
  USE mo_linked_list,   ONLY: SURFACE, HYBRID, NETCDF
  USE mo_memory_base,   ONLY: new_stream, add_stream_element,  &
                              default_stream_setting, add_stream_reference, &
                              delete_stream, t_stream
  USE mo_time_control,  ONLY: delta_time, lstart

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: socol_scav              ! calculate scavenging of different species
  PUBLIC :: cond_scav               ! calculate cloud and precipitation parameters for scavenging
  PUBLIC :: init_scav               ! initialize scav
  PUBLIC :: construct_stream_scav   ! construct the stream
  PUBLIC :: destruct_stream_scav    ! destruct the stream
  PUBLIC :: init_stream_scav        ! initialize the stream
  PUBLIC :: accumulate_stream_scav  ! accumulate stream elements

  TYPE (t_stream), PUBLIC, POINTER :: scav     !the scav stream

  INTEGER, PUBLIC, PARAMETER :: nspec_scav = 5

  INTEGER, PUBLIC, DIMENSION(nspec_scav) :: idx_spec

  REAL(dp), PUBLIC, DIMENSION(nspec_scav) :: zdg, & ! gas-phase diffusivity of species
                                             dht, & ! temperature depedence of henry coefficients
                                             khcp   ! Henry coefficient representing solubility [M/atm]
                                                    ! (i.e. aqueous-phase composition divided by gas-phase composition) 

  REAL(dp), PUBLIC, POINTER :: wd_ch3o2h(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_ch2o(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_h2o2(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_hno3(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_o3(:,:)

  REAL(dp), PUBLIC, POINTER :: wd_ch3o2h_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_ch2o_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_h2o2_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_hno3_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wd_o3_d(:,:)

  REAL(dp), PUBLIC, POINTER :: wdcv_ch3o2h(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_ch2o(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_h2o2(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_hno3(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_o3(:,:)

  REAL(dp), PUBLIC, POINTER :: wdcv_ch3o2h_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_ch2o_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_h2o2_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_hno3_d(:,:)
  REAL(dp), PUBLIC, POINTER :: wdcv_o3_d(:,:)

  REAL(dp), PUBLIC, POINTER :: wohno3(:,:,:)
  REAL(dp), PUBLIC, POINTER :: wohno3_d(:,:,:)  

  DATA zdg /0.131_dp, 0.164_dp, 0.184_dp, 0.136_dp, 0.164_dp/
  DATA dht /5322._dp, 6425._dp, 6338._dp, 8694._dp, 2560._dp/
  DATA khcp /3.0E2_dp, 7.0E3_dp, 1.E5_dp, 2.1E5_dp, 1.2E-02_dp/

!-------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE init_scav

    ! Initializes idx_spec
    ! called from call_submodels.f90/call_request_tracer

    ! Andrea Stenke, ETH Zurich, April 2010

    USE mo_socol_tracers, ONLY: idt_ch3o2h, idt_ch2o, idt_h2o2, idt_hno3, idt_o3

    idx_spec(1) = idt_ch3o2h
    idx_spec(2) = idt_ch2o   
    idx_spec(3) = idt_h2o2 
    idx_spec(4) = idt_hno3 
    idx_spec(5) = idt_o3  

  END SUBROUTINE init_scav

!-------------------------------------------------------------------------------
  SUBROUTINE socol_scav(krow, kproma, kbdim, klev, ktrac, p, ph, ptvm1, pxtm1, pxtte, &
       ptm1, paclc, zmratep, zfprec, zfevap, zmlwc, xn3depcv, cvdprec, kconbot, totliqcv)

    ! Description:
    !
    ! Calculation of scavenging of soluble gases
    ! based on Roelofs and Lelieveld, JGR, 100 (D10), 20983-20998, 1995
    ! original code CLSCAV.F ECHAM4.L39(DLR)/CHEM

    ! *socol_scav* is called from chem2
    
    ! Andrea Stenke, ETH Zurich, April 2010
  
   
  USE mo_constants,               ONLY: api, argas, avo, amd, g, rd
  USE mo_kind,                    ONLY: dp
  USE mo_time_control,            ONLY: delta_time, time_step_len
  USE mo_socol_tracers,           ONLY: trac_chemspec, idt_ch3o2h, idt_ch2o, idt_h2o2, &
                                        idt_hno3, idt_o3

  IMPLICIT NONE

  ! Subroutine arguments:
  INTEGER, INTENT(in) :: krow, kproma, kbdim, klev, ktrac
  INTEGER, INTENT(in) :: kconbot(kbdim) 
  REAL(dp), INTENT(in) :: p(kbdim,klev)                        ! pressure at full levels
  REAL(dp), INTENT(in) :: ph(kbdim,klev+1)                     ! pressure at half levels 
  REAL(dp), INTENT(in) :: ptvm1(kbdim,klev)                    ! virtual temperature
  REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac)              ! mixing ratios of chemical species
  REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)           ! tendencies of chemical species
  REAL(dp), INTENT(in) :: ptm1(kbdim,klev)                     ! temperature
  REAL(dp), INTENT(inout) :: paclc(kbdim,klev)                    ! cloud cover
  REAL(dp), INTENT(in) :: xn3depcv(kbdim,klev,nspec_scav), cvdprec(kbdim,klev)
  REAL(dp), INTENT(in) :: zfprec(kbdim,klev)                   ! Flux of precip [kg/m猖]
  REAL(dp), INTENT(in) :: zfevap(kbdim,klev)                   ! total flux evaporating within one timestep [kg/m猖]
  REAL(dp), INTENT(inout) :: zmratep(kbdim,klev)               ! precip formation rate [kg/m連]
  REAL(dp), INTENT(in) :: zmlwc(kbdim,klev)                    ! liquid water content
  REAL(dp), INTENT(in) :: totliqcv(kbdim,klev)                   ! liquid water content, convective

  ! Local variables:

  INTEGER, DIMENSION(kbdim) :: kcltop
  INTEGER :: jl, jk, jt 
                                      
  REAL(dp), PARAMETER :: zmoltom2s=avo*1000._dp/amd
  REAL(dp), PARAMETER :: zfac=avo/amd
  REAL(dp), PARAMETER :: cvcover = 0.05_dp ! ECHAM-assumed convective cloud-cover
  REAL(dp), PARAMETER :: dgair = 0.133_dp, dghno3 = 0.136_dp
  REAL(dp), DIMENSION(kbdim,klev) :: zmtof
                               
  REAL(dp), DIMENSION(kbdim,nspec_scav) :: wd, wdcv
  REAL(dp), DIMENSION(kbdim,klev) :: zrho  ! air density
  REAL(dp), DIMENSION(kbdim,klev) :: hno3_old
  REAL(dp), DIMENSION(kbdim,klev,ktrac) :: zm
  REAL(dp), DIMENSION(kbdim) :: raincv, oldcov
  REAL(dp) :: zdp ! delta_p
  REAL(dp) :: beta, betaex ! scavenging coefficient
  REAL(dp) :: zkg          ! mass transfer coefficient
  REAL(dp) :: rdrad, &     ! average cloud drop size
              rn,    &     ! cloud droplet number concentration
              cl2rain, zinflux, zclear, zeop, rflx, rlwc, rdvol, ru, &
              znre, znsc, znsh, &
              zmxcov, washfrac, dwd, &
              zl, d, dl, fac1, &
              totliq, henper, prefracl

  ! Executable statements:

  oldcov(:) = 0._dp
  kcltop(:) = klev
  raincv(:) = 0._dp
  wd(:,:) = 0._dp
  wdcv(:,:) = 0._dp    
                                                                           
  DO jl = 1, kproma
     DO jk = 1, klev
        
!!$        zrho(jl,jk) = p(jl,jk) / (rd * ptvm1(jl,jk)) ! air density
        zrho(jl,jk) = p(jl,jk) * amd * 1.E-6_dp / (argas * ptvm1(jl,jk)) ! air density
        zdp = ph(jl,jk+1) - ph(jl,jk) ! delta_p
        zmtof(jl,jk) = zdp / (delta_time * g * zrho(jl,jk) * 1.E-3_dp)  ! mass-to-flux (cm3 m-2 s-1)  
        zm(jl,jk,:) = (pxtm1(jl,jk,:) + pxtte(jl,jk,:) * delta_time) * zfac * zrho(jl,jk) ! VMR -> molec/cm3

        hno3_old(jl,jk) = zm(jl,jk,idt_hno3) 

        IF (paclc(jl,jk) .LT. 1.E-5_dp .OR. zmlwc(jl,jk) .LT. 1.E-10_dp) THEN
           zmratep(jl,jk) = 0._dp
           paclc(jl,jk) = 0._dp
        END IF

     END DO
  END DO

  DO jk = 1, klev                                            
     
     ! 1. wet deposition due to stratiform precipitation 

     !   Below cloud scavenging:                                                
     !   adds to total wet deposition flux                                     
     !   Note that ZFPREC includes the additional flux in the cloud!       
     
     DO jl = 1, kproma  
                                         
        cl2rain = zmratep(jl,jk) * delta_time * zmtof(jl,jk) * 1.E-6_dp                     
        zinflux = MAX(zfprec(jl,jk) - cl2rain, 0._dp) 
                         
        IF (zinflux.LT.1.E-15_dp) kcltop(jl) = klev                     
        IF ((zinflux - zfevap(jl,jk)).GT.1.E-15_dp .AND. jk.GT.kcltop(jl)) THEN  
           
           ! calculate rain parameters (Kumar, 1985)                               
           rflx = zinflux / oldcov(jl) * 3600._dp ! Intensity of rain entering the level [mm/h]  
           rlwc = 72._dp * rflx**0.88_dp ! effective liquid water content               
           rdrad = 0.3659_dp * rflx**0.21_dp ! effective rain drop radius  
           rdvol = 4._dp * api / 3._dp * rdrad**3._dp ! rain drop volume         
           rn = rlwc / rdvol                                                      
           ru = 9.58_dp * (1._dp - EXP(-(rdrad / 0.885_dp)**1.147_dp))           
            
           DO jt = 1, nspec_scav

              ! calculate mass transfer coefficient(Durham et. al., 1981)       
              znre = SQRT(20._dp * rdrad * rd / dgair)                                       
              znsc = dgair / zdg(jt)       
!!$              ZNSC=DGAIR/DGHNO3
              znsh = 1._dp + 0.3_dp * znre * (znsc**(1._dp / 3._dp))                                   
              zkg = 10._dp * zdg(jt) / rdrad * znsh                                      
              
              ! calculate washout (Kumar, 1985); account for cloud overlap           
              beta = 4._dp * api * rdrad**2._dp * rn * zkg *1.E-8_dp                                
              betaex = EXP(-1._dp * beta * delta_time) - 1._dp

              ! take into account solubility of gas
              totliq = zmlwc(jl,jk) * zrho(jl,jk) * 1.E3_dp ! totliq [l/m設, zmlwc [kg/kg], zrho [g/cm設
              henper = (4.088E-2_dp / khcp(jt)) * & ! T-dependence of henry coef
                   EXP(-dht(jt) * (1._dp/ptm1(jl,jk) - 1._dp/298._dp)) 
              prefracl = 1.E-3_dp / (henper / totliq + 1.E-3_dp)
!!$
              d = zm(jl,jk,idx_spec(jt)) * betaex * prefracl   
!!$              d = zm(jl,jk,idx_spec(jt)) * betaex    
                       
              zmxcov = MIN(oldcov(jl), paclc(jl,jk))                             
              washfrac = oldcov(jl) - zmxcov
              d = d * washfrac

              zm(jl,jk,idx_spec(jt)) = zm(jl,jk,idx_spec(jt)) + d

              dwd = d * zmtof(jl,jk)
              wd(jl,jt) = wd(jl,jt) - dwd

           END DO ! jt
        ENDIF
     END DO ! jl
                                                                 
     ! Re-evaporation of incoming flux                                        
                                                                          
      DO jl = 1, kproma  
         cl2rain = zmratep(jl,jk) * delta_time * zmtof(jl,jk) * 1.E-6_dp                   
         zinflux = MAX(zfprec(jl,jk) - cl2rain, 0._dp)         
         zclear = 1._dp - paclc(JL,JK) 
                                         
         IF (zclear.GT.1.E-20_dp .AND. zinflux.GT.1.E-20_dp) THEN 
                     
            zeop = zfevap(jl,jk) / zinflux
            zeop = MAX(0._dp,zeop)
            zeop = MIN(1._dp,zeop)  

            DO jt = 1, nspec_scav

               d = zeop*wd(jl,jt) / zmtof(jl,jk)
               wd(jl,jt) = (1._dp - zeop) * wd(jl,jt)

               zm(jl,jk,idx_spec(jt)) = zm(jl,jk,idx_spec(jt)) + d
               
            END DO ! jt
         ENDIF
      END DO ! jl
                                     
      ! in-cloud precipitation formation
      ! Transfer from cloud to rain water:                                     
      ! adds to wet deposition flux                                           
                                                                          
      DO jl = 1, kproma
                             
         IF (zmratep(jl,jk) .GT. 1.E-20_dp) THEN  
            
            IF(jk.LT.kcltop(jl)) kcltop(jl) = jk
            
            DO jt = 1, nspec_scav
               
               zl = paclc(jl,jk) * zm(jl,jk,idx_spec(jt))
               fac1 = zl * zmratep(jl,jk) / (zrho(jl,jk) * zmlwc(jl,jk))

               ! take into account solubility of gas
               totliq = zmlwc(jl,jk) * zrho(jl,jk) * 1.E3_dp ! totliq [l/m設, zmlwc [kg/kg], zrho [g/cm設
               henper = (4.088E-2_dp / khcp(jt)) * & ! T-dependence of henry coef
                    EXP(-dht(jt) * (1._dp/ptm1(jl,jk) - 1._dp/298._dp)) 
               prefracl = 1.E-3_dp / (henper / totliq + 1.E-3_dp)

               dl = fac1 * delta_time * 1.E-3_dp * prefracl
!!$               dl = fac1 * delta_time * 1.E-3_dp

               dl = MAX(0._dp, dl)
               dl = MIN(zm(jl,jk,idx_spec(jt)), dl)
               wd(jl,jt) = wd(jl,jt) + dl * zmtof(jl,jk)

               zm(jl,jk,idx_spec(jt)) = zm(jl,jk,idx_spec(jt)) - dl  
              
            END DO ! jt
        
            ! Set cover for gridbox below                                          
            oldcov(jl) = paclc(jl,jk)
            
         ENDIF
      END DO ! jl
      
      ! 2. wet deposition due to convective precipitation 
  
      ! Below cloud scavenging                                                  
                                                                          
      DO jl = 1, kproma  
         IF (kconbot(jl).GT.0) THEN            
            IF (jk.GT.kconbot(jl) .AND. raincv(jl).GT.0._dp) THEN
               
               ! calculate rain parameters
               rflx = raincv(jl) / cvcover*3600._dp
               rlwc = 72._dp * rflx**0.88_dp ! effective liquid water content                
               rdrad = 0.3659_dp * rflx**0.21_dp ! effective rain drop radius               
               rdvol = 4._dp * api / 3._dp * rdrad**3._dp ! rain drop volume             
               rn = rlwc / rdvol                                                      
               ru = 9.58_dp * (1._dp - EXP(-(rdrad / 0.885_dp)**1.147_dp))        
            
               DO jt = 1, nspec_scav

                  ! calculate mass transfer coefficient (Durham et. al., 1981)       
                  znre = SQRT(20._dp * rdrad * rd / dgair)                                       
                  znsc = dgair / zdg(jt) 
!!$                  ZNSC = DGAIR / DGHNO3                            
                  znsh = 1._dp + 0.3_dp * znre * (znsc**(1._dp / 3._dp))    
                  zkg = 10._dp * zdg(jt) / rdrad * znsh  
                                        
                  beta = 4._dp * api * rdrad**2._dp * rn * zkg *1.E-8_dp                                
                  betaex = EXP(-1._dp * beta * delta_time) - 1._dp

                  ! take into account solubility of gas
                  totliq = zmlwc(jl,jk) * zrho(jl,jk) * 1.E3_dp ! totliq [l/m設, zmlwc [kg/kg], zrho [g/cm設
                  henper = (4.088E-2_dp / khcp(jt)) * & ! T-dependence of henry coef
                       EXP(-dht(jt) * (1._dp/ptm1(jl,jk) - 1._dp/298._dp)) 
                  prefracl = 1.E-3_dp / (henper / totliq + 1.E-3_dp)
                  
                  d = zm(jl,jk,idx_spec(jt)) * betaex * prefracl    
!!$                  d = zm(jl,jk,idx_spec(jt)) * betaex

                  d = d * cvcover
                  zm(jl,jk,idx_spec(jt)) = zm(jl,jk,idx_spec(jt)) + d
                  dwd = d * zmtof(jl,jk)
                  wdcv(jl,jt) = wdcv(jl,jt) - dwd
     
               END DO
            ENDIF          
         ENDIF
      END DO
      
      !  Below cloud evaporation                                                 
      
      DO jl = 1, kproma                                               
         IF (kconbot(jl).GT.0 .AND. raincv(jl).GT.0._dp) THEN
            IF (jk.GT.kconbot(jl) .AND. cvdprec(jl,jk).LT.0._dp) THEN   
            
               zeop = -1._dp * cvdprec(jl,jk) / raincv(jl)
               zeop = MAX(0._dp,zeop)
               zeop = MIN(1._dp,zeop)  

               DO jt = 1, nspec_scav
                  
                  d = zeop * wdcv(jl,jt) / zmtof(jl,jk)
                  wdcv(jl,jt) = (1._dp - zeop) * wdcv(jl,jt)
                  zm(jl,jk,idx_spec(jt)) = zm(jl,jk,idx_spec(jt)) + d

               END DO               
            ENDIF
         ENDIF
      END DO
                                                           
      ! Add up the in-cloud scavenging calculated in Tiedtke                                     
                     
      DO jl = 1, kproma                                                     
         IF (kconbot(jl).GT.0 ) THEN    
            
            raincv(jl) = raincv(jl) + cvdprec(jl,jk)
                               
            IF (jk.LE.kconbot(jl)) THEN                                        
               
               DO jt = 1, nspec_scav

                  wdcv(jl,jt) = wdcv(jl,jt) + xn3depcv(jl,jk,jt) * zmoltom2s

               END DO                   
            ENDIF
         ENDIF
      END DO

      ! store HNO3 washout
      do jl = 1, kproma
         
         wohno3_d(jl,jk,krow) = hno3_old(jl,jk) - zm(jl,jk,idt_hno3)

      end do

      ! Update Tendencies, molc/cm3 -> VMR

      DO jl = 1, kproma
         DO jt = 1, nspec_scav

            IF (idx_spec(jt) == idt_ch3o2h) THEN             
               pxtte(jl,jk,idt_ch3o2h) = (zm(jl,jk,idt_ch3o2h) / &
                    (zfac * zrho(jl,jk)) - pxtm1(jl,jk,idt_ch3o2h)) / delta_time
            ELSEIF (idx_spec(jt) == idt_ch2o) THEN             
               pxtte(jl,jk,idt_ch2o) = (zm(jl,jk,idt_ch2o) / &
                    (zfac * zrho(jl,jk)) - pxtm1(jl,jk,idt_ch2o)) / delta_time
            ELSEIF (idx_spec(jt) == idt_h2o2) THEN             
               pxtte(jl,jk,idt_h2o2) = (zm(jl,jk,idt_h2o2) / &
                    (zfac * zrho(jl,jk)) - pxtm1(jl,jk,idt_h2o2)) / delta_time
            ELSEIF (idx_spec(jt) == idt_hno3) THEN             
               pxtte(jl,jk,idt_hno3) = (zm(jl,jk,idt_hno3) / &
                    (zfac * zrho(jl,jk)) - pxtm1(jl,jk,idt_hno3)) / delta_time
            ELSEIF (idx_spec(jt) == idt_o3) THEN             
               pxtte(jl,jk,idt_o3) = ( zm(jl,jk,idt_o3) / &
                    (zfac * zrho(jl,jk)) - pxtm1(jl,jk,idt_o3) ) / delta_time
            END IF

         ENDDO
      END DO
                                   
   END DO ! jk

   ! Diagnostics

   DO jl = 1, kproma
      DO jt = 1, nspec_scav
         
         IF (idx_spec(jt) == idt_ch3o2h) THEN
            wd_ch3o2h_d(jl,krow) = wd(jl,jt) / zfac ! kg/m猖 
            wdcv_ch3o2h_d(jl,krow) = wdcv(jl,jt) / zfac
         ELSEIF (idx_spec(jt) == idt_ch2o) THEN
            wd_ch2o_d(jl,krow) = wd(jl,jt) / zfac 
            wdcv_ch2o_d(jl,krow) = wdcv(jl,jt) / zfac 
         ELSEIF (idx_spec(jt) == idt_h2o2) THEN
            wd_h2o2_d(jl,krow) = wd(jl,jt) / zfac
            wdcv_h2o2_d(jl,krow) = wdcv(jl,jt) / zfac    
         ELSEIF (idx_spec(jt) == idt_hno3) THEN
            wd_hno3_d(jl,krow) = wd(jl,jt) / zfac          
            wdcv_hno3_d(jl,krow) = wdcv(jl,jt) / zfac 
         ELSEIF (idx_spec(jt) == idt_o3) THEN
            wd_o3_d(jl,krow) = wd(jl,jt) / zfac
            wdcv_o3_d(jl,krow) = wdcv(jl,jt) / zfac 
         END IF
      END DO
   END DO

 END SUBROUTINE socol_scav

!-------------------------------------------------------------------------------
SUBROUTINE cond_scav   ( kproma, kbdim, KTDIA,KLEV,KLEVP1 & 
!!$C----------------------------------------------------------------------   
!!$C - INPUT  2D .                                                           
     , KLAB,     PACLCACM, PAPHM1, PAPHP1 &          
     , PAPM1,    PAPP1,    PGEOM1, PQM1,   PTM1 &        
     , PXLM1,    PXIM1,    PXTEC &                                  
!!$C - INPUT  1D .                                                           
     , LALAND,   KTYPE &                               
     , PACLCOVM, pxlvim, pxivim, PAPRLM, PQVIM &             
!!$C - OUTPUT 2D .                                                           
     , PACLC,    PACLCAC &                               
!!$C - OUTPUT 1D .                                                           
     , PACLCOV,  pxlvi, pxivi, PAPRL,  PQVI, PSSFL &    
!!$C - INPUT/OUTPUT 2D .                                                     
     , PTTE,     PQTE, pxlte, pxite &                    
!!$C - INPUT/OUTPUT 1D .                                                     
     , PAPRS,    PRSFL &                                  
     , zmratep, zfprec, zfevap, zmlwc &                                          
     , tropo, krow)                                                 
!!$C                                                                         
!!$C**** *COND* - COMPUTES LARGE-SCALE WATER PHASE CHANGES AND CLOUD COVER.  
!!$C                                                                         
!!$C                                                                         
!!$C     SUBJECT.                                                            
!!$C     --------                                                            
!!$C                                                                         
!!$C          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE THREE     
!!$C     PROGNOSTIC VARIABLES T, Q, AND X (CLOUD WATER) DUE TO WATER PHASE   
!!$C     CHANGES (CONDENSATION, EVAPORATION OF FALLING PRECIPITATION IN      
!!$C     UNSATURATED LAYERS AND MELTING OR FREEZING OF THE FALLING WATER)    
!!$C     AND PRECIPITAION FORMATION (COALESCENSE, SEDIMENTATION).            
!!$C                                                                         
!!$C     RAIN, SNOWFALL, SURFACE FLUXES, CLOUD WATER AND CLOUD COVER         
!!$C     LATER TO BE USED FOR SOIL PROCESSES AND RADIATION ARE STORED.       
!!$C                                                                         
!!$C                                                                         
!!$C**   INTERFACE.                                                          
!!$C     ----------                                                          
!!$C                                                                         
!!$C     *CALL* *COND*                                                       
!!$C                                                                         
!!$C     INPUT ARGUMENTS.                                                    
!!$C     ----- ----------                                                    
!!$C  - 2D                                                                   
!!$C  KLAB     : CONVECTION FLAG (0: NO CONVECTION, 1:  ....                 
!!$C  PACLCACM : CLOUD COVER, ACCUMULATED (OLD VALUE)                        
!!$C  PAPHM1   : PRESSURE AT HALF LEVELS (T-DT)                              
!!$C  PAPHP1   : PRESSURE AT HALF LEVELS (T+DT)                              
!!$C  PAPM1    : PRESSURE AT FULL LEVELS (T-DT)                              
!!$C  PAPP1    : PRESSURE AT FULL LEVELS (T+DT)                              
!!$C  PGEOM1   : GEOPOTENTIAL AT FULL LEVELS (T-DT)                          
!!$C  PQM1     : SPECIFIC HUMIDITY (T-DT)                                    
!!$C  PTM1     : TEMPERATURE (T-DT)                                          
!!$C  PXM1     : CLOUD WATER (T-DT)                                          
!!$C  PXTEC    : TENDENCY OF DETRAINED CONVECTIVE CLOUD WATER                
!!$C  - 1D                                                                   
!!$C  KTYPE    : TYPE OF CONVECTION                                          
!!$C  LALAND   : LOGICAL MASK FOR LAND                                       
!!$C  PACLCOVM : TOTAL CLOUD COVER, ACCUMULATED (OLD VALUE)                  
!!$C  PALWCVIM : VERTICALLY INTEGRATED CLOUD WATER, ACCUMULATED (OLD VALUE)  
!!$C  PAPRLM   : STRATIFORM PRECIPITAION, ACCUMULATED (OLD VALUE)            
!!$C  PQVIM    : VERTICALLY INTEGRATED SPEC. HUMIDITY, ACCUMULATED (OLD VALUE)  
!!$C                                                                         
!!$C     OUTPUT ARGUMENTS.                                                   
!!$C     ------ ----------                                                   
!!$C  - 2D                                                                   
!!$C  PACLC    : CLOUD COVER                                                 
!!$C  PACLCAC  : CLOUD COVER, ACCUMULATED                                    
!!$C  - 1D                                                                   
!!$C  PACLCOV  : TOTAL CLOUD COVER                                           
!!$C  PALWCVI  : VERTICALLY INTEGRATED CLOUD WATER, ACCUMULATED              
!!$C  PAPRL    : STRATIFORM PRECIPITATION, ACCUMULATED                       
!!$C  PQVI     : VERTICALLY INTEGRATED SPEC. HUMIDITY, ACCUMULATED           
!!$C  PSSFL    : SURFACE SNOW FLUX                                           
!!$C                                                                         
!!$C     INPUT/OUTPUT ARGUMENTS.                                             
!!$C     ------------ ----------                                             
!!$C  - 2D                                                                   
!!$C  PTTE     : TENDENCY OF TEMPERATURE                                     
!!$C  PQTE     : TENDENCY OF SPECIFIC HUMIDITY                               
!!$C  PXTE     : TENDENCY OF CLOUD WATER                                     
!!$C  - 1D                                                                   
!!$C  PAPRS    : SNOW FALL, ACCUMULATED                                      
!!$C  PRSFL    : SURFACE RAIN FLUX                                           
!!$C                                                                         
!!$C                                                                         
!!$C     EXTERNALS.                                                          
!!$C     ----------                                                          
!!$C                                                                         
!!$C     METHOD.                                                             
!!$C     -------                                                             
!!$C                                                                         
!!$C          SEE REFERENCES                                                 
!!$C                                                                         
!!$C     REFERENCES.                                                         
!!$C     ----------                                                          
!!$C                                                                         
!!$C     1. LARGE SCALE PHASE CHANGES' PART OF THE MODEL'S DOCUMENTATION     
!!$C     2. ROECKNER, E. AND U. SCHLESE (1985), ECMWF-WORKSHOP ON            
!!$C        "CLOUD COVER PARAMETERIZATION IN NUMERICAL MODELS", PP87-108.    
!!$C     3. ROECKNER ET AL. (1991) ECMWF/WRCP WORKSHOP ON                    
!!$C        "CLOUDS, RADIATIVE TRANSFER AND THE HYDROLOGICAL CYCLE",         
!!$C         199-222, ECMWF, READING, U.K.                                   
!!$C     4. SUNDQVIST, H. (1978), QJRMS, 104, 677-690.                       
!!$C                                                                         
!!$C     AUTHOR.                                                             
!!$C     -------                                                             
!!$C     U.SCHLESE     DKRZ-HAMBURG   NOV-92                                 
!!$C                                                                         
!!$C     MODIFICATIONS.                                                      
!!$C     --------------                                                      
!!$C                                                                         
  USE mo_kind,           ONLY : dp
  USE mo_constants,      ONLY : cpd, vtmpc2, g, rd, alv, als, rv, api    &
       , vtmpc1, rhoh2o, c2es, tmelt, c4les, c5les, c4ies, c5ies
  USE mo_convect_tables, ONLY : lookuperror, lookupoverflow, jptlucu1    &
       , jptlucu2, tlucua, tlucub, tlucuaw
  USE mo_time_control,   ONLY : delta_time, time_step_len
  USE mo_geoloc,         ONLY : philat_2d
  !
  IMPLICIT NONE                              

  !     
  INTEGER :: kbdim, klevp1, KLEVM1, klev2, klevp2, klev2p1, klev, kproma, ktdia, krow
  INTEGER :: jl, jk, jb   

  INTEGER KLAB(kbdim,KLEV), KTYPE(kbdim), ktropo(kbdim)                  
  REAL(dp) PACLCACM(kbdim,KLEV), PAPHM1(kbdim,KLEVP1), PAPHP1(kbdim,KLEVP1) &
       ,PAPM1(kbdim,KLEV),    PAPP1(kbdim,KLEV),    PQM1(kbdim,KLEV) & 
       ,PGEOM1(kbdim,KLEV) &                                           
       ,PTM1(kbdim,KLEV),     PXM1(kbdim,KLEV),     PXTEC(kbdim,KLEV) &
       ,PXLM1(kbdim,KLEV),    PXIM1(kbdim,KLEV) 
  REAL(dp) PACLCOVM(kbdim),      PALWCVIM(kbdim),      PAPRLM(kbdim) &       
       ,PQVIM(kbdim)                                                    
  REAL(dp) PACLC(kbdim,KLEV),    PACLCAC(kbdim,KLEV)                        
  REAL(dp) PACLCOV(kbdim),       PALWCVI(kbdim),       PAPRL(kbdim) &       
       ,PQVI(kbdim),         PSSFL(kbdim)
  REAL(dp) pxlvi(kbdim), pxlvim(kbdim), pxivi(kbdim), pxivim(kbdim)                             
  REAL(dp) PTTE(kbdim,KLEV),     PQTE(kbdim,KLEV),     PXTE(kbdim,KLEV)
  REAL(dp) pxlte(kbdim,klev),    pxite(kbdim,klev)     
  REAL(dp) PAPRS(kbdim),         PRSFL(kbdim),         tropo(kbdim)                               
  !                                                                         
  !    TEMPORARY ARRAYS                                                     
  !                      
  
  INTEGER, PARAMETER :: NEXC=4, NEXL=4
  INTEGER jswitch
  INTEGER INVB(kbdim) 

  REAL(dp), PARAMETER :: ZEPFLM=1.E-24_dp & ! Security Parameters                                      
       , ZEPSEC=1.E-12_dp &                                                      
       , ZEPQP1=0._dp &                                                          
       , ZEPZWP=1.E-20_dp &                                                      
       , ZEPCLC=1.E-2_dp &                                                       
       , ZDNOPR=1.E4_dp

  REAL(dp), PARAMETER ::   ZDELP=1.E4_dp &                                                        
       , ZVTFAC=3.29_dp &
       , ZRS=0.99_dp &                                       
       , ZZEPCONS=0.015_dp &           !Original L39 value: 0.012                                       
       , ZRHSC=0.4_dp &                                                           
       , ZSATSC=0.8_dp &
       , ZRTC=0.6_dp &                 !Original L39 value: 0.7
       , ZRTL=0.6_dp &                 !Original L39 value: 0.7 
       , ZTMSC=0.70_dp &               !Original L39 value: 0.48                                   
       , ZCAA=0.0059_dp &                                                     
       , ZCAB=0.003102_dp &                                                                     
       , ZWRSC=0.5E-3_dp &                                                       
       , ZWRL=0.5E-3_dp &                                                        
       , ZWRS=0.3E-3_dp &                            
       , ZCLOIA=1.0E+02_dp 

  REAL(dp) :: ztmst, zqtmst, zdiagt, zdiagw, zcons, zconsi, zcons1, zcons2, zmeltp2, z2eS
  REAL(dp) :: zsnmlt, zdt, zqr, zcor, it, zzqt, zcvm5, zcvm4, zlat, zrcp, zztc, zs
  REAL(dp) :: zwmvn, zwi, zwpvn, zxetmst, zcpn, zwma, zwpa, zqpa, zqsm1, zcnull, zc1, rdcpd
  REAL(dp) :: zclcm1, zwmax, zwps, ztc, zl, zrho, zzfac, zalpha, zwpni, zaci, zqq, zwpvnl
  REAL(dp) :: zdacl, zdaci, zwrc, zwpvni, zacl, zpcc, zzdr, zclear, zzrh, zzep, zfrac
  REAL(dp) :: zrflo, zsflo, zdqdt, zcvml, zdtdt, zdxdt, zclp1, zalwcold, zdxcor, zdthdp

  REAL(dp):: ZALWC(kbdim,klev)    ,ZBAS(kbdim,klev)     ,ZCLCP1(kbdim,klev) &    
       ,ZSAT(kbdim,klev) &                                          
       ,ZCONRAT(kbdim)       ,ZDAC(kbdim) &                               
       ,ZDP(kbdim)           ,ZEP(kbdim),           ZEPR(kbdim) &          
       ,ZEPS(kbdim)          ,ZFLN(kbdim),          ZGEOH(kbdim,klevp1) &   
       ,ZLSDCP(kbdim,klev),   ZLVDCP(kbdim,klev) &                     
       ,ZQCON(kbdim,klev)    ,ZQP1(kbdim,klev) &                          
       ,ZQSP1(kbdim,klev)    ,ZRAIN(kbdim),         ZRFL(kbdim) &           
       ,ZRFLN(kbdim)         ,ZRHC(kbdim,klev),     ZSFL(kbdim) &          
       ,ZSFLN(kbdim)         ,ZSNOW(kbdim),         ZSNMT(kbdim) &         
       ,ZTP1(kbdim,klev)     ,ZTOP(kbdim,klev) &                         
       ,ZDTHMIN(kbdim)       ,ZTHETA(kbdim,klev)                       
  !                                                                         
  LOGICAL LO1,LO,LOLL                                                 
  LOGICAL LALAND(kbdim)  

  ! eth_as_scav
  REAL(dp):: zmratep(kbdim,klev), & ! precipitation formation rate [kg/m3s]
       zfprec(kbdim,klev),  & ! flux of precip [kg/m2s]
       zfevap(kbdim,klev),  & ! evaporation of precip [kg/m2s]  
       zclcover(kbdim,klev), &                          
       zmlwc(kbdim,klev),RQSP1(kbdim,klev)                                 

  ! Executive Statements

!!$C                                                                         
!!$C*    COMPUTATIONAL CONSTANTS.                                            
!!$C     ------------- ----------                                            
!!$C                                                                                    
  ZTMST=time_step_len                                
  ZQTMST=1._dp/ZTMST                                                     
  ZDIAGT=0.5_dp*time_step_len                                                    
  ZDIAGW=ZDIAGT/RHOH2O                                                
  !                                                                         
  ZCONS1=CPD*VTMPC2                                                   
  ZCONS2=1._dp/(ZTMST*G)                                                 
  ZMELTP2=TMELT+2._dp                                                    
  Z2ES=C2ES                                                           
  ZCNULL=1.5E-4_dp*ZTMST      !Original L39 value: 2.0E-4*ZTMST             
  ZCONSI=ZVTFAC*G*ZTMST*ZTMSC/ZDELP
  ZC1=2._dp*ZTMST 
  KLEVM1=KLEV-1                                                       
  KLEV2=KLEV/2                                                        
  KLEV2P1=KLEV2+1                                                     
  RDCPD=RD/CPD

  palwcvi(1:kproma) = pxlvi(1:kproma) + pxivi(1:kproma)
  palwcvim(1:kproma) = pxlvim(1:kproma) + pxivim(1:kproma)
  pxm1(1:kproma,:) = pxlm1(1:kproma,:) + pxim1(1:kproma,:)
  pxte(1:kproma,:) = pxlte(1:kproma,:) + pxite(1:kproma,:)

  ! Initialize arrays
  zmratep = 0._dp
  zfprec = 0._dp
  zfevap = 0._dp
  zclcover = 0._dp
  zmlwc = 0._dp                                                  
!!$C                                                                         
!!$C     ------------------------------------------------------------------  
!!$C                                                                         
!!$C*         2.  TOP BOUNDARY CONDITIONS AND QUANTITIES NEEDED FOR          
!!$C*             --- -------- ---------- --- ---------- ------ ---          
!!$C*                CONDENSATION AND PRECIPITATION CALCULATIONS.            
!!$C*                ------------ --- ------------- ------------             
!!$C                                                                                                     
!!$C                                                                         
!!$C*         2.1     SET TO ZERO PRECIPITATION FLUXES AT THE TOP.           
!!$C                                                                         
  DO JL=1, kproma                                         
     ZRFL(JL)=0._dp                                                     
     ZSFL(JL)=0._dp                                                         
     ZRAIN(JL)=0._dp                                                       
     ZSNOW(JL)=0._dp                                                       
     ZEPR(JL)=0._dp                                                        
     ZEPS(JL)=0._dp                                                         
     ZSNMT(JL)=0._dp                                                      
     ZRFLN(JL)=0._dp                                                        
     ZSFLN(JL)=0._dp                                                        
     ZDTHMIN(JL)=0._dp                                                     
     INVB(JL)=1                                                          
  END DO
  ZS=0._dp                                                               
!!$C                                                                         
!!$C*        2.2      CALCULATE POTENTIAL TEMPERATURES                       
!!$C                                                                                                    
  DO  JK=KLEV,KLEV2,-1                                             
     DO  jl = 1, kproma                                               
        ZTHETA(JL,JK)=PTM1(JL,JK)*(1.E5_dp/PAPM1(JL,JK))**RDCPD               
     END DO
  END DO
!!$C                                                                         
!!$C*         2.3    CHECK FOR PRESENCE OF LOW-LEVEL INVERSION               
!!$C                  ( SEA POINTS ONLY )                                    
!!$C                                                                                                 
  DO JK=KLEV,KLEV2P1,-1                                           
     DO  jl = 1, kproma                                               
        IF((.NOT.LALAND(JL)) .AND. KTYPE(JL).EQ.0) THEN                     
           ZDTHDP=(ZTHETA(JL,JK)-ZTHETA(JL,JK-1))*ZCLOIA/(PAPM1(JL,JK) &       
                -PAPM1(JL,JK-1))                                           
           LO=ZDTHDP.LT.ZDTHMIN(JL)
           IF (LO) THEN
              ZDTHMIN(JL) = ZDTHDP
              INVB(JL) = JK
           ELSE
              ZDTHMIN(JL) = ZDTHMIN(JL)
              INVB(JL) = INVB(JL)
           END IF                                    
        END IF
     END DO
  END DO
!!$C                                                                         
!!$C*      2.3.1   DETERMINE LAYER ABOVE THE TROPOPAUSE             
!!$C                                                                                        
  ZLAT=philat_2d(1,krow)                                             
!!$C                                                                        
!!$C HOT TOWERS MUST NOT BE FORBIDDEN!                                      
!!$C COND IS ALLOWED TO WORK IN THE troposphere as usual,                   
!!$c in the stratosphere only at latitudes between 55S and 55N                                               
!!$C                                                                        
  JSWITCH=0                                                          
  IF (ABS(ZLAT).GT.55._dp*API/180._dp) THEN                                
     JSWITCH=1                                                        
     DO JK=1,KLEV                                                 
        DO jl = 1, kproma                                          
           IF ((PAPHP1(JL,JK)  .LE.TROPO(JL)).AND. &            
                (PAPHP1(JL,JK+1).GE.TROPO(JL))) ktropo(jl)=jk       
        END DO
     END DO
  ENDIF

  DO JK=KTDIA,KLEV                                                
     DO jl = 1, kproma                                               
!!$C                                                                         
!!$C*     2.4   T, Q AND QS PROVISIONAL VALUES AT T+DT                       
!!$C*            EFFECTIVES L FOR VAPORISATION AND SUBLIMATION               
!!$C                                                                         
        ZTP1(JL,JK)=PTM1(JL,JK)+ZTMST*PTTE(JL,JK)                           
        ZQP1(JL,JK)=PQM1(JL,JK)+ZTMST*PQTE(JL,JK)                           
        ZQP1(JL,JK)=MAX(ZQP1(JL,JK),ZEPQP1)                               

        ZRCP=1._dp/(CPD+ZCONS1*ZQP1(JL,JK))                                    
        ZLVDCP(JL,JK)=ALV*ZRCP                                              
        ZLSDCP(JL,JK)=ALS*ZRCP                                              

        ZZTC=TMELT-ZTP1(JL,JK)  
        IF (zztc .lt. 0._dp) THEN
           ZCVM4=C4LES
           ZCVM5=C5LES*ZLVDCP(JL,JK)
        ELSE
           ZCVM4=C4IES
           ZCVM5=C5IES*ZLSDCP(JL,JK)
        END IF          
        ZZQT=1._dp/(ZTP1(JL,JK)-ZCVM4)                                         
        IT=ZTP1(JL,JK)*1000._dp                                                
        ZQSP1(JL,JK)=TLUCUA(IT)/PAPP1(JL,JK)                                
        ZQSP1(JL,JK)=MIN(ZQSP1(JL,JK),0.5)                                    
        ZCOR=1._dp/(1._dp-VTMPC1*ZQSP1(JL,JK))                                    
        ZQSP1(JL,JK)=ZQSP1(JL,JK)*ZCOR                                      
        ZQCON(JL,JK)=1._dp/(1._dp+ZCVM5*ZQSP1(JL,JK)*ZCOR*ZZQT**2)                
!!$C                                                                         
!!$C*     2.5    THRESHOLD RELATIVE HUMIDITY                                 
!!$C                                                                         
        JB=INVB(JL)                                                         
        LO=(JB.GT.KLEV-8.AND. JB.LT.KLEVM1)                                    
        IF(KLAB(JL,JK).EQ.2) THEN                                           
           ZRHC(JL,JK)=ZRTC+(ZRS-ZRTC)*EXP(1._dp- &                              
                (PAPHP1(JL,KLEVP1)/PAPP1(JL,JK))**NEXC)                         
           ZSAT(JL,JK)=1._dp                                                    
        ELSE                                                                
           IF ((JK.GE.JB+1).AND.LO) THEN                                      
              ZRHC(JL,JK)=ZRHSC                                               
              ZSAT(JL,JK)=ZSATSC                                              
           ELSE                                                              
              ZRHC(JL,JK)=ZRTL+(ZRS-ZRTL)*EXP(1._dp- &                      
                   (PAPHP1(JL,KLEVP1)/PAPP1(JL,JK))**NEXL)                       
              ZSAT(JL,JK)=1._dp                                                  
           END IF
        ENDIF
!!$C                                                                         
!!$C*    2.6  CLOUD COVER AT T+DT                                            
!!$C                                                                         
        ZQR=ZQSP1(JL,JK)*ZSAT(JL,JK)*ZRHC(JL,JK)                            
        ZCLCP1(JL,JK)=(ZQP1(JL,JK)-ZQR)/(ZQSP1(JL,JK)*ZSAT(JL,JK)-ZQR)      
        ZCLCP1(JL,JK)=MAX(ZCLCP1(JL,JK),0._dp)                               
        ZCLCP1(JL,JK)=MIN(ZCLCP1(JL,JK),1._dp)                               
        ZCLCP1(JL,JK)=1._dp-SQRT(1._dp-ZCLCP1(JL,JK))                             
        !                                                                         
        IF((JK.LT.KTROPO(JL)).AND.(JSWITCH.EQ.1))THEN     
           ZCLCP1(JL,JK)=0._dp                                                
        ENDIF
        ZCLCOVER(JL,JK)=ZCLCP1(JL,JK)                                        
     END DO
  END DO
!!$C                                                                         
!!$C                                                                         
!!$C*      2.8 CLOUD TOPS/BASES AND THICKNESS.                               
!!$C                                                                                                     
!!$C                                                                         
!!$C      2.81 GEOPOTENTIAL AT HALF LEVELS                                   
!!$C                                                                         
  DO  JK=KTDIA,KLEV-1                                              
     DO jl = 1, kproma                                               
        ZGEOH(JL,JK+1)=(PGEOM1(JL,JK+1)+PGEOM1(JL,JK))*0.5                  
     END DO
  END DO

  DO jl = 1, kproma                                               
     ZGEOH(JL,1)=2._dp*PGEOM1(JL,1)-ZGEOH(JL,2)                             
     ZGEOH(JL,KLEVP1)=0._dp                                                 
  END DO
!!$C                                                                         
!!$C       2.82  CLOUD TOPS                                                  
!!$C                                                                         
  DO jl = 1, kproma                                               
     IF(ZCLCP1(JL,1).GE.ZEPCLC) THEN                                     
        ZTOP(JL,1)=ZGEOH(JL,1)                                            
     ELSE                                                                
        ZTOP(JL,1)=-1._dp                                                    
     ENDIF
  END DO
  !                                                                         
  DO JK=KTDIA+1,KLEV                                              
     DO  jl = 1, kproma                                               
        IF(ZCLCP1(JL,JK).GE.ZEPCLC.AND.ZTOP(JL,JK-1).LT.0._dp) THEN            
           ZTOP(JL,JK)=ZGEOH(JL,JK)                                          
        ELSEIF(ZCLCP1(JL,JK).GE.ZEPCLC) THEN                                
           ZTOP(JL,JK)=ZTOP(JL,JK-1)                                         
        ELSE                                                                
           ZTOP(JL,JK)=-1._dp                                                   
        ENDIF
     END DO
  END DO
!!$C                                                                         
!!$C      2.83 CLOUD BASES                                                   
!!$C                                                                         
  DO jl = 1, kproma                                               
     IF(ZCLCP1(JL,KLEV).GE.ZEPCLC) THEN                                  
        ZBAS(JL,KLEV)=ZGEOH(JL,KLEVP1)                                    
     ELSE                                                                
        ZBAS(JL,KLEV)=-1._dp                                                 
     ENDIF
  END DO
  !                                                                        
  DO JK=KLEV-1,KTDIA,-1                                           
     DO jl = 1, kproma                                               
        IF(ZCLCP1(JL,JK).GE.ZEPCLC.AND.ZBAS(JL,JK+1).LT.0._dp) THEN            
           ZBAS(JL,JK)=ZGEOH(JL,JK+1)                                        
        ELSEIF(ZCLCP1(JL,JK).GE.ZEPCLC) THEN                                
           ZBAS(JL,JK)=ZBAS(JL,JK+1)                                         
        ELSE                                                                
           ZBAS(JL,JK)=-1._dp                                                   
        ENDIF
     END DO
  END DO
!!$C                                                                         
!!$C***                                                                      
  DO JK=KTDIA,KLEV                                                
!!$C***                                                                      
!!$C                                                                         
!!$C     ------------------------------------------------------------------  
!!$C                                                                         
!!$C*         3.   SNOW MELT AND SATURATION SPECIFIC HUMIDITIES AT T-DT.     
!!$C               ---- ---- --- ---------- -------- ---------- -- ----      
!!$C                                                                              
     !                                                                                             
     DO jl = 1, kproma                                               
!!$C                                                                         
!!$C*         3.1   MELTING OF INCOMING SNOW                                 
!!$C                                                                         
!!$C                                                                         
        ZDP(JL)=PAPHP1(JL,JK+1)-PAPHP1(JL,JK)                               
        ZCONS=ZCONS2*ZDP(JL)/(ZLSDCP(JL,JK)-ZLVDCP(JL,JK))                  
        ZSNMLT=MIN(ZSFL(JL),ZCONS*MAX(0._dp,(ZTP1(JL,JK)-ZMELTP2)))        
        ZSNMT(JL)=ZSNMT(JL)+ZSNMLT                                          
        ZRFLN(JL)=ZRFL(JL)+ZSNMLT                                           
        ZSFLN(JL)=ZSFL(JL)-ZSNMLT                                           
        ZDT=-ZSNMLT/ZCONS                                                   
        ZTP1(JL,JK)=ZTP1(JL,JK)+ZDT                                         
        PTTE(JL,JK)=PTTE(JL,JK)+ZDT*ZQTMST                                  
!!$C                                                                         
!!$C*    3.3   SATURATION SPECIFIC HUMIDITY FOR T-DT.                        
!!$C                                                                         
        IT=PTM1(JL,JK)*1000._dp                                                
        ZQSM1=TLUCUA(IT)/PAPM1(JL,JK)                                       
        ZQSM1=MIN(ZQSM1,0.5)                                                  
        ZCOR=1._dp/(1._dp-VTMPC1*ZQSM1)                                           
        ZQSM1=ZQSM1*ZCOR                                                    
!!$C                                                                         
!!$C     ----------------------------------------------------------------    
!!$C                                                                         
!!$C*        4.   PROVISIONAL VALUES OF CLOUD COVER, HUMIDITY                
!!$C              ------------ ------ -- ----- -----  --------               
!!$C                        AND CLOUD WATER.                                 
!!$C                        --- ----- -----                                  
!!$C                                                                         
!!$C                                                                         
!!$C*        4.1  CLOUD COVER AT  T-DT                                       
!!$C                                                                         
!!$C                                                                         
        ZQR=ZQSM1*ZSAT(JL,JK)*ZRHC(JL,JK)                                   
        ZCLCM1=(PQM1(JL,JK)-ZQR)/(ZQSM1*ZSAT(JL,JK)-ZQR)                    
        ZCLCM1=MAX(ZCLCM1,0._dp)                                             
        ZCLCM1=MIN(ZCLCM1,1._dp)                                             
        ZCLCM1=1._dp-SQRT(1._dp-ZCLCM1)                                           
        !                                                                         
        IF ((JK.LT.KTROPO(JL)).AND.(JSWITCH.EQ.1)) THEN     
           ZCLCM1=0._dp                                                          
        ENDIF
!!$C*        4.3  HUMIDITY AND LIQUID WATER IN THE CLOUDY AND CLOUD-FREE     
!!$C              PART FOR T-1.                                              
!!$C                                                                         
        IF(ZCLCM1.GT.1.E-10) THEN                                           
           ZQPA=MAX(ZQSM1,PQM1(JL,JK))                                      
        ELSE                                                                
           ZQPA=PQM1(JL,JK)                                                   
        ENDIF
        IF(ZCLCP1(JL,JK).GT.1.E-10.AND.PXM1(JL,JK).GT.0._dp) THEN              
           ZWPA=PXM1(JL,JK)/ZCLCP1(JL,JK)                                     
           ZWMA=0._dp                                                            
        ELSE                                                                
           ZWPA=PXM1(JL,JK)                                                   
           ZWMA=ZWPA                                                          
        ENDIF
!!$C                                                                         
!!$C     -----------------------------------------------------------------   
!!$C                                                                         
!!$C*         5.   CONDENSATION/EVAPORATION RATE.                            
!!$C               ------------ ----------- ----                             
!!$C                                                                         
        ZCPN=(ZQPA+PQTE(JL,JK)*ZTMST-ZQSP1(JL,JK))*ZQCON(JL,JK)             
        ZXETMST=PXTE(JL,JK)*ZTMST                                           
        ZWPVN=ZWPA+ZXETMST                                                  
        ZWMVN=ZWMA+ZXETMST                                                  
        LO=ZWPVN.GE.0._dp                                                      
        LO1=LO.AND.-ZCPN.GT.ZWPVN 
        IF (LO1) THEN
           ZCPN = -ZWPVN
        ELSE
           ZCPN = ZCPN
        END IF                                      
        ZWI=-ZWPVN                                                          
        ZWMAX=MAX(ZCPN,0._dp)                                                
        ZWPS=ZWPVN+ZCPN
        IF (LO) THEN
           ZWPVN = ZWPS
           ZCPN = ZCPN
        ELSE
           ZWPVN = ZWMAX
           ZCPN = ZWPVN+ZWI
        END IF                                                                               
        ZCONRAT(JL)=ZCLCP1(JL,JK)*ZCPN+(1._dp-ZCLCP1(JL,JK))*(-ZWMVN)          
!!$C                                                                         
!!$C     -----------------------------------------------------------------   
!!$C                                                                         
!!$C*       6.  CLOUD PHYSICS AND PRECIPITATION FLUXES.                      
!!$C            ----- ------- --- ------------- ------                       
!!$C                                                                         
!!$C                                                                         
!!$C*       6.1  AUTOCONVERSION AND ACCRETION OF CLOUD DROPLETS IN WARM      
!!$C*            CLOUDS AND SEDIMENTATION OF ICE CRYSTALS IN COLD CLOUDS     
!!$C                                                                         
        ZTC=ZTP1(JL,JK)-TMELT
        IF (-ZTC .LT. 0._dp) THEN
           ZL=ZLVDCP(JL,JK)
        ELSE
           ZL=ZLSDCP(JL,JK)
        END IF
        IF(ZCONRAT(JL).GT.0._dp) THEN                                          
           ZCONRAT(JL)=MIN(ZCONRAT(JL),ZQP1(JL,JK))                        
        ENDIF
        ZTP1(JL,JK)=ZTP1(JL,JK)+ZL*ZCONRAT(JL)                              
        ZQP1(JL,JK)=ZQP1(JL,JK)-   ZCONRAT(JL)                              
        ZWPVN=ZWPVN+PXTEC(JL,JK)*ZTMST                                      
        IF(ZWPVN*ZCLCP1(JL,JK).GT.ZEPZWP) THEN  
           IF (LALAND(JL)) THEN
              ZWRC=ZWRL
           ELSE
              ZWRC=ZWRS
           END IF
           LO=ZTOP(JL,JK)-ZBAS(JL,JK).LT.ZDNOPR 
           IF (LO) THEN
              ZWRC=ZWRSC
           ELSE
              ZWRC=ZWRC
           END IF
           ZTC=ZTP1(JL,JK)-TMELT                                             
           ZRHO=PAPP1(JL,JK)/(RD*ZTP1(JL,JK))  
           IF (ZTC .LT. 0._dp) THEN
              ZZFAC=1._dp
           ELSE
              ZZFAC=0._dp
           END IF
           ZALPHA=(1._dp-ZZFAC)+ZZFAC*(ZCAA+(1._dp-ZCAA)*EXP(-ZCAB*ZTC**2))      
           ZWPVNI=(1._dp-ZALPHA)*ZWPVN                                          
           ZACI=ZCONSI*(ZWPVNI*ZRHO)**1.16                                   
           ZWPVNL=ZALPHA*ZWPVN                                               
           ZQQ=(ZWPVNL/ZWRC)**2                                              
           ZQQ=MAX(ZQQ,1.E-6)                                              
           ZQQ=MIN(ZQQ,100._dp)                                               
           ZACL=ZCNULL*ZWPVNL*(1._dp-EXP(-ZQQ))                               
           ZPCC=ZWPVNL*(ZRFLN(JL)+ZSFLN(JL))*ZC1                             
           ZPCC=MAX(ZPCC,0._dp)                                               
           ZDACL=MIN((ZACL+ZPCC),ZWPVNL)                                   
           ZDACI=MIN(ZWPVNI,MAX(0.1*ZWPVNI,ZACI))                        
           ZDAC(JL)=(ZDACL+ZDACI)*ZCLCP1(JL,JK)                              
!!$C - - -                                                                 
!!$C     ZMRATEP : PRECIPITATION FORMATION RATE [KG/(M3*SEC)]              
!!$C     ZMLWC : LIQUID WATER CONTENT [KGH2O/KGAIR]                        
!!$C - - -                                                                 
           ZMRATEP(JL,JK)=ZDACL*ZRHO/ZTMST                                   
           ZMLWC(JL,JK)=ZWPVNL                                               
        ELSE                                                                
           ZDAC(JL)=ZWPVN*ZCLCP1(JL,JK)                                      
        ENDIF
        !                                                                         
        ZEP(JL)=0._dp                                                          
!!$C                                                                         
!!$C*      6.3   NEW PRECIPITATION FLUXES AFTER CONSIDERATION OF             
!!$C*            COALESCENCE AND SEDIMENTATION PROCESSES.                    
!!$C                                                                         
        ZZDR=MAX(0._dp,ZCONS2*ZDP(JL)*ZDAC(JL))                              
        IF (ZTP1(JL,JK).GT.TMELT) THEN                                      
           ZRFLN(JL)=ZRFLN(JL)+ZZDR                                           
           ZRAIN(JL)=ZRAIN(JL)+ZZDR                                           
        ELSE                                                                
           ZSFLN(JL)=ZSFLN(JL)+ZZDR                                           
           ZSNOW(JL)=ZSNOW(JL)+ZZDR                                           
        ENDIF
        ZFLN(JL)=ZRFLN(JL)+ZSFLN(JL)     
     END DO
!!$C                                                                         
!!$C*         6.5     AT TOP LAYER OR IF THERE IS NO PRECIPITATION AVOID     
!!$C*                 EVAPORATION.                                           
!!$C                                                                                                       
     !                                                                      
     IF (ZS.EQ.0._dp.AND.JK.GT.1) THEN                                      
        ZS=SUM(ZFLN)                                            
     END IF
     IF (JK.GT.1.AND.ZS.NE.0._dp) THEN                                      
!!$C***                                                                      
!!$C                                                                         
!!$C     ------------------------------------------------------------------  
!!$C                                                                         
!!$C*         7.     EVAPORATION OF PRECIPITATION.                           
!!$C*                ----------- -- -------------                            
!!$C                                                                                                        
        DO jl = 1, kproma                                               
           ZCLEAR=1._dp-ZCLCP1(JL,JK)                                             
           ZZRH=ZRHC(JL,JK)*ZSAT(JL,JK) &                                   
                +ZCLCP1(JL,JK)*(1._dp-ZRHC(JL,JK)*ZSAT(JL,JK))                     
           ZZEP=MIN(0._dp,(MIN(ZQP1(JL,JK), &                                  
                ZQSP1(JL,JK)*ZZRH)-ZQSP1(JL,JK)))                        
           ZZEP=ZZEPCONS*ZZEP                                                  
           ZZEP=MAX(-(ZRFL(JL)+ZSFL(JL)),ZZEP*ZQCON(JL,JK)*ZCONS2 &       
                *ZDP(JL))       
!!$C - - -                                                                  
!!$C     ZFPREC : FLUX OF PRECIPITAION [KG/(M2*SEC)]                        
!!$C     ZFEVAP : TOTALER FLUSS, DER IN EINEM ZEITSCHRITT                   
!!$C              VERDUNSTET [KG/(M2*SEC)]                                  
!!$C - - -                                                                  
           ZFPREC(JL,JK)=ZRFLN(JL)                                            
           ZFEVAP(JL,JK)=-ZZEP*ZCLEAR                                         
           ZFEVAP(JL,JK)=ZFEVAP(JL,JK)*ZRFLN(JL)/ &                             
                MAX(ZRFLN(JL)+ZSFLN(JL),ZEPFLM)                      
           ZFLN(JL)=ZFLN(JL)+ZZEP*ZCLEAR                                       
           ZEP(JL)=ZZEP/(ZCONS2*ZDP(JL))*ZCLEAR                                
           ZFRAC=ZFLN(JL)/MAX(ZRFLN(JL)+ZSFLN(JL),ZEPFLM)                    
           ZRFLO=ZRFLN(JL)                                                     
           ZSFLO=ZSFLN(JL)                                                     
           ZRFLN(JL)=ZRFLN(JL)*ZFRAC                                           
           ZSFLN(JL)=ZSFLN(JL)*ZFRAC                                           
           ZEPR(JL)=ZEPR(JL)+(ZRFLO-ZRFLN(JL))                                 
           ZEPS(JL)=ZEPS(JL)+(ZSFLO-ZSFLN(JL))                                 
        END DO
!!$C                                                                         
!!$C     ------------------------------------------------------------------  
!!$C                                                                         
!!$C*      8.    INCREMENTATION OF T,Q AND X TENDENCIES AND FLUXES' SWAP.    
!!$C             -------------- -- - - --- - ---------- --- ------- ----     
!!$C                                                            
!!$C                                                                         
!!$C*         8.1     RESUME COMPUTATIONS FOR THE TOP LAYER.                 
!!$C                                                            
!!$C***                                                                      
     END IF
!!$C***                                                                      
!!$C                                                                         
!!$C*         8.2     MODIFICATION OF THE T,Q AND X TENDENCIES,              
!!$C*                 CLOUD COVER FOR DIAGNOSTICS AND RADIATION.             
!!$C                                                                                                    
     DO jl = 1, kproma                                               
        ZDQDT=-(ZCONRAT(JL)+ZEP(JL))*ZQTMST                                 
        ZZTC=ZTP1(JL,JK)-TMELT
        IF (-ZZTC .LT. 0._dp) THEN
           ZCVML=ZLVDCP(JL,JK)  
        ELSE
           ZCVML=ZLSDCP(JL,JK)
        END IF
        ZDTDT=-ZCVML*ZDQDT                                                 
        ZDXDT=(ZCONRAT(JL)-ZDAC(JL))*ZQTMST                                
        PXTE(JL,JK)=PXTE(JL,JK)+PXTEC(JL,JK)+ZDXDT                         
        !                                                                        
        IT=ZTP1(JL,JK)*1000._dp                                              
        ZQSP1(JL,JK)=TLUCUA(IT)/PAPP1(JL,JK)                               
        ZQSP1(JL,JK)=MIN(ZQSP1(JL,JK),0.5)                                   
        ZCOR=1._dp/(1._dp-VTMPC1*ZQSP1(JL,JK))                                   
        ZQSP1(JL,JK)=ZQSP1(JL,JK)*ZCOR                                     
        !                                                                        
        ZQR=ZQSP1(JL,JK)*ZSAT(JL,JK)*ZRHC(JL,JK)                           
        ZCLP1=(ZQP1(JL,JK)-ZQR)/(ZQSP1(JL,JK)*ZSAT(JL,JK)-ZQR)             
        ZCLP1=MAX(ZCLP1,0._dp)                                              
        ZCLP1=MIN(ZCLP1,1._dp)                                             
        ZCLP1=1._dp-SQRT(1._dp-ZCLP1)                                            
        !                                                                         
        IF((JK.LT.KTROPO(JL)).AND.(JSWITCH.EQ.1))THEN       
           ZCLP1=0._dp                                                             
        ENDIF
        ZALWC(JL,JK)=PXM1(JL,JK)+ZTMST*PXTE(JL,JK)                         
        ZALWCOLD=ZALWC(JL,JK)                                              
        LO=(ZALWC(JL,JK).LT.ZEPZWP.OR.ZCLP1.LT.ZEPCLC) 
        IF (LO) THEN
           ZALWC(JL,JK)=0._dp
        ELSE
           ZALWC(JL,JK)=ZALWC(JL,JK)
        END IF
        ZDXCOR=(ZALWC(JL,JK)-ZALWCOLD)/ZTMST                               
        PXTE(JL,JK)=PXTE(JL,JK)+ZDXCOR                                    
!!$        PQTE(JL,JK)=PQTE(JL,JK)-ZDXCOR+ZDQDT                              
!!$        PTTE(JL,JK)=PTTE(JL,JK)+ZCVML*ZDXCOR+ZDTDT
        IF (LO) THEN
           PACLC(JL,JK)=0._dp
        ELSE
           PACLC(JL,JK)=ZCLP1
        END IF
        PACLCAC(JL,JK)=PACLCACM(JL,JK)+PACLC(JL,JK)*ZDIAGT                   
     END DO
!!$C                                                                         
!!$C*         8.4     SWAP OF THE FLUXES AND END OF THE VERTICAL LOOP.       
!!$C                                                                                                      
     DO jl = 1, kproma                                               
        ZRFL(JL)=ZRFLN(JL)                                                  
        ZSFL(JL)=ZSFLN(JL)                                                  
     END DO
!!$      C***                                                                      
  END DO
!!$C***                                                                      
!!$C                                                                         
!!$C*         8.5     SURFACE FLUXES.                                        
!!$C                                                                                        
!!$  DO jl = 1, kproma                                               
!!$     PRSFL(JL)=PRSFL(JL)+ZRFL(JL)                                        
!!$     PSSFL(JL)=ZSFL(JL)                                                  
!!$     PAPRL(JL)=PAPRLM(JL)+ZDIAGW*(PRSFL(JL)+PSSFL(JL))                   
!!$     PAPRS(JL)=PAPRS(JL)+ZDIAGW*PSSFL(JL)                                
!!$  END DO
!!$C                                                                         
!!$C*          8.6    ACCUMULATED TOTAL CLOUDCOVER                           
!!$C                                                                                                      
!!$  DO jl = 1, kproma                                               
!!$     PACLCOV(JL)=1._dp-PACLC(JL,1)                                          
!!$  END DO
!!$  DO JK=2,KLEV                                                    
!!$     DO jl = 1, kproma                                               
!!$        PACLCOV(JL)=PACLCOV(JL)*(1._dp-MAX(PACLC(JL,JK),PACLC(JL,JK-1))) &    
!!$             /(1._dp-MIN(PACLC(JL,JK-1),1._dp-ZEPSEC))          
!!$     END DO
!!$  END DO
!!$  DO jl = 1, kproma                                               
!!$     PACLCOV(JL)=1._dp-PACLCOV(JL)                                          
!!$     PACLCOV(JL)=PACLCOVM(JL)+ZDIAGT*PACLCOV(JL)                         
!!$  END DO
!!$C                                                                         
!!$C*          8.7     VERTICALLY INTEGRATED HUMIDITY AND CLOUD WATER        
!!$C                                                                                                     
!!$  DO jl = 1, kproma                                               
!!$     PQVI(JL)=0._dp                                                         
!!$     PALWCVI(JL)=0._dp                                                      
!!$  END DO
!!$  !                                                                         
!!$  DO JK=KTDIA,KLEV                                                
!!$     DO jl = 1, kproma                                               
!!$        PQVI(JL)=PQVI(JL)+PQM1(JL,JK)*(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))       
!!$        PALWCVI(JL)=PALWCVI(JL)+ZALWC(JL,JK) &                              
!!$             *(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))                        
!!$     END DO
!!$  END DO
!!$  !                                                                         
!!$  DO jl = 1, kproma                                               
!!$     PQVI(JL)=PQVIM(JL)+ZDIAGT*PQVI(JL)/G                                
!!$     PALWCVI(JL)=PALWCVIM(JL)+ZDIAGT*PALWCVI(JL)/G                       
!!$  END DO                      
  !                                                                        
  RETURN                                                             
END SUBROUTINE cond_scav

!-------------------------------------------------------------------------------
  SUBROUTINE construct_stream_scav

    ! Allocates output streams

    ! *construct_stream_ch4* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (scav ,'scav', filetype=NETCDF)      
  
  ! add entries for geopotential and log surface pressure (used by the
  ! afterburner):
  CALL add_stream_reference (scav, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
  CALL add_stream_reference (scav, 'lsp'     ,'sp'    ,lpost=.TRUE.)

  CALL default_stream_setting (scav, lrerun    = .TRUE. , &
                                       leveltype = SURFACE, &
                                       table     = 199,     &
                                       laccu     = .TRUE.)

  CALL add_stream_element (scav, 'WD_CH3O2H', wd_ch3o2h, lpost=.TRUE., &
                           longname='Large-scale wet deposition CH3O2H', &
                           units='kg/m2/s', code=100)

  CALL add_stream_element (scav, 'WD_CH2O', wd_ch2o, lpost=.TRUE., &
                           longname='Large-scale wet deposition CH2O', &
                           units='kg/m2/s', code=100)

  CALL add_stream_element (scav, 'WD_H2O2', wd_h2o2, lpost=.TRUE., &
                           longname='Large-scale wet deposition H2O2', &
                           units='kg/m2/s', code=100)

  CALL add_stream_element (scav, 'WD_HNO3', wd_hno3, lpost=.TRUE., &
                           longname='Large-scale wet deposition HNO3', &
                           units='kg/m2/s', code=100)

  CALL add_stream_element (scav, 'WD_O3', wd_o3, lpost=.TRUE., &
                           longname='Large-scale wet deposition O3', &
                           units='kg/m2/s', code=100)


  CALL add_stream_element (scav, 'WDCV_CH3O2H', wdcv_ch3o2h, lpost=.TRUE., &
                           longname='Convective wet deposition CH3O2H', &
                           units='kg/m2/s', code=101)

  CALL add_stream_element (scav, 'WDCV_CH2O', wdcv_ch2o, lpost=.TRUE., &
                           longname='Convective wet deposition CH2O', &
                           units='kg/m2/s', code=101)

  CALL add_stream_element (scav, 'WDCV_H2O2', wdcv_h2o2, lpost=.TRUE., &
                           longname='Convective wet deposition H2O2', &
                           units='kg/m2/s', code=101)

  CALL add_stream_element (scav, 'WDCV_HNO3', wdcv_hno3, lpost=.TRUE., &
                           longname='Convective wet deposition HNO3', &
                           units='kg/m2/s', code=101)

  CALL add_stream_element (scav, 'WDCV_O3', wdcv_o3, lpost=.TRUE., &
                           longname='Convective wet deposition O3', &
                           units='kg/m2/s', code=101)


  CALL add_stream_element (scav, 'WD_CH3O2H_d', wd_ch3o2h_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Large-scale wet deposition CH3O2H', &
                           units='kg/m2/s', code=106)

  CALL add_stream_element (scav, 'WD_CH2O_d', wd_ch2o_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Large-scale wet deposition CH2O', &
                           units='kg/m2/s', code=106)

  CALL add_stream_element (scav, 'WD_H2O2_d', wd_h2o2_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Large-scale wet deposition H2O2', &
                           units='kg/m2/s', code=106)

  CALL add_stream_element (scav, 'WD_HNO3_d', wd_hno3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Large-scale wet deposition HNO3', &
                           units='kg/m2/s', code=106)

  CALL add_stream_element (scav, 'WD_O3_d', wd_o3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Large-scale wet deposition O3', &
                           units='kg/m2/s', code=106)

  CALL add_stream_element (scav, 'WDCV_CH3O2H_d', wdcv_ch3o2h_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Convective wet deposition CH3O2H', &
                           units='kg/m2/s', code=107)

  CALL add_stream_element (scav, 'WDCV_CH2O_d', wdcv_ch2o_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Convective wet deposition CH2O', &
                           units='kg/m2/s', code=107)

  CALL add_stream_element (scav, 'WDCV_H2O2_d', wdcv_h2o2_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Convective wet deposition H2O2', &
                           units='kg/m2/s', code=107)

  CALL add_stream_element (scav, 'WDCV_HNO3_d', wdcv_hno3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Convective wet deposition HNO3', &
                           units='kg/m2/s', code=107)

  CALL add_stream_element (scav, 'WDCV_O3_d', wdcv_o3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Convective wet deposition O3', &
                           units='kg/m2/s', code=107)

  !
  CALL add_stream_element (scav, 'WOHNO3', wohno3, lpost=.TRUE., &
                           longname='Washout HNO3', leveltype = HYBRID, &
                           units='cm-3 s-1', code=107)

  CALL add_stream_element (scav, 'WOHNO3_d', wohno3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Washout HNO3', &
                           leveltype = HYBRID, units='cm-3 s-1', code=107)

END SUBROUTINE construct_stream_scav

!---------------------------------------------------------------------------------

SUBROUTINE destruct_stream_scav

  ! Deallocates memory.

  ! *destruct_stream_* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(scav)

END SUBROUTINE destruct_stream_scav

!---------------------------------------------------------------------------------

SUBROUTINE init_stream_scav

  ! Initializes streams with zero.

  ! *init_stream_* is called from *call_init_tracers*,
  ! src/call_submodels.f90.

  IF (lstart) THEN     ! use restart fields otherwise
     wd_ch3o2h(:,:) = 0._dp
     wd_ch2o(:,:) = 0._dp
     wd_h2o2(:,:) = 0._dp
     wd_hno3(:,:) = 0._dp
     wd_o3(:,:) = 0._dp
     
     wd_ch3o2h_d(:,:) = 0._dp
     wd_ch2o_d(:,:) = 0._dp
     wd_h2o2_d(:,:) = 0._dp
     wd_hno3_d(:,:) = 0._dp
     wd_o3_d(:,:) = 0._dp
     
     wdcv_ch3o2h(:,:) = 0._dp
     wdcv_ch2o(:,:) = 0._dp
     wdcv_h2o2(:,:) = 0._dp
     wdcv_hno3(:,:) = 0._dp
     wdcv_o3(:,:) = 0._dp
     
     wdcv_ch3o2h_d(:,:) = 0._dp
     wdcv_ch2o_d(:,:) = 0._dp
     wdcv_h2o2_d(:,:) = 0._dp
     wdcv_hno3_d(:,:) = 0._dp
     wdcv_o3_d(:,:) = 0._dp

     wohno3(:,:,:) = 0._dp
     wohno3_d(:,:,:) = 0._dp   

  ENDIF
 
END SUBROUTINE init_stream_scav

!---------------------------------------------------------------------------------

SUBROUTINE accumulate_stream_scav (kproma, krow)


  ! This subroutine accumulates the current value of a variable at every time 
  ! step, such that finally monthly streams can be calculated.

  ! *accumulate_stream_* is called from *call_diagn*, 
  ! src/call_submodels

  ! Scalar arguments
  INTEGER, INTENT(in) :: kproma, krow


  ! Executable statements:

  wd_ch3o2h(1:kproma,krow) = wd_ch3o2h(1:kproma,krow) + delta_time * &
       wd_ch3o2h_d(1:kproma,krow)

  wd_ch2o(1:kproma,krow) = wd_ch2o(1:kproma,krow) + delta_time * &
       wd_ch2o_d(1:kproma,krow)

  wd_h2o2(1:kproma,krow) = wd_h2o2(1:kproma,krow) + delta_time * &
       wd_h2o2_d(1:kproma,krow)

  wd_hno3(1:kproma,krow) = wd_hno3(1:kproma,krow) + delta_time * &
       wd_hno3_d(1:kproma,krow)

  wd_o3(1:kproma,krow) = wd_o3(1:kproma,krow) + delta_time * &
       wd_o3_d(1:kproma,krow)

  wdcv_ch3o2h(1:kproma,krow) = wdcv_ch3o2h(1:kproma,krow) + delta_time * &
       wdcv_ch3o2h_d(1:kproma,krow)

  wdcv_ch2o(1:kproma,krow) = wdcv_ch2o(1:kproma,krow) + delta_time * &
       wdcv_ch2o_d(1:kproma,krow)

  wdcv_h2o2(1:kproma,krow) = wdcv_h2o2(1:kproma,krow) + delta_time * &
       wdcv_h2o2_d(1:kproma,krow)

  wdcv_hno3(1:kproma,krow) = wdcv_hno3(1:kproma,krow) + delta_time * &
       wdcv_hno3_d(1:kproma,krow)

  wdcv_o3(1:kproma,krow) = wdcv_o3(1:kproma,krow) + delta_time * &
       wdcv_o3_d(1:kproma,krow)

  wohno3(1:kproma,:,krow) = wohno3(1:kproma,:,krow) + wohno3_d(1:kproma,:,krow)

END SUBROUTINE accumulate_stream_scav
          
END MODULE mo_socol_scav

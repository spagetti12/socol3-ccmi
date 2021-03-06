SUBROUTINE socol_sediment (krow, xt) 

  ! Description:
  !
  ! Calculation of sedimentation of 
  ! - HNO3
  ! - PSCs
  
  ! *socol_sediment* is called from *mezon*.
  
  ! Martin Schraner, ETH Zurich, April 2009
   
  USE mo_constants,               ONLY: api
  USE mo_kind,                    ONLY: dp
  USE mo_socol_constants,         ONLY: eps_mr, sigmanat, rhonat, &
                                        sigmaice, rhoice
  USE mo_socol_dimensions,        ONLY: nbdim, nproma, nlev, ntrac
  USE mo_socol_gcmfields,         ONLY: t, prest
  USE mo_socol_grid_calculations, ONLY: zlevb, dens
  USE mo_socol_namelist,          ONLY: hetnat_rmode, hetice_ndens, lscav, lsynth
  USE mo_socol_streams,           ONLY: sedpsc1_d, sedpsc2_d
  USE mo_socol_time_control,      ONLY: time_step_len_chem
  USE mo_socol_tracers,           ONLY: idt_h2o, idt_psc1, idt_psc2, idt_hno3, idt_h2o2
  USE mo_socol_synth_tracers,     ONLY: idt_NH_50W

  IMPLICIT NONE

  ! Subroutine arguments:
  INTEGER, INTENT(in) :: krow  ! Index of current grid point block
  REAL(dp), INTENT(inout) :: xt(nbdim,nlev,ntrac)  ! Chemical species


  ! Executable statements:

  ! eth_as_scav
  ! now handled in socol_scav
  ! 1. Rainout and washout of HNO3:
  IF (.not. lscav) CALL sediment_hno3

  ! 2. Sedimentation of PSC I and PSC II:
  CALL sediment_pscs

CONTAINS

  SUBROUTINE sediment_hno3

    ! Description:

    ! Calculation of rainout and washout of HNO3 for layers below 160 hPa.

    ! Local variables:
    REAL(dp), PARAMETER :: upbound_sediment_hno3 = 160._dp   ! hPa

    ! Executable statements:

    WHERE (prest(1:nproma,:) .GE. upbound_sediment_hno3) &
         xt(1:nproma,:,idt_hno3) = xt(1:nproma,:,idt_hno3) / &
         (1._dp+time_step_len_chem*4.0E-06_dp)

    WHERE (prest(1:nproma,:) .GE. upbound_sediment_hno3) &
         xt(1:nproma,:,idt_h2o2) = xt(1:nproma,:,idt_h2o2) / &
         (1._dp+time_step_len_chem*4.0E-06_dp)

    if (lsynth) &
    WHERE (prest(1:nproma,:) .GE. upbound_sediment_hno3) &
         xt(1:nproma,:,idt_NH_50W) = xt(1:nproma,:,idt_NH_50W) / &
         (1._dp+time_step_len_chem*4.0E-06_dp)
    
  END SUBROUTINE sediment_hno3

    
  SUBROUTINE sediment_pscs

    ! Description:
  
    ! Calculation of sedimentation of PSCs.
    !
    ! Modified sedimentation code with number of NAT fixed, number of ice fixed,
    ! radius of NAT fixed (-> PSC2 fixed (or zero)), radius of ice calculated
    ! from PSC2-concentration. Sedimentation velocity dependent on radius and 
    ! temperature.
    !
    ! *socol_sediment* is called from *mezon*.
    !
    ! May, 1998: original code
    ! Martin Schraner, Beiping Luo, ETH Zurich, August 2005: Modified code 
    !                                                        (!MSHET)
    ! Martin Schraner, ETH Zurich, April 2009: Modified for SOCOLvs3.0
    
    ! Local variables:
    INTEGER :: jl, jk
    REAL(dp) :: dz, nnat, rmodeice, algsignat, algsigice, sedm, tc, eta, c1, &
         facnat1, facnat2, facnat3, facnat4, facice1, facice2, facice3, &
         facice4
    REAL(dp) :: ymi(nlev), ymn(nlev)
    LOGICAL :: lpscs, lfirst
    
 
    ! Executable statements:

    lfirst = .TRUE.
  
    ! Initialize output streams:
    sedpsc1_d(:,:,krow) = 0.0_dp
    sedpsc2_d(:,:,krow) = 0.0_dp

    DO jl = 1, nproma
       ! Do we have some PSC-I and/or PSC-II in the column?
       lpscs = .FALSE.
  
       DO jk = 1, nlev
          IF (xt(jl,jk,idt_psc1) .GT. eps_mr .OR. &
               xt(jl,jk,idt_psc2) .GT. eps_mr) THEN
             lpscs = .TRUE.
             EXIT
          ENDIF             
       ENDDO
  
       ! Exit if there aren't any PSCs in the column:
       IF (.NOT. lpscs) CYCLE

       ! Calculate some constants when this point is passed for the first time:
       IF (lfirst) THEN
          algsignat=LOG(sigmanat)
          facnat1 = 3.0_dp/4.0_dp*117.0_dp*1.661E-24_dp/(api*hetnat_rmode* &
               hetnat_rmode*hetnat_rmode*rhonat* &
               EXP(9._dp/2._dp*algsignat*algsignat))
          facnat2 = 2.0_dp*980.0_dp*1.62_dp/9.0_dp
          facnat3 = 4._dp/3._dp*api*rhonat*(hetnat_rmode**5)*time_step_len_chem* &
               EXP(25._dp/2._dp*algsignat*algsignat)
          facnat4 = 117._dp*1.661E-24_dp*100._dp
          algsigice=LOG(sigmaice)
          facice1 = 3.0_dp/4.0_dp*18._dp*1.661E-24_dp/(hetice_ndens*api*rhoice* &
               EXP(9._dp/2._dp*algsigice*algsigice))
          facice2 = 2.0_dp*980.0_dp*0.91_dp/9.0_dp
          facice3 = hetice_ndens*4._dp/3._dp*api*rhoice*time_step_len_chem* &
               EXP(25._dp/2._dp*algsigice*algsigice)
          facice4 = 18._dp*1.661E-24_dp*100._dp
          lfirst = .FALSE.
       ENDIF

       ! Calculation of gas-phase HNO3 and H2O:
       xt(jl,:,idt_hno3) = xt(jl,:,idt_hno3) - xt(jl,:,idt_psc1)
       xt(jl,:,idt_h2o)  = xt(jl,:,idt_h2o)  - xt(jl,:,idt_psc2)
       
       DO jk = 1, nlev
          ! Layer vertical size [m]:
          dz = (zlevb(jl,jk)-zlevb(jl,jk+1))*1000._dp
          
          ! Preparations for sedimentation-velocity:
          tc = t(jl,jk)-273.15_dp    !Temperature in Celsius
          eta = (1.718_dp+.0049_dp*tc-1.2E-5_dp*tc*tc)*1.0E-4_dp
          
          ! Do we have PSC-I?
          IF (xt(jl,jk,idt_psc1) .GT. eps_mr) THEN
             
             ! Calculation of number density of NAT:
             nnat = xt(jl,jk,idt_psc1)*dens(jl,jk)*facnat1
             
             ! Coefficient for sedimentation-velocity of NAT
             ! (sedvel [cm/s] = c1 * radius [cm] ^ 2):
             c1 = facnat2/eta
             
             ! Mass of NAT sedimented below from level jk (per timestep per 
             ! cm2) (=sedvel*mass*time step) [g/cm2]:
             sedm = nnat*c1*facnat3
             
             ! Number of NAT molecules sedimented below from level jk (per time 
             ! step per cm2), devided by depth of the level [cm] (= loss of NAT 
             ! molcules per time step per cm3) [molecules/cm3]
             ymn(jk) = sedm/(facnat4*dz)           
          ELSE
             ymn(jk) = 0.0_dp
          ENDIF
          
          ! Do we have PSC-II?
          IF (xt(jl,jk,idt_psc2) .GT. eps_mr) THEN

             ! Mode radius of ice [cm]:
             rmodeice = (xt(jl,jk,idt_psc2)*dens(jl,jk)* &
                  facice1)**(1.0_dp/3.0_dp)
             
             ! Coefficient for sedimentation-velocity of ICE
             ! (sedvel [cm/s] = c1 * radius [cm] ^ 2):
             c1 = facice2/eta
             
             ! Mass of ICE sedimented below from level jk (per timestep per 
             ! cm2) (=sedvel*mass*timestep) [g/cm2]:
             sedm = rmodeice**5*c1*facice3
             
             ! Number of NAT molecules sedimented below from level jk (per time
             ! step per cm2), devided by depth of the level [cm] (= loss of NAT 
             ! molcules per timestep per cm3) [molecules/cm3]
             ymi(jk) = sedm/(facice4*dz)           
          ELSE          
             ymi(jk) = 0.0_dp
          ENDIF 
       
       ENDDO
       
       ! Recalculation of the mixing ratio after sedimentation:    
       ! First layer: 
       sedpsc1_d(jl,1,krow) = -ymn(1)           
       xt(jl,1,idt_psc1) = xt(jl,1,idt_psc1) + sedpsc1_d(jl,1,krow)/dens(jl,1)
  
       sedpsc2_d(jl,1,krow) = -ymi(1)
       xt(jl,1,idt_psc2) = xt(jl,1,idt_psc2) + sedpsc2_d(jl,1,krow)/dens(jl,1)
  
       ! Other layers:
       DO jk = 2, nlev              
          sedpsc1_d(jl,jk,krow) = (ymn(jk-1) - ymn(jk))
          xt(jl,jk,idt_psc1) = &
               xt(jl,jk,idt_psc1) + sedpsc1_d(jl,jk,krow)/dens(jl,jk)
  
          sedpsc2_d(jl,jk,krow) = (ymi(jk-1) - ymi(jk))
          xt(jl,jk,idt_psc2) = &
               xt(jl,jk,idt_psc2) + sedpsc2_d(jl,jk,krow)/dens(jl,jk)    
       ENDDO
          
       ! Back calculation of gas-phase water and HNO3:     
       xt(jl,:,idt_hno3) = xt(jl,:,idt_hno3) + xt(jl,:,idt_psc1)
       xt(jl,:,idt_h2o)  = xt(jl,:,idt_h2o)  + xt(jl,:,idt_psc2)  
       
    ENDDO
   
  END SUBROUTINE sediment_pscs
  
END SUBROUTINE socol_sediment

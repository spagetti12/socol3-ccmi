SUBROUTINE mezon (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, ptte, &
     ptm1, pxtm1, pxtte, pqm1, pqte, ktrpwmo, krow, kglat)

  ! Description:
 
  ! Driver of chemistry module MEZON.
  !
  ! *mezon* is called from *call_chem1*, src/call_submodels.f90.   
  !
  ! Martin Schraner, ETH Zurich, April 2009

  USE mo_constants,               ONLY: amw, amd
  USE mo_exception,               ONLY: message
  USE mo_kind,                    ONLY: dp
  USE mo_socol_ch4,               ONLY: ch4_emiss_wetland ! eth_as_ch4
  USE mo_socol_isotope,           ONLY: idt_12c, idt_13c
  USE mo_socol_constants,         ONLY: eps_mr
  USE mo_socol_dimensions,        ONLY: nbdim, nproma, nlev, nlevp1, ntrac
  USE mo_socol_gcmfields,         ONLY: t, p, ph, prest, presb
  USE mo_socol_grid_calculations, ONLY: calculate_boundaries_altitude, &
                                        calculate_solar_angles, calculate_dens, &
                                        calculate_hetpsc_lowlevind
  USE mo_socol_namelist,          ONLY: lchem, lh2o_coupl, &
       spe, gcr, eep,  lch4_wetland, lch4_isotope, & ! eth_as_ch4
       lext_oh, lo3orig, lsynth , lnudg_2h
  USE mo_socol_time_control,      ONLY: l_trigchem, time_step_len_chem
  USE mo_socol_tracers,           ONLY: trac_chemspec, n_trac_chemspec, xtte_chem, &
                                        pos_socol_tracers, idt_h2o, idt_ch4, idt_oh, idt_no ! eth_as_ch4
  USE mo_time_control,            ONLY: delta_time, time_step_len, &
                                        get_month_len, current_date, get_date_components
  USE mo_socol_ionization,        ONLY: socol_gcr_spe_ionization
  ! eth_as_ohclim+
  USE mo_ohclim,                  ONLY: ohclim
  USE mo_mpi,                     ONLY: p_pe
  ! eth_as_ohclim-
  USE mo_socol_ch4,               ONLY: interpolate_ch4_clim, ch4_surface_flux ! eth_as_ch4
  USE mo_socol_o3orig,            ONLY: o3orig_physc
  USE mo_socol_synth_tracers, ONLY: idt_NH_50W
  USE mo_socol_streams,            ONLY: init_stream_nudg_2h

  IMPLICIT NONE

  ! Subroutine arguments:
  INTEGER, INTENT(in) :: krow, kglat, kproma, kbdim, klev, klevp1, ktrac, &
       ktrpwmo(kbdim)
  REAL(dp), INTENT(in) :: paphp1(kbdim,klevp1), papp1(kbdim,klev), &
         ptte(kbdim,klev), ptm1(kbdim,klev)
  REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac), &
       pqm1(kbdim,klev), pqte(kbdim,klev)

  ! Local variables:
  REAL(dp) :: pxt_bchem(kbdim,klev,ktrac), pxt_achem(kbdim,klev,ktrac), dnox(krow,kproma,klev), dhox(krow,kproma,klev)
  REAL(dp), DIMENSION(kbdim,klev)   :: zoh
  INTEGER  :: jl, jk
  INTEGER  :: icurrentday, icurrentmonth, icurrentyear

!WB  CHARACTER(100) :: kroww
!WB  CHARACTER(200) :: kglatt

  ! External subroutines:
  EXTERNAL interpolate_socol_bcond, set_socol_bcond_mixrat, &
       set_socol_bcond_fluxes, chem, socol_sediment, socol_col_o3


  ! Executable statements:

  IF (lchem) THEN

     ! 1. Save dimensions and some fields calculated by the GCM with the names 
     !    used by MEZON (see 6.):

     ! 1.1 Dimensions:
     nbdim  = kbdim
     nproma = kproma
     nlev   = klev
     nlevp1 = klevp1
     ntrac  = ktrac
     
     ! 1.2 Fields calculated by the GCM:
     
     ! Temperature [K]:
     t(1:kproma,:)     = ptm1(1:kproma,:) + ptte(1:kproma,:)*time_step_len
     ! Pressure at full levels [Pa]:
     p(1:kproma,:)     = papp1(1:kproma,:)
     ! Pressure at half levels [Pa]:
     ph(1:kproma,:)    = paphp1(1:kproma,:)
     ! Pressure at full levels [hPa]:
     prest(1:kproma,:) = papp1(1:kproma,:)*0.01_dp
     ! Pressure at half levels [hPa]: 
     presb(1:kproma,:) = paphp1(1:kproma,:)*0.01_dp


     ! 2. Convert water vapour from g/g (used in GCM) to mixing ratio (used in
     !    CTM):
     IF (lh2o_coupl) THEN
        pxtm1(1:kproma,:,idt_h2o) = amd/amw * pqm1(1:kproma,:)
        pxtte(1:kproma,:,idt_h2o) = amd/amw * pqte(1:kproma,:)
     ENDIF

     ! 3. Calculate grid related quantities:

     ! 3.1 Altitude of model layers at model boundaries:
     CALL calculate_boundaries_altitude (krow, kproma, kbdim, klev, klevp1, &
          p, ph, t)
     
     ! 3.2 Density of air molecules:
     CALL calculate_dens (kproma, kbdim, klev, p, t)

     ! 3.3 Solar angles:
     IF (l_trigchem) CALL calculate_solar_angles (krow, kproma, kbdim, klev)

     ! 3.4 Pressure level indices of the lower boundary of NAT and ice PSCs:
     CALL calculate_hetpsc_lowlevind (kproma, kbdim, klev, p, ktrpwmo)

     ! 3.5 Calculate CH4 emissions from wetlands
     
     IF (lch4_wetland) &
          CALL ch4_emiss_wetland(krow, kproma, kbdim) ! eth_as_ch4

     ! 4. Set boundary conditions: (a) before call of MEZON:

     ! 4.1 Boundary conditions given as fluxes:
     CALL set_socol_bcond_fluxes (krow, kproma, kbdim, klev, klevp1, ktrac, &
          ph, pxtm1, pxtte, ktrpwmo, p)

     CALL get_date_components (current_date,  month=icurrentmonth, year=icurrentyear,  &
            day=icurrentday)

     ! 4.2 Mixing ratio boundary conditions:
     CALL set_socol_bcond_mixrat (kproma, kbdim, klev, ktrac, p, pxtm1, &
          pxtte, pqm1, pqte, ktrpwmo, icurrentday, krow)  

     ! 4.3 Call the EEP flux calculation based on the Ap index
     !     parameterization by Baumgaertner (2008), including SPE
     !     influence on HOx & NOx
     IF ((spe) .OR. (eep) .OR. (gcr)) THEN
        CALL socol_gcr_spe_ionization (krow,klev,kbdim,ktrac,prest,kproma,pxtte,pxtm1) ! JGA
     endif

     ! eth_as_ohclim+
     ! 4.3 External OH distribution
     IF (lext_oh) THEN
        zoh(1:kproma,:) = ohclim(krow,kproma,kbdim,klev,ph,p)
        DO jk = 1, klev
           WHERE (p(1:kproma,jk) .gt. 10000._dp) 
              pxtm1(1:kproma,jk,idt_oh) = zoh(1:kproma,jk)
              pxtte(1:kproma,jk,idt_oh) = 0._dp
           END WHERE
        END DO
     END IF
     ! eth_as_ohclim-

     ! Call o3 orig 
     IF (lo3orig) CALL o3orig_physc(1,kproma,kbdim,krow,ktrac,pxtm1,pxtte,p)
     
     IF (l_trigchem) THEN   ! only at chemical time steps     
     
     ! 5. Input / output of MEZON subroutines:
     
        DO jk = 1, n_trac_chemspec   !FORALL (jk = 1:n_trac_chemspec)
           pxt_bchem(1:nproma,:,trac_chemspec(jk)) = &
                pxtm1(1:nproma,:,trac_chemspec(jk)) + &
                pxtte(1:kproma,:,trac_chemspec(jk))*time_step_len
        END DO   !END FORALL

        IF (lsynth) THEN
           pxt_bchem(1:nproma,:,idt_NH_50W) = &
             pxtm1(1:nproma,:,idt_NH_50W) + &
             pxtte(1:kproma,:,idt_NH_50W)*time_step_len
        END IF

        ! eth_as_ch4+
        IF (lch4_isotope) THEN
           pxt_bchem(1:nproma,:,idt_12c) = &
                pxtm1(1:nproma,:,idt_12c) + &
                pxtte(1:kproma,:,idt_12c)*time_step_len
           
           pxt_bchem(1:nproma,:,idt_13c) = &
                pxtm1(1:nproma,:,idt_13c) + &
             pxtte(1:kproma,:,idt_13c)*time_step_len
        END IF
        ! eth_as_ch4-

        ! Set negative gas concentrations to zero:
        WHERE (pxt_bchem .LT. eps_mr) pxt_bchem = 0.0_dp
        
        pxt_achem(1:nproma,:,:) = pxt_bchem(1:nproma,:,:)

!!$     ! Call o3 orig 
!!$     IF (lo3orig) CALL o3orig_physc(1,kproma,kbdim,krow,ktrac,pxtm1,pxtte,p)
   
     ! 6. Call MEZON subroutines:

     
     ! 6.1 Calculate chemistry:
     CALL chem (krow, pxt_achem, ph)

     ! 6.2 Calculate sedimentation of PSCs and HNO3:
     CALL socol_sediment (krow, pxt_achem)
     
     ! 6.3 Calculate ozone column:
     CALL socol_col_o3 (krow, pxt_achem)
     
     ! 7. Update tracer tendencies:
     
     ! 7.1 Calculate tendency due to chemistry:
     DO jk = 1, n_trac_chemspec   !FORALL (jk = 1:n_trac_chemspec)
        xtte_chem(1:nproma,:,trac_chemspec(jk),krow) = &
             (pxt_achem(1:nproma,:,trac_chemspec(jk)) - &
             pxt_bchem(1:nproma,:,trac_chemspec(jk)))/time_step_len_chem
     END DO   !END FORALL



     IF (lsynth) THEN
        xtte_chem(1:nproma,:,idt_NH_50W,krow) = &
             (pxt_achem(1:nproma,:,idt_NH_50W) - &
             pxt_bchem(1:nproma,:,idt_NH_50W))/time_step_len_chem
     END IF
     
     ! eth_as_ch4+
     IF (lch4_isotope) THEN
        xtte_chem(1:nproma,:,idt_12c,krow) = &
             (pxt_achem(1:nproma,:,idt_12c) - &
             pxt_bchem(1:nproma,:,idt_12c))/time_step_len_chem

        xtte_chem(1:nproma,:,idt_13c,krow) = &
             (pxt_achem(1:nproma,:,idt_13c) - &
             pxt_bchem(1:nproma,:,idt_13c))/time_step_len_chem
     END IF
     ! eth_as_ch4-
     
     
  ENDIF ! l_trigchem

  
     ! 7.2 Add chemical tendencies to pxtte:
     DO jk = 1, n_trac_chemspec   !FORALL (jk = 1:n_trac_chemspec)
        pxtte(1:kproma,:,trac_chemspec(jk)) = &
             pxtte(1:kproma,:,trac_chemspec(jk)) + &
             xtte_chem(1:kproma,:,trac_chemspec(jk),krow)
     END DO   !END FORALL

     IF (lsynth) pxtte(1:kproma,:,idt_NH_50W) = &
          pxtte(1:kproma,:,idt_NH_50W) + &
          xtte_chem(1:kproma,:,idt_NH_50W,krow)

     ! Call o3 orig
!     IF (lo3orig) CALL o3orig_physc(2,kproma,kbdim,krow,ktrac,pxtm1,pxtte,p)

     ! eth_as_ch4+
     IF (lch4_isotope) THEN  
        pxtte(1:kproma,:,idt_12c) = pxtte(1:kproma,:,idt_12c) + &
             xtte_chem(1:kproma,:,idt_12c,krow)
        
        pxtte(1:kproma,:,idt_13c) = pxtte(1:kproma,:,idt_13c) + &
             xtte_chem(1:kproma,:,idt_13c,krow) 
     END IF

     ! eth_as_ch4-   

     ! 8. Set boundary conditions: (b) after call of MEZON:

     ! 8.1 Mixing ratio boundary conditions:
     CALL set_socol_bcond_mixrat (kproma, kbdim, klev, ktrac, p, pxtm1, &
          pxtte, pqm1, pqte, ktrpwmo, icurrentday, krow)

     ! eth_as_ohclim+
     ! 8.3 External OH distribution
     IF (lext_oh) THEN
        DO jk = 1, klev
           WHERE (p(1:kproma,jk) .gt. 10000._dp) 
              pxtm1(1:kproma,jk,idt_oh) = zoh(1:kproma,jk)
              pxtte(1:kproma,jk,idt_oh) = 0._dp
           END WHERE
        END DO
     END IF
     ! eth_as_ohclim-

     ! 9. Check tendencies (force that chemical species remain positive):
     CALL pos_socol_tracers(kproma, kbdim, klev, ktrac, krow, pxtm1, pxtte)

     If (lnudg_2h) CALL init_stream_nudg_2h (nproma,krow,kbdim,ktrac,klev,pxtm1,pxtte,t)



     ! 10. Convert water vapour from mixing ratio (used in CTM) to g/g 
     !     (used in GCM):
     IF (lh2o_coupl) pqte(1:kproma,:) = amw/amd * pxtte(1:kproma,:,idt_h2o)

     ! 11. CH4 nudging

     CALL interpolate_ch4_clim(kglat) ! Interpolate CMDL CH4 climatology to current ts

     CALL ch4_surface_flux(krow, kproma, kbdim, klev, ktrac, pxtm1, pxtte)  
  ENDIF

END SUBROUTINE mezon

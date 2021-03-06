  SUBROUTINE read_socol_bcond_y(yr) 
    
    ! *read_socol_bcond_y* calls subroutines to read the SOCOL boundary
    ! conditions for the current model year.
    ! To be read after (re-)start and at the beginning of every year.

    ! *read_socol_bcond_y* is called from *call_read_bcond*,
    ! src/call_submodels.f90.
    !
    ! M. Schraner, ETH Zurich, November 2008

    USE mo_socol_ghg_ods,         ONLY: read_socol_ghg_ods_global, &
                                        read_socol_odsmem_global
    USE mo_socol_namelist
    USE mo_socol_sun,             ONLY: read_socol_sun_data

    IMPLICIT NONE

    INTEGER :: yr

    ! Executable statements:

    ! Greenhouse gases and ODSs (globally constant fields):
    CALL read_socol_ghg_ods_global(yr)

    ! Preparation of conversion factor for ODS families:
    IF (lchem .OR. lch4_nocoupl_3d .OR. ln2o_nocoupl_3d .OR. &
         lodscl_nocoupl_3d) CALL read_socol_odsmem_global(yr)

    ! Solar irradiance of SW spectral bands, Schumann-Runge-bands, 
    ! Lyman alpha, and NO-photolysis parameters:
    CALL read_socol_sun_data(yr)

  END SUBROUTINE read_socol_bcond_y
!==============================================================================
  SUBROUTINE read_socol_bcond_m(yr,mo,day)
    
    ! *read_socol_bcond_m* calls subroutines to read the SOCOL boundary
    ! conditions for the current model month.
    ! To be read after (re-)start and at the beginning of every month.

    ! *read_socol_bcond_m* is called from *call_read_bcond*,
    ! src/call_submodels.f90.
    !
    ! M. Schraner, ETH Zurich, November 2008

    USE mo_radiation,             ONLY: iaero
    USE mo_socol_co_nox,          ONLY: read_socol_co_nox
    USE mo_socol_nmvoc,           ONLY: read_socol_nmvoc   ! eth_as_nmvoc
    USE mo_socol_ch4,             ONLY: read_socol_ch4     ! eth_as_ch4
    USE mo_socol_ghg_ods,         ONLY: read_socol_ghg_3d
    USE mo_socol_namelist
    USE mo_socol_qbo,             ONLY: read_socol_qbo
    USE mo_socol_sun,             ONLY: read_socol_photolysis
    USE mo_socol_strataerosols,   ONLY: read_socol_strataerosols
    USE mo_socol_tropoaerosols,   ONLY: read_socol_tropoaerosols
    USE mo_socol_ionization,      ONLY: read_socol_ionization_data
    USE mo_socol_photo,           ONLY: read_socol_photo   ! eth_as_tropchem

    IMPLICIT NONE

    INTEGER :: yr, mo, day
    LOGICAL :: newmonth

    ! Executable statements:

    ! Greenhouse gases from former simulation (3d-fields):
    ! (Initialization of "normal" GHG/ODS boundary conditions is called from
    ! *read_socol_bcond_y*.)
    IF ((.NOT. lch4_coupl .AND. lch4_nocoupl_3d) .OR. &
         (.NOT. ln2o_coupl .AND. ln2o_nocoupl_3d) .OR. &
         (.NOT. lodscl_coupl .AND. lodscl_nocoupl_3d)) &
         CALL read_socol_ghg_3d(yr,mo)

    ! Stratospheric aerosols:
    CALL read_socol_strataerosols(yr,mo)

    ! Tropospheric aerosols:
    IF (iaero .EQ. 10) CALL read_socol_tropoaerosols(yr,mo)

    ! CO and NOx emissions:
    IF (lchem) CALL read_socol_co_nox(yr,mo)

    ! eth_as_nmvoc
    ! NMVOC emissions:
    IF (lchem) CALL read_socol_nmvoc(yr,mo)

    ! eth_as_ch4
    ! CH4 emissions:
    IF (lch4_flux) CALL read_socol_ch4(yr,mo)

    ! Photolysis rates:
    IF (lchem) CALL read_socol_photolysis(yr,mo)

    ! QBO nudging:
    CALL read_socol_qbo(yr,mo)

    ! GCRs and Ap-Index which is needed for the HOx/NOx production
    ! through ionization of the Mesosphere  ! JGAFUPSOL

    IF (gcr .OR. eep .OR. spe) THEN
       newmonth = .TRUE.
       CALL read_socol_ionization_data(yr,mo,day,newmonth)
    ENDIF

    ! eth_as_tropchem
!!$    ! Precalculated photolysis rates:
!!$    IF (lchem) CALL read_socol_photo(yr,mo)

  END SUBROUTINE read_socol_bcond_m
!==============================================================================
  SUBROUTINE read_socol_bcond_d(yr,mo,dy)

    ! *read_socol_bcond_d* calls subroutines to read the SOCOL boundary
    ! conditions for the current model month.
    ! To be read after (re-)start and at the beginning of every month.

    ! *read_socol_bcond_d* is called from *call_read_bcond*,
    ! src/call_submodels.f90.
    !
    ! M. Schraner, ETH Zurich, November 2008

    USE mo_socol_ionization,      ONLY: read_socol_ionization_data
    USE mo_socol_namelist,        ONLY: gcr, spe, eep

    IMPLICIT NONE

    INTEGER :: yr, mo, dy
    LOGICAL :: newmonth

    newmonth = .TRUE.
    ! Executable statements:

    ! Get only the SPE/EEP info which is needed for the HOx/NOx production
    ! through ionization of the Mesosphere  ! JGASOCOL

   IF (gcr .OR. eep .OR. spe) THEN
       newmonth=.FALSE.
       CALL read_socol_ionization_data(yr,mo,dy,newmonth)
   ENDIF

  END SUBROUTINE read_socol_bcond_d
!==============================================================================
  SUBROUTINE interpolate_socol_bcond_global

    ! *interpolate_socol_bcond_global* calls subroutines to interpolate 
    ! the globally constant boundary conditions of SOCOL 
    ! to the current time step.

    ! *interpolate_socol_bcond_global* is called from
    ! *call_read_bcond*, src/call_submodels.f90.
    !
    ! M. Schraner, ETH Zurich, November 2008

    USE mo_socol_ghg_ods,         ONLY: interpolate_socol_ghg_ods_glob
    USE mo_socol_namelist
    USE mo_socol_sun,             ONLY: interpolate_socol_sun_data, &
                                        interpolate_socol_photolysis
    USE mo_socol_time_control,    ONLY: l_trigchem
    USE mo_time_control,          ONLY: l_trigrad

    IMPLICIT NONE

    ! Executable statements:

    ! Greenhouse gases and ODS's (1d fields):
    IF (l_trigrad .OR. lchem) CALL interpolate_socol_ghg_ods_glob

    ! Solar irradiance of SW spectral bands, Schumann-Runge-bands, 
    ! Lyman alpha, and NO-photolysis parameters:
    IF (l_trigrad .OR. l_trigchem) CALL interpolate_socol_sun_data

    ! Photolysis rates:
    IF (l_trigchem) CALL interpolate_socol_photolysis 

  END SUBROUTINE interpolate_socol_bcond_global
!==============================================================================
  SUBROUTINE interpolate_socol_bcond (krow, kproma, kbdim, klev, klevp1, p, &
       ph, t, ptropo, ktrpwmop1)

    ! *interpolate_socol_bcond* calls subroutines to interpolate 
    ! the boundary conditions of SOCOL (except for globally constant boundary
    ! conditions) to the current time step (and to the current pressure grid, 
    ! if necessary).

    ! *interpolate_socol_bcond* is called from *physc*.
    !
    ! M. Schraner, ETH Zurich, November 2008

    USE mo_kind,                    ONLY: dp
    USE mo_radiation,               ONLY: iaero
    USE mo_socol_co_nox,            ONLY: interpolate_socol_co_nox
    USE mo_socol_nmvoc,             ONLY: interpolate_socol_nmvoc      ! eth_as_nmvoc
    USE mo_socol_ch4,               ONLY: interpolate_socol_ch4        ! eth_as_ch4
    USE mo_socol_ghg_ods,           ONLY: interpolate_socol_ghg_ods_3d, &
                                          interpolate_socol_ods_fa2mem, &
                                          interpolate_socol_co2trac
    USE mo_socol_grid_calculations, ONLY: calculate_zaetr_socol, &
                                          calculate_boundaries_altitude
    USE mo_socol_namelist
    USE mo_socol_strataerosols,     ONLY: interpolate_socol_strataerosols
    USE mo_socol_time_control,      ONLY: l_trigchem
    USE mo_socol_tropoaerosols,     ONLY: interpolate_socol_tropoaerosols
    USE mo_time_control,            ONLY: l_trigrad
    USE mo_socol_photo,           ONLY: interpolate_socol_photo  ! eth_as_tropchem

    IMPLICIT NONE

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim, klev, klevp1
    INTEGER, INTENT(in)  :: ktrpwmop1(kbdim) ! Pressure index of uppermost 
                            ! level with level center below tropopause
    REAL(dp), INTENT(in) :: p(kbdim,klev), ph(kbdim,klevp1) 
                            ! Pressure at full/half levels [Pa]
    REAL(dp), INTENT(in) :: t(kbdim,klev)    ! Temperature [K]
    REAL(dp), INTENT(in) :: ptropo(kbdim)    ! Tropopause height [Pa]


    ! Executable statements:

    IF (l_trigrad .OR. lchem) THEN

       IF (l_trigrad .OR. l_trigchem) THEN
          ! 1. Preparations:
          
          ! 1.1 Calculate stratosphere/ troposphere flag with continuos  
          !     transition below tropopause ("zaetr"):
          CALL calculate_zaetr_socol(kproma, kbdim, klev, ktrpwmop1)
       
          ! 1.2 Calculate altitude of model layers at model boundaries:
          CALL calculate_boundaries_altitude(krow, kproma, kbdim, klev, &
               klevp1, p, ph, t)
       ENDIF
       
       ! 2. Linear interpolation to chemistry / radiation time step:
       
       ! 2.1 Greenhouse gases from 3d fields (if not coupled to chemistry 
       !     module):
       IF (l_trigrad .AND. (lch4_nocoupl_3d .OR. ln2o_nocoupl_3d .OR. &
            lodscl_nocoupl_3d)) &
            CALL interpolate_socol_ghg_ods_3d(krow, kproma, kbdim, klev, p)
       
       ! 2.2 Conversion factor for ODS families:
       IF (l_trigchem .OR. (l_trigrad .AND. (lch4_nocoupl_3d .OR. &
            lodscl_nocoupl_3d))) &
            CALL interpolate_socol_ods_fa2mem(krow, kproma, kbdim, klev, p, &
            ptropo)
       
       ! 2.3 Stratospheric aerosols:
       IF (l_trigrad .OR. l_trigchem) &
            CALL interpolate_socol_strataerosols(krow, kproma, kbdim, klev)
       
       ! 2.4 Tropospheric aerosols:
       IF (l_trigrad .AND. iaero .EQ. 10) &
            CALL interpolate_socol_tropoaerosols(krow, kproma, kbdim)
       
       ! 2.5 CO and NOx emissions:
       IF (lchem) CALL interpolate_socol_co_nox(krow, kproma, kbdim)

       ! 2.6 CO2-tracer
       IF (lpco2) CALL interpolate_socol_co2trac(krow, kproma, kbdim)

       ! eth_as_ch4
       ! 2.6 CH4 emissions:
       IF (lch4_flux) CALL interpolate_socol_ch4(krow, kproma, kbdim)   

       ! eth_as_nmvoc
       ! 2.7 NMVOC emissions:
       IF (lchem) CALL interpolate_socol_nmvoc(krow, kproma, kbdim)  

       ! eth_as_tropchem
!!$       ! 2.8 Precalculated photolysis rates:
!!$       IF (l_trigchem) CALL interpolate_socol_photo(krow, kproma, kbdim, klev)  

    ENDIF
       
  END SUBROUTINE interpolate_socol_bcond
!==============================================================================
  SUBROUTINE set_socol_bcond_mixrat(kproma, kbdim, klev, ktrac, p, pxtm1, &
     pxtte, pqm1, pqte, ktrpwmo, day, krow) 

    ! Sets mixing ratio boundary conditions.

    ! *set_socol_bcond_mixgrat* is called from *socol_mezon*.
    !
    ! M. Schraner, ETH Zurich, January 2009

    USE mo_constants,             ONLY: amd, amw
    USE mo_kind,                  ONLY: dp
    USE mo_socol_ghg_ods,         ONLY: co2_bcond_chem, ch4_bcond_chem, &
                                        n2o_bcond_chem, odscls_bcond_chem, &
                                        odscll_bcond_chem, odsbr_bcond_chem, &
                                        hcl_bcond_chem, clno3_bcond_chem, &
                                        hbr_bcond_chem, co2_chem, &
                                        co2_trac, &
                                        odsmem_bcond_chem
    USE mo_socol_namelist,        ONLY: lh2o_coupl, nlevabtropo_h2ogcm, lpco2, lch4_flux, &
                                        lheppa
    USE mo_socol_time_control,    ONLY: l_trigchem
    USE mo_socol_tracers,         ONLY: idt_ch4, idt_n2o, idt_odscls, &
                                        idt_odscll, idt_odsbr, idt_hcl, &
                                        idt_clno3, idt_hbr, idt_h2o, &
                                        idt_o3, idt_h2o2, idt_ch3o2h, &
                                        idt_f11, idt_f12, idt_cbrf3, idt_cfc113, idt_cfc114, &
                                        idt_cfc115, idt_ccl4, idt_ch3ccl3, idt_hcfc22, idt_hcfc141b, &
                                        idt_hcfc142b, idt_h1211, idt_ch3br, idt_ch3cl, idt_hcfc21, &
                                        idt_hcfc123, idt_h2402, idt_chbr3, idt_ch2br2, &
                                        idt_co2, idt_linearage, idt_idealage, idt_no, idt_no2, &
                                        idt_n, idt_no3, idt_n2o5, &
                                        idt_hno3,idt_hno4
    USE mo_socol_co_nox,            ONLY: nox_heppa

    USE mo_time_control,          ONLY: current_date, next_date, start_date, get_date_components

    IMPLICIT NONE

    ! Subroutine arguments:
    INTEGER, INTENT(in)     :: kproma, kbdim, klev, ktrac, krow
    REAL(dp), INTENT(in) :: p(kbdim,klev)   ! Pressure of full levels
    REAL(dp), INTENT(in) :: pqm1(kbdim,klev), pqte(kbdim,klev)
    REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)
    INTEGER, INTENT(in)  :: ktrpwmo(kbdim)   ! Pressure index of pressure level
                                             ! containing the tropopause 

    ! Local variables:
    INTEGER :: i, jl, jk, klev0
    INTEGER :: monlen, day

    INTEGER :: icurrentyear, icurrentmonth, icurrentday
    INTEGER :: istartyear, istartmonth, istartday
    INTEGER :: inextyear, inextmonth, inextday
    
    LOGICAL :: lnewday

    ! Local parameters:
    ! Number of lowermost levels where mixing ratios are prescribed by lower 
    ! boundary condition:
    INTEGER, PARAMETER :: nlev_lowbound_bcond_mixrat = 6

    ! Upper limits for O3, H2O2, and CH3O2H:
    REAL(dp), PARAMETER :: o3_upperlimit_chem     = 1.5E-05_dp
    REAL(dp), PARAMETER :: h2o2_upperlimit_chem   = 8.0E-09_dp
    REAL(dp), PARAMETER :: ch3o2h_upperlimit_chem = 2.5E-09_dp

!WB    CHARACTER(100) :: kroww
!WB    CHARACTER(100) :: kpromaa

    ! Executable statements:

    ! Lower boundary conditions:
    klev0 = klev-nlev_lowbound_bcond_mixrat+1

    IF (.not. lch4_flux) &
         pxtm1(1:kproma,klev0:klev,idt_ch4)    = ch4_bcond_chem    ! eth_as_ch4
    pxtm1(1:kproma,klev0:klev,idt_n2o)    = n2o_bcond_chem
    pxtm1(1:kproma,klev0:klev,idt_odscls) = odscls_bcond_chem
    pxtm1(1:kproma,klev0:klev,idt_odscll) = odscll_bcond_chem
    pxtm1(1:kproma,klev0:klev,idt_odsbr)  = odsbr_bcond_chem
    pxtm1(1:kproma,klev0:klev,idt_hcl)    = hcl_bcond_chem
    pxtm1(1:kproma,klev0:klev,idt_clno3)  = clno3_bcond_chem
    pxtm1(1:kproma,klev0:klev,idt_hbr)    = hbr_bcond_chem
    pxtm1(1:kproma,klev0:klev,idt_f11)      = odsmem_bcond_chem(1)
    pxtm1(1:kproma,klev0:klev,idt_f12)      = odsmem_bcond_chem(2)
    pxtm1(1:kproma,klev0:klev,idt_cbrf3)    = odsmem_bcond_chem(12)
    pxtm1(1:kproma,klev0:klev,idt_cfc113)   = odsmem_bcond_chem(3)
    pxtm1(1:kproma,klev0:klev,idt_cfc114)   = odsmem_bcond_chem(4)
    pxtm1(1:kproma,klev0:klev,idt_cfc115)   = odsmem_bcond_chem(5)
    pxtm1(1:kproma,klev0:klev,idt_ccl4)     = odsmem_bcond_chem(6)
    pxtm1(1:kproma,klev0:klev,idt_ch3ccl3)  = odsmem_bcond_chem(7)
    pxtm1(1:kproma,klev0:klev,idt_hcfc22)   = odsmem_bcond_chem(8)
    pxtm1(1:kproma,klev0:klev,idt_hcfc141b) = odsmem_bcond_chem(9)
    pxtm1(1:kproma,klev0:klev,idt_hcfc142b) = odsmem_bcond_chem(10)
    pxtm1(1:kproma,klev0:klev,idt_h1211)    = odsmem_bcond_chem(11)
    pxtm1(1:kproma,klev0:klev,idt_ch3br)    = odsmem_bcond_chem(13)
    pxtm1(1:kproma,klev0:klev,idt_ch3cl)    = odsmem_bcond_chem(14)
    pxtm1(1:kproma,klev0:klev,idt_hcfc21)   = odsmem_bcond_chem(15)
    pxtm1(1:kproma,klev0:klev,idt_hcfc123)  = odsmem_bcond_chem(16)
    pxtm1(1:kproma,klev0:klev,idt_h2402)    = odsmem_bcond_chem(17)
    pxtm1(1:kproma,klev0:klev,idt_chbr3)    = odsmem_bcond_chem(18)
    pxtm1(1:kproma,klev0:klev,idt_ch2br2)   = odsmem_bcond_chem(19)

    IF (lpco2) THEN
       do i = klev0, klev
          pxtm1(1:kproma,i,idt_co2) = co2_trac(1:kproma)
       end do
    END IF

    ! Set tendencies of above species to zero in the lowermost levels:
    IF (.not. lch4_flux) &
         pxtte(1:kproma,klev0:klev,idt_ch4)    = 0._dp             ! eth_as_ch4
    pxtte(1:kproma,klev0:klev,idt_n2o)    = 0._dp
    pxtte(1:kproma,klev0:klev,idt_odscls) = 0._dp 
    pxtte(1:kproma,klev0:klev,idt_odscll) = 0._dp 
    pxtte(1:kproma,klev0:klev,idt_odsbr)  = 0._dp 
    pxtte(1:kproma,klev0:klev,idt_hcl)    = 0._dp 
    pxtte(1:kproma,klev0:klev,idt_clno3)  = 0._dp 
    pxtte(1:kproma,klev0:klev,idt_hbr)    = 0._dp 
    pxtte(1:kproma,klev0:klev,idt_f11)      = 0._dp
    pxtte(1:kproma,klev0:klev,idt_f12)      = 0._dp
    pxtte(1:kproma,klev0:klev,idt_cbrf3)    = 0._dp
    pxtte(1:kproma,klev0:klev,idt_cfc113)   = 0._dp
    pxtte(1:kproma,klev0:klev,idt_cfc114)   = 0._dp
    pxtte(1:kproma,klev0:klev,idt_cfc115)   = 0._dp
    pxtte(1:kproma,klev0:klev,idt_ccl4)     = 0._dp
    pxtte(1:kproma,klev0:klev,idt_ch3ccl3)  = 0._dp
    pxtte(1:kproma,klev0:klev,idt_hcfc22)   = 0._dp
    pxtte(1:kproma,klev0:klev,idt_hcfc141b) = 0._dp
    pxtte(1:kproma,klev0:klev,idt_hcfc142b) = 0._dp
    pxtte(1:kproma,klev0:klev,idt_h1211)    = 0._dp
    pxtte(1:kproma,klev0:klev,idt_ch3br)    = 0._dp
    pxtte(1:kproma,klev0:klev,idt_ch3cl)    = 0._dp
    pxtte(1:kproma,klev0:klev,idt_hcfc21)   = 0._dp
    pxtte(1:kproma,klev0:klev,idt_hcfc123)  = 0._dp
    pxtte(1:kproma,klev0:klev,idt_h2402)    = 0._dp
    pxtte(1:kproma,klev0:klev,idt_chbr3)    = 0._dp
    pxtte(1:kproma,klev0:klev,idt_h2402)    = 0._dp

    IF (lpco2) pxtte(1:kproma,klev0:klev,idt_co2)    = 0._dp

    ! idealized age tracers
    CALL get_date_components (start_date, month=istartmonth, year=istartyear, &
         day=istartday)

    CALL get_date_components (current_date, month=icurrentmonth, year=icurrentyear, &
         day=icurrentday)

    CALL get_date_components (next_date, month=inextmonth, year=inextyear, &
         day=inextday)

    ! linearly increasing age tracer
    ! surface mixing ratio increases by 1 every day

    lnewday=icurrentday/=inextday

    IF (lnewday) THEN
       pxtm1(1:kproma,klev,idt_linearage)    = pxtm1(1:kproma,klev,idt_linearage) + 1._dp
    ENDIF

    pxtte(1:kproma,klev,idt_linearage)       = 0._dp

    ! "ideal age" tracer

    pxtm1(1:kproma,klev,idt_idealage)       = 0._dp  ! surface mixing ratio equals 0 all time
    pxtte(1:kproma,klev,idt_idealage)       = 0._dp
    
    IF (lnewday) THEN
       pxtm1(1:kproma,1:klev-1,idt_idealage)    = pxtm1(1:kproma,1:klev-1,idt_idealage) + 1._dp  ! mixing ratio increases by 1 every day
    ENDIF
    
    ! Boundary condition for CO2 (globally constant value):
    IF (l_trigchem) co2_chem = co2_bcond_chem

    ! Replace water vapour of the CTM by the one calculated by the GCM for 
    ! levels below seplevind_h2o (!! only in case that lcoupl_h2o=.FALSE. !!):
    IF (.NOT. lh2o_coupl) THEN
       DO jl = 1, kproma
          klev0 = ktrpwmo(jl) - nlevabtropo_h2ogcm
          pxtm1(jl,klev0:klev,idt_h2o) = &
               amd/amw * pqm1(jl,klev0:klev)  ! g/g -> mixing ratio
          pxtte(jl,klev0:klev,idt_h2o) = &
               amd/amw * pqte(jl,klev0:klev)  ! g/g -> mixing ratio
       ENDDO
    ENDIF

    ! Upper limits for certain species:
    WHERE(pxtm1(1:kproma,:,idt_o3) .GT. o3_upperlimit_chem)
       pxtm1(1:kproma,:,idt_o3) = o3_upperlimit_chem
       pxtte(1:kproma,:,idt_o3) = 0._dp
    ENDWHERE
    WHERE(pxtm1(1:kproma,:,idt_h2o2) .GT. h2o2_upperlimit_chem)
       pxtm1(1:kproma,:,idt_h2o2) = h2o2_upperlimit_chem
       pxtte(1:kproma,:,idt_h2o2) = 0._dp
    ENDWHERE
    WHERE(pxtm1(1:kproma,:,idt_ch3o2h) .GT. ch3o2h_upperlimit_chem)
       pxtm1(1:kproma,:,idt_ch3o2h) = ch3o2h_upperlimit_chem
       pxtte(1:kproma,:,idt_ch3o2h) = 0._dp
    ENDWHERE


  If (lheppa) then
     IF (lheppa.and.icurrentday<=inextday) then

!WB     WRITE(kroww,  '(i8)') krow
!WB          kroww = TRIM(adjustl(kroww))
!WB     WRITE(kpromaa,  '(i8)') kproma
!WB          kpromaa = TRIM(adjustl(kpromaa))

!WB     write(*,*)shape(pxtm1)
!WB     write(*,*)idt_no2
!WB     write(*,*)idt_no
!WB     open(1338,file='no_bc_chk_pxtm1b_'//TRIM(kroww)//'_'//TRIM(kpromaa)//'.txt')
!WB     write(1338,*)pxtm1
!WB     close(1338)

!WB     write(*,*)shape(nox_heppa)
!WB     open(1339,file='no_bc_chk_nox_heppa_'//TRIM(kroww)//'_'//TRIM(kpromaa)//'.txt')
!WB     write(1339,*)nox_heppa
!WB     close(1339)

       IF ((day==1).and.(icurrentmonth==11)) THEN
       else

       DO jl=1,kproma
!WB       pxtm1(jl,1,idt_no2)  = pxtm1(jl,2,idt_no2)*nox_heppa(jl,day,krow)/ &
!WB                        (pxtm1(jl,2,idt_no2)+pxtm1(jl,2,idt_no))
!WB       pxtm1(jl,1,idt_no)  = pxtm1(jl,2,idt_no)*nox_heppa(jl,day,krow)/ &
!WB                        (pxtm1(jl,2,idt_no2)+pxtm1(jl,2,idt_no))

         pxtm1(jl,1,idt_n)  = 0._dp!
         pxtm1(jl,1,idt_n2o5)  = 0._dp!
         pxtm1(jl,1,idt_no3)  = 0._dp!
         pxtm1(jl,1,idt_hno3)  = 0._dp!
         pxtm1(jl,1,idt_hno4)  = 0._dp!
         pxtm1(jl,1,idt_no2)  = 0._dp!
         pxtm1(jl,1,idt_no)   = nox_heppa(jl,day,krow)

       ENDDO

!WB     write(*,*)shape(pxtte)
!WB     open(1341,file='no_bc_chk_pxtteb_'//TRIM(kroww)//'_'//TRIM(kpromaa)//'.txt')
!WB     write(1341,*)pxtte
!WB     close(1341)

        pxtte(1:kproma,1,idt_n)    = 0._dp!
        pxtte(1:kproma,1,idt_n2o5)    = 0._dp!
        pxtte(1:kproma,1,idt_no3)    = 0._dp!
        pxtte(1:kproma,1,idt_hno3)    = 0._dp!
        pxtte(1:kproma,1,idt_hno4)    = 0._dp!
        pxtte(1:kproma,1,idt_no2)    = 0._dp!
        pxtte(1:kproma,1,idt_no)     = 0._dp!


!WB     write(*,*)shape(pxtm1)
!WB     open(1340,file='no_bc_chk_pxtm1a_'//TRIM(kroww)//'_'//TRIM(kpromaa)//'.txt')
!WB     write(1340,*)pxtm1
!WB     close(1340)

!WB     write(*,*)shape(pxtte)
!WB     open(1342,file='no_bc_chk_pxttea_'//TRIM(kroww)//'_'//TRIM(kpromaa)//'.txt')
!WB     write(1342,*)pxtte
!WB     close(1342)

      endif
     endif
    endif

  END SUBROUTINE set_socol_bcond_mixrat
!===============================================================================
  SUBROUTINE set_socol_bcond_fluxes(krow, kproma, kbdim, klev, klevp1, ktrac,  &
       ph, pxtm1, pxtte, ktrpwmo, papp1) 

    ! Prescribes boundary conditions given as fluxes. 
    ! - The surface emissions are added uniformely to the 
    !   <nlev_lowbound_bcond_fluxes> (see below) lowermost levels. 
    ! - The removal of species by deposition acts uniformly to the
    !   <nlev_lowbound_bcond_fluxes> (see below) lowermost levels. 
    ! - The lighning and aircraft emissions are added locally to the
    !   corresponding grid boxes.

    ! *set_socol_bcond_fluxes* is called from *mezon*.
    !
    ! M. Schraner, ETH Zurich, February 2009 

    USE mo_constants,               ONLY: g, amd
    USE mo_kind,                    ONLY: dp
    USE mo_exception,               ONLY: finish, message, message_text ! eth_as_ch4
    USE mo_memory_g3b,              ONLY: slm  !land-sea mask (1: land, 0: sea)
    USE mo_socol_constants,         ONLY: atomweight, amco, amno, amno2, amch4, &
                                          nofracnox_bcond, no2fracnox_bcond, &
                                          amc5h8, amch2o, amhcooh, amch3cooh ! eth_as_tropchem
    USE mo_socol_co_nox,            ONLY: co_emiss_surf, nox_emiss_surf, &
                                          nox_emiss_aircr, nox_lightn_col, &
                                          nox_lightn_profile, &
                                          k0_aircr, k1_aircr, k0_lightn, &
                                          k1_lightn
    USE mo_socol_nmvoc,             ONLY: nmvoc_emiss_biogen, nmvoc_emiss_bb, &
                                          nmvoc_emiss_anthrop, & ! eth_as_nmvoc
                                          c5h8_emiss, ch2o_emiss, hcooh_emiss, & ! eth_as_tropchem
                                          ch3cooh_emiss ! eth_as_tropchem
    USE mo_socol_ch4,               ONLY: ch4_emiss_ship, ch4_emiss_grassfire, & ! eth_as_ch4
                                          ch4_emiss_forestfire, ch4_emiss_agr, &
                                          ch4_emiss_awb, ch4_emiss_dom, &
                                          ch4_emiss_ene, ch4_emiss_ind, &
                                          ch4_emiss_tra, ch4_emiss_wst, &
                                          ch4_emiss_slv, ch4_emiss_wet, &
                                          ch4_emiss_wetlands, ch4_emiss_bb, &
                                          ch4_emiss_avi, ch4_emiss_air, &        
                                          ch4_emiss_trans, ch4_emiss_intship, &
                                          ch4_emiss_resi, ch4_emiss_fuels, &      
                                          ch4_emiss_agri, ch4_emiss_waste, &      
                                          ch4_emiss_wstwat, ch4_emiss_other   
    USE mo_socol_deposition,        ONLY: wcodep, wnodep, wno2dep, wo3dep, &
                                          whno3dep, wh2o2dep, &
                                          wh2dep ! eth_as_20100527
    USE mo_socol_grid_calculations, ONLY: zlevb, dens
    USE mo_socol_tracers,           ONLY: idt_co, idt_no, idt_no2, &
                                          idt_o3, idt_hno3, idt_h2o2, &
                                          idt_h2, &  ! eth_as_20100527
                                          idt_ch4, &  ! eth_as_ch4
                                          idt_c5h8, idt_ch2o, idt_hcooh, idt_ch3cooh ! eth_as_tropchem
    USE mo_time_control,            ONLY: time_step_len
    USE mo_socol_namelist,          ONLY: lch4_flux, lch4_wetland, lch4_isotope, &  ! eth_as_ch4
                                          lch4_ipcc, ch4emissfac, &
                                          lnmvoc, lnmvoc_accmip, coemissfac, & ! eth_as_co
                                          interactivelnox,&
                                          lsynth !CCMI tracers
    USE mo_socol_isotope                                                         ! eth_as_ch4
    ! eth_as_nmvoc+
    USE mo_socol_streams,           ONLY: elnox, lnox_d
    ! eth_as_nmvoc-
    USE mo_time_control,            ONLY: time_step_len,lstart
    USE mo_geoloc,                  ONLY: philon_2d, philat_2d
    USE mo_socol_lnox_main,         ONLY: PLT, lnt1, llnox_calc
    USE mo_socol_synth_tracers,     ONLY: synth_trac_flux  !CCMI tracers

    IMPLICIT NONE

    ! Subroutine arguments:
    INTEGER, INTENT(in)     :: krow, kproma, kbdim, klev, klevp1, ktrac, ktrpwmo(kbdim)
    REAL(dp), INTENT(in)    :: ph(kbdim,klevp1) 
                               ! Pressure at half levels [Pa]
    REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac)
    REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)    
    REAL(dp), INTENT(in) :: papp1(kbdim,klev) 
    ! Local parameters:
    ! Number of lowermost levels where mixing ratios are prescribed by lower 
    ! boundary condition:
    INTEGER, PARAMETER :: nlev_lowbound_bcond_fluxes = 1 ! 6 ! eth_as_20110112 

    ! Local variables:
    INTEGER  :: jl, jk, klev0, jkp, i, profile_type
    real(dp) :: zplt(16), pltm(klev), zkmh(klev) 
    REAL(dp) :: fac, amnox, scale_ch4
    REAL(dp) :: zanthrop, zbb, zwet ! eth_as_ch4
    REAL(dp), ALLOCATABLE :: airmolec_layer(:,:)
    REAL(dp) :: airmolec_col(kbdim), dz_layer(kbdim), &
         gasmolec_layer(kbdim), dz_col(kbdim), emiss(kbdim), &
         depos(kbdim), deposvel(kbdim), &
         emiss_anthrop(kbdim), emiss_bb(kbdim), emiss_wet(kbdim) ! eth_as_ch4
    REAL(dp) :: aircr(kbdim,k0_aircr:k1_aircr)
    REAL(dp) :: lightn(kbdim,klev)


    ! Executable statements:

    ! 1. Preparations:

    ! Upper boundary index:
    klev0 = klev-nlev_lowbound_bcond_fluxes+1

    ! Allocate memory:
    ALLOCATE(airmolec_layer(kproma,klev-nlev_lowbound_bcond_fluxes+1:klev))

    ! Calculate airmolec_col (total number of air molecules in a column of the  
    ! lowermost layers [molec/m^2]), airmolec_layer (number of air 
    ! molecules in individual layer [molec/m^2]) and dz_col (total thickness 
    ! of lowermost layers [m]):  
    airmolec_col(:) = 0._dp
    dz_col(:) = 0._dp
    DO jk = klev0, klev
       dz_layer(:) = (zlevb(1:kproma,jk)-zlevb(1:kproma,jk+1))*1000._dp 
                              ! thickness [m]
       airmolec_layer(:,jk) = dz_layer(:)*dens(1:kproma,jk)*1.E06_dp
                              ! *1.E06: [1/cm^3] -> [1/m^3]
       airmolec_col(:) = airmolec_col(:) + airmolec_layer(:,jk)
       dz_col(:) = dz_col(:) + dz_layer(:)
    ENDDO


    ! 2. CO:

    ! 2.1 CO surface emissions:

    fac = 1.E-03_dp/(atomweight*amco) ! [g(CO)/m^2/s] -> [CO molec/m^2/s]
    emiss(:) = coemissfac*co_emiss_surf(1:kproma)*fac/airmolec_col(:) ! eth_as_co


    ! CCMI synthetic tracers:
    if (lsynth) CALL synth_trac_flux(krow, kproma, kbdim, klev, ktrac, pxtm1, pxtte, papp1, ktrpwmo, airmolec_col)

    ! eth_as_nmvoc+
    IF (lnmvoc) THEN
       IF (lnmvoc_accmip) THEN

          fac = 1._dp/(atomweight*amco) ! [kg/m^2/s] -> [molec/m^2/s]
          ! 2.1.1 add contribution from NMVOC emissions to CO emissions
          ! Emissions are NMVOC emissions, scaling factors transfer NMVOC to right amount of CO
          emiss(:) = emiss(:) &
               + (0.83_dp * nmvoc_emiss_biogen(1:kproma) & 
               + 0.31_dp * nmvoc_emiss_bb(1:kproma) & 
               + 1._dp * nmvoc_emiss_anthrop(1:kproma)) & 
               * fac/airmolec_col(:)      

!!$          ! write CO contribution from NMVOC to output stream, kg(CO)/m2/s
!!$          do jl = 1, kproma
!!$             co_nmvocbio_d(jl,krow) = 0.83_dp * nmvoc_emiss_biogen(jl)
!!$             co_nmvocbb_d(jl,krow) =  0.31_dp * nmvoc_emiss_bb(jl)
!!$             co_nmvocanthrop_d(jl,krow) = 1._dp * nmvoc_emiss_anthrop(jl) 
!!$          end do
  
       ELSE ! old GEIA and GFED emissions, AS 29.05.2013

          fac = 1.E-03_dp/(atomweight*amco) ! [g/m^2/s] -> [molec/m^2/s]
          ! 2.1.1 add contribution from NMVOC emissions to CO emissions
          emiss(:) = emiss(:) &
               + ( 0.16_dp * nmvoc_emiss_biogen(1:kproma) * amco/12.011_dp & ! g(C) -> g(CO)
               + 0.43_dp * nmvoc_emiss_bb(1:kproma) & ! emissions already g(CO)
               + 0.43_dp * nmvoc_emiss_anthrop(1:kproma) * amco/12.011_dp) & ! g(C) -> g(CO)
               * fac/airmolec_col(:)
          
!!$          ! write CO contribution from NMVOC to output stream, kg(CO)/m2/s
!!$          do jl = 1, kproma
!!$             co_nmvocbio_d(jl,krow) = 1.e-3 * 0.16_dp * nmvoc_emiss_biogen(jl) * amco/12.011_dp
!!$             co_nmvocbb_d(jl,krow) =  1.e-3 * 0.43_dp * nmvoc_emiss_bb(jl)
!!$             co_nmvocanthrop_d(jl,krow) = 1.e-3 * 0.43_dp * nmvoc_emiss_anthrop(jl) * amco/12.011_dp
!!$          end do
       END IF
    END IF
    ! eth_as_nmvoc-

    ! 2.2 CO deposition:

    ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
    DO jl = 1, kproma   !FORALL (jl = 1:kproma) 
       deposvel(jl) = slm(jl,krow)*wcodep(1) + (1._dp-slm(jl,krow))*wcodep(2)
    END DO

    DO jk = klev0, klev
       ! Number of CO molecules in the lowermost layers [CO molec/m^2]:
       gasmolec_layer(:) = airmolec_layer(:,jk) * &
            (pxtm1(1:kproma,jk,idt_co)+time_step_len*pxtte(1:kproma,jk,idt_co))

       ! Calculate the ratio of CO molecules removed by deposition in 
       ! a column of the lowermost layers. The CO molecules are removed from 
       ! the individual layers proportional to the number of CO molecules in 
       ! the layer. 
       ! (1/airmolec_layer: [CO molec/m^2/s] -> [mixing ratio/m^2/s]): 
       depos(:) = deposvel(:)/dz_col(:)*gasmolec_layer(:)/ &
            airmolec_layer(:,jk)

       ! 2.3 Total tendency:

       pxtte(1:kproma,jk,idt_co) = &
            pxtte(1:kproma,jk,idt_co) + emiss(:) - depos(:)
    END DO


    ! 3. NO and NO2:

    ! 3.1 NOx emissions from aircrafts and NOx production from lightning:
    
    ! 3.1.1 Molecular weight of NOx:
    amnox = nofracnox_bcond*amno + no2fracnox_bcond*amno2

    ! 3.1.2 NOx emissions from aircrafts [NOx mixing ratio/s]:

    ! [NOx molec/cm^3/s] -> [NOx mixing ratio/s]:
    aircr(:,k0_aircr:k1_aircr) = nox_emiss_aircr(1:kproma,k0_aircr:k1_aircr) / &
         dens(1:kproma,k0_aircr:k1_aircr)
    
    ! 3.1.3 NOx production from lightning [NOx mixing ratio/s]:
    
    ! (a) NOx production from lightning [kg(NOx)/kg(air)/s] in layer k:
    !     ([kg(NOx)/m^2/s] in layer jk) / ([kg(air)/m^2] in layer jk) 

    do i=1,16
       zplt(i)=(i-1)*1._dp  ! z-axis of lightning profile (km), 1: surface -> 16: top
    enddo

    do jkp  = 1, kproma
       ! Profiles interpolation
       do jk=1,klev
          zkmh(klev+1-jk)=(zlevb(jkp,jk)+zlevb(jkp,jk+1))/2._dp ! altitude of full levels (km), 1: surface -> klev: top
       enddo

       profile_type = 0

       if(abs(philat_2d(jkp,krow)).le.40._dp .and. abs(philat_2d(jkp,krow)).gt.20._dp) profile_type=1
       if(abs(philat_2d(jkp,krow)).gt.40._dp) profile_type=2
       if(abs(philat_2d(jkp,krow)).le.20._dp .and.(slm(jkp,krow).ge.0.1_dp)) profile_type=3
       if(abs(philat_2d(jkp,krow)).le.20._dp .and.(slm(jkp,krow).lt.0.1_dp)) profile_type=4

       call lnt1(plt(:,profile_type),zplt,16,pltm,zkmh,klev)

       do jk=1,klev
          pltm(klev+1-jk)=pltm(klev+1-jk)*(zlevb(jkp,jk)-zlevb(jkp,jk+1)) ! interpolated profile, weighted by layer thickness, 1: surface -> klev: top
       enddo

       pltm(:)=pltm(:)/sum(pltm) ! fract. instead of %

       IF (.NOT.llnox_calc .or. .NOT.interactivelnox) THEN
          DO jk = 1, klev
             IF(pltm(klev+1-jk) .gt. 1.e-10_dp) then      
                lightn(jkp,jk) = & 
!!$                  nox_lightn_col(jkp)*pltm(jk)*1.E-03_dp / &
!!$                ((ph(jkp,jk+1)-ph(jkp,jk))*g)
                  nox_lightn_col(jkp)*pltm(klev+1-jk)*1.E-03_dp / &
                  ((ph(jkp,jk+1)-ph(jkp,jk))/g) ! bug fix dp*g -> dp/g

              ! additional output kg/m2/s
              elnox(jkp,jk,krow) = nox_lightn_col(jkp)*pltm(klev+1-jk)*1.E-03_dp

            ELSE
              lightn(jkp,jk) = 0._dp
              elnox(jkp,jk,krow) = 0._dp  
            END IF
          ENDDO
       ENDIF
       
       IF (llnox_calc .and. interactivelnox) then
          DO jk = 1, klev
            IF(pltm(klev+1-jk) .gt. 1.e-10_dp) THEN
              lightn(jkp,jk) = &
!!$                  lnox(jkp)*pltm(klev+1-jk)*1.E-03_dp / &
!!$                  ((ph(jkp,jk+1)-ph(jkp,jk))/g) ! bug fix dp*g -> dp/g
                  lnox_d(jkp,krow)*pltm(klev+1-jk)*1.E-03_dp / &
                  ((ph(jkp,jk+1)-ph(jkp,jk))/g) ! bug fix dp*g -> dp/g

              ! additional output kg/m2/s
              elnox(jkp,jk,krow) = lnox_d(jkp,krow)*pltm(klev+1-jk)*1.E-03_dp

            ELSE
              lightn(jkp,jk) = 0._dp
              elnox(jkp,jk,krow) = 0._dp           
            END IF  
          ENDDO
       ENDIF

    ENDDO ! jkp

    ! (b) [kg(NOx)/kg(air)/s] -> [NOx mixing ratio/s]
    lightn(:,1:klev) = lightn(:,1:klev)*amd/amnox
    
    ! 3.2 NO:

    ! 3.2.1 NO emissions from aircrafts:
    pxtte(1:kproma,k0_aircr:k1_aircr,idt_no) = &
         pxtte(1:kproma,k0_aircr:k1_aircr,idt_no) + &
         nofracnox_bcond*aircr(:,k0_aircr:k1_aircr)

    ! 3.2.2 NO production from lightning:
    pxtte(1:kproma,1:klev,idt_no) = &
         pxtte(1:kproma,1:klev,idt_no) + &
         nofracnox_bcond*lightn(:,1:klev)

    ! 3.2.3 NO surface emissions:

    fac = 1.E-03_dp/(atomweight*amno) ! [g(NO)/m^2/s] -> [NO molec/m^2/s]
    emiss(:) = nofracnox_bcond*nox_emiss_surf(1:kproma)*fac/airmolec_col(:)

    ! 3.2.4 NO deposition:

    ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
    DO jl = 1, kproma   !FORALL (jl = 1:kproma)
       deposvel(jl) = slm(jl,krow)*wnodep(1) + (1._dp-slm(jl,krow))*wnodep(2)
    END DO

    DO jk = klev0, klev
       ! Number of NO molecules in the lowermost layers [NO molec/m^2]:
       gasmolec_layer(:) = airmolec_layer(:,jk) * &
            (pxtm1(1:kproma,jk,idt_no)+time_step_len*pxtte(1:kproma,jk,idt_no))

       ! Calculate the ratio of NO molecules removed by deposition in 
       ! a column of the lowermost layers. The NO molecules are removed from 
       ! the individual layers proportional to the number of NO molecules in 
       ! the layer. 
       ! (1/airmolec_layer: [NO molec/m^2/s] -> [mixing ratio/m^2/s]): 
       depos(:) = deposvel(:)/dz_col(:)*gasmolec_layer(:)/ &
            airmolec_layer(:,jk)

       ! 3.2.5 Total tendency:

       pxtte(1:kproma,jk,idt_no) = &
            pxtte(1:kproma,jk,idt_no) + emiss(:) - depos(:)
    END DO

    ! 3.3 NO2:

    ! 3.3.1 NO2 emissions from aircrafts:
    pxtte(1:kproma,k0_aircr:k1_aircr,idt_no2) = &
         pxtte(1:kproma,k0_aircr:k1_aircr,idt_no2) + &
         no2fracnox_bcond*aircr(:,k0_aircr:k1_aircr)

    ! 3.3.2 NO2 production from lightning:
    pxtte(1:kproma,1:klev,idt_no2) = &
         pxtte(1:kproma,1:klev,idt_no2) + &
         no2fracnox_bcond*lightn(:,1:klev)

    ! 3.3.3 NO2 surface emissions:

    fac = 1.E-03_dp/(atomweight*amno2) ! [g(NO2)/m^2/s] -> [NO2 molec/m^2/s]
    emiss(:) = no2fracnox_bcond*nox_emiss_surf(1:kproma)*fac/airmolec_col(:)

    ! 3.3.4 NO2 deposition:

    ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
    DO jl = 1, kproma   !FORALL (jl = 1:kproma) 
       deposvel(jl) = slm(jl,krow)*wno2dep(1) + (1._dp-slm(jl,krow))*wno2dep(2)
    END DO

    DO jk = klev0, klev
       ! Number of NO2 molecules in the lowermost layers [NO2 molec/m^2]:
       gasmolec_layer(:) = airmolec_layer(:,jk) * &
            (pxtm1(1:kproma,jk,idt_no2)+&
            time_step_len*pxtte(1:kproma,jk,idt_no2))

       ! Calculate the ratio of NO2 molecules removed by deposition in 
       ! a column of the lowermost layers. The NO2 molecules are removed from 
       ! the individual layers proportional to the number of NO2 molecules in 
       ! the layer. 
       ! (1/airmolec_layer: [NO2 molec/m^2/s] -> [mixing ratio/m^2/s]): 
       depos(:) = deposvel(:)/dz_col(:)*gasmolec_layer(:)/ &
            airmolec_layer(:,jk)

       ! 3.3.5 Total tendency:

       pxtte(1:kproma,jk,idt_no2) = &
            pxtte(1:kproma,jk,idt_no2) + emiss(:) - depos(:)
    END DO


    ! 4. HNO3:

    ! 4.1 HNO3 deposition:

    ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
    DO jl = 1, kproma   !FORALL (jl = 1:kproma) 
       deposvel(jl) = &
            slm(jl,krow)*whno3dep(1) + (1._dp-slm(jl,krow))*whno3dep(2)
    END DO

    DO jk = klev0, klev
       ! Number of HNO3 molecules in the lowermost layers [HNO3 molec/m^2]:
       gasmolec_layer(:) = airmolec_layer(:,jk) * &
            (pxtm1(1:kproma,jk,idt_hno3)+ &
            time_step_len*pxtte(1:kproma,jk,idt_hno3))

       ! Calculate the ratio of HNO3 molecules removed by deposition in 
       ! a column of the lowermost layers. The HNO3 molecules are removed from 
       ! the individual layers proportional to the number of HNO3 molecules in 
       ! the layer. 
       ! (1/airmolec_layer: [HNO3 molec/m^2/s] -> [mixing ratio/m^2/s]): 
       depos(:) = deposvel(:)/dz_col(:)*gasmolec_layer(:)/ &
            airmolec_layer(:,jk)

       ! 4.2 Tendency:

       pxtte(1:kproma,jk,idt_hno3) = pxtte(1:kproma,jk,idt_hno3) - depos(:)
    END DO


    ! 5. O3:

    ! 5.1 O3 deposition:

    ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
    DO jl = 1, kproma   !FORALL (jl = 1:kproma) 
       deposvel(jl) = slm(jl,krow)*wo3dep(1) + (1._dp-slm(jl,krow))*wo3dep(2)
    END DO

    DO jk = klev0, klev
       ! Number of O3 molecules in the lowermost layers [O3 molec/m^2]:
       gasmolec_layer(:) = airmolec_layer(:,jk) * &
            (pxtm1(1:kproma,jk,idt_o3)+ &
            time_step_len*pxtte(1:kproma,jk,idt_o3))

       ! Calculate the ratio of O3 molecules removed by deposition in 
       ! a column of the lowermost layers. The O3 molecules are removed from 
       ! the individual layers proportional to the number of O3 molecules in 
       ! the layer. 
       ! (1/airmolec_layer: [O3 molec/m^2/s] -> [mixing ratio/m^2/s]): 
       depos(:) = deposvel(:)/dz_col(:)*gasmolec_layer(:)/ &
            airmolec_layer(:,jk)

       ! 5.2 Tendency:

       pxtte(1:kproma,jk,idt_o3) = pxtte(1:kproma,jk,idt_o3) - depos(:)
    END DO


    ! 6. H2O2:

    ! 6.1 H2O2 deposition:

    ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
    DO jl = 1, kproma   !FORALL (jl = 1:kproma) 
       deposvel(jl) = &
            slm(jl,krow)*wh2o2dep(1) + (1._dp-slm(jl,krow))*wh2o2dep(2)
    END DO

    DO jk = klev0, klev
       ! Number of H2O2 molecules in the lowermost layers [H2O2 molec/m^2]:
       gasmolec_layer(:) = airmolec_layer(:,jk) * &
            (pxtm1(1:kproma,jk,idt_h2o2)+ &
            time_step_len*pxtte(1:kproma,jk,idt_h2o2))

       ! Calculate the ratio of H2O2 molecules removed by deposition in 
       ! a column of the lowermost layers. The H2O2 molecules are removed from 
       ! the individual layers proportional to the number of H2O2 molecules in 
       ! the layer. 
       ! (1/airmolec_layer: [H2O2 molec/m^2/s] -> [mixing ratio/m^2/s]): 
       depos(:) = deposvel(:)/dz_col(:)*gasmolec_layer(:)/ &
            airmolec_layer(:,jk)

       ! 6.2 Tendency:

       pxtte(1:kproma,jk,idt_h2o2) = pxtte(1:kproma,jk,idt_h2o2) - depos(:)
    END DO

    ! eth_as_ch4

    IF ( lch4_flux) THEN
       ! 7. CH4:

       ! 7.1 CH4 surface emissions:
 
       scale_ch4 = ch4emissfac         
       fac = scale_ch4 * 1._dp/(atomweight*amch4) ! [kg(CH4)/m^2/s] -> [CH4 molec/m^2/s]
                                                ! Factor 1.31, 20110525
                                                ! Factor 1.32, 20110603 

       IF ( lch4_ipcc) THEN
          IF ( lch4_wetland) THEN
             emiss(:) = (ch4_emiss_ship(1:kproma) + ch4_emiss_grassfire(1:kproma) + &
                  ch4_emiss_forestfire(1:kproma) + ch4_emiss_agr(1:kproma) + &
                  ch4_emiss_awb(1:kproma) + ch4_emiss_dom(1:kproma) + &
                  ch4_emiss_ene(1:kproma) + ch4_emiss_ind(1:kproma) + &
                  ch4_emiss_tra(1:kproma) + ch4_emiss_wst(1:kproma) + &
                  ch4_emiss_slv(1:kproma) + ch4_emiss_wet(1:kproma))*fac/airmolec_col(:)
          ELSE
             emiss(:) = (ch4_emiss_ship(1:kproma) + ch4_emiss_grassfire(1:kproma) + &
                  ch4_emiss_forestfire(1:kproma) + ch4_emiss_agr(1:kproma) + &
                  ch4_emiss_awb(1:kproma) + ch4_emiss_dom(1:kproma) + &
                  ch4_emiss_ene(1:kproma) + ch4_emiss_ind(1:kproma) + &
                  ch4_emiss_tra(1:kproma) + ch4_emiss_wst(1:kproma) + &
                  ch4_emiss_slv(1:kproma))*fac/airmolec_col(:) 
          END IF
       ELSE
          emiss(:) = (ch4_emiss_wetlands(1:kproma) + ch4_emiss_bb(1:kproma) + &
                ch4_emiss_ene(1:kproma) + ch4_emiss_ind(1:kproma) + &
                ch4_emiss_avi(1:kproma) + ch4_emiss_air(1:kproma) + &
                ch4_emiss_intship(1:kproma) + ch4_emiss_ship(1:kproma) + &
                ch4_emiss_trans(1:kproma) + ch4_emiss_resi(1:kproma) + &              
                ch4_emiss_fuels(1:kproma) + ch4_emiss_agri(1:kproma) + &      
                ch4_emiss_waste(1:kproma) + ch4_emiss_wstwat(1:kproma) + &    
                ch4_emiss_other(1:kproma))*fac/airmolec_col(:)  
       END IF
          

       ! 7.2 CH4 deposition:
       ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
       ! no deposition of CH4 (Hauglustaine et al., 1994)
       
       ! 7.3 Total tendency:
       
       DO jk = klev0, klev
          
          pxtte(1:kproma,jk,idt_ch4) = &
               pxtte(1:kproma,jk,idt_ch4) + emiss(:)
       END DO


       ! 7.4 Isotopic structure of methane emissions
       
       IF ( lch4_isotope ) THEN

          fac = scale_ch4 * 1._dp/(atomweight*amch4) ! [kg(CH4)/m^2/s] -> [CH4 molec/m^2/s]

          ! anthropogenic emissions
          IF (lch4_ipcc) THEN
             emiss_anthrop(:) = (ch4_emiss_ship(1:kproma) + ch4_emiss_agr(1:kproma) + &
                  ch4_emiss_awb(1:kproma) +  ch4_emiss_dom(1:kproma) + &
                  ch4_emiss_ene(1:kproma) +  ch4_emiss_ind(1:kproma) + &
                  ch4_emiss_tra(1:kproma) +  ch4_emiss_wst(1:kproma) + &
                  ch4_emiss_slv(1:kproma))*fac/airmolec_col(:)
          ELSE
             emiss_anthrop(:) = (ch4_emiss_ene(1:kproma) + ch4_emiss_ind(1:kproma) + &
                  ch4_emiss_avi(1:kproma) + ch4_emiss_air(1:kproma) + &
                  ch4_emiss_intship(1:kproma) + ch4_emiss_ship(1:kproma) + &
                  ch4_emiss_trans(1:kproma) + ch4_emiss_resi(1:kproma) + &              
                  ch4_emiss_fuels(1:kproma) + ch4_emiss_agri(1:kproma) + &      
                  ch4_emiss_waste(1:kproma) + ch4_emiss_wstwat(1:kproma) + &    
                  ch4_emiss_other(1:kproma))*fac/airmolec_col(:)
          END IF

          zanthrop = (delta13C_anthrop/1000._dp + 1._dp) * r_ref

          ! biomass burning emissions
          IF (lch4_ipcc) THEN
             emiss_bb(:) = (ch4_emiss_grassfire(1:kproma) + ch4_emiss_forestfire(1:kproma)) * &
                  fac/airmolec_col(:)
          ELSE
             emiss_bb(:) = ch4_emiss_bb(1:kproma) * fac/airmolec_col(:)
          END IF

          zbb = (delta13C_wildfire/1000._dp + 1._dp) * r_ref

          ! prescribed wetland emissions
          
          IF (.not. lch4_ipcc) THEN
             emiss_wet(:) = ch4_emiss_wetlands(1:kproma)*fac/airmolec_col(:)
             zwet = (delta13C_wetlands/1000._dp + 1._dp) * r_ref
          END IF

          IF (lch4_ipcc) THEN
             ! wetland emissions
             IF ( lch4_wetland ) THEN
                
                emiss_wet(:) = ch4_emiss_wet(1:kproma)*fac/airmolec_col(:)
                
                zwet = (delta13C_wetlands/1000._dp + 1._dp) * r_ref
                
                DO jk = klev0, klev
                   
                   pxtte(1:kproma,jk,idt_12c) = pxtte(1:kproma,jk,idt_12c) + &
                        emiss_anthrop(:)/(1._dp + zanthrop) + emiss_bb(:)/(1._dp + zbb) + &
                        emiss_wet(:)/(1._dp + zwet)
                   
                   pxtte(1:kproma,jk,idt_13c) = pxtte(1:kproma,jk,idt_13c) + &
                        zanthrop*emiss_anthrop(:)/(1._dp + zanthrop) + zbb*emiss_bb(:)/(1._dp + zbb) + &
                        zwet*emiss_wet(:)/(1._dp + zwet)
                   
                END DO
                
             ELSE
                
                DO jk = klev0, klev
                   
                   pxtte(1:kproma,jk,idt_12c) = pxtte(1:kproma,jk,idt_12c) + &
                        emiss_anthrop(:)/(1._dp + zanthrop) + emiss_bb(:)/(1._dp + zbb)
                   
                   pxtte(1:kproma,jk,idt_13c) = pxtte(1:kproma,jk,idt_13c) + &
                        zanthrop*emiss_anthrop(:)/(1._dp + zanthrop) + zbb*emiss_bb(:)/(1._dp + zbb)
                   
!!$                DO jl = 1, kproma
!!$                   IF ( (pxtte(jl,jk,idt_12c) + pxtte(jl,jk,idt_13c)) &
!!$                        .ne. pxtte(jl,jk,idt_ch4)) THEN
!!$                      WRITE(message_text,*) &
!!$                           'CH4 isotopes: 12C-(CH4) + 13C-(CH4) /= CH4!', pxtte(jl,jk,idt_12c), pxtte(jl,jk,idt_13c), pxtte(jl,jk,idt_ch4)
!!$                      CALL message('',TRIM(message_text))
!!$                      CALL finish('socol_boundary_conditions','Run terminated')
!!$                   ENDIF
!!$                END DO
                   
                END DO

             END IF

          ELSE

             DO jk = klev0, klev
                   
                pxtte(1:kproma,jk,idt_12c) = pxtte(1:kproma,jk,idt_12c) + &
                     emiss_anthrop(:)/(1._dp + zanthrop) + emiss_bb(:)/(1._dp + zbb) + &
                     emiss_wet(:)/(1._dp + zwet)
                   
                pxtte(1:kproma,jk,idt_13c) = pxtte(1:kproma,jk,idt_13c) + &
                     zanthrop*emiss_anthrop(:)/(1._dp + zanthrop) + zbb*emiss_bb(:)/(1._dp + zbb) + &
                     zwet*emiss_wet(:)/(1._dp + zwet)
                   
             END DO

          END IF
       
       END IF
 
    END IF

    ! eth_as_20100527+

    ! 8. H2:

    ! 8.1 H2 deposition:

    ! Deposition velocity (over land (1) and/or sea (2)) [m/s]:
    DO jl = 1, kproma   !FORALL (jl = 1:kproma)
       deposvel(jl) = &
            slm(jl,krow)*wh2dep(1) + (1._dp-slm(jl,krow))*wh2dep(2)
    END DO

    DO jk = klev0, klev
       ! Number of H2 molecules in the lowermost layers [H2 molec/m^2]:
       gasmolec_layer(:) = airmolec_layer(:,jk) * &
            (pxtm1(1:kproma,jk,idt_h2)+ &
            time_step_len*pxtte(1:kproma,jk,idt_h2))

       ! Calculate the ratio of H2 molecules removed by deposition in
       ! a column of the lowermost layers. The H2 molecules are removed from
       ! the individual layers proportional to the number of H2 molecules in
       ! the layer.
       ! (1/airmolec_layer: [H2 molec/m^2/s] -> [mixing ratio/m^2/s]):
       depos(:) = deposvel(:)/dz_col(:)*gasmolec_layer(:)/ &
            airmolec_layer(:,jk)

       ! 8.2 Tendency:

       pxtte(1:kproma,jk,idt_h2) = pxtte(1:kproma,jk,idt_h2) - depos(:)
    END DO

    ! eth_as_20100527-
!!$
    ! eth_as_tropchem+

    ! 9. Species for isoprene chemistry

    ! 9.1 C5H8 surface emissions:

    fac = 1._dp/(atomweight*amc5h8) ! [kg/m^2/s] -> [molec/m^2/s]
    emiss(:) = c5h8_emiss(1:kproma)*fac/airmolec_col(:) ! * 2._dp ! * 1.2_dp ! eth_as_20110512
    
    DO jk = klev0, klev
       pxtte(1:kproma,jk,idt_c5h8) = &
            pxtte(1:kproma,jk,idt_c5h8) + emiss(:)
    END DO
    
    ! 9.2 HCHO surface emissions:

    fac = 1._dp/(atomweight*amch2o) ! [kg/m^2/s] -> [molec/m^2/s]
    emiss(:) = ch2o_emiss(1:kproma)*fac/airmolec_col(:)
    
    DO jk = klev0, klev
       pxtte(1:kproma,jk,idt_ch2o) = &
            pxtte(1:kproma,jk,idt_ch2o) + emiss(:)
    END DO    

    ! 9.3 HCOOH surface emissions:

    fac = 1._dp/(atomweight*amhcooh) ! [kg/m^2/s] -> [molec/m^2/s]
    !!$ emiss(:) = hcooh_emiss(1:kproma)*fac/airmolec_col(:)
    emiss(:) = 0._dp    

    DO jk = klev0, klev
       pxtte(1:kproma,jk,idt_hcooh) = &
            pxtte(1:kproma,jk,idt_hcooh) + emiss(:)
    END DO    

    ! 9.4 CH3COOH surface emissions:

    fac = 1._dp/(atomweight*amch3cooh) ! [kg/m^2/s] -> [molec/m^2/s]
    emiss(:) = ch3cooh_emiss(1:kproma)*fac/airmolec_col(:)
    
    DO jk = klev0, klev
       pxtte(1:kproma,jk,idt_ch3cooh) = &
            pxtte(1:kproma,jk,idt_ch3cooh) + emiss(:)
    END DO  

    ! eth_as_tropchem-

    ! 10. Deallocate memory:
    DEALLOCATE(airmolec_layer)

  END SUBROUTINE set_socol_bcond_fluxes
!===============================================================================
  SUBROUTINE cleanup_socol_bcond

    ! Calls subroutines to deallocate memory from boundary condition fields 
    ! (allocatable module variables).

    ! *cleanup_socol_bcond* is called from *call_free_submodel_memory*, 
    ! src/call_submodels.f90.
    !
    ! M. Schraner, ETH Zurich, January 2009

    USE mo_radiation,             ONLY: iaero
    USE mo_socol_co_nox,          ONLY: cleanup_socol_co_nox
    USE mo_socol_ghg_ods,         ONLY: cleanup_socol_ghg_ods
    USE mo_socol_nmvoc,           ONLY: cleanup_socol_nmvoc      ! eth_as_nmvoc
    USE mo_socol_ch4,             ONLY: cleanup_socol_ch4        ! eth_as_ch4
    USE mo_socol_namelist
    USE mo_socol_qbo,             ONLY: cleanup_socol_qbo
    USE mo_socol_strataerosols,   ONLY: cleanup_socol_strataerosols
    USE mo_socol_tropoaerosols,   ONLY: cleanup_socol_tropoaerosols
    USE mo_socol_photo,           ONLY: cleanup_socol_photo      ! eth_as_tropchem

    IMPLICIT NONE
    
    ! Greenhouse gases and ODSs:
    CALL cleanup_socol_ghg_ods

    ! Stratospheric aerosols:
    CALL cleanup_socol_strataerosols

    ! Tropospheric aerosols:
    IF (iaero .EQ. 10) CALL cleanup_socol_tropoaerosols

    ! CO and NOx emissions:
    IF (lchem) CALL cleanup_socol_co_nox

    ! eth_as_nmvoc
    ! NMVOC emissions:
    IF (lchem) CALL cleanup_socol_nmvoc

    ! eth_as_ch4
    ! CH4 emissions:
    IF (lch4_flux) CALL cleanup_socol_ch4

    ! QBO nudging:
    IF (lqbonudg) CALL cleanup_socol_qbo

    ! eth_as_tropchem
!!$    ! Photolysis rates
!!$    call cleanup_socol_photo 

  END SUBROUTINE cleanup_socol_bcond
!===============================================================================

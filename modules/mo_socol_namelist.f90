MODULE mo_socol_namelist 

  ! Control variables for SOCOL (namelist *socolctl*)
  ! Contains subroutine to read namelist
  !
  ! Martin Schraner, ETH Zurich, October 2008

  USE mo_kind,               ONLY: dp

  IMPLICIT NONE

  LOGICAL :: lchem               = .TRUE.  ! .TRUE. for calling chemistry
                                           ! module MEZON; 
                                           ! .FALSE. for running only GCM
  REAL(dp):: delta_time_chem     = 7200._dp ! time interval between two calls of
                                           ! chemistry (MEZON) [s]
  LOGICAL :: lfamily_correct_adv = .TRUE.  ! family correction of Cly, Bry,
                                           ! and NOy after advection
  LOGICAL :: lco2_var_rad        = .TRUE.  ! concentrations changing in time
                                           ! (radiation)
  LOGICAL :: lco2_var_chem       = .TRUE.  ! concentrations changing in time
                                           ! (chemistry)
  LOGICAL :: lch4_var_rad        = .TRUE.  ! concentrations changing in time
                                           ! (radiation)
  LOGICAL :: lch4_var_chem       = .TRUE.  ! concentrations changing in time
                                           ! (chemistry)
  LOGICAL :: lch4_coupl          = .TRUE.  ! coupling of radiation and 
                                           ! chemistry module
  LOGICAL :: lch4_nocoupl_3d     = .FALSE. ! no coupling, but precribed 3d-data
                                           ! files in radiation module from 
                                           ! former coupled simulation
  LOGICAL :: lch4_flux           = .FALSE. ! methane fluxes instead of
                                           ! prescribed concentrations  ! eth_as_ch4
  LOGICAL :: lsurfemch4_var_chem = .FALSE.  ! CH4 surface emissions changing in
                                           ! time (chemistry)           ! eth_as_ch4
  LOGICAL :: lch4_wetland        = .FALSE. ! CH4 emissions from wetlands ! eth_as_ch4
  LOGICAL :: lch4_isotope        = .FALSE. ! CH4 isotopic composition diagnostic ! eth_as_ch4
  LOGICAL :: lch4_ipcc           = .TRUE.  ! IPCC CH4 emission inventory ! eth_as_ch4
  LOGICAL :: ln2o_var_rad        = .TRUE.  ! concentrations changing in time
                                           ! (radiation)
  LOGICAL :: ln2o_var_chem       = .TRUE.  ! concentrations changing in time
                                           ! (chemistry)
  LOGICAL :: ln2o_coupl          = .TRUE.  ! coupling of radiation and
                                           ! chemistry module
  LOGICAL :: ln2o_nocoupl_3d     = .FALSE. ! no coupling, but prescribed 3d-data
                                           ! files in radiation module from 
                                           ! former coupled simulation
  LOGICAL :: lodscl_var_rad      = .TRUE.  ! concentrations changing in time
                                           ! (radiation)
  LOGICAL :: lodscl_var_chem     = .TRUE.  ! concentrations changing in time
                                           ! (chemistry)
  LOGICAL :: lodscl_coupl        = .TRUE.  ! coupling of radiation and
                                           ! chemistry module
  LOGICAL :: lodscl_nocoupl_3d   = .FALSE. ! no coupling, but precribed 3d-data
                                           ! files in radiation module from 
                                           ! former coupled simulation
  LOGICAL :: lodsbr_var_chem     = .TRUE.  ! concentrations changing in time
                                           ! (chemistry)
  LOGICAL :: lo3_coupl           = .TRUE.  ! coupling of radiation and
                                           ! chemistry module
  LOGICAL :: lh2o_coupl          = .TRUE.  ! coupling of radiation and
                                           ! chemistry module
  INTEGER :: nlevabtropo_h2ogcm  = 3       ! Tropopause level
                                           ! + nlevabtropo_seph2o
                                           ! = Pressure level, below which 
                                           !   GCM field is used for the 
                                           !   CTM. 
                                           ! !! Only necessary in case that
                                           !    lh2o_coupl=.FALSE. and
                                           !    lchem=.TRUE. !!
  LOGICAL :: laircrnox_var_chem  = .TRUE.  ! NOx aircraft emissions changing in
                                           ! time (chemistry)
  LOGICAL :: lsurfemnox_var_chem = .TRUE.  ! NOx surface emissions changing in
                                           ! time (chemistry)
  LOGICAL :: lsurfemco_var_chem  = .TRUE.  ! CO surface emissions changing in
                                           ! time (chemistry)
  LOGICAL :: lnmvoc              = .TRUE.  ! additional CO from NMVOC  ! eth_as_nmvoc

  LOGICAL :: lnmvoc_accmip       = .FALSE. ! ACCMIP NMHC emissions  !eth_as_nmvoc

  LOGICAL :: lsurfemnmvoc_var_chem  = .TRUE.  ! NMVOC surface emissions changing in
                                           ! time (chemistry) ! eth_as_nmvoc
  LOGICAL :: lsurfemc5h8_var_chem  = .TRUE.  ! C5H8 surface emissions changing in
                                           ! time (chemistry) ! eth_as
  LOGICAL :: lsurfemch2o_var_chem  = .TRUE.  ! CH2O surface emissions changing in
                                           ! time (chemistry) ! eth_as
  LOGICAL :: lsurfemch3cooh_var_chem  = .TRUE.  ! CH3COOH surface emissions changing in
                                           ! time (chemistry) ! eth_as
  LOGICAL :: lstrataer_var_rad   = .TRUE.  ! annual and monthly changing 
                                           ! concentrations of stratospheric 
                                           ! aerosols (background and 
                                           ! volcanoes) for radiation module
  LOGICAL :: lstrataer_var_chem  = .TRUE.  ! ditto for chemistry module
  LOGICAL :: lstrataer_bg_rad    = .TRUE.  ! stratospheric background aerosol 
                                           ! climatology for radiation module
  LOGICAL :: lstrataer_bg_chem   = .TRUE.  ! ditto for chemistry module
  LOGICAL :: ltropaer_var_rad   = .TRUE.  ! annual and monthly changing 
                                           ! concentrations of tropospheric 
                                           ! aerosols for radiation module
  LOGICAL :: lsolarvar           = .TRUE.  ! annual and monthly changing solar
                                           ! constant (radiation and chemistry)
  LOGICAL :: lsrb_lya_heating    = .TRUE.  ! Schumman-Runge and Lyman alpha 
                                           ! heating parameterization
  LOGICAL :: gcr                 = .TRUE. ! Activate galactic cosmic rays param
                                            ! from Usoskin et al.
  LOGICAL :: eep                 = .TRUE. ! Activate NOx influx as a function of Ap index
                                            ! (Baumgaertner et al., 2009)
  LOGICAL :: spe                 = .TRUE. ! Activate HOx and NOx production by SPEs
  LOGICAL :: deltaecorr          = .FALSE. ! Delta-E correction factor for photolysis
                                            ! of N2O, N2O5, X11 and X12
  LOGICAL :: lext_oh             = .FALSE. ! prescribe external (tropospheric) OH distribution for
                                           ! chemical calculations ! eth_as_ohclim
  LOGICAL :: lscav               = .FALSE. ! use scavenging module ! eth_as_scav
  LOGICAL :: lo3orig             = .FALSE.  ! use O3 transport diagnostic
  LOGICAL :: interactivelnox     = .FALSE. ! .true. for interactive lightning  !eth_ts
  LOGICAL :: lsynth              = .TRUE.  !   .true. for synthetic tracers   !CCMI tracers
  LOGICAL :: lnudg_2h            = .TRUE.  ! 2 hour output for CAO
  LOGICAL :: lheppa              = .TRUE. ! Heppa nox nudging

  REAL(dp):: sun_irrad_const(6)  = (/ 2.0494015_dp, 175.44220_dp, &
       448.35958_dp, 439.98852_dp, 249.60072_dp, 50.655947_dp/)
                                           ! Solar irradiance for 6 SW spectral
                                           ! bands if lsolarvar=.FALSE.
                                           ! (Default: mean over 1977-1998)
  REAL(dp):: sun_srb_const(12) = (/ &
       0.14128786_dp,  0.13687350_dp,  0.11948691_dp,  0.097687500_dp, &
       0.075445000_dp, 0.059813582_dp, 0.053836942_dp, 0.060293764_dp, &
       0.076832954_dp, 0.098558728_dp, 0.11871395_dp,  0.13763564_dp  /)
                                           ! Monthly values of Schumann-Runge 
                                           ! parameter if lsolarvar=.FALSE.
                                           ! (Default: mean over 1977-1998)
  REAL(dp):: sun_lya_const(12) = (/ &
       1.8168500_dp, 1.8146227_dp, 1.6960455_dp, 1.5721864_dp, &
       1.4345545_dp, 1.3454818_dp, 1.3084364_dp, 1.3585182_dp, &
       1.4669546_dp, 1.5910318_dp, 1.6758000_dp, 1.8349682_dp /)
                                           ! Monthly values of Lyman alpha
                                           ! parameter if lsolarvar=.FALSE.
                                           ! (Default: mean over 1977-1998)
  REAL(dp):: sun_nopho_const(12) = (/ &
       5.6252500E-06_dp, 5.6074546E-06_dp, 5.5348955E-06_dp, 5.4445590E-06_dp, &
       5.3521727E-06_dp, 5.2874045E-06_dp, 5.2625500E-06_dp, 5.2895319E-06_dp, &
       5.3583409E-06_dp, 5.4484682E-06_dp, 5.5315227E-06_dp, 5.6108273E-06_dp /)
                                           ! Monthly values of NO-photolysis 
                                           ! parameter if lsolarvar=.FALSE.
                                           ! (Default: mean over 1977-1998)
  LOGICAL :: lphotfull           = .FALSE. ! Full photolysis calculation
  LOGICAL :: lsphericdependphot  = .FALSE. ! Spherical dependency for photolysis
                                           ! lookup-table 
  LOGICAL :: lhetchem            = .TRUE.  ! Heterogeneous chemistry
  REAL(dp):: hetnat_lowlev       = 0._dp! Lower boundary for PSCs [hPa]
                                           ! 0.: level below tropopause
  REAL(dp):: hetnat_uplev        = 11._dp  ! Upper boundary for PSCs [hPa]
  REAL(dp):: hetnat_north        = 0._dp   ! Bound for PSCs on north hemis [deg]
  REAL(dp):: hetnat_south        = 0._dp   ! Bound for PSCs on south hemis [deg]
  REAL(dp):: hetice_lowlev       = 130._dp   ! Lower boundary for PSCs [hPa]
                                           ! 0.: level below tropopause
  REAL(dp):: hetice_uplev        = 11._dp  ! Upper boundary for PSCs [hPa]
  REAL(dp):: hetice_north        = 50._dp  ! Bound for PSCs on north hemis [deg]
  REAL(dp):: hetice_south        = -50._dp ! Bound for PSCs on south hemis [deg]
  REAL(dp):: hetnat_rmode        = 5.0E-04_dp ! Mode radius of NAT particles 
                                              ! [cm]
  REAL(dp):: hetice_ndens        = 1.0E-02_dp ! Number of ice particles [1/cm3]
  LOGICAL :: lqbonudg            = .FALSE.  ! QBO nudging
  REAL(dp):: tauqbonudg          = 7._dp   ! Damping time of QBO nuding [days]
                                           ! (7 days: medium damping)
  INTEGER :: cyear                         ! Boundary conditions are fixed 
                                           ! at year <cyear>
  REAL(dp):: co2fac              = 1._dp   ! CO2 factor (for ensembles)
  LOGICAL :: lpco2               = .FALSE. ! passive CO2 tracer for age of air
  REAL(dp):: vini_pco2           = 0._dp   ! initialization value for passive CO2 tracer
  REAL(dp):: coemissfac          = 1.5_dp  ! scaling factor for direct CO emissions ! eth_as_co
  REAL(dp):: ch4emissfac         = 1._dp   ! scaling factor for CH4 emissions ! eth_as_ch4
!===============================================================================
! switch for boundary condition initialisation (not part of the namelist):
  LOGICAL :: linit_socol_bcond   = .TRUE.      

CONTAINS

  SUBROUTINE init_socol

    ! Reads namelist *socolctl*
    ! *init_socol* is called from *call_init_submodels*, src/call_submodels.f90

    USE mo_doctor,             ONLY: nout, nerr
    USE mo_exception,          ONLY: finish
    USE mo_mpi,                ONLY: p_io, p_parallel, p_parallel_io, p_bcast
    USE mo_namelist,           ONLY: position_nml, nnml, POSITIONED, MISSING, &
                                      LENGTH_ERROR, READ_ERROR
    ! eth_as_ohclim+
    USE mo_ohclim         ,ONLY: su_ohclim                     
    ! eth_as_ohclim-

    IMPLICIT NONE
 
    ! Local scalars: 

    INTEGER           :: ierr  ! error return value from position_nml

    INCLUDE 'socolctl.inc'

    ! Executable statements:

    cyear = -9999

    IF (p_parallel_io) THEN
      CALL position_nml ('SOCOLCTL', status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
        READ (nnml, socolctl)
      CASE(MISSING)
         CALL finish ('init_socol', &
              'namelist socolctl is missing in namelist.echam')
      CASE(LENGTH_ERROR)
         CALL finish ('init_socol', &
              'namelist socolctl is longer than in the definition of socolctl.in')
      CASE(READ_ERROR)
         CALL finish ('init_socol', 'general read error when reading namelist')
      END SELECT
    ENDIF

    IF (p_parallel) THEN
      CALL p_bcast (lchem,              p_io)
      CALL p_bcast (delta_time_chem,    p_io)
      CALL p_bcast (lfamily_correct_adv,p_io)
      CALL p_bcast (lco2_var_rad,       p_io)
      CALL p_bcast (lco2_var_chem,      p_io)
      CALL p_bcast (lch4_var_rad,       p_io)
      CALL p_bcast (lch4_var_chem,      p_io)
      CALL p_bcast (lch4_coupl,         p_io)
      CALL p_bcast (lch4_nocoupl_3d,    p_io)
      CALL p_bcast (lch4_flux,          p_io)  ! eth_as_ch4
      CALL p_bcast (lsurfemch4_var_chem,p_io)  ! eth_as_ch4
      CALL p_bcast (lch4_wetland,       p_io)  ! eth_as_ch4
      CALL p_bcast (lch4_isotope,       p_io)  ! eth_as_ch4
      CALL p_bcast (lch4_ipcc,          p_io)  ! eth_as_ch4
      CALL p_bcast (ln2o_var_rad,       p_io)
      CALL p_bcast (ln2o_var_chem,      p_io)
      CALL p_bcast (ln2o_coupl,         p_io)
      CALL p_bcast (ln2o_nocoupl_3d,    p_io)
      CALL p_bcast (lodscl_var_rad,     p_io)
      CALL p_bcast (lodscl_var_chem,    p_io)
      CALL p_bcast (lodscl_coupl,       p_io)
      CALL p_bcast (lodscl_nocoupl_3d,  p_io)
      CALL p_bcast (lodsbr_var_chem,    p_io)
      CALL p_bcast (lo3_coupl,          p_io)
      CALL p_bcast (lh2o_coupl,         p_io)
      CALL p_bcast (nlevabtropo_h2ogcm, p_io)
      CALL p_bcast (laircrnox_var_chem, p_io)
      CALL p_bcast (lsurfemnox_var_chem,p_io)
      CALL p_bcast (lsurfemco_var_chem, p_io)
      CALL p_bcast (lnmvoc,             p_io) ! eth_as_nmvoc
      
      print*, 'test p_bcast1'
      CALL p_bcast (lnmvoc_accmip,      p_io) ! eth_as_nmvoc
      print*, 'test p_bcast2'
      
      CALL p_bcast (lsurfemnmvoc_var_chem, p_io) ! eth_as_nmvoc
      CALL p_bcast (lsurfemc5h8_var_chem, p_io) ! eth_as
      CALL p_bcast (lsurfemch2o_var_chem, p_io) ! eth_as
      CALL p_bcast (lsurfemch3cooh_var_chem, p_io) ! eth_as
      CALL p_bcast (lstrataer_var_rad,  p_io)
      CALL p_bcast (lstrataer_var_chem, p_io)
      CALL p_bcast (lstrataer_bg_rad,   p_io)
      CALL p_bcast (lstrataer_bg_chem,  p_io)
      CALL p_bcast (ltropaer_var_rad,   p_io)
      CALL p_bcast (lsolarvar,          p_io)
      CALL p_bcast (lsrb_lya_heating,   p_io)
      CALL p_bcast (deltaecorr,   p_io)
      CALL p_bcast (sun_irrad_const,    p_io)
      CALL p_bcast (sun_srb_const,      p_io)
      CALL p_bcast (sun_lya_const,      p_io)
      CALL p_bcast (sun_nopho_const,    p_io)
      CALL p_bcast (lphotfull,          p_io)
      CALL p_bcast (lsphericdependphot, p_io)
      CALL p_bcast (lhetchem,           p_io)
      CALL p_bcast (hetnat_lowlev,      p_io)
      CALL p_bcast (hetnat_uplev,       p_io)
      CALL p_bcast (hetnat_north,       p_io)
      CALL p_bcast (hetnat_south,       p_io)
      CALL p_bcast (hetice_lowlev,      p_io)
      CALL p_bcast (hetice_uplev,       p_io)
      CALL p_bcast (hetice_north,       p_io)
      CALL p_bcast (hetice_south,       p_io)
      CALL p_bcast (hetnat_rmode,       p_io)
      CALL p_bcast (hetice_ndens,       p_io)
      CALL p_bcast (lqbonudg,           p_io)
      CALL p_bcast (tauqbonudg,         p_io)
      CALL p_bcast (cyear,              p_io)
      CALL p_bcast (co2fac,             p_io)
      CALL p_bcast (lpco2,              p_io)
      CALL p_bcast (vini_pco2,          p_io)
      CALL p_bcast (gcr,                p_io)
      CALL p_bcast (eep,                p_io)
      CALL p_bcast (spe,                p_io)
      CALL p_bcast (coemissfac,         p_io) ! eth_as_co
      CALL p_bcast (ch4emissfac,        p_io) ! eth_as_ch4
      CALL p_bcast (lext_oh,            p_io) ! eth_as_ohclim
      CALL p_bcast (lscav,              p_io) ! eth_as_scav
      CALL p_bcast (lo3orig,            p_io) !O3 transport diagnostic
      CALL p_bcast (interactivelnox, p_io)    !eth_ts
      CALL p_bcast (lsynth, p_io)    !CCMI tracers
      CALL p_bcast (linit_socol_bcond, p_io) 
      CALL p_bcast (lnudg_2h,           p_io)
      CALL p_bcast (lheppa,             p_io)
    ENDIF

    IF (.NOT. lchem) THEN
       lch4_coupl   = .FALSE.
       ln2o_coupl   = .FALSE.
       lodscl_coupl = .FALSE.
       lo3_coupl    = .FALSE.
       lh2o_coupl   = .FALSE.
       lo3orig      = .FALSE.
    ENDIF

    IF ((lch4_var_rad .AND. .NOT. lch4_var_chem) .OR. &
         (.NOT. lch4_var_rad .AND. lch4_var_chem)) lch4_coupl       = .FALSE.
    IF ((ln2o_var_rad .AND. .NOT. ln2o_var_chem) .OR. &
         (.NOT. ln2o_var_rad .AND. ln2o_var_chem)) ln2o_coupl       = .FALSE.
    IF ((lodscl_var_rad .AND. .NOT. lodscl_var_chem) .OR. &
         (.NOT. lodscl_var_rad .AND. lodscl_var_chem)) lodscl_coupl = .FALSE.

    IF (lch4_coupl) lch4_nocoupl_3d           = .FALSE.
    IF (ln2o_coupl) ln2o_nocoupl_3d           = .FALSE.
    IF (lodscl_coupl) lodscl_nocoupl_3d       = .FALSE.

    IF (lstrataer_var_rad) lstrataer_bg_rad   = .FALSE. 
    IF (lstrataer_var_chem) lstrataer_bg_chem = .FALSE.

    ! eth_as_ch4+
    IF ( .not. lch4_flux) lch4_wetland = .FALSE.
    IF ( .not. lch4_flux) lch4_isotope = .FALSE.
    IF ( .not. lch4_flux) lch4_ipcc = .FALSE.  
    IF ( .not. lch4_flux) lsurfemch4_var_chem = .FALSE. 
    ! eth_as_ch4-

    ! eth_as_nmvoc+
    IF (.not. lnmvoc)  lsurfemnmvoc_var_chem = .FALSE.
    ! eth_as_nmvoc-

    IF ((cyear .EQ. -9999) .AND. &
         (.NOT. lco2_var_rad .OR. .NOT. lco2_var_chem .OR. &
         .NOT. lch4_var_rad .OR. .NOT. lch4_var_chem .OR. &
         .NOT. ln2o_var_rad .OR. .NOT. ln2o_var_chem .OR. &
         .NOT. lodscl_var_rad .OR. .NOT. lodscl_var_chem .OR. &
         .NOT. lodsbr_var_chem.OR. .NOT. laircrnox_var_chem .OR. &
         .NOT. lsurfemnox_var_chem .OR. .NOT. lsurfemco_var_chem .OR. &
         .NOT. lstrataer_var_rad .OR. .NOT. lstrataer_var_chem .OR. &
         .NOT. lsolarvar)) &
         CALL finish ('init_socol', &
         'cyear must be provided in namelist *SOCOLCTL*')

   IF (SIZE(sun_irrad_const) .NE. 6) &
        CALL finish ('init_socol', 'sun_irrad_const must have length 6')
   IF (SIZE(sun_srb_const) .NE. 12) &
        CALL finish ('init_socol', 'sun_srb_const must have length 12')
   IF (SIZE(sun_lya_const) .NE. 12) &
        CALL finish ('init_socol', 'sun_lya_const must have length 12')
   IF (SIZE(sun_nopho_const) .NE. 12) &
        CALL finish ('init_socol', 'sun_nopho_const must have length 12')
         
  END SUBROUTINE init_socol

END MODULE mo_socol_namelist

MODULE mo_socol_sun

  ! Description:

  ! Organises time depending solar data, i.e.
  ! a) irradiance of 6 SW spectral bands
  ! b) Parameters for Schumann-Runge bands, Lyman-alpha parameterization
  ! c) Parameter for NO-photolysis
  ! The module contains subroutines for reading data sets, interpolation of 
  ! to the current time step and parameterizations for additional heating by
  ! the Schumann-Runge bands and Lyman-alpha.

  ! M. Schraner, ETH ZÅ¸rich, January 2009

  USE mo_constants,          ONLY: amd, g, amo3, api
  USE mo_control,            ONLY: nlev
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_geoloc,             ONLY: amu0_x
  USE mo_kind,               ONLY: dp
  USE mo_mpi,                ONLY: p_io, p_parallel, p_parallel_io, p_bcast
  USE mo_radiation,          ONLY: solc
  USE mo_socol_constants,    ONLY: amo2, mro2, avoinv, cp 
  USE mo_socol_interpo,      ONLY: wgt1_chem, wgt2_chem, wgt1_rad, wgt2_rad, &
                                   yw1_chem, yw2_chem, yw1_rad, yw2_rad, &
                                   m3w1_chem, m3w2_chem
  USE mo_socol_namelist
  USE mo_socol_readfile,     ONLY: socol_ascii_yearmo_tab_y, socol_read_netcdf
  USE mo_socol_time_control, ONLY: l_trigchem
  USE mo_sw,                 ONLY: nsw, rsun          ! number of SW bands
  USE mo_time_control,       ONLY: l_trigrad
  USE mo_socol_tracers,      ONLY: idt_o3
  USE mo_memory_g1a,         ONLY: xtm1

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  ! Irradiance of SW spectral bands (monthly values):
  REAL(dp) :: sun_irrad_y(nsw,0:13)

  ! Schumann-Runge-bands parameter (monthly values):
  REAL(dp) :: sun_srb_y(0:13)

  ! Lyman alpha parameter (monthly values):
  REAL(dp) :: sun_lya_y(0:13)

  !Hartley parameter (monthly values):
  REAL(dp) :: sun_har_y(0:13)
    
  ! Huggins parameters (monthly values):
    REAL(dp) :: sun_hug1_y(0:13)
    REAL(dp) :: sun_hug2_y(0:13)
    
  ! NO-photolysis parameter (monthly values):
  REAL(dp) :: sun_nopho_y(0:13)

  ! Ditto, interpolated to current time step:
  REAL(dp), PUBLIC :: sun_nopho
  REAL(dp) :: sun_srb
  REAL(dp) :: sun_lya
  REAL(dp) :: sun_har
  REAL(dp) :: sun_hug1
  REAL(dp) :: sun_hug2

  ! Dimensions of photolysis rates:
  INTEGER, PARAMETER, PUBLIC :: n_photolytab_o2 = 40
  INTEGER, PARAMETER, PUBLIC :: n_photolytab_o3 = 44
  INTEGER, PARAMETER, PUBLIC :: nphotolysis = 54 ! number of different 
                                                  ! photolysis reactions

  ! Number of O2/O3 molecules passed by a solar beam [molec/cm^2] in
  ! look up table of the photolysis rates:
  REAL(dp), PUBLIC :: xo2(n_photolytab_o2), xo3(n_photolytab_o3)

  ! Names of photolysis reactions:
  CHARACTER(31) :: photolysis_name(nphotolysis)

  ! Photolysis rates of preceeding, current and following month:
  REAL(dp) :: photolysis_m3(n_photolytab_o2,n_photolytab_o3,nphotolysis,0:2)

  ! Photolysis rates of current time step:
  REAL(dp), DIMENSION(n_photolytab_o2,n_photolytab_o3), PUBLIC :: &
       do2a,       do3p,       do3d,       dno,        dno2,       dhno3,     &
       dno3a,      dno3b,      dn2o5,      dn2o,       dhno4,      dclno3a,   &
       dh2o2,      df11,       df12,       dhocl,      dh2o,       dco2,      &
       dcl2,       dcl2o2,     dch2oa,     dch2ob,     dch3o2h,    dbro,      &
       dbrno3a,    dbrno3b,    dbrcl,      dhobr,      dcbrf3,     do2b,      &
       dch4,       dclno3b,    dhcl,       dcfc113,    dcfc114,    dcfc115,   &
       dccl4,      dch3ccl3,   dhcfc22,    dhcfc141b,  dhcfc142b,  dh1211,    &
       dch3br,     dch3cl,     dhcfc21,    dhcfc123,   dh2402,     dchbr3,    &
       dch2br2,  &
       dpan, dmacr, dhac, dmgly, dpaa

  ! Photolysis rates for current time step and current grid point box. Values
  ! calculated by subroutine *socol_dis* (deduced from do2a, do3p,...).
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
       tjo2a,      tjo3p,      tjo3d,      tjno,       tjno2,      tjhno3,    &
       tjno3a,     tjno3b,     tjn2o5,     tjn2o,      tjhno4,     tjclno3a,  &
       tjh2o2,     tjf11,      tjf12,      tjhocl,     tjh2o,      tjco2,     &
       tjcl2,      tjcl2o2,    tjch2oa,    tjch2ob,    tjch3o2h,   tjbro,     &
       tjbrno3a,   tjbrno3b,   tjbrcl,     tjhobr,     tjcbrf3,    tjo2b,     &
       tjch4,      tjclno3b,   tjhcl,      tjcfc113,   tjcfc114,   tjcfc115,  &
       tjccl4,     tjch3ccl3,  tjhcfc22,   tjhcfc141b, tjhcfc142b, tjh1211,   &
       tjch3br,    tjch3cl,    tjhcfc21,   tjhcfc123,  tjh2402,    tjchbr3,   &
       tjch2br2,                                                              &
       tjpan, tjmacr, tjhac, tjmgly, tjpaa

  ! Delta-E corrections of photolysis rates of N2O5, N2O, F11, and F12 [-]:
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: xn2o5, xn2o, xf11, xf12

  LOGICAL, SAVE, PRIVATE :: lnot_used = .TRUE.

  PUBLIC :: allocate_socol_photolysis, read_socol_sun_data, &
       read_socol_photolysis, interpolate_socol_sun_data, &
       interpolate_socol_photolysis, srb_lya_heating, &
       cleanup_socol_photolysis

  ! Intrinsic functions: 
  INTRINSIC :: SUM, SQRT, EXP

CONTAINS

  SUBROUTINE allocate_socol_photolysis

    ! Allocates allocatable module variables.
    
    ! *allocate_socol_photolysis* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90

    IF (lnot_used) THEN
       
       ALLOCATE(tjo2a(nlev))           ; tjo2a(:)      = 0.0_dp
       ALLOCATE(tjo3p(nlev))           ; tjo3p(:)      = 0.0_dp
       ALLOCATE(tjo3d(nlev))           ; tjo3d(:)      = 0.0_dp
       ALLOCATE(tjno(nlev))            ; tjno(:)       = 0.0_dp
       ALLOCATE(tjno2(nlev))           ; tjno2(:)      = 0.0_dp
       ALLOCATE(tjhno3(nlev))          ; tjhno3(:)     = 0.0_dp

       ALLOCATE(tjno3a(nlev))          ; tjno3a(:)     = 0.0_dp
       ALLOCATE(tjno3b(nlev))          ; tjno3b(:)     = 0.0_dp
       ALLOCATE(tjn2o5(nlev))          ; tjn2o5(:)     = 0.0_dp
       ALLOCATE(tjn2o(nlev))           ; tjn2o(:)      = 0.0_dp
       ALLOCATE(tjhno4(nlev))          ; tjhno4(:)     = 0.0_dp
       ALLOCATE(tjclno3a(nlev))        ; tjclno3a(:)   = 0.0_dp

       ALLOCATE(tjh2o2(nlev))          ; tjh2o2(:)     = 0.0_dp
       ALLOCATE(tjf11(nlev))           ; tjf11(:)      = 0.0_dp
       ALLOCATE(tjf12(nlev))           ; tjf12(:)      = 0.0_dp
       ALLOCATE(tjhocl(nlev))          ; tjhocl(:)     = 0.0_dp
       ALLOCATE(tjh2o(nlev))           ; tjh2o(:)      = 0.0_dp
       ALLOCATE(tjco2(nlev))           ; tjco2(:)      = 0.0_dp

       ALLOCATE(tjcl2(nlev))           ; tjcl2(:)      = 0.0_dp
       ALLOCATE(tjcl2o2(nlev))         ; tjcl2o2(:)    = 0.0_dp
       ALLOCATE(tjch2oa(nlev))         ; tjch2oa(:)    = 0.0_dp
       ALLOCATE(tjch2ob(nlev))         ; tjch2ob(:)    = 0.0_dp
       ALLOCATE(tjch3o2h(nlev))        ; tjch3o2h(:)   = 0.0_dp
       ALLOCATE(tjbro(nlev))           ; tjbro(:)      = 0.0_dp

       ALLOCATE(tjbrno3a(nlev))        ; tjbrno3a(:)   = 0.0_dp
       ALLOCATE(tjbrno3b(nlev))        ; tjbrno3b(:)   = 0.0_dp
       ALLOCATE(tjbrcl(nlev))          ; tjbrcl(:)     = 0.0_dp
       ALLOCATE(tjhobr(nlev))          ; tjhobr(:)     = 0.0_dp
       ALLOCATE(tjcbrf3(nlev))         ; tjcbrf3(:)    = 0.0_dp
       ALLOCATE(tjo2b(nlev))           ; tjo2b(:)      = 0.0_dp

       ALLOCATE(tjch4(nlev))           ; tjch4(:)      = 0.0_dp
       ALLOCATE(tjclno3b(nlev))        ; tjclno3b(:)   = 0.0_dp
       ALLOCATE(tjhcl(nlev))           ; tjhcl(:)      = 0.0_dp
       ALLOCATE(tjcfc113(nlev))        ; tjcfc113(:)   = 0.0_dp
       ALLOCATE(tjcfc114(nlev))        ; tjcfc114(:)   = 0.0_dp
       ALLOCATE(tjcfc115(nlev))        ; tjcfc115(:)   = 0.0_dp

       ALLOCATE(tjccl4(nlev))          ; tjccl4(:)     = 0.0_dp
       ALLOCATE(tjch3ccl3(nlev))       ; tjch3ccl3(:)  = 0.0_dp
       ALLOCATE(tjhcfc22(nlev))        ; tjhcfc22(:)   = 0.0_dp
       ALLOCATE(tjhcfc141b(nlev))      ; tjhcfc141b(:) = 0.0_dp
       ALLOCATE(tjhcfc142b(nlev))      ; tjhcfc142b(:) = 0.0_dp
       ALLOCATE(tjh1211(nlev))         ; tjh1211(:)    = 0.0_dp

       ALLOCATE(tjch3br(nlev))         ; tjch3br(:)    = 0.0_dp
       ALLOCATE(tjch3cl(nlev))         ; tjch3cl(:)    = 0.0_dp
       ALLOCATE(tjhcfc21(nlev))        ; tjhcfc21(:)   = 0.0_dp
       ALLOCATE(tjhcfc123(nlev))       ; tjhcfc123(:)  = 0.0_dp
       ALLOCATE(tjh2402(nlev))         ; tjh2402(:)    = 0.0_dp
       ALLOCATE(tjchbr3(nlev))         ; tjchbr3(:)    = 0.0_dp
       ALLOCATE(tjch2br2(nlev))        ; tjch2br2(:)   = 0.0_dp

       ALLOCATE(tjpan(nlev))           ; tjpan(:)  = 0.0_dp
       ALLOCATE(tjmacr(nlev))         ; tjmacr(:) = 0.0_dp
       ALLOCATE(tjhac(nlev))           ; tjhac(:) = 0.0_dp
       ALLOCATE(tjmgly(nlev))         ; tjmgly(:) = 0.0_dp
       ALLOCATE(tjpaa(nlev))           ; tjpaa(:) = 0.0_dp

       lnot_used = .FALSE.
       
    ENDIF

  END SUBROUTINE allocate_socol_photolysis


  SUBROUTINE read_socol_sun_data(yr)

    ! Reads monthly values of prescribed solar irradiance of the SW spectral
    ! bands and Schumann-Runge-bands, Lyman alpha, and NO-photolysis parameters
    ! for months of current year plus December of preceeding / January of 
    ! following year.

    ! *read_socol_sun_data* is called from *read_socol_ghg_ods_global*.

    ! Subroutine arguments:
    INTEGER, intent(in) :: yr

    ! Local variables:
    INTEGER, PARAMETER :: nsw_data=6
    REAL(dp) :: sun_par_y(6,0:13)
    CHARACTER(31) :: fn 

    IF (lsolarvar) THEN

       ! Solar irradiance:
       IF (nsw_data .NE. nsw) THEN
          WRITE(message_text,*) 'Number of SW bands in data file: ', nsw_data
          CALL message('',TRIM(message_text))
          WRITE(message_text,*) 'Number of SW bands in ECHAM5: ', nsw
          CALL message('',TRIM(message_text))
          CALL finish('read_socol_sun_data','Run terminated')    
       ENDIF

       sun_irrad_y(:,:) = socol_ascii_yearmo_tab_y(nsw_data,yr,'sun_irrad_y')
       WRITE (message_text,*) &
            'Reading annually and monthly changing solar irradiance ', &
            'from file sun_irrad_y'
       CALL message('',TRIM(message_text))
             
       ! Schumann-Runge-bands, Lyman alpha, Hartley, Huggins and NO-photolysis parameter:
       sun_par_y(:,:) = socol_ascii_yearmo_tab_y(6,yr,'sun_par_y')     

       IF (lsrb_lya_heating) THEN
          sun_lya_y(:) = sun_par_y(1,:)
          sun_srb_y(:) = sun_par_y(2,:)
          sun_har_y(:)  = sun_par_y(3,:)
          sun_hug1_y(:) = sun_par_y(4,:)/250.0_dp  ! from w/m**2 to w/m**2/angstrom
          sun_hug2_y(:) = sun_par_y(5,:)/550.0_dp  ! from w/m**2 to w/m**2/angstrom
       ENDIF

       IF (lchem) sun_nopho_y(:) = sun_par_y(6,:)

       WRITE (message_text,*) &
            'Reading annually and monthly changing solar irradiance ', &
            'from file sun_irrad_y'
       CALL message('',TRIM(message_text))

       IF (lchem) THEN
          WRITE (message_text,*) &
               'Reading annually and monthly changing Schumann-Runge-bands, ', &
               'Lyman alpha, and NO-photolysis parameters from file sun_par_y'
       ELSE
          WRITE (message_text,*) &
               'Reading annually and monthly changing Schumann-Runge-bands, ', &
               'and Lyman alpha from file sun_par_y'
       ENDIF
       CALL message('',TRIM(message_text))

    ELSE   ! if lsolvar=.false.

       ! Schumann-Runge-bands, Lyman alpha, and NO-photolysis parameter:
       IF (lsrb_lya_heating) THEN
          sun_srb_y(1:12) = sun_srb_const(1:12)
          sun_srb_y(0) = sun_srb_const(12)
          sun_srb_y(13) = sun_srb_const(1)
          sun_lya_y(1:12) = sun_lya_const(1:12)
          sun_lya_y(0) = sun_lya_const(12)
          sun_lya_y(13) = sun_lya_const(1)
       ENDIF
       IF (lchem) THEN
          sun_nopho_y(1:12) = sun_nopho_const(1:12)
          sun_nopho_y(0) = sun_nopho_const(12)
          sun_nopho_y(13) = sun_nopho_const(1)
       ENDIF

       WRITE (message_text,*) &
            'Annually and monthly constant solar irradiance'
       CALL message('',TRIM(message_text))
       IF (lchem) THEN
          WRITE (message_text,*) &
               'Using annually constant Schumann-Runge-bands, and Lyman ', &
               'alpha'
       ELSE
          WRITE (message_text,*) &
               'Annually constant Schumann-Runge-bands, and Lyman alpha'
       ENDIF
       CALL message('',TRIM(message_text))


    ENDIF

    ! Delta-E corrections of photolysis rates of N2O5, N2O, F11, and F12 [-]:
    ! Caution, the corrections are factors - so for a "neutral" influence they
    ! have to be set to one, not zero.
    IF (lchem) THEN
       IF (deltaecorr) THEN
          fn = 'photolysis_delta_e_corr'
          CALL socol_read_netcdf(TRIM(fn), 'XN2O5', 'LEV', data1d=xn2o5)
          CALL socol_read_netcdf(TRIM(fn), 'XN2O',  'LEV', data1d=xn2o )
          CALL socol_read_netcdf(TRIM(fn), 'XF11',  'LEV', data1d=xf11 )
          CALL socol_read_netcdf(TRIM(fn), 'XF12',  'LEV', data1d=xf12 )
          WRITE (message_text,*) &
                  'Using Delta-E correction of photolysis rates.'
       ELSE
          ALLOCATE(xn2o5(nlev))        ; xn2o5(:)    = 1.0_dp
          ALLOCATE(xn2o(nlev))         ; xn2o(:)    = 1.0_dp
          ALLOCATE(xf11(nlev))         ; xf11(:)    = 1.0_dp
          ALLOCATE(xf12(nlev))         ; xf12(:)    = 1.0_dp
          WRITE (message_text,*) &
                  'Delta-E corrections are DEACTIVATED'
       ENDIF
       CALL message('',TRIM(message_text))
    ENDIF

  END SUBROUTINE read_socol_sun_data


  SUBROUTINE read_socol_photolysis(yr,mo)
    
    ! Reads photolysis rates for preceeding, current and following month from
    ! pre-calculated lookup-tables.

    ! *read_socol_photolysis* is called from *read_socol_bcond_m*, 
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo

    ! Local variables:
    INTEGER :: n_photolytab_o2_file, n_photolytab_o3_file, i, k
    REAL(dp), ALLOCATABLE :: photolysis_data(:,:,:)
    REAL(dp) :: xit(4)
    CHARACTER(31), PARAMETER :: fnvar = 'photolysis', fnmean = 'photolysis_mean'
    CHARACTER(31) :: varname


    ! Executable statements:

    ! Names of photolysis reactions (as in nc-file):
    photolysis_name(1) = 'do2a'
    photolysis_name(2) = 'do3p'
    photolysis_name(3) = 'do3d'
    photolysis_name(4) = 'dno'
    photolysis_name(5) = 'dno2'
  
    photolysis_name(6) = 'dhno3'
    photolysis_name(7) = 'dno3a'
    photolysis_name(8) = 'dno3b'
    photolysis_name(9) = 'dn2o5'
    photolysis_name(10) = 'dn2o'
  
    photolysis_name(11) = 'dhno4'
    photolysis_name(12) = 'dclno3a'
    photolysis_name(13) = 'dh2o2'
    photolysis_name(14) = 'df11'
    photolysis_name(15) = 'df12'
  
    photolysis_name(16) = 'dhocl'
    photolysis_name(17) = 'dh2o'
    photolysis_name(18) = 'dco2'
    photolysis_name(19) = 'dcl2'
    photolysis_name(20) = 'dcl2o2'
  
    photolysis_name(21) = 'dch2oa'
    photolysis_name(22) = 'dch2ob'
    photolysis_name(23) = 'dch3o2h'
    photolysis_name(24) = 'dbro'
    photolysis_name(25) = 'dbrno3a'
  
    photolysis_name(26) = 'dbrno3b'
    photolysis_name(27) = 'dbrcl'
    photolysis_name(28) = 'dhobr'
    photolysis_name(29) = 'dcbrf3'
    photolysis_name(30) = 'do2b'
  
    photolysis_name(31) = 'dch4'
    photolysis_name(32) = 'dclno3b'
    photolysis_name(33) = 'dhcl'
    photolysis_name(34) = 'dcfc113'
    photolysis_name(35) = 'dcfc114'
  
    photolysis_name(36) = 'dcfc115'
    photolysis_name(37) = 'dccl4'
    photolysis_name(38) = 'dch3ccl3'
    photolysis_name(39) = 'dhcfc22'
    photolysis_name(40) = 'dhcfc141b'
  
    photolysis_name(41) = 'dhcfc142b'
    photolysis_name(42) = 'dh1211'
    photolysis_name(43) = 'dch3br'
    photolysis_name(44) = 'dch3cl'
    photolysis_name(45) = 'dhcfc21'
  
    photolysis_name(46) = 'dhcfc123'
    photolysis_name(47) = 'dh2402'
    photolysis_name(48) = 'dchbr3'
    photolysis_name(49) = 'dch2br2'

    photolysis_name(50) = 'dpan'
    photolysis_name(51) = 'dmacr'
    photolysis_name(52) = 'dhac'
    photolysis_name(53) = 'dmgly'
    photolysis_name(54) = 'dpaa' 

    DO i = 1, nphotolysis

       ! Read photolysis rates from file:
       varname = TRIM(photolysis_name(i))
       IF (lsolarvar) THEN
          CALL socol_read_netcdf(fnvar, varname, '-', yr=yr, mo=mo, &
               varname_longname='-', extradimname='tabo2', &
               extradim2name='tabo3', data3d=photolysis_data, &
               nextradim=n_photolytab_o2_file, nextradim2=n_photolytab_o3_file)
       ELSE
          CALL socol_read_netcdf(fnmean, varname, '-', mo=mo, &
               varname_longname='-', extradimname='tabo2', &
               extradim2name='tabo3', data3d=photolysis_data, &
               nextradim=n_photolytab_o2_file, nextradim2=n_photolytab_o3_file)
       ENDIF

       ! Check dimensions:
       IF (i .EQ. 1) THEN
          IF (n_photolytab_o2_file .NE. n_photolytab_o2) THEN
             WRITE(message_text,*) 'Number of O2 tables in data file: ', &
                  n_photolytab_o2_file
             CALL message('',TRIM(message_text))
             WRITE(message_text,*) 'Number of O2 tables in ECHAM5: ', &
                  n_photolytab_o2
             CALL message('',TRIM(message_text))
             CALL finish('read_strataerosols','Run terminated')      
          ENDIF
          IF (n_photolytab_o3_file .NE. n_photolytab_o3) THEN
             WRITE(message_text,*) 'Number of O3 tables in data file: ', &
                  n_photolytab_o3_file
             CALL message('',TRIM(message_text))
             WRITE(message_text,*) 'Number of O3 tables in ECHAM5: ', &
                  n_photolytab_o3
             CALL message('',TRIM(message_text))
             CALL finish('read_strataerosols','Run terminated')      
          ENDIF
       ENDIF

       ! Allocate data to photolysis_m3:
       photolysis_m3(:,:,i,:) = photolysis_data(:,:,:)

    ENDDO

    ! Number of O2/O3 molecules passed by a solar beam [molec/cm^2] in
    ! look up table of the photolysis rates:
    xit(1) = 1.0_dp
    xit(2) = 2.0_dp
    xit(3) = 4.0_dp
    xit(4) = 7.0_dp    
    DO i = 1,4
       DO k = 1, 10
          xo2(i+(k-1)*4) = xit(i)*10.0_dp**REAL(k+17,dp) !40
       ENDDO
       DO k = 1, 11
          xo3(i+(k-1)*4) = xit(i)*10.0_dp**REAL(k+11,dp) !44
       ENDDO
    ENDDO    

    ! Print message:
    IF (p_parallel_io) THEN
       IF (lsolarvar) THEN
          WRITE (message_text,*) &
               'Reading annually and monthly changing precalulated ', &
               'photolysis rates' 
       ELSE
          WRITE (message_text,*) &
               'Reading annually constant monthly changing precalulated ', &
               'photolysis rates'
       ENDIF
       CALL message('',TRIM(message_text))
    ENDIF

  END SUBROUTINE read_socol_photolysis


  SUBROUTINE interpolate_socol_sun_data

    ! Interpolates monthly stratospheric solar data to the current time step.

    ! *interpolate_socol_sun_data* is called from 
    ! *interpolate_socol_bcond_global*, src/socol_boundary_conditions.

    REAL(dp) :: sun_irrad(nsw)

    IF (l_trigrad) THEN
       IF (lsolarvar) THEN
          sun_irrad(:) = &
               wgt1_rad*sun_irrad_y(:,yw1_rad) + wgt2_rad*sun_irrad_y(:,yw2_rad)
       ELSE
          sun_irrad(:) = sun_irrad_const(:)
       ENDIF
       IF (lsrb_lya_heating) THEN
          sun_srb = &
               wgt1_rad*sun_srb_y(yw1_rad) + wgt2_rad*sun_srb_y(yw2_rad)

          sun_lya = &
               wgt1_rad*sun_lya_y(yw1_rad) + wgt2_rad*sun_lya_y(yw2_rad)
         
          sun_har = &
               wgt1_rad*sun_har_y(yw1_rad) + wgt2_rad*sun_har_y(yw2_rad)
          
          sun_hug1 = &
               wgt1_rad*sun_hug1_y(yw1_rad) + wgt2_rad*sun_hug1_y(yw2_rad)
          
          sun_hug2 = &
               wgt1_rad*sun_hug2_y(yw1_rad) + wgt2_rad*sun_hug2_y(yw2_rad)
                         
       ENDIF
       
       ! Total irradiance (solar constant):
       solc = SUM(sun_irrad)
       ! Ratios:
       rsun = sun_irrad/solc 
    ENDIF

    IF (l_trigchem) THEN
       sun_nopho = &
            wgt1_chem*sun_nopho_y(yw1_chem) + wgt2_chem*sun_nopho_y(yw2_chem)
    ENDIF

  END SUBROUTINE interpolate_socol_sun_data


  SUBROUTINE interpolate_socol_photolysis

    ! Interpolates monthly stratospheric photolysis rates to the current 
    ! time step.

    ! *interpolate_socol_photolysis* is called from 
    ! *interpolate_socol_bcond_global*, src/socol_boundary_conditions.

    do2a(:,:)      = wgt1_chem * photolysis_m3(:,:,1,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,1,m3w2_chem)
    do3p(:,:)      = wgt1_chem * photolysis_m3(:,:,2,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,2,m3w2_chem)
    do3d(:,:)      = wgt1_chem * photolysis_m3(:,:,3,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,3,m3w2_chem)
    dno(:,:)       = wgt1_chem * photolysis_m3(:,:,4,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,4,m3w2_chem)
    dno2(:,:)      = wgt1_chem * photolysis_m3(:,:,5,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,5,m3w2_chem)

    dhno3(:,:)     = wgt1_chem * photolysis_m3(:,:,6,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,6,m3w2_chem) 
    dno3a(:,:)     = wgt1_chem * photolysis_m3(:,:,7,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,7,m3w2_chem)
    dno3b(:,:)     = wgt1_chem * photolysis_m3(:,:,8,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,8,m3w2_chem)
    dn2o5(:,:)     = wgt1_chem * photolysis_m3(:,:,9,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,9,m3w2_chem)
    dn2o(:,:)      = wgt1_chem * photolysis_m3(:,:,10,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,10,m3w2_chem)

    dhno4(:,:)     = wgt1_chem * photolysis_m3(:,:,11,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,11,m3w2_chem)
    dclno3a(:,:)   = wgt1_chem * photolysis_m3(:,:,12,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,12,m3w2_chem)
    dh2o2(:,:)     = wgt1_chem * photolysis_m3(:,:,13,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,13,m3w2_chem)
    df11(:,:)      = wgt1_chem * photolysis_m3(:,:,14,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,14,m3w2_chem)
    df12(:,:)      = wgt1_chem * photolysis_m3(:,:,15,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,15,m3w2_chem)

    dhocl(:,:)     = wgt1_chem * photolysis_m3(:,:,16,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,16,m3w2_chem)
    dh2o(:,:)      = wgt1_chem * photolysis_m3(:,:,17,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,17,m3w2_chem)
    dco2(:,:)      = wgt1_chem * photolysis_m3(:,:,18,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,18,m3w2_chem)
    dcl2(:,:)      = wgt1_chem * photolysis_m3(:,:,19,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,19,m3w2_chem)
    dcl2o2(:,:)    = wgt1_chem * photolysis_m3(:,:,20,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,20,m3w2_chem)

    dch2oa(:,:)    = wgt1_chem * photolysis_m3(:,:,21,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,21,m3w2_chem)
    dch2ob(:,:)    = wgt1_chem * photolysis_m3(:,:,22,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,22,m3w2_chem)
    dch3o2h(:,:)   = wgt1_chem * photolysis_m3(:,:,23,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,23,m3w2_chem)
    dbro(:,:)      = wgt1_chem * photolysis_m3(:,:,24,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,24,m3w2_chem)
    dbrno3a(:,:)   = wgt1_chem * photolysis_m3(:,:,25,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,25,m3w2_chem)

    dbrno3b(:,:)   = wgt1_chem * photolysis_m3(:,:,26,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,26,m3w2_chem)
    dbrcl(:,:)     = wgt1_chem * photolysis_m3(:,:,27,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,27,m3w2_chem)
    dhobr(:,:)     = wgt1_chem * photolysis_m3(:,:,28,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,28,m3w2_chem)
    dcbrf3(:,:)    = wgt1_chem * photolysis_m3(:,:,29,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,29,m3w2_chem)
    do2b(:,:)      = wgt1_chem * photolysis_m3(:,:,30,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,30,m3w2_chem)
    
    dch4(:,:)      = wgt1_chem * photolysis_m3(:,:,31,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,31,m3w2_chem)
    dclno3b(:,:)   = wgt1_chem * photolysis_m3(:,:,32,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,32,m3w2_chem)
    dhcl(:,:)      = wgt1_chem * photolysis_m3(:,:,33,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,33,m3w2_chem)
    dcfc113(:,:)   = wgt1_chem * photolysis_m3(:,:,34,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,34,m3w2_chem)
    dcfc114(:,:)   = wgt1_chem * photolysis_m3(:,:,35,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,35,m3w2_chem)

    dcfc115(:,:)   = wgt1_chem * photolysis_m3(:,:,36,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,36,m3w2_chem)
    dccl4(:,:)     = wgt1_chem * photolysis_m3(:,:,37,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,37,m3w2_chem)
    dch3ccl3(:,:)  = wgt1_chem * photolysis_m3(:,:,38,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,38,m3w2_chem)
    dhcfc22(:,:)   = wgt1_chem * photolysis_m3(:,:,39,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,39,m3w2_chem)
    dhcfc141b(:,:) = wgt1_chem * photolysis_m3(:,:,40,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,40,m3w2_chem)

    dhcfc142b(:,:) = wgt1_chem * photolysis_m3(:,:,41,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,41,m3w2_chem)
    dh1211(:,:)    = wgt1_chem * photolysis_m3(:,:,42,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,42,m3w2_chem)
    dch3br(:,:)    = wgt1_chem * photolysis_m3(:,:,43,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,43,m3w2_chem)
    dch3cl(:,:)    = wgt1_chem * photolysis_m3(:,:,44,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,44,m3w2_chem)
    dhcfc21(:,:)   = wgt1_chem * photolysis_m3(:,:,45,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,45,m3w2_chem)

    dhcfc123(:,:)  = wgt1_chem * photolysis_m3(:,:,46,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,46,m3w2_chem)
    dh2402(:,:)    = wgt1_chem * photolysis_m3(:,:,47,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,47,m3w2_chem)
    dchbr3(:,:)    = wgt1_chem * photolysis_m3(:,:,48,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,48,m3w2_chem)
    dch2br2(:,:)   = wgt1_chem * photolysis_m3(:,:,49,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,49,m3w2_chem)

    dpan(:,:) = wgt1_chem * photolysis_m3(:,:,50,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,50,m3w2_chem)
    dmacr(:,:) = wgt1_chem * photolysis_m3(:,:,51,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,51,m3w2_chem)
    dhac(:,:) = wgt1_chem * photolysis_m3(:,:,52,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,52,m3w2_chem)
    dmgly(:,:) = wgt1_chem * photolysis_m3(:,:,53,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,53,m3w2_chem)  
    dpaa(:,:) = wgt1_chem * photolysis_m3(:,:,54,m3w1_chem) +      &
         wgt2_chem * photolysis_m3(:,:,54,m3w2_chem)  

  END SUBROUTINE interpolate_socol_photolysis


    SUBROUTINE srb_lya_heating (krow, kproma, kbdim, klev, papm1, paphm1, ptte)

    ! Calculates heating by Schumann Runge bands (175-200 nm) and by Lyman 
    ! alpha line (121.6 nm) not taken into account by MA-ECHAM5 and 
	! extra heating due to solar variability in Hartley and Huggins bands.
    ! Parameterization based on Zhu(1994) and Nicolet(1985).

    ! *srb_lya_heating* is called from *physc*.

    ! Subroutine arguments:
    INTEGER, INTENT(in)     :: krow, kproma, kbdim, klev
	REAL(dp) ,INTENT(IN)    :: paphm1(kbdim,klev+1) ! half level pressure
    REAL(dp) ,INTENT(IN)    :: papm1(kbdim,klev)    ! full level pressure
    REAL(dp), INTENT(inout) :: ptte(kbdim,klev)

    ! Local constants and variables:
    INTEGER :: jl, jk
    REAL(dp), PARAMETER :: srb_lya_heating_lowbound = 1100_dp ! Apply parame-
                                                   ! terization for levels
                                                   ! above 1100 Pa (30 km)
    REAL(dp), PARAMETER :: o3_heating_lowbound = 20000_dp ! Apply parame-
                                                       ! terization for levels
                                                       ! above 20000 Pa (15 km)
                                                       
    REAL(dp), PARAMETER :: vclo2_above = 2.4E19_dp ! Part from above the model 
                                                   ! top
    REAL(dp), PARAMETER :: vclo3_above = 6.0E13_dp ! Part from above the model 
                                                   ! top
    REAL(dp), PARAMETER :: srb_a = 0.143_dp        ! Schumann-Runge constant
    REAL(dp), PARAMETER :: srb_b = 9.64E08_dp      ! Schumann-Runge constant
    REAL(dp), PARAMETER :: sigma_lya = 1.E-20_dp   ! Absorption cross section
                                                   ! of O2 for Lyman alpha line 

    REAL(dp) :: srb_heating, lya_heating, har_heating, hug_heating                
    REAL(dp) :: vclo2(klev), vclo3(klev), vmr_o3(kbdim,klev), yy1, yy2, vmr_o3c
    
    REAL(dp) :: vmr_o2, heat_cap, grav
    
    vmr_o2   = mro2 ! oxygen VMR
    heat_cap = cp   ! air heat capacity [J/kg/K]
    grav     = g    ! gravity accel [m/s**2] 

    vmr_o3(:,:)=xtm1(:,:,idt_o3,krow)

    DO jl = 1, kproma

       ! Calculations only for positive zenith angles:

       IF (amu0_x(jl,krow) .GT. 0.0_dp) THEN         

           vclo2(1) = vclo2_above + vmr_o2*(papm1(jl,1)-paphm1(jl,1))/grav/amu0_x(jl,krow)*1.0e-4_dp/(amo2*avoinv)         !slant mol/cm**2
           vclo3(1) = vclo3_above + vmr_o3(jl,1)*(papm1(jl,1)-paphm1(jl,1))/grav/amu0_x(jl,krow)*1.5e-4_dp/(amo3*avoinv)       !slant mol/cm**2 ???

          DO jk=2,klev
       ! Calculations of heating due to oxygen absorptio only above ~30 km:

             IF (papm1(jl,jk) .LE. srb_lya_heating_lowbound) THEN
                vclo2(jk) = vclo2(jk-1) + vmr_o2*(papm1(jl,jk)-papm1(jl,jk-1))/grav/amu0_x(jl,krow)*1e-4_dp/(amo2*avoinv)  !slant mol/cm**2
             ENDIF
             
             IF (papm1(jl,jk) .LE. O3_heating_lowbound) THEN
                vmr_o3c = (vmr_o3(jl,jk-1)+vmr_o3(jl,jk))*0.5_dp  ! o3 VMR on half-integer levels
                vclo3(jk) = vclo3(jk-1) + vmr_o3c*(papm1(jl,jk)-papm1(jl,jk-1))/grav/amu0_x(jl,krow)*1.5e-4_dp/(amo2*avoinv)  !slant mol/cm**2 ???
             ENDIF
          ENDDO
          
          
          DO jk=1,klev
          
             IF (papm1(jl,jk) .LE. srb_lya_heating_lowbound) THEN
                lya_heating = vmr_o2*sun_lya*(1.725E-18_dp/(vclo2(jk)**0.1175_dp))*1e-4_dp/(avoinv*amd)/heat_cap*exp(-2.115E-18_dp*vclo2(jk)**0.8855_dp)  ! K/s
                yy1 = sqrt(1.0_dp+4.0_dp*2.07e-20_dp*vclo2(jk)/api/0.0152_dp)
                yy2 = (vclo2(1)/vclo2(jk))**0.3_dp*2.07e-24_dp*2.6_dp


                ! Heating by Schumann Runge bands:
                srb_heating = vmr_o2*sun_srb*yy2/(avoinv*amd)/heat_cap/yy1*exp((-api*0.0152_dp/2._dp)*(yy1-1.0_dp))   ! ???
             
                ! Total heating [K/s]:
                ptte(jl,jk) = ptte(jl,jk) + srb_heating + lya_heating
             ENDIF ! Layer limitations
             
             IF (papm1(jl,jk) .LE. O3_heating_lowbound) THEN
                har_heating = sun_har*vmr_o3(jl,jk)*8.7e-22_dp*0.9_dp/(avoinv*amd)/heat_cap*exp(-8.7e-18_dp*vclo3(jk)) !??? some times sigma is in cm**2

                yy1 = 2.144693348e-17_dp           !exp(-0.01273*3015.0)
                yy2 = 3.10716734e-16_dp            !exp(-0.01273*2805.0)
                hug_heating = vmr_o3(jl,jk)/(vclo3(jk)*1e4_dp)/(avoinv*amd)/heat_cap*0.6_dp*&
                (sun_hug1+(sun_hug2-sun_hug1)*exp(-1.15e-2_dp*vclo3(jk)*yy1)-sun_hug2*exp(-1.15e-2_dp*vclo3(jk)*yy2))  !???
                ptte(jl,jk) = ptte(jl,jk) + har_heating + hug_heating
             ENDIF ! layer limitations
          ENDDO ! klev
       ENDIF   ! daytime

    ENDDO  !kproma


  END SUBROUTINE srb_lya_heating


  SUBROUTINE cleanup_socol_photolysis

    ! Deallocates allocatated module variables.
    
    ! *cleanup_socol_gcmfields* is called from *call_free_submodel_memory*,
    ! src/socol_photolysis.f90.

    IF (.NOT. lnot_used) THEN
       
       DEALLOCATE(tjo2a)
       DEALLOCATE(tjo3p)
       DEALLOCATE(tjo3d)
       DEALLOCATE(tjno)
       DEALLOCATE(tjno2)
       DEALLOCATE(tjhno3)
       
       DEALLOCATE(tjno3a)
       DEALLOCATE(tjno3b)
       DEALLOCATE(tjn2o5)
       DEALLOCATE(tjn2o)
       DEALLOCATE(tjhno4)
       DEALLOCATE(tjclno3a)
       
       DEALLOCATE(tjh2o2)
       DEALLOCATE(tjf11)
       DEALLOCATE(tjf12)
       DEALLOCATE(tjhocl)
       DEALLOCATE(tjh2o)
       DEALLOCATE(tjco2)
       
       DEALLOCATE(tjcl2)
       DEALLOCATE(tjcl2o2)
       DEALLOCATE(tjch2oa)
       DEALLOCATE(tjch2ob)
       DEALLOCATE(tjch3o2h)
       DEALLOCATE(tjbro)

       DEALLOCATE(tjbrno3a)
       DEALLOCATE(tjbrno3b)
       DEALLOCATE(tjbrcl)
       DEALLOCATE(tjhobr)
       DEALLOCATE(tjcbrf3)
       DEALLOCATE(tjo2b)
       
       DEALLOCATE(tjch4)
       DEALLOCATE(tjclno3b)
       DEALLOCATE(tjhcl)
       DEALLOCATE(tjcfc113)
       DEALLOCATE(tjcfc114)
       DEALLOCATE(tjcfc115)
       
       DEALLOCATE(tjccl4)
       DEALLOCATE(tjch3ccl3)
       DEALLOCATE(tjhcfc22)
       DEALLOCATE(tjhcfc141b)
       DEALLOCATE(tjhcfc142b)
       DEALLOCATE(tjh1211)
       
       DEALLOCATE(tjch3br)
       DEALLOCATE(tjch3cl)
       DEALLOCATE(tjhcfc21)
       DEALLOCATE(tjhcfc123)
       DEALLOCATE(tjh2402)
       DEALLOCATE(tjchbr3)     
       DEALLOCATE(tjch2br2)

       DEALLOCATE(tjpan)
       DEALLOCATE(tjmacr)
       DEALLOCATE(tjhac)
       DEALLOCATE(tjmgly)
       DEALLOCATE(tjpaa)

       lnot_used = .TRUE.
       
    ENDIF

    IF (ALLOCATED(xn2o5)) DEALLOCATE(xn2o5)
    IF (ALLOCATED(xn2o)) DEALLOCATE(xn2o)
    IF (ALLOCATED(xf11)) DEALLOCATE(xf11)
    IF (ALLOCATED(xf12)) DEALLOCATE(xf12)

  END SUBROUTINE cleanup_socol_photolysis

END MODULE mo_socol_sun

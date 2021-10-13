MODULE mo_socol_tropoaerosols

  ! Description:

  ! Organises tropospheric aerosols (reading data set; interpolation to
  ! current time step). Tropospheric aerosols are used for radiation. 
  ! Aerosol-cloud interaction and heterogeneous chemistry are not taken into
  ! account. The climatological data set of aerosol mass mixing ratio 
  ! distributions is based on Lohmann et al, 1999, JGR , and is representative 
  ! for 1980. As optical parameters the ones from GADS (already implemented
  ! in ECHAM5, see *mo_aero_gads*) are used. 
  !
  ! Aerosols types:
  !
  ! ======================================================================
  ! |    SPECIES        | RELATIVE HUMIDITY * | GADS AEROSOL-TYPE |
  ! ======================================================================
  ! MSA                        MODEL              11 = SUSO
  ! SULFATE                    MODEL              11 = SUSO
  ! DUST  0-1 mu               0%                 8 = MIAM
  ! DUST  1-2 mu               0%                 9 = MICM
  ! SEA SALT 0-1 mu            MODEL              5 = SSAM
  ! SEA SALT 1-10 mu           MODEL              6 = SSCM
  ! BLACK CARBON PHOBIC        0%                 3 = SOOT
  ! BLACK CARBON PHYLIC        MODEL              3 = SOOT
  ! ORGANIC CARBON PHO         0%                 1 = INSO
  ! ORGANIC CARBON PHY         MODEL              2 = WASO
  ! ======================================================================
  ! ======================================================================
  ! * MODEL: calculated by the model; 0%: fixed with 0%.
  !
  ! All fields are stored from top to bottom.
  ! 
  ! NB: This aerosol climatologogy is used if iaero is set to 10 in the 
  !     namelist *RADCTL*. Besides lgadsrh should be set to .TRUE. in the 
  !     namelist *RADCTL* to enable aerosol optical properties depending
  !     relative humidity.
  !
  ! M. Schraner, ETH Zürich, December 2008

  USE mo_aero_gads,               ONLY: fcvaer
  USE mo_constants,               ONLY: g
  USE mo_control,                 ONLY: nlon, ngl, nlev
  USE mo_decomposition,           ONLY: lc => local_decomposition, &
                                        global_decomposition
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_io
  USE mo_kind,                    ONLY: dp
  USE mo_mpi,                     ONLY: p_io, p_parallel, p_parallel_io, &
                                        p_bcast
  USE mo_radiation,               ONLY: ndfaer, newaer
  USE mo_socol_grid_calculations, ONLY: zaetr_socol
  USE mo_socol_interpo,           ONLY: wgt1_rad, wgt2_rad, m3w1_rad, m3w2_rad
  USE mo_socol_namelist
  USE mo_socol_readfile,          ONLY: socol_read_netcdf
  USE mo_time_control,            ONLY: l_trigrad
  USE mo_transpose,               ONLY: scatter_gp

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:
  
  ! Corresponding GADS aerosol types:

  INTEGER, PARAMETER, PUBLIC, DIMENSION(12) ::  &
       ndfaer_socol_tropoaer=(/11,11,8,9,5,6,3,3,1,2,0,0/)

  ! Vertical dimension of data set: first / last corresponding ECHAM-index:
  INTEGER :: k0_tropoaer, k1_tropoaer

  ! Preceeding, current and following month:
  REAL(dp), ALLOCATABLE :: mr_tropoaer_m3(:,:,:,:,:)  ! Aerosol mass mixing 
                                                      ! ratios [kg/kg]
  ! Interpolated to current time step:
  REAL(dp), ALLOCATABLE :: mr_tropoaer(:,:,:)         ! ditto

  PUBLIC :: read_socol_tropoaerosols, interpolate_socol_tropoaerosols, &
       cleanup_socol_tropoaerosols, socol_tropoaero

  ! Intrinsic functions: 
  INTRINSIC :: MAX, MIN

CONTAINS

  SUBROUTINE read_socol_tropoaerosols(yr,mo)

    ! Reads climatology of tropospheric aerosol mixing ratios [kg/kg].

    ! *read_socol_tropoaerosols* is called from *read_socol_bcond_m*,
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: mo, yr

    ! Local variables:
    INTEGER :: jt
    INTEGER, PARAMETER :: nvar = 10  ! Number of variables in data set
    REAL(dp), ALLOCATABLE :: tropoaer_data(:,:,:,:)
    CHARACTER(25), PARAMETER :: fntrop = 'tropoaer', fntropclim = 'tropoaerclim'
    CHARACTER(31) :: varname(nvar)
    CHARACTER(50) :: varname_longname(nvar)
    

    ! Executable statements:

    ! Check number of tropospheric aerosols:
    IF (newaer .NE. nvar) THEN
       WRITE(message_text,*) 'Unexpected number of tropospheric aerosols:', &
            newaer
       CALL message('',TRIM(message_text))
       CALL finish('read_socol_tropoaerosols','Run terminated')
    END IF

    ! Variable names:
    varname(1) = 'MSA'
    varname(2) = 'SULFATE'
    varname(3) = 'DUSTACCU'
    varname(4) = 'DUSTCOARSE'
    varname(5) = 'SEASALTACCU'
    varname(6) = 'SEASALTCOARSE'
    varname(7) = 'BCPHOBIC'
    varname(8) = 'BCPHILIC'
    varname(9) = 'OCPHOBIC'
    varname(10) = 'OCPHILIC'

    varname_longname(1) = 'methane sulfonate'
    varname_longname(2) = 'sulfate'
    varname_longname(3) = 'dust (0-1 um)'
    varname_longname(4) = 'dust (1-2 um)'
    varname_longname(5) = 'sea salt (0-1 um)'
    varname_longname(6) = 'sea salt (1-10 um)'
    varname_longname(7) = 'black carbon (hydrophobic)'
    varname_longname(8) = 'black carbon (hydrophilic)'
    varname_longname(9) = 'organic carbon (hydrophobic)'
    varname_longname(10) = 'organic carbon (hydrophilic)'

    ! Read tropospheric aerosols from data file:
    DO jt = 1, nvar

       IF (ltropaer_var_rad) THEN
          CALL socol_read_netcdf(fntrop, TRIM(varname(jt)), 'LONLATLEV', &
               yr=yr, mo=mo, varname_longname=TRIM(varname_longname(jt)), &
               data4d=tropoaer_data, lowlevind=k0_tropoaer, &
               uplevind=k1_tropoaer)
       ELSE
          CALL socol_read_netcdf(fntropclim, TRIM(varname(jt)), 'LONLATLEV', &
               mo=mo, varname_longname=TRIM(varname_longname(jt)), &
               data4d=tropoaer_data, lowlevind=k0_tropoaer, &
               uplevind=k1_tropoaer)
       ENDIF
       
       ! Allocate memory (if necessary):
       IF (.NOT. ALLOCATED(mr_tropoaer_m3)) &
            ALLOCATE(mr_tropoaer_m3(lc%nproma,k0_tropoaer:k1_tropoaer, &
            newaer,lc%ngpblks,0:2))
       
       ! Allocate data to mr_tropoaer_m3:
       mr_tropoaer_m3(:,:,jt,:,:) = tropoaer_data(:,:,:,:)
    END DO

    ! Deallocate memory:
    IF (ALLOCATED(tropoaer_data)) DEALLOCATE(tropoaer_data)

  END SUBROUTINE read_socol_tropoaerosols


  SUBROUTINE interpolate_socol_tropoaerosols(krow, kproma, kbdim)

    ! Interpolates monthly stratospheric aerosol data to the current time step.

    ! *interpolate_socol_tropoaerosols* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: krow, kproma, kbdim

    IF (.NOT. ALLOCATED(mr_tropoaer)) &
         ALLOCATE(mr_tropoaer(kbdim,k0_tropoaer:k1_tropoaer,newaer))
    
    mr_tropoaer(1:kproma,:,:) = &
         wgt1_rad*mr_tropoaer_m3(1:kproma,:,:,krow,m3w1_rad) + &
         wgt2_rad*mr_tropoaer_m3(1:kproma,:,:,krow,m3w2_rad)

  END SUBROUTINE interpolate_socol_tropoaerosols

  
  FUNCTION socol_tropoaero(kproma, kbdim, klev, pdp)

    ! Multiplication of aerosol mass mixing ratios by 
    ! (a) fcvaer*delta p/10/g (transformation [kg/kg] -> [part/cm2], i.e.
    !     aerosol amount per model layer. (Extinction coefficients of GADS
    !     aerosols have the unit [cm2/part], see modules/mo_aero_gads.f90.)) 
    ! (b) (1-zaetr_socol) (avoid tropospheric aerosols in the stratosphere).

    ! *socol_tropoaero* is called from *radiation*.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: kproma, kbdim, klev
    REAL(dp), INTENT(in) :: pdp(kbdim,klev)
    REAL(dp) :: socol_tropoaero(kproma,klev,newaer)

    ! Local variables:
    INTEGER :: jaer
    REAL(dp) :: inv10g

    inv10g = 0.1_dp/g

    socol_tropoaero(1:kproma,1:k0_tropoaer-1,:) = 0._dp
    IF (k1_tropoaer .LT. klev) &
         socol_tropoaero(1:kproma,k1_tropoaer+1:klev,:) = 0._dp  

    DO jaer = 1, newaer
       socol_tropoaero(1:kproma,k0_tropoaer:k1_tropoaer,jaer) = &
            fcvaer(ndfaer(jaer))*pdp(1:kproma,k0_tropoaer:k1_tropoaer) * &
            inv10g * &
            (1._dp-zaetr_socol(1:kproma,k0_tropoaer:k1_tropoaer)) * &
            mr_tropoaer(1:kproma,k0_tropoaer:k1_tropoaer,jaer)
    ENDDO

    ! MSA and sulfate data must be multiplied by 3 
    ! (Ulrike Lohmann, pers. communication, 2005): 
    socol_tropoaero(1:kproma,k0_tropoaer:k1_tropoaer,1:2) = &
         socol_tropoaero(1:kproma,k0_tropoaer:k1_tropoaer,1:2)*3._dp

  END FUNCTION socol_tropoaero


  SUBROUTINE cleanup_socol_tropoaerosols

    ! Deallocates module variables.
    
    ! *cleanup_socol_tropoaerosols* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(mr_tropoaer_m3)) DEALLOCATE(mr_tropoaer_m3)
    IF (ALLOCATED(mr_tropoaer)) DEALLOCATE(mr_tropoaer)

  END SUBROUTINE cleanup_socol_tropoaerosols

END MODULE mo_socol_tropoaerosols



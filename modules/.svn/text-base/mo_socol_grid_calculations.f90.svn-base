MODULE mo_socol_grid_calculations

  ! This module contains module variables of grid related quantities as
  ! - model levels for standard surface pressure
  ! - heigth of the model layers above sea level
  ! - air density of the grid boxes
  ! - solar angles
  ! - stratosphere/troposphere flag of the grid boxes
  ! - pressure level indices of the lower boundary of NAT and ice PSCs
  ! Besides it contains subroutines to determine these quantities.
  !
  ! M. Schraner, ETH Zurich, January 2009

  USE mo_constants,           ONLY: api, avo, amd, g, rd
  USE mo_control,             ONLY: nlon, ngl, nlev, nlevp1, nvclev, vct
  USE mo_geoloc,              ONLY: philon_2d, philat_2d
  USE mo_hyb,                 ONLY: nlmsgl, nlmslp, nplvp1, nplvp2 
  USE mo_io,                  ONLY: vlon, vlat
  USE mo_kind,                ONLY: dp
  USE mo_socol_constants,     ONLY: api180, psurf_stand
  USE mo_socol_namelist,      ONLY: hetnat_lowlev, hetice_lowlev
  USE mo_socol_time_control,  ONLY: l_trigchem, chemistry_date_yr, &
                                    chemistry_date_dyofyr, &
                                    chemistry_date_daypart, ndays_currentyear
  USE mo_time_control,        ONLY: l_trigrad
 
  IMPLICIT NONE

  PRIVATE

  ! Module variables:

  ! Pressure of model layers for standard surface pressure (half levels) [hPa]:
  REAL(dp), ALLOCATABLE, PUBLIC :: pres_stand_socol(:)

  ! Ditto for full levels [hPa]:
  REAL(dp), ALLOCATABLE, PUBLIC :: presf_stand_socol(:)  

  ! Altitude of model layers at model box boundaries [km] (no elevations
  ! included):
  REAL(dp), ALLOCATABLE, PUBLIC :: zlevb(:,:)

  ! Altitude of model layers at model box boundaries [km] for photolysis rates 
  ! calculations (no elevations included; shifted pressure indices):
  REAL(dp), ALLOCATABLE, PUBLIC :: zetb(:,:)

  ! Density of air molecules [molec/cm^3]:
  REAL(dp), ALLOCATABLE, PUBLIC :: dens(:,:)

  ! Stratosphere/ troposphere flag with continuos transition below tropopause:
  REAL(dp), ALLOCATABLE, PUBLIC :: zaetr_socol(:,:)

  ! Pressure level indices of the lower boundary of NAT and ice PSCs:
  INTEGER, ALLOCATABLE, PUBLIC :: hetice_lowlevind(:), hetnat_lowlevind(:)

  ! Cosine of zenith_angle:
  REAL(dp), ALLOCATABLE, PUBLIC :: cosz(:)  ! Cosine of zenith_angle

  ! Epsilon used for comparison of data and model grids:
  REAL(dp), PARAMETER, PUBLIC :: dim_eps_rel=0.001_dp
  REAL(dp), PARAMETER, PUBLIC :: dim_eps_abs=0.01_dp

  PUBLIC :: check_dim_socol, calculate_pres_stand_socol, &
       calculate_boundaries_altitude, calculate_solar_angles, calculate_dens, &
       calculate_zaetr_socol, calculate_hetpsc_lowlevind, &
       socol_interpolate_index, cleanup_socol_grid_calc

  ! Intrinsic functions:
  INTRINSIC :: ABS, MERGE, SUM

CONTAINS

  SUBROUTINE check_dim_socol(dimname, dim_f, dimerror, lowlevind, uplevind)

    ! Checks for dimension <dimname> if grid of the data file (vector <dim_f>)
    ! agrees with the corresponding model grid. For a vertical grid, 
    ! <lowlevind> and <uplevind> are determined.

    ! INPUT:
    !
    ! <dimname>     : (STRING): name of dimension. Valid are 'lon', lat', and 
    !                 'lev'
    ! <dim_f>       : (REAL(:)): <dimname> grid of data file (to be compared 
    !                 with model grid)
    !
    ! OUTPUT:
    !
    ! <dimerror>    : (INTEGER): dimerror is 0, if grids agree, and > 0 else.
    ! <lowlevind>, <uplevind> : (INTEGER, OPTIONAL): If the data array depends
    !                 on the pressure dimension (<dimname>='lev'), <lowlevind> /
    !                 <uplevind> indicate the lower / upper bounds of the 
    !                 vertical dimension of the output array such that the 
    !                 vertical indices of the model output array agree with the
    !                 ones of ECHAM5.

    ! Subroutine arguments:
    REAL(dp), INTENT(in), DIMENSION(:) :: dim_f
    CHARACTER(*), INTENT(in)           :: dimname
    INTEGER, INTENT(out)               :: dimerror
    INTEGER, OPTIONAL, INTENT(out)     :: lowlevind, uplevind

    ! Local variables:
    INTEGER  :: jk, ndim_f, ndim_e
    INTEGER, ALLOCATABLE  :: dimerrorv(:)
    REAL(dp), ALLOCATABLE :: dim_e(:)
    LOGICAL, ALLOCATABLE  :: ll(:)


    ! Executable statements:

    ndim_f=SIZE(dim_f,1)

    dimerror=0

    ! ECHAM5 grid:
    SELECT CASE (dimname)
    CASE ('lon')
       ndim_e=nlon
       ALLOCATE(dim_e(ndim_e))
       dim_e=vlon
    CASE ('lat')
       ndim_e=ngl
       ALLOCATE(dim_e(ndim_e))
       dim_e=vlat
    CASE ('lev')
       ! Determine uplevind, lowlevind:
       ! Upper index:
       lowlevind=-1
       DO jk=1, nlev
          IF (ABS(dim_f(1)-presf_stand_socol(jk))/presf_stand_socol(jk) &
               .LE. dim_eps_rel) THEN
             lowlevind=jk
             EXIT
          ENDIF
       ENDDO
       IF (lowlevind .EQ. -1) dimerror=1 
        
       ! Lower index:
       uplevind = lowlevind+ndim_f-1

       ndim_e=uplevind-lowlevind+1
       ALLOCATE(dim_e(ndim_e))
       dim_e(:)=presf_stand_socol(lowlevind:uplevind)
    CASE DEFAULT
       dimerror = 9999
    END SELECT

    IF (ndim_f .NE. ndim_e) dimerror=1
       
    IF (dimerror .EQ. 0) THEN
       ALLOCATE(ll(ndim_e))
       ALLOCATE(dimerrorv(ndim_e))

       IF (dimname .EQ. 'lev') THEN
          ! relative error:
          ll(:) = (ABS(dim_f(:)-dim_e(:))/dim_e(:) .GT. dim_eps_rel)
       ELSE
          ! absolute error:
          ll(:) = (ABS(dim_f(:)-dim_e(:)) .GT. dim_eps_abs)
       ENDIF
       dimerrorv(:) = MERGE(1,0,ll(:))

       dimerror = SUM(dimerrorv)
    ENDIF

    ! Deallocate memory:
    IF (ALLOCATED(dimerrorv)) DEALLOCATE (dimerrorv)
    IF (ALLOCATED(dim_e)) DEALLOCATE (dim_e)
    IF (ALLOCATED(ll)) DEALLOCATE (ll)

  END SUBROUTINE check_dim_socol
    

  SUBROUTINE calculate_pres_stand_socol
      
    ! Calculates pressure at model levels for standard surface pressure in 
    ! hPa.

    ! *calculate_pres_stand_socol* is called from *check_dim_netcdf*, 
    ! modules/mo_socol_readfile.f90.
    
    ! Local variables:
    INTEGER :: jk
    
    ! Executable statements:

    ! Return if pres_stand_socol and presf_stand_socol are already calculated:
    IF (ALLOCATED(pres_stand_socol) .AND. ALLOCATED(presf_stand_socol)) RETURN
    
    ! Allocate memory:
    IF (.NOT. ALLOCATED(pres_stand_socol)) ALLOCATE(pres_stand_socol(nlevp1))
    IF (.NOT. ALLOCATED(presf_stand_socol)) ALLOCATE(presf_stand_socol(nlev))
    
    ! Pressure level values:
    DO jk = 1, nplvp1
       pres_stand_socol(jk) = vct(jk)
    ENDDO
    
    ! Hybrid-level values:
    DO jk = nplvp2, nlmsgl
       pres_stand_socol(jk) = vct(jk) + vct(jk+nvclev)*psurf_stand
    ENDDO
    
    ! Sigma-level values:
    DO jk = nlmslp, nlevp1
       pres_stand_socol(jk) = vct(jk+nvclev)*psurf_stand
    ENDDO
    
    ! [Pa] -> [hPa]:
    pres_stand_socol = 0.01_dp*pres_stand_socol
    
    ! Full levels:
    DO jk = 1, nlev
       presf_stand_socol(jk) = &
            (pres_stand_socol(jk)+pres_stand_socol(jk+1))*0.5_dp
    ENDDO
    
  END SUBROUTINE calculate_pres_stand_socol
    
    
  SUBROUTINE calculate_boundaries_altitude (krow, kproma, kbdim, klev, &
       klevp1, p, ph, t)
    
    ! Calculates altitude of boundaries according to barometric formula
    ! dh=-dp/(rho(h)*g).

    ! *calculate_boundaries_altitude* is called from *physc* and from
    ! *mezon*, src/socol_mezon.f90.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim, klev, klevp1
    REAL(dp), INTENT(in) :: p(kbdim,klev), ph(kbdim,klevp1) 
    ! Pressure at full/half levels [Pa]
    REAL(dp), INTENT(in) :: t(kbdim,klev)    ! Temperature [K]
    
    ! Local variables:
    INTEGER  :: jl, jk
    REAL(dp) :: rho


    ! Executable statements:
    
    IF (.NOT. ALLOCATED(zlevb)) ALLOCATE(zlevb(kbdim,klevp1))
    IF (.NOT. ALLOCATED(zetb)) ALLOCATE(zetb(kbdim,klevp1+1))
    
    zlevb(1:kproma,klevp1) = 0._dp
    
    DO jl = 1, kproma
       DO jk = klev, 1, -1
          ! ideal gas equation: pV = nRT = nM(R/M)T = mR* T -> rho=p/(R* T),
          ! 
          ! where: n  number of moles
          !        R  universal gas constant
          !        R* = R/M gas constant of (dry) air
          !        M  molecular mass of air
          !        m  mass of air:
          rho = p(jl,jk)/(rd*t(jl,jk))
          
          zlevb(jl,jk) = zlevb(jl,jk+1) + (ph(jl,jk+1)-ph(jl,jk))/(rho*g)
       ENDDO
    ENDDO
    
    zlevb(1:kproma,:) = 0.001_dp*zlevb(1:kproma,:)   ! [m] -> [km]

    ! zetb (! shifted pressure indices !):
    zetb(1:kproma,1) = 200._dp   ! [km]
    zetb(1:kproma,2:klevp1+1) = zlevb(1:kproma,1:klevp1)
    
  END SUBROUTINE calculate_boundaries_altitude


  SUBROUTINE calculate_solar_angles (krow, kproma, kbdim, klev)
    
    ! Calculates solar angles for chemistry date.

    ! *calculate_solar_angles* is called from *mezon*, src/socol_mezon.f90.
    
    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim, klev

    ! Local variables:
    REAL(dp) :: rot, vernequinox, decl, cosf(kbdim)
    REAL, PARAMETER    :: api2 = 2._dp*api
    

    ! Executable statements:
    
    IF (.NOT. ALLOCATED(cosz)) ALLOCATE(cosz(kbdim))

    ! Earth's rotation angle of chemistry_date [radians]:
    rot = api2*chemistry_date_daypart
   
    ! Day of vernal equinox for current year (same parameterization as in
    ! *prerad*):
    vernequinox = 78.41_dp - 0.0078_dp*REAL(chemistry_date_yr-1987,dp) + &
         0.25_dp*MOD(chemistry_date_yr,4)

    ! Solar declination of chemistry_date:
    decl = 0.409_dp*SIN(api2*(chemistry_date_dyofyr - vernequinox)/ &
         ndays_currentyear)

    ! Cosine of hourly angle:
    cosf(1:kproma) = SIN(philon_2d(1:kproma,krow)*api180)*SIN(rot) - &
         COS(philon_2d(1:kproma,krow)*api180)*COS(rot)

    ! Cosine of zenith angle:
    cosz(1:kproma) = SIN(philat_2d(1:kproma,krow)*api180)*SIN(decl) + &
         COS(philat_2d(1:kproma,krow)*api180)*COS(decl)*cosf(1:kproma)

  END SUBROUTINE calculate_solar_angles
      

  SUBROUTINE calculate_dens (kproma, kbdim, klev, p, t)

    ! Calculates density of air molecules [molec/cm^3] according to ideal
    ! gas equation:
    !
    ! dens = N/V = (Na n)/V = (Na p)/(RT) = Na p/(R* MT)
    !
    ! where: N    number of molecules
    !        n    number of moles
    !        Na   number of molecules per mole (Avogadro constant)
    !        R    universal gas constant
    !        R* = R/M  gas constant of (dry) air
    !        M    molecular mass of an air molecule

    ! *calculate_dens* is called from *mezon*, src/socol_mezon.f90.
    
    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: kproma, kbdim, klev
    REAL(dp), INTENT(in) :: p(kbdim,klev) ! Pressure  [Pa]
    REAL(dp), INTENT(in) :: t(kbdim,klev) ! Temperature [K]

    ! Local variable:
    REAL(dp) :: fac

    
    ! Executable statements:
    
    IF (.NOT. ALLOCATED(dens)) ALLOCATE(dens(kbdim,klev))
    
    fac = avo/(rd*amd*1.E-03_dp)*1.E-06_dp  !1.E-03: M [g] -> M [kg]
    !1.E-06: [1/m^3] -> [1/cm^3]
    
    dens(1:kproma,:) = fac*p(1:kproma,:)/t(1:kproma,:)
    
  END SUBROUTINE calculate_dens


  SUBROUTINE calculate_zaetr_socol (kproma, kbdim, klev, ktrpwmop1)

    ! Calculates stratosphere/troposphere flag zaetr_socol: 
    ! zaetr = 0. in troposphere; zaetr = 1. in stratosphere; smooth
    ! transition in the three upper most tropospheric layers

    ! *calculate_zaetr_socol* is calculated from 
    ! *interpolate_socol_bcond*, src/socol_boundary_conditions.f90
    
    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: kproma, kbdim, klev
    INTEGER, INTENT(in)  :: ktrpwmop1(kbdim) ! Pressure index of uppermost 
                            ! level with level center below tropopause
    
    ! Local variables:
    INTEGER :: jl, jk
    

    ! Executable statements:

    IF (.NOT. ALLOCATED(zaetr_socol)) ALLOCATE(zaetr_socol(kbdim,klev))
    
    DO jl = 1, kproma
       zaetr_socol(jl,1:ktrpwmop1(jl)-1) = 1._dp      ! stratosphere
       
       ! Transition from stratospheric to tropospheric aerosol in the 
       ! three upper most tropospheric layers:
       zaetr_socol(jl,ktrpwmop1(jl)) = 0.6666666666666667_dp
       zaetr_socol(jl,ktrpwmop1(jl)+1) = 0.3333333333333333_dp
       
       zaetr_socol(jl,ktrpwmop1(jl)+2:klev) = 0._dp   ! troposphere
    ENDDO
    
  END SUBROUTINE calculate_zaetr_socol


  SUBROUTINE calculate_hetpsc_lowlevind (kproma, kbdim, klev, p, ktrpwmo)

    ! Determine the pressure level indices of the lower boundary of NAT and ice 
    ! PSCs.

    ! *calculate_hetpsc_lowlevind* is calculated from *mezon*.
    
    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: kproma, kbdim, klev
    REAL(dp), INTENT(in) :: p(kbdim,klev)    ! Pressure at full levels [Pa]
    INTEGER, INTENT(in)  :: ktrpwmo(kbdim)   ! Pressure index of pressure level
                                             ! containing the tropopause  
    ! Local variables:
    REAL(dp) :: lowlev_pa
    INTEGER :: jl, jk


    ! Executable statements:

    IF (.NOT. ALLOCATED(hetnat_lowlevind)) ALLOCATE(hetnat_lowlevind(kbdim))
    IF (.NOT. ALLOCATED(hetice_lowlevind)) ALLOCATE(hetice_lowlevind(kbdim))

    ! NAT:

    ! hPa -> Pa:
    lowlev_pa = 100._dp*hetnat_lowlev

    IF (lowlev_pa .GE. 1.E-06) THEN
       DO jl = 1, kproma
          DO jk = klev, 1, -1
             IF (p(jl,jk) .LT. lowlev_pa) EXIT
          ENDDO
          hetnat_lowlevind(jl) = jk
       ENDDO
    ELSE   ! layer containing tropopause
       hetnat_lowlevind(1:kproma) = ktrpwmo(1:kproma)
    ENDIF

    ! Ice:

    ! hPa -> Pa:
    lowlev_pa = 100._dp*hetice_lowlev

    IF (lowlev_pa .GE. 1.E-06) THEN
       DO jl = 1, kproma
          DO jk = klev, 1, -1
             IF (p(jl,jk) .LT. lowlev_pa) EXIT
          ENDDO
          hetice_lowlevind(jl) = jk
       ENDDO
    ELSE   ! layer containing tropopause
       hetice_lowlevind(1:kproma) = ktrpwmo(1:kproma)
    ENDIF

  END SUBROUTINE calculate_hetpsc_lowlevind


  SUBROUTINE socol_interpolate_index(g1, g2, k2first, k2last, k2index)
    
    ! This subroutine initializes the module variables which
    ! are used by the interpolation subroutines.
    
    ! Input:
    ! ------
    ! g1      : coordinates of grid 1
    ! g2      : coordinates of grid 2
    
    ! Output:
    ! -------
    ! k2first : index of first coordinate of grid 2
    !           in the range of grid 1
    ! k2last  : index of last coordinate of grid 2
    !           in the range of grid 1
    ! k2index : contains for each coordinate of grid2
    !           in the range of grid1 the index of the
    !           nearest bigger coordinate of grid1
    
    ! Data on grid 1 may be interpolated to the index
    ! range [k2first,k2last] of grid 2.

    ! *socol_interpolate_index* is called from *read_socol_qbo*.
    
    ! Authors: Marco Giorgetta, MPI for Meteorology, Hamburg, October 1999
    !          Martin Schraner, ETH Zurich, February 2009

    USE mo_kind,      ONLY: dp
    
    IMPLICIT NONE
    
    ! Subroutine arguments:
    REAL(dp), DIMENSION(:), INTENT(in) :: g1, g2
    INTEGER, INTENT(out)               :: k2first, k2last
    INTEGER, ALLOCATABLE, INTENT(out)  :: k2index(:)
   
    ! Local variables:
    INTEGER  :: ndim1, ndim2, k1, k2


    ! Executable statements:

    ndim1=SIZE(g1)
    ndim2=SIZE(g2)

    IF (.NOT. ALLOCATED(k2index)) ALLOCATE(k2index(ndim2))
    
    ! Find index of first element of grid 2
    ! which is in the range of grid 1
    ! -------------------------------------
    k2first=0
    k2=1
    DO
       IF (g2(k2) .GE. g1(1)) THEN
          k2first=k2
       ELSE
          k2=k2+1
       ENDIF
       IF (k2first .NE. 0) EXIT
    ENDDO

    ! Find index of last element of grid 2
    ! which is in the range of grid 1
    ! ------------------------------------
    k2last=0
    k2=ndim2
    DO
       IF (g2(k2) .LE. g1(ndim1)) THEN
          k2last=k2
       ELSE
          k2=k2-1
       ENDIF
       IF (k2last .NE. 0) EXIT
    ENDDO
    
    ! Find indices k2index for elements k2first
    ! to k2last of grid 2
    ! -----------------------------------------
    k2index(:)=0
    k2=k2first
    k1=2
    DO
       IF (g2(k2) .LE. g1(k1)) THEN
          k2index(k2)=k1
          k2=k2+1
       ELSE
          k1=k1+1
       ENDIF
       IF (k2 .GT. k2last) EXIT
    ENDDO
    
  END SUBROUTINE socol_interpolate_index

  
  SUBROUTINE cleanup_socol_grid_calc

    ! Deallocates module variables.
    
    ! *cleanup_socol_grid_calc* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(pres_stand_socol)) DEALLOCATE(pres_stand_socol)
    IF (ALLOCATED(presf_stand_socol)) DEALLOCATE(presf_stand_socol)
    IF (ALLOCATED(zlevb)) DEALLOCATE(zlevb)
    IF (ALLOCATED(zaetr_socol)) DEALLOCATE(zaetr_socol)

  END SUBROUTINE cleanup_socol_grid_calc
  
END MODULE mo_socol_grid_calculations

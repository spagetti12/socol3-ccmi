MODULE mo_socol_co_nox

  ! Description:

  ! Organises CO/NOx emissions at surface, by aircrafts and NOx produced by
  ! lightning.
  !
  ! M. Schraner, ETH Zurich, February 2009

  USE mo_kind,               ONLY: dp
  USE mo_socol_interpo,      ONLY: wgt1_chem, wgt2_chem, m3w1_chem, m3w2_chem
  USE mo_socol_namelist
  USE mo_socol_readfile,     ONLY: socol_read_netcdf
  USE mo_time_control,       ONLY: get_month_len

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  ! Pressure dimension of data set: first / last corresponding ECHAM-index:
  INTEGER, PUBLIC :: k0_aircr, k1_aircr, k0_lightn, k1_lightn

  ! Preceeding, current and following month:
  REAL(dp), ALLOCATABLE :: co_emiss_surf_m3(:,:,:)     ! CO surface emissions
                                                       ! [g/m^2/s]
  REAL(dp), ALLOCATABLE :: nox_emiss_surf_m3(:,:,:)    ! NOx surface emissions
                                                       ! [g/m^2/s]
  REAL(dp), ALLOCATABLE :: nox_emiss_aircr_m3(:,:,:,:) ! NOx aircraft emissions
                                                       ! [molec/cm^3/s]
  REAL(dp), ALLOCATABLE :: nox_lightn_col_m3(:,:,:)    ! NOx produced by 
                                                       ! lightning in a column
                                                       ! [g/m^2/s] 
  ! Interpolated to current time step:
  REAL(dp), ALLOCATABLE, PUBLIC :: co_emiss_surf(:)    ! [g/m^2/s]
  REAL(dp), ALLOCATABLE, PUBLIC :: nox_emiss_surf(:)   ! [g/m^2/s]
  REAL(dp), ALLOCATABLE, PUBLIC :: nox_emiss_aircr(:,:)  ! [molec/cm^3/s]
  REAL(dp), ALLOCATABLE, PUBLIC :: nox_lightn_col(:)   ! [g/m^2/s] 

  REAL(dp), ALLOCATABLE, PUBLIC :: nox_heppa(:,:,:)    ! Heppa nox [mol/mol]

  ! Vertical lightning profile:
  REAL(dp), ALLOCATABLE, PUBLIC :: nox_lightn_profile(:)

  PUBLIC :: read_socol_co_nox, interpolate_socol_co_nox, cleanup_socol_co_nox

CONTAINS

  SUBROUTINE read_socol_co_nox(yr,mo)

    ! Reads CO and NOx emissions for current, preceeding and following month.

    ! *read_socol_co_nox* is called from *read_socol_bcond_m*,
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo
    INTEGER :: monlen
    CHARACTER(LEN=100) :: year, month

    ! Local variables:
    CHARACTER(25), PARAMETER :: fn_surf_aircr0 = 'co_nox_emiss_surf_aircr', &
                                fn_lightn = 'nox_lightning'
    CHARACTER(31) :: varname, fn_heppa
    CHARACTER(50) :: varname_longname


    ! Executable statements:

    ! CO surface emissions:
    varname = 'CO_EMISS_SURF'
    varname_longname = 'CO surface emissions'
    IF (lsurfemco_var_chem) THEN
       CALL socol_read_netcdf(fn_surf_aircr0, varname, 'LONLAT', yr=yr, mo=mo, &
            varname_longname=varname_longname, data3d=co_emiss_surf_m3)
    ELSE
       CALL socol_read_netcdf(fn_surf_aircr0, varname, 'LONLAT', mo=mo, &
            varname_longname=varname_longname, data3d=co_emiss_surf_m3)
    ENDIF

    ! NOx surface emissions:
    varname = 'NOX_EMISS_SURF'
    varname_longname = 'NOx surface emissions'
    IF (lsurfemnox_var_chem) THEN
       CALL socol_read_netcdf(fn_surf_aircr0, varname, 'LONLAT', yr=yr, mo=mo, &
            varname_longname=varname_longname, data3d=nox_emiss_surf_m3)
    ELSE
       CALL socol_read_netcdf(fn_surf_aircr0, varname, 'LONLAT', mo=mo, &
            varname_longname=varname_longname, data3d=nox_emiss_surf_m3)
    ENDIF

    ! NOx aircraft emissions:
    varname = 'NOX_EMISS_AIRCR'
    varname_longname = 'NOx aircraft emissions'
    IF (laircrnox_var_chem) THEN
       CALL socol_read_netcdf(fn_surf_aircr0, varname, 'LONLATLEV', &
            yr=yr, mo=mo, varname_longname=varname_longname, &
            data4d=nox_emiss_aircr_m3, lowlevind=k0_aircr, uplevind=k1_aircr)
    ELSE
       CALL socol_read_netcdf(fn_surf_aircr0, varname, 'LONLATLEV', &
            mo=mo, varname_longname=varname_longname, &
            data4d=nox_emiss_aircr_m3, lowlevind=k0_aircr, uplevind=k1_aircr)
    ENDIF
 
    ! NOx production from lightning:
    varname = 'NOX_LIGHTN_COL'
    varname_longname = 'NOx production from lightning'
    CALL socol_read_netcdf(fn_lightn, varname, 'LONLAT', mo=mo, &
         varname_longname=varname_longname, data3d=nox_lightn_col_m3)

    varname = 'NOX_LIGHTN_PROFILE'
    varname_longname = 'vertical lightning profile'
    CALL socol_read_netcdf(fn_lightn, varname, 'LEV', &
         varname_longname=varname_longname, data1d=nox_lightn_profile, &
         lowlevind=k0_lightn, uplevind=k1_lightn)


    if ((yr.eq.2008).and.(mo.ge.10)) lheppa=.TRUE.
    if ((yr.eq.2009).and.(mo.le.5))  lheppa=.TRUE.

    If (lheppa) then        !Heppa nox  

      varname  = 'vmr'
      monlen=get_month_len(yr,mo)
      WRITE(year,  '(i8)') yr
          year = TRIM(adjustl(year))
      WRITE(month, '(i8)') mo
          month = TRIM(adjustl(month))
      fn_heppa = 'nox_heppa'//TRIM(year)//TRIM(month)
      if (mo.le.9) fn_heppa = 'nox_heppa'//TRIM(year)//'0'//TRIM(month)
      varname_longname = 'HEPPA uppermost layer NOx'
      CALL socol_read_netcdf(fn_heppa, varname, 'lonlat', &
           varname_longname=varname_longname, data3d=nox_heppa,nextradim=monlen,&
           extradimname='time')


!WB!----------------TB REMOVED----------------
!WB!     open(1333,file='no_heppa_reprint.txt') 
!WB!     write(1333,*)nox_heppa
!WB!     close(1333)
!WB!     open(1333,file='no_heppa_details.txt')
!WB!     write(1333,*)fn_heppa
!WB!     write(1333,*)SHAPE(nox_heppa)
!WB!     write(1333,*)monlen
!WB!     close(1333)
!WB!     stop
!WB!----------------TB REMOVED----------------


    endif

  END SUBROUTINE read_socol_co_nox


  SUBROUTINE interpolate_socol_co_nox(krow, kproma, kbdim)

    ! Interpolates monthly CO and NOx emissions to the current time step.

    ! *interpolate_socol_co_nox* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: krow, kproma, kbdim

    IF (.NOT. ALLOCATED(co_emiss_surf)) &
         ALLOCATE(co_emiss_surf(kbdim))
    IF (.NOT. ALLOCATED(nox_emiss_surf)) &
         ALLOCATE(nox_emiss_surf(kbdim))
    IF (.NOT. ALLOCATED(nox_emiss_aircr)) &
         ALLOCATE(nox_emiss_aircr(kbdim,k0_aircr:k1_aircr))
    IF (.NOT. ALLOCATED(nox_lightn_col)) &
         ALLOCATE(nox_lightn_col(kbdim))
           
    co_emiss_surf(1:kproma) = &
         wgt1_chem*co_emiss_surf_m3(1:kproma,krow,m3w1_chem) + &
         wgt2_chem*co_emiss_surf_m3(1:kproma,krow,m3w2_chem)
    nox_emiss_surf(1:kproma) = &
         wgt1_chem*nox_emiss_surf_m3(1:kproma,krow,m3w1_chem) + &
         wgt2_chem*nox_emiss_surf_m3(1:kproma,krow,m3w2_chem)
    nox_emiss_aircr(1:kproma,:) = &
         wgt1_chem*nox_emiss_aircr_m3(1:kproma,:,krow,m3w1_chem) + &
         wgt2_chem*nox_emiss_aircr_m3(1:kproma,:,krow,m3w2_chem)
    nox_lightn_col(1:kproma) = &
         wgt1_chem*nox_lightn_col_m3(1:kproma,krow,m3w1_chem) + &
         wgt2_chem*nox_lightn_col_m3(1:kproma,krow,m3w2_chem)

  END SUBROUTINE interpolate_socol_co_nox


  SUBROUTINE cleanup_socol_co_nox

    ! Deallocates module variables.
    
    ! *cleanup_socol_co_nox* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(co_emiss_surf_m3)) DEALLOCATE(co_emiss_surf_m3)
    IF (ALLOCATED(nox_emiss_surf_m3)) DEALLOCATE(nox_emiss_surf_m3)
    IF (ALLOCATED(nox_emiss_aircr_m3)) DEALLOCATE(nox_emiss_aircr_m3)
    IF (ALLOCATED(nox_lightn_col_m3)) DEALLOCATE(nox_lightn_col_m3)
    IF (ALLOCATED(co_emiss_surf)) DEALLOCATE(co_emiss_surf)
    IF (ALLOCATED(nox_emiss_surf)) DEALLOCATE(nox_emiss_surf)
    IF (ALLOCATED(nox_emiss_aircr)) DEALLOCATE(nox_emiss_aircr)
    IF (ALLOCATED(nox_lightn_col)) DEALLOCATE(nox_lightn_col)
    IF (ALLOCATED(nox_lightn_profile)) DEALLOCATE(nox_lightn_profile)
    if (lnudg_2h) then
      IF (ALLOCATED(nox_heppa)) DEALLOCATE(nox_heppa)
    endif
  END SUBROUTINE cleanup_socol_co_nox

END MODULE mo_socol_co_nox

  

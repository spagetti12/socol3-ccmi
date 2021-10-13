MODULE mo_socol_nmvoc

  ! Description:

  ! Organises NMVOC emissions from various sources at surface
  !
  ! A. Stenke, ETH Zurich, March 2010
  !
  ! A. Stenke, ETH Zurich, July 2010:
  ! Emissions for isoprene chemistry included

  USE mo_kind,               ONLY: dp
  USE mo_socol_interpo,      ONLY: wgt1_chem, wgt2_chem, m3w1_chem, m3w2_chem
  USE mo_socol_namelist
  USE mo_socol_readfile,     ONLY: socol_read_netcdf

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  ! Preceeding, current and following month:
  REAL(dp), ALLOCATABLE :: nmvoc_emiss_biogen_m3(:,:,:), & ! NMVOC biogen emissions, [g/m^2/s]
                           nmvoc_emiss_bb_m3(:,:,:), &     ! NMVOC biomass burning emissions, [g/m^2/s]
                           nmvoc_emiss_anthrop_m3(:,:,:)   ! NMVOC anthropogenic emissions, [g/m^2/s]

  REAL(dp), ALLOCATABLE :: c5h8_emiss_m3(:,:,:),  &        ! emissions for isoprene scheme
                           ch2o_emiss_m3(:,:,:),  &
                           hcooh_emiss_m3(:,:,:), &
                           ch3cooh_emiss_m3(:,:,:)

  ! Interpolated to current time step:
  REAL(dp), ALLOCATABLE, PUBLIC :: nmvoc_emiss_biogen(:), & ! NMVOC biogen emissions, [g/m^2/s]
                                   nmvoc_emiss_bb(:), &     ! NMVOC biomass burning emissions, [g/m^2/s]
                                   nmvoc_emiss_anthrop(:)   ! NMVOC anthropogenic emissions, [g/m^2/s]

  REAL(dp), ALLOCATABLE, PUBLIC :: c5h8_emiss(:),  &        ! emissions for isoprene scheme
                                   ch2o_emiss(:),  &
                                   hcooh_emiss(:), &
                                   ch3cooh_emiss(:)

  PUBLIC :: read_socol_nmvoc, interpolate_socol_nmvoc, cleanup_socol_nmvoc

CONTAINS

  SUBROUTINE read_socol_nmvoc(yr,mo)

    ! Reads NMVOC emissions for current, preceeding and following month.

    ! *read_socol_nmvoc* is called from *read_socol_bcond_m*,
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo

    ! Local variables:
    CHARACTER(25), PARAMETER :: fn_nmvoc_biogen = 'nmvoc_emiss_biogenic', &
                                fn_nmvoc_bb = 'nmvoc_emiss_bb', &
                                fn_nmvoc_anthrop = 'nmvoc_emiss_anthrop', &
                                fn_c5h8 = 'c5h8_emissions', &
                                fn_ch2o = 'ch2o_emissions', &
                                fn_hcooh = 'hcooh_emissions', &
                                fn_ch3cooh = 'ch3cooh_emissions'

    CHARACTER(31) :: varname
    CHARACTER(50) :: varname_longname


    ! Executable statements:

    IF (lnmvoc) THEN
       ! NMVOC biogen emissions:
       varname = 'NMVOC'
       varname_longname = 'biogenic emissions of NMHC (MEGANv2)'
       IF (lsurfemnmvoc_var_chem) THEN
          CALL socol_read_netcdf( fn_nmvoc_biogen, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=nmvoc_emiss_biogen_m3)
       ELSE
          CALL socol_read_netcdf( fn_nmvoc_biogen, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=nmvoc_emiss_biogen_m3)
       ENDIF
       
       ! NMVOC biomass burning emissions:
       varname = 'bb_emission_flux'
       varname_longname = 'biomass burning emissions of NMHC (ACCMIP)'
       IF (lsurfemnmvoc_var_chem) THEN
          CALL socol_read_netcdf( fn_nmvoc_bb, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=nmvoc_emiss_bb_m3)
       ELSE
          CALL socol_read_netcdf( fn_nmvoc_bb, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=nmvoc_emiss_bb_m3)
       ENDIF
       
       ! NMVOC anthropogenic emissions:
       varname = 'nmvoc'
       varname_longname = 'anthropogenic emissions of NMHC (ACCMIP)'
       IF (lsurfemnmvoc_var_chem) THEN
          CALL socol_read_netcdf( fn_nmvoc_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=nmvoc_emiss_anthrop_m3)
       ELSE
          CALL socol_read_netcdf( fn_nmvoc_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=nmvoc_emiss_anthrop_m3)
       ENDIF
    END IF

    ! C5H8 emissions
    varname = 'grassfire'
    varname_longname = 'biomass burning (ACCMIP) and biogenic (MEGANv2) emissions of isoprene (C5H8)'
    IF (lsurfemc5h8_var_chem) THEN
      CALL socol_read_netcdf( fn_c5h8, varname, 'lonlat', yr=yr, mo=mo, &
         varname_longname=varname_longname, data3d=c5h8_emiss_m3)
    ELSE
      CALL socol_read_netcdf( fn_c5h8, varname, 'lonlat', mo=mo, &
         varname_longname=varname_longname, data3d=c5h8_emiss_m3)
    ENDIF

    ! CH2O emissions
    varname = 'grassfire'
    varname_longname = 'biomass burning, anthropogenic (ACCMIP) and biogenic (MEGANv2) emissions of formaldehyde (CH2O)'
    IF (lsurfemch2o_var_chem) THEN
      CALL socol_read_netcdf( fn_ch2o, varname, 'lonlat',  yr=yr, mo=mo, &
         varname_longname=varname_longname, data3d=ch2o_emiss_m3)
    ELSE
      CALL socol_read_netcdf( fn_ch2o, varname, 'lonlat', mo=mo, &
         varname_longname=varname_longname, data3d=ch2o_emiss_m3)
    ENDIF

!!$    ! HCOOH emissions
!!$    varname = 'grassfire'
!!$    varname_longname = 'Grassland Fire Emissions of Formic Acid (HCOOH)'
!!$    CALL socol_read_netcdf( fn_hcooh, varname, 'lonlat', mo=mo, &
!!$         varname_longname=varname_longname, data3d=hcooh_emiss_m3)

    ! CH3COOH emissions
    varname = 'grassfire'
    varname_longname = 'forest, savanna and grassland fire emissions of acetic acid (CH3COOH)'
    IF (lsurfemch3cooh_var_chem) THEN
      CALL socol_read_netcdf( fn_ch3cooh, varname, 'lonlat', yr=yr, mo=mo, &
         varname_longname=varname_longname, data3d=ch3cooh_emiss_m3)
    ELSE
      CALL socol_read_netcdf( fn_ch3cooh, varname, 'lonlat', mo=mo, &
         varname_longname=varname_longname, data3d=ch3cooh_emiss_m3)
    ENDIF

  END SUBROUTINE read_socol_nmvoc


  SUBROUTINE interpolate_socol_nmvoc(krow, kproma, kbdim)

    ! Interpolates monthly NMVOC emissions to the current time step.

    ! *interpolate_socol_nmvoc* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: krow, kproma, kbdim

    IF (.NOT. ALLOCATED(nmvoc_emiss_biogen)) &
         ALLOCATE(nmvoc_emiss_biogen(kbdim))
    IF (.NOT. ALLOCATED(nmvoc_emiss_bb)) &
         ALLOCATE(nmvoc_emiss_bb(kbdim))
    IF (.NOT. ALLOCATED(nmvoc_emiss_anthrop)) &
         ALLOCATE(nmvoc_emiss_anthrop(kbdim))

    IF (.NOT. ALLOCATED(c5h8_emiss)) &
         ALLOCATE(c5h8_emiss(kbdim))
    IF (.NOT. ALLOCATED(ch2o_emiss)) &
         ALLOCATE(ch2o_emiss(kbdim))
    IF (.NOT. ALLOCATED(hcooh_emiss)) &
         ALLOCATE(hcooh_emiss(kbdim))
    IF (.NOT. ALLOCATED(ch3cooh_emiss)) &
         ALLOCATE(ch3cooh_emiss(kbdim))

    IF (lnmvoc) THEN         
       nmvoc_emiss_biogen(1:kproma) = &
            wgt1_chem*nmvoc_emiss_biogen_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*nmvoc_emiss_biogen_m3(1:kproma,krow,m3w2_chem)

       nmvoc_emiss_bb(1:kproma) = &
            wgt1_chem*nmvoc_emiss_bb_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*nmvoc_emiss_bb_m3(1:kproma,krow,m3w2_chem)
       
       nmvoc_emiss_anthrop(1:kproma) = &
            wgt1_chem*nmvoc_emiss_anthrop_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*nmvoc_emiss_anthrop_m3(1:kproma,krow,m3w2_chem)
    END IF
    
    c5h8_emiss(1:kproma) = &
         wgt1_chem*c5h8_emiss_m3(1:kproma,krow,m3w1_chem) + &
         wgt2_chem*c5h8_emiss_m3(1:kproma,krow,m3w2_chem)

    ch2o_emiss(1:kproma) = &
         wgt1_chem*ch2o_emiss_m3(1:kproma,krow,m3w1_chem) + &
         wgt2_chem*ch2o_emiss_m3(1:kproma,krow,m3w2_chem)

!!$    hcooh_emiss(1:kproma) = &
!!$         wgt1_chem*hcooh_emiss_m3(1:kproma,krow,m3w1_chem) + &
!!$         wgt2_chem*hcooh_emiss_m3(1:kproma,krow,m3w2_chem)

    ch3cooh_emiss(1:kproma) = &
         wgt1_chem*ch3cooh_emiss_m3(1:kproma,krow,m3w1_chem) + &
         wgt2_chem*ch3cooh_emiss_m3(1:kproma,krow,m3w2_chem)

  END SUBROUTINE interpolate_socol_nmvoc


  SUBROUTINE cleanup_socol_nmvoc

    ! Deallocates module variables.
    
    ! *cleanup_socol_nmvoc* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(nmvoc_emiss_biogen_m3))  DEALLOCATE(nmvoc_emiss_biogen_m3)
    IF (ALLOCATED(nmvoc_emiss_bb_m3))      DEALLOCATE(nmvoc_emiss_bb_m3)
    IF (ALLOCATED(nmvoc_emiss_anthrop_m3)) DEALLOCATE(nmvoc_emiss_anthrop_m3)

    IF (ALLOCATED(nmvoc_emiss_biogen))  DEALLOCATE(nmvoc_emiss_biogen)
    IF (ALLOCATED(nmvoc_emiss_bb))      DEALLOCATE(nmvoc_emiss_bb)
    IF (ALLOCATED(nmvoc_emiss_anthrop)) DEALLOCATE(nmvoc_emiss_anthrop)

    IF (ALLOCATED(c5h8_emiss_m3))     DEALLOCATE(c5h8_emiss_m3)
    IF (ALLOCATED(ch2o_emiss_m3))     DEALLOCATE(ch2o_emiss_m3)
    IF (ALLOCATED(hcooh_emiss_m3))    DEALLOCATE(hcooh_emiss_m3)
    IF (ALLOCATED(ch3cooh_emiss_m3))  DEALLOCATE(ch3cooh_emiss_m3)

    IF (ALLOCATED(c5h8_emiss))     DEALLOCATE(c5h8_emiss)
    IF (ALLOCATED(ch2o_emiss))     DEALLOCATE(ch2o_emiss)
    IF (ALLOCATED(hcooh_emiss))    DEALLOCATE(hcooh_emiss)
    IF (ALLOCATED(ch3cooh_emiss))  DEALLOCATE(ch3cooh_emiss)

  END SUBROUTINE cleanup_socol_nmvoc

END MODULE mo_socol_nmvoc

  

MODULE mo_socol_ch4

  ! Description:

  ! Organises CH4 emissions at surface.
  !
  ! Andrea Stenke, ETH Zurich, November 2009

  USE mo_kind,               ONLY: dp
  USE mo_socol_interpo,      ONLY: wgt1_chem, wgt2_chem, m3w1_chem, m3w2_chem, &
                                   yw1_chem, yw2_chem
  USE mo_socol_namelist
  USE mo_socol_readfile,     ONLY: socol_read_netcdf
  USE mo_decomposition,      ONLY: lc => local_decomposition, global_decomposition

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  REAL(dp), ALLOCATABLE, PUBLIC :: wetland_frac(:,:)           ! (nlon,ngl) in global coordinates
  REAL(dp), ALLOCATABLE, PUBLIC :: npp(:,:)                    ! (nlon,ngl) in global coordinates 
  REAL(dp), ALLOCATABLE, PUBLIC :: tsurf(:,:)                  ! (nlon,ngl) in global coordinates 

  REAL(dp), ALLOCATABLE, PUBLIC :: ch4_cmdl_m12(:,:)

  REAL(dp) :: ch4_cmdl

  LOGICAL, SAVE :: lnot_used2 = .TRUE.

 
  ! Preceeding, current and following month:
  ! IPCC emission inventories
  REAL(dp), ALLOCATABLE :: ch4_emiss_ship_m3(:,:,:), &         ! CH4 ship emissions [kg/m^2/s]
                           ch4_emiss_grassfire_m3(:,:,:), &    ! CH4 grassfire emissions [kg/m^2/s]
                           ch4_emiss_forestfire_m3(:,:,:), &   ! CH4 forestfire emissions [kg/m^2/s]
                           ch4_emiss_agr_m3(:,:,:), &          ! CH4 agriculture emissions [kg/m^2/s]
                           ch4_emiss_awb_m3(:,:,:), &          ! CH4 agricultural waste burning emissions [kg/m^2/s]
                           ch4_emiss_dom_m3(:,:,:), &          ! CH4 domestic emissions [kg/m^2/s]
                           ch4_emiss_ene_m3(:,:,:), &          ! CH4 energy emissions [kg/m^2/s]
                           ch4_emiss_ind_m3(:,:,:), &          ! CH4 industry emissions [kg/m^2/s]
                           ch4_emiss_tra_m3(:,:,:), &          ! CH4 traffic emissions [kg/m^2/s]
                           ch4_emiss_wst_m3(:,:,:), &          ! CH4 waste emissions [kg/m^2/s]
                           ch4_emiss_slv_m3(:,:,:)             ! CH4 solvent emissions [kg/m^2/s]

  ! Alternative emission inventories
  REAL(dp), ALLOCATABLE :: ch4_emiss_wetlands_m3(:,:,:), &     ! CH4 wetland emissions [kg/m^2/s]
                           ch4_emiss_bb_m3(:,:,:), &           ! CH4 fire emissions [kg/m^2/s]
                           ch4_emiss_avi_m3(:,:,:), &          ! CH4 aviation emissions [kg/m^2/s]
                           ch4_emiss_air_m3(:,:,:), &          ! CH4 aircraft emissions [kg/m^2/s]
                           ch4_emiss_trans_m3(:,:,:), &        ! CH4 transport emissions [kg/m^2/s]
                           ch4_emiss_intship_m3(:,:,:), &      ! CH4 international shipping emissions [kg/m^2/s]
                           ch4_emiss_resi_m3(:,:,:), &         ! CH4 residential emissions [kg/m^2/s]
                           ch4_emiss_fuels_m3(:,:,:), &        ! CH4 fuels emissions [kg/m^2/s]
                           ch4_emiss_agri_m3(:,:,:), &         ! CH4 agriculture emissions [kg/m^2/s]
                           ch4_emiss_waste_m3(:,:,:), &        ! CH4 waste emissions [kg/m^2/s]
                           ch4_emiss_wstwat_m3(:,:,:), &       ! CH4 waste water emissions [kg/m^2/s]
                           ch4_emiss_other_m3(:,:,:)           ! CH4 other emissions [kg/m^2/s]

  ! Interpolated to current time step:
  ! IPCC emission inventories    
  REAL(dp), ALLOCATABLE, PUBLIC :: ch4_emiss_ship(:), &        ! [kg/m^2/s]
                                   ch4_emiss_grassfire(:), &
                                   ch4_emiss_forestfire(:), &
                                   ch4_emiss_agr(:), &
                                   ch4_emiss_awb(:), &
                                   ch4_emiss_dom(:), &
                                   ch4_emiss_ene(:), &
                                   ch4_emiss_ind(:), &
                                   ch4_emiss_tra(:), &
                                   ch4_emiss_wst(:), &
                                   ch4_emiss_slv(:)

  ! Alternative emission inventories
  REAL(dp), ALLOCATABLE, PUBLIC :: ch4_emiss_wetlands(:), &     ! [kg/m^2/s]
                                   ch4_emiss_bb(:), &  
                                   ch4_emiss_avi(:), &        
                                   ch4_emiss_air(:), &        
                                   ch4_emiss_trans(:), &      
                                   ch4_emiss_intship(:), &
                                   ch4_emiss_resi(:), &       
                                   ch4_emiss_fuels(:), &      
                                   ch4_emiss_agri(:), &       
                                   ch4_emiss_waste(:), &      
                                   ch4_emiss_wstwat(:), &     
                                   ch4_emiss_other(:)          

  REAL(dp), ALLOCATABLE, PUBLIC :: ch4_emiss_wet(:)            ! [kg/m^2/s]      

  PUBLIC :: read_socol_ch4, read_wetland, read_npp, read_tsurf, interpolate_socol_ch4, &
       cleanup_socol_ch4, ch4_emiss_wetland, ch4_surface_flux, read_ch4_clim, &
       interpolate_ch4_clim, allocate_ch4_clim, cleanup_ch4_clim

CONTAINS

  SUBROUTINE read_socol_ch4(yr,mo)

    ! Reads CH4 emissions for current, preceeding and following month.

    ! *read_socol_ch4* is called from *read_socol_bcond_m*,
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo

    ! Local variables:
    CHARACTER(25), PARAMETER :: fn_ch4_anthrop = 'ch4_emiss_surf_anthrop', &
                                fn_ch4_ship = 'ch4_emiss_surf_ship', &
                                fn_ch4_bb = 'ch4_emiss_surf_bb', & 
                                fn_ch4_wet = 'ch4_emiss_surf_wet'                               
    CHARACTER(31) :: varname
    CHARACTER(50) :: varname_longname


    ! Executable statements:

    if (lch4_ipcc) then ! read in IPCC emission inventory 
       ! CH4 agricultural emissions:
       varname = 'emiss_agr'
       varname_longname = 'Agricultural sector emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_agr_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_agr_m3)
       ENDIF

       ! CH4 agricultural waste burning emissions:
       varname = 'emiss_awb'
       varname_longname = 'Agricultural waste burning emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_awb_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_awb_m3)
       ENDIF
       
       ! CH4 domestic combustion emissions:
       varname = 'emiss_dom'
       varname_longname = 'Residential and commercial combustion emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_dom_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_dom_m3)
       ENDIF
       
       ! CH4 energy production emissions:
       varname = 'emiss_ene'
       varname_longname = 'Energy production and distribution emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ene_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ene_m3)
       ENDIF
       
       ! CH4 industry emissions:
       varname = 'emiss_ind'
       varname_longname = 'Industrial processes and combustion emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ind_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ind_m3)
       ENDIF
       
       ! CH4 land transport emissions:
       varname = 'emiss_tra'
       varname_longname = 'Land transport emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_tra_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_tra_m3)
       ENDIF
       
       ! CH4 waste treatment emissions:
       varname = 'emiss_wst'
       varname_longname = 'Waste treatment and disposal emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_wst_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_wst_m3)
       ENDIF
       
       ! CH4 solvent production emissions:
       IF (yr .gt. 2000) THEN ! emissions from solvent production only in future emission data base
          varname = 'emiss_slv'
          varname_longname = 'Solvent production and use emissions for CH4'
          IF (lsurfemch4_var_chem) THEN
             CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
                  varname_longname=varname_longname, data3d=ch4_emiss_slv_m3)
          ELSE
             CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
                  varname_longname=varname_longname, data3d=ch4_emiss_slv_m3)
          ENDIF
       END IF

       ! CH4 ship emissions:
       varname = 'emiss_shp'
       varname_longname = 'Shipping emissions for CH4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_ship, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ship_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_ship, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ship_m3)
       ENDIF

       ! CH4 grassfire emissions:
       varname = 'grassfire'
       varname_longname = 'Grassland Fire Emissions of Methane (CH4)'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_bb, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_grassfire_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_bb, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_grassfire_m3)
       ENDIF
       
       ! CH4 forestfire emissions:
       varname = 'forestfire'
       varname_longname = 'Forest Fire Emissions of Methane (CH4)'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_bb, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_forestfire_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_bb, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_forestfire_m3)
       ENDIF

    ELSE ! alternative emission inventories

       ! CH4 wetland emissions
       varname = 'methane'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_wet, varname, 'lonlat', yr=yr, mo=mo, &
               data3d=ch4_emiss_wetlands_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_wet, varname, 'lonlat', mo=mo, &
               data3d=ch4_emiss_wetlands_m3)
       ENDIF

       ! CH4 fire emissions
       varname = 'ch4_bb'
       varname_longname = 'aggregated fire emissions for species ch4'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_bb, varname, 'lonlat', yr=yr, mo=mo, &
               data3d=ch4_emiss_bb_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_bb, varname, 'lonlat', mo=mo, &
               data3d=ch4_emiss_bb_m3)
       ENDIF

       ! CH4 energy
       varname = 'ch4_ene'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ene_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ene_m3)
       ENDIF  

       ! CH4 industry
       varname = 'ch4_ind'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ind_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ind_m3)
       ENDIF         

       ! CH4 aviation
       varname = 'ch4_avi'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_avi_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_avi_m3)
       ENDIF 

       ! CH4 aircraft
       varname = 'ch4_air'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_air_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_air_m3)
       ENDIF         

       ! CH4 transport
       varname = 'ch4_trans'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_trans_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_trans_m3)
       ENDIF   

       ! CH4 international shipping
       varname = 'ch4_intship'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_intship_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_intship_m3)
       ENDIF    

       ! CH4 shipping
       varname = 'ch4_ship'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ship_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_ship_m3)
       ENDIF   

       ! CH4 residential
       varname = 'ch4_resi'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_resi_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_resi_m3)
       ENDIF   
     
       ! CH4 fuels
       varname = 'ch4_fuels'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_fuels_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_fuels_m3)
       ENDIF   

       ! CH4 agriculture
       varname = 'ch4_agri'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_agri_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_agri_m3)
       ENDIF   

       ! CH4 waste
       varname = 'ch4_waste'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_waste_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_waste_m3)
       ENDIF

       ! CH4 wastewater
       varname = 'ch4_wstwat'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_wstwat_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_wstwat_m3)
       ENDIF 

       ! CH4 other
       varname = 'ch4_other'
       varname_longname = 'Emissions of CH4 -'
       IF (lsurfemch4_var_chem) THEN
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', yr=yr, mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_other_m3)
       ELSE
          CALL socol_read_netcdf(fn_ch4_anthrop, varname, 'lonlat', mo=mo, &
               varname_longname=varname_longname, data3d=ch4_emiss_other_m3)
       ENDIF 

    END IF

  END SUBROUTINE read_socol_ch4

  !-----------------------------------------------------------------------------------

  SUBROUTINE read_wetland

    ! A. Stenke, ETH, November 2009
    ! reads in wetland fraction for the calculation
    ! of methane emissions from wetlands
    ! [Percent of gridcell covered by wetland]
    ! *read_wetland* is called from stepon

    USE mo_control,       ONLY: nwlm
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_wetl(:,:)

    LOGICAL       :: lex
    INTEGER       :: nvarid

    IF (p_parallel_io) THEN

       CALL message('','CH4 emissions: Read in wetland map.')
       INQUIRE (nwlm, exist=lex)
       WRITE(message_text,*) 'lex: ', lex
       CALL message('read_wetland',message_text)
       IF (lex) THEN
          wetlnc1%format = NETCDF
          CALL IO_open_unit (nwlm, wetlnc1, IO_READ)
       ELSE
          CALL finish ('read_wetland', 'Could not open wetland-map file')
       ENDIF

    ENDIF

    !     Allocate memory for wetland_frac per PE

    IF (.NOT. ALLOCATED(wetland_frac)) ALLOCATE (wetland_frac(lc%nproma, lc%ngpblks))
    
    !     Read wetland-file
    IF (p_parallel_io) THEN
       
       !     Allocate memory for wetland_frac global fields
       
       ALLOCATE (zin(lc%nlon,lc%nlat,8))
       
       CALL IO_INQ_VARID (wetlnc1%file_id, 'wetland', nvarid)
       CALL IO_GET_VAR_DOUBLE (wetlnc1%file_id, nvarid, zin)

    END IF

    NULLIFY (gl_wetl)
    IF (p_pe == p_io) gl_wetl => zin(:,:,1)
    CALL scatter_gp (gl_wetl, wetland_frac(:,:), global_decomposition)
    
    IF (p_parallel_io) THEN
       DEALLOCATE (zin)
       
       !    Close file(s)
       
       CALL IO_close(wetlnc1)
       
    END IF
    
  END SUBROUTINE read_wetland

  !-----------------------------------------------------------------------------------

  SUBROUTINE read_npp

    ! A. Stenke, ETH, November 2009
    ! reads in climatological NPP as calculated by LPJ
    ! [kg(C)/m²]
    ! *read_npp* is called from stepon

    USE mo_control,       ONLY: nnpp
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:)
    REAL(dp), POINTER :: gl_npp(:,:)

    LOGICAL       :: lex
    INTEGER       :: nvarid

    IF (p_parallel_io) THEN

       CALL message('','CH4 emissions: Read in net primary production.')
       INQUIRE (nnpp, exist=lex)
       WRITE(message_text,*) 'lex: ', lex
       CALL message('read_npp',message_text)
       IF (lex) THEN
          nppnc1%format = NETCDF
          CALL IO_open_unit (nnpp, nppnc1, IO_READ)
       ELSE
          CALL finish ('read_npp', 'Could not open npp file')
       ENDIF

    ENDIF

    !     Allocate memory for npp per PE

    IF (.NOT. ALLOCATED(npp)) ALLOCATE (npp(lc%nproma, lc%ngpblks))
    
    !     Read npp-file
    IF (p_parallel_io) THEN
       
       !     Allocate memory for npp global fields
       
       ALLOCATE (zin(lc%nlon,lc%nlat))
       
!!$       CALL IO_INQ_VARID (nppnc1%file_id, 'NPP', nvarid)
       CALL IO_INQ_VARID (nppnc1%file_id, 'HR', nvarid)
       CALL IO_GET_VAR_DOUBLE (nppnc1%file_id, nvarid, zin(:,:))

    END IF

    NULLIFY (gl_npp)
    IF (p_pe == p_io) gl_npp => zin(:,:)
    CALL scatter_gp (gl_npp, npp(:,:), global_decomposition)

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)
       
       !    Close file(s)
       
       CALL IO_close(nppnc1)
       
    END IF

  END SUBROUTINE read_npp

  !-----------------------------------------------------------------------------------

  SUBROUTINE read_tsurf

    ! A. Stenke, ETH, November 2009
    ! reads in climatological annual mean surface temperature
    ! [K]
    ! *read_tsurf* is called from stepon

    USE mo_control,       ONLY: ntsurf
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:)
    REAL(dp), POINTER :: gl_tsurf(:,:)

    LOGICAL       :: lex
    INTEGER       :: nvarid

    IF (p_parallel_io) THEN

       CALL message('','CH4 emissions: Read in annual mean surface temperature.')
       INQUIRE (ntsurf, exist=lex)
       WRITE(message_text,*) 'lex: ', lex
       CALL message('read_tsurf',message_text)
       IF (lex) THEN
          tsnc1%format = NETCDF
          CALL IO_open_unit (ntsurf, tsnc1, IO_READ)
       ELSE
          CALL finish ('read_tsruf', 'Could not open surface temperature file')
       ENDIF

    ENDIF

    !     Allocate memory for tsurf per PE

    IF (.NOT. ALLOCATED(tsurf)) ALLOCATE (tsurf(lc%nproma, lc%ngpblks))
    
    !     Read npp-file
    IF (p_parallel_io) THEN
       
       !     Allocate memory for annual mean Tsurf
       
       ALLOCATE (zin(lc%nlon,lc%nlat))
       
       CALL IO_INQ_VARID (tsnc1%file_id, 'tslm1', nvarid)
       CALL IO_GET_VAR_DOUBLE (tsnc1%file_id, nvarid, zin(:,:))

    END IF

    NULLIFY (gl_tsurf)
    IF (p_pe == p_io) gl_tsurf => zin(:,:)
    CALL scatter_gp (gl_tsurf, tsurf(:,:), global_decomposition)

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)
       
       !    Close file(s)
       
       CALL IO_close(tsnc1)
       
    END IF

  END SUBROUTINE read_tsurf

  !-----------------------------------------------------------------------------------

  SUBROUTINE ch4_emiss_wetland(krow, kproma, kbdim)

    USE mo_geoloc, ONLY : gboxarea_2d 
    USE mo_socol_ch4_streams, ONLY: ch4emswet_d

    ! Subroutine arguments:
    INTEGER, INTENT(in)     :: krow, kproma, kbdim

    ! local parameters:
    INTEGER :: jl

    REAL(dp), PARAMETER :: ef = 0.005_dp ! emission factor
    REAL(dp), PARAMETER :: ms = 0.19_dp  ! moisture factor

    REAL(dp) :: tpeat, & ! low-emitting peatland
                tflood   ! high-emitting floodplain wetlands 
    REAL(dp) :: pl       ! factor

    IF (.NOT. ALLOCATED(ch4_emiss_wet)) &
         ALLOCATE(ch4_emiss_wet(kbdim))

    DO jl = 1, kproma

       IF (npp(jl,krow) .gt. 0._dp .and. wetland_frac(jl,krow) .gt. 0._dp) THEN 
 
          pl = EXP((tsurf(jl,krow) - 303._dp) / 8._dp)

          tflood = npp(jl,krow) * ms & 
                  * wetland_frac(jl,krow) * 0.01_dp & ! % -> fraction
                  / 31536000._dp ! -> kg(CH4)/m²/s  ???
          
          tpeat = npp(jl,krow) * ef & 
                  * wetland_frac(jl,krow) * 0.01_dp & ! % -> fraction
                  / 31536000._dp ! -> kg(CH4)/m²/s  ???

          ch4_emiss_wet(jl) = pl * tflood + (1. - pl) * tpeat

       ELSE
          
          ch4_emiss_wet(jl) = 0._dp
          
       END IF

       ch4emswet_d(jl,krow) = ch4_emiss_wet(jl)

    END DO

  END SUBROUTINE ch4_emiss_wetland

  !-----------------------------------------------------------------------------------

  SUBROUTINE interpolate_socol_ch4(krow, kproma, kbdim)

    ! Interpolates monthly CH4 emissions to the current time step.

    ! *interpolate_socol_ch4* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: krow, kproma, kbdim

    IF (.NOT. ALLOCATED(ch4_emiss_ship)) &
         ALLOCATE(ch4_emiss_ship(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_grassfire)) &
         ALLOCATE(ch4_emiss_grassfire(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_forestfire)) &
         ALLOCATE(ch4_emiss_forestfire(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_agr)) &
         ALLOCATE(ch4_emiss_agr(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_awb)) &
         ALLOCATE(ch4_emiss_awb(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_dom)) &
         ALLOCATE(ch4_emiss_dom(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_ene)) &
         ALLOCATE(ch4_emiss_ene(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_ind)) &
         ALLOCATE(ch4_emiss_ind(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_tra)) &
         ALLOCATE(ch4_emiss_tra(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_wst)) &
         ALLOCATE(ch4_emiss_wst(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_slv)) &
         ALLOCATE(ch4_emiss_slv(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_wetlands)) &
         ALLOCATE(ch4_emiss_wetlands(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_bb)) &
         ALLOCATE(ch4_emiss_bb(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_avi)) &
         ALLOCATE(ch4_emiss_avi(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_air)) &
         ALLOCATE(ch4_emiss_air(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_trans)) &
         ALLOCATE(ch4_emiss_trans(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_intship)) &
         ALLOCATE(ch4_emiss_intship(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_resi)) &
         ALLOCATE(ch4_emiss_resi(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_fuels)) &
         ALLOCATE(ch4_emiss_fuels(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_agri)) &
         ALLOCATE(ch4_emiss_agri(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_waste)) &
         ALLOCATE(ch4_emiss_waste(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_wstwat)) &
         ALLOCATE(ch4_emiss_wstwat(kbdim))

    IF (.NOT. ALLOCATED(ch4_emiss_other)) &
         ALLOCATE(ch4_emiss_other(kbdim))
           
    IF (ALLOCATED(ch4_emiss_ship_m3)) THEN
       ch4_emiss_ship(1:kproma) = &
            wgt1_chem*ch4_emiss_ship_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_ship_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_ship(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_grassfire_m3)) THEN
       ch4_emiss_grassfire(1:kproma) = &
            wgt1_chem*ch4_emiss_grassfire_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_grassfire_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_grassfire(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_forestfire_m3)) THEN
       ch4_emiss_forestfire(1:kproma) = &
            wgt1_chem*ch4_emiss_forestfire_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_forestfire_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_forestfire(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_agr_m3)) THEN
       ch4_emiss_agr(1:kproma) = &
            wgt1_chem*ch4_emiss_agr_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_agr_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_agr(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_awb_m3)) THEN
       ch4_emiss_awb(1:kproma) = &
            wgt1_chem*ch4_emiss_awb_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_awb_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_awb(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_dom_m3)) THEN
       ch4_emiss_dom(1:kproma) = &
            wgt1_chem*ch4_emiss_dom_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_dom_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_dom(1:kproma) = 0._dp
    END IF
       
    IF (ALLOCATED(ch4_emiss_ene_m3)) THEN
       ch4_emiss_ene(1:kproma) = &
            wgt1_chem*ch4_emiss_ene_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_ene_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_ene(1:kproma) = 0._dp
    END IF
    
    IF (ALLOCATED(ch4_emiss_ind_m3)) THEN
       ch4_emiss_ind(1:kproma) = &
            wgt1_chem*ch4_emiss_ind_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_ind_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_ind(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_tra_m3)) THEN
       ch4_emiss_tra(1:kproma) = &
            wgt1_chem*ch4_emiss_tra_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_tra_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_tra(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_wst_m3)) THEN
       ch4_emiss_wst(1:kproma) = &
            wgt1_chem*ch4_emiss_wst_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_wst_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_wst(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_slv_m3)) THEN  ! slv not presented in all emission files
       ch4_emiss_slv(1:kproma) = &
            wgt1_chem*ch4_emiss_slv_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_slv_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_slv(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_avi_m3)) THEN
       ch4_emiss_avi(1:kproma) = &
            wgt1_chem*ch4_emiss_avi_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_avi_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_avi(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_air_m3)) THEN
       ch4_emiss_air(1:kproma) = &
            wgt1_chem*ch4_emiss_air_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_air_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_air(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_trans_m3)) THEN
       ch4_emiss_trans(1:kproma) = &
            wgt1_chem*ch4_emiss_trans_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_trans_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_trans(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_intship_m3)) THEN
       ch4_emiss_intship(1:kproma) = &
            wgt1_chem*ch4_emiss_intship_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_intship_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_intship(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_resi_m3)) THEN
       ch4_emiss_resi(1:kproma) = &
            wgt1_chem*ch4_emiss_resi_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_resi_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_resi(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_fuels_m3)) THEN
       ch4_emiss_fuels(1:kproma) = &
            wgt1_chem*ch4_emiss_fuels_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_fuels_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_fuels(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_agri_m3)) THEN
       ch4_emiss_agri(1:kproma) = &
            wgt1_chem*ch4_emiss_agri_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_agri_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_agri(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_waste_m3)) THEN
       ch4_emiss_waste(1:kproma) = &
            wgt1_chem*ch4_emiss_waste_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_waste_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_waste(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_wstwat_m3)) THEN
       ch4_emiss_wstwat(1:kproma) = &
            wgt1_chem*ch4_emiss_wstwat_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_wstwat_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_wstwat(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_other_m3)) THEN
       ch4_emiss_other(1:kproma) = &
            wgt1_chem*ch4_emiss_other_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_other_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_other(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_wetlands_m3)) THEN
       ch4_emiss_wetlands(1:kproma) = &
            wgt1_chem*ch4_emiss_wetlands_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_wetlands_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_wetlands(1:kproma) = 0._dp
    END IF

    IF (ALLOCATED(ch4_emiss_bb_m3)) THEN
       ch4_emiss_bb(1:kproma) = &
            wgt1_chem*ch4_emiss_bb_m3(1:kproma,krow,m3w1_chem) + &
            wgt2_chem*ch4_emiss_bb_m3(1:kproma,krow,m3w2_chem)
    ELSE
       ch4_emiss_bb(1:kproma) = 0._dp
    END IF

  END SUBROUTINE interpolate_socol_ch4

  !-----------------------------------------------------------------------------------

  SUBROUTINE cleanup_socol_ch4

    ! Deallocates module variables.
    
    ! *cleanup_socol_ch4* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(ch4_emiss_ship_m3)) DEALLOCATE(ch4_emiss_ship_m3)
    IF (ALLOCATED(ch4_emiss_grassfire_m3)) DEALLOCATE(ch4_emiss_grassfire_m3)
    IF (ALLOCATED(ch4_emiss_forestfire_m3)) DEALLOCATE(ch4_emiss_forestfire_m3)
    IF (ALLOCATED(ch4_emiss_agr_m3)) DEALLOCATE(ch4_emiss_agr_m3)
    IF (ALLOCATED(ch4_emiss_awb_m3)) DEALLOCATE(ch4_emiss_awb_m3)
    IF (ALLOCATED(ch4_emiss_dom_m3)) DEALLOCATE(ch4_emiss_dom_m3)
    IF (ALLOCATED(ch4_emiss_ene_m3)) DEALLOCATE(ch4_emiss_ene_m3)
    IF (ALLOCATED(ch4_emiss_ind_m3)) DEALLOCATE(ch4_emiss_ind_m3)
    IF (ALLOCATED(ch4_emiss_tra_m3)) DEALLOCATE(ch4_emiss_tra_m3)
    IF (ALLOCATED(ch4_emiss_wst_m3)) DEALLOCATE(ch4_emiss_wst_m3)
    IF (ALLOCATED(ch4_emiss_slv_m3)) DEALLOCATE(ch4_emiss_slv_m3)

    IF (ALLOCATED(ch4_emiss_wetlands_m3)) DEALLOCATE(ch4_emiss_wetlands_m3)
    IF (ALLOCATED(ch4_emiss_bb_m3)) DEALLOCATE(ch4_emiss_bb_m3)
    IF (ALLOCATED(ch4_emiss_avi_m3)) DEALLOCATE(ch4_emiss_avi_m3)
    IF (ALLOCATED(ch4_emiss_air_m3)) DEALLOCATE(ch4_emiss_air_m3)
    IF (ALLOCATED(ch4_emiss_trans_m3)) DEALLOCATE(ch4_emiss_trans_m3)
    IF (ALLOCATED(ch4_emiss_intship_m3)) DEALLOCATE(ch4_emiss_intship_m3)
    IF (ALLOCATED(ch4_emiss_resi_m3)) DEALLOCATE(ch4_emiss_resi_m3)
    IF (ALLOCATED(ch4_emiss_fuels_m3)) DEALLOCATE(ch4_emiss_fuels_m3)
    IF (ALLOCATED(ch4_emiss_agri_m3)) DEALLOCATE(ch4_emiss_agri_m3)
    IF (ALLOCATED(ch4_emiss_waste_m3)) DEALLOCATE(ch4_emiss_waste_m3)
    IF (ALLOCATED(ch4_emiss_wstwat_m3)) DEALLOCATE(ch4_emiss_wstwat_m3)
    IF (ALLOCATED(ch4_emiss_other_m3)) DEALLOCATE(ch4_emiss_other_m3)

    IF (ALLOCATED(ch4_emiss_ship)) DEALLOCATE(ch4_emiss_ship)
    IF (ALLOCATED(ch4_emiss_grassfire)) DEALLOCATE(ch4_emiss_grassfire)
    IF (ALLOCATED(ch4_emiss_forestfire)) DEALLOCATE(ch4_emiss_forestfire)
    IF (ALLOCATED(ch4_emiss_agr)) DEALLOCATE(ch4_emiss_agr)
    IF (ALLOCATED(ch4_emiss_awb)) DEALLOCATE(ch4_emiss_awb)
    IF (ALLOCATED(ch4_emiss_dom)) DEALLOCATE(ch4_emiss_dom)
    IF (ALLOCATED(ch4_emiss_ene)) DEALLOCATE(ch4_emiss_ene)
    IF (ALLOCATED(ch4_emiss_ind)) DEALLOCATE(ch4_emiss_ind)
    IF (ALLOCATED(ch4_emiss_tra)) DEALLOCATE(ch4_emiss_tra)
    IF (ALLOCATED(ch4_emiss_wst)) DEALLOCATE(ch4_emiss_wst)
    IF (ALLOCATED(ch4_emiss_slv)) DEALLOCATE(ch4_emiss_slv)

    IF (ALLOCATED(ch4_emiss_wetlands)) DEALLOCATE(ch4_emiss_wetlands)
    IF (ALLOCATED(ch4_emiss_bb)) DEALLOCATE(ch4_emiss_bb)
    IF (ALLOCATED(ch4_emiss_avi)) DEALLOCATE(ch4_emiss_avi)
    IF (ALLOCATED(ch4_emiss_air)) DEALLOCATE(ch4_emiss_air)
    IF (ALLOCATED(ch4_emiss_trans)) DEALLOCATE(ch4_emiss_trans)
    IF (ALLOCATED(ch4_emiss_intship)) DEALLOCATE(ch4_emiss_intship)
    IF (ALLOCATED(ch4_emiss_resi)) DEALLOCATE(ch4_emiss_resi)
    IF (ALLOCATED(ch4_emiss_fuels)) DEALLOCATE(ch4_emiss_fuels)
    IF (ALLOCATED(ch4_emiss_agri)) DEALLOCATE(ch4_emiss_agri)
    IF (ALLOCATED(ch4_emiss_waste)) DEALLOCATE(ch4_emiss_waste)
    IF (ALLOCATED(ch4_emiss_wstwat)) DEALLOCATE(ch4_emiss_wstwat)
    IF (ALLOCATED(ch4_emiss_other)) DEALLOCATE(ch4_emiss_other)

    IF (lch4_wetland) THEN
       IF (ALLOCATED(ch4_emiss_wet)) DEALLOCATE(ch4_emiss_wet)
    END IF

  END SUBROUTINE cleanup_socol_ch4

  !-----------------------------------------------------------------------------------

  SUBROUTINE ch4_surface_flux(krow, kproma, kbdim, klev, ktrac, pxtm1, pxtte)

    ! SR constrain the calculated surface layer CH4 concentrations to
    ! a read-in CMDL CH4 climatology. 
    ! The calculated CH4 surface emission flux is copied to the output stream.

    ! *ch4_surface_flux* is called from *mezon* 
    !
    ! A. Stenke, ETH Zurich, April 2011

    USE mo_exception,               ONLY: finish
    USE mo_kind,                    ONLY: dp
    USE mo_socol_grid_calculations, ONLY: dens, zlevb
    USE mo_socol_tracers,           ONLY: idt_ch4
    USE mo_time_control,            ONLY: time_step_len
    USE mo_socol_ch4_streams,       ONLY: ch4_emflux_d
    
    ! Subroutine arguments:
    INTEGER, INTENT(in)     :: krow, kproma, kbdim, klev, ktrac

    REAL(dp), INTENT(in)    :: pxtm1(kbdim,klev,ktrac)
    REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)

    ! Local variables:
    INTEGER  :: jl
    REAL, PARAMETER :: trelax_ch4 = 21600._dp ! 6 hours
    REAL(dp) :: zxtm1, zdz



    ! calculate CH4 tendency resulting from the difference of the calculated and 
    ! "fixed" surface layer volume mixing ratio
    ! 
    ! the strength of the forcing towards the perscribed CH4 surface layer vmr 
    ! is determined by using a relaxation coefficient *trelax_ch4*;
    ! currently *trelax_ch4* is set to 6 hours; the relaxation is used to avoid
    ! numerical problems in the polar regions probably due to advection problems...

    DO jl = 1, kproma

       zxtm1 = ch4_cmdl

       IF (zxtm1 .lt. 1.e-20_dp) &
            CALL finish('ch4_surface_flux','Initial CH4 volume mixing ratio < 0 !')


       ! calculate surface emission flux from the difference of the calculated and 
       ! "fixed" surface layer volume mixing ratio
       ! * dens -> molec(CH4)/cm3
       ! * 1.e6 -> molec(CH4)/m3
       ! * zdz  -> molec(CH4)/m2
       ! / time_step_len [s] -> [molec/m2/s]

       ! thickness of surface layer

       zdz = (zlevb(jl,klev)-zlevb(jl,klev+1))*1000._dp ! [m]

       ch4_emflux_d(jl,krow) = (zxtm1 - &
            (pxtm1(jl,klev,idt_ch4) + pxtte(jl,klev,idt_ch4) * time_step_len)) * &
            dens(jl,klev) * 1.e6_dp * zdz / time_step_len

       ! update CH4 tendency
       !!pxtte(jl,klev,idt_ch4) = pxtte(jl,klev,idt_ch4) - &
       !!     (pxtm1(jl,klev,idt_ch4) + pxtte(jl,klev,idt_ch4) * time_step_len - &
       !!     zxtm1) / trelax_ch4

    END DO
  


  END SUBROUTINE ch4_surface_flux
  
  !-----------------------------------------------------------------------------------

  SUBROUTINE allocate_ch4_clim

    USE mo_control,       ONLY: ngl
    
    ! Allocates module variables.
    
    ! *allocate_ch4_clim* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90
    
    IF (lnot_used2) THEN

       ALLOCATE (ch4_cmdl_m12(ngl, 0:13))
       
       lnot_used2 = .FALSE.
       
    ENDIF
    
  END SUBROUTINE allocate_ch4_clim
  !----------------------------------------------------------------------------------- 
  
  SUBROUTINE cleanup_ch4_clim
    
    ! Deallocates module variables.
    
    ! *cleanup_ch4_clim* is called from *call_free_submodel_memory*,
    ! src/socol_submodels.f90.
    
    IF (.NOT. lnot_used2) THEN
       
       DEALLOCATE(ch4_cmdl_m12)
       
       lnot_used2 = .TRUE.
       
    ENDIF
    
  END SUBROUTINE cleanup_ch4_clim
  !-----------------------------------------------------------------------------------

  SUBROUTINE read_ch4_clim

    ! A. Stenke, ETH, May 2011
    ! reads in zonal mean surface CH4 climatology, CMDL
    ! [mol/mol]
    ! *read_ch4_clim* is called from stepon

    USE mo_control,       ONLY: nch4, ngl
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io, p_parallel, p_bcast, p_io

    REAL(dp), ALLOCATABLE :: zin(:,:,:,:)

    LOGICAL       :: lex
    INTEGER       :: nvarid

    IF (p_parallel_io) THEN

       CALL message('','CH4 emissions: Read in CMDL CH4 climatology.')
       INQUIRE (nch4, exist=lex)
       WRITE(message_text,*) 'lex: ', lex
       CALL message('read_ch4_clim',message_text)
       IF (lex) THEN
          ch4nc1%format = NETCDF
          CALL IO_open_unit (nch4, ch4nc1, IO_READ)
       ELSE
          CALL finish ('read_ch4_clim', 'Could not open CMDL CH4 climatology file')
       ENDIF

    ENDIF
    
    ! Read CMDL CH4 climatology
    IF (p_parallel_io) THEN
       
       ! Allocate memory
       
       ALLOCATE (zin(1,ngl,1,12))
       
       CALL IO_INQ_VARID (ch4nc1%file_id, 'CH4', nvarid)
       CALL IO_GET_VAR_DOUBLE (ch4nc1%file_id, nvarid, zin(:,:,:,:))
       
       ch4_cmdl_m12(:,1:12) = zin(1,:,1,:)
       ch4_cmdl_m12(:,0) = zin(1,:,1,12)
       ch4_cmdl_m12(:,13) = zin(1,:,1,1)

    END IF


    IF(p_parallel) CALL p_bcast(ch4_cmdl_m12, p_io)

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)
       
       !    Close file(s)
       
       CALL IO_close(ch4nc1)
       
    END IF

  END SUBROUTINE read_ch4_clim

  !-----------------------------------------------------------------------------------

  SUBROUTINE interpolate_ch4_clim(jglat)

    ! Interpolates monthly CH4 climatology to the current time step.

    ! *interpolate_ch4_clim* is called from *mezon* (src/socol_mezon.f90)

    USE mo_gaussgrid,     ONLY: gl_gw
    USE mo_control,       ONLY: ngl
    USE mo_socol_ghg_ods, ONLY: ch4_bcond_chem
    

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: jglat
    INTEGER :: jg

    REAL(dp) :: zch4mean, zscale

    ! calculate global mean of CH4 climatology

    zch4mean = 0._dp
    
    ! gl_gw: N->S
    ! ch4_cmdl: N->S

    do jg = 1, ngl
       zch4mean = zch4mean + gl_gw(jg)*(wgt1_chem*ch4_cmdl_m12(jg,yw1_chem) + &
            wgt2_chem*ch4_cmdl_m12(jg,yw2_chem))
    end do

    zscale = ch4_bcond_chem/zch4mean

    !!! jglat = lc% glat(krow) ! N->S latitude index

    ch4_cmdl = &
         wgt1_chem*ch4_cmdl_m12(jglat,yw1_chem) + &
         wgt2_chem*ch4_cmdl_m12(jglat,yw2_chem)

    ch4_cmdl = ch4_cmdl * zscale
    

  END SUBROUTINE interpolate_ch4_clim

  !-----------------------------------------------------------------------------------

END MODULE mo_socol_ch4

  

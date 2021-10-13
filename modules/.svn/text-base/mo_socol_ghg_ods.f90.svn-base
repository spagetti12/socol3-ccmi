MODULE mo_socol_ghg_ods

  ! Description:

  ! *mo_socol_ghg_ods* contains subroutines to prepare greenhouse gases
  ! and ODS's for radiation and chemistry module. Besides it contains 
  ! subroutines calculating the conversion factors to determine the  
  ! individual family members from the ODS families.
  !
  ! M. Schraner, ETH Zurich, November 2008

  USE mo_control,            ONLY: nlon, ngl
  USE mo_constants,          ONLY: amd, amco2, amch4, amn2o, amo3, amw
  USE mo_decomposition,      ONLY: lc => local_decomposition, &
                                   global_decomposition
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_geoloc,             ONLY: philat_2d
  USE mo_greenhouse_gases,   ONLY: ghg_no_cfc !number of ODS's for radiation
  USE mo_io
  USE mo_kind,               ONLY: dp
  USE mo_mpi,                ONLY: p_io, p_parallel, p_parallel_io, p_bcast
  USE mo_netcdf,             ONLY: FILE_INFO, io_get_vara_double
  USE mo_radiation,          ONLY: ico2, ich4, in2o, icfc, io3, ighg
  USE mo_socol_interpo,      ONLY: wgt1_chem, wgt2_chem, wgt1_rad, wgt2_rad, &
                                   m3w1_chem, m3w2_chem, m3w1_rad, m3w2_rad, &
                                   yw1_chem, yw2_chem, yw1_rad, yw2_rad
  USE mo_socol_namelist
  USE mo_socol_time_control, ONLY: l_trigchem
  USE mo_socol_tracers,      ONLY: idt_ch4, idt_n2o, idt_odscls, idt_odscll, &
                                   idt_o3, &
                                   idt_f11, idt_f12, idt_cfc113, idt_cfc114, &
                                   idt_cfc115, idt_ccl4, idt_ch3ccl3, idt_hcfc22, idt_hcfc141b, &
                                   idt_hcfc142b, idt_h1211, idt_ch3cl, idt_hcfc21, &
                                   idt_hcfc123
  
  USE mo_time_control,       ONLY: get_time_step, l_trigrad  
  USE mo_transpose,          ONLY: scatter_gp    

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  ! Parameters:  

  INTEGER, PARAMETER :: nodsyb = 5      ! number of years (before present) that
                                        ! are used to determine concentrations
                                        ! of family members in CFC-11 and CFC-12
                                        ! family (depending on age of air)
  INTEGER, PARAMETER :: nods_tot = 24   ! total number of ODS's (chemistry 
                                        ! and/or radiation)
  INTEGER, PARAMETER :: nods_tab  = 19  ! number of ODS's in dataset (for which 
                                        ! values exist and which are treated 
                                        ! by the model; all other ODS species
                                        ! are set to zero)
  INTEGER, PARAMETER, PUBLIC :: nods_chem = 20  ! number of ODS's for chemistry

  ! Radiation:

  ! Vertical levels of 3d boundary conditions:
  INTEGER               :: nflev
  REAL(dp), ALLOCATABLE :: flev(:)

  ! Monthly values of boundary conditions:
  REAL(dp)              :: co2_y_rad(0:13), ch4_y_rad(0:13), &
                           n2o_y_rad(0:13), cfc_y_rad(ghg_no_cfc,0:13)

  ! Boundary conditions at current time step:
  REAL(dp), PUBLIC      :: co2_bcond_rad, ch4_bcond_rad, n2o_bcond_rad, &
                           cfc_bcond_rad(ghg_no_cfc)

  ! Monthly values of 3d boundary conditions:
  REAL(dp), ALLOCATABLE :: ch4_3dflev_m3_rad(:,:,:,:), &
                           n2o_3dflev_m3_rad(:,:,:,:), &
                           odscls_3dflev_m3_rad(:,:,:,:), &
                           odscll_3dflev_m3_rad(:,:,:,:)

  ! 3d boundary conditions at current time step:
  REAL(dp), ALLOCATABLE :: odscls_3d_bcond_rad(:,:), &
                           odscll_3d_bcond_rad(:,:)
  REAL(dp), ALLOCATABLE, PUBLIC :: ch4_3d_bcond_rad(:,:), &
                           n2o_3d_bcond_rad(:,:), cfc_3d_bcond_rad(:,:,:)

  ! Conversion factors ODS family -> ODS family members at current time step:
  REAL(dp), ALLOCATABLE, PUBLIC :: odscl_fa2mem_rad(:,:,:) 
   
                        
  ! Chemistry:

  ! Monthly values of (lower) boundary conditions:
  REAL(dp)              :: co2_y_chem(0:13), ch4_y_chem(0:13), &
                           n2o_y_chem(0:13), odscls_y_chem(0:13), &
                           odscll_y_chem(0:13), odsbr_y_chem(0:13)

  REAL(dp)              :: odsmem_y_tab(nods_tab,0:13) 

  ! (Lower) Boundary conditions at current time step:
  REAL(dp), PUBLIC      :: co2_bcond_chem, ch4_bcond_chem, n2o_bcond_chem, &
                           odscls_bcond_chem, odscll_bcond_chem, &
                           odsbr_bcond_chem

  REAL(dp), PUBLIC      :: odsmem_bcond_chem(nods_tab)

  ! Lower boundary conditions independent on time:
  REAL(dp), PARAMETER, PUBLIC :: hcl_bcond_chem = 5.0E-12_dp, &
                           clno3_bcond_chem = 1.0E-14_dp, &
                           hbr_bcond_chem = 1.0E-13_dp
                           
  ! CO2 at current time step [mol/mol]:
  REAL(dp), PUBLIC      :: co2_chem

  ! Number of Cl/Br atoms of ODSs:
  INTEGER               :: odsncl(nods_tot), odsnbr(nods_tot)
  INTEGER, PUBLIC       :: odsnat_chem(nods_chem)

  ! Monthly values of conversion factors ODS families -> ODS members:
  REAL(dp)              :: odscl_fa2mem_yb_y(nodsyb,nods_tot,0:13), &
                           odsbr_fa2mem_yb_y(nodsyb,nods_tot,0:13), &
                           odscl_fa2mem_yb_y_const(nodsyb,nods_tot), &
                           odsbr_fa2mem_yb_y_const(nodsyb,nods_tot)

  ! Conversion factors ODS families -> ODS members at current time step:
  REAL(dp), ALLOCATABLE, PUBLIC :: ods_fa2mem_chem(:,:,:)

  ! CO2-Tracer, monthly values of lower boundary condition
  REAL, ALLOCATABLE             :: co2_y_trac(:,:,:)

  ! CO2-Tracer, lower boundary condition at current time step
  REAL(dp), ALLOCATABLE, PUBLIC :: co2_trac(:)

  ! Public entities:

  PUBLIC :: read_socol_ghg_ods_global, read_socol_ghg_3d, &
       read_socol_odsmem_global, interpolate_socol_ghg_ods_glob, &
       interpolate_socol_ghg_ods_3d, interpolate_socol_ods_fa2mem, &
       cleanup_socol_ghg_ods, ch4_chemrad, n2o_chemrad, cfc_chemrad, &
       o3_chemrad, ghg_message, interpolate_socol_co2trac
      
CONTAINS

  SUBROUTINE read_socol_ghg_ods_global(yr)

    ! Calls subroutines to read monthly values of prescribed greenhouse gases
    ! and ODS's of the current year (globally constant values).

    ! *read_socol_ghg_ods_global* is called from *read_socol_bcond_y*, 
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr

    CALL read_socol_co2_global(yr)
    CALL read_socol_ch4_global(yr)
    CALL read_socol_n2o_global(yr)
    CALL read_socol_odscl_global(yr)
    CALL read_socol_odsbr_global(yr)

    ! Set flags for ozone and water vapour coupling between radiation
    ! and chemistry:
    IF (lo3_coupl) THEN
       CALL message('',' Coupling of chemistry with radiation for O3')
       CALL message('','  -> reading O3 (radiation) from O3 (chemistry)')
       IF (io3 .NE. 0) io3 = 10
    ELSE
       IF (lchem) &
            CALL message('',' No coupling of chemistry with radiation for O3')
    ENDIF
    IF (lchem .AND. .NOT. lh2o_coupl) &
         CALL message('',' No coupling of chemistry with radiation for H2O')
    ighg = 0  ! no ECHAM5 GHG scenarios

  END SUBROUTINE read_socol_ghg_ods_global


  SUBROUTINE read_socol_ghg_3d(yr,mo)

    ! Calls subroutines to read monthly values of prescribed CH4 for 
    ! preceeding, current, and following month (3d-files from former SOCOL 
    ! simulation).

    ! *read_socol_ghg_ods_3d* is called from *read_socol_bcond_m*, 
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo

    IF (.NOT. lch4_coupl .AND. lch4_nocoupl_3d) CALL read_socol_ch4_3d(yr,mo)
    IF (.NOT. ln2o_coupl .AND. ln2o_nocoupl_3d) CALL read_socol_n2o_3d(yr,mo)
    IF (.NOT. lodscl_coupl .AND. lodscl_nocoupl_3d) &
         CALL read_socol_odscl_3d(yr,mo)

  END SUBROUTINE read_socol_ghg_3d


  SUBROUTINE read_socol_co2_global(yr)

    ! Reads monthly values of prescribed CO2 of the current year.

    ! *read_socol_co2_global* is called from *read_socol_ghg_ods_global*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr

    ! Radiation:

    CALL read_ghg_ods_monthfile('co2_y   ', 'CO2   ', '(radiation) ', &
         lco2_var_rad, yr, cyear, co2_y_rad)
    IF (ico2 .NE. 0) ico2 = 11

    ! Mixing ratio -> g/g:
    co2_y_rad(:) = co2_y_rad(:)*1.0e-06_dp*amco2/amd*co2fac

    ! Chemistry:

    IF (lchem) THEN
       CALL read_ghg_ods_monthfile('co2_y   ', 'CO2   ', '(chemistry) ', &
            lco2_var_chem, yr, cyear, co2_y_chem)
       co2_y_chem(:) = co2_y_chem(:)*1.0e-06_dp*co2fac
    ENDIF

  END SUBROUTINE read_socol_co2_global


  SUBROUTINE read_socol_ch4_global(yr)

    ! Reads monthly values of prescribed CH4 of the current year.

    ! *read_socol_ch4_global* is called from *read_socol_ghg_ods_global*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr

    ! Radiation:

    IF (lch4_coupl) THEN
       CALL message('',' Coupling of chemistry with radiation for CH4')
       CALL message('','  -> reading CH4 (radiation) from CH4 (chemistry)')
       IF (ich4 .NE. 0) ich4 = 10
    ELSE
       IF (lchem) &
            CALL message('',' No coupling of chemistry with radiation for CH4')

       IF (.NOT. lch4_nocoupl_3d) THEN
          CALL read_ghg_ods_monthfile('ch4_y   ', 'CH4   ', '(radiation) ', &
               lch4_var_rad, yr, cyear, ch4_y_rad)

          ! Mixing ratio -> g/g:
          ch4_y_rad(:) = ch4_y_rad(:)*1.0e-09_dp*amch4/amd
          IF (ich4 .NE. 0) ich4 = 11
       ENDIF
    ENDIF

    ! Chemistry:

    IF (lchem) THEN
       CALL read_ghg_ods_monthfile('ch4_y   ', 'CH4   ', '(chemistry) ', &
            lch4_var_chem, yr, cyear, ch4_y_chem)
       ch4_y_chem(:) = ch4_y_chem(:)*1.0e-09_dp
    ENDIF

  END SUBROUTINE read_socol_ch4_global


  SUBROUTINE read_socol_ch4_3d(yr,mo)

    ! Reads monthly values of prescribed CH4 for preceeding, current, and
    ! following month (3d-files from former SOCOL simulation).

    ! *read_socol_ch4_3d* is called from *read_socol_ghg_3d*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo
    
    IF (.NOT. ALLOCATED(ch4_3dflev_m3_rad)) &
         ALLOCATE(ch4_3dflev_m3_rad(lc%nproma,lc%nlev,lc%ngpblks,0:2))
    
    CALL read_ghg_3dfile('CH4   ', 'CH4   ', &
         lch4_var_rad, yr, mo, cyear, ch4_3dflev_m3_rad)
    
    ! Mixing ratio -> g/g:
    ch4_3dflev_m3_rad(:,:,:,:) = &
         ch4_3dflev_m3_rad(:,:,:,:)*amch4/amd
    IF (ich4 .NE. 0) ich4 = 12
    
  END SUBROUTINE read_socol_ch4_3d


  SUBROUTINE read_socol_n2o_global(yr)

    ! Reads monthly values of prescribed N2O of the current year.

    ! *read_socol_n2o_global* is called from *read_socol_ghg_ods_global*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr

    ! Radiation:

    IF (ln2o_coupl) THEN
       CALL message('',' Coupling of chemistry with radiation for N2O')
       CALL message('','  -> reading N2O (radiation) from N2O (chemistry)')
       IF (in2o .NE. 0) in2o = 10
    ELSE
       IF (lchem) &
            CALL message('',' No coupling of chemistry with radiation for N2O')

       IF (.NOT. ln2o_nocoupl_3d) THEN
          CALL read_ghg_ods_monthfile('n2o_y   ', 'N2O   ', '(radiation) ', &
               ln2o_var_rad, yr, cyear, n2o_y_rad)

          ! Mixing ratio -> g/g:
          n2o_y_rad(:) = n2o_y_rad(:)*1.0e-09_dp*amn2o/amd
          IF (in2o .NE. 0) in2o = 11
       ENDIF
    ENDIF

    ! Chemistry:

    IF (lchem) THEN
       CALL read_ghg_ods_monthfile('n2o_y   ', 'N2O   ', '(chemistry) ', &
            ln2o_var_chem, yr, cyear, n2o_y_chem)
       n2o_y_chem(:) = n2o_y_chem(:)*1.0e-09_dp
    ENDIF

  END SUBROUTINE read_socol_n2o_global


  SUBROUTINE read_socol_n2o_3d(yr,mo)

    ! Reads monthly values of prescribed N2O for preceeding, current, and
    ! following month (3d-files from former SOCOL simulation).

    ! *read_socol_n2o_3d* is called from *read_socol_ghg_3d*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo
    
    IF (.NOT. ALLOCATED(n2o_3dflev_m3_rad)) &
         ALLOCATE(n2o_3dflev_m3_rad(lc%nproma,lc%nlev,lc%ngpblks,0:2))
    
    CALL read_ghg_3dfile('N2O   ', 'N2O   ', &
         ln2o_var_rad, yr, mo, cyear, n2o_3dflev_m3_rad)
    
    ! Mixing ratio -> g/g:
    n2o_3dflev_m3_rad(:,:,:,:) = &
         n2o_3dflev_m3_rad(:,:,:,:)*amn2o/amd
    IF (in2o .NE. 0) in2o = 12
    
  END SUBROUTINE read_socol_n2o_3d


  SUBROUTINE read_socol_odscl_global(yr)

    ! Reads monthly values of prescribed chlorine containing ODSs of the 
    ! current year.

    ! *read_socol_odscl_global* is called from *read_socol_ghg_ods_global*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr

    ! Local variables:
!!$    REAL(dp) :: odsmem_y_tab(nods_tab,0:13)  

    ! Radiation:

    IF (lodscl_coupl) THEN
       CALL message &
            ('',' Coupling of chemistry with radiation for ODSCLS/ODSCLL')
       CALL message &
            ('','  -> reading ODSCLS/ODSLCLL (radiation) from ODSCLS/ODSCLL (chemistry)')
       IF (icfc .NE. 0) icfc = 10
    ELSE
       IF (lchem) CALL message &
            ('',' No coupling of chemistry with radiation for ODSCLS/ODSCLL')

       IF (.NOT. lodscl_nocoupl_3d) THEN
          CALL read_ghg_ods_monthfile('odsmem_y', 'ODSCL ', '(radiation) ', &
               lodscl_var_rad, yr, cyear, ods_y=odsmem_y_tab)

          ! Allocate corresponding indices of ECHAM5 radiation code:
          ! NB: In contrast to ECHAM4, only extinction cofficients of short
          !     and long lived CFCs are distinguished in ECHAM5.
          cfc_y_rad(1,:) = odsmem_y_tab(1,:) + &  ! CFC-11
                           odsmem_y_tab(8,:) + &  ! HCFC-22
                           odsmem_y_tab(16,:)+ &  ! HCFC-123
                           odsmem_y_tab(9,:) + &  ! HCFC-141b
                           odsmem_y_tab(10,:)+ &  ! HCFC-142b
                           odsmem_y_tab(6,:) + &  ! CCl4
                           odsmem_y_tab(7,:) + &  ! CH3CCl3
                           odsmem_y_tab(11,:)+ &  ! H-1211
                           odsmem_y_tab(14,:)+ &  ! CH3Cl
                           odsmem_y_tab(15,:)     ! HCFC-21
                           
          cfc_y_rad(2,:) = odsmem_y_tab(2,:) + &  ! CFC-12
                           odsmem_y_tab(3,:) + &  ! CFC-113
                           odsmem_y_tab(4,:) + &  ! CFC-114
                           odsmem_y_tab(5,:)      ! CFC-115
          ! Scale CFCs only, keep the volume mixing ratio:
          cfc_y_rad(:,:) = cfc_y_rad(:,:)*1.0e-12_dp
          IF (icfc .NE. 0) icfc = 11
       ENDIF
    ENDIF

    ! Chemistry:

    IF (lchem) THEN
       CALL read_ghg_ods_monthfile('odscls_y', 'ODSCLS', '(chemistry) ', &
            lodscl_var_chem, yr, cyear, odscls_y_chem)
       CALL read_ghg_ods_monthfile('odscll_y', 'ODSCLL', '(chemistry) ', &
            lodscl_var_chem, yr, cyear, odscll_y_chem)
       odscls_y_chem(:) = odscls_y_chem(:)*1.0e-12_dp
       odscll_y_chem(:) = odscll_y_chem(:)*1.0e-12_dp
    ENDIF

    CALL read_ghg_ods_monthfile('odsmem_y', 'ODSCL ', '(radiation) ', &
         lodscl_var_rad, yr, cyear, ods_y=odsmem_y_tab)

  END SUBROUTINE read_socol_odscl_global


  SUBROUTINE read_socol_odscl_3d(yr,mo)
    ! Reads monthly values of prescribed chlorine containing ODSs for
    ! preceeding, current, and following month (3d-files from former SOCOL 
    ! simulation).
    
    ! *read_socol_odscl_3d* is called from *read_socol_ghg_3d*.
    
    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo
    
    IF (.NOT. ALLOCATED(odscls_3dflev_m3_rad)) &
         ALLOCATE(odscls_3dflev_m3_rad(lc%nproma,lc%nlev,lc%ngpblks,0:2))
    IF (.NOT. ALLOCATED(odscll_3dflev_m3_rad)) &
         ALLOCATE(odscll_3dflev_m3_rad(lc%nproma,lc%nlev,lc%ngpblks,0:2))
    CALL read_ghg_3dfile('ODSCLS', 'ODSCLS', &
         lodscl_var_rad, yr, mo, cyear, odscls_3dflev_m3_rad)
    CALL read_ghg_3dfile('ODSCLL', 'ODSCLS', &
         lodscl_var_rad, yr, mo, cyear, odscll_3dflev_m3_rad)
    
    odscls_3dflev_m3_rad(:,:,:,:) = odscls_3dflev_m3_rad(:,:,:,:)
    odscll_3dflev_m3_rad(:,:,:,:) = odscll_3dflev_m3_rad(:,:,:,:)
    IF (icfc .NE. 0) icfc = 12

  END SUBROUTINE read_socol_odscl_3d


  SUBROUTINE read_socol_odsbr_global(yr)

    ! Reads monthly values of prescribed bromine containing ODS's of the 
    ! current year.

    ! *read_socol_odsbr_global* is called from *read_socol_ghg_ods_global*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr

    ! Chemistry:

    IF (lchem) CALL read_ghg_ods_monthfile('odsbr_y ', &
         'ODSBR ', '(chemistry) ', lodsbr_var_chem, yr, cyear, odsbr_y_chem)
    odsbr_y_chem(:) = odsbr_y_chem(:)*1.0e-12_dp

  END SUBROUTINE read_socol_odsbr_global


  SUBROUTINE interpolate_socol_ghg_ods_glob

    ! Interpolates globally constant monthly GHG/ODS boundary condition 
    ! fields to the current time step.
    !
    ! *interpolate_socol_ghg_ods_glob* is called from 
    ! *interpolate_socol_bcond_global*, src/socol_boundary_conditions.90.

    ! Radiation:
    IF (l_trigrad) THEN
       co2_bcond_rad = wgt1_rad*co2_y_rad(yw1_rad) + &
            wgt2_rad*co2_y_rad(yw2_rad)
       IF (.NOT. lch4_coupl .AND. .NOT. lch4_nocoupl_3d) &
            ch4_bcond_rad = wgt1_rad*ch4_y_rad(yw1_rad) + &
            wgt2_rad*ch4_y_rad(yw2_rad)
       IF (.NOT. ln2o_coupl .AND. .NOT. ln2o_nocoupl_3d) &
            n2o_bcond_rad = wgt1_rad*n2o_y_rad(yw1_rad) + &
            wgt2_rad*n2o_y_rad(yw2_rad)
       IF (.NOT. lodscl_coupl .AND. .NOT. lodscl_nocoupl_3d) &
            cfc_bcond_rad(:) = wgt1_rad*cfc_y_rad(:,yw1_rad) + &
            wgt2_rad*cfc_y_rad(:,yw2_rad)
    ENDIF

    ! Chemistry:
    IF (l_trigchem) co2_bcond_chem = wgt1_chem*co2_y_chem(yw1_chem) + &
         wgt2_chem*co2_y_chem(yw2_chem)
    ch4_bcond_chem = wgt1_chem*ch4_y_chem(yw1_chem) + &
         wgt2_chem*ch4_y_chem(yw2_chem)
    n2o_bcond_chem = wgt1_chem*n2o_y_chem(yw1_chem) + &
         wgt2_chem*n2o_y_chem(yw2_chem)
    odscls_bcond_chem = wgt1_chem*odscls_y_chem(yw1_chem) + &
         wgt2_chem*odscls_y_chem(yw2_chem)
    odscll_bcond_chem = wgt1_chem*odscll_y_chem(yw1_chem) + &
         wgt2_chem*odscll_y_chem(yw2_chem)
    odsbr_bcond_chem  = wgt1_chem*odsbr_y_chem(yw1_chem) + &
         wgt2_chem*odsbr_y_chem(yw2_chem)

    odsmem_bcond_chem(:) = wgt1_chem*odsmem_y_tab(:,yw1_chem) + &
         wgt2_chem*odsmem_y_tab(:,yw2_chem)
    odsmem_bcond_chem(:) = odsmem_bcond_chem(:)*1.0e-12_dp ! volume mixing ratio

  END SUBROUTINE interpolate_socol_ghg_ods_glob


  SUBROUTINE interpolate_socol_ghg_ods_3d(krow, kproma, kbdim, klev, p)

    ! Interpolates monthly GHG/ODS boundary condition 3d-fields to the current
    ! time step.
    !
    ! *interpolate_socol_ghg_ods_3d* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim, klev
    REAL(dp), INTENT(in) :: p(kbdim,klev)

    ! CH4:
    IF (l_trigrad .AND. lch4_nocoupl_3d) THEN
       IF (.NOT. ALLOCATED(ch4_3d_bcond_rad)) & 
            ALLOCATE(ch4_3d_bcond_rad(kbdim,klev))  
       CALL interpolate_3dflev(krow, kproma, kbdim, klev, p, &
            ch4_3dflev_m3_rad(:,:,krow,:), ch4_3d_bcond_rad)
    ENDIF

    ! N2O:
    IF (l_trigrad .AND. ln2o_nocoupl_3d) THEN
       IF (.NOT. ALLOCATED(n2o_3d_bcond_rad)) &
            ALLOCATE(n2o_3d_bcond_rad(kbdim,klev))  
       CALL interpolate_3dflev(krow, kproma, kbdim, klev, p, &
            n2o_3dflev_m3_rad(:,:,krow,:), n2o_3d_bcond_rad)
    ENDIF

    ! CFCs:
    IF (lodscl_nocoupl_3d) THEN
       IF (.NOT. ALLOCATED(odscls_3d_bcond_rad)) &
            ALLOCATE(odscls_3d_bcond_rad(kbdim,klev))
       IF (.NOT. ALLOCATED(odscll_3d_bcond_rad)) &
            ALLOCATE(odscll_3d_bcond_rad(kbdim,klev))
       IF (.NOT. ALLOCATED(cfc_3d_bcond_rad)) THEN
          ALLOCATE(cfc_3d_bcond_rad(kbdim,klev,ghg_no_cfc))
          CALL interpolate_3dflev(krow, kproma, kbdim, klev, p, &
               odscls_3dflev_m3_rad(:,:,krow,:), odscls_3d_bcond_rad)
          CALL interpolate_3dflev(krow, kproma, kbdim, klev, p, &
               odscll_3dflev_m3_rad(:,:,krow,:), odscll_3d_bcond_rad)
          ! Allocate corresponding indices of ECHAM5 radiation code:
          cfc_3d_bcond_rad(1:kproma,:,1) = odscls_3d_bcond_rad(1:kproma,:) * &
               (odscl_fa2mem_rad(1:kproma,:,1) + &  ! CFC-11
                odscl_fa2mem_rad(1:kproma,:,6) + &  ! HCFC-22
                odscl_fa2mem_rad(1:kproma,:,7) + &  ! HCFC-123
                odscl_fa2mem_rad(1:kproma,:,11)+ &  ! HCFC-141b
                odscl_fa2mem_rad(1:kproma,:,12)+ &  ! HCFC-142b
                odscl_fa2mem_rad(1:kproma,:,15)+ &  ! CCl4
                odscl_fa2mem_rad(1:kproma,:,16)+ &  ! CH3CCl3
                odscl_fa2mem_rad(1:kproma,:,17)+ &  ! H-1211
                odscl_fa2mem_rad(1:kproma,:,20)+ &  ! CH3Cl
                odscl_fa2mem_rad(1:kproma,:,21))    ! HCFC-21
          cfc_3d_bcond_rad(1:kproma,:,2) = odscls_3d_bcond_rad(1:kproma,:) * &
               (odscl_fa2mem_rad(1:kproma,:,2) + &  ! CFC-12
                odscl_fa2mem_rad(1:kproma,:,3) + &  ! CFC-113
                odscl_fa2mem_rad(1:kproma,:,4) + &  ! CFC-114
                odscl_fa2mem_rad(1:kproma,:,5))     ! CFC-115
       ENDIF
    ENDIF

  END SUBROUTINE interpolate_socol_ghg_ods_3d

  SUBROUTINE interpolate_socol_co2trac(krow, kproma, kbdim)

    ! Interpolates monthly lower boundary condition of CO2-Tracer to the current
    ! time step.
    !
    ! *interpolate_socol_co2trac* is called from 
    ! *interpolate_socol_bcond*, src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim

    !CO2-Tracer

    IF (.NOT. ALLOCATED(co2_trac)) &
         ALLOCATE(co2_trac(kbdim))

    co2_trac(1:kproma) = wgt1_chem*co2_y_trac(1:kproma,krow,yw1_chem) + &
         wgt2_chem*co2_y_trac(1:kproma,krow,yw2_chem)

  END SUBROUTINE interpolate_socol_co2trac

  SUBROUTINE read_socol_odsmem_global(yr)

    ! Reads monthly values of individual members of ODSCLS, ODSCLL and ODSBR 
    ! families of the current year and the nodsyb preceeding years as well as
    ! the chemical lifetimes and the number of Cl/Br-atoms. Besides calls
    ! subroutine *calculate_ods_fa2mem* to prepare conversion factors 
    ! odscl_fa2mem/ odsbr_fa2mem used for calculation the individual ODS 
    ! members from their family.

    ! *read_socol_odsmem_global* is called from
    ! *read_socol_bcond_y*, src/socol_boundary_conditions.f90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr

    ! Local variables:
    LOGICAL :: lvar
    INTEGER :: jj, jj1, jj2
    INTEGER :: n_cl(nods_tab), n_br(nods_tab)
    REAL(dp) :: ltods(nods_tab), odsmem_yb(nods_tab,-nodsyb*12+1:13)
    REAL(dp) :: odslifetime(nods_tot)

    jj1=2
    jj2=1

    IF (lodscl_var_rad .OR. lodscl_var_chem .OR. lodsbr_var_chem) jj1=1
    IF (.NOT. lodscl_var_rad .OR. .NOT. lodscl_var_chem .OR. &
         .NOT. lodsbr_var_chem) jj2=2

    DO jj = jj1, jj2

       IF (jj .EQ. 1) lvar = .TRUE.
       IF (jj .EQ. 2) lvar = .FALSE.

       ! Read values of ODS family members odsmem_yb [pptv] for all months back
       ! to nodsyb years before present as well as number of Cl/Br atoms
       ! n_cl, n_br and chemical lifetimes ltods:
       CALL read_ghg_ods_monthfile('odsmem_y', 'ODS   ', '(fa members)', &
            lvar, yr, cyear, n_cl=n_cl, n_br=n_br, ltods=ltods, &
            odsmem_yb=odsmem_yb)
       
       ! Save number of Cl-atoms (new order of indices):
       odsncl(1)  = n_cl(1)   ! CFC-11
       odsncl(2)  = n_cl(2)   ! CFC-12
       odsncl(3)  = n_cl(3)   ! CFC-113
       odsncl(4)  = n_cl(4)   ! CFC-114
       odsncl(5)  = n_cl(5)   ! CFC-115
       odsncl(6)  = n_cl(8)   ! HCFC-22
       odsncl(7)  = n_cl(16)  ! HCFC-123
       odsncl(8)  = 1         ! HCFC-124
       odsncl(9)  = 0         ! HCFC-125
       odsncl(10) = 0         ! HCFC-134a
       odsncl(11) = n_cl(9)   ! HCFC-141b
       odsncl(12) = n_cl(10)  ! HCFC-142b
       odsncl(13) = 4         ! HCFC-143a
       odsncl(14) = 3         ! HCFC-152a
       odsncl(15) = n_cl(6)   ! CCl4
       odsncl(16) = n_cl(7)   ! CH3CCl3
       odsncl(17) = n_cl(11)  ! H-1211
       odsncl(18) = n_cl(12)  ! H-1301
       odsncl(19) = n_cl(13)  ! CH3Br
       odsncl(20) = n_cl(14)  ! CH3Cl
       odsncl(21) = n_cl(15)  ! HCFC-21
       odsncl(22) = n_cl(17)  ! H-2402
       odsncl(23) = n_cl(18)  ! CHBr3
       odsncl(24) = n_cl(19)  ! CH2Br2

       ! Save number of Br-atoms (new order of indices):      
       odsnbr(1)  = n_br(1)   ! CFC-11
       odsnbr(2)  = n_br(2)   ! CFC-12
       odsnbr(3)  = n_br(3)   ! CFC-113
       odsnbr(4)  = n_br(4)   ! CFC-114
       odsnbr(5)  = n_br(5)   ! CFC-115
       odsnbr(6)  = n_br(8)   ! HCFC-22
       odsnbr(7)  = n_br(16)  ! HCFC-123
       odsnbr(8)  = 0         ! HCFC-124
       odsnbr(9)  = 0         ! HCFC-125
       odsnbr(10) = 0         ! HCFC-134a
       odsnbr(11) = n_br(9)   ! HCFC-141b
       odsnbr(12) = n_br(10)  ! HCFC-142b
       odsnbr(13) = 0         ! HCFC-143a
       odsnbr(14) = 0         ! HCFC-152a
       odsnbr(15) = n_br(6)   ! CCl4
       odsnbr(16) = n_br(7)   ! CH3CCl3
       odsnbr(17) = n_br(11)  ! H-1211
       odsnbr(18) = n_br(12)  ! H-1301
       odsnbr(19) = n_br(13)  ! CH3Br
       odsnbr(20) = n_br(14)  ! CH3Cl
       odsnbr(21) = n_br(15)  ! HCFC-21
       odsnbr(22) = n_br(17)  ! H-2402
       odsnbr(23) = n_br(18)  ! CHBr3
       odsnbr(24) = n_br(19)  ! CH2Br2
 
       ! Save ODS lifetime (new order of indices):      
       odslifetime(1)  = ltods(1)   ! CFC-11
       odslifetime(2)  = ltods(2)   ! CFC-12
       odslifetime(3)  = ltods(3)   ! CFC-113
       odslifetime(4)  = ltods(4)   ! CFC-114
       odslifetime(5)  = ltods(5)   ! CFC-115
       odslifetime(6)  = ltods(8)   ! HCFC-22
       odslifetime(7)  = ltods(16)  ! HCFC-123
       odslifetime(8)  = -9999._dp  ! HCFC-124
       odslifetime(9)  = -9999._dp  ! HCFC-125
       odslifetime(10) = -9999._dp  ! HCFC-134a
       odslifetime(11) = ltods(9)   ! HCFC-141b
       odslifetime(12) = ltods(10)  ! HCFC-142b
       odslifetime(13) = -9999._dp  ! HCFC-143a
       odslifetime(14) = -9999._dp  ! HCFC-152a
       odslifetime(15) = ltods(6)   ! CCl4
       odslifetime(16) = ltods(7)   ! CH3CCl3
       odslifetime(17) = ltods(11)  ! H-1211
       odslifetime(18) = ltods(12)  ! H-1301
       odslifetime(19) = ltods(13)  ! CH3Br
       odslifetime(20) = ltods(14)  ! CH3Cl
       odslifetime(21) = ltods(15)  ! HCFC-21
       odslifetime(22) = ltods(17)  ! H-2402
       odslifetime(23) = ltods(18)  ! CHBr3
       odslifetime(24) = ltods(19)  ! CH2Br2

       ! Calculation of conversion factors odscl_fa2mem_yb_y, odsbr_fa2mem_yb_y
       ! for ages of air from present back to nodsyb years before present. 
       ! If ODS's are kept constant in time (lvar=.FALSE.), 
       ! odscl_fa2mem_yb_y_const and odsbr_fa2mem_yb_y_const are
       ! calculated instead. 
       CALL calculate_ods_fa2mem(odsmem_yb, odslifetime, 'ODSCLS', lvar)
       CALL calculate_ods_fa2mem(odsmem_yb, odslifetime, 'ODSCLL', lvar)
       CALL calculate_ods_fa2mem(odsmem_yb, odslifetime, 'ODSBR', lvar)

    ENDDO

    ! Save odsncl, odsnbr with order of indices as in MEZON (-> odsnat_chem):
    IF (lchem) THEN
       odsnat_chem(1)  = odsncl(1)     !CFC11
       odsnat_chem(2)  = odsncl(2)     !CFC12
       odsnat_chem(4)  = odsncl(3)     !CFC113
       odsnat_chem(5)  = odsncl(4)     !CFC114
       odsnat_chem(6)  = odsncl(5)     !CFC115
       odsnat_chem(7)  = odsncl(15)    !CCL4
       odsnat_chem(8)  = odsncl(16)    !CH3CCL3
       odsnat_chem(9)  = odsncl(6)     !HCFC22
       odsnat_chem(10) = odsncl(11)    !HCFC141B
       odsnat_chem(11) = odsncl(12)    !HCFC142B
       odsnat_chem(14) = odsncl(20)    !CH3CL
       odsnat_chem(15) = odsncl(21)    !HCFC21
       odsnat_chem(16) = odsncl(7)     !HCFC123
          
       odsnat_chem(3)  = odsnbr(18)    !H1301
       odsnat_chem(12) = odsnbr(17)    !H1211
       odsnat_chem(13) = odsnbr(19)    !CH3Br
       odsnat_chem(17) = odsnbr(22)    !H2402
       odsnat_chem(18) = odsnbr(23)    !CHBR3
       odsnat_chem(19) = odsnbr(24)    !CH2BR2
       odsnat_chem(20) = odsncl(17)    !H1211 (Cl)
    ENDIF
       
  END SUBROUTINE read_socol_odsmem_global


  SUBROUTINE calculate_ods_fa2mem(odsmem_yb, odslifetime, ods_name, lvar)

    ! Calculates conversion factors odscl_fa2mem_yb_y, odsbr_fa2mem_yb_y
    ! for ages of air from present back to nodsyb years before present.
    ! If ODS's are kept constant in time (lvar=.FALSE.), 
    ! odscl_fa2mem_yb_y_const and odsbr_fa2mem_yb_y_const are
    ! calculated instead.

    ! *calculate_ods_fa2mem* is called from *read_socol_odsmem*. 

    ! Subroutine arguments:
    CHARACTER (*), INTENT(in) :: ods_name
    LOGICAL, INTENT(in)       :: lvar
    REAL(dp), INTENT(in)      :: odsmem_yb(nods_tab,-nodsyb*12+1:13), &
                                 odslifetime(nods_tot)

    ! Local variables:
    LOGICAL  :: lpo, lex0, lex1, lex2
    INTEGER  :: i, i1, i2, ii, iii, ipo, iipo, jc, m
    REAL(dp) :: odsexponatsu, mai
    REAL(dp) :: odsmem_me(nods_tab)
    INTEGER  :: odsnat(nods_tot) 
    REAL(dp) :: odsmem(nods_tot), odsexpo(nods_tot), odspr(nods_tot), &
                ods_fa2mem_yb_y(nodsyb,nods_tot,0:13)
                
    ipo = -1
    iipo = -1
    lpo = .FALSE.

    IF (lvar) THEN
       i1 = 0
       i2 = 13
    ELSE
       i1 = 12
       i2 = 12
    ENDIF

    DO ii = 1, nodsyb  !loop over years
       DO i = i1, i2      !loop over months

          IF (lvar) THEN
             iii = ii
          ELSE
             iii = 1
          ENDIF

          ! Mean of ODS family members [pptv] over preceeding 12 months:
          odsmem_me(:) = 0._dp
          DO jc =1, nods_tab
             DO m = 1, 12
                odsmem_me(jc) = odsmem_me(jc)+odsmem_yb(jc,-iii*12+i+m)
             ENDDO
          ENDDO
          odsmem_me(:) = odsmem_me(:)/12._dp

          ! Save ODS family members [pptv] for month/ year of the loop (new 
          ! order of indices):
          SELECT CASE (TRIM(ods_name))

          CASE ('ODSCLS')
             odsmem(1)  = odsmem_me(1)    ! CFC-11
             odsmem(2)  = 0._dp           ! CFC-12
             odsmem(3)  = 0._dp           ! CFC-113
             odsmem(4)  = 0._dp           ! CFC-114
             odsmem(5)  = 0._dp           ! CFC-115
             odsmem(6)  = odsmem_me(8)    ! HCFC-22
             odsmem(7)  = odsmem_me(16)   ! HCFC-123
             odsmem(8)  = 0._dp           ! HCFC-124
             odsmem(9)  = 0._dp           ! HCFC-125
             odsmem(10) = 0._dp           ! HCFC-134a
             odsmem(11) = odsmem_me(9)    ! HCFC-141b
             odsmem(12) = odsmem_me(10)   ! HCFC-142b
             odsmem(13) = 0._dp           ! HCFC-143a
             odsmem(14) = 0._dp           ! HCFC-152a
             odsmem(15) = odsmem_me(6)    ! CCl4
             odsmem(16) = odsmem_me(7)    ! CH3CCl3
             odsmem(17) = odsmem_me(11)   ! H-1211
             odsmem(18) = 0._dp           ! H-1301
             odsmem(19) = 0._dp           ! CH3Br
             odsmem(20) = odsmem_me(14)   ! CH3Cl
             odsmem(21) = odsmem_me(15)   ! HCFC-21
             odsmem(22) = 0._dp           ! H-2402
             odsmem(23) = 0._dp           ! CHBr3
             odsmem(24) = 0._dp           ! CH2Br2

             odsnat(:) = odsncl(:)        ! number of Cl atoms
             
          CASE ('ODSCLL')             
             odsmem(1)  = 0._dp           ! CFC-11
             odsmem(2)  = odsmem_me(2)    ! CFC-12
             odsmem(3)  = odsmem_me(3)    ! CFC-113
             odsmem(4)  = odsmem_me(4)    ! CFC-114
             odsmem(5)  = odsmem_me(5)    ! CFC-115
             odsmem(6)  = 0._dp           ! HCFC-22
             odsmem(7)  = 0._dp           ! HCFC-123
             odsmem(8)  = 0._dp           ! HCFC-124
             odsmem(9)  = 0._dp           ! HCFC-125
             odsmem(10) = 0._dp           ! HCFC-134a
             odsmem(11) = 0._dp           ! HCFC-141b
             odsmem(12) = 0._dp           ! HCFC-142b
             odsmem(13) = 0._dp           ! HCFC-143a
             odsmem(14) = 0._dp           ! HCFC-152a
             odsmem(15) = 0._dp           ! CCl4
             odsmem(16) = 0._dp           ! CH3CCl3
             odsmem(17) = 0._dp           ! H-1211
             odsmem(18) = 0._dp           ! H-1301
             odsmem(19) = 0._dp           ! CH3Br
             odsmem(20) = 0._dp           ! CH3Cl
             odsmem(21) = 0._dp           ! HCFC-21
             odsmem(22) = 0._dp           ! H-2402
             odsmem(23) = 0._dp           ! CHBr3
             odsmem(24) = 0._dp           ! CH2Br2

             odsnat(:) = odsncl(:)        ! number of Cl atoms
             
          CASE ('ODSBR')     
             odsmem(1)  = 0._dp           ! CFC-11
             odsmem(2)  = 0._dp           ! CFC-12
             odsmem(3)  = 0._dp           ! CFC-113
             odsmem(4)  = 0._dp           ! CFC-114
             odsmem(5)  = 0._dp           ! CFC-115
             odsmem(6)  = 0._dp           ! HCFC-22
             odsmem(7)  = 0._dp           ! HCFC-123
             odsmem(8)  = 0._dp           ! HCFC-124
             odsmem(9)  = 0._dp           ! HCFC-125
             odsmem(10) = 0._dp           ! HCFC-134a
             odsmem(11) = 0._dp           ! HCFC-141b
             odsmem(12) = 0._dp           ! HCFC-142b
             odsmem(13) = 0._dp           ! HCFC-143a
             odsmem(14) = 0._dp           ! HCFC-152a
             odsmem(15) = 0._dp           ! CCl4
             odsmem(16) = 0._dp           ! CH3CCl3
             odsmem(17) = odsmem_me(11)   ! H-1211
             odsmem(18) = odsmem_me(12)   ! H-1301
             odsmem(19) = odsmem_me(13)   ! CH3Br
             odsmem(20) = 0._dp           ! CH3Cl
             odsmem(21) = 0._dp           ! HCFC-21
             odsmem(22) = odsmem_me(17)   ! H-2402
             odsmem(23) = odsmem_me(18)   ! CHBr3
             odsmem(24) = odsmem_me(19)   ! CH2Br2

             odsnat(:) = odsnbr(:)        !number of Br atoms

          END SELECT

          mai = REAL(ii,dp) - 0.5_dp

          ! Calculate 1/n(i)*P(i) (n(i) number of Cl/Br atoms, P(i) as in 
          ! Schraner et al. 2008):
          DO jc = 1, nods_tot
             IF (odslifetime(jc) .LT. -9998._dp .OR. &
                  odsmem(jc) .LT. EPSILON(1.0_dp)) THEN  !empty values
                odsexpo(jc) = 0._dp
             ELSE
                odsexpo(jc) = odsmem(jc)*EXP(-mai/odslifetime(jc))
             ENDIF
          ENDDO

          odsexponatsu = 0._dp
          DO jc = 1, nods_tot
             odsexponatsu = odsexponatsu + odsexpo(jc)*REAL(odsnat(jc),dp)
          ENDDO

          IF (odsexponatsu .GE. EPSILON(1.0_dp)) THEN
             
             IF (iipo .EQ. -1 .AND. (i .GE. i1+1 .OR. ii .GE. 2)) THEN
                ! ODSCLS/L/BR > 0 for current year/month. But in all year(s)/
                ! month(s) read in before, ODSCLS/L/BR was = 0.
                ! -> Take the current value of ods_fa2mem_yb_y to calculate
                ! the distribution between family and family members for
                ! all year(s)/month(s) read in before.
                lpo = .TRUE.
             END IF
             
             ipo = i
             iipo = ii
                
             DO jc = 1, nods_tot
                ods_fa2mem_yb_y(ii,jc,i) = odsexpo(jc)/odsexponatsu
             ENDDO

             IF (lpo) THEN
                ! ODSCLS/L/BR > 0 for current year/month. But in all year(s)/
                ! month(s) read in before, ODSCLS/L/BR was = 0.
                ! -> Take the current value of ods_fa2mem_yb_y to calculate
                ! the distribution between family and family members for
                ! all year(s)/month(s) read in before.
                DO jc = 1, nods_tot
                   IF (ii .GE. 2) THEN 
                      ods_fa2mem_yb_y(0:ii-1,jc,:) = ods_fa2mem_yb_y(ii,jc,i)
                   END IF
                   IF (i .GE. 1) THEN
                      ods_fa2mem_yb_y(ii,jc,0:i-1) = ods_fa2mem_yb_y(ii,jc,i)
                   END IF
                   lpo = .FALSE.
                END DO
             ENDIF

          ELSE
             IF (iipo .EQ. -1) THEN
                ods_fa2mem_yb_y(ii,:,i) = 0._dp
             ELSE
                ! ODSCLS/L/BR = 0 for current month of year. But there are
                ! year(s)/month(s), where ODSCLS/L/BR > 0.
                ! -> Take this value of ods_fa2mem_yb_y to calculate
                ! the distribution between family and family members for 
                ! current month of year. 
                ods_fa2mem_yb_y(ii,:,i) = ods_fa2mem_yb_y(iipo,:,ipo)
             ENDIF
          ENDIF

       ENDDO
    ENDDO

    SELECT CASE (TRIM(ods_name))

    CASE ('ODSCLS')
       IF (lvar) THEN
          odscl_fa2mem_yb_y(:,:,:) = ods_fa2mem_yb_y(:,:,:)
       ELSE
          odscl_fa2mem_yb_y_const(:,:) = ods_fa2mem_yb_y(:,:,12)
       ENDIF
    CASE ('ODSCLL')
       IF (lvar) THEN
          DO ii = 1, nodsyb        
             DO i = 0, 13  
                DO jc = 1, nods_tot
                   odscl_fa2mem_yb_y(ii,jc,i) = MAX(ods_fa2mem_yb_y(ii,jc,i), &
                        odscl_fa2mem_yb_y(ii,jc,i))
                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO ii = 1, nodsyb        
             DO jc = 1, nods_tot
                odscl_fa2mem_yb_y_const(ii,jc) = &
                     MAX(ods_fa2mem_yb_y(ii,jc,12), &
                     odscl_fa2mem_yb_y_const(ii,jc))
             ENDDO
          ENDDO
       ENDIF
    CASE ('ODSBR')
       IF (lvar) THEN
          odsbr_fa2mem_yb_y(:,:,:) = ods_fa2mem_yb_y(:,:,:)
       ELSE
          odsbr_fa2mem_yb_y_const(:,:) = ods_fa2mem_yb_y(:,:,12)
       ENDIF
    END SELECT

  END SUBROUTINE calculate_ods_fa2mem


  SUBROUTINE interpolate_socol_ods_fa2mem(krow, kproma, kbdim, klev, p, ptropo)

    ! Calculates 3d conversion factors odscl_fa2mem and odsbr_fa2mem for 
    ! current time step:
    ! 1. Interpolation of odscl_fa2mem_yb_y and odsbr_fa2mem_yb_y to current
    !    time step (-> odscl_fa2mem_yb and odsbr_fa2mem_yb).
    ! 2. Allocation of indices of odscl_fa2mem_yb and odsbr_fa2mem_yb 
    !    according to age of air (Manzini and Feichter, 1999)
    !    (-> odscl_fa2mem, odsbr_fa2mem).
    ! 3. Order of indices as in MEZON (-> ods_fa2mem_chem).
    !
    ! *interpolate_socol_ods_fa2mem* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.f90.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim, klev
    REAL(dp), INTENT(in) :: p(kbdim,klev), ptropo(kbdim)

    ! Local variables:
    INTEGER  :: jl, jk, jj, jj1, jj2
    REAL(dp) :: lat
    REAL(dp) :: odscl_fa2mem_yb(nodsyb,nods_tot), &
                odsbr_fa2mem_yb(nodsyb,nods_tot)
    REAL(dp) :: odscl_fa2mem(kbdim,klev,nods_tot), &
                odsbr_fa2mem(kbdim,klev,nods_tot), & 
                odscl_fa2memconst(kbdim,klev,nods_tot), &
                odsbr_fa2memconst(kbdim,klev,nods_tot)

    ! Allocate memory:
    IF (.NOT. ALLOCATED(odscl_fa2mem_rad)) &
         ALLOCATE(odscl_fa2mem_rad(kbdim,klev,nods_tot))
    IF (lchem) THEN
       IF (.NOT. ALLOCATED(ods_fa2mem_chem)) &
         ALLOCATE(ods_fa2mem_chem(kbdim,klev,nods_chem))
    ENDIF

    jj1=2
    jj2=1

    IF (lodscl_var_rad .OR. lodscl_var_chem .OR. lodsbr_var_chem) jj1=1
    IF (.NOT. lodscl_var_rad .OR. .NOT. lodscl_var_chem .OR. &
         .NOT. lodsbr_var_chem) jj2=2

    DO jj = jj2, jj1, -1

       ! Interpolation to current time step:
       IF (jj .EQ. 1) THEN
          odscl_fa2mem_yb(:,:) = (wgt1_chem*odscl_fa2mem_yb_y(:,:,yw1_chem) + &
               wgt2_chem*odscl_fa2mem_yb_y(:,:,yw2_chem))
          odsbr_fa2mem_yb(:,:) = (wgt1_chem*odsbr_fa2mem_yb_y(:,:,yw1_chem) + &
               wgt2_chem*odsbr_fa2mem_yb_y(:,:,yw2_chem))
       ENDIF
       IF (jj .EQ. 2) THEN
          odscl_fa2mem_yb(:,:) = odscl_fa2mem_yb_y_const(:,:)
          odsbr_fa2mem_yb(:,:) = odsbr_fa2mem_yb_y_const(:,:)
       ENDIF

       ! Determination of conversion factor according to pressure level and 
       ! latitude of grid box (parameterisation of age of air by M. Schraner
       ! after Manzini and Feichter, JGR, Vol. 104, 1999):
       DO jl = 1, kproma
          
          lat = philat_2d(jl,krow)  !latitude in degrees

          IF (ABS(lat) .LE. 30._dp) THEN
             ! 1.  0° <= latitude <= 30°

             DO jk =1, klev
              
                IF (p(jl,jk) .GE. 5000._dp) THEN
                   ! 1.a  0 <= age of air < 1 year 
                   odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(1,:)
                   odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(1,:)
                ELSE
                   IF (p(jl,jk) .GE. 3000._dp) THEN
                      ! 1.b  1 <= age of air < 2 years 
                      odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(2,:)
                      odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(2,:)
                   ELSE
                      IF (p(jl,jk) .GE. 1500._dp) THEN
                         ! 1.c  2 <= age of air < 3 years
                         odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(3,:)
                         odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(3,:)
                      ELSE
                         IF (p(jl,jk) .GE. 50._dp) THEN
                            ! 1.c  3 <= age of air < 4 years 
                            odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(4,:)
                            odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(4,:)
                         ELSE
                            ! 1.d  age of air >= 4 years
                            odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(5,:)
                            odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(5,:)
                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF

             ENDDO

          ELSE

             IF (ABS(lat) .LE. 60._dp) THEN
                ! 2.  30° < latitude <= 60°

                DO jk =1, klev
              
                   IF (p(jl,jk) .GE. ptropo(jl)) THEN
                      ! 2.a  0 <= age of air < 1 year 
                      odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(1,:)
                      odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(1,:)
                   ELSE
                      IF (p(jl,jk) .GE. 6000._dp) THEN
                         ! 2.b  1 <= age of air < 2 years 
                         odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(2,:)
                         odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(2,:)
                      ELSE
                         IF (p(jl,jk) .GE. 4000._dp) THEN
                            ! 2.c  2 <= age of air < 3 years
                            odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(3,:)
                            odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(3,:)
                         ELSE
                            IF (p(jl,jk) .GE. 1000._dp) THEN
                               ! 2.c  3 <= age of air < 4 years 
                               odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(4,:)
                               odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(4,:)
                            ELSE
                               ! 2.d  age of air >= 4 years
                               odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(5,:)
                               odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(5,:)
                            ENDIF
                         ENDIF
                      ENDIF
                   ENDIF
                   
                ENDDO

             ELSE

                ! 3.  60° < latitude <= 90°

                DO jk =1, klev
              
                   IF (p(jl,jk) .GE. ptropo(jl)) THEN
                      ! 3.a  0 <= age of air < 1 year 
                      odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(1,:)
                      odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(1,:)
                   ELSE
                      IF (p(jl,jk) .GE. 10000._dp) THEN
                         ! 3.b  1 <= age of air < 2 years 
                         odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(2,:)
                         odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(2,:)
                      ELSE
                         IF (p(jl,jk) .GE. 5000._dp) THEN
                            ! 3.c  2 <= age of air < 3 years
                            odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(3,:)
                            odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(3,:)
                         ELSE
                            IF (p(jl,jk) .GE. 2000._dp) THEN
                               ! 3.c  3 <= age of air < 4 years 
                               odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(4,:)
                               odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(4,:)
                            ELSE
                               ! 3.d  age of air >= 4 years
                               odscl_fa2mem(jl,jk,:) = odscl_fa2mem_yb(5,:)
                               odsbr_fa2mem(jl,jk,:) = odsbr_fa2mem_yb(5,:)
                            ENDIF
                         ENDIF
                      ENDIF
                   ENDIF
                   
                ENDDO

             ENDIF
          ENDIF

       ENDDO

       IF (jj .EQ. 2) THEN
          odscl_fa2memconst(1:kproma,:,:) = odscl_fa2mem(1:kproma,:,:)
          odsbr_fa2memconst(1:kproma,:,:) = odsbr_fa2mem(1:kproma,:,:)
       ENDIF

    ENDDO

    IF (l_trigrad) THEN
       IF (lodscl_var_rad) THEN
          odscl_fa2mem_rad(1:kproma,:,:) = odscl_fa2mem(1:kproma,:,:)
       ELSE
          odscl_fa2mem_rad(1:kproma,:,:) = odscl_fa2memconst(1:kproma,:,:)
       ENDIF
    ENDIF

    ! Save odscl_fa2mem and odsbr_fa2mem with order of indices as
    ! in MEZON (-> ods_fa2mem_chem):
    IF (l_trigchem) THEN
       IF (lodscl_var_chem) THEN 
          ods_fa2mem_chem(1:kproma,:,1)  = &
               odscl_fa2mem(1:kproma,:,1)       !CFC11
          ods_fa2mem_chem(1:kproma,:,2)  = &
               odscl_fa2mem(1:kproma,:,2)       !CFC12
          ods_fa2mem_chem(1:kproma,:,4)  = &
               odscl_fa2mem(1:kproma,:,3)       !CFC113
          ods_fa2mem_chem(1:kproma,:,5)  = &
               odscl_fa2mem(1:kproma,:,4)       !CFC114
          ods_fa2mem_chem(1:kproma,:,6)  = &
               odscl_fa2mem(1:kproma,:,5)       !CFC115
          ods_fa2mem_chem(1:kproma,:,7)  = &
               odscl_fa2mem(1:kproma,:,15)      !CCL4
          ods_fa2mem_chem(1:kproma,:,8)  = &
               odscl_fa2mem(1:kproma,:,16)      !CH3CCL3
          ods_fa2mem_chem(1:kproma,:,9)  = &
               odscl_fa2mem(1:kproma,:,6)       !HCFC22
          ods_fa2mem_chem(1:kproma,:,10) = &
               odscl_fa2mem(1:kproma,:,11)      !HCFC141B
          ods_fa2mem_chem(1:kproma,:,11) = &
               odscl_fa2mem(1:kproma,:,12)      !HCFC142B
          ods_fa2mem_chem(1:kproma,:,20) = &
               odscl_fa2mem(1:kproma,:,17)      !H1211 (Cl)
          ods_fa2mem_chem(1:kproma,:,14) = &
               odscl_fa2mem(1:kproma,:,20)      !CH3CL
          ods_fa2mem_chem(1:kproma,:,15) = &
               odscl_fa2mem(1:kproma,:,21)      !HCFC21
          ods_fa2mem_chem(1:kproma,:,16) = &
               odscl_fa2mem(1:kproma,:,7)       !HCFC123
       ELSE
          ods_fa2mem_chem(1:kproma,:,1)  = &
               odscl_fa2memconst(1:kproma,:,1)  !CFC11
          ods_fa2mem_chem(1:kproma,:,2)  = &
               odscl_fa2memconst(1:kproma,:,2)  !CFC12
          ods_fa2mem_chem(1:kproma,:,4)  = &
               odscl_fa2memconst(1:kproma,:,3)  !CFC113
          ods_fa2mem_chem(1:kproma,:,5)  = &
               odscl_fa2memconst(1:kproma,:,4)  !CFC114
          ods_fa2mem_chem(1:kproma,:,6)  = &
               odscl_fa2memconst(1:kproma,:,5)  !CFC115
          ods_fa2mem_chem(1:kproma,:,7)  = &
               odscl_fa2memconst(1:kproma,:,15) !CCL4
          ods_fa2mem_chem(1:kproma,:,8)  = &
               odscl_fa2memconst(1:kproma,:,16) !CH3CCL3
          ods_fa2mem_chem(1:kproma,:,9)  = &
               odscl_fa2memconst(1:kproma,:,6)  !HCFC22
          ods_fa2mem_chem(1:kproma,:,10) = &
               odscl_fa2memconst(1:kproma,:,11) !HCFC141B
          ods_fa2mem_chem(1:kproma,:,11) = &
               odscl_fa2memconst(1:kproma,:,12) !HCFC142B
          ods_fa2mem_chem(1:kproma,:,20) = &
               odscl_fa2memconst(1:kproma,:,17) !H1211 (Cl)
          ods_fa2mem_chem(1:kproma,:,14) = &
               odscl_fa2memconst(1:kproma,:,20) !CH3CL
          ods_fa2mem_chem(1:kproma,:,15) = &
               odscl_fa2memconst(1:kproma,:,21) !HCFC21
          ods_fa2mem_chem(1:kproma,:,16) = &
               odscl_fa2memconst(1:kproma,:,7)  !HCFC123
       ENDIF

       IF (lodsbr_var_chem) THEN
          ods_fa2mem_chem(1:kproma,:,3)  = &
               odsbr_fa2mem(1:kproma,:,18)      !H1301
          ods_fa2mem_chem(1:kproma,:,12) = &
               odsbr_fa2mem(1:kproma,:,17)      !H1211 (Br)
          ods_fa2mem_chem(1:kproma,:,13) = &
               odsbr_fa2mem(1:kproma,:,19)      !CH3Br
          ods_fa2mem_chem(1:kproma,:,17) = &
               odsbr_fa2mem(1:kproma,:,22)      !H2402
          ods_fa2mem_chem(1:kproma,:,18) = &
               odsbr_fa2mem(1:kproma,:,23)      !CHBr3
          ods_fa2mem_chem(1:kproma,:,19) = &
               odsbr_fa2mem(1:kproma,:,24)      !CH2Br2
       ELSE
          ods_fa2mem_chem(1:kproma,:,3)  = &
               odsbr_fa2memconst(1:kproma,:,18) !H1301
          ods_fa2mem_chem(1:kproma,:,12) = &
               odsbr_fa2memconst(1:kproma,:,17) !H1211 (Br)
          ods_fa2mem_chem(1:kproma,:,13) = &
               odsbr_fa2memconst(1:kproma,:,19) !CH3Br
          ods_fa2mem_chem(1:kproma,:,17) = &
               odsbr_fa2memconst(1:kproma,:,22) !H2402
          ods_fa2mem_chem(1:kproma,:,18) = &
               odsbr_fa2memconst(1:kproma,:,23) !CHBr3
          ods_fa2mem_chem(1:kproma,:,19) = &
               odsbr_fa2memconst(1:kproma,:,24) !CH2Br2
       ENDIF
    ENDIF

  END SUBROUTINE interpolate_socol_ods_fa2mem


  SUBROUTINE read_ghg_ods_monthfile(filename, ghg_name, radchem_name, lvar, &
       yr, cy, ghg_y, ods_y, n_cl, n_br, ltods, odsmem_yb)

    ! Reads GHG/ODS monthly file for months of current year plus December of
    ! preceeding as well as January of following year.

    ! *read_ghg_ods_monthfile* is called from *read_socol_co2_global*,
    ! *read_socol_ch4_global*, *read_socol_n2o_global*, 
    ! *read_socol_odscl_global*, *read_socol_odsbr_global*, and 
    ! *read_socol_odsmem_global*.

    ! Subroutine arguments:
    CHARACTER (*), INTENT(in) :: filename
    CHARACTER (*), INTENT(in) :: ghg_name
    CHARACTER (*), INTENT(in) :: radchem_name
    LOGICAL, INTENT(in) :: lvar
    INTEGER, INTENT(in) :: yr, cy
    INTEGER, OPTIONAL, INTENT(out)  :: n_cl(nods_tab), n_br(nods_tab)
    REAL(dp), OPTIONAL, INTENT(out) :: ghg_y(0:13), ods_y(nods_tab,0:13), &
         ltods(nods_tab), odsmem_yb(nods_tab,-nodsyb*12+1:13)

    ! Local variables:
    CHARACTER (7) :: notused
    LOGICAL :: lods, lodsyb
    INTEGER  :: i, l, m, ii, i0, i1, ii0, m0, m1, yearmo, yearmo0, yearmo1, &
         y0, y1, y2, ioerror
    INTEGER, PARAMETER :: lun = 20
    INTEGER :: veci(nods_tab)
    REAL(dp) :: vecr(nods_tab)
    CHARACTER(4) :: yr0_st, yr1_st, yr2_st
    CHARACTER(2) :: mo0_st, mo1_st

    ioerror = 0
    lods = .FALSE.
    lodsyb = .FALSE.

    IF (PRESENT(ods_y) .OR. PRESENT(odsmem_yb)) lods = .TRUE.
    IF (PRESENT(odsmem_yb)) lodsyb = .TRUE. 

    ! Determine year:
    IF (lvar) THEN
       y1 = yr*100
    ELSE
       y1 = cy*100
    ENDIF

    y0 = y1-100
    y2 = y1+100       

    IF (lodsyb .AND. lvar) THEN
       yearmo0 = (yr-nodsyb)*100+1
       i0  = -nodsyb
       ii0 = -nodsyb*12+1
    ELSE
       yearmo0 = y0+12
       i0  = 1
       ii0 = 0
    ENDIF

    ! Open file:
    IF (p_parallel_io) THEN
       OPEN (lun, ERR=110, FILE=TRIM(filename), STATUS='old', &
            FORM='formatted', ACTION='read')
       GOTO 111
110    ioerror = 1
111    CONTINUE
    ENDIF

    IF (p_parallel) CALL p_bcast (ioerror, p_io)   
    IF (ioerror .EQ. 1) THEN
       WRITE(message_text,*) 'Could not open ', TRIM(filename)
       CALL message('',TRIM(message_text))
       CALL finish('read_ghg_ods_monthfile','Run terminated')
    ENDIF

    ! Read file:
    IF (p_parallel_io) THEN
       IF (lods) THEN  !ODS-File
          DO i = 1, 5
             READ (lun,'(A1)') notused    ! comment
          ENDDO
          READ (lun, *) veci           ! # Cl
          IF (PRESENT(n_cl)) n_cl(:) = veci(:)
          READ (lun,'(A1)') notused    ! comment
          READ (lun,'(A1)') notused    ! comment
          READ (lun, *) veci           ! # Br
          IF (PRESENT(n_br)) n_br(:) = veci(:)
          READ (lun,'(A1)') notused    ! comment
          READ (lun,'(A1)') notused    ! comment
          READ (lun, *) vecr           ! lifetime
          IF (PRESENT(ltods)) ltods(:) = vecr(:)
          READ (lun,'(A1)') notused    ! comment          
       ENDIF
       DO
          IF (lods) THEN
             IF (lodsyb) THEN
                READ (lun,*, END=120) yearmo, odsmem_yb(:,ii0)
             ELSE
                READ (lun,*, END=120) yearmo, ods_y(:,0)
             ENDIF
          ELSE
             READ (lun,*, END=120) yearmo, ghg_y(0)
          ENDIF
          IF (yearmo .EQ. yearmo0) EXIT
       ENDDO
120    CONTINUE
    ENDIF
    
    IF (p_parallel) CALL p_bcast (yearmo, p_io) 
    IF (lvar .AND. yearmo .NE. yearmo0) THEN
       WRITE(yr0_st, '(i4)') yearmo0/100
       IF (yearmo0-yearmo0/100 .LE. 9) THEN
          WRITE(mo0_st, '(i1)') yearmo0-yearmo0/100
       ELSE
          WRITE(mo0_st, '(i2)') yearmo0-yearmo0/100
       ENDIF
       WRITE(message_text,*) &
            'Year ', yr0_st, ', month ', TRIM(mo0_st), ' not found in file <', &
            TRIM(filename),'>'
       CALL message('',TRIM(message_text))
       CALL finish ('read_ghg_ods_monthfile', 'Run terminated.')
    ENDIF

    IF (p_parallel_io) THEN
       IF (lodsyb) THEN

          !Read file: 
          DO i = i0, 1
             IF (i .EQ. -nodsyb) THEN 
                m0 = 2
             ELSE 
                m0 = 1
             ENDIF
             IF (i .EQ. 1 .AND. lvar) THEN 
                m1 = 1
             ELSE 
                m1 = 12
             ENDIF
             DO m = m0, m1
                IF (lvar) THEN
                   yearmo1 = (yr+i)*100+m
                   ii = i*12+m
                ELSE
                   yearmo1 = y1+m
                   ii = m
                ENDIF
            
                READ (lun,*, END=121) yearmo, odsmem_yb(:,ii)

121             IF (yearmo .NE. yearmo1) THEN
                   ioerror = 1
                   EXIT
                ENDIF
             ENDDO
          ENDDO

       ELSE
          IF (lvar) THEN
             i1 = 13
          ELSE
             i1 = 12
          ENDIF

          !Read file: 
          DO  i=1, i1
             IF (i .EQ. 13) THEN
                yearmo1 = y2+1
             ELSE
                yearmo1 = y1+i
             ENDIF
             IF (lods) THEN
                READ (lun,*, END=123) yearmo, ods_y(:,i)
             ELSE
                READ (lun,*, END=123) yearmo, ghg_y(i)
             ENDIF

123          IF (yearmo .NE. yearmo1) THEN
                ioerror = 1
                EXIT
             ENDIF
          ENDDO

       ENDIF
    ENDIF

    IF (p_parallel) CALL p_bcast (ioerror, p_io) 
    IF (ioerror .GE. 1) THEN
       IF (p_parallel) THEN
          CALL p_bcast (yearmo, p_io)
          CALL p_bcast (i, p_io)
       ENDIF
       IF (ioerror .EQ. 1) THEN
          WRITE(yr1_st, '(i4)') yearmo1/100
          IF (yearmo1-yearmo1/100 .LE. 9) THEN
             WRITE(mo1_st, '(i1)') yearmo1-yearmo1/100
          ELSE
             WRITE(mo1_st, '(i2)') yearmo1-yearmo1/100
          ENDIF
          WRITE(message_text,*) &
               'Year ', yr1_st, ', month ', TRIM(mo1_st), &
               ' not found in file <', TRIM(filename),'>'
          CALL message('',TRIM(message_text))
          CALL finish ('read_ghg_ods_monthfile', 'Run terminated.')
       ENDIF
    ENDIF

    IF (p_parallel_io .AND. .NOT. lvar) THEN
       ! Mean of monthly values:
       IF (lodsyb) THEN
          DO l=1, nods_tab
             odsmem_yb(l,:) = SUM(odsmem_yb(l,1:12))/12._dp
          ENDDO
       ELSE
          IF (lods) THEN
             DO l=1, nods_tab
                ods_y(l,:) = SUM(ods_y(l,1:12))/12._dp
             ENDDO
          ELSE
             ghg_y(:) = SUM(ghg_y(1:12))/12._dp
          ENDIF
       ENDIF
    ENDIF

    IF (p_parallel_io) THEN
       ! Print message:
          WRITE (message_text,*) 'Reading ', TRIM(ghg_name), ' ', &
               TRIM(radchem_name), ' from file ', TRIM(filename)
          CALL message('',TRIM(message_text))
          IF (lvar) THEN
             WRITE(yr0_st, '(i4)') y0/100
             WRITE(yr1_st, '(i4)') y1/100
             WRITE(yr2_st, '(i4)') y2/100
             WRITE (message_text,*) ' for years ', yr0_st, ', ', yr1_st, &
                  ', ', yr2_st
             CALL message('',TRIM(message_text))
             CALL message('', &
                  '  (annually and monthly changing, globally constant values)')
          ELSE
             WRITE(yr1_st, '(i4)') y1/100
             WRITE (message_text,*) ' for year ', yr1_st
             CALL message('',TRIM(message_text))
             CALL message('',&
                  '  (annually and monthly constant, globally constant value)')
          ENDIF

       ! Close file:
       CLOSE (lun)
    ENDIF

    ! Distribute arrays:
    IF (p_parallel) THEN
       IF (lods) THEN 
          IF (lodsyb) THEN
             CALL p_bcast (odsmem_yb, p_io)
          ELSE
             CALL p_bcast (ods_y, p_io)
          ENDIF
          ! Only arrays of type REAL allowed for *b_cast*:
          IF (PRESENT(n_cl)) THEN
             vecr = REAL(n_cl,dp)
             CALL p_bcast (vecr, p_io)
             n_cl = NINT(vecr)
          ENDIF
          IF (PRESENT(n_br)) THEN
             vecr = REAL(n_br,dp)
             CALL p_bcast (vecr, p_io)
             n_br = NINT(vecr)
          ENDIF
          IF (PRESENT(ltods)) CALL p_bcast (ltods, p_io)
       ELSE
          CALL p_bcast (ghg_y, p_io)
       ENDIF
    ENDIF
        
  END SUBROUTINE read_ghg_ods_monthfile


  SUBROUTINE read_ghg_3dfile(ghg_varname, ghg_name, lvar, yr, mo, cy, &
       ghg_3dflev_m3_loc)

    ! Reads 3d GHG netcdf-files for months of current year plus December of
    ! preceeding as well as January of following year. Output on pressure
    ! grid as in the netcdf-file.

    ! *read_ghg_3dfile* is called from *read_socol_ch4_3d*,
    ! *read_socol_n2o_3d*, and *read_socol_odscl_3d*.
    
    ! Subroutine arguments:
    CHARACTER (*), INTENT(in) :: ghg_varname, ghg_name
    LOGICAL, INTENT(in) :: lvar
    INTEGER, INTENT(in) :: yr, mo, cy
    REAL(dp), INTENT(out), ALLOCATABLE :: ghg_3dflev_m3_loc(:,:,:,:)

    ! Local variables:
    CHARACTER (8)  :: fn0, fn1, fn2
    LOGICAL  :: lex0, lex1, lex2
    INTEGER  :: start(4), count(4)
    INTEGER  :: i, jk, y0, y1, y2, ioerror, io_nlon, io_ngl
    TYPE (FILE_INFO) :: gg3dnc0, gg3dnc1, gg3dnc2 
    REAL(dp), ALLOCATABLE, TARGET :: zin (:,:,:,:)
    REAL(dp), POINTER :: ghg_3dflev_m3_gl(:,:,:,:)
    CHARACTER(4)  :: yr0_st, yr1_st, yr2_st
    CHARACTER(2)  :: mo0_st, mo1_st, mo2_st


    ! Executable statements:

    ioerror = 0

    ! Initial value for open/close status of gg3dnc0, gg3dnc1, and gg3dnc2:
    gg3dnc0%opened = .FALSE.
    gg3dnc1%opened = .FALSE.
    gg3dnc2%opened = .FALSE.

    ! Determine year:
    IF (lvar) THEN
       y1 = yr
       y0 = y1-1
       y2 = y1+1
    ELSE
       y1 = cy
    ENDIF

    ! Open file(s):
    IF (p_parallel_io) THEN
       WRITE (fn1, '("gg3d", i4)') y1
       INQUIRE (FILE=fn1, EXIST=lex1)
       IF (lex1) CALL IO_open (fn1, gg3dnc1, IO_READ)
       IF (lvar) THEN
          WRITE (fn0, '("gg3d", i4)') y0
          WRITE (fn2, '("gg3d", i4)') y2
          IF (mo .EQ. 1) THEN
             INQUIRE (FILE=fn0, EXIST=lex0)
             IF (lex0) CALL IO_open (fn0, gg3dnc0, IO_READ)
          ENDIF
          IF (mo .EQ. 12) THEN
             INQUIRE (FILE=fn2, EXIST=lex2)
             IF (lex2) CALL IO_open (fn2, gg3dnc2, IO_READ)
          ENDIF
       ELSE
          gg3dnc0 = gg3dnc1 
          gg3dnc2 = gg3dnc1
       ENDIF
    ENDIF

    ! Terminate program if file(s) could not be opened:
    IF (p_parallel) CALL p_bcast (lex1, p_io)
    IF (.NOT. lex1) THEN
       WRITE (message_text,*) 'Could not open file <',fn1,'>'
       CALL message('',TRIM(message_text))
       CALL finish ('read_ghg_3dfile', 'Run terminated.')
    ENDIF
    IF (lvar) THEN
       IF (mo .EQ. 1) THEN
          IF (p_parallel) CALL p_bcast (lex0, p_io)
          IF (.NOT. lex0) THEN
             WRITE (message_text,*) 'Could not open file <',fn1,'>'
             CALL message('',TRIM(message_text))
             CALL finish ('read_ghg_3dfile', 'Run terminated.')
          ENDIF
       ENDIF
       IF (mo .EQ. 12) THEN
          IF (p_parallel) CALL p_bcast (lex2, p_io)
          IF (.NOT. lex2) THEN
             WRITE (message_text,*) 'Could not open file <',fn2,'>'
             CALL message('',TRIM(message_text))
             CALL finish ('read_ghg_3dfile', 'Run terminated.')
          ENDIF
       ENDIF 
    ENDIF

    ! Read file(s):
    IF (p_parallel_io) THEN

       ! 1. Current year:
       !Check resolution:
       CALL IO_inq_dimid  (gg3dnc1%file_id, 'lat', io_var_id)
       CALL IO_inq_dimlen (gg3dnc1%file_id, io_var_id, io_ngl)
       CALL IO_inq_dimid  (gg3dnc1%file_id, 'lon', io_var_id)
       CALL IO_inq_dimlen (gg3dnc1%file_id, io_var_id, io_nlon)
       CALL IO_inq_dimid  (gg3dnc1%file_id, 'p', io_var_id)
       CALL IO_inq_dimlen (gg3dnc1%file_id, io_var_id, nflev)
             
       IF  (io_ngl /= ngl .OR. io_nlon /= nlon) THEN
          WRITE (message_text,*) 'Unexpected resolution of 3d ', &
               TRIM(ghg_name),'-field: ', io_nlon, io_ngl
          CALL message('',TRIM(message_text))
          CALL finish ('read_ghg_3dfile', 'unexpected resolution')
       ENDIF

       ALLOCATE(zin(nlon,nflev,ngl,0:2)) 
       IF (.NOT. ALLOCATED(flev)) ALLOCATE(flev(nflev))
       
       CALL IO_inq_varid (gg3dnc1%file_id, 'p', io_var_id)
       CALL IO_get_var_double (gg3dnc1%file_id, io_var_id, flev)
       
       ! Read values:
       CALL IO_inq_varid (gg3dnc1%file_id, TRIM(ghg_varname), io_var_id)

       ! Preceeding month:
       IF (mo .NE. 1) THEN
          DO jk = 1, nflev
             count(:) = (/ nlon, ngl,  1, 1 /)
             start(:) = (/    1,   1, jk, mo-1 /)   
             CALL IO_get_vara_double (gg3dnc1%file_id, io_var_id, &
               start, count, zin(:,jk,:,0))
          ENDDO
       ENDIF

       ! Current month:
       DO jk = 1, nflev
          count(:) = (/ nlon, ngl,  1, 1 /)
          start(:) = (/    1,   1, jk, mo /)   
          CALL IO_get_vara_double (gg3dnc1%file_id, io_var_id, &
               start, count, zin(:,jk,:,1))
       ENDDO

       ! Following month:
       IF (mo .NE. 12) THEN
          DO jk = 1, nflev
             count(:) = (/ nlon, ngl,  1, 1 /)
             start(:) = (/    1,   1, jk, mo+1 /)   
             CALL IO_get_vara_double (gg3dnc1%file_id, io_var_id, &
                  start, count, zin(:,jk,:,2))
          ENDDO
       ENDIF 
       
       ! 2. December of preceeding year (current month is January):
       IF (mo .EQ. 1) THEN
          ! Read values:
          CALL IO_inq_varid (gg3dnc0%file_id, TRIM(ghg_varname), io_var_id)
          DO jk = 1, nflev
             count(:) = (/ nlon, ngl,  1,  1 /)
             start(:) = (/    1,   1, jk, 12 /)
             CALL IO_get_vara_double (gg3dnc0%file_id, io_var_id, &
                  start, count, zin(:,jk,:,0))
          ENDDO
       ENDIF
    
       ! 3. January of following year (current month is December):
       IF (mo .EQ. 12) THEN     
          ! Read values:
          CALL IO_inq_varid (gg3dnc2%file_id, TRIM(ghg_varname), io_var_id)
          DO jk = 1, nflev
             count(:) = (/ nlon, ngl,  1,  1 /)
             start(:) = (/    1,   1, jk,  1 /)
             CALL IO_get_vara_double (gg3dnc2%file_id, io_var_id, &
                  start, count, zin(:,jk,:,2))
          ENDDO
       ENDIF

       ! Close file(s):
       CALL IO_close (gg3dnc1)
       IF (lvar .AND. mo .EQ. 1) CALL IO_close (gg3dnc0)
       IF (lvar .AND. mo .EQ. 12) CALL IO_close (gg3dnc2)

       ! Print message:
       IF (mo-1 .LE. 9) THEN
          WRITE(mo0_st, '(i1)') mo-1
       ELSE
          WRITE(mo0_st, '(i2)') mo-1
       ENDIF
       IF (mo .LE. 9) THEN
          WRITE(mo1_st, '(i1)') mo
       ELSE
          WRITE(mo1_st, '(i2)') mo
       ENDIF
       IF (mo+1 .LE. 9) THEN
          WRITE(mo2_st, '(i1)') mo+1
       ELSE
          WRITE(mo2_st, '(i2)') mo+1
       ENDIF

       IF (lvar) THEN
          WRITE(yr0_st, '(i4)') y0
          WRITE(yr1_st, '(i4)') y1
          WRITE(yr2_st, '(i4)') y2
          IF (mo .NE. 1 .AND. mo .NE. 12) THEN 
             WRITE (message_text,*) 'Reading 3-dim ', TRIM(ghg_name), &
                  '-fields for radiation from file ', fn1
             CALL message('',TRIM(message_text))
             WRITE (message_text,*) ' for year ', yr1_st, ', months ', &
                  TRIM(mo0_st), ', ', TRIM(mo1_st), ', ', TRIM(mo2_st)
          ELSE
             IF (mo .EQ. 1) THEN
                WRITE (message_text,*) 'Reading 3-dim ', TRIM(ghg_name), &
                     '-fields for radiation from files ', fn0, ', ', fn1
                CALL message('',TRIM(message_text))
                WRITE (message_text,*) ' for year ', yr0_st, &
                     ', month 12 and for year ', yr1_st, ', months ', &
                     TRIM(mo1_st), ', ', TRIM(mo2_st)
             ELSE
                WRITE (message_text,*) 'Reading 3-dim ', TRIM(ghg_name), &
                     '-fields for radiation from files ', fn1, ', ', fn2
                CALL message('',TRIM(message_text))
                WRITE (message_text,*) ' for year ', yr1_st, &
                     ', months ', TRIM(mo0_st), ', ', TRIM(mo1_st), &
                     ', and for year ', yr2_st, ', month 1'
             ENDIF
          ENDIF
          CALL message('',TRIM(message_text))
          WRITE (message_text,*) '  (annually and monthly changing 3-d-values)'
          CALL message('',TRIM(message_text))
       ELSE
          WRITE (message_text,*) 'Reading 3-dim ', TRIM(ghg_name), &
               '-fields for radiation from file ', fn1
          CALL message('',TRIM(message_text))
          IF (mo .NE. 1 .AND. mo .NE. 12) THEN
             WRITE (message_text,*) ' for months ', TRIM(mo0_st), ', ', &
                  TRIM(mo1_st), ', ', TRIM(mo2_st)
          ELSE
             IF (mo .EQ. 1) THEN
                WRITE (message_text,*) ' for months ', TRIM(mo1_st), ', ', &
                     TRIM(mo2_st), ', 12'
             ELSE
                WRITE (message_text,*) ' for months 1, ', TRIM(mo0_st), ', ', &
                     TRIM(mo1_st)
             ENDIF
          ENDIF
          CALL message('',TRIM(message_text))
          WRITE (message_text,*) '  (annually constant 3-d-values)'
          CALL message('',TRIM(message_text))
       ENDIF
    ENDIF

    ! Distribute pressure grid:
    IF (p_parallel) THEN
       CALL p_bcast (nflev, p_io)
       CALL p_bcast (flev, p_io)
    ENDIF

    ! Allocate memory for local fields:
    IF (.NOT. ALLOCATED(ghg_3dflev_m3_loc)) &
         ALLOCATE(ghg_3dflev_m3_loc(lc%nproma,lc%nlev,lc%ngpblks,0:2))

    NULLIFY (ghg_3dflev_m3_gl)
    DO i = 0, 2
       IF (p_parallel_io) ghg_3dflev_m3_gl => zin(:,:,:,i:i)

       ! Scatter the field over the processors:
       CALL scatter_gp (ghg_3dflev_m3_gl, ghg_3dflev_m3_loc(:,:,:,i:i), &
            global_decomposition)
    ENDDO

    ! Deallocate memory:
    IF (p_parallel_io) DEALLOCATE (zin)

    ! hPa -> Pa
    flev(:) = 100._dp * flev(:)
        
  END SUBROUTINE read_ghg_3dfile


  SUBROUTINE calculate_fixlev_interpolation(kproma, kbdim, klev, p, fixlev, &
       nfixlev, k0, zk, zkp1)

    ! Prepares interpolation from fixlev-pressure grid to model pressure grid
    ! (calculated in *interpolate_3dflev*).

    ! *calculate_fixlev_interpolation* is called from *interpolate_3dflev*.

    ! Subroutine arguments:
    INTEGER, INTENT(in)   :: kproma, kbdim, klev, nfixlev
    REAL(dp), INTENT(in)  :: fixlev(nfixlev)     !pressure level
    REAL(dp), INTENT(in)  :: p(kbdim,klev)
    INTEGER, INTENT(out)  :: k0(kbdim,nfixlev)
    REAL(dp), INTENT(out) :: zk(kbdim,nfixlev), zkp1(kbdim,nfixlev)

    ! Local variables:
    LOGICAL :: k_flag   
    INTEGER :: k, jk, jl

    DO jk=1, klev
       k_flag = .TRUE.
       DO jl=1, kproma
          IF (p(jl,jk) .GE. fixlev(1) .AND. &
               p(jl,jk) .LE. fixlev(nfixlev)) THEN
             DO k=1, nfixlev-1
                IF (p(jl,jk) .GE. fixlev(k) .AND. &
                     p(jl,jk) .LT. fixlev(k+1) .AND. k_flag) THEN
                   k0(jl,jk) = k
                   zkp1(jl,jk) = fixlev(k+1) - p(jl,jk)
                   zk(jl,jk) = p(jl,jk) - fixlev(k)
                   k_flag = .FALSE.
                ENDIF
             ENDDO
          ELSE
             IF (p(jl,jk) .LT. fixlev(1)) THEN
                k0(jl,jk) = nfixlev
                zk(jl,jk) = 0._dp
                zkp1(jl,jk) = 1._dp
             ELSE
                k0(jl,jk) = 1
                zk(jl,jk) = 1._dp
                zkp1(jl,jk) = 0._dp
             ENDIF
          ENDIF
       ENDDO
    ENDDO
            
  END SUBROUTINE calculate_fixlev_interpolation


  SUBROUTINE interpolate_3dflev(krow, kproma, kbdim, klev, p, ghg_3dflev_m3, &
       ghg_3d)

    ! Interpolates to current time step and to model pressure grid.

    ! *interpolate_3dflev* is called from *interpolate_socol_ghg_ods_3d*.

    ! Subroutine arguments:
    INTEGER, INTENT(in)   :: krow, kproma, kbdim, klev
    REAL(dp), INTENT(in)  :: p(kbdim,klev)
    REAL(dp), INTENT(in)  :: ghg_3dflev_m3(kbdim,nflev,0:2)
    REAL(dp), INTENT(out) :: ghg_3d(kbdim,klev)

    ! Local variables:
    INTEGER  :: k, jk, jl
    INTEGER  :: k0(kbdim,nflev)
    REAL(dp) :: zk(kbdim,nflev), zkp1(kbdim,nflev)
    REAL(dp) :: ghg_3dflev(kbdim,nflev)

    ! Interpolation to current time step:
    ghg_3dflev(1:kproma,:) = wgt1_rad*ghg_3dflev_m3(1:kproma,:,m3w1_rad) + &
         wgt2_rad*ghg_3dflev_m3(1:kproma,:,m3w2_rad)
    
    ! Interpolation from flev-pressure grid to model pressure grid:
    CALL calculate_fixlev_interpolation(kproma, kbdim, klev, p, flev, nflev, &
         k0, zk, zkp1)

    DO jk=1, klev
       DO jl=1, kproma
          ghg_3d(jl,jk) = (zk(jl,jk)*ghg_3dflev(jl,k0(jl,jk)+1) &
               + zkp1(jl,jk)*ghg_3dflev(jl,k0(jl,jk))) &
               / (zk(jl,jk)+zkp1(jl,jk))
       ENDDO
    ENDDO

  END SUBROUTINE interpolate_3dflev


  FUNCTION ch4_chemrad(krow, kproma, klev)

    ! Conversion from chemistry module to radiation module for CH4.

    ! *ch4_chemrad* is called from *radiation*.

    USE mo_memory_g1a,         ONLY: xtm1
    
    INTEGER, INTENT(in) :: krow, kproma, klev
    REAL(dp) :: ch4_chemrad(kproma,klev) 

    ! Mixing ratio -> g/g:
    ch4_chemrad(:,:) =  amch4/amd*xtm1(1:kproma,:,idt_ch4,krow)
        
    
  END FUNCTION ch4_chemrad


  FUNCTION n2o_chemrad(krow, kproma, klev)

    ! Conversion from chemistry module to radiation module for N2O.

    ! *n2o_chemrad* is called from *radiation*.

    USE mo_memory_g1a,         ONLY: xtm1

    INTEGER, INTENT(in) :: krow, kproma, klev
    REAL(dp) :: n2o_chemrad(kproma,klev) 

    ! Mixing ratio -> g/g:
    n2o_chemrad(:,:) = amn2o/amd*xtm1(1:kproma,:,idt_n2o,krow)

  END FUNCTION n2o_chemrad


  FUNCTION cfc_chemrad(krow, kproma, klev)

    ! Conversion from chemistry module to radiation module for CFCs.

    ! *cfc_chemrad* is called from *radiation*.

    USE mo_memory_g1a,         ONLY: xtm1
    
    INTEGER, INTENT(in) :: krow, kproma, klev
    REAL(dp) :: cfc_chemrad(kproma,klev,ghg_no_cfc)

!!$    cfc_chemrad(:,:,1) = &
!!$                         (odscl_fa2mem_rad(1:kproma,:,1) + &   ! CFC-11
!!$                          odscl_fa2mem_rad(1:kproma,:,6) + &   ! HCFC-22
!!$                          odscl_fa2mem_rad(1:kproma,:,7) + &   ! HCFC-123
!!$                          odscl_fa2mem_rad(1:kproma,:,11)+ &   ! HCFC-141b
!!$                          odscl_fa2mem_rad(1:kproma,:,12)+ &   ! HCFC-142b
!!$                          odscl_fa2mem_rad(1:kproma,:,15)+ &   ! CCl4
!!$                          odscl_fa2mem_rad(1:kproma,:,16)+ &   ! CH3CCl3
!!$                          odscl_fa2mem_rad(1:kproma,:,17)+ &   ! H-1211
!!$                          odscl_fa2mem_rad(1:kproma,:,20)+ &   ! CH3Cl
!!$                          odscl_fa2mem_rad(1:kproma,:,21)) * & ! HCFC-21
!!$                          xtm1(1:kproma,:,idt_odscls,krow)
!!$    cfc_chemrad(:,:,2) = &
!!$                         (odscl_fa2mem_rad(1:kproma,:,2) + &   ! CFC-12
!!$                          odscl_fa2mem_rad(1:kproma,:,3) + &   ! CFC-113
!!$                          odscl_fa2mem_rad(1:kproma,:,4) + &   ! CFC-114
!!$                          odscl_fa2mem_rad(1:kproma,:,5)) * &  ! CFC-115
!!$                          xtm1(1:kproma,:,idt_odscll,krow)

    cfc_chemrad(:,:,1) = xtm1(1:kproma,:,idt_f11,krow) + &
                         xtm1(1:kproma,:,idt_hcfc22,krow) + &
                         xtm1(1:kproma,:,idt_hcfc123,krow) + &
                         xtm1(1:kproma,:,idt_hcfc141b,krow) + &
                         xtm1(1:kproma,:,idt_hcfc142b,krow) + &
                         xtm1(1:kproma,:,idt_ccl4,krow) + &
                         xtm1(1:kproma,:,idt_ch3ccl3,krow) + &
                         xtm1(1:kproma,:,idt_h1211,krow) + &
                         xtm1(1:kproma,:,idt_ch3cl,krow) + &
                         xtm1(1:kproma,:,idt_hcfc21,krow)

    cfc_chemrad(:,:,2) = xtm1(1:kproma,:,idt_f12,krow) + &
                         xtm1(1:kproma,:,idt_cfc113,krow) + &
                         xtm1(1:kproma,:,idt_cfc114,krow) + &
                         xtm1(1:kproma,:,idt_cfc115,krow)

  END FUNCTION cfc_chemrad


  FUNCTION o3_chemrad(krow, kproma, klev)

    ! Conversion from chemistry module to radiation module for O3.

    ! *o3_chemrad* is called from *radiation*.

    USE mo_memory_g1a,         ONLY: xtm1

    INTEGER, INTENT(in)  :: krow, kproma, klev
    REAL(dp) :: o3_chemrad(kproma,klev) 

    ! Mixing ratio -> g/g:
    o3_chemrad(:,:) = amo3/amd*xtm1(1:kproma,:,idt_o3,krow)

  END FUNCTION o3_chemrad


  SUBROUTINE ghg_message(co2, ch4, n2o, cfcs)

    ! Prints GHG concentrations at planetary boundary layer used in radiation
    ! module.

    ! *ghg_message* is called from *radiation*.

    REAL(dp), INTENT(in) :: co2, ch4, n2o, cfcs(ghg_no_cfc)

    WRITE (message_text,'(a,f9.4,a,f8.2,a,f8.3)') &
         'Lower boundary conditions: CO2 = ', co2*1.0e06_dp*amd/amco2, &
         ' CH4 = ', ch4*1.0e09_dp*amd/amch4, ' N2O = ', n2o*1.0e09_dp*amd/amn2o
    CALL message('', TRIM(message_text))
    WRITE (message_text,'(a,8f7.2,/,7x,8f7.2)') &
         ' CFC = ', cfcs(1:ghg_no_cfc)*1.0e12_dp
    CALL message('', TRIM(message_text))

  END SUBROUTINE ghg_message


  SUBROUTINE cleanup_socol_ghg_ods

    ! Deallocates module variables.
    
    ! *cleanup_socol_ghg_ods* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(flev)) DEALLOCATE(flev)
    IF (ALLOCATED(ch4_3dflev_m3_rad)) DEALLOCATE(ch4_3dflev_m3_rad)
    IF (ALLOCATED(n2o_3dflev_m3_rad)) DEALLOCATE(n2o_3dflev_m3_rad)
    IF (ALLOCATED(odscls_3dflev_m3_rad)) DEALLOCATE(odscls_3dflev_m3_rad)
    IF (ALLOCATED(odscll_3dflev_m3_rad)) DEALLOCATE(odscll_3dflev_m3_rad)
    IF (ALLOCATED(odscls_3d_bcond_rad)) DEALLOCATE(odscls_3d_bcond_rad)
    IF (ALLOCATED(odscll_3d_bcond_rad)) DEALLOCATE(odscll_3d_bcond_rad)
    IF (ALLOCATED(ch4_3d_bcond_rad)) DEALLOCATE(ch4_3d_bcond_rad)
    IF (ALLOCATED(n2o_3d_bcond_rad)) DEALLOCATE(n2o_3d_bcond_rad)
    IF (ALLOCATED(cfc_3d_bcond_rad)) DEALLOCATE(cfc_3d_bcond_rad)
    IF (ALLOCATED(odscl_fa2mem_rad)) DEALLOCATE(odscl_fa2mem_rad)
    IF (ALLOCATED(ods_fa2mem_chem)) DEALLOCATE(ods_fa2mem_chem)

  END SUBROUTINE cleanup_socol_ghg_ods

END MODULE mo_socol_ghg_ods

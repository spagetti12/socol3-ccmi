MODULE mo_socol_ch4_streams 

  ! *mo_socol_ch4_streams* contains subroutines to define streams of 
  ! different methane sinks and of wetland methane emissions
  !
  ! Andrea Stenke, ETH Zurich, December 2009

  USE mo_kind,          ONLY: dp      
  USE mo_linked_list,   ONLY: HYBRID, NETCDF
  USE mo_memory_base,   ONLY: new_stream, add_stream_element,  &
                              default_stream_setting, add_stream_reference, &
                              delete_stream, t_stream
  USE mo_time_control,  ONLY: delta_time, lstart

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_stream_ch4   ! construct the ch4 stream
  PUBLIC :: destruct_stream_ch4    ! destruct the ch4 stream
  PUBLIC :: init_stream_ch4        ! initialize the ch4 stream
  PUBLIC :: accumulate_stream_ch4  ! accumulate stream elements

  TYPE (t_stream), PUBLIC, POINTER :: ch4     !the ch4 stream

  REAL(dp), PUBLIC, POINTER :: ch4descl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4desod_1(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4desod_2(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4desoh(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4dest01(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4emswet(:,:)
  REAL(dp), PUBLIC, POINTER :: ch4_emflux(:,:)

  REAL(dp), PUBLIC, POINTER :: ch4descl_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4desod_1_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4desod_2_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4desoh_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4dest01_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4emswet_d(:,:)
  REAL(dp), PUBLIC, POINTER :: ch4_emflux_d(:,:)

  ! eth_as_oh+
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: &
       ohprod_rn3, ohdest_rn8, ohprod_rn11, &  
       ohdest_rn12, ohdest_rn13, ohdest_rn21, &  
       ohdest_rn24a, ohprod_rn27, ohdest_rn31, &  
       ohprod_rh1, ohdest_rh2, ohdest_rh3, &   
       ohdest_rh4, ohdest_rh13, ohdest_rh11, &  
       ohprod_rh5, ohprod_rh6, ohprod_rh8, &   
       ohdest_rh12, ohprod_rh12a, ohprod_rh9, &   
       ohprod_rh14a, ohdest_rh15, ohprod_rh16, &  
       ohprod_rh18, ohdest_rc5a, ohdest_rc5b, &  
       ohprod_rc9, ohdest_rc10, ohprod_rc13a, & 
       ohdest_rc6, ohprod_rc26, ohdest_rb14, &  
       ohprod_rb15, ohprod_rb16, ohprod_rb17, &  
       ohprod_rb22, ohprod_rch1, ohdest_rch3, &  
       ohdest_rch7, ohprod_rch7a, ohdest_rch14, & 
       ohprod_rch15, ohprod_rt06, ohdest_rt07, &  
       ohprod_rt08, ohprod_rt11, ohprod_rt12, &  
       ohprod_rt17, ohdest_rt18, ohdest_rt19, &  
       ohprod_rt20, ohdest_rt26, ohdest_ods10, & 
       ohdest_ods11, ohdest_ods13, ohdest_ods14, & 
       ohdest_ods15, ohdest_ods16, ohdest_ods18, & 
       ohdest_ods19, ohdest_rti01, ohprod_rti02, &
       ohprod_rti08, ohdest_rti09, ohdest_rti10, &
       ohprod_rti11, ohprod_rti17, ohdest_rti18, &
       ohdest_rti19, ohdest_rti20, ohdest_rti22, &
       ohdest_rti31, ohdest_rti33, ohdest_rti34, &
       ohprod_rti35, ohprod_rti39, ohprod_rti44

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: &
       ohprod_rn3_d, ohdest_rn8_d, ohprod_rn11_d, &  
       ohdest_rn12_d, ohdest_rn13_d, ohdest_rn21_d, &  
       ohdest_rn24a_d, ohprod_rn27_d, ohdest_rn31_d, &  
       ohprod_rh1_d, ohdest_rh2_d, ohdest_rh3_d, &   
       ohdest_rh4_d, ohdest_rh13_d, ohdest_rh11_d, &  
       ohprod_rh5_d, ohprod_rh6_d, ohprod_rh8_d, &   
       ohdest_rh12_d, ohprod_rh12a_d, ohprod_rh9_d, &   
       ohprod_rh14a_d, ohdest_rh15_d, ohprod_rh16_d, &  
       ohprod_rh18_d, ohdest_rc5a_d, ohdest_rc5b_d, &  
       ohprod_rc9_d, ohdest_rc10_d, ohprod_rc13a_d, & 
       ohdest_rc6_d, ohprod_rc26_d, ohdest_rb14_d, &  
       ohprod_rb15_d, ohprod_rb16_d, ohprod_rb17_d, &  
       ohprod_rb22_d, ohprod_rch1_d, ohdest_rch3_d, &  
       ohdest_rch7_d, ohprod_rch7a_d, ohdest_rch14_d, & 
       ohprod_rch15_d, ohprod_rt06_d, ohdest_rt07_d, &  
       ohprod_rt08_d, ohprod_rt11_d, ohprod_rt12_d, &  
       ohprod_rt17_d, ohdest_rt18_d, ohdest_rt19_d, &  
       ohprod_rt20_d, ohdest_rt26_d, ohdest_ods10_d, & 
       ohdest_ods11_d, ohdest_ods13_d, ohdest_ods14_d, & 
       ohdest_ods15_d, ohdest_ods16_d, ohdest_ods18_d, & 
       ohdest_ods19_d, ohdest_rti01_d, ohprod_rti02_d, &
       ohprod_rti08_d, ohdest_rti09_d, ohdest_rti10_d, &
       ohprod_rti11_d, ohprod_rti17_d, ohdest_rti18_d, &
       ohdest_rti19_d, ohdest_rti20_d, ohdest_rti22_d, &
       ohdest_rti31_d, ohdest_rti33_d, ohdest_rti34_d, &
       ohprod_rti35_d, ohprod_rti39_d, ohprod_rti44_d

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: &
       ho2dest_rn3, ho2prod_rn20, &
       ho2dest_rn22, ho2prod_rn22a, & 
       ho2dest_rn27, ho2prod_rh3, &   
       ho2dest_rh4, ho2dest_rh5, & 
       ho2dest_rh6, ho2dest_rh7, & 
       ho2dest_rh7a, ho2prod_rh12, & 
       ho2prod_rh12a, ho2prod_rh10, &
       ho2dest_rh14, ho2dest_rh14a, &
       ho2prod_rc5a, ho2dest_rc11, &
       ho2dest_rc13, ho2dest_rc13a, & 
       ho2prod_rc14, ho2dest_rb02, &
       ho2prod_rb05, ho2dest_rb13, &
       ho2prod_rch6, ho2prod_rch10, &
       ho2dest_rch13, ho2prod_rch21, & 
       ho2prod_rt07, ho2dest_rt08, &  
       ho2dest_rt09, ho2prod_rt10, &  
       ho2prod_rt19, ho2prod_rti02, &
       ho2prod_rti04, ho2dest_rti06, &
       ho2prod_rti07, ho2prod_rti11, &
       ho2prod_rti12, ho2dest_rti13, &
       ho2prod_rti14, ho2prod_rti19, &
       ho2dest_rti23, ho2dest_rti24, &
       ho2prod_rti27, ho2prod_rti35, &
       ho2prod_rti36, ho2prod_rti37, &
       ho2prod_rti39, ho2prod_rti40, &
       ho2prod_rti41, ho2prod_rti42

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: &
       ho2dest_rn3_d, ho2prod_rn20_d, &
       ho2dest_rn22_d, ho2prod_rn22a_d, & 
       ho2dest_rn27_d, ho2prod_rh3_d, &   
       ho2dest_rh4_d, ho2dest_rh5_d, & 
       ho2dest_rh6_d, ho2dest_rh7_d, & 
       ho2dest_rh7a_d, ho2prod_rh12_d, & 
       ho2prod_rh12a_d, ho2prod_rh10_d, &
       ho2dest_rh14_d, ho2dest_rh14a_d, &
       ho2prod_rc5a_d, ho2dest_rc11_d, &
       ho2dest_rc13_d, ho2dest_rc13a_d, & 
       ho2prod_rc14_d, ho2dest_rb02_d, &
       ho2prod_rb05_d, ho2dest_rb13_d, &
       ho2prod_rch6_d, ho2prod_rch10_d, &
       ho2dest_rch13_d, ho2prod_rch21_d, & 
       ho2prod_rt07_d, ho2dest_rt08_d, &  
       ho2dest_rt09_d, ho2prod_rt10_d, &  
       ho2prod_rt19_d, ho2prod_rti02_d, &
       ho2prod_rti04_d, ho2dest_rti06_d, &
       ho2prod_rti07_d, ho2prod_rti11_d, &
       ho2prod_rti12_d, ho2dest_rti13_d, &
       ho2prod_rti14_d, ho2prod_rti19_d, &
       ho2dest_rti23_d, ho2dest_rti24_d, &
       ho2prod_rti27_d, ho2prod_rti35_d, &
       ho2prod_rti36_d, ho2prod_rti37_d, &
       ho2prod_rti39_d, ho2prod_rti40_d, &
       ho2prod_rti41_d, ho2prod_rti42_d

  ! eth_as_oh-

!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
  SUBROUTINE construct_stream_ch4

    ! Allocates output streams

    ! *construct_stream_ch4* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (ch4 ,'ch4', filetype=NETCDF)      
  
  ! add entries for geopotential and log surface pressure (used by the
  ! afterburner):
  CALL add_stream_reference (ch4, 'geosp'   ,'g3b'   ,lpost=.FALSE.)
  CALL add_stream_reference (ch4, 'lsp'     ,'sp'    ,lpost=.FALSE.)

  CALL default_stream_setting (ch4, lrerun    = .TRUE. , &
                                       leveltype = HYBRID , &
                                       table     = 199,     &
                                       laccu     = .TRUE.)

  CALL add_stream_element (ch4, 'CH4DESCL', ch4descl, lpost=.FALSE., &
                           longname='CH4+Cl->HCl+CH3', &
                           units='mol/cm3/s', code=100)
  CALL add_stream_element (ch4, 'CH4DESOD1', ch4desod_1, lpost=.FALSE., &
                           longname='CH4+O1D->OH+CH3', &
                           units='mol/cm3/s', code=101)
  CALL add_stream_element (ch4, 'CH4DESOD2', ch4desod_2, lpost=.FALSE., &
                           longname='CH4+O1D->CH2O+H2', &
                           units='mol/cm3/s', code=102)
  CALL add_stream_element (ch4, 'CH4DESOH', ch4desoh,lpost=.FALSE., &
                           longname='CH4+OH->H2O+CH3', &
                           units='mol/cm3/s', code=103)
  CALL add_stream_element (ch4, 'CH4DEST01', ch4dest01,lpost=.FALSE., &
                           longname='CH4->CH3+H', &
                           units='mol/cm3/s', code=104)
  CALL add_stream_element (ch4, 'CH4EMISWET', ch4emswet,lpost=.FALSE., &
                           longname='CH4 emissions from wetlands', &
                           units='kg(CH4)/m2/s', code=105)
  CALL add_stream_element (ch4, 'CH4_EMFLUX', ch4_emflux,lpost=.FALSE., &
                           longname='CH4 surface emission flux', &
                           units='molec/m2/s', code=105)

  CALL add_stream_element (ch4, 'CH4DESCL_d', ch4descl_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4+Cl->HCl+CH3', &
                           units='mol/cm3/s', code=106)
  CALL add_stream_element (ch4, 'CH4DESOD1_d', ch4desod_1_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4+O1D->OH+CH3', &
                           units='mol/cm3/s', code=107)
  CALL add_stream_element (ch4, 'CH4DESOD2_d', ch4desod_2_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4+O1D->CH2O+H2', &
                           units='mol/cm3/s', code=108)
  CALL add_stream_element (ch4, 'CH4DESOH_d', ch4desoh_d,lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4+OH->H2O+CH3', &
                           units='mol/cm3/s', code=109)
  CALL add_stream_element (ch4, 'CH4DEST01_d', ch4dest01_d,lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4->CH3+H', &
                           units='mol/cm3/s', code=110)
  CALL add_stream_element (ch4, 'CH4EMISWET_d', ch4emswet_d,lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4 emissions from wetlands', &
                           units='kg(CH4)/m2/s', code=111)
  CALL add_stream_element (ch4, 'CH4_EMFLUX_d', ch4_emflux_d,lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4 surface emission flux', &
                           units='molec/m2/s', code=111)

  ! eth_as_oh+
  CALL add_stream_element (ch4, 'OHPROD_RN3', ohprod_rn3, lpost=.FALSE., &
                           longname='NO+HO2->OH+NO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RN3_d', ohprod_rn3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO+HO2->OH+NO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RN8', ohdest_rn8, lpost=.FALSE., &
                           longname='NO2+OH+MU->HNO3+MU', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RN8_d', ohdest_rn8_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO2+OH+MU->HNO3+MU', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RN11', ohprod_rn11, lpost=.FALSE., &
                           longname='HNO3->NO2+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RN11_d', ohprod_rn11_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HNO3->NO2+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RN12', ohdest_rn12, lpost=.FALSE., &
                           longname='HNO3+OH->NO3+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RN12_d', ohdest_rn12_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HNO3+OH->NO3+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RN13', ohdest_rn13, lpost=.FALSE., &
                           longname='HNO3+OH+MU->NO3+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RN13_d', ohdest_rn13_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HNO3+OH+MU->NO3+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RN21', ohdest_rn21, lpost=.FALSE., &
                           longname='HNO4+OH->H2O+O2+NO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RN21_d', ohdest_rn21_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HNO4+OH->H2O+O2+NO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RN24a', ohdest_rn24a, lpost=.FALSE., &
                           longname='CLNO3+OH->HOCL+NO3', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RN24a_d', ohdest_rn24a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CLNO3+OH->HOCL+NO3', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RN27', ohprod_rn27, lpost=.FALSE., &
                           longname='N+HO2->NO+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RN27_d', ohprod_rn27_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='N+HO2->NO+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RN31', ohdest_rn31, lpost=.FALSE., &
                           longname='N+OH->NO+H', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RN31_d', ohdest_rn31_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='N+OH->NO+H', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH1', ohprod_rh1, lpost=.FALSE., &
                           longname='H2O+OD->OH+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH1_d', ohprod_rh1_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2O+OD->OH+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RH2', ohdest_rh2, lpost=.FALSE., &
                           longname='OH+O->O2+H', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RH2_d', ohdest_rh2_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+O->O2+H', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RH3', ohdest_rh3, lpost=.FALSE., &
                           longname='OH+O3->HO2+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RH3_d', ohdest_rh3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+O3->HO2+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RH4', ohdest_rh4, lpost=.FALSE., &
                           longname='OH+HO2->H2O+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RH4_d', ohdest_rh4_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+HO2->H2O+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RH13', ohdest_rh13, lpost=.FALSE., &
                           longname='OH+OH->H2O+O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RH13_d', ohdest_rh13_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+OH->H2O+O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RH11', ohdest_rh11, lpost=.FALSE., &
                           longname='CO+OH->H+CO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RH11_d', ohdest_rh11_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CO+OH->H+CO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH5', ohprod_rh5, lpost=.FALSE., &
                           longname='HO2+O->OH+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH5_d', ohprod_rh5_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HO2+O->OH+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH6', ohprod_rh6, lpost=.FALSE., &
                           longname='HO2+O3->OH+O2+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH6_d', ohprod_rh6_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HO2+O3->OH+O2+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH8', ohprod_rh8, lpost=.FALSE., &
                           longname='H2O2->OH+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH8_d', ohprod_rh8_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2O2->OH+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RH12', ohdest_rh12, lpost=.FALSE., &
                           longname='H2O2+OH->HO2+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RH12_d', ohdest_rh12_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2O2+OH->HO2+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH12a', ohprod_rh12a, lpost=.FALSE., &
                           longname='H2O2+O->HO2+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH12a_d', ohprod_rh12a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2O2+O->HO2+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH9', ohprod_rh9, lpost=.FALSE., &
                           longname='H+O3->OH+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH9_d', ohprod_rh9_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H+O3->OH+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH14a', ohprod_rh14a, lpost=.FALSE., &
                           longname='H+HO2->OH+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH14a_d', ohprod_rh14a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H+HO2->OH+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RH15', ohdest_rh15, lpost=.FALSE., &
                           longname='OH+H2->H+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RH15_d', ohdest_rh15_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+H2->H+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH16', ohprod_rh16, lpost=.FALSE., &
                           longname='H2+OD->OH+H', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH16_d', ohprod_rh16_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2+OD->OH+H', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RH18', ohprod_rh18, lpost=.FALSE., &
                           longname='H2O->H+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RH18_d', ohprod_rh18_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2O->H+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RC5a', ohdest_rc5a, lpost=.FALSE., &
                           longname='CLO+OH->CL+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RC5a_d', ohdest_rc5a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CLO+OH->CL+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RC5b', ohdest_rc5b, lpost=.FALSE., &
                           longname='CLO+OH->HCL+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RC5b_d', ohdest_rc5b_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CLO+OH->HCL+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RC9', ohprod_rc9, lpost=.FALSE., &
                           longname='HOCL->CL+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RC9_d', ohprod_rc9_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HOCL->CL+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RC10', ohdest_rc10, lpost=.FALSE., &
                           longname='HOCL+OH->CLO+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RC10_d', ohdest_rc10_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HOCL+OH->CLO+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RC13a', ohprod_rc13a, lpost=.FALSE., &
                           longname='CL+HO2->CLO +OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RC13a_d', ohprod_rc13a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CL+HO2->CLO +OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RC6', ohdest_rc6, lpost=.FALSE., &
                           longname='HCL+OH->CL+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RC6_d', ohdest_rc6_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HCL+OH->CL+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RC26', ohprod_rc26, lpost=.FALSE., &
                           longname='HCL+O->CL+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RC26_d', ohprod_rc26_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HCL+O->CL+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RB14', ohdest_rb14, lpost=.FALSE., &
                           longname='HBR+OH->BR+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RB14_d', ohdest_rb14_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HBR+OH->BR+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RB15', ohprod_rb15, lpost=.FALSE., &
                           longname='HBR+O->BR+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RB15_d', ohprod_rb15_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HBR+O->BR+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RB16', ohprod_rb16, lpost=.FALSE., &
                           longname='HBR+OD->BR+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RB16_d', ohprod_rb16_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HBR+OD->BR+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RB17', ohprod_rb17, lpost=.FALSE., &
                           longname='HOBR+O->BRO+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RB17_d', ohprod_rb17_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HOBR+O->BRO+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RB22', ohprod_rb22, lpost=.FALSE., &
                           longname='HOBR->BR+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RB22_d', ohprod_rb22_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HOBR->BR+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RCH1', ohprod_rch1, lpost=.FALSE., &
                           longname='CH4+OD->OH+CH3', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RCH1_d', ohprod_rch1_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4+OD->OH+CH3', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RCH3', ohdest_rch3, lpost=.FALSE., &
                           longname='CH4+OH->H2O+CH3', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RCH3_d', ohdest_rch3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH4+OH->H2O+CH3', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RCH7', ohdest_rch7, lpost=.FALSE., &
                           longname='CH2O+OH->HCO+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RCH7_d', ohdest_rch7_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH2O+OH->HCO+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RCH7a', ohprod_rch7a, lpost=.FALSE., &
                           longname='CH2O+O->HCO+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RCH7a_d', ohprod_rch7a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH2O+O->HCO+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RCH14', ohdest_rch14, lpost=.FALSE., &
                           longname='CH3O2H+OH->CH3O2+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RCH14_d', ohdest_rch14_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3O2H+OH->CH3O2+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RCH15', ohprod_rch15, lpost=.FALSE., &
                           longname='CH3O2H->CH3O+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RCH15_d', ohprod_rch15_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3O2H->CH3O+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RT06', ohprod_rt06, lpost=.FALSE., &
                           longname='H2+O->OH+H', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RT06_d', ohprod_rt06_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2+O->OH+H', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RT07', ohdest_rt07, lpost=.FALSE., &
                           longname='NO3+OH->NO2+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RT07_d', ohdest_rt07_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO3+OH->NO2+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RT08', ohprod_rt08, lpost=.FALSE., &
                           longname='NO3+HO2->NO2+OH+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RT08_d', ohprod_rt08_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO3+HO2->NO2+OH+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RT11', ohprod_rt11, lpost=.FALSE., &
                           longname='O+HOCL->OH+CLO', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RT11_d', ohprod_rt11_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='O+HOCL->OH+CLO', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RT12', ohprod_rt12, lpost=.FALSE., &
                           longname='CL+HOCL->OH+CL2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RT12_d', ohprod_rt12_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CL+HOCL->OH+CL2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RT17', ohprod_rt17, lpost=.FALSE., &
                           longname='O+HCL->OH+CL', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RT17_d', ohprod_rt17_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='O+HCL->OH+CL', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RT18', ohdest_rt18, lpost=.FALSE., &
                           longname='OH+CL2->HOCL+CL', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RT18_d', ohdest_rt18_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+CL2->HOCL+CL', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RT19', ohdest_rt19, lpost=.FALSE., &
                           longname='BRO+OH->HO2+BR', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RT19_d', ohdest_rt19_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='BRO+OH->HO2+BR', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RT20', ohprod_rt20, lpost=.FALSE., &
                           longname='H+NO2->OH+NO', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RT20_d', ohprod_rt20_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H+NO2->OH+NO', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RT26', ohdest_rt26, lpost=.FALSE., &
                           longname='OH+OH+MU->H2O2+MU', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RT26_d', ohdest_rt26_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+OH+MU->H2O2+MU', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS10', ohdest_ods10, lpost=.FALSE., &
                           longname='HCFC141B->2*CL', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS10_d', ohdest_ods10_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HCFC141B->2*CL', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS11', ohdest_ods11, lpost=.FALSE., &
                           longname='HCFC142B->CL', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS11_d', ohdest_ods11_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HCFC142B->CL', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS13', ohdest_ods13, lpost=.FALSE., &
                           longname='CH3BR->BR', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS13_d', ohdest_ods13_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3BR->BR', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS14', ohdest_ods14, lpost=.FALSE., &
                           longname='CH3CL->CL', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS14_d', ohdest_ods14_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3CL->CL', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS15', ohdest_ods15, lpost=.FALSE., &
                           longname='HCFC21->2*CL', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS15_d', ohdest_ods15_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HCFC21->2*CL', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS16', ohdest_ods16, lpost=.FALSE., &
                           longname='HCFC123->2*CL', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS16_d', ohdest_ods16_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HCFC123->2*CL', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS18', ohdest_ods18, lpost=.FALSE., &
                           longname='CHBR3->3*BR', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS18_d', ohdest_ods18_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CHBR3->3*BR', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_ODS19', ohdest_ods19, lpost=.FALSE., &
                           longname='CH2BR2->2*BR', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_ODS19_d', ohdest_ods19_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH2BR2->2*BR', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI01', ohdest_rti01, lpost=.FALSE., &
                           longname='C5H8+OH->ISO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI01_d', ohdest_rti01_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='C5H8+OH->ISO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RTI02', ohprod_rti02, lpost=.FALSE., &
                           longname='C5H8+O3->OH+HO+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RTI02_d', ohprod_rti02_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='C5H8+O3->OH+HO+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RTI08', ohprod_rti08, lpost=.FALSE., &
                           longname='ISO2H+OH->MACR+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RTI08_d', ohprod_rti08_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='ISO2H+OH->MACR+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI09', ohdest_rti09, lpost=.FALSE., &
                           longname='ISON+OH->NALD', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI09_d', ohdest_rti09_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='ISON+OH->NALD', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI10', ohdest_rti10, lpost=.FALSE., &
                           longname='MACR+OH->MACRO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI10_d', ohdest_rti10_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACR+OH->MACRO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RTI11', ohprod_rti11, lpost=.FALSE., &
                           longname='MACR+O3->OH+HO2+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RTI11_d', ohprod_rti11_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACR+O3->OH+HO2+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RTI17', ohprod_rti17, lpost=.FALSE., &
                           longname='MPAN+OH->NO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RTI17_d', ohprod_rti17_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MPAN+OH->NO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI18', ohdest_rti18, lpost=.FALSE., &
                           longname='MACRO2H+OH = MACRO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI18_d', ohdest_rti18_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACRO2H+OH = MACRO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI19', ohdest_rti19, lpost=.FALSE., &
                           longname='HACET+OH->MGLY+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI19_d', ohdest_rti19_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HACET+OH->MGLY+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI20', ohdest_rti20, lpost=.FALSE., &
                           longname='MGLY+OH->CH3CO3+CO', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI20_d', ohdest_rti20_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MGLY+OH->CH3CO3+CO', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI22', ohdest_rti22, lpost=.FALSE., &
                           longname='NALD+OH->CH2O+CO+NO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI22_d', ohdest_rti22_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NALD+OH->CH2O+CO+NO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI31', ohdest_rti31, lpost=.FALSE., &
                           longname='PAN+OH->CH2O+NO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI31_d', ohdest_rti31_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='PAN+OH->CH2O+NO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI33', ohdest_rti33, lpost=.FALSE., &
                           longname='CH3CO3H+OH->CH3CO3', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI33_d', ohdest_rti33_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3CO3H+OH->CH3CO3', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHDEST_RTI34', ohdest_rti34, lpost=.FALSE., &
                           longname='CH3COOH+OH->CH3O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHDEST_RTI34_d', ohdest_rti34_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3COOH+OH->CH3O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RTI35', ohprod_rti35, lpost=.FALSE., &
                           longname='ISO2H->OH+MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RTI35_d', ohprod_rti35_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='ISO2H->OH+MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RTI39', ohprod_rti39, lpost=.FALSE., &
                           longname='MACRO2H->OH+HO2+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RTI39_d', ohprod_rti39_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACRO2H->OH+HO2+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'OHPROD_RTI44', ohprod_rti44, lpost=.FALSE., &
                           longname='CH3CO3H->CH3O2+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'OHPROD_RTI44_d', ohprod_rti44_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3CO3H->CH3O2+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RN3', ho2dest_rn3, lpost=.FALSE., &
                           longname='NO+HO2->OH+NO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RN3_d', ho2dest_rn3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO+HO2->OH+NO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RN22', ho2dest_rn22, lpost=.FALSE., &
                           longname='HO2+NO2+MU->HNO4+MU', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RN22_d', ho2dest_rn22_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HO2+NO2+MU->HNO4+MU', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RN27', ho2dest_rn27, lpost=.FALSE., &
                           longname='N+HO2->NO+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RN27_d', ho2dest_rn27_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='N+HO2->NO+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RH4', ho2dest_rh4, lpost=.FALSE., &
                           longname='OH+HO2->H2O+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RH4_d', ho2dest_rh4_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+HO2->H2O+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RH5', ho2dest_rh5, lpost=.FALSE., &
                           longname='HO2+O->OH+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RH5_d', ho2dest_rh5_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HO2+O->OH+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RH6', ho2dest_rh6, lpost=.FALSE., &
                           longname='HO2+O3->OH+O2+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RH6_d', ho2dest_rh6_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HO2+O3->OH+O2+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RH7', ho2dest_rh7, lpost=.FALSE., &
                           longname='HO2+HO2->H2O2+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RH7_d', ho2dest_rh7_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HO2+HO2->H2O2+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RH7a', ho2dest_rh7a, lpost=.FALSE., &
                           longname='HO2+HO2->H2O2+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RH7a_d', ho2dest_rh7a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HO2+HO2->H2O2+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RH14', ho2dest_rh14, lpost=.FALSE., &
                           longname='H+HO2->H2+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RH14_d', ho2dest_rh14_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H+HO2->H2+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RH14a', ho2dest_rh14a, lpost=.FALSE., &
                           longname='H+HO2->OH+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RH14a_d', ho2dest_rh14a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H+HO2->OH+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RC11', ho2dest_rc11, lpost=.FALSE., &
                           longname='CLO+HO2->HOCL+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RC11_d', ho2dest_rc11_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CLO+HO2->HOCL+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RC13', ho2dest_rc13, lpost=.FALSE., &
                           longname='CL+HO2->HCL+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RC13_d', ho2dest_rc13_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CL+HO2->HCL+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RC13a', ho2dest_rc13a, lpost=.FALSE., &
                           longname='CL+HO2->CLO+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RC13a_d', ho2dest_rc13a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CL+HO2->CLO+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RB02', ho2dest_rb02, lpost=.FALSE., &
                           longname='BR+HO2->HBR+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RB02_d', ho2dest_rb02_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='BR+HO2->HBR+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RB13', ho2dest_rb13, lpost=.FALSE., &
                           longname='BRO+HO2->HOBR+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RB13_d', ho2dest_rb13_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='BRO+HO2->HOBR+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RCH13', ho2dest_rch13, lpost=.FALSE., &
                           longname='CH3O2+HO2->CH3O2H+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RCH13_d', ho2dest_rch13_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3O2+HO2->CH3O2H+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RT08', ho2dest_rt08, lpost=.FALSE., &
                           longname='NO3+HO2->NO2+OH+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RT08_d', ho2dest_rt08_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO3+HO2->NO2+OH+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RT09', ho2dest_rt09, lpost=.FALSE., &
                           longname='NO3+HO2->HNO3+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RT09_d', ho2dest_rt09_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO3+HO2->HNO3+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RN20', ho2prod_rn20, lpost=.FALSE., &
                           longname='HNO4->NO2+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RN20_d', ho2prod_rn20_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HNO4->NO2+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RN22a', ho2prod_rn22a, lpost=.FALSE., &
                           longname='HNO4+MU+MU->HO2+NO2+MU', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RN22a_d', ho2prod_rn22a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HNO4+MU+MU->HO2+NO2+MU', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RH3', ho2prod_rh3, lpost=.FALSE., &
                           longname='OH+O3->HO2+O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RH3_d', ho2prod_rh3_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='OH+O3->HO2+O2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RH12', ho2prod_rh12, lpost=.FALSE., &
                           longname='H2O2+OH->HO2+H2O', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RH12_d', ho2prod_rh12_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2O2+OH->HO2+H2O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RH12a', ho2prod_rh12a, lpost=.FALSE., &
                           longname='H2O2+O->HO2+OH', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RH12a_d', ho2prod_rh12a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H2O2+O->HO2+OH', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RH10', ho2prod_rh10, lpost=.FALSE., &
                           longname='H+O2+MU->HO2+MU', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RH10_d', ho2prod_rh10_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='H+O2+MU->HO2+MU', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RC5a', ho2prod_rc5a, lpost=.FALSE., &
                           longname='CLO+OH->CL+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RC5a_d', ho2prod_rc5a_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CLO+OH->CL+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RC14', ho2prod_rc14, lpost=.FALSE., &
                           longname='CL+H2O2->HCL+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RC14_d', ho2prod_rc14_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CL+H2O2->HCL+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RB05', ho2prod_rb05, lpost=.FALSE., &
                           longname='BR+H2O2->HBR+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RB05_d', ho2prod_rb05_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='BR+H2O2->HBR+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RCH6', ho2prod_rch6, lpost=.FALSE., &
                           longname='CH3O+O2->CH2O+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RCH6_d', ho2prod_rch6_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3O+O2->CH2O+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RCH10', ho2prod_rch10, lpost=.FALSE., &
                           longname='HCO+O2->CO+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RCH10_d', ho2prod_rch10_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HCO+O2->CO+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RCH21', ho2prod_rch21, lpost=.FALSE., &
                           longname='CH3O2+CH3O2->CH2O+CH3O+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RCH21_d', ho2prod_rch21_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3O2+CH3O2->CH2O+CH3O+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RT07', ho2prod_rt07, lpost=.FALSE., &
                           longname='NO3+OH->NO2+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RT07_d', ho2prod_rt07_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NO3+OH->NO2+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RT10', ho2prod_rt10, lpost=.FALSE., &
                           longname='CH2O+NO3->CO+HO2+HNO3', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RT10_d', ho2prod_rt10_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH2O+NO3->CO+HO2+HNO3', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RT19', ho2prod_rt19, lpost=.FALSE., &
                           longname='BRO+OH->HO2+BR', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RT19_d', ho2prod_rt19_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='BRO+OH->HO2+BR', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI02', ho2prod_rti02, lpost=.FALSE., &
                           longname='C5H8+O3->OH+HO2+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI02_d', ho2prod_rti02_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='C5H8+O3->OH+HO2+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI04', ho2prod_rti04, lpost=.FALSE., &
                           longname='ISO2+NO->NO2+MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI04_d', ho2prod_rti04_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='ISO2+NO->NO2+MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RTI06', ho2dest_rti06, lpost=.FALSE., &
                           longname='ISO2+HO2->ISO2H', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RTI06_d', ho2dest_rti06_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='ISO2+HO2->ISO2H', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI07', ho2prod_rti07, lpost=.FALSE., &
                           longname='2ISO2->2MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI07_d', ho2prod_rti07_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='2ISO2->2MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI11', ho2prod_rti11, lpost=.FALSE., &
                           longname='MACR+O3->HO2+OH+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI11_d', ho2prod_rti11_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACR+O3->HO2+OH+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI12', ho2prod_rti12, lpost=.FALSE., &
                           longname='MACRO2+NO->HO2+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI12_d', ho2prod_rti12_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACRO2+NO->HO2+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RTI13', ho2dest_rti13, lpost=.FALSE., &
                           longname='MACRO2+HO2->MACRO2H', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RTI13_d', ho2dest_rti13_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACRO2+HO2->MACRO2H', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI14', ho2prod_rti14, lpost=.FALSE., &
                           longname='2MACRO2->HO2+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI14_d', ho2prod_rti14_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='2MACRO2->HO2+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI19', ho2prod_rti19, lpost=.FALSE., &
                           longname='HACET+OH->MGLY+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI19_d', ho2prod_rti19_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HACET+OH->MGLY+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RTI23', ho2dest_rti23, lpost=.FALSE., &
                           longname='CH3CO3+HO2->CH3CO3H', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RTI23_d', ho2dest_rti23_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3CO3+HO2->CH3CO3H', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2DEST_RTI24', ho2dest_rti24, lpost=.FALSE., &
                           longname='CH3CO3+HO2->CH3COOH+O3', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2ODEST_RTI24_d', ho2dest_rti24_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3CO3+HO2->CH3COOH+O3', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI27', ho2prod_rti27, lpost=.FALSE., &
                           longname='CH3CO3+CH3O2->CH2O+HO2+CH3O2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI27_d', ho2prod_rti27_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='CH3CO3+CH3O2->CH2O+HO2+CH3O', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI35', ho2prod_rti35, lpost=.FALSE., &
                           longname='ISO2H->OH+MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI35_d', ho2prod_rti35_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='ISO2H->OH+MACR+CH2O+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI36', ho2prod_rti36, lpost=.FALSE., &
                           longname='ISON->NO2+MACR+CH2O+HO', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI36_d', ho2prod_rti36_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='ISON->NO2+MACR+CH2O+HO', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI37', ho2prod_rti37, lpost=.FALSE., &
                           longname='MACR->CH3CO3+CH2O+CO+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI37_d', ho2prod_rti37_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACR->CH3CO3+CH2O+CO+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI39', ho2prod_rti39, lpost=.FALSE., &
                           longname='MACRO2H->OH+HO2+', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI39_d', ho2prod_rti39_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MACRO2H->OH+HO2+', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI40', ho2prod_rti40, lpost=.FALSE., &
                           longname='HACET->CH3CO3+CH2O+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI40_d', ho2prod_rti40_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='HACET->CH3CO3+CH2O+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI41', ho2prod_rti41, lpost=.FALSE., &
                           longname='MGLY->CH3CO3+CO+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI41_d', ho2prod_rti41_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='MGLY->CH3CO3+CO+HO2', &
                           units='mol/cm3/s', code=106)

  CALL add_stream_element (ch4, 'HO2PROD_RTI42', ho2prod_rti42, lpost=.FALSE., &
                           longname='NALD->CH2O+CO+NO2+HO2', &
                           units='mol/cm3/s', code=100)

  CALL add_stream_element (ch4, 'H2OPROD_RTI42_d', ho2prod_rti42_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='NALD->CH2O+CO+NO2+HO2', &
                           units='mol/cm3/s', code=106)
  ! eth_as_oh-

END SUBROUTINE construct_stream_ch4

!---------------------------------------------------------------------------------

SUBROUTINE destruct_stream_ch4

  ! Deallocates memory.

  ! *destruct_stream_ch4* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(ch4)

END SUBROUTINE destruct_stream_ch4

!---------------------------------------------------------------------------------

SUBROUTINE init_stream_ch4

  ! Initializes streams with zero.

  ! *init_stream_ch4* is called from *call_init_tracers*,
  ! src/call_submodels.f90.

  IF (lstart) THEN     ! use restart fields otherwise
     ch4descl(:,:,:)   = 0._dp
     ch4desod_1(:,:,:) = 0._dp
     ch4desod_2(:,:,:) = 0._dp
     ch4desoh(:,:,:)   = 0._dp
     ch4dest01(:,:,:)  = 0._dp
     ch4emswet(:,:)    = 0._dp
     ch4_emflux(:,:)   = 0._dp

     ch4descl_d(:,:,:)   = 0._dp
     ch4desod_1_d(:,:,:) = 0._dp
     ch4desod_2_d(:,:,:) = 0._dp
     ch4desoh_d(:,:,:)   = 0._dp
     ch4dest01_d(:,:,:)  = 0._dp
     ch4emswet_d(:,:)    = 0._dp
     ch4_emflux_d(:,:)   = 0._dp

     ! eth_as_oh+
     ohprod_rn3(:,:,:)   = 0._dp 
     ohdest_rn8(:,:,:)   = 0._dp 
     ohprod_rn11(:,:,:)  = 0._dp  
     ohdest_rn12(:,:,:)  = 0._dp
     ohdest_rn13(:,:,:)  = 0._dp
     ohdest_rn21(:,:,:)  = 0._dp  
     ohdest_rn24a(:,:,:) = 0._dp
     ohprod_rn27(:,:,:)  = 0._dp
     ohdest_rn31(:,:,:)  = 0._dp
     ohprod_rh1(:,:,:)   = 0._dp
     ohdest_rh2(:,:,:)   = 0._dp
     ohdest_rh3(:,:,:)   = 0._dp   
     ohdest_rh4(:,:,:)   = 0._dp
     ohdest_rh13(:,:,:)  = 0._dp 
     ohdest_rh11(:,:,:)  = 0._dp   
     ohprod_rh5(:,:,:)   = 0._dp 
     ohprod_rh6(:,:,:)   = 0._dp 
     ohprod_rh8(:,:,:)   = 0._dp  
     ohdest_rh12(:,:,:)  = 0._dp
     ohprod_rh12a(:,:,:) = 0._dp
     ohprod_rh9(:,:,:)   = 0._dp   
     ohprod_rh14a(:,:,:) = 0._dp
     ohdest_rh15(:,:,:)  = 0._dp
     ohprod_rh16(:,:,:)  = 0._dp  
     ohprod_rh18(:,:,:)  = 0._dp
     ohdest_rc5a(:,:,:)  = 0._dp 
     ohdest_rc5b(:,:,:)  = 0._dp  
     ohprod_rc9(:,:,:)   = 0._dp
     ohdest_rc10(:,:,:)  = 0._dp 
     ohprod_rc13a(:,:,:) = 0._dp 
     ohdest_rc6(:,:,:)   = 0._dp 
     ohprod_rc26(:,:,:)  = 0._dp 
     ohdest_rb14(:,:,:)  = 0._dp  
     ohprod_rb15(:,:,:)  = 0._dp
     ohprod_rb16(:,:,:)  = 0._dp
     ohprod_rb17(:,:,:)  = 0._dp  
     ohprod_rb22(:,:,:)  = 0._dp
     ohprod_rch1(:,:,:)  = 0._dp
     ohdest_rch3(:,:,:)  = 0._dp  
     ohdest_rch7(:,:,:)  = 0._dp
     ohprod_rch7a(:,:,:) = 0._dp
     ohdest_rch14(:,:,:) = 0._dp 
     ohprod_rch15(:,:,:) = 0._dp
     ohprod_rt06(:,:,:)  = 0._dp
     ohdest_rt07(:,:,:)  = 0._dp  
     ohprod_rt08(:,:,:)  = 0._dp 
     ohprod_rt11(:,:,:)  = 0._dp 
     ohprod_rt12(:,:,:)  = 0._dp  
     ohprod_rt17(:,:,:)  = 0._dp
     ohdest_rt18(:,:,:)  = 0._dp 
     ohdest_rt19(:,:,:)  = 0._dp  
     ohprod_rt20(:,:,:)  = 0._dp
     ohdest_rt26(:,:,:)  = 0._dp
     ohdest_ods10(:,:,:) = 0._dp 
     ohdest_ods11(:,:,:) = 0._dp 
     ohdest_ods13(:,:,:) = 0._dp 
     ohdest_ods14(:,:,:) = 0._dp
     ohdest_ods15(:,:,:) = 0._dp
     ohdest_ods16(:,:,:) = 0._dp
     ohdest_ods18(:,:,:) = 0._dp 
     ohdest_ods19(:,:,:) = 0._dp
     ohdest_rti01(:,:,:) = 0._dp
     ohprod_rti02(:,:,:) = 0._dp
     ohprod_rti08(:,:,:) = 0._dp
     ohdest_rti09(:,:,:) = 0._dp
     ohdest_rti10(:,:,:) = 0._dp
     ohprod_rti11(:,:,:) = 0._dp
     ohprod_rti17(:,:,:) = 0._dp
     ohdest_rti18(:,:,:) = 0._dp
     ohdest_rti19(:,:,:) = 0._dp
     ohdest_rti20(:,:,:) = 0._dp
     ohdest_rti22(:,:,:) = 0._dp
     ohdest_rti31(:,:,:) = 0._dp
     ohdest_rti33(:,:,:) = 0._dp
     ohdest_rti34(:,:,:) = 0._dp
     ohprod_rti35(:,:,:) = 0._dp
     ohprod_rti39(:,:,:) = 0._dp
     ohprod_rti44(:,:,:) = 0._dp

     ohprod_rn3_d(:,:,:)   = 0._dp 
     ohdest_rn8_d(:,:,:)   = 0._dp 
     ohprod_rn11_d(:,:,:)  = 0._dp  
     ohdest_rn12_d(:,:,:)  = 0._dp
     ohdest_rn13_d(:,:,:)  = 0._dp
     ohdest_rn21_d(:,:,:)  = 0._dp  
     ohdest_rn24a_d(:,:,:) = 0._dp
     ohprod_rn27_d(:,:,:)  = 0._dp
     ohdest_rn31_d(:,:,:)  = 0._dp
     ohprod_rh1_d(:,:,:)   = 0._dp
     ohdest_rh2_d(:,:,:)   = 0._dp
     ohdest_rh3_d(:,:,:)   = 0._dp   
     ohdest_rh4_d(:,:,:)   = 0._dp
     ohdest_rh13_d(:,:,:)  = 0._dp 
     ohdest_rh11_d(:,:,:)  = 0._dp   
     ohprod_rh5_d(:,:,:)   = 0._dp 
     ohprod_rh6_d(:,:,:)   = 0._dp 
     ohprod_rh8_d(:,:,:)   = 0._dp  
     ohdest_rh12_d(:,:,:)  = 0._dp
     ohprod_rh12a_d(:,:,:) = 0._dp
     ohprod_rh9_d(:,:,:)   = 0._dp   
     ohprod_rh14a_d(:,:,:) = 0._dp
     ohdest_rh15_d(:,:,:)  = 0._dp
     ohprod_rh16_d(:,:,:)  = 0._dp  
     ohprod_rh18_d(:,:,:)  = 0._dp
     ohdest_rc5a_d(:,:,:)  = 0._dp 
     ohdest_rc5b_d(:,:,:)  = 0._dp  
     ohprod_rc9_d(:,:,:)   = 0._dp
     ohdest_rc10_d(:,:,:)  = 0._dp 
     ohprod_rc13a_d(:,:,:) = 0._dp 
     ohdest_rc6_d(:,:,:)   = 0._dp 
     ohprod_rc26_d(:,:,:)  = 0._dp 
     ohdest_rb14_d(:,:,:)  = 0._dp  
     ohprod_rb15_d(:,:,:)  = 0._dp
     ohprod_rb16_d(:,:,:)  = 0._dp
     ohprod_rb17_d(:,:,:)  = 0._dp  
     ohprod_rb22_d(:,:,:)  = 0._dp
     ohprod_rch1_d(:,:,:)  = 0._dp
     ohdest_rch3_d(:,:,:)  = 0._dp  
     ohdest_rch7_d(:,:,:)  = 0._dp
     ohprod_rch7a_d(:,:,:) = 0._dp
     ohdest_rch14_d(:,:,:) = 0._dp 
     ohprod_rch15_d(:,:,:) = 0._dp
     ohprod_rt06_d(:,:,:)  = 0._dp
     ohdest_rt07_d(:,:,:)  = 0._dp  
     ohprod_rt08_d(:,:,:)  = 0._dp 
     ohprod_rt11_d(:,:,:)  = 0._dp 
     ohprod_rt12_d(:,:,:)  = 0._dp  
     ohprod_rt17_d(:,:,:)  = 0._dp
     ohdest_rt18_d(:,:,:)  = 0._dp 
     ohdest_rt19_d(:,:,:)  = 0._dp  
     ohprod_rt20_d(:,:,:)  = 0._dp
     ohdest_rt26_d(:,:,:)  = 0._dp
     ohdest_ods10_d(:,:,:) = 0._dp 
     ohdest_ods11_d(:,:,:) = 0._dp 
     ohdest_ods13_d(:,:,:) = 0._dp 
     ohdest_ods14_d(:,:,:) = 0._dp
     ohdest_ods15_d(:,:,:) = 0._dp
     ohdest_ods16_d(:,:,:) = 0._dp
     ohdest_ods18_d(:,:,:) = 0._dp 
     ohdest_ods19_d(:,:,:) = 0._dp
     ohdest_rti01_d(:,:,:) = 0._dp
     ohprod_rti02_d(:,:,:) = 0._dp
     ohprod_rti08_d(:,:,:) = 0._dp
     ohdest_rti09_d(:,:,:) = 0._dp
     ohdest_rti10_d(:,:,:) = 0._dp
     ohprod_rti11_d(:,:,:) = 0._dp
     ohprod_rti17_d(:,:,:) = 0._dp
     ohdest_rti18_d(:,:,:) = 0._dp
     ohdest_rti19_d(:,:,:) = 0._dp
     ohdest_rti20_d(:,:,:) = 0._dp
     ohdest_rti22_d(:,:,:) = 0._dp
     ohdest_rti31_d(:,:,:) = 0._dp
     ohdest_rti33_d(:,:,:) = 0._dp
     ohdest_rti34_d(:,:,:) = 0._dp
     ohprod_rti35_d(:,:,:) = 0._dp
     ohprod_rti39_d(:,:,:) = 0._dp
     ohprod_rti44_d(:,:,:) = 0._dp

     ho2dest_rn3(:,:,:) = 0._dp 
     ho2prod_rn20(:,:,:) = 0._dp
     ho2dest_rn22(:,:,:) = 0._dp 
     ho2prod_rn22a(:,:,:) = 0._dp 
     ho2dest_rn27(:,:,:) = 0._dp 
     ho2prod_rh3(:,:,:) = 0._dp  
     ho2dest_rh4(:,:,:) = 0._dp 
     ho2dest_rh5(:,:,:) = 0._dp
     ho2dest_rh6(:,:,:) = 0._dp 
     ho2dest_rh7(:,:,:) = 0._dp
     ho2dest_rh7a(:,:,:) = 0._dp 
     ho2prod_rh12(:,:,:) = 0._dp
     ho2prod_rh12a(:,:,:) = 0._dp 
     ho2prod_rh10(:,:,:) = 0._dp
     ho2dest_rh14(:,:,:) = 0._dp 
     ho2dest_rh14a(:,:,:) = 0._dp
     ho2prod_rc5a(:,:,:) = 0._dp 
     ho2dest_rc11(:,:,:) = 0._dp
     ho2dest_rc13(:,:,:) = 0._dp 
     ho2dest_rc13a(:,:,:) = 0._dp 
     ho2prod_rc14(:,:,:) = 0._dp 
     ho2dest_rb02(:,:,:) = 0._dp
     ho2prod_rb05(:,:,:) = 0._dp 
     ho2dest_rb13(:,:,:) = 0._dp
     ho2prod_rch6(:,:,:) = 0._dp 
     ho2prod_rch10(:,:,:) = 0._dp
     ho2dest_rch13(:,:,:) = 0._dp 
     ho2prod_rch21(:,:,:) = 0._dp
     ho2prod_rt07(:,:,:) = 0._dp 
     ho2dest_rt08(:,:,:) = 0._dp 
     ho2dest_rt09(:,:,:) = 0._dp 
     ho2prod_rt10(:,:,:) = 0._dp 
     ho2prod_rt19(:,:,:) = 0._dp 
     ho2prod_rti02(:,:,:) = 0._dp
     ho2prod_rti04(:,:,:) = 0._dp
     ho2dest_rti06(:,:,:) = 0._dp
     ho2prod_rti07(:,:,:) = 0._dp
     ho2prod_rti11(:,:,:) = 0._dp
     ho2prod_rti12(:,:,:) = 0._dp
     ho2dest_rti13(:,:,:) = 0._dp
     ho2prod_rti14(:,:,:) = 0._dp
     ho2prod_rti19(:,:,:) = 0._dp
     ho2dest_rti23(:,:,:) = 0._dp
     ho2dest_rti24(:,:,:) = 0._dp
     ho2prod_rti27(:,:,:) = 0._dp
     ho2prod_rti35(:,:,:) = 0._dp
     ho2prod_rti36(:,:,:) = 0._dp
     ho2prod_rti37(:,:,:) = 0._dp
     ho2prod_rti39(:,:,:) = 0._dp
     ho2prod_rti40(:,:,:) = 0._dp
     ho2prod_rti41(:,:,:) = 0._dp
     ho2prod_rti42(:,:,:) = 0._dp

     ho2dest_rn3_d(:,:,:) = 0._dp 
     ho2prod_rn20_d(:,:,:) = 0._dp
     ho2dest_rn22_d(:,:,:) = 0._dp 
     ho2prod_rn22a_d(:,:,:) = 0._dp 
     ho2dest_rn27_d(:,:,:) = 0._dp 
     ho2prod_rh3_d(:,:,:) = 0._dp  
     ho2dest_rh4_d(:,:,:) = 0._dp 
     ho2dest_rh5_d(:,:,:) = 0._dp
     ho2dest_rh6_d(:,:,:) = 0._dp 
     ho2dest_rh7_d(:,:,:) = 0._dp
     ho2dest_rh7a_d(:,:,:) = 0._dp 
     ho2prod_rh12_d(:,:,:) = 0._dp
     ho2prod_rh12a_d(:,:,:) = 0._dp 
     ho2prod_rh10_d(:,:,:) = 0._dp
     ho2dest_rh14_d(:,:,:) = 0._dp 
     ho2dest_rh14a_d(:,:,:) = 0._dp
     ho2prod_rc5a_d(:,:,:) = 0._dp 
     ho2dest_rc11_d(:,:,:) = 0._dp
     ho2dest_rc13_d(:,:,:) = 0._dp 
     ho2dest_rc13a_d(:,:,:) = 0._dp 
     ho2prod_rc14_d(:,:,:) = 0._dp 
     ho2dest_rb02_d(:,:,:) = 0._dp
     ho2prod_rb05_d(:,:,:) = 0._dp 
     ho2dest_rb13_d(:,:,:) = 0._dp
     ho2prod_rch6_d(:,:,:) = 0._dp 
     ho2prod_rch10_d(:,:,:) = 0._dp
     ho2dest_rch13_d(:,:,:) = 0._dp 
     ho2prod_rch21_d(:,:,:) = 0._dp
     ho2prod_rt07_d(:,:,:) = 0._dp 
     ho2dest_rt08_d(:,:,:) = 0._dp 
     ho2dest_rt09_d(:,:,:) = 0._dp 
     ho2prod_rt10_d(:,:,:) = 0._dp 
     ho2prod_rt19_d(:,:,:) = 0._dp 
     ho2prod_rti02_d(:,:,:) = 0._dp
     ho2prod_rti04_d(:,:,:) = 0._dp
     ho2dest_rti06_d(:,:,:) = 0._dp
     ho2prod_rti07_d(:,:,:) = 0._dp
     ho2prod_rti11_d(:,:,:) = 0._dp
     ho2prod_rti12_d(:,:,:) = 0._dp
     ho2dest_rti13_d(:,:,:) = 0._dp
     ho2prod_rti14_d(:,:,:) = 0._dp
     ho2prod_rti19_d(:,:,:) = 0._dp
     ho2dest_rti23_d(:,:,:) = 0._dp
     ho2dest_rti24_d(:,:,:) = 0._dp
     ho2prod_rti27_d(:,:,:) = 0._dp
     ho2prod_rti35_d(:,:,:) = 0._dp
     ho2prod_rti36_d(:,:,:) = 0._dp
     ho2prod_rti37_d(:,:,:) = 0._dp
     ho2prod_rti39_d(:,:,:) = 0._dp
     ho2prod_rti40_d(:,:,:) = 0._dp
     ho2prod_rti41_d(:,:,:) = 0._dp
     ho2prod_rti42_d(:,:,:) = 0._dp

     ! eth_as_oh-
  ENDIF
 
END SUBROUTINE init_stream_ch4

!---------------------------------------------------------------------------------

SUBROUTINE accumulate_stream_ch4 (kproma, krow)

  USE mo_socol_namelist, ONLY: lch4_wetland

  ! This subroutine accumulates the current value of a variable at every time 
  ! step, such that finally monthly streams can be calculated.

  ! *accumulate_stream_ch4* is called from *call_diagn*, 
  ! src/call_submodels

  ! Scalar arguments
  INTEGER, INTENT(in) :: kproma, krow


  ! Executable statements:

  ch4descl(1:kproma,:,krow)  = ch4descl(1:kproma,:,krow) + delta_time * &
       ch4descl_d(1:kproma,:,krow)
  ch4desod_1(1:kproma,:,krow)  = ch4desod_1(1:kproma,:,krow) + delta_time * &
       ch4desod_1_d(1:kproma,:,krow)
  ch4desod_2(1:kproma,:,krow)  = ch4desod_2(1:kproma,:,krow) + delta_time * &
       ch4desod_2_d(1:kproma,:,krow)
  ch4desoh(1:kproma,:,krow)  = ch4desoh(1:kproma,:,krow) + delta_time * &
       ch4desoh_d(1:kproma,:,krow)
  ch4dest01(1:kproma,:,krow)  = ch4dest01(1:kproma,:,krow) + delta_time * &
       ch4dest01_d(1:kproma,:,krow)

  IF (lch4_wetland) &
       ch4emswet(1:kproma,krow) = ch4emswet(1:kproma,krow) + delta_time * &
       ch4emswet_d(1:kproma,krow)

  ch4_emflux(1:kproma,krow) = ch4_emflux(1:kproma,krow) + delta_time * &
       ch4_emflux_d(1:kproma,krow)

  ! eth_as_oh+
  ohprod_rn3(1:kproma,:,krow) = ohprod_rn3(1:kproma,:,krow) + delta_time * &
       ohprod_rn3_d(1:kproma,:,krow)
  ohdest_rn8(1:kproma,:,krow) = ohdest_rn8(1:kproma,:,krow) + delta_time * &
       ohdest_rn8_d(1:kproma,:,krow)
  ohprod_rn11(1:kproma,:,krow) = ohprod_rn11(1:kproma,:,krow) + delta_time * &
       ohprod_rn11_d(1:kproma,:,krow)
  ohdest_rn12(1:kproma,:,krow) = ohdest_rn12(1:kproma,:,krow) + delta_time * &
       ohdest_rn12_d(1:kproma,:,krow)
  ohdest_rn13(1:kproma,:,krow) = ohdest_rn13(1:kproma,:,krow) + delta_time * &
       ohdest_rn13_d(1:kproma,:,krow)
  ohdest_rn21(1:kproma,:,krow) = ohdest_rn21(1:kproma,:,krow) + delta_time * &
       ohdest_rn21_d(1:kproma,:,krow)
  ohdest_rn24a(1:kproma,:,krow) = ohdest_rn24a(1:kproma,:,krow) + delta_time * &
       ohdest_rn24a_d(1:kproma,:,krow)
  ohprod_rn27(1:kproma,:,krow) = ohprod_rn27(1:kproma,:,krow) + delta_time * &
       ohprod_rn27_d(1:kproma,:,krow)
  ohdest_rn31(1:kproma,:,krow) = ohdest_rn31(1:kproma,:,krow) + delta_time * &
       ohdest_rn31_d(1:kproma,:,krow)
  ohprod_rh1(1:kproma,:,krow) = ohprod_rh1(1:kproma,:,krow) + delta_time * &
       ohprod_rh1_d(1:kproma,:,krow)
  ohdest_rh2(1:kproma,:,krow) = ohdest_rh2(1:kproma,:,krow) + delta_time * &
       ohdest_rh2_d(1:kproma,:,krow)
  ohdest_rh3(1:kproma,:,krow) = ohdest_rh3(1:kproma,:,krow) + delta_time * &
       ohdest_rh3_d(1:kproma,:,krow)
  ohdest_rh4(1:kproma,:,krow) = ohdest_rh4(1:kproma,:,krow) + delta_time * &
       ohdest_rh4_d(1:kproma,:,krow)
  ohdest_rh13(1:kproma,:,krow) = ohdest_rh13(1:kproma,:,krow) + delta_time * &
       ohdest_rh13_d(1:kproma,:,krow)
  ohdest_rh11(1:kproma,:,krow) = ohdest_rh11(1:kproma,:,krow) + delta_time * &
       ohdest_rh11_d(1:kproma,:,krow)
  ohprod_rh5(1:kproma,:,krow) = ohprod_rh5(1:kproma,:,krow) + delta_time * &
       ohprod_rh5_d(1:kproma,:,krow)
  ohprod_rh6(1:kproma,:,krow) = ohprod_rh6(1:kproma,:,krow) + delta_time * &
       ohprod_rh6_d(1:kproma,:,krow)
  ohprod_rh8(1:kproma,:,krow) = ohprod_rh8(1:kproma,:,krow) + delta_time * &
       ohprod_rh8_d(1:kproma,:,krow)
  ohdest_rh12(1:kproma,:,krow) = ohdest_rh12(1:kproma,:,krow) + delta_time * &
       ohdest_rh12_d(1:kproma,:,krow)
  ohprod_rh12a(1:kproma,:,krow) = ohprod_rh12a(1:kproma,:,krow) + delta_time * &
       ohprod_rh12a_d(1:kproma,:,krow)
  ohprod_rh9(1:kproma,:,krow) = ohprod_rh9(1:kproma,:,krow) + delta_time * &
       ohprod_rh9_d(1:kproma,:,krow)
  ohprod_rh14a(1:kproma,:,krow) = ohprod_rh14a(1:kproma,:,krow) + delta_time * &
       ohprod_rh14a_d(1:kproma,:,krow)
  ohdest_rh15(1:kproma,:,krow) = ohdest_rh15(1:kproma,:,krow) + delta_time * &
       ohdest_rh15_d(1:kproma,:,krow)
  ohprod_rh16(1:kproma,:,krow) = ohprod_rh16(1:kproma,:,krow) + delta_time * &
       ohprod_rh16_d(1:kproma,:,krow)
  ohprod_rh18(1:kproma,:,krow) = ohprod_rh18(1:kproma,:,krow) + delta_time * &
       ohprod_rh18_d(1:kproma,:,krow)
  ohdest_rc5a(1:kproma,:,krow) = ohdest_rc5a(1:kproma,:,krow) + delta_time * &
       ohdest_rc5a_d(1:kproma,:,krow)
  ohdest_rc5b(1:kproma,:,krow) = ohdest_rc5b(1:kproma,:,krow) + delta_time * &
       ohdest_rc5b_d(1:kproma,:,krow)
  ohprod_rc9(1:kproma,:,krow) = ohprod_rc9(1:kproma,:,krow) + delta_time * &
       ohprod_rc9_d(1:kproma,:,krow)
  ohdest_rc10(1:kproma,:,krow) = ohdest_rc10(1:kproma,:,krow) + delta_time * &
       ohdest_rc10_d(1:kproma,:,krow)
  ohprod_rc13a(1:kproma,:,krow) = ohprod_rc13a(1:kproma,:,krow) + delta_time * &
       ohprod_rc13a_d(1:kproma,:,krow)
  ohdest_rc6(1:kproma,:,krow) = ohdest_rc6(1:kproma,:,krow) + delta_time * &
       ohdest_rc6_d(1:kproma,:,krow)
  ohprod_rc26(1:kproma,:,krow) = ohprod_rc26(1:kproma,:,krow) + delta_time * &
       ohprod_rc26_d(1:kproma,:,krow)
  ohdest_rb14(1:kproma,:,krow) = ohdest_rb14(1:kproma,:,krow) + delta_time * &
       ohdest_rb14_d(1:kproma,:,krow)
  ohprod_rb15(1:kproma,:,krow) = ohprod_rb15(1:kproma,:,krow) + delta_time * &
       ohprod_rb15_d(1:kproma,:,krow)
  ohprod_rb16(1:kproma,:,krow) = ohprod_rb16(1:kproma,:,krow) + delta_time * &
       ohprod_rb16_d(1:kproma,:,krow)
  ohprod_rb17(1:kproma,:,krow) = ohprod_rb17(1:kproma,:,krow) + delta_time * &
       ohprod_rb17_d(1:kproma,:,krow)
  ohprod_rb22(1:kproma,:,krow) = ohprod_rb22(1:kproma,:,krow) + delta_time * &
       ohprod_rb22_d(1:kproma,:,krow)
  ohprod_rch1(1:kproma,:,krow) = ohprod_rch1(1:kproma,:,krow) + delta_time * &
       ohprod_rch1_d(1:kproma,:,krow)
  ohdest_rch3(1:kproma,:,krow) = ohdest_rch3(1:kproma,:,krow) + delta_time * &
       ohdest_rch3_d(1:kproma,:,krow)
  ohdest_rch7(1:kproma,:,krow) = ohdest_rch7(1:kproma,:,krow) + delta_time * &
       ohdest_rch7_d(1:kproma,:,krow)
  ohprod_rch7a(1:kproma,:,krow) = ohprod_rch7a(1:kproma,:,krow) + delta_time * &
       ohprod_rch7a_d(1:kproma,:,krow)
  ohdest_rch14(1:kproma,:,krow) = ohdest_rch14(1:kproma,:,krow) + delta_time * &
       ohdest_rch14_d(1:kproma,:,krow)
  ohprod_rch15(1:kproma,:,krow) = ohprod_rch15(1:kproma,:,krow) + delta_time * &
       ohprod_rch15_d(1:kproma,:,krow)
  ohprod_rt06(1:kproma,:,krow) = ohprod_rt06(1:kproma,:,krow) + delta_time * &
       ohprod_rt06_d(1:kproma,:,krow)
  ohdest_rt07(1:kproma,:,krow) = ohdest_rt07(1:kproma,:,krow) + delta_time * &
       ohdest_rt07_d(1:kproma,:,krow)
  ohprod_rt08(1:kproma,:,krow) = ohprod_rt08(1:kproma,:,krow) + delta_time * &
       ohprod_rt08_d(1:kproma,:,krow)
  ohprod_rt11(1:kproma,:,krow) = ohprod_rt11(1:kproma,:,krow) + delta_time * &
       ohprod_rt11_d(1:kproma,:,krow)
  ohprod_rt12(1:kproma,:,krow) = ohprod_rt12(1:kproma,:,krow) + delta_time * &
       ohprod_rt12_d(1:kproma,:,krow)
  ohprod_rt17(1:kproma,:,krow) = ohprod_rt17(1:kproma,:,krow) + delta_time * &
       ohprod_rt17_d(1:kproma,:,krow)
  ohdest_rt18(1:kproma,:,krow) = ohdest_rt18(1:kproma,:,krow) + delta_time * &
       ohdest_rt18_d(1:kproma,:,krow)
  ohdest_rt19(1:kproma,:,krow) = ohdest_rt19(1:kproma,:,krow) + delta_time * &
       ohdest_rt19_d(1:kproma,:,krow)
  ohprod_rt20(1:kproma,:,krow) = ohprod_rt20(1:kproma,:,krow) + delta_time * &
       ohprod_rt20_d(1:kproma,:,krow)
  ohdest_rt26(1:kproma,:,krow) = ohdest_rt26(1:kproma,:,krow) + delta_time * &
       ohdest_rt26_d(1:kproma,:,krow)
  ohdest_ods10(1:kproma,:,krow) = ohdest_ods10(1:kproma,:,krow) + delta_time * &
       ohdest_ods10_d(1:kproma,:,krow)
  ohdest_ods11(1:kproma,:,krow) = ohdest_ods11(1:kproma,:,krow) + delta_time * &
       ohdest_ods11_d(1:kproma,:,krow)
  ohdest_ods13(1:kproma,:,krow) = ohdest_ods13(1:kproma,:,krow) + delta_time * &
       ohdest_ods13_d(1:kproma,:,krow)
  ohdest_ods14(1:kproma,:,krow) = ohdest_ods14(1:kproma,:,krow) + delta_time * &
       ohdest_ods14_d(1:kproma,:,krow)
  ohdest_ods15(1:kproma,:,krow) = ohdest_ods15(1:kproma,:,krow) + delta_time * &
       ohdest_ods15_d(1:kproma,:,krow)
  ohdest_ods16(1:kproma,:,krow) = ohdest_ods16(1:kproma,:,krow) + delta_time * &
       ohdest_ods16_d(1:kproma,:,krow)
  ohdest_ods18(1:kproma,:,krow) = ohdest_ods18(1:kproma,:,krow) + delta_time * &
       ohdest_ods18_d(1:kproma,:,krow)
  ohdest_ods19(1:kproma,:,krow) = ohdest_ods19(1:kproma,:,krow) + delta_time * &
       ohdest_ods19_d(1:kproma,:,krow)
  ohdest_rti01(1:kproma,:,krow) = ohdest_rti01(1:kproma,:,krow) + delta_time * &
       ohdest_rti01_d(1:kproma,:,krow)
  ohprod_rti02(1:kproma,:,krow) = ohprod_rti02(1:kproma,:,krow) + delta_time * &
       ohprod_rti02_d(1:kproma,:,krow)
  ohprod_rti08(1:kproma,:,krow) = ohprod_rti08(1:kproma,:,krow) + delta_time * &
       ohprod_rti08_d(1:kproma,:,krow)
  ohdest_rti09(1:kproma,:,krow) = ohdest_rti09(1:kproma,:,krow) + delta_time * &
       ohdest_rti09_d(1:kproma,:,krow)
  ohdest_rti10(1:kproma,:,krow) = ohdest_rti10(1:kproma,:,krow) + delta_time * &
       ohdest_rti10_d(1:kproma,:,krow)
  ohprod_rti11(1:kproma,:,krow) = ohprod_rti11(1:kproma,:,krow) + delta_time * &
       ohprod_rti11_d(1:kproma,:,krow)
  ohprod_rti17(1:kproma,:,krow) = ohprod_rti17(1:kproma,:,krow) + delta_time * &
       ohprod_rti17_d(1:kproma,:,krow) 
  ohdest_rti18(1:kproma,:,krow) = ohdest_rti18(1:kproma,:,krow) + delta_time * &
       ohdest_rti18_d(1:kproma,:,krow)
  ohdest_rti19(1:kproma,:,krow) = ohdest_rti19(1:kproma,:,krow) + delta_time * &
       ohdest_rti19_d(1:kproma,:,krow)
  ohdest_rti20(1:kproma,:,krow) = ohdest_rti20(1:kproma,:,krow) + delta_time * &
       ohdest_rti20_d(1:kproma,:,krow)
  ohdest_rti22(1:kproma,:,krow) = ohdest_rti22(1:kproma,:,krow) + delta_time * &
       ohdest_rti22_d(1:kproma,:,krow)
  ohdest_rti31(1:kproma,:,krow) = ohdest_rti31(1:kproma,:,krow) + delta_time * &
       ohdest_rti31_d(1:kproma,:,krow)
  ohdest_rti33(1:kproma,:,krow) = ohdest_rti33(1:kproma,:,krow) + delta_time * &
       ohdest_rti33_d(1:kproma,:,krow)
  ohdest_rti34(1:kproma,:,krow) = ohdest_rti34(1:kproma,:,krow) + delta_time * &
       ohdest_rti34_d(1:kproma,:,krow)
  ohprod_rti35(1:kproma,:,krow) = ohprod_rti35(1:kproma,:,krow) + delta_time * &
       ohprod_rti35_d(1:kproma,:,krow) 
  ohprod_rti39(1:kproma,:,krow) = ohprod_rti39(1:kproma,:,krow) + delta_time * &
       ohprod_rti39_d(1:kproma,:,krow)
  ohprod_rti44(1:kproma,:,krow) = ohprod_rti44(1:kproma,:,krow) + delta_time * &
       ohprod_rti44_d(1:kproma,:,krow)

  ho2dest_rn3(1:kproma,:,krow) = ho2dest_rn3(1:kproma,:,krow) + delta_time * &
       ho2dest_rn3_d(1:kproma,:,krow)
  ho2prod_rn20(1:kproma,:,krow) = ho2prod_rn20(1:kproma,:,krow) + delta_time * &
       ho2prod_rn20_d(1:kproma,:,krow)
  ho2dest_rn22(1:kproma,:,krow) = ho2dest_rn22(1:kproma,:,krow) + delta_time * &
       ho2dest_rn22_d(1:kproma,:,krow)
  ho2prod_rn22a(1:kproma,:,krow) = ho2prod_rn22a(1:kproma,:,krow) + delta_time * &
       ho2prod_rn22a_d(1:kproma,:,krow)
  ho2dest_rn27(1:kproma,:,krow) = ho2dest_rn27(1:kproma,:,krow) + delta_time * &
       ho2dest_rn27_d(1:kproma,:,krow)
  ho2prod_rh3(1:kproma,:,krow) = ho2prod_rh3(1:kproma,:,krow) + delta_time * &
       ho2prod_rh3_d(1:kproma,:,krow)
  ho2dest_rh4(1:kproma,:,krow) = ho2dest_rh4(1:kproma,:,krow) + delta_time * &
       ho2dest_rh4_d(1:kproma,:,krow)
  ho2dest_rh5(1:kproma,:,krow) = ho2dest_rh5(1:kproma,:,krow) + delta_time * &
       ho2dest_rh5_d(1:kproma,:,krow)
  ho2dest_rh6(1:kproma,:,krow) = ho2dest_rh6(1:kproma,:,krow) + delta_time * &
       ho2dest_rh6_d(1:kproma,:,krow)
  ho2dest_rh7(1:kproma,:,krow) = ho2dest_rh7(1:kproma,:,krow) + delta_time * &
       ho2dest_rh7_d(1:kproma,:,krow)
  ho2dest_rh7a(1:kproma,:,krow) = ho2dest_rh7a(1:kproma,:,krow) + delta_time * &
       ho2dest_rh7a_d(1:kproma,:,krow)
  ho2prod_rh12(1:kproma,:,krow) = ho2prod_rh12(1:kproma,:,krow) + delta_time * &
       ho2prod_rh12_d(1:kproma,:,krow)
  ho2prod_rh12a(1:kproma,:,krow) = ho2prod_rh12a(1:kproma,:,krow) + delta_time * &
       ho2prod_rh12a_d(1:kproma,:,krow)
  ho2prod_rh10(1:kproma,:,krow) = ho2prod_rh10(1:kproma,:,krow) + delta_time * &
       ho2prod_rh10_d(1:kproma,:,krow)
  ho2dest_rh14(1:kproma,:,krow) = ho2dest_rh14(1:kproma,:,krow) + delta_time * &
       ho2dest_rh14_d(1:kproma,:,krow)
  ho2dest_rh14a(1:kproma,:,krow) = ho2dest_rh14a(1:kproma,:,krow) + delta_time * &
       ho2dest_rh14a_d(1:kproma,:,krow)
  ho2prod_rc5a(1:kproma,:,krow) = ho2prod_rc5a(1:kproma,:,krow) + delta_time * &
       ho2prod_rc5a_d(1:kproma,:,krow)
  ho2dest_rc11(1:kproma,:,krow) = ho2dest_rc11(1:kproma,:,krow) + delta_time * &
       ho2dest_rc11_d(1:kproma,:,krow)
  ho2dest_rc13(1:kproma,:,krow) = ho2dest_rc13(1:kproma,:,krow) + delta_time * &
       ho2dest_rc13_d(1:kproma,:,krow)
  ho2dest_rc13a(1:kproma,:,krow) = ho2dest_rc13a(1:kproma,:,krow) + delta_time * &
       ho2dest_rc13a_d(1:kproma,:,krow)
  ho2prod_rc14(1:kproma,:,krow) = ho2prod_rc14(1:kproma,:,krow) + delta_time * &
       ho2prod_rc14_d(1:kproma,:,krow)
  ho2dest_rb02(1:kproma,:,krow) = ho2dest_rb02(1:kproma,:,krow) + delta_time * &
       ho2dest_rb02_d(1:kproma,:,krow)
  ho2prod_rb05(1:kproma,:,krow) = ho2prod_rb05(1:kproma,:,krow) + delta_time * &
       ho2prod_rb05_d(1:kproma,:,krow)
  ho2dest_rb13(1:kproma,:,krow) = ho2dest_rb13(1:kproma,:,krow) + delta_time * &
       ho2dest_rb13_d(1:kproma,:,krow)
  ho2prod_rch6(1:kproma,:,krow) = ho2prod_rch6(1:kproma,:,krow) + delta_time * &
       ho2prod_rch6_d(1:kproma,:,krow)
  ho2prod_rch10(1:kproma,:,krow) = ho2prod_rch10(1:kproma,:,krow) + delta_time * &
       ho2prod_rch10_d(1:kproma,:,krow)
  ho2dest_rch13(1:kproma,:,krow) = ho2dest_rch13(1:kproma,:,krow) + delta_time * &
       ho2dest_rch13_d(1:kproma,:,krow)
  ho2prod_rch21(1:kproma,:,krow) = ho2prod_rch21(1:kproma,:,krow) + delta_time * &
       ho2prod_rch21_d(1:kproma,:,krow)
  ho2prod_rt07(1:kproma,:,krow) = ho2prod_rt07(1:kproma,:,krow) + delta_time * &
       ho2prod_rt07_d(1:kproma,:,krow)
  ho2dest_rt08(1:kproma,:,krow) = ho2dest_rt08(1:kproma,:,krow) + delta_time * &
       ho2dest_rt08_d(1:kproma,:,krow)
  ho2dest_rt09(1:kproma,:,krow) = ho2dest_rt09(1:kproma,:,krow) + delta_time * &
       ho2dest_rt09_d(1:kproma,:,krow)
  ho2prod_rt10(1:kproma,:,krow) = ho2prod_rt10(1:kproma,:,krow) + delta_time * &
       ho2prod_rt10_d(1:kproma,:,krow)
  ho2prod_rt19(1:kproma,:,krow) = ho2prod_rt19(1:kproma,:,krow) + delta_time * &
       ho2prod_rt19_d(1:kproma,:,krow)
  ho2prod_rti02(1:kproma,:,krow) = ho2prod_rti02(1:kproma,:,krow) + delta_time * &
       ho2prod_rti02_d(1:kproma,:,krow)
  ho2prod_rti04(1:kproma,:,krow) = ho2prod_rti04(1:kproma,:,krow) + delta_time * &
       ho2prod_rti04_d(1:kproma,:,krow)
  ho2dest_rti06(1:kproma,:,krow) = ho2dest_rti06(1:kproma,:,krow) + delta_time * &
       ho2dest_rti06_d(1:kproma,:,krow)    
  ho2prod_rti07(1:kproma,:,krow) = ho2prod_rti07(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti07_d(1:kproma,:,krow)
  ho2prod_rti11(1:kproma,:,krow) = ho2prod_rti11(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti11_d(1:kproma,:,krow)
  ho2prod_rti12(1:kproma,:,krow) = ho2prod_rti12(1:kproma,:,krow) + delta_time * &   
       ho2prod_rti12_d(1:kproma,:,krow)
  ho2dest_rti13(1:kproma,:,krow) = ho2dest_rti13(1:kproma,:,krow) + delta_time * &
       ho2dest_rti13_d(1:kproma,:,krow)
  ho2prod_rti14(1:kproma,:,krow) = ho2prod_rti14(1:kproma,:,krow) + delta_time * &  
       ho2prod_rti14_d(1:kproma,:,krow)
  ho2prod_rti19(1:kproma,:,krow) = ho2prod_rti19(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti19_d(1:kproma,:,krow)
  ho2dest_rti23(1:kproma,:,krow) = ho2dest_rti23(1:kproma,:,krow) + delta_time * &
       ho2dest_rti23_d(1:kproma,:,krow)
  ho2dest_rti24(1:kproma,:,krow) = ho2dest_rti24(1:kproma,:,krow) + delta_time * &
       ho2dest_rti24_d(1:kproma,:,krow)
  ho2prod_rti27(1:kproma,:,krow) = ho2prod_rti27(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti27_d(1:kproma,:,krow)
  ho2prod_rti35(1:kproma,:,krow) = ho2prod_rti35(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti35_d(1:kproma,:,krow)
  ho2prod_rti36(1:kproma,:,krow) = ho2prod_rti36(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti36_d(1:kproma,:,krow)
  ho2prod_rti37(1:kproma,:,krow) = ho2prod_rti37(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti37_d(1:kproma,:,krow)
  ho2prod_rti39(1:kproma,:,krow) = ho2prod_rti39(1:kproma,:,krow) + delta_time * &  
       ho2prod_rti39_d(1:kproma,:,krow)
  ho2prod_rti40(1:kproma,:,krow) = ho2prod_rti40(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti40_d(1:kproma,:,krow)
  ho2prod_rti41(1:kproma,:,krow) = ho2prod_rti41(1:kproma,:,krow) + delta_time * & 
       ho2prod_rti41_d(1:kproma,:,krow)  
  ho2prod_rti42(1:kproma,:,krow) = ho2prod_rti42(1:kproma,:,krow) + delta_time * &  
       ho2prod_rti42_d(1:kproma,:,krow)
  ! eth_as_oh-

END SUBROUTINE accumulate_stream_ch4
          
END MODULE mo_socol_ch4_streams

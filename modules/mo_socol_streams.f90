MODULE mo_socol_streams

  ! *mo_socol_streams* contains subroutines to define streams of chemical 
  ! species as monthly means for SOCOL.
  !
  ! Martin Schraner, ETH Zurich, September 2008 

  USE mo_doctor,        ONLY: nout
  USE mo_kind,          ONLY: dp      
  USE mo_linked_list,   ONLY: HYBRID, NETCDF
  USE mo_memory_base,   ONLY: new_stream, add_stream_element,  &
                              default_stream_setting, add_stream_reference, &
                              delete_stream, t_stream
  USE mo_socol_namelist, ONLY: interactivelnox, &! eth_as_03062013
                               lsynth  !CCMI tracers
  USE mo_socol_tracers
  USE mo_socol_o3orig, ONLY: n_trac_orig, itrac_o3orig
  USE mo_time_control,  ONLY: l_trigrad, trigrad, delta_time, lstart
  USE mo_time_event,    ONLY: io_time_event, TIME_INC_SECONDS, &
                              TIME_INC_MINUTES, TIME_INC_HOURS , &
                              TIME_INC_DAYS, TIME_INC_MONTHS, TRIG_LAST
  USE mo_socol_synth_tracers, ONLY: idt_NH_5, idt_NH_50,  idt_NH_50W, & !CCMI tracers
                                    idt_ST80_25, idt_CO_25, idt_CO_50

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_stream_chem_m  ! construct the chem_m stream
  PUBLIC :: construct_stream_chem2   ! construct the chem2 stream
  PUBLIC :: construct_stream_chem2d   ! construct the chem2 strea
  PUBLIC :: destruct_stream_chem_m   ! destruct the chem_m stream
  PUBLIC :: destruct_stream_chem2    ! destruct the chem2 stream
  PUBLIC :: destruct_stream_chem2d   ! destruct the chem2 stream
  PUBLIC :: init_stream_chem_m       ! initialize the chem_m stream
  PUBLIC :: init_stream_chem2        ! initialize the chem2 stream
  PUBLIC :: init_stream_chem2d        ! initialize the chem2d stream
  PUBLIC :: accumulate_stream_chem_m ! accumulation of streams
  PUBLIC :: accumulate_stream_chem2 ! accumulation of streams
  
  PUBLIC :: construct_stream_nudg_2h  ! construct the nudg_2h stream
  PUBLIC :: destruct_stream_nudg_2h   ! destruct the nudg_2h stream
  PUBLIC :: init_stream_nudg_2h       ! initialize the nudg_2h stream

  TYPE (t_stream), PUBLIC, POINTER :: chem_m    !the chem_m stream
  TYPE (t_stream), PUBLIC, POINTER :: chem2     !the chem2 stream
  TYPE (t_stream), PUBLIC, POINTER :: chem2d     !the chem2 stream
  TYPE (t_stream), PUBLIC, POINTER :: nudg_2h   !the nudg_2h stream
  
  ! chem_m (=chem1) stream
  REAL(dp), PUBLIC, POINTER :: o3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: o_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: no_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: no2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hno3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: no3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: n2o5_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hno4_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ho2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: clno3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: clo_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: n_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: oh_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: h_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cl_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hocl_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: n2o_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: co_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hcl_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: h2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: h2o2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: h2o_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cl2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cl2o2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3o2h_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch2o_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: bro_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: brno3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: brcl_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hbr_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hobr_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3o2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3o_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hco_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: br_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: od_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3co3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: pan_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: c5h8_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: f11_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: f12_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cbrf3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cfc113_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cfc114_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cfc115_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccl4_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3ccl3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hcfc22_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hcfc141b_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hcfc142b_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: h1211_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3br_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch3cl_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hcfc21_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hcfc123_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: h2402_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: chbr3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch2br2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: linage_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: idealage_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: NH_5_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: NH_50_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: NH_50W_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ST80_25_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: CO_25_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: CO_50_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: O3ONPBL_1_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: O3ONMBL_2_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3OTRBL_3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: O3OSMBL_4_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3OSPBL_5_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3ONPFT_6_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3ONMFT_7_m(:,:,:)   
  REAL(dp), PUBLIC, POINTER :: O3OTRFT_8_m(:,:,:)   
  REAL(dp), PUBLIC, POINTER :: O3OSMFT_9_m(:,:,:)   
  REAL(dp), PUBLIC, POINTER :: O3OSPFT_10_m(:,:,:)  
  REAL(dp), PUBLIC, POINTER :: O3ONPLS_11_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3ONMLS_12_m(:,:,:)  
  REAL(dp), PUBLIC, POINTER :: O3OTRLS_13_m(:,:,:)  
  REAL(dp), PUBLIC, POINTER :: O3OTRMS_14_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3OSMLS_15_m(:,:,:)  
  REAL(dp), PUBLIC, POINTER :: O3OSPLS_16_m(:,:,:)  
  REAL(dp), PUBLIC, POINTER :: O3ONPUS_17_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3ONMUS_18_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3OTRUS_19_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3OSMUS_20_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: O3OSPUS_21_m(:,:,:)  

 ! chem2 stream
  REAL(dp), PUBLIC, POINTER :: totoz_m(:,:)
  REAL(dp), PUBLIC, POINTER :: sadsts_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: sadpsc1_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: sadpsc2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: jo2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: jo3d_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: jno2_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: jcl2o2_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: lfp_m(:,:) 
  REAL(dp), PUBLIC, POINTER :: lnox_m(:,:) 
  REAL(dp), PUBLIC, POINTER :: elnox_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: cu_uvelo_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: xupdr_m(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: rxn_1_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_2_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_3_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_4_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_5_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_6_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_7_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_8_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_9_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_10_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_11_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_12_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_13_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_14_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_15_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_16_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_17_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_18_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_19_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_20_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_21_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_22_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_23_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_24_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_25_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_26_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_27_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_28_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_29_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_30_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_31_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_32_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_33_m(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_34_m(:,:,:)

!!!!
  REAL(dp), PUBLIC, POINTER :: totoz_d(:,:)
  REAL(dp), PUBLIC, POINTER :: sadsts_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: sadpsc1_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: sadpsc2_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: jo2(:,:,:)
  REAL(dp), PUBLIC, POINTER :: jo3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: jno2(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: jcl2o2(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: lfp_d(:,:)
  REAL(dp), PUBLIC, POINTER :: lnox_d(:,:)
  REAL(dp), PUBLIC, POINTER :: elnox(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: cu_uvelo(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: xupdr(:,:,:) 
  REAL(dp), PUBLIC, POINTER :: rxn_1(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_2(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_3(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_4(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_5(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_6(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_7(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_8(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_9(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_10(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_11(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_12(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_13(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_14(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_15(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_16(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_17(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_18(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_19(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_20(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_21(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_22(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_23(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_24(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_25(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_26(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_27(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_28(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_29(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_30(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_31(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_32(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_33(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rxn_34(:,:,:)

!!! 
  REAL(dp), PUBLIC, POINTER :: famfixcl_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: famfixbr_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: famfixn_d(:,:,:)

  REAL(dp), PUBLIC, POINTER :: sedpsc1_d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: sedpsc2_d(:,:,:)
  
  REAL(dp), PUBLIC, POINTER :: t_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: o3_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: no_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: no2_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hno3_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: n2o5_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: n2o_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ch4_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: co_h(:,:,:)
  REAL(dp), PUBLIC, POINTER :: h2o_h(:,:,:)


CONTAINS

!-------------------------------------------------------------------------------
  SUBROUTINE construct_stream_nudg_2h

    ! Allocates output streams for 2 hourly output for CAO Moscow; 

    ! *construct_stream_nudg_2h* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (nudg_2h ,'nudg_2h', filetype=NETCDF, &  !output as netcdf-file
                   interval=io_time_event(2,TIME_INC_HOURS,TRIG_LAST,0))      


  CALL default_stream_setting (nudg_2h, lrerun    = .FALSE. , &
                                       leveltype = HYBRID , &
                                       table     = 495,     &
                                       laccu     = .FALSE.,&
                                       contnorest = .FALSE.)

  CALL add_stream_element (nudg_2h, 'te_h',     t_h,      lpost=.TRUE.,  &
                           longname='Temperature',                &
                           units='K',       code=9 )
  CALL add_stream_element (nudg_2h, 'O3_h',     o3_h,     lpost=.TRUE.,  &
                           longname='O3 volume mixing ratio',     &
                           units='mol/mol', code=11)
  CALL add_stream_element (nudg_2h, 'NO_h',     no_h,     lpost=.TRUE.,  &
                           longname='NO volume mixing ratio',     &
                           units='mol/mol', code=13)
  CALL add_stream_element (nudg_2h, 'NO2_h',    no2_h,    lpost=.TRUE.,  &
                           longname='NO2 volume mixing ratio',    &
                           units='mol/mol', code=14)
  CALL add_stream_element (nudg_2h, 'HNO3_h',   hno3_h,   lpost=.TRUE.,  &
                           longname='HNO3 volume mixing ratio',   &
                           units='mol/mol', code=15)
  CALL add_stream_element (nudg_2h, 'N2O5_h',   n2o5_h,   lpost=.TRUE.,  &
                           longname='N2O5 volume mixing ratio',   &
                           units='mol/mol', code=17)
  CALL add_stream_element (nudg_2h, 'N2O_h',    n2o_h,    lpost=.TRUE.,  &
                           longname='N2O volume mixing ratio',    &
                           units='mol/mol', code=29)
  CALL add_stream_element (nudg_2h, 'CH4_h',    ch4_h,    lpost=.TRUE.,  &
                           longname='CH4 volume mixing ratio',    &
                           units='mol/mol', code=30)
  CALL add_stream_element (nudg_2h, 'CO_h',     co_h,     lpost=.TRUE.,  &
                           longname='CO volume mixing ratio',     &
                           units='mol/mol', code=31)
  CALL add_stream_element (nudg_2h, 'H2O_h',    h2o_h,    lpost=.TRUE.,  &
                           longname='H2O volume mixing ratio',    &
                           units='mol/mol', code=35)

  CALL default_stream_setting ( nudg_2h, lrerun = .FALSE.,contnorest = .FALSE.)
 
END SUBROUTINE construct_stream_nudg_2h




  SUBROUTINE construct_stream_chem_m

    ! Allocates output streams for monthly mean output of MEZON; 
    ! besides calculates time interval between two chemical time steps 
    ! (=radiation time steps).

    ! *construct_stream_chem_m* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (chem_m ,'chem1', filetype=NETCDF, &  !output as netcdf-file
                   interval=io_time_event(1,TIME_INC_DAYS,TRIG_LAST,0))      


  CALL default_stream_setting (chem_m, lrerun    = .FALSE. , &
                                       leveltype = HYBRID , &
                                       table     = 199,     &
                                       laccu     = .TRUE.)

  CALL add_stream_reference (chem_m, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
  CALL add_stream_reference (chem_m, 'lsp'     ,'sp'    ,lpost=.TRUE.)
  CALL add_stream_reference (chem_m, 'aps'     ,'g3b'    ,lpost=.TRUE.)


  CALL add_stream_element (chem_m, 'O3',     o3_m,     lpost=.TRUE.,  &
                           longname='O3 volume mixing ratio',     &
                           units='mol/mol', code=11)
  CALL add_stream_element (chem_m, 'O' ,     o_m,      lpost=.TRUE.,  &
                           longname='O volume mixing ratio',      &
                           units='mol/mol', code=12)
  CALL add_stream_element (chem_m, 'NO',     no_m,     lpost=.TRUE.,  &
                           longname='NO volume mixing ratio',     &
                           units='mol/mol', code=13)
  CALL add_stream_element (chem_m, 'NO2',    no2_m,    lpost=.TRUE.,  &
                           longname='NO2 volume mixing ratio',    &
                           units='mol/mol', code=14)
  CALL add_stream_element (chem_m, 'HNO3',   hno3_m,   lpost=.FALSE.,  &
                           longname='HNO3 volume mixing ratio',   &
                           units='mol/mol', code=15)
  CALL add_stream_element (chem_m, 'NO3',    no3_m,    lpost=.FALSE.,  &
                           longname='NO3 volume mixing ratio',    &
                           units='mol/mol', code=16)
  CALL add_stream_element (chem_m, 'N2O5',   n2o5_m,   lpost=.FALSE.,  &
                           longname='N2O5 volume mixing ratio',   &
                           units='mol/mol', code=17)
  CALL add_stream_element (chem_m, 'HNO4',   hno4_m,   lpost=.FALSE.,  &
                           longname='HNO4 volume mixing ratio',   &
                           units='mol/mol', code=18)
  CALL add_stream_element (chem_m, 'HO2',    ho2_m,    lpost=.TRUE.,  &
                           longname='HO2 volume mixing ratio',    &
                           units='mol/mol', code=19)
  CALL add_stream_element (chem_m, 'ClNO3',  clno3_m,  lpost=.FALSE.,  &
                           longname='ClNO3 volume mixing ratio',  &
                           units='mol/mol', code=20)
  CALL add_stream_element (chem_m, 'ClO',    clo_m,    lpost=.FALSE.,  &
                           longname='ClO volume mixing ratio',    &
                           units='mol/mol', code=21)
  CALL add_stream_element (chem_m, 'N',      n_m,      lpost=.FALSE.,  &
                           longname='N volume mixing ratio',      &
                           units='mol/mol', code=22)
  CALL add_stream_element (chem_m, 'OH',     oh_m,     lpost=.TRUE.,  &
                           longname='OH volume mixing ratio',     &
                           units='mol/mol', code=23)
  CALL add_stream_element (chem_m, 'H',      h_m,      lpost=.FALSE.,  &
                           longname='H volume mixing ratio',      &
                           units='mol/mol', code=24)
  CALL add_stream_element (chem_m, 'Cl',     cl_m,     lpost=.FALSE.,  &
                           longname='Cl volume mixing ratio',     &
                           units='mol/mol', code=25)
  CALL add_stream_element (chem_m, 'HOCl',   hocl_m,   lpost=.FALSE.,  &
                           longname='HOCl volume mixing ratio',   &
                           units='mol/mol', code=26)

  CALL add_stream_element (chem_m, 'N2O',    n2o_m,    lpost=.TRUE.,  &
                           longname='N2O volume mixing ratio',    &
                           units='mol/mol', code=29)
  CALL add_stream_element (chem_m, 'CH4',    ch4_m,    lpost=.FALSE.,  &
                           longname='CH4 volume mixing ratio',    &
                           units='mol/mol', code=30)
  CALL add_stream_element (chem_m, 'CO',     co_m,     lpost=.FALSE.,  &
                           longname='CO volume mixing ratio',     &
                           units='mol/mol', code=31)
  CALL add_stream_element (chem_m, 'HCl',    hcl_m,    lpost=.FALSE.,  &
                           longname='HCl volume mixing ratio',    &
                           units='mol/mol', code=32)
  CALL add_stream_element (chem_m, 'H2',     h2_m,     lpost=.FALSE.,  &
                           longname='H2 volume mixing ratio',     &
                           units='mol/mol', code=33)
  CALL add_stream_element (chem_m, 'H2O2',   h2o2_m,   lpost=.FALSE.,  &
                           longname='H2O2 volume mixing ratio',   &
                           units='mol/mol', code=34)

  CALL add_stream_element (chem_m, 'H2O',    h2o_m,    lpost=.TRUE.,  &
                           longname='H2O volume mixing ratio',    &
                           units='mol/mol', code=35)

  CALL add_stream_element (chem_m, 'Cl2',    cl2_m,    lpost=.FALSE.,  &
                           longname='Cl2 volume mixing ratio',    &
                           units='mol/mol', code=36)
  CALL add_stream_element (chem_m, 'Cl2O2',  cl2o2_m,  lpost=.FALSE.,  &
                           longname='Cl2O2 volume mixing ratio',  &
                           units='mol/mol', code=37)

  CALL add_stream_element (chem_m, 'CH3O2H', ch3o2h_m, lpost=.FALSE.,  &
                           longname='CH3O2H volume mixing ratio', &
                           units='mol/mol', code=40)
  CALL add_stream_element (chem_m, 'CH2O',   ch2o_m,   lpost=.FALSE.,  &
                           longname='cH2O volume mixing ratio',   &
                           units='mol/mol', code=41)
  CALL add_stream_element (chem_m, 'BrO',    bro_m,    lpost=.FALSE.,  &
                           longname='BrO volume mixing ratio',    &
                           units='mol/mol', code=42)
  CALL add_stream_element (chem_m, 'BrNO3',  brno3_m,  lpost=.FALSE.,  &
                           longname='BrNO3 volume mixing ratio',  &
                           units='mol/mol', code=43)
  CALL add_stream_element (chem_m, 'BrCl',   brcl_m,   lpost=.FALSE.,  &
                           longname='BrCl volume mixing ratio',   &
                           units='mol/mol', code=44)
  CALL add_stream_element (chem_m, 'HBr',    hbr_m,    lpost=.FALSE.,  &
                           longname='HBr volume mixing ratio',    &
                           units='mol/mol', code=45)
  CALL add_stream_element (chem_m, 'HOBr',   hobr_m,   lpost=.FALSE.,  &
                           longname='HOBr volume mixing ratio',   &
                           units='mol/mol', code=46)

  CALL add_stream_element (chem_m, 'Br',     br_m,     lpost=.FALSE.,  &
                           longname='Br volume mixing ratio',     &
                           units='mol/mol', code=52)
  CALL add_stream_element (chem_m, 'OD',     od_m,     lpost=.FALSE.,  &
                           longname='OD volume mixing ratio',     &
                           units='mol/mol', code=53)
  CALL add_stream_element (chem_m, 'CH3CO3',     ch3co3_m,     lpost=.FALSE.,  &
                           longname='Peroxyacetyl radical volume mixing ratio',     &
                           units='mol/mol', code=54)
  CALL add_stream_element (chem_m, 'PAN', pan_m, lpost=.FALSE.,  &
                           longname='Peroxyacetylnitrate volume mixing ratio',&
                           units='mol/mol', code=55)
  CALL add_stream_element (chem_m, 'C5H8', c5h8_m, lpost=.FALSE.,  &
                           longname='Isoprene volume mixing ratio',&
                           units='mol/mol', code=66)
  CALL add_stream_element (chem_m, 'CFC11',f11_m,lpost=.FALSE., &
                           longname='CFC11 volume mixing ratio',           &
                           units='mol/mol', code=70)
  CALL add_stream_element (chem_m, 'CFC12',f12_m,lpost=.FALSE., &
                           longname='CFC12 volume mixing ratio',           &
                           units='mol/mol', code=71)
  CALL add_stream_element (chem_m, 'CBRF3',cbrf3_m,lpost=.FALSE., &
                           longname='CBRF3 volume mixing ratio',           &
                           units='mol/mol', code=72)
  CALL add_stream_element (chem_m, 'CFC113',cfc113_m,lpost=.FALSE., &
                           longname='CFC113 volume mixing ratio',           &
                           units='mol/mol', code=73)
  CALL add_stream_element (chem_m, 'CFC114',cfc114_m,lpost=.FALSE., &
                           longname='CFC114 volume mixing ratio',           &
                           units='mol/mol', code=74)
  CALL add_stream_element (chem_m, 'CFC115',cfc115_m,lpost=.FALSE., &
                           longname='CFC115 volume mixing ratio',           &
                           units='mol/mol', code=75)
  CALL add_stream_element (chem_m, 'CCL4',ccl4_m,lpost=.FALSE., &
                           longname='CCL4 volume mixing ratio',           &
                           units='mol/mol', code=76)
  CALL add_stream_element (chem_m, 'CH3CCL3',ch3ccl3_m,lpost=.FALSE., &
                           longname='CH3CCL3 volume mixing ratio',           &
                           units='mol/mol', code=77)
  CALL add_stream_element (chem_m, 'HCFC22',hcfc22_m,lpost=.FALSE., &
                           longname='HCFC22 volume mixing ratio',           &
                           units='mol/mol', code=78)
  CALL add_stream_element (chem_m, 'HCFC141b',hcfc141b_m,lpost=.FALSE., &
                           longname='HCFC141b volume mixing ratio',           &
                           units='mol/mol', code=79)
  CALL add_stream_element (chem_m, 'HCFC142b',hcfc142b_m,lpost=.FALSE., &
                           longname='HCFC142b volume mixing ratio',           &
                           units='mol/mol', code=80)
  CALL add_stream_element (chem_m, 'H1211',h1211_m,lpost=.FALSE., &
                           longname='H1211 volume mixing ratio',           &
                           units='mol/mol', code=81)
  CALL add_stream_element (chem_m, 'CH3BR',ch3br_m,lpost=.FALSE., &
                           longname='CH3BR volume mixing ratio',           &
                           units='mol/mol', code=82)
  CALL add_stream_element (chem_m, 'CH3CL',ch3cl_m,lpost=.FALSE., &
                           longname='CH3CL volume mixing ratio',           &
                           units='mol/mol', code=83)
  CALL add_stream_element (chem_m, 'HCFC21',hcfc21_m,lpost=.FALSE., &
                           longname='HCFC21 volume mixing ratio',           &
                           units='mol/mol', code=84)
  CALL add_stream_element (chem_m, 'HCFC123',hcfc123_m,lpost=.FALSE., &
                           longname='HCFC123 volume mixing ratio',           &
                           units='mol/mol', code=85)
  CALL add_stream_element (chem_m, 'H2402',h2402_m,lpost=.FALSE., &
                           longname='H2402 volume mixing ratio',           &
                           units='mol/mol', code=86)
  CALL add_stream_element (chem_m, 'CHBR3',chbr3_m,lpost=.FALSE., &
                           longname='CHBR3 volume mixing ratio',           &
                           units='mol/mol', code=87)
  CALL add_stream_element (chem_m, 'CH2BR2',ch2br2_m,lpost=.FALSE., &
                           longname='CH2BR2 volume mixing ratio',           &
                           units='mol/mol', code=88)
  CALL add_stream_element (chem_m, 'LinearAge',linage_m,lpost=.FALSE., &
                           longname='Linearly increasing age tracer volume mixing ratio',           &
                           units='mol/mol', code=89)
  CALL add_stream_element (chem_m, 'IdealAge',idealage_m,lpost=.FALSE., &
                           longname='Ideal age tracer volume mixing ratio',           &
                           units='mol/mol', code=90)
  !CCMI synthetic tracers
  if (lsynth) then
    CALL add_stream_element (chem_m, 'NH_5',NH_5_m,lpost=.FALSE., &
                           longname='CCMI tracer 30-50?,N 5-day decay',      &
                           units='mol/mol', code=91)
    CALL add_stream_element (chem_m, 'NH_50',NH_50_m,lpost=.FALSE., &
                           longname='CCMI tracer 30-50?,N 50-day decay', &
                           units='mol/mol', code=92)
    CALL add_stream_element (chem_m, 'NH_50W',NH_50W_m,lpost=.FALSE., &
                           longname='CCMI tracer 3050?N, 50-day decay + wet deposition', &
                           units='mol/mol', code=93)
    CALL add_stream_element (chem_m, 'ST80_25',ST80_25_m,lpost=.FALSE., &
                           longname='CCMI tracer above 80hPa, 25-day troposphere decay',       &
                           units='mol/mol', code=94)
    CALL add_stream_element (chem_m, 'CO_25',CO_25_m,lpost=.FALSE., &
                           longname='CCMI tracer anthropogenic CO, 25-day decay',      &
                           units='mol/mol', code=95)
    CALL add_stream_element (chem_m, 'CO_50',CO_50_m,lpost=.FALSE., &
                           longname='CCMI tracer anthropogenic CO, 50-day decay',      &
                           units='mol/mol', code=96)
  endif

   CALL add_stream_element (chem_m, 'O3ONPBL_1',O3ONPBL_1_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONPBL_1',      &
                           units='mol/mol', code=101)
   CALL add_stream_element (chem_m, 'O3ONMBL_2',O3ONMBL_2_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONMBL_2',      &
                           units='mol/mol', code=102)
   CALL add_stream_element (chem_m, 'O3OTRBL_3',O3OTRBL_3_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OTRBL_3',      &
                           units='mol/mol', code=103)
   CALL add_stream_element (chem_m, 'O3OSMBL_4',O3OSMBL_4_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSMBL_4',      &
                           units='mol/mol', code=104)
   CALL add_stream_element (chem_m, 'O3OSPBL_5',O3OSPBL_5_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSPBL_5',      &
                           units='mol/mol', code=105)
   CALL add_stream_element (chem_m, 'O3ONPFT_6',O3ONPFT_6_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONPFT_6',      &
                           units='mol/mol', code=106)
   CALL add_stream_element (chem_m, 'O3ONMFT_7',O3ONMFT_7_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONMFT_7',      &
                           units='mol/mol', code=107)
   CALL add_stream_element (chem_m, 'O3OTRFT_8',O3OTRFT_8_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OTRFT_8',      &
                           units='mol/mol', code=108)
   CALL add_stream_element (chem_m, 'O3OSMFT_9',O3OSMFT_9_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSMFT_9',      &
                           units='mol/mol', code=109)
   CALL add_stream_element (chem_m, 'O3OSPFT_10',O3OSPFT_10_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSPFT_10',      &
                           units='mol/mol', code=110)
   CALL add_stream_element (chem_m, 'O3ONPLS_11',O3ONPLS_11_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONPLS_11',      &
                           units='mol/mol', code=111)
   CALL add_stream_element (chem_m, 'O3ONMLS_12',O3ONMLS_12_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONMLS_12',      &
                           units='mol/mol', code=112)
   CALL add_stream_element (chem_m, 'O3OTRLS_13',O3OTRLS_13_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OTRLS_13',      &
                           units='mol/mol', code=113)
   CALL add_stream_element (chem_m, 'O3OTRMS_14',O3OTRMS_14_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OTRMS_14',      &
                           units='mol/mol', code=114)
   CALL add_stream_element (chem_m, 'O3OSMLS_15',O3OSMLS_15_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSMLS_15',      &
                           units='mol/mol', code=115)
   CALL add_stream_element (chem_m, 'O3OSPLS_16',O3OSPLS_16_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSPLS_16',      &
                           units='mol/mol', code=116)
   CALL add_stream_element (chem_m, 'O3ONPUS_17',O3ONPUS_17_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONPUS_17',      &
                           units='mol/mol', code=117)
   CALL add_stream_element (chem_m, 'O3ONMUS_18',O3ONMUS_18_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3ONMUS_18',      &
                           units='mol/mol', code=118)
   CALL add_stream_element (chem_m, 'O3OTRUS_19',O3OTRUS_19_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OTRUS_19',      &
                           units='mol/mol', code=119)
   CALL add_stream_element (chem_m, 'O3OSMUS_20',O3OSMUS_20_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSMUS_20',      &
                           units='mol/mol', code=120)
   CALL add_stream_element (chem_m, 'O3OSPUS_21',O3OSPUS_21_m,lpost=.FALSE., &
                           longname='O3orig tracer - O3OSPUS_21',      &
                           units='mol/mol', code=121)

END SUBROUTINE construct_stream_chem_m


  SUBROUTINE construct_stream_chem2

    ! Allocates output streams for 12 h output of MEZON

    ! *construct_stream_chem2* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (chem2 ,'chem2', filetype=NETCDF, &  !output as netcdf-file
                   interval=io_time_event(1,TIME_INC_DAYS,TRIG_LAST,0))      


  CALL default_stream_setting (chem2, lrerun    = .FALSE. , &
                                       leveltype = HYBRID , &
                                       table     = 199,     &
                                       laccu     = .TRUE.)

  CALL add_stream_reference (chem2, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
  CALL add_stream_reference (chem2, 'lsp'     ,'sp'    ,lpost=.TRUE.)
  CALL add_stream_reference (chem2, 'aps'     ,'g3b'    ,lpost=.TRUE.)

  CALL add_stream_element (chem2, 'TOTOZ',   totoz_m,  lpost=.TRUE.,  &
                           longname='Ozone column',               &
                           units='DU',      code=54)

  CALL add_stream_element (chem2, 'SADSTS', sadsts_m,   lpost=.FALSE., &
                           longname='STS surface area density', &
                           units='um2/cm3', code=100)
  CALL add_stream_element (chem2, 'SADPSC1', sadpsc1_m, lpost=.FALSE., &
                           longname='PSC1 surface area density', &
                           units='um2/cm3', code=101)
  CALL add_stream_element (chem2, 'SADPSC2', sadpsc2_m, lpost=.FALSE., &
                           longname='PSC2 surface area density', &
                           units='um2/cm3', code=102)

   CALL add_stream_element (chem2, 'JO2', jo2_m, lpost=.TRUE.,  &
       longname='J-values O2->O+O',&
       units='1/s', code=110) 
   CALL add_stream_element (chem2, 'JO3d', jo3d_m, lpost=.TRUE.,  &
       longname='J-values  O3->O2+O(1D)',&
       units='1/s', code=111) 
   CALL add_stream_element (chem2, 'JNO2', jno2_m, lpost=.FALSE.,  &
       longname='J-values NO2->NO+O',&
       units='1/s', code=112) 
   CALL add_stream_element (chem2, 'JCl2O2', jcl2o2_m, lpost=.FALSE.,  &
       longname='J-values Cl2O2->2Cl',&
       units='1/s', code=113) 

   CALL add_stream_element (chem2, 'lfp',lfp_m,lpost=.FALSE.,  &
        longname='Lightning flash frequency by Price',&
        units='1/s', code=114)
   CALL add_stream_element (chem2, 'lnox',lnox_m,lpost=.FALSE.,  &
        longname='NOx from lightning',&
        units='g(NOx)/m**2/s', code=115)
   CALL add_stream_element (chem2, 'lnox_emiss',elnox_m,lpost=.FALSE.,  &
        longname='NOx from lightning',&
        units='kg(NOx)/m**2/s', code=116)
   CALL add_stream_element (chem2, 'cu_uvelo',cu_uvelo_m,lpost=.FALSE.,  &
        longname='convective updraft velocity',&
        units='m/s', code=117)
   CALL add_stream_element (chem2, 'xupdr',xupdr_m,lpost=.FALSE.,  &
        longname='convective mass flux',&
        units='kg/m2/s', code=118)

  CALL add_stream_element (chem2, 'RXN_1',rxn_1_m,lpost=.FALSE.,       &
                           longname='Ozone-depleting NOx cycles',   &
                           units='mol/cm3/s', code=119)
  CALL add_stream_element (chem2, 'RXN_2',rxn_2_m,lpost=.FALSE.,  &
                           longname='Ozone-depleting HOx cycles',   &
                           units='mol/cm3/s', code=120)
  CALL add_stream_element (chem2, 'RXN_3',rxn_3_m,lpost=.FALSE.,  &
                           longname='Ozone-depleting Cl cycles',    &
                           units='mol/cm3/s', code=121)
  CALL add_stream_element (chem2, 'RXN_4',rxn_4_m,lpost=.FALSE.,  &
                           longname='Ozone-depleting Br cycles',    &
                           units='mol/cm3/s', code=122)
  CALL add_stream_element (chem2, 'RXN_5',rxn_5_m,lpost=.FALSE.,  &
                           longname='Net chemical ozone production',&
                           units='mol/cm3/s', code=123)

  CALL add_stream_element (chem2, 'RXN_6',rxn_6_m,lpost=.FALSE.,  &
                           longname='Sum of the HO2+NO and RO2+NO reactions', &
                           units='mol/cm3/s', code=124)
  CALL add_stream_element (chem2, 'RXN_7',rxn_7_m,lpost=.FALSE.,  &
                           longname='Sum of OD+H2O, O3+HO2, O3+OH and C5H8+O3', &
                           units='mol/cm3/s', code=125)
  CALL add_stream_element (chem2, 'RXN_8',rxn_8_m,lpost=.FALSE.,  &
                           longname='NO2+hv->NO+O',                 &
                           units='mol/cm3/s', code=126)
  CALL add_stream_element (chem2, 'RXN_9',rxn_9_m,lpost=.FALSE.,  &
                           longname='OD+H2O->OH+OH',                 &
                           units='mol/cm3/s', code=127)
  CALL add_stream_element (chem2, 'RXN_10',rxn_10_m,lpost=.FALSE., &
                           longname='Total OH loss',                 &
                           units='mol/cm3/s', code=128)

  CALL add_stream_element (chem2, 'RXN_11',rxn_11_m,lpost=.FALSE.,&
                           longname='CO+OH+M->H+CO2+M', &
                           units='mol/cm3/s', code=156)
  CALL add_stream_element (chem2, 'RXN_12',rxn_12_m,lpost=.FALSE.,&
                           longname='CH4+OH->CH3+H2O',               &
                           units='mol/cm3/s', code=130)
  CALL add_stream_element (chem2, 'RXN_13',rxn_13_m,lpost=.FALSE.,&
                           longname='CH4+CL->CH3+HCL',               &
                           units='mol/cm3/s', code=131)
  CALL add_stream_element (chem2, 'RXN_14',rxn_14_m,lpost=.FALSE.,&
                           longname='H2O2 production',               &
                           units='mol/cm3/s', code=132)
  CALL add_stream_element (chem2, 'RXN_15',rxn_15_m,lpost=.FALSE., &
                           longname='Production of all hydrogen peroxides', &
                           units='mol/cm3/s', code=133)

  CALL add_stream_element (chem2, 'RXN_16',rxn_16_m,lpost=.FALSE., &
                           longname='HNO3 production',               &
                           units='mol/cm3/s', code=159)
  CALL add_stream_element (chem2, 'RXN_17',rxn_17_m,lpost=.FALSE., &
                           longname='Sum of RO2+NO reactions',       &
                           units='mol/cm3/s', code=135)
  CALL add_stream_element (chem2, 'RXN_18',rxn_18_m,lpost=.FALSE., &
                           longname='Sum of RO2+HO2 reactions',      &
                           units='mol/cm3/s', code=136)
  CALL add_stream_element (chem2, 'RXN_19',rxn_19_m,lpost=.FALSE., &
                           longname='Sum of RO2+RO2 reactions',      &
                           units='mol/cm3/s', code=137)
  CALL add_stream_element (chem2, 'RXN_20',rxn_20_m,lpost=.FALSE., &
                           longname='CH3CO3+NO2+M->PAN+M', &
                           units='mol/cm3/s', code=138)

  CALL add_stream_element (chem2, 'RXN_21',rxn_21_m,lpost=.FALSE., &
                           longname='O2 photolysis',                 &
                           units='mol/cm3/s', code=139)
  CALL add_stream_element (chem2, 'RXN_22',rxn_22_m,lpost=.FALSE., &
                           longname='Cl2O2 photolysis',              &
                           units='mol/cm3/s', code=140)
  CALL add_stream_element (chem2, 'RXN_23',rxn_23_m,lpost=.FALSE., &
                           longname='O3+hv->O(1D)',      &
                           units='mol/cm3/s', code=141)
  CALL add_stream_element (chem2, 'RXN_24',rxn_24_m,lpost=.FALSE., &
                           longname='Total O(1D) production',      &
                           units='mol/cm3/s', code=142)
  CALL add_stream_element (chem2, 'RXN_25',rxn_25_m,lpost=.FALSE., &
                           longname='Total ozone production', &
                           units='mol/cm3/s', code=143)

  CALL add_stream_element (chem2, 'RXN_26',rxn_26_m,lpost=.FALSE., &
                           longname='Total ozone loss',               &
                           units='mol/cm3/s', code=144)
  CALL add_stream_element (chem2, 'RXN_27',rxn_27_m,lpost=.FALSE., &
                           longname='O3+OH->HO2+O2',       &
                           units='mol/cm3/s', code=145)
  CALL add_stream_element (chem2, 'RXN_28',rxn_28_m,lpost=.FALSE., &
                           longname='O3+HO2->OH+O2',      &
                           units='mol/cm3/s', code=146)
  CALL add_stream_element (chem2, 'RXN_29',rxn_29_m,lpost=.FALSE., &
                           longname='HO2+NO->NO2+OH',      &
                           units='mol/cm3/s', code=147)
  CALL add_stream_element (chem2, 'RXN_30',rxn_30_m,lpost=.FALSE., &
                           longname='CH3O2+NO', &
                           units='mol/cm3/s', code=148)

  CALL add_stream_element (chem2, 'RXN_31',rxn_31_m,lpost=.FALSE., &
                           longname='O3+isoprene',               &
                           units='mol/cm3/s', code=149)
  CALL add_stream_element (chem2, 'RXN_32',rxn_32_m,lpost=.FALSE., &
                           longname='Total CH4 loss',       &
                           units='mol/cm3/s', code=150)
  CALL add_stream_element (chem2, 'RXN_33',rxn_33_m,lpost=.FALSE., &
                           longname='Total CO loss',      &
                           units='mol/cm3/s', code=151)
  CALL add_stream_element (chem2, 'RXN_34',rxn_34_m,lpost=.FALSE., &
                           longname='Total OH production',      &
                           units='mol/cm3/s', code=157)

END SUBROUTINE construct_stream_chem2

  SUBROUTINE construct_stream_chem2d


    ! *construct_stream_chem2* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (chem2d ,'chem2d')      
 
  CALL default_stream_setting (chem2, lrerun    = .TRUE. , &
                                       leveltype = HYBRID , &
                                       table     = 199,     &
                                       laccu     = .FALSE.)

  CALL add_stream_element (chem2d, 'TOTOZ',   totoz_d,  lpost=.FALSE.,  &
                           longname='Ozone column',               &
                           units='DU',      code=54)

  CALL add_stream_element (chem2d, 'FAMFIXCL',famfixcl_d,lpost=.FALSE., &
                           longname='Cly-family fixer',           &
                           units='-', code=95)
  CALL add_stream_element (chem2d, 'FAMFIXBR',famfixbr_d,lpost=.FALSE., &
                           longname='Bry-family fixer',           &
                           units='-', code=96)
  CALL add_stream_element (chem2d, 'FAMFIXN', famfixn_d, lpost=.FALSE., &
                           longname='NOy-family fixer',           &
                           units='-', code=97)

  CALL add_stream_element (chem2d, 'SEDPSC1', sedpsc1_d, lpost=.FALSE., &
                           longname='PSC1 sedimentation: HNO3 loss/gain', &
                           units='molec/cm3/s', code=98)
  CALL add_stream_element (chem2d, 'SEDPSC2', sedpsc2_d, lpost=.FALSE., &
                           longname='PSC2 sedimentation: H2O loss/gain', &
                           units='molec/cm3/s', code=99)
  
  CALL add_stream_element (chem2d, 'SADSTS', sadsts_d,   lpost=.FALSE., &
                           longname='STS surface area density', &
                           units='um2/cm3', code=100)
  CALL add_stream_element (chem2d, 'SADPSC1', sadpsc1_d, lpost=.FALSE., &
                           longname='PSC1 surface area density', &
                           units='um2/cm3', code=101)
  CALL add_stream_element (chem2d, 'SADPSC2', sadpsc2_d, lpost=.FALSE., &
                           longname='PSC2 surface area density', &
                           units='um2/cm3', code=102)
 
   CALL add_stream_element (chem2d, 'JO2', jo2, lpost=.FALSE.,  &
       longname='J-values O2->O+O',&
       units='1/s', code=110) 
   CALL add_stream_element (chem2d, 'JO3d', jo3d, lpost=.FALSE.,  &
       longname='J-values  O3->O2+O(1D)',&
       units='1/s', code=111) 
   CALL add_stream_element (chem2d, 'JNO2', jno2, lpost=.FALSE.,  &
       longname='J-values NO2->NO+O',&
       units='1/s', code=112) 
   CALL add_stream_element (chem2d, 'JCl2O2', jcl2o2, lpost=.FALSE.,  &
       longname='J-values Cl2O2->2Cl',&
       units='1/s', code=113) 

   CALL add_stream_element (chem2d, 'lfp',lfp_d,lpost=.FALSE.,  &
        longname='Lightning flash frequency by Price',&
        units='1/s', code=114)
   CALL add_stream_element (chem2d, 'lnox',lnox_d,lpost=.FALSE.,  &
        longname='NOx from lightning',&
        units='g(NOx)/m**2/s', code=115)
   CALL add_stream_element (chem2d, 'lnox_emiss',elnox,lpost=.FALSE.,  &
        longname='NOx from lightning',&
        units='kg(NOx)/m**2/s', code=116)
   CALL add_stream_element (chem2d, 'cu_uvelo',cu_uvelo,lpost=.FALSE.,  &
        longname='convective updraft velocity',&
        units='m/s', code=117)
   CALL add_stream_element (chem2d, 'xupdr',xupdr,lpost=.FALSE.,  &
        longname='convective mass flux',&
        units='kg/m2/s', code=118)

  CALL add_stream_element (chem2d, 'RXN_1',rxn_1,lpost=.FALSE.,       &
                           longname='Ozone-depleting NOx cycles',   &
                           units='mol/cm3/s', code=119)
  CALL add_stream_element (chem2d, 'RXN_2',rxn_2,lpost=.FALSE.,  &
                           longname='Ozone-depleting HOx cycles',   &
                           units='mol/cm3/s', code=120)
  CALL add_stream_element (chem2d, 'RXN_3',rxn_3,lpost=.FALSE.,  &
                           longname='Ozone-depleting Cl cycles',    &
                           units='mol/cm3/s', code=121)
  CALL add_stream_element (chem2d, 'RXN_4',rxn_4,lpost=.FALSE.,  &
                           longname='Ozone-depleting Br cycles',    &
                           units='mol/cm3/s', code=122)
  CALL add_stream_element (chem2d, 'RXN_5',rxn_5,lpost=.FALSE.,  &
                           longname='Net chemical ozone production',&
                           units='mol/cm3/s', code=123)

  CALL add_stream_element (chem2d, 'RXN_6',rxn_6,lpost=.FALSE.,  &
                           longname='Sum of the HO2+NO and RO2+NO reactions', &
                           units='mol/cm3/s', code=124)
  CALL add_stream_element (chem2d, 'RXN_7',rxn_7,lpost=.FALSE.,  &
                           longname='Sum of OD+H2O, O3+HO2, O3+OH and C5H8+O3', &
                           units='mol/cm3/s', code=125)
  CALL add_stream_element (chem2d, 'RXN_8',rxn_8,lpost=.FALSE.,  &
                           longname='NO2+hv->NO+O',                 &
                           units='mol/cm3/s', code=126)
  CALL add_stream_element (chem2d, 'RXN_9',rxn_9,lpost=.FALSE.,  &
                           longname='OD+H2O->OH+OH',                 &
                           units='mol/cm3/s', code=127)
  CALL add_stream_element (chem2d, 'RXN_10',rxn_10,lpost=.FALSE., &
                           longname='Total OH loss',                 &
                           units='mol/cm3/s', code=128)

  CALL add_stream_element (chem2d, 'RXN_11',rxn_11,lpost=.FALSE.,&
                           longname='CO+OH+M->H+CO2+M', &
                           units='mol/cm3/s', code=156)
  CALL add_stream_element (chem2d, 'RXN_12',rxn_12,lpost=.FALSE.,&
                           longname='CH4+OH->CH3+H2O',               &
                           units='mol/cm3/s', code=130)
  CALL add_stream_element (chem2d, 'RXN_13',rxn_13,lpost=.FALSE.,&
                           longname='CH4+CL->CH3+HCL',               &
                           units='mol/cm3/s', code=131)
  CALL add_stream_element (chem2d, 'RXN_14',rxn_14,lpost=.FALSE.,&
                           longname='H2O2 production',               &
                           units='mol/cm3/s', code=132)
  CALL add_stream_element (chem2d, 'RXN_15',rxn_15,lpost=.FALSE., &
                           longname='Production of all hydrogen peroxides', &
                           units='mol/cm3/s', code=133)

  CALL add_stream_element (chem2d, 'RXN_16',rxn_16,lpost=.FALSE., &
                           longname='HNO3 production',               &
                           units='mol/cm3/s', code=159)
  CALL add_stream_element (chem2d, 'RXN_17',rxn_17,lpost=.FALSE., &
                           longname='Sum of RO2+NO reactions',       &
                           units='mol/cm3/s', code=135)
  CALL add_stream_element (chem2d, 'RXN_18',rxn_18,lpost=.FALSE., &
                           longname='Sum of RO2+HO2 reactions',      &
                           units='mol/cm3/s', code=136)
  CALL add_stream_element (chem2d, 'RXN_19',rxn_19,lpost=.FALSE., &
                           longname='Sum of RO2+RO2 reactions',      &
                           units='mol/cm3/s', code=137)
  CALL add_stream_element (chem2d, 'RXN_20',rxn_20,lpost=.FALSE., &
                           longname='CH3CO3+NO2+M->PAN+M', &
                           units='mol/cm3/s', code=138)

  CALL add_stream_element (chem2d, 'RXN_21',rxn_21,lpost=.FALSE., &
                           longname='O2 photolysis',                 &
                           units='mol/cm3/s', code=139)
  CALL add_stream_element (chem2d, 'RXN_22',rxn_22,lpost=.FALSE., &
                           longname='Cl2O2 photolysis',              &
                           units='mol/cm3/s', code=140)
  CALL add_stream_element (chem2d, 'RXN_23',rxn_23,lpost=.FALSE., &
                           longname='O3+hv->O(1D)',      &
                           units='mol/cm3/s', code=141)
  CALL add_stream_element (chem2d, 'RXN_24',rxn_24,lpost=.FALSE., &
                           longname='Total O(1D) production',      &
                           units='mol/cm3/s', code=142)
  CALL add_stream_element (chem2d, 'RXN_25',rxn_25,lpost=.FALSE., &
                           longname='Total ozone production', &
                           units='mol/cm3/s', code=143)

  CALL add_stream_element (chem2d, 'RXN_26',rxn_26,lpost=.FALSE., &
                           longname='Total ozone loss',               &
                           units='mol/cm3/s', code=144)
  CALL add_stream_element (chem2d, 'RXN_27',rxn_27,lpost=.FALSE., &
                           longname='O3+OH->HO2+O2',       &
                           units='mol/cm3/s', code=145)
  CALL add_stream_element (chem2d, 'RXN_28',rxn_28,lpost=.FALSE., &
                           longname='O3+HO2->OH+O2',      &
                           units='mol/cm3/s', code=146)
  CALL add_stream_element (chem2d, 'RXN_29',rxn_29,lpost=.FALSE., &
                           longname='HO2+NO->NO2+OH',      &
                           units='mol/cm3/s', code=147)
  CALL add_stream_element (chem2d, 'RXN_30',rxn_30,lpost=.FALSE., &
                           longname='CH3O2+NO', &
                           units='mol/cm3/s', code=148)

  CALL add_stream_element (chem2d, 'RXN_31',rxn_31,lpost=.FALSE., &
                           longname='O3+isoprene',               &
                           units='mol/cm3/s', code=149)
  CALL add_stream_element (chem2d, 'RXN_32',rxn_32,lpost=.FALSE., &
                           longname='Total CH4 loss',       &
                           units='mol/cm3/s', code=150)
  CALL add_stream_element (chem2d, 'RXN_33',rxn_33,lpost=.FALSE., &
                           longname='Total CO loss',      &
                           units='mol/cm3/s', code=151)
  CALL add_stream_element (chem2d, 'RXN_34',rxn_34,lpost=.FALSE., &
                           longname='Total OH production',      &
                           units='mol/cm3/s', code=157)

END SUBROUTINE construct_stream_chem2d


SUBROUTINE destruct_stream_chem_m

  ! Deallocates memory.

  ! *destruct_stream_chem_m* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(chem_m)

END SUBROUTINE destruct_stream_chem_m


SUBROUTINE destruct_stream_chem2

  ! Deallocates memory.

  ! *destruct_stream_chem2* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(chem2)

END SUBROUTINE destruct_stream_chem2


SUBROUTINE destruct_stream_nudg_2h

  ! Deallocates memory.

  ! *destruct_stream_nudg_2h* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(nudg_2h)

END SUBROUTINE destruct_stream_nudg_2h


SUBROUTINE destruct_stream_chem2d

  ! Deallocates memory.

  ! *destruct_stream_chem2d* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(chem2d)

END SUBROUTINE destruct_stream_chem2d



SUBROUTINE init_stream_nudg_2h(kproma,krow,kbdim,ktrac,klev,pxtm1,pxtte,t)

  ! Initializes streams with zero.

  ! *init_stream_nudg_2h* is called from *mezon*,
  ! src/socol_mezon.f90.

  ! Scalar arguments
  INTEGER, INTENT(in) :: kproma, krow, kbdim, ktrac,klev

  ! Array arguments
  REAL(dp), INTENT(in) :: t(kbdim,klev)
  REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)

  t_h(1:kproma,:,krow)        = t(1:kproma,:)
  o3_h(1:kproma,:,krow)       = pxtm1(1:kproma,:,idt_o3)     + delta_time * pxtte(1:kproma,:,idt_o3)
  no_h(1:kproma,:,krow)       =  pxtm1(1:kproma,:,idt_no)     + delta_time * pxtte(1:kproma,:,idt_no)
  no2_h(1:kproma,:,krow)      =  pxtm1(1:kproma,:,idt_no2)     + delta_time * pxtte(1:kproma,:,idt_no2)
  hno3_h(1:kproma,:,krow)     =  pxtm1(1:kproma,:,idt_hno3)     + delta_time * pxtte(1:kproma,:,idt_hno3)
  n2o5_h(1:kproma,:,krow)     =  pxtm1(1:kproma,:,idt_n2o5)     + delta_time * pxtte(1:kproma,:,idt_n2o5)
  n2o_h(1:kproma,:,krow)      =  pxtm1(1:kproma,:,idt_n2o)     + delta_time * pxtte(1:kproma,:,idt_n2o)
  ch4_h(1:kproma,:,krow)      = pxtm1(1:kproma,:,idt_ch4)     + delta_time * pxtte(1:kproma,:,idt_ch4)
  co_h(1:kproma,:,krow)       =  pxtm1(1:kproma,:,idt_co)     + delta_time * pxtte(1:kproma,:,idt_co)
  h2o_h(1:kproma,:,krow)      =  pxtm1(1:kproma,:,idt_h2o)     + delta_time * pxtte(1:kproma,:,idt_h2o)

END SUBROUTINE init_stream_nudg_2h



SUBROUTINE init_stream_chem_m

  ! Initializes streams with zero.

  ! *init_stream_chem_m* is called from *call_init_tracers*,
  ! src/call_submodels.f90.

  o3_m(:,:,:)       = 0._dp
  o_m(:,:,:)        = 0._dp
  no_m(:,:,:)       = 0._dp
  no2_m(:,:,:)      = 0._dp
  hno3_m(:,:,:)     = 0._dp

  no3_m(:,:,:)      = 0._dp
  n2o5_m(:,:,:)     = 0._dp
  hno4_m(:,:,:)     = 0._dp
  ho2_m(:,:,:)      = 0._dp
  clno3_m(:,:,:)    = 0._dp

  clo_m(:,:,:)      = 0._dp
  n_m(:,:,:)        = 0._dp
  oh_m(:,:,:)       = 0._dp
  h_m(:,:,:)        = 0._dp
  cl_m(:,:,:)       = 0._dp

  hocl_m(:,:,:)     = 0._dp
  n2o_m(:,:,:)      = 0._dp
  ch4_m(:,:,:)      = 0._dp

  co_m(:,:,:)       = 0._dp
  hcl_m(:,:,:)      = 0._dp
  h2_m(:,:,:)       = 0._dp
  h2o2_m(:,:,:)     = 0._dp
  h2o_m(:,:,:)      = 0._dp

  cl2_m(:,:,:)      = 0._dp
  cl2o2_m(:,:,:)    = 0._dp
  ch3o2h_m(:,:,:)   = 0._dp

  ch2o_m(:,:,:)     = 0._dp
  bro_m(:,:,:)      = 0._dp
  brno3_m(:,:,:)    = 0._dp
  brcl_m(:,:,:)     = 0._dp
  hbr_m(:,:,:)      = 0._dp

  hobr_m(:,:,:)     = 0._dp
  ch3_m(:,:,:)      = 0._dp
  ch3o2_m(:,:,:)    = 0._dp
  ch3o_m(:,:,:)     = 0._dp

  hco_m(:,:,:)      = 0._dp
  br_m(:,:,:)       = 0._dp
  od_m(:,:,:)       = 0._dp

  ch3co3_m(:,:,:)  = 0._dp
  pan_m(:,:,:)     = 0._dp
  c5h8_m(:,:,:)    = 0._dp

  f11_m(:,:,:)    = 0._dp
  f12_m(:,:,:)    = 0._dp
  cbrf3_m(:,:,:)    = 0._dp
  cfc113_m(:,:,:)    = 0._dp
  cfc114_m(:,:,:)    = 0._dp
  cfc115_m(:,:,:)    = 0._dp
  ccl4_m(:,:,:)    = 0._dp
  ch3ccl3_m(:,:,:)    = 0._dp
  hcfc22_m(:,:,:)    = 0._dp
  hcfc141b_m(:,:,:)    = 0._dp
  hcfc142b_m(:,:,:)    = 0._dp
  h1211_m(:,:,:)    = 0._dp
  ch3br_m(:,:,:)    = 0._dp
  ch3cl_m(:,:,:)    = 0._dp
  hcfc21_m(:,:,:)    = 0._dp
  hcfc123_m(:,:,:)    = 0._dp
  h2402_m(:,:,:)    = 0._dp
  chbr3_m(:,:,:)    = 0._dp
  ch2br2_m(:,:,:)    = 0._dp
  linage_m(:,:,:)    = 0._dp
  idealage_m(:,:,:)    = 0._dp
 
  !CCMI synthetic tracers
  NH_5_m(:,:,:) = 0._dp
  NH_50_m(:,:,:) = 0._dp
  NH_50W_m(:,:,:) = 0._dp
  ST80_25_m(:,:,:) = 0._dp
  CO_25_m(:,:,:) = 0._dp
  CO_50_m(:,:,:) = 0._dp

  O3ONPBL_1_m(:,:,:) = 0._dp
  O3ONMBL_2_m(:,:,:) = 0._dp
  O3OTRBL_3_m(:,:,:) = 0._dp
  O3OSMBL_4_m(:,:,:) = 0._dp
  O3OSPBL_5_m(:,:,:) = 0._dp
  O3ONPFT_6_m(:,:,:) = 0._dp
  O3ONMFT_7_m(:,:,:) = 0._dp 
  O3OTRFT_8_m(:,:,:) = 0._dp
  O3OSMFT_9_m(:,:,:) = 0._dp 
  O3OSPFT_10_m(:,:,:) = 0._dp 
  O3ONPLS_11_m(:,:,:) = 0._dp
  O3ONMLS_12_m(:,:,:) = 0._dp 
  O3OTRLS_13_m(:,:,:)  = 0._dp
  O3OTRMS_14_m(:,:,:) = 0._dp
  O3OSMLS_15_m(:,:,:) = 0._dp
  O3OSPLS_16_m(:,:,:) = 0._dp
  O3ONPUS_17_m(:,:,:) = 0._dp
  O3ONMUS_18_m(:,:,:) = 0._dp
  O3OTRUS_19_m(:,:,:) = 0._dp
  O3OSMUS_20_m(:,:,:) = 0._dp
  O3OSPUS_21_m(:,:,:) = 0._dp

END SUBROUTINE init_stream_chem_m

SUBROUTINE init_stream_chem2

  ! Initializes streams with zero.
  
  ! *init_stream_chem2* is called from *call_init_tracers*,
  ! src/call_submodels.f90.
  
  totoz_m(:,:) = 0._dp
  sadsts_m(:,:,:) = 0._dp
  sadpsc1_m(:,:,:) = 0._dp
  sadpsc2_m(:,:,:) = 0._dp
  jo2_m(:,:,:) = 0._dp
  jo3d_m(:,:,:) = 0._dp
  jno2_m(:,:,:)  = 0._dp
  jcl2o2_m(:,:,:)  = 0._dp
  lfp_m(:,:)  = 0._dp
  lnox_m(:,:)  = 0._dp
  elnox_m(:,:,:)  = 0._dp
  cu_uvelo_m(:,:,:)  = 0._dp
  xupdr_m(:,:,:)  = 0._dp
  rxn_1_m(:,:,:) = 0._dp
  rxn_2_m(:,:,:) = 0._dp
  rxn_3_m(:,:,:) = 0._dp
  rxn_4_m(:,:,:) = 0._dp
  rxn_5_m(:,:,:) = 0._dp
  rxn_6_m(:,:,:) = 0._dp
  rxn_7_m(:,:,:) = 0._dp
  rxn_8_m(:,:,:) = 0._dp
  rxn_9_m(:,:,:) = 0._dp
  rxn_10_m(:,:,:) = 0._dp
  rxn_11_m(:,:,:) = 0._dp
  rxn_12_m(:,:,:) = 0._dp
  rxn_13_m(:,:,:) = 0._dp
  rxn_14_m(:,:,:) = 0._dp
  rxn_15_m(:,:,:) = 0._dp
  rxn_16_m(:,:,:) = 0._dp
  rxn_17_m(:,:,:) = 0._dp
  rxn_18_m(:,:,:) = 0._dp
  rxn_19_m(:,:,:) = 0._dp
  rxn_20_m(:,:,:) = 0._dp
  rxn_21_m(:,:,:) = 0._dp
  rxn_22_m(:,:,:) = 0._dp
  rxn_23_m(:,:,:) = 0._dp
  rxn_24_m(:,:,:) = 0._dp
  rxn_25_m(:,:,:) = 0._dp
  rxn_26_m(:,:,:) = 0._dp
  rxn_27_m(:,:,:) = 0._dp
  rxn_28_m(:,:,:) = 0._dp
  rxn_29_m(:,:,:) = 0._dp
  rxn_30_m(:,:,:) = 0._dp
  rxn_31_m(:,:,:) = 0._dp
  rxn_32_m(:,:,:) = 0._dp
  rxn_33_m(:,:,:) = 0._dp
  rxn_34_m(:,:,:) = 0._dp

END SUBROUTINE init_stream_chem2

SUBROUTINE init_stream_chem2d

  ! Initializes streams with zero.

  ! *init_stream_chem2* is called from *call_init_tracers*,
  ! src/call_submodels.f90.

  IF (lstart) THEN     ! use restart fields otherwise
     totoz_d(:,:)      = 0._dp
  ENDIF

  lfp_d(:,:) = 0._dp
  lnox_d(:,:) = 0._dp
  elnox(:,:,:) = 0._dp
  xupdr(:,:,:) = 0._dp
  cu_uvelo(:,:,:) = 0._dp

  IF (lstart) THEN     ! use restart fields otherwise
     famfixcl_d(:,:,:) = 0._dp
     famfixbr_d(:,:,:) = 0._dp
     famfixn_d(:,:,:)  = 0._dp

     sedpsc1_d(:,:,:)  = 0._dp
     sedpsc2_d(:,:,:)  = 0._dp

     sadsts_d(:,:,:)   = 0._dp
     sadpsc1_d(:,:,:)  = 0._dp
     sadpsc2_d(:,:,:)  = 0._dp

     jo3d(:,:,:) = 0._dp
     jno2(:,:,:) = 0._dp
     jo2(:,:,:) = 0._dp
     jcl2o2(:,:,:) = 0._dp

     rxn_1(:,:,:) = 0._dp
     rxn_2(:,:,:) = 0._dp
     rxn_3(:,:,:) = 0._dp
     rxn_4(:,:,:) = 0._dp
     rxn_5(:,:,:) = 0._dp

     rxn_6(:,:,:) = 0._dp
     rxn_7(:,:,:) = 0._dp
     rxn_8(:,:,:) = 0._dp
     rxn_9(:,:,:) = 0._dp
     rxn_10(:,:,:) = 0._dp

     rxn_11(:,:,:) = 0._dp
     rxn_12(:,:,:) = 0._dp
     rxn_13(:,:,:) = 0._dp
     rxn_14(:,:,:) = 0._dp 
     rxn_15(:,:,:) = 0._dp

     rxn_16(:,:,:) = 0._dp
     rxn_17(:,:,:) = 0._dp
     rxn_18(:,:,:) = 0._dp
     rxn_19(:,:,:) = 0._dp
     rxn_20(:,:,:) = 0._dp

     rxn_21(:,:,:) = 0._dp
     rxn_22(:,:,:) = 0._dp
     rxn_23(:,:,:) = 0._dp
     rxn_24(:,:,:) = 0._dp
     rxn_25(:,:,:) = 0._dp

     rxn_26(:,:,:) = 0._dp
     rxn_27(:,:,:) = 0._dp
     rxn_28(:,:,:) = 0._dp
     rxn_29(:,:,:) = 0._dp
     rxn_30(:,:,:) = 0._dp

     rxn_31(:,:,:) = 0._dp
     rxn_32(:,:,:) = 0._dp
     rxn_33(:,:,:) = 0._dp
     rxn_34(:,:,:) = 0._dp
  ENDIF
 
END SUBROUTINE init_stream_chem2d


SUBROUTINE accumulate_stream_chem_m (kproma, kbdim, klev, ktrac, krow, &
     pxtm1, pxtte)

  ! This subroutine accumulates the current value of a variable at every time 
  ! step, such that finally monthly streams can be calculated.

  ! *accumulate_stream_chem_m* is called from *call_diagn*, 
  ! src/call_submodels

  ! Scalar arguments
  INTEGER, INTENT(in) :: kproma, kbdim, klev, ktrac, krow

  ! Array arguments
  REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)
  
  ! Executable statements:
  o3_m(1:kproma,:,krow)       = o3_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_o3)     + delta_time * pxtte(1:kproma,:,idt_o3))
  o_m(1:kproma,:,krow)        = o_m(1:kproma,:,krow)        + delta_time * &
       (pxtm1(1:kproma,:,idt_o)      + delta_time * pxtte(1:kproma,:,idt_o))
  no_m(1:kproma,:,krow)       = no_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_no)     + delta_time * pxtte(1:kproma,:,idt_no))
  no2_m(1:kproma,:,krow)      = no2_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_no2)    + delta_time * pxtte(1:kproma,:,idt_no2))
  hno3_m(1:kproma,:,krow)     = hno3_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_hno3)   + delta_time * pxtte(1:kproma,:,idt_hno3))

  no3_m(1:kproma,:,krow)      = no3_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_no3)    + delta_time * pxtte(1:kproma,:,idt_no3))
  n2o5_m(1:kproma,:,krow)     = n2o5_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_n2o5)   + delta_time * pxtte(1:kproma,:,idt_n2o5))
  hno4_m(1:kproma,:,krow)     = hno4_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_hno4)   + delta_time * pxtte(1:kproma,:,idt_hno4))
  ho2_m(1:kproma,:,krow)      = ho2_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_ho2)    + delta_time * pxtte(1:kproma,:,idt_ho2))
  clno3_m(1:kproma,:,krow)    = clno3_m(1:kproma,:,krow)    + delta_time * &
       (pxtm1(1:kproma,:,idt_clno3)  + delta_time * pxtte(1:kproma,:,idt_clno3))

  clo_m(1:kproma,:,krow)      = clo_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_clo)    + delta_time * pxtte(1:kproma,:,idt_clo))
  n_m(1:kproma,:,krow)        = n_m(1:kproma,:,krow)        + delta_time * &
       (pxtm1(1:kproma,:,idt_n)      + delta_time * pxtte(1:kproma,:,idt_n))
  oh_m(1:kproma,:,krow)       = oh_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_oh)     + delta_time * pxtte(1:kproma,:,idt_oh))
  h_m(1:kproma,:,krow)        = h_m(1:kproma,:,krow)        + delta_time * &
       (pxtm1(1:kproma,:,idt_h)      + delta_time * pxtte(1:kproma,:,idt_h))
  cl_m(1:kproma,:,krow)       = cl_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_cl)     + delta_time * pxtte(1:kproma,:,idt_cl))

  hocl_m(1:kproma,:,krow)     = hocl_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_hocl)   + delta_time * pxtte(1:kproma,:,idt_hocl))
  n2o_m(1:kproma,:,krow)      = n2o_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_n2o)    + delta_time * pxtte(1:kproma,:,idt_n2o))
  ch4_m(1:kproma,:,krow)      = ch4_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_ch4)    + delta_time * pxtte(1:kproma,:,idt_ch4))

  co_m(1:kproma,:,krow)       = co_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_co)     + delta_time * pxtte(1:kproma,:,idt_co))
  hcl_m(1:kproma,:,krow)      = hcl_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_hcl)    + delta_time * pxtte(1:kproma,:,idt_hcl))
  h2_m(1:kproma,:,krow)       = h2_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_h2)     + delta_time * pxtte(1:kproma,:,idt_h2))
  h2o2_m(1:kproma,:,krow)     = h2o2_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_h2o2)   + delta_time * pxtte(1:kproma,:,idt_h2o2))
  h2o_m(1:kproma,:,krow)      = h2o_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_h2o)    + delta_time * pxtte(1:kproma,:,idt_h2o))

  cl2_m(1:kproma,:,krow)      = cl2_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_cl2)    + delta_time * pxtte(1:kproma,:,idt_cl2))
  cl2o2_m(1:kproma,:,krow)    = cl2o2_m(1:kproma,:,krow)    + delta_time * &
       (pxtm1(1:kproma,:,idt_cl2o2)  + delta_time * pxtte(1:kproma,:,idt_cl2o2))
  ch3o2h_m(1:kproma,:,krow)   = ch3o2h_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3o2h) + delta_time * pxtte(1:kproma,:,idt_ch3o2h))

  ch2o_m(1:kproma,:,krow)     = ch2o_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_ch2o)   + delta_time * pxtte(1:kproma,:,idt_ch2o))
  bro_m(1:kproma,:,krow)      = bro_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_bro)    + delta_time * pxtte(1:kproma,:,idt_bro))
  brno3_m(1:kproma,:,krow)    = brno3_m(1:kproma,:,krow)    + delta_time * &
       (pxtm1(1:kproma,:,idt_brno3)  + delta_time * pxtte(1:kproma,:,idt_brno3))
  brcl_m(1:kproma,:,krow)     = brcl_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_brcl)   + delta_time * pxtte(1:kproma,:,idt_brcl))
  hbr_m(1:kproma,:,krow)      = hbr_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_hbr)    + delta_time * pxtte(1:kproma,:,idt_hbr))

  hobr_m(1:kproma,:,krow)     = hobr_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_hobr)   + delta_time * pxtte(1:kproma,:,idt_hobr))
  ch3_m(1:kproma,:,krow)      = ch3_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3)    + delta_time * pxtte(1:kproma,:,idt_ch3))
  ch3o2_m(1:kproma,:,krow)    = ch3o2_m(1:kproma,:,krow)    + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3o2)  + delta_time * pxtte(1:kproma,:,idt_ch3o2))
  ch3o_m(1:kproma,:,krow)     = ch3o_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3o)   + delta_time * pxtte(1:kproma,:,idt_ch3o))

  hco_m(1:kproma,:,krow)      = hco_m(1:kproma,:,krow)      + delta_time * &
       (pxtm1(1:kproma,:,idt_hco)    + delta_time * pxtte(1:kproma,:,idt_hco))
  br_m(1:kproma,:,krow)       = br_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_br)     + delta_time * pxtte(1:kproma,:,idt_br))
  od_m(1:kproma,:,krow)       = od_m(1:kproma,:,krow)       + delta_time * &
       (pxtm1(1:kproma,:,idt_od)     + delta_time * pxtte(1:kproma,:,idt_od))

  ch3co3_m(1:kproma,:,krow) = ch3co3_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3co3) + delta_time * pxtte(1:kproma,:,idt_ch3co3))
  pan_m(1:kproma,:,krow) = pan_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_pan) + delta_time * pxtte(1:kproma,:,idt_pan))
  c5h8_m(1:kproma,:,krow) = c5h8_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_c5h8) + delta_time * pxtte(1:kproma,:,idt_c5h8))

  f11_m(1:kproma,:,krow)   = f11_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_f11)   + delta_time * pxtte(1:kproma,:,idt_f11))
  f12_m(1:kproma,:,krow)   = f12_m(1:kproma,:,krow)     + delta_time * &
       (pxtm1(1:kproma,:,idt_f12)   + delta_time * pxtte(1:kproma,:,idt_f12))
  cbrf3_m(1:kproma,:,krow) = cbrf3_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_cbrf3)   + delta_time * pxtte(1:kproma,:,idt_cbrf3))
  cfc113_m(1:kproma,:,krow) = cfc113_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_cfc113)   + delta_time * pxtte(1:kproma,:,idt_cfc113))
  cfc114_m(1:kproma,:,krow) = cfc114_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_cfc114)   + delta_time * pxtte(1:kproma,:,idt_cfc114))
  cfc115_m(1:kproma,:,krow) = cfc115_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_cfc115)   + delta_time * pxtte(1:kproma,:,idt_cfc115))
  ccl4_m(1:kproma,:,krow) = ccl4_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_ccl4)   + delta_time * pxtte(1:kproma,:,idt_ccl4))
  ch3ccl3_m(1:kproma,:,krow) = ch3ccl3_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3ccl3)   + delta_time * pxtte(1:kproma,:,idt_ch3ccl3))
  hcfc22_m(1:kproma,:,krow) = hcfc22_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_hcfc22)   + delta_time * pxtte(1:kproma,:,idt_hcfc22))
  hcfc141b_m(1:kproma,:,krow) = hcfc141b_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_hcfc141b)   + delta_time * pxtte(1:kproma,:,idt_hcfc141b))    
  hcfc142b_m(1:kproma,:,krow) = hcfc142b_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_hcfc142b)   + delta_time * pxtte(1:kproma,:,idt_hcfc142b))    
  h1211_m(1:kproma,:,krow) = h1211_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_h1211)   + delta_time * pxtte(1:kproma,:,idt_h1211))
  ch3br_m(1:kproma,:,krow) = ch3br_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3br)   + delta_time * pxtte(1:kproma,:,idt_ch3br))
  ch3cl_m(1:kproma,:,krow) = ch3cl_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_ch3cl)   + delta_time * pxtte(1:kproma,:,idt_ch3cl))
  hcfc21_m(1:kproma,:,krow) = hcfc21_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_hcfc21)   + delta_time * pxtte(1:kproma,:,idt_hcfc21))
  hcfc123_m(1:kproma,:,krow) = hcfc123_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_hcfc123)   + delta_time * pxtte(1:kproma,:,idt_hcfc123))
  h2402_m(1:kproma,:,krow) = h2402_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_h2402)   + delta_time * pxtte(1:kproma,:,idt_h2402))
  chbr3_m(1:kproma,:,krow) = chbr3_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_chbr3)   + delta_time * pxtte(1:kproma,:,idt_chbr3))
  ch2br2_m(1:kproma,:,krow) = ch2br2_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_ch2br2)   + delta_time * pxtte(1:kproma,:,idt_ch2br2))

  linage_m(1:kproma,:,krow) = linage_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_linearage)   + delta_time * pxtte(1:kproma,:,idt_linearage))
  idealage_m(1:kproma,:,krow) = idealage_m(1:kproma,:,krow)   + delta_time * &
       (pxtm1(1:kproma,:,idt_idealage)   + delta_time * pxtte(1:kproma,:,idt_idealage))

  !CCMI synthetic tracers
  if (lsynth) then
    NH_5_m(1:kproma,:,krow) = NH_5_m(1:kproma,:,krow) + delta_time * &
         (pxtm1(1:kproma,:,idt_NH_5) + delta_time * pxtte(1:kproma,:,idt_NH_5))

    NH_50_m(1:kproma,:,krow) = NH_50_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_NH_50) + delta_time * pxtte(1:kproma,:,idt_NH_50))

    NH_50W_m(1:kproma,:,krow) = NH_50W_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_NH_50W) + delta_time * pxtte(1:kproma,:,idt_NH_50W))

    ST80_25_m(1:kproma,:,krow) = ST80_25_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_ST80_25) + delta_time * pxtte(1:kproma,:,idt_ST80_25))

    CO_25_m(1:kproma,:,krow) = CO_25_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_CO_25) + delta_time * pxtte(1:kproma,:,idt_CO_25))

    CO_50_m(1:kproma,:,krow) = CO_50_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,idt_CO_50) + delta_time * pxtte(1:kproma,:,idt_CO_50))
 end if

 O3ONPBL_1_m(1:kproma,:,krow) = O3ONPBL_1_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(1)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(1)))

 O3ONMBL_2_m(1:kproma,:,krow) = O3ONMBL_2_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(2)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(2)))

 O3OTRBL_3_m(1:kproma,:,krow) =  O3OTRBL_3_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(3)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(3)))

 O3OSMBL_4_m(1:kproma,:,krow) = O3OSMBL_4_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(4)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(4)))

 O3OSPBL_5_m(1:kproma,:,krow) = O3OSPBL_5_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(5)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(5)))

 O3ONPFT_6_m(1:kproma,:,krow) = O3ONPFT_6_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(6)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(6)))
 
 O3ONMFT_7_m(1:kproma,:,krow) = O3ONMFT_7_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(7)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(7)))

 O3OTRFT_8_m(1:kproma,:,krow) = O3OTRFT_8_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(8)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(8)))

 O3OSMFT_9_m(1:kproma,:,krow) = O3OSMFT_9_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(9)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(9)))

 O3OSPFT_10_m(1:kproma,:,krow) = O3OSPFT_10_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(10)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(10)))

 O3ONPLS_11_m(1:kproma,:,krow) = O3ONPLS_11_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(11)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(11)))
 
 O3ONMLS_12_m(1:kproma,:,krow) = O3ONMLS_12_m(1:kproma,:,krow)  + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(12)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(12)))

 O3OTRLS_13_m(1:kproma,:,krow) = O3OTRLS_13_m(1:kproma,:,krow)  + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(13)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(13)))

 O3OTRMS_14_m(1:kproma,:,krow) = O3OTRMS_14_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(14)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(14)))

 O3OSMLS_15_m(1:kproma,:,krow) = O3OSMLS_15_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(15)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(15)))  

 O3OSPLS_16_m(1:kproma,:,krow) = O3OSPLS_16_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(16)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(16)))

 O3ONPUS_17_m(1:kproma,:,krow) = O3ONPUS_17_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(17)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(17)))

 O3ONMUS_18_m(1:kproma,:,krow) = O3ONMUS_18_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(18)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(18)))

 O3OTRUS_19_m(1:kproma,:,krow) = O3OTRUS_19_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(19)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(19)))

 O3OSMUS_20_m(1:kproma,:,krow) = O3OSMUS_20_m(1:kproma,:,krow)  + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(20)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(20)))

 O3OSPUS_21_m(1:kproma,:,krow) = O3OSPUS_21_m(1:kproma,:,krow) + delta_time * &
       (pxtm1(1:kproma,:,itrac_o3orig(21)) + delta_time * pxtte(1:kproma,:,itrac_o3orig(21)))

END SUBROUTINE accumulate_stream_chem_m

SUBROUTINE accumulate_stream_chem2 (kproma, kbdim, klev, ktrac, krow)

  ! This subroutine accumulates the current value of a variable at every time 
  ! step, such that finally monthly streams can be calculated.

  ! *accumulate_stream_chem2* is called from *call_diagn*, 
  ! src/call_submodels

  ! Scalar arguments
  INTEGER, INTENT(in) :: kproma, kbdim, klev, ktrac, krow
  
  ! Executable statements:
  
  totoz_m(:,krow)  = totoz_m(:,krow)  + delta_time * totoz_d(:,krow)
  
  sadsts_m(1:kproma,:,krow)   = sadsts_m(1:kproma,:,krow)   + delta_time * &
       sadsts_d(1:kproma,:,krow)
  sadpsc1_m(1:kproma,:,krow)  = sadpsc1_m(1:kproma,:,krow)  + delta_time * &
       sadpsc1_d(1:kproma,:,krow)
  sadpsc2_m(1:kproma,:,krow)  = sadpsc2_m(1:kproma,:,krow)  + delta_time * &
       sadpsc2_d(1:kproma,:,krow)

  jo2_m(1:kproma,:,krow) = jo2_m(1:kproma,:,krow) + delta_time * &
       jo2(1:kproma,:,krow)
  jo3d_m(1:kproma,:,krow) = jo3d_m(1:kproma,:,krow) + delta_time * &
       jo3d(1:kproma,:,krow)
  jno2_m(1:kproma,:,krow) = jno2_m(1:kproma,:,krow) + delta_time * &
       jno2(1:kproma,:,krow)
  jcl2o2_m(1:kproma,:,krow) = jcl2o2_m(1:kproma,:,krow) + delta_time * &
       jcl2o2(1:kproma,:,krow)

  lfp_m(1:kproma,krow) = lfp_m(1:kproma,krow) + delta_time * &
       lfp_d(1:kproma,krow)
  lnox_m(1:kproma,krow) = lnox_m(1:kproma,krow) + delta_time * &
       lnox_d(1:kproma,krow)
  elnox_m(1:kproma,:,krow) = elnox_m(1:kproma,:,krow) + delta_time * &
       elnox(1:kproma,:,krow)
  cu_uvelo_m(1:kproma,:,krow) = cu_uvelo_m(1:kproma,:,krow) + delta_time * &
       cu_uvelo(1:kproma,:,krow)
  xupdr_m(1:kproma,:,krow) = xupdr_m(1:kproma,:,krow) + delta_time * &
       xupdr(1:kproma,:,krow)
  
  rxn_1_m(1:kproma,:,krow) = rxn_1_m(1:kproma,:,krow) + delta_time * &
       rxn_1(1:kproma,:,krow)
  rxn_2_m(1:kproma,:,krow) = rxn_2_m(1:kproma,:,krow) + delta_time * &
       rxn_2(1:kproma,:,krow)
  rxn_3_m(1:kproma,:,krow) = rxn_3_m(1:kproma,:,krow) + delta_time * &
       rxn_3(1:kproma,:,krow)
  rxn_4_m(1:kproma,:,krow) = rxn_4_m(1:kproma,:,krow) + delta_time * &
       rxn_4(1:kproma,:,krow)
  rxn_5_m(1:kproma,:,krow) = rxn_5_m(1:kproma,:,krow) + delta_time * &
       rxn_5(1:kproma,:,krow)

  rxn_6_m(1:kproma,:,krow) = rxn_6_m(1:kproma,:,krow) + delta_time * &
       rxn_6(1:kproma,:,krow)
  rxn_7_m(1:kproma,:,krow) = rxn_7_m(1:kproma,:,krow) + delta_time * &
       rxn_7(1:kproma,:,krow)
  rxn_8_m(1:kproma,:,krow) = rxn_8_m(1:kproma,:,krow) + delta_time * &
       rxn_8(1:kproma,:,krow)
  rxn_9_m(1:kproma,:,krow) = rxn_9_m(1:kproma,:,krow) + delta_time * &
       rxn_9(1:kproma,:,krow)
  rxn_10_m(1:kproma,:,krow) = rxn_10_m(1:kproma,:,krow) + delta_time * &
       rxn_10(1:kproma,:,krow)

  rxn_11_m(1:kproma,:,krow) = rxn_11_m(1:kproma,:,krow) + delta_time * &
       rxn_11(1:kproma,:,krow)
  rxn_12_m(1:kproma,:,krow) = rxn_12_m(1:kproma,:,krow) + delta_time * &
       rxn_12(1:kproma,:,krow)
  rxn_13_m(1:kproma,:,krow) = rxn_13_m(1:kproma,:,krow) + delta_time * &
       rxn_13(1:kproma,:,krow)
  rxn_14_m(1:kproma,:,krow) = rxn_14_m(1:kproma,:,krow) + delta_time * &
       rxn_14(1:kproma,:,krow)
  rxn_15_m(1:kproma,:,krow) = rxn_15_m(1:kproma,:,krow) + delta_time * &
       rxn_15(1:kproma,:,krow)

  rxn_16_m(1:kproma,:,krow) = rxn_16_m(1:kproma,:,krow) + delta_time * &
       rxn_16(1:kproma,:,krow)
  rxn_17_m(1:kproma,:,krow) = rxn_17_m(1:kproma,:,krow) + delta_time * &
       rxn_17(1:kproma,:,krow)
  rxn_18_m(1:kproma,:,krow) = rxn_18_m(1:kproma,:,krow) + delta_time * &
       rxn_18(1:kproma,:,krow)
  rxn_19_m(1:kproma,:,krow) = rxn_19_m(1:kproma,:,krow) + delta_time * &
       rxn_19(1:kproma,:,krow)
  rxn_20_m(1:kproma,:,krow) = rxn_20_m(1:kproma,:,krow) + delta_time * &
       rxn_20(1:kproma,:,krow)
  rxn_21_m(1:kproma,:,krow) = rxn_21_m(1:kproma,:,krow) + delta_time * &
       rxn_21(1:kproma,:,krow)
  rxn_22_m(1:kproma,:,krow) = rxn_22_m(1:kproma,:,krow) + delta_time * &
       rxn_22(1:kproma,:,krow)
  rxn_23_m(1:kproma,:,krow) = rxn_23_m(1:kproma,:,krow) + delta_time * &
       rxn_23(1:kproma,:,krow)
  rxn_24_m(1:kproma,:,krow) = rxn_24_m(1:kproma,:,krow) + delta_time * &
       rxn_24(1:kproma,:,krow)
  rxn_25_m(1:kproma,:,krow) = rxn_25_m(1:kproma,:,krow) + delta_time * &
       rxn_25(1:kproma,:,krow)

  rxn_26_m(1:kproma,:,krow) = rxn_26_m(1:kproma,:,krow) + delta_time * &
       rxn_26(1:kproma,:,krow)
  rxn_27_m(1:kproma,:,krow) = rxn_27_m(1:kproma,:,krow) + delta_time * &
       rxn_27(1:kproma,:,krow)
  rxn_28_m(1:kproma,:,krow) = rxn_28_m(1:kproma,:,krow) + delta_time * &
       rxn_28(1:kproma,:,krow)
  rxn_29_m(1:kproma,:,krow) = rxn_29_m(1:kproma,:,krow) + delta_time * &
       rxn_29(1:kproma,:,krow)
  rxn_30_m(1:kproma,:,krow) = rxn_30_m(1:kproma,:,krow) + delta_time * &
       rxn_30(1:kproma,:,krow)

  rxn_31_m(1:kproma,:,krow) = rxn_31_m(1:kproma,:,krow) + delta_time * &
       rxn_31(1:kproma,:,krow)
  rxn_32_m(1:kproma,:,krow) = rxn_32_m(1:kproma,:,krow) + delta_time * &
       rxn_32(1:kproma,:,krow)
  rxn_33_m(1:kproma,:,krow) = rxn_33_m(1:kproma,:,krow) + delta_time * &
       rxn_33(1:kproma,:,krow)
  rxn_34_m(1:kproma,:,krow) = rxn_34_m(1:kproma,:,krow) + delta_time * &
       rxn_34(1:kproma,:,krow)

END SUBROUTINE accumulate_stream_chem2
          
END MODULE mo_socol_streams

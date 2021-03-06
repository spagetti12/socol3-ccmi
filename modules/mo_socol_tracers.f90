MODULE mo_socol_tracers

  ! Contains subroutines to handle chemical species (tracers) for MEZON.

  ! Martin Schraner, ETH Zurich, September 2008

  USE mo_advection,         ONLY: iadvec   ! selected advection scheme 
  USE mo_constants,         ONLY: amw, amd
  USE mo_control,           ONLY: nlon, ngl, nlev
  USE mo_decomposition,     ONLY: dcl=>local_decomposition, &
                                   dcg=>global_decomposition
  USE mo_doctor,            ONLY: nout
  USE mo_exception,         ONLY: finish
  USE mo_kind,              ONLY: dp
  USE mo_memory_gl,         ONLY: xt, q    ! tracer field array (t)
  USE mo_memory_g1a,        ONLY: xtm1     ! tracer field array (t-1) 
  USE mo_mpi,               ONLY: p_bcast, p_io, p_parallel, p_parallel_io
  USE mo_semi_impl,         ONLY: eps
  USE mo_socol_constants,   ONLY: eps_mr
  USE mo_socol_namelist,    ONLY: lh2o_coupl, lpco2, vini_pco2
  USE mo_socol_readfile,    ONLY: socol_read_netcdf
  USE mo_time_control,      ONLY: lresume, time_step_len, lstart
  USE mo_tracdef,           ONLY: trlist, t_trlist, t_trinfo
  USE mo_tracer,            ONLY: ntrac, & ! number of tracers
                                  new_tracer, & ! request resources for a tracer
                                  GAS, & ! phase indicators
                                  OFF, ON, &
                                  trlist, RESTART, INITIAL, CONSTANT
  USE mo_transpose,         ONLY: scatter_gp

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  ! Number of chemical species:
  INTEGER, PUBLIC :: n_trac_chemspec

  ! Vector with indices in order of initialization file:
  INTEGER, ALLOCATABLE, PUBLIC :: trac_chemspec(:)

  ! Tracer indices (chemical species) [mass mixing ratio]:
  INTEGER, PUBLIC :: idt_o3,  idt_o,  idt_od,  idt_no,  idt_no2,  idt_hno3,   &
       idt_no3,    idt_n2o5,  idt_hno4,   idt_ho2,    idt_clno3,  idt_clo,    &
       idt_n,      idt_oh,    idt_h,      idt_cl,     idt_hocl,   idt_odscls, &
       idt_odscll, idt_n2o,   idt_ch4,    idt_co,     idt_hcl,    idt_h2,     &
       idt_h2o2,   idt_h2o,   idt_cl2,    idt_cl2o2,  idt_psc1,   idt_psc2,   &
       idt_ch3,    idt_ch3o2, idt_ch3o,   idt_ch2o,   idt_hco,    idt_ch3o2h, &
       idt_br,     idt_bro,   idt_hbr,    idt_hobr,   idt_brno3,  idt_brcl,   &
       idt_odsbr,                                                             &
       idt_cly_adv, idt_bry_adv, idt_noy_adv, idt_co2, &
       idt_ch3co3, idt_pan, idt_ch3co3h, idt_ch3cooh, idt_nald, idt_hcooh, & ! eth_as_tropchem
       idt_hacet, idt_mgly, idt_macr, idt_macro2, idt_macro2h, idt_mpan,   & ! eth_as_tropchem
       idt_c5h8, idt_iso2, idt_iso2h, idt_ison, &                               ! eth_as_tropchem
        idt_f11, idt_f12, idt_cbrf3, idt_cfc113, idt_cfc114, &
       idt_cfc115, idt_ccl4, idt_ch3ccl3, idt_hcfc22, idt_hcfc141b, &
       idt_hcfc142b, idt_h1211, idt_ch3br, idt_ch3cl, idt_hcfc21, &
       idt_hcfc123, idt_h2402, idt_chbr3, idt_ch2br2, &
       idt_linearage, idt_idealage
      
  ! Tracer tendencies due to chemistry:
  REAL(dp), ALLOCATABLE, PUBLIC :: xtte_chem(:,:,:,:)

  PUBLIC :: request_socol_tracers, init_socol_tracers, pos_socol_tracers, &
       destruct_socol_tracers

CONTAINS

  SUBROUTINE request_socol_tracers

    ! Defines chemical species and families (only used for advection scheme 
    ! corrections) as tracers.

    ! *request_socol_tracers* is called from *call_request_tracer*,
    ! src/call_submodels.f90.

    INTEGER :: trac_tran  ! switch for transport in ECHAM5
    INTEGER :: trac_vdiff ! switch for vertical diffusion in ECHAM5
    INTEGER :: trac_conv  ! switch for convection in ECHAM5
    
     ! Executable statements:

    IF (lpco2) THEN
       n_trac_chemspec = 81
    ELSE
       n_trac_chemspec = 80
    END IF

    ALLOCATE(trac_chemspec(n_trac_chemspec))

    trac_tran=iadvec
    trac_vdiff=ON
    trac_conv=ON


    trac_tran=iadvec
    trac_vdiff=ON
    trac_conv=ON
    
    ! 1. Chemical species (n_trac_chemspec):
    
    CALL new_tracer('O3',     'MEZON', idx=idt_o3,     nwrite=OFF,  &
         longname='O3 volume mixing ratio',     units='mol/mol', &
         code=11, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('O',      'MEZON', idx=idt_o,      nwrite=OFF, &
         longname='O volume mixing ratio',      units='mol/mol', &
         code=12, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('NO',     'MEZON', idx=idt_no,     nwrite=OFF, &
         longname='NO volume mixing ratio',     units='mol/mol', &
         code=13, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('NO2',    'MEZON', idx=idt_no2,    nwrite=OFF, &
         longname='NO2 volume mixing ratio',    units='mol/mol', &
         code=14, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('HNO3',   'MEZON', idx=idt_hno3,   nwrite=OFF, &
         longname='HNO3 volume mixing ratio',   units='mol/mol', &
         code=15, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    CALL new_tracer('NO3',    'MEZON', idx=idt_no3,    nwrite=OFF, &
         longname='NO3 volume mixing ratio',    units='mol/mol', &
         code=16, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('N2O5',   'MEZON', idx=idt_n2o5,   nwrite=OFF, &
         longname='N2O5 volume mixing ratio',   units='mol/mol', &
         code=17, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('HNO4',   'MEZON', idx=idt_hno4,   nwrite=OFF, &
         longname='HNO4 volume mixing ratio',   units='mol/mol', &
         code=18, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('HO2',    'MEZON', idx=idt_ho2,    nwrite=OFF, &
         longname='HO2 volume mixing ratio',    units='mol/mol', &
         code=19, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('CLNO3',  'MEZON', idx=idt_clno3,  nwrite=OFF, &
         longname='ClNO3 volume mixing ratio',  units='mol/mol', &
         code=20, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    CALL new_tracer('CLO',    'MEZON', idx=idt_clo,    nwrite=OFF, &
         longname='ClO volume mixing ratio',    units='mol/mol', &
         code=21, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('N',      'MEZON', idx=idt_n,      nwrite=OFF, &
         longname='N volume mixing ratio',      units='mol/mol', &
         code=22, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('OH',     'MEZON', idx=idt_oh,     nwrite=OFF, &
         longname='OH volume mixing ratio',     units='mol/mol', &
         code=23, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('H',      'MEZON', idx=idt_h,      nwrite=OFF, &
         longname='H volume mixing ratio',      units='mol/mol', &
         code=24, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('CL',     'MEZON', idx=idt_cl,     nwrite=OFF, &
         longname='Cl volume mixing ratio',     units='mol/mol', &
         code=25, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    CALL new_tracer('HOCL',   'MEZON', idx=idt_hocl,   nwrite=OFF, &
         longname='HOCl volume mixing ratio',   units='mol/mol', &
         code=26, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('ODSCLS', 'MEZON', idx=idt_odscls, nwrite=OFF, &
         longname='ODSCLS volume mixing ratio', units='mol/mol', &
         code=27, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('ODSCLL', 'MEZON', idx=idt_odscll, nwrite=OFF, &
         longname='ODSCLL volume mixing ratio', units='mol/mol', &
         code=28, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('N2O',    'MEZON', idx=idt_n2o,    nwrite=OFF, &
         longname='N2O volume mixing ratio',    units='mol/mol', &
         code=29, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('CH4',    'MEZON', idx=idt_ch4,    nwrite=OFF, &
         longname='CH4 volume mixing ratio',    units='mol/mol', &
         code=30, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    CALL new_tracer('CO',     'MEZON', idx=idt_co,     nwrite=OFF, &
         longname='Co volume mixing ratio',     units='mol/mol', &
         code=31, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('HCL',    'MEZON', idx=idt_hcl,    nwrite=OFF, &
         longname='HCl volume mixing ratio',    units='mol/mol', &
         code=32, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('H2',     'MEZON', idx=idt_h2,     nwrite=OFF, &
         longname='H2 volume mixing ratio',     units='mol/mol', &
         code=33, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('H2O2',   'MEZON', idx=idt_h2o2,   nwrite=OFF, &
         longname='H2O2 volume mixing ratio',   units='mol/mol', &
         code=34, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    IF (lh2o_coupl) THEN    ! H2O CTM = H2O GCM -> no advection necessary!
       CALL new_tracer('H2O',    'MEZON', idx=idt_h2o,    nwrite=OFF, &
            longname='H2O volume mixing ratio',    units='mol/mol', &
            code=35, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
            ntran=OFF, nvdiff=OFF, nconv=OFF, nphase=GAS)
    ELSE
       CALL new_tracer('H2O',    'MEZON', idx=idt_h2o,    nwrite=OFF, &
            longname='H2O volume mixing ratio',    units='mol/mol', &
            code=35, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
            ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    ENDIF
    
    CALL new_tracer('CL2',    'MEZON', idx=idt_cl2,    nwrite=OFF, &
         longname='Cl2 volume mixing ratio',    units='mol/mol', &
         code=36, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('CL2O2',  'MEZON', idx=idt_cl2o2,  nwrite=OFF, &
         longname='Cl2O2 volume mixing ratio',  units='mol/mol', &
         code=37, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('PSC1',   'MEZON', idx=idt_psc1,   nwrite=OFF, &
         longname='PSC1 volume mixing ratio',   units='mol/mol', &
         code=38, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=0, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    !PSC1 is included in H2O for transport (and output)  (see chem/em_chemini)
    CALL new_tracer('PSC2',   'MEZON', idx=idt_psc2,   nwrite=OFF, &
         longname='PSC2 volume mixing ratio',   units='mol/mol', &
         code=39, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=0, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    !PSC2 is included in H2O for transport (and output)  (see chem/em_chemini)
    CALL new_tracer('CH3O2H', 'MEZON', idx=idt_ch3o2h, nwrite=OFF, &
         longname='CH3O2H volume mixing ratio', units='mol/mol', &
         code=40, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    CALL new_tracer('CH2O',   'MEZON', idx=idt_ch2o,   nwrite=OFF, &
         longname='CH2O volume mixing ratio',   units='mol/mol', &
         code=41, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('BRO',    'MEZON', idx=idt_bro,    nwrite=OFF, &
         longname='BrO volume mixing ratio',    units='mol/mol', &
         code=42, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('BRNO3',  'MEZON', idx=idt_brno3,  nwrite=OFF, &
         longname='BrNO3 volume mixing ratio',  units='mol/mol', &
         code=43, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('BRCL',   'MEZON', idx=idt_brcl,   nwrite=OFF, &
         longname='BrCl volume mixing ratio',   units='mol/mol', &
         code=44, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('HBR',    'MEZON', idx=idt_hbr,    nwrite=OFF, &
         longname='HBr volume mixing ratio',    units='mol/mol', &
         code=45, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    CALL new_tracer('HOBR',   'MEZON', idx=idt_hobr,   nwrite=OFF, &
         longname='HOBr volume mixing ratio',   units='mol/mol', &
         code=46, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('ODSBR',  'MEZON', idx=idt_odsbr,  nwrite=OFF, &
         longname='ODSBR volume mixing ratio',  units='mol/mol', &
         code=47, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('CH3',    'MEZON', idx=idt_ch3,    nwrite=OFF, &
         longname='CH3 volume mixing ratio',    units='mol/mol', &
         code=48, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('CH3O2',  'MEZON', idx=idt_ch3o2,  nwrite=OFF, &
         longname='CH3O2 volume mixing ratio',  units='mol/mol', &
         code=49, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('CH3O',   'MEZON', idx=idt_ch3o,   nwrite=OFF, &
         longname='CH3O volume mixing ratio',   units='mol/mol', &
         code=50, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    CALL new_tracer('HCO',    'MEZON', idx=idt_hco,    nwrite=OFF, &
         longname='HCO volume mixing ratio',    units='mol/mol', &
         code=51, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('BR',     'MEZON', idx=idt_br,     nwrite=OFF,  &
         longname='Br volume mixing ratio',     units='mol/mol', &
         code=52, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('OD',     'MEZON', idx=idt_od,     nwrite=OFF, &
         longname='OD volume mixing ratio',     units='mol/mol', &
         code=53, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    ! eth_as_tropchem+
    ! Tracer for Isoprene scheme
    
    ! C1 and C2 compounds
    CALL new_tracer('CH3CO3',     'MEZON', idx=idt_ch3co3,     nwrite=OFF, &
         longname='Peroxyacetyl radical, volume mixing ratio',     units='mol/mol', &
         code=54, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('PAN',     'MEZON', idx=idt_pan,     nwrite=OFF, &
         longname='Peroxyacetylnitrate, volume mixing ratio',     units='mol/mol', &
         code=55, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('CH3CO3H',     'MEZON', idx=idt_ch3co3h,     nwrite=OFF, &
         longname='Peroxyacetic acid, volume mixing ratio',     units='mol/mol', &
         code=56, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS) 

    CALL new_tracer('CH3COOH',     'MEZON', idx=idt_ch3cooh,     nwrite=OFF, &
         longname='Acetic acid, volume mixing ratio',     units='mol/mol', &
         code=57, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('NALD',     'MEZON', idx=idt_nald,     nwrite=OFF, &
         longname='Nitrooxyacetaldehyde, volume mixing ratio',     units='mol/mol', &
         code=58, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('HCOOH',     'MEZON', idx=idt_hcooh,     nwrite=OFF, &
         longname='Formic acid, volume mixing ratio',     units='mol/mol', &
         code=59, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    ! C3 compounds
    CALL new_tracer('HACET',     'MEZON', idx=idt_hacet,     nwrite=OFF, &
         longname='Hydroxyacetone and other C3 ketones, volume mixing ratio',     units='mol/mol', &
         code=60, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('MGLY',     'MEZON', idx=idt_mgly,     nwrite=OFF, &
         longname='Methylglyoxal and other C3 aldehydes, volume mixing ratio',     units='mol/mol', &
         code=61, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    ! C4 compounds
    CALL new_tracer('MACR',     'MEZON', idx=idt_macr,     nwrite=OFF, &
         longname='Methacrolein and other C4 carbonyls, volume mixing ratio',     units='mol/mol', &
         code=62, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('MACRO2',     'MEZON', idx=idt_macro2,     nwrite=OFF, &
         longname='Peroxy radicals from MACR+OH, volume mixing ratio',     units='mol/mol', &
         code=63, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('MACRO2H',     'MEZON', idx=idt_macro2h,     nwrite=OFF, &
         longname='Hydroperoxides from MACRO2+HO2, volume mixing ratio',     units='mol/mol', &
         code=64, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('MPAN',     'MEZON', idx=idt_mpan,     nwrite=OFF, &
         longname='Peroxymethacryloylnitrate, volume mixing ratio',     units='mol/mol', &
         code=65, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    ! C5 compounds
    CALL new_tracer('C5H8',     'MEZON', idx=idt_c5h8,     nwrite=OFF, &
         longname='Isoprene, volume mixing ratio',     units='mol/mol', &
         code=66, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('ISO2',     'MEZON', idx=idt_iso2,     nwrite=OFF, &
         longname='Peroxy radicals from C5H8+OH, volume mixing ratio',     units='mol/mol', &
         code=67, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('ISO2H',     'MEZON', idx=idt_iso2h,     nwrite=OFF, &
         longname='Beta-Hydroxyhydroperoxides, volume mixing ratio',     units='mol/mol', &
         code=68, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('ISON',     'MEZON', idx=idt_ison,     nwrite=OFF, &
         longname='Beta-hydroxyalkylnitrates, volume mixing ratio',     units='mol/mol', &
         code=69, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    ! eth_as_tropchem-

    CALL new_tracer('CFC11',     'MEZON', idx=idt_f11,     nwrite=OFF, &
         longname='CFC11 volume mixing ratio',     units='mol/mol', &
         code=70, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('CFC12',     'MEZON', idx=idt_f12,     nwrite=OFF, &
         longname='CFC12 volume mixing ratio',     units='mol/mol', &
         code=71, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CBRF3',     'MEZON', idx=idt_cbrf3,     nwrite=OFF, &
         longname='CBRF3 volume mixing ratio',     units='mol/mol', &
         code=72, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CFC113',     'MEZON', idx=idt_cfc113,     nwrite=OFF, &
         longname='CFC113 volume mixing ratio',     units='mol/mol', &
         code=73, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CFC114',     'MEZON', idx=idt_cfc114,     nwrite=OFF, &
         longname='CFC114 volume mixing ratio',     units='mol/mol', &
         code=74, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CFC115',     'MEZON', idx=idt_cfc115,     nwrite=OFF, &
         longname='CFC115 volume mixing ratio',     units='mol/mol', &
         code=75, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CCL4',     'MEZON', idx=idt_ccl4,     nwrite=OFF, &
         longname='CCL4 volume mixing ratio',     units='mol/mol', &
         code=76, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CH3CCL3',     'MEZON', idx=idt_ch3ccl3,  nwrite=OFF, &
         longname='CH3CCL3 volume mixing ratio',     units='mol/mol', &
         code=77, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('HCFC22',     'MEZON', idx=idt_hcfc22,    nwrite=OFF, &
         longname='HCFC22 volume mixing ratio',     units='mol/mol', &
         code=78, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('HCFC141b',     'MEZON', idx=idt_hcfc141b, nwrite=OFF, &
         longname='HCFC141b volume mixing ratio',     units='mol/mol', &
         code=79, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('HCFC142b',     'MEZON', idx=idt_hcfc142b, nwrite=OFF, &
         longname='HCFC142b volume mixing ratio',     units='mol/mol', &
         code=80, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('H1211',     'MEZON', idx=idt_h1211,     nwrite=OFF, &
         longname='H1211 volume mixing ratio',     units='mol/mol', &
         code=81, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CH3BR',     'MEZON', idx=idt_ch3br,     nwrite=OFF, &
         longname='CH3BR volume mixing ratio',     units='mol/mol', &
         code=82, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CH3CL',     'MEZON', idx=idt_ch3cl,     nwrite=OFF, &
         longname='CH3CL volume mixing ratio',     units='mol/mol', &
         code=83, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('HCFC21',     'MEZON', idx=idt_hcfc21,   nwrite=OFF, &
         longname='HCFC21 volume mixing ratio',     units='mol/mol', &
         code=84, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('HCFC123',     'MEZON', idx=idt_hcfc123, nwrite=OFF, &
         longname='HCFC123 volume mixing ratio',     units='mol/mol', &
         code=85, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('H2402',     'MEZON', idx=idt_h2402,     nwrite=OFF, &
         longname='H2402 volume mixing ratio',     units='mol/mol', &
         code=86, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CHBR3',     'MEZON', idx=idt_chbr3,     nwrite=OFF, &
         longname='CHBR3 volume mixing ratio',     units='mol/mol', &
         code=87, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('CH2BR2',     'MEZON', idx=idt_ch2br2,   nwrite=OFF, &
         longname='CH2BR2 volume mixing ratio',     units='mol/mol', &
         code=88, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('LinearAge',     'MEZON', idx=idt_linearage,     nwrite=OFF, &
          longname='Linearly increasing age tracer volume mixing ratio',     units='mol/mol', &
          code=89, table=131, ninit=RESTART+CONSTANT, vini=0._dp, nrerun=ON, &
          ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

     CALL new_tracer('IdealAge',     'MEZON', idx=idt_idealage,     nwrite=OFF, &
          longname='Ideal age tracer volume mixing ratio',     units='mol/mol', &
          code=90, table=131, ninit=RESTART+CONSTANT, vini=1._dp, nrerun=ON, &
          ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    IF (lpco2) &
         CALL new_tracer('CO2',     'MEZON', idx=idt_co2,   nwrite=OFF, &
         longname='CO2 volume mixing ratio',     units='ppm', &
         code=91, table=131, ninit=RESTART+CONSTANT, vini=0.0_dp, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
    ! Vector with tracer indices in order of initialisation file:
    trac_chemspec( 1) = idt_o3
    trac_chemspec( 2) = idt_o
    trac_chemspec( 3) = idt_no
    trac_chemspec( 4) = idt_no2
    trac_chemspec( 5) = idt_no3
    
    trac_chemspec( 6) = idt_n2o5
    trac_chemspec( 7) = idt_hocl
    trac_chemspec( 8) = idt_oh
    trac_chemspec( 9) = idt_ho2
    trac_chemspec(10) = idt_hno3
    
    trac_chemspec(11) = idt_h2o2
    trac_chemspec(12) = idt_cl
    trac_chemspec(13) = idt_clo
    trac_chemspec(14) = idt_hcl
    trac_chemspec(15) = idt_clno3
    
    trac_chemspec(16) = idt_n
    trac_chemspec(17) = idt_hno4
    trac_chemspec(18) = idt_n2o
    trac_chemspec(19) = idt_h2o
    trac_chemspec(20) = idt_ch4
    
    trac_chemspec(21) = idt_odscls
    trac_chemspec(22) = idt_odscll
    trac_chemspec(23) = idt_co
    trac_chemspec(24) = idt_h
    trac_chemspec(25) = idt_cl2o2
    
    trac_chemspec(26) = idt_cl2
    trac_chemspec(27) = idt_h2
    trac_chemspec(28) = idt_psc1
    trac_chemspec(29) = idt_psc2
    trac_chemspec(30) = idt_ch3
    
    trac_chemspec(31) = idt_ch3o2
    trac_chemspec(32) = idt_ch3o
    trac_chemspec(33) = idt_ch2o
    trac_chemspec(34) = idt_hco
    trac_chemspec(35) = idt_ch3o2h
    
    trac_chemspec(36) = idt_br
    trac_chemspec(37) = idt_bro
    trac_chemspec(38) = idt_brno3
    trac_chemspec(39) = idt_hbr
    trac_chemspec(40) = idt_brcl
    
    trac_chemspec(41) = idt_hobr
    trac_chemspec(42) = idt_odsbr
    trac_chemspec(43) = idt_od

    ! eth_as_tropchem+
    trac_chemspec(44) = idt_c5h8
    trac_chemspec(45) = idt_iso2
    trac_chemspec(46) = idt_iso2h
    trac_chemspec(47) = idt_ison
    trac_chemspec(48) = idt_macr
    trac_chemspec(49) = idt_macro2
    trac_chemspec(50) = idt_macro2h
    trac_chemspec(51) = idt_mpan
    trac_chemspec(52) = idt_hacet
    trac_chemspec(53) = idt_mgly
    trac_chemspec(54) = idt_ch3co3
    trac_chemspec(55) = idt_pan
    trac_chemspec(56) = idt_ch3co3h
    trac_chemspec(57) = idt_ch3cooh
    trac_chemspec(58) = idt_nald
    trac_chemspec(59) = idt_hcooh
    ! eth_as_tropchem-

    trac_chemspec(60) = idt_f11
    trac_chemspec(61) = idt_f12

    trac_chemspec(62) = idt_cbrf3
    trac_chemspec(63) = idt_cfc113
    trac_chemspec(64) = idt_cfc114
    trac_chemspec(65) = idt_cfc115
    trac_chemspec(66) = idt_ccl4

    trac_chemspec(67) = idt_ch3ccl3
    trac_chemspec(68) = idt_hcfc22
    trac_chemspec(69) = idt_hcfc141b
    trac_chemspec(70) = idt_hcfc142b
    trac_chemspec(71) = idt_h1211

    trac_chemspec(72) = idt_ch3br
    trac_chemspec(73) = idt_ch3cl
    trac_chemspec(74) = idt_hcfc21
    trac_chemspec(75) = idt_hcfc123
    trac_chemspec(76) = idt_h2402

    trac_chemspec(77) = idt_chbr3
    trac_chemspec(78) = idt_ch2br2
    trac_chemspec(79) = idt_linearage
    trac_chemspec(80) = idt_idealage

    IF (lpco2) trac_chemspec(n_trac_chemspec) = idt_co2

    ! 2. Families (only used for advection scheme corrections):

    CALL new_tracer('CLY_adv', 'MEZON', idx=idt_cly_adv, nwrite=OFF,  &
         longname='Cly (advection scheme)',     units='mol/mol', &
         code=97, table=131, ninit=CONSTANT, nrerun=OFF, vini=0.0_dp, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('BRY_adv', 'MEZON', idx=idt_bry_adv, nwrite=OFF,  &
         longname='Bry (advection scheme)',     units='mol/mol', &
         code=98, table=131, ninit=CONSTANT, nrerun=OFF, vini=0.0_dp, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    CALL new_tracer('NOY_adv', 'MEZON', idx=idt_noy_adv, nwrite=OFF,  &
         longname='NOy (advection scheme)',     units='mol/mol', &
         code=99, table=131, ninit=CONSTANT, nrerun=OFF, vini=0.0_dp, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    
  END SUBROUTINE request_socol_tracers
  
  
  SUBROUTINE init_socol_tracers
    
    ! Reads initial file of chemistry module MEZON if no restart file and 
    ! initialize tracer tendencies due to chemistry.
 
    ! *init_socol_tracers* is called from *call_init_tracers*, 
    ! src/call_submodels.f90.   

    IMPLICIT NONE

    INTEGER :: jt
    REAL(dp), ALLOCATABLE :: chemfield(:,:,:)
 
    
    ! Executable statements:

    DO jt = 1, n_trac_chemspec
       IF (trlist% ti(trac_chemspec(jt))% init .EQ. 0 .AND. &
            IAND (trlist% ti(trac_chemspec(jt))% ninit, INITIAL) /= 0) THEN

          ! Read initialization fields (each variable individually):
          CALL socol_read_netcdf('chem_initial', &
               TRIM(trlist% ti(trac_chemspec(jt))% fullname), 'LONLATLEV', &
               varname_longname='-', data3d=chemfield)

          ! Allocate field to xt:
          xt(:,:,trac_chemspec(jt),:) = chemfield(:,:,:)

          IF (lresume) xtm1(:,:,trac_chemspec(jt),:) = (1._dp - eps) &
               * xt(:,:,trac_chemspec(jt),:)
          trlist% ti(trac_chemspec(jt))% init = INITIAL

       ENDIF
    ENDDO

    ! Initialize chemical tendencies:
    ALLOCATE(xtte_chem(dcl%nproma,dcl%nlev,ntrac,dcl%ngpblks))
    xtte_chem(:,:,:,:) = 0._dp

    ! Deallocate memory:
    IF (ALLOCATED(chemfield)) DEALLOCATE(chemfield)

    ! Copy initial field of H2O (CTM) to q (GCM) at first time step:
    IF (lstart .AND. lh2o_coupl) q = amw/amd * xt(:,:,idt_h2o,:)

  END SUBROUTINE init_socol_tracers


  SUBROUTINE pos_socol_tracers (kproma, kbdim, klev, ktrac, krow, pxtm1, pxtte)

    ! Forces that chemical tracers remain positive.

    ! *pos_socol_tracers* is called from called from *mezon*, 
    ! src/socol_mezon.f90, and from *call_chem2*, src/call_submodels.f90.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: kproma, kbdim, klev, ktrac, krow
    REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac)
    REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)

    ! Local variables:
    INTEGER  :: jl, jk, jt, cor(kbdim,klev,ktrac), corn(kbdim,klev,ktrac), &
         cormax(ktrac), corfam(kbdim,klev), zero(kbdim,klev), &
         one(kbdim,klev), cornfamb, cornfama
    REAL(dp) :: pxtteb(kbdim,klev,ktrac), pxttec(ktrac), wg(ktrac), wgsum, &
         famteb, famtea, famdiff, sum
    LOGICAL  :: lcor(kbdim,klev)


    ! Executable statements:

    cormax(:) = 0

    pxtteb(1:kproma,:,:) = pxtte(1:kproma,:,:)

    DO jt = 1, n_trac_chemspec

       if (trac_chemspec(jt) .eq. idt_linearage .or. &
            trac_chemspec(jt) .eq. idt_idealage) cycle

       ! Check if tracer remains positive at t+2*dt:
       lcor(1:kproma,:) =  pxtm1(1:kproma,:,trac_chemspec(jt)) + &
            pxtte(1:kproma,:,trac_chemspec(jt))*time_step_len .LT. 0._dp
       cor(1:kproma,:,trac_chemspec(jt)) = MERGE(1, 0, lcor(1:kproma,:))
       cormax(trac_chemspec(jt)) = MAXVAL(cor(1:kproma,:,trac_chemspec(jt)))

       ! Corrections necessary:
       IF (cormax(trac_chemspec(jt)) .EQ. 1) THEN
          DO jl = 1, kproma
             DO jk = 1, klev
                IF (lcor(jl,jk)) THEN

!!$                   ! Print warning:
!!$                   sum = pxtm1(jl,jk,trac_chemspec(jt)) + &
!!$                        time_step_len*pxtte(jl,jk,trac_chemspec(jt))
!!$                   IF (-sum .GE. eps_mr) THEN
!!$                      WRITE(*,*) ' Problems reported in ', &
!!$                           trlist%ti(trac_chemspec(jt))%fullname, jl, krow, jk
!!$                      WRITE(*,*) ' pxtm1+time_step_len*pxtte = ', sum
!!$                   ENDIF

                   pxtte(jl,jk,trac_chemspec(jt)) = &
                        -pxtm1(jl,jk,trac_chemspec(jt))/time_step_len
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    
    ENDDO

    ! If corrections of members of the Cly, Bry and/or NOy families were
    ! necessary, fix the corresponding family:
    IF (MAXVAL(cormax) .EQ. 1) THEN

       ! Cly:
       IF (cormax(idt_cl) + cormax(idt_clo) + cormax(idt_hocl)+ &
            cormax(idt_cl2) + cormax(idt_cl2o2) + cormax(idt_clno3) + &
            cormax(idt_hcl) + cormax(idt_brcl) .GE. 1) THEN

          corfam(1:kproma,:) = cor(1:kproma,:,idt_cl) + &
               cor(1:kproma,:,idt_clo) + cor(1:kproma,:,idt_hocl) + &
               cor(1:kproma,:,idt_cl2) + cor(1:kproma,:,idt_cl2o2) + &
               cor(1:kproma,:,idt_clno3) + cor(1:kproma,:,idt_hcl) + &
               cor(1:kproma,:,idt_brcl)

          DO jl = 1, kproma
             DO jk = 1, klev
                IF (corfam(jl,jk) .GE. 1) THEN

                   ! Tendency of family before correction:
                   famteb = pxtteb(jl,jk,idt_cl) + pxtteb(jl,jk,idt_clo) + &
                        pxtteb(jl,jk,idt_hocl) + &
                        2._dp*pxtteb(jl,jk,idt_cl2) + &
                        2._dp*pxtteb(jl,jk,idt_cl2o2) + &
                        pxtteb(jl,jk,idt_clno3) + &
                        pxtteb(jl,jk,idt_hcl) + pxtteb(jl,jk,idt_brcl)
                   
                   ! Tendency of family after correction:
                   famtea = pxtte(jl,jk,idt_cl) + pxtte(jl,jk,idt_clo) + &
                        pxtte(jl,jk,idt_hocl) + &
                        2._dp*pxtte(jl,jk,idt_cl2) + &
                        2._dp*pxtte(jl,jk,idt_cl2o2) + &
                        pxtte(jl,jk,idt_clno3) + &
                        pxtte(jl,jk,idt_hcl) + pxtte(jl,jk,idt_brcl)

                   ! Error (to be corrected):
                   famdiff = famtea -famteb

                   ! No further corrections if error is too small:
                   IF (ABS(famdiff) .LE. eps_mr) CYCLE

                   corn(1:kproma,:,idt_cl) = 1 - cor(1:kproma,:,idt_cl)
                   corn(1:kproma,:,idt_clo) = 1 - cor(1:kproma,:,idt_clo)
                   corn(1:kproma,:,idt_hocl) = 1 - cor(1:kproma,:,idt_hocl)
                   corn(1:kproma,:,idt_cl2) = 1 - cor(1:kproma,:,idt_cl2)
                   corn(1:kproma,:,idt_cl2o2) = 1 - cor(1:kproma,:,idt_cl2o2)
                   corn(1:kproma,:,idt_hcl) = 1 - cor(1:kproma,:,idt_hcl)
       
                   DO                      
                      ! Weights (only uncorrected values; BrCl and ClNO3 
                      ! excluded, because they belong to more than one family):
                      wg(idt_cl) = ABS(pxtte(jl,jk,idt_cl))* &
                           REAL(corn(jl,jk,idt_cl),dp)
                      wg(idt_clo) = ABS(pxtte(jl,jk,idt_clo))* &
                           REAL(corn(jl,jk,idt_clo),dp)
                      wg(idt_hocl) = ABS(pxtte(jl,jk,idt_hocl))* &
                           REAL(corn(jl,jk,idt_hocl),dp)
                      wg(idt_cl2) = 2._dp*ABS(pxtte(jl,jk,idt_cl2))* &
                           REAL(corn(jl,jk,idt_cl2),dp)
                      wg(idt_cl2o2) = 2._dp*ABS(pxtte(jl,jk,idt_cl2o2))* &
                           REAL(corn(jl,jk,idt_cl2o2),dp)
                      wg(idt_hcl) = ABS(pxtte(jl,jk,idt_hcl))* &
                           REAL(corn(jl,jk,idt_hcl),dp)

                      wgsum = wg(idt_cl) + wg(idt_clo) + wg(idt_hocl) + &
                           wg(idt_cl2) + wg(idt_cl2o2) + wg(idt_hcl)

                      ! Exit endless loop if wgsum is too small:
                      IF (wgsum .LE. eps_mr) THEN
                         cornfama = 0
                         EXIT
                      ENDIF

                      ! Corrections:
                      pxttec(idt_cl) = pxtte(jl,jk,idt_cl) - &
                           wg(idt_cl)/wgsum*famdiff
                      pxttec(idt_clo) = pxtte(jl,jk,idt_clo) - &
                           wg(idt_clo)/wgsum*famdiff
                      pxttec(idt_hocl) = pxtte(jl,jk,idt_hocl) - &
                           wg(idt_hocl)/wgsum*famdiff
                      pxttec(idt_cl2) = pxtte(jl,jk,idt_cl2) - &
                           wg(idt_cl2)/wgsum*famdiff
                      pxttec(idt_cl2o2) = pxtte(jl,jk,idt_cl2o2) - &
                           wg(idt_cl2o2)/wgsum*famdiff
                      pxttec(idt_hcl) = pxtte(jl,jk,idt_hcl) - &
                           wg(idt_hcl)/wgsum*famdiff

                      cornfamb = corn(jl,jk,idt_cl) + corn(jl,jk,idt_clo) + &
                           corn(jl,jk,idt_hocl) + corn(jl,jk,idt_cl2) + &
                           corn(jl,jk,idt_cl2o2) + corn(jl,jk,idt_hcl)

                      ! Do tracers get negative after corrections above?
                      IF (pxtm1(jl,jk,idt_cl) + pxttec(idt_cl)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_cl) = 0
                      IF (pxtm1(jl,jk,idt_clo) + pxttec(idt_clo)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_clo) = 0
                      IF (pxtm1(jl,jk,idt_hocl) + pxttec(idt_hocl)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_hocl) = 0
                      IF (pxtm1(jl,jk,idt_cl2) + pxttec(idt_cl2)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_cl2) = 0
                      IF (pxtm1(jl,jk,idt_cl2o2) + pxttec(idt_cl2o2)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_cl2o2) = 0
                      IF (pxtm1(jl,jk,idt_hcl) + pxttec(idt_hcl)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_hcl) = 0

                      cornfama = corn(jl,jk,idt_cl) + corn(jl,jk,idt_clo) + &
                           corn(jl,jk,idt_hocl) + corn(jl,jk,idt_cl2) + &
                           corn(jl,jk,idt_cl2o2) + corn(jl,jk,idt_hcl)

                      ! Exit endless loop if no further corrections necessary:
                      IF (cornfamb .EQ. cornfama .OR. cornfama .EQ. 0) EXIT

                   END DO

                   ! Correct tendencies:
                   IF (cornfama .NE. 0) THEN
                      pxtte(jl,jk,idt_cl) = pxttec(idt_cl)
                      pxtte(jl,jk,idt_clo) = pxttec(idt_clo)
                      pxtte(jl,jk,idt_hocl) = pxttec(idt_hocl)
                      pxtte(jl,jk,idt_cl2) = pxttec(idt_cl2)
                      pxtte(jl,jk,idt_cl2o2) = pxttec(idt_cl2o2)
                      pxtte(jl,jk,idt_hcl) = pxttec(idt_hcl)
                   ENDIF

                END IF
             END DO
          END DO
       END IF

       ! Bry:
       IF (cormax(idt_br) + cormax(idt_bro) + cormax(idt_hbr) + &
            cormax(idt_hobr) + cormax(idt_brno3) + cormax(idt_brcl) .GE. 1) THEN

          corfam(1:kproma,:) = cor(1:kproma,:,idt_br) + &
               cor(1:kproma,:,idt_bro) + cor(1:kproma,:,idt_hbr) + &
               cor(1:kproma,:,idt_hobr) + cor(1:kproma,:,idt_brno3) + &
               cor(1:kproma,:,idt_brcl)

          DO jl = 1, kproma
             DO jk = 1, klev
                IF (corfam(jl,jk) .GE. 1) THEN

                   ! Tendency of family before correction:
                   famteb = pxtteb(jl,jk,idt_br) + pxtteb(jl,jk,idt_bro) + &
                        pxtteb(jl,jk,idt_hbr) + pxtteb(jl,jk,idt_hobr) + &
                        pxtteb(jl,jk,idt_brno3) + pxtteb(jl,jk,idt_brcl)
                   
                   ! Tendency of family after correction:
                   famtea = pxtte(jl,jk,idt_br) + pxtte(jl,jk,idt_bro) + &
                        pxtte(jl,jk,idt_hbr) + pxtte(jl,jk,idt_hobr) + &
                        pxtte(jl,jk,idt_brno3) + pxtte(jl,jk,idt_brcl)

                   ! Error (to be corrected):
                   famdiff = famtea -famteb
  
                   ! No further corrections if error is too small:
                   IF (ABS(famdiff) .LE. eps_mr) CYCLE

                   corn(1:kproma,:,idt_br) = 1 - cor(1:kproma,:,idt_br)
                   corn(1:kproma,:,idt_bro) = 1 - cor(1:kproma,:,idt_bro)
                   corn(1:kproma,:,idt_hbr) = 1 - cor(1:kproma,:,idt_hbr)
                   corn(1:kproma,:,idt_hobr) = 1 - cor(1:kproma,:,idt_hobr)

                   DO
                      ! Weights (only uncorrected values; BrCl and BrNO3 
                      ! exbruded, because they belong to more than one family):
                      wg(idt_br) = ABS(pxtte(jl,jk,idt_br))* &
                           REAL(corn(jl,jk,idt_br),dp)
                      wg(idt_bro) = ABS(pxtte(jl,jk,idt_bro))* &
                           REAL(corn(jl,jk,idt_bro),dp)
                      wg(idt_hbr) = ABS(pxtte(jl,jk,idt_hbr))* &
                           REAL(corn(jl,jk,idt_hbr),dp)
                      wg(idt_hobr) = 2._dp*ABS(pxtte(jl,jk,idt_hobr))* &
                           REAL(corn(jl,jk,idt_hobr),dp)

                      wgsum = wg(idt_br) + wg(idt_bro) + wg(idt_hbr) + &
                           wg(idt_hobr)

                      ! Exit endless loop if wgsum is too small:
                      IF (wgsum .LE. eps_mr) THEN
                         cornfama = 0
                         EXIT
                      ENDIF

                      ! Corrections:
                      pxttec(idt_br) = pxtte(jl,jk,idt_br) - &
                           wg(idt_br)/wgsum*famdiff
                      pxttec(idt_bro) = pxtte(jl,jk,idt_bro) - &
                           wg(idt_bro)/wgsum*famdiff
                      pxttec(idt_hbr) = pxtte(jl,jk,idt_hbr) - &
                           wg(idt_hbr)/wgsum*famdiff
                      pxttec(idt_hobr) = pxtte(jl,jk,idt_hobr) - &
                           wg(idt_hobr)/wgsum*famdiff

                      cornfamb = corn(jl,jk,idt_br) + corn(jl,jk,idt_bro) + &
                           corn(jl,jk,idt_hbr) + corn(jl,jk,idt_hobr)
                      
                      ! Do tracers get negative after corrections above?
                      IF (pxtm1(jl,jk,idt_br) + pxttec(idt_br)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_br) = 0
                      IF (pxtm1(jl,jk,idt_bro) + pxttec(idt_bro)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_bro) = 0
                      IF (pxtm1(jl,jk,idt_hbr) + pxttec(idt_hbr)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_hbr) = 0
                      IF (pxtm1(jl,jk,idt_hobr) + pxttec(idt_hobr)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_hobr) = 0

                      cornfama = corn(jl,jk,idt_br) + corn(jl,jk,idt_bro) + &
                           corn(jl,jk,idt_hbr) + corn(jl,jk,idt_hobr)

                      ! Exit endless loop if no further corrections necessary:
                      IF (cornfamb .EQ. cornfama .OR. cornfama .EQ. 0) EXIT

                   END DO

                   ! Correct tendencies:
                   IF (cornfama .NE. 0) THEN
                      pxtte(jl,jk,idt_br) = pxttec(idt_br)
                      pxtte(jl,jk,idt_bro) = pxttec(idt_bro)
                      pxtte(jl,jk,idt_hbr) = pxttec(idt_hbr)
                      pxtte(jl,jk,idt_hobr) = pxttec(idt_hobr)
                   ENDIF

                END IF
             END DO
          END DO
       END IF

       ! NOy:
       IF (cormax(idt_n) + cormax(idt_no) + cormax(idt_no2)+ &
            cormax(idt_no3) + cormax(idt_n2o5) + cormax(idt_hno3) + &
            cormax(idt_hno4) + cormax(idt_clno3) + cormax(idt_brno3) .GE. 1) &
            THEN

          corfam(1:kproma,:) = cor(1:kproma,:,idt_n) + &
               cor(1:kproma,:,idt_no) + cor(1:kproma,:,idt_no2) + &
               cor(1:kproma,:,idt_no3) + cor(1:kproma,:,idt_n2o5) + &
               cor(1:kproma,:,idt_hno3) + cor(1:kproma,:,idt_hno4) + &
               cor(1:kproma,:,idt_clno3) + cor(1:kproma,:,idt_brno3) 

          DO jl = 1, kproma
             DO jk = 1, klev
                IF (corfam(jl,jk) .GE. 1) THEN

                   ! Tendency of family before correction:
                   famteb = pxtteb(jl,jk,idt_n) + pxtteb(jl,jk,idt_no) + &
                        pxtteb(jl,jk,idt_no2) + pxtteb(jl,jk,idt_no3) + &
                        2._dp*pxtteb(jl,jk,idt_n2o5) + &
                        pxtteb(jl,jk,idt_hno3) + pxtteb(jl,jk,idt_hno4) + &
                        pxtteb(jl,jk,idt_clno3) + pxtteb(jl,jk,idt_brno3)
                   
                   ! Tendency of family after correction:
                   famtea = pxtte(jl,jk,idt_n) + pxtte(jl,jk,idt_no) + &
                        pxtte(jl,jk,idt_no2) + pxtte(jl,jk,idt_no3) + &
                        2._dp*pxtte(jl,jk,idt_n2o5) + &
                        pxtte(jl,jk,idt_hno3) + pxtte(jl,jk,idt_hno4) + &
                        pxtte(jl,jk,idt_clno3) + pxtte(jl,jk,idt_brno3)

                   ! Error (to be corrected):
                   famdiff = famtea -famteb
  
                   ! No further corrections if error is too small:
                   IF (ABS(famdiff) .LE. eps_mr) CYCLE

                   corn(1:kproma,:,idt_n) = 1 - cor(1:kproma,:,idt_n)
                   corn(1:kproma,:,idt_no) = 1 - cor(1:kproma,:,idt_no)
                   corn(1:kproma,:,idt_no2) = 1 - cor(1:kproma,:,idt_no2)
                   corn(1:kproma,:,idt_no3) = 1 - cor(1:kproma,:,idt_no3)
                   corn(1:kproma,:,idt_n2o5) = 1 - cor(1:kproma,:,idt_n2o5)
                   corn(1:kproma,:,idt_hno3) = 1 - cor(1:kproma,:,idt_hno3)
                   corn(1:kproma,:,idt_hno4) = 1 - cor(1:kproma,:,idt_hno4)

                   DO
                      ! Weights (only uncorrected values; ClNO3 and BrNO3 
                      ! exnuded, because they belong to more than one family):
                      wg(idt_n) = ABS(pxtte(jl,jk,idt_n))* &
                           REAL(corn(jl,jk,idt_n),dp)
                      wg(idt_no) = ABS(pxtte(jl,jk,idt_no))* &
                           REAL(corn(jl,jk,idt_no),dp)
                      wg(idt_no2) = ABS(pxtte(jl,jk,idt_no2))* &
                           REAL(corn(jl,jk,idt_no2),dp)
                      wg(idt_no3) = ABS(pxtte(jl,jk,idt_no3))* &
                           REAL(corn(jl,jk,idt_no3),dp)
                      wg(idt_n2o5) = 2._dp*ABS(pxtte(jl,jk,idt_n2o5))* &
                           REAL(corn(jl,jk,idt_n2o5),dp)
                      wg(idt_hno3) = ABS(pxtte(jl,jk,idt_hno3))* &
                           REAL(corn(jl,jk,idt_hno3),dp)
                      wg(idt_hno4) = ABS(pxtte(jl,jk,idt_hno4))* &
                           REAL(corn(jl,jk,idt_hno4),dp)

                      wgsum = wg(idt_n) + wg(idt_no) + wg(idt_no2) + &
                           wg(idt_no3) + wg(idt_n2o5) + wg(idt_hno3) + &
                           wg(idt_hno4)

                      ! Exit endless loop if wgsum is too small:
                      IF (wgsum .LE. eps_mr) THEN
                         cornfama = 0
                         EXIT
                      ENDIF

                      ! Corrections:
                      pxttec(idt_n) = pxtte(jl,jk,idt_n) - &
                           wg(idt_n)/wgsum*famdiff
                      pxttec(idt_no) = pxtte(jl,jk,idt_no) - &
                           wg(idt_no)/wgsum*famdiff
                      pxttec(idt_no2) = pxtte(jl,jk,idt_no2) - &
                           wg(idt_no2)/wgsum*famdiff
                      pxttec(idt_no3) = pxtte(jl,jk,idt_no3) - &
                           wg(idt_no3)/wgsum*famdiff
                      pxttec(idt_n2o5) = pxtte(jl,jk,idt_n2o5) - &
                           wg(idt_n2o5)/wgsum*famdiff
                      pxttec(idt_hno3) = pxtte(jl,jk,idt_hno3) - &
                           wg(idt_hno3)/wgsum*famdiff
                      pxttec(idt_hno4) = pxtte(jl,jk,idt_hno4) - &
                           wg(idt_hno4)/wgsum*famdiff

                      cornfamb = corn(jl,jk,idt_n) + corn(jl,jk,idt_no) + &
                           corn(jl,jk,idt_no2) + corn(jl,jk,idt_no3) + &
                           corn(jl,jk,idt_n2o5) + corn(jl,jk,idt_hno3) + &
                           corn(jl,jk,idt_hno4)
                      
                      ! Do tracers get negative after corrections above?
                      IF (pxtm1(jl,jk,idt_n) + pxttec(idt_n)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_n) = 0
                      IF (pxtm1(jl,jk,idt_no) + pxttec(idt_no)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_no) = 0
                      IF (pxtm1(jl,jk,idt_no2) + pxttec(idt_no2)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_no2) = 0
                      IF (pxtm1(jl,jk,idt_no3) + pxttec(idt_no3)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_no3) = 0
                      IF (pxtm1(jl,jk,idt_n2o5) + pxttec(idt_n2o5)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_n2o5) = 0
                      IF (pxtm1(jl,jk,idt_hno3) + pxttec(idt_hno3)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_hno3) = 0
                      IF (pxtm1(jl,jk,idt_hno4) + pxttec(idt_hno4)* &
                           time_step_len .LT. 0._dp) corn(jl,jk,idt_hno4) = 0

                      cornfama = corn(jl,jk,idt_n) + corn(jl,jk,idt_no) + &
                           corn(jl,jk,idt_no2) + corn(jl,jk,idt_no3) + &
                           corn(jl,jk,idt_n2o5) + corn(jl,jk,idt_hno3) + &
                           corn(jl,jk,idt_hno4)

                      ! Exit endless loop if no further corrections necessary:
                      IF (cornfamb .EQ. cornfama .OR. cornfama .EQ. 0) EXIT

                   END DO

                   ! Correct tendencies:
                   IF (cornfama .NE. 0) THEN
                      pxtte(jl,jk,idt_n) = pxttec(idt_n)
                      pxtte(jl,jk,idt_no) = pxttec(idt_no)
                      pxtte(jl,jk,idt_no2) = pxttec(idt_no2)
                      pxtte(jl,jk,idt_no3) = pxttec(idt_no3)
                      pxtte(jl,jk,idt_n2o5) = pxttec(idt_n2o5)
                      pxtte(jl,jk,idt_hno3) = pxttec(idt_hno3)
                      pxtte(jl,jk,idt_hno4) = pxttec(idt_hno4)
                   ENDIF

                END IF
             END DO
          END DO
       END IF

    END IF
 
  END SUBROUTINE pos_socol_tracers

 
  SUBROUTINE destruct_socol_tracers
    
    ! Deallocates module variables.
    
    ! *destruct_socol_tracers* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(xtte_chem)) DEALLOCATE(xtte_chem)

  END SUBROUTINE destruct_socol_tracers

END MODULE mo_socol_tracers

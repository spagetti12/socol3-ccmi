 SUBROUTINE chem(krow, xt, ph)
 
   ! Description:
   !
   ! Chemistry package Ver. 3.1 (May.,1999)
   ! Main features:
   ! 1. Twenty five species
   ! 2. Look-up table parameterization for photodissociation rates
   ! 3. Numerical scheme: Implicit Newton-Raphson scheme
   !
   ! *chem* is called from *mezon*.
   !
   ! Eugene Rozanov, PMOD/WRC Davos, original code
   ! Martin Schraner, ETH Zurich, April 2009: Modifications for SOCOLvs3.0
 
   USE mo_exception,               ONLY: finish, message, message_text
   USE mo_kind,                    ONLY: dp
   USE mo_socol_chem,              ONLY: fico, cofi
   USE mo_socol_constants,         ONLY: eps_mr
   USE mo_socol_dimensions,        ONLY: nbdim, nproma, nlev, ntrac, &
                                         ndiso, nhiso, ngasch
   USE mo_socol_gcmfields,         ONLY: t, prest
   USE mo_socol_ghg_ods,           ONLY: ods_fa2mem_chem, odsnat_chem
   USE mo_socol_grid_calculations, ONLY: dens, zetb
   USE mo_socol_hetero,            ONLY: het_new
   USE mo_socol_namelist,          ONLY: lphotfull, lhetchem, lch4_isotope, lo3orig 
   USE mo_socol_strataerosols,     ONLY: sad_strataer, rm_strataer, av_strataer
   USE mo_socol_streams,           ONLY: rxn_1, rxn_2, rxn_3, rxn_4, &
        rxn_5, rxn_6, rxn_7, rxn_8, &
        rxn_9, rxn_10, rxn_11, rxn_12, &
        rxn_13, rxn_14, rxn_15, rxn_16, &
        rxn_17, rxn_18, rxn_19, rxn_20, &
        rxn_21, rxn_22, rxn_23, rxn_24, &
        rxn_25, rxn_26, rxn_27, rxn_28, &
        rxn_29, rxn_30, rxn_31, rxn_32, &
        rxn_33, rxn_34, &
        jo2, jo3d, jno2, jcl2o2 !!$, & 
!!$        jpan, jpaa, jmgly

   USE mo_socol_sun,               ONLY: tjo2a,      tjo3p,      tjo3d,      &
                                         tjno,       tjno2,      tjhno3,     &
                                         tjno3a,     tjno3b,     tjn2o5,     &
                                         tjn2o,      tjhno4,     tjclno3a,   &
                                         tjh2o2,     tjf11,      tjf12,      &
                                         tjhocl,     tjh2o,      tjco2,      &
                                         tjcl2,      tjcl2o2,    tjch2oa,    &
                                         tjch2ob,    tjch3o2h,   tjbro,      &
                                         tjbrno3a,   tjbrno3b,   tjbrcl,     &
                                         tjhobr,     tjcbrf3,    tjo2b,      &
                                         tjch4,      tjclno3b,   tjhcl,      &
                                         tjcfc113,   tjcfc114,   tjcfc115,   &
                                         tjccl4,     tjch3ccl3,  tjhcfc22,   &
                                         tjhcfc141b, tjhcfc142b, tjh1211,    &
                                         tjch3br,    tjch3cl,    tjhcfc21,   &
                                         tjhcfc123,  tjh2402,    tjchbr3,    &
                                         tjch2br2,                           &
                                         xn2o5, xn2o, xf11, xf12,            &
                                         tjpan, tjmacr, tjhac, tjmgly, tjpaa

   USE mo_socol_isotope,           ONLY: idt_12c, idt_13c ! eth_as_ch4
   USE mo_socol_photo,             ONLY: pan_jval, ch3co3h_jval, mgly_jval
   USE mo_socol_time_control,      ONLY: time_step_len_chem
   USE mo_socol_tracers,           ONLY: idt_O3      ,idt_OD      ,idt_O       ,&
                                         idt_NO      ,idt_HO2     ,idt_CLO     ,&
                                         idt_NO2     ,idt_OH      ,idt_NO3     ,&
                                         idt_N2O5    ,idt_HNO3    ,idt_HNO4    ,&
                                         idt_CLNO3   ,idt_CL      ,idt_N       ,&
                                         idt_H2O2    ,idt_H       ,idt_HOCL    ,&
                                         idt_CL2     ,idt_CL2O2   ,idt_HCL     ,&
                                         idt_BR      ,idt_CH2O    ,idt_BRO     ,&
                                         idt_HBR     ,idt_HOBR    ,idt_BRNO3   ,&
                                         idt_BRCL    ,idt_CH3     ,idt_CH3O2   ,&
                                         idt_CH3O    ,idt_HCO     ,idt_CH3O2H  ,&
                                         idt_C5H8    ,idt_ISO2    ,idt_ISO2H   ,&
                                         idt_ISON    ,idt_MACR    ,idt_MACRO2  ,&
                                         idt_MPAN    ,idt_MACRO2H ,idt_HACET   ,&
                                         idt_MGLY    ,idt_NALD    ,idt_CH3CO3  ,&
                                         idt_PAN     ,idt_CH3CO3H ,idt_CH3COOH ,&
                                         idt_HCOOH   ,&
                                         idt_H2O     ,idt_N2O     ,idt_CH4     ,&
                                         idt_CO      ,idt_H2      ,&
                                         idt_odscls  ,idt_odscll  ,idt_odsbr   ,&
                                         idt_psc1    ,idt_psc2, &
                                         idt_f11, idt_f12, idt_cbrf3, idt_cfc113, idt_cfc114, &
                                         idt_cfc115, idt_ccl4, idt_ch3ccl3, idt_hcfc22, idt_hcfc141b, &
                                         idt_hcfc142b, idt_h1211, idt_ch3br, idt_ch3cl, idt_hcfc21, &
                                         idt_hcfc123, idt_h2402, idt_chbr3, idt_ch2br2

   USE mo_socol_photo,              ONLY: cloud_mod, cloud_mul_d ! eth_as_photo
   USE mo_geoloc,                   ONLY: amu0_x
   USE mo_memory_g1a,               ONLY: xlm1, xim1
   USE mo_memory_g3b,               ONLY: albedo, aclc
   USE mo_socol_ch4_streams,        ONLY: ch4descl_d, ch4desod_1_d, ch4desod_2_d, &
                                          ch4desoh_d, ch4dest01_d, &
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
                                          ohdest_ods19_d, &
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
                                          ho2prod_rt19_d
   USE mo_socol_o3orig,             ONLY: o3losst, o3prod

   IMPLICIT NONE
 
   ! Subroutine arguments:
   INTEGER, INTENT(in) :: krow  ! Index of current grid point block
   REAL(dp), INTENT(in) :: ph(nbdim,nlev+1) ! eth_as_photo
   REAL(dp), INTENT(inout) :: xt(nbdim,nlev,ntrac)      ! Input / output:
                                ! chemical species before / after chemistry
                                ! with time step length = time_step_len_chem
 
   ! Local variables:
   INTEGER :: jl, jk, it, k, nakop, i
 
   REAL(dp) :: xxd, pw, tw, sadw, r_mean, avol, tmww
   REAL(dp) :: brsuma, brsumb, clsuma, clsumb, brsumaa, brsumbb, clsumaa, &
        clsumbb, noysuma, noysumb, clprodods, brprodods, noyloss
   REAL(dp) :: gasch94_cl, gasch94_b
   REAL(dp) :: diso(ndiso), hiso(nhiso), gasch(ngasch)
 
   LOGICAL :: kluerrb, kluerrc, kluerrn, lexitloop
 
   ! eth_as_photo+

   logical  ::  zagtz(nbdim)             ! zenith angle > 0 flag array
   real(dp) ::  secant  
   real(dp), dimension(nlev) :: cld_line, &  ! vertical cloud array
        lwc_line, &            ! vertical lwc array
        fac1, &                
        eff_alb, &             ! effective albedo from cloud modifications
        cld_mult               ! clould multiplier
   real(dp):: clouds(nbdim,nlev) ! cloud fraction
   
   ! eth_as_photo-

   ! Local parameters:
   INTEGER, PARAMETER :: nakopmax = 128   ! Maximal number to devide chemical
                                          ! time step for Newton-Raphson  !MSv3
   INTEGER, PARAMETER :: nuit = 100       ! Maximal number of iterations for
                                          ! Newton-Raphson (subroutine
                                          ! *newraph*)
   REAL(dp), PARAMETER :: eeps = 0.05_dp  ! Maximal relative error allowed in
                                          ! Newton-Raphson
   ! External subroutines:
   EXTERNAL dis, newraph
 
 
   ! Executable statements:
 
   ! eth_as_photo
   ! attention: amu0_x is cosine of zenith angle
   zagtz(:) = amu0_x(:,krow) > 0._dp   

   DO jl = 1, nproma
 
      IF (lphotfull) THEN
 
         ! Calculate photolysis rates (lphotofull=.TRUE.)
         ! or get photolysis rates
         ! from look-up table (lphotofull=.FALSE.):
 
         ! CALL dispho
 
         WRITE(message_text,*) &
              'Calculation of full photolysis not implemented yet!', &
              'Set lphotfull=.FALSE. in namelist SOCOLCTL!'
         CALL message('',TRIM(message_text))
         CALL finish('socol_chem','Run terminated')
 
      ELSE
 
         ! Get photolyis rates from look up table according to the oxygen and
         ! ozone amount in the air column above a grid box:
         CALL dis(jl, krow, xt(jl,:,idt_o3))
 
      ENDIF

      ! eth_as_photo+
      
      cld_mult(:) = 1._dp
      
      if( zagtz(jl) ) then

         secant = 1._dp / amu0_x(jl,krow)
         
         if( secant <= 50._dp ) then
            do jk = 1, nlev
               fac1(jk)     = ph(jl,jk+1) - ph(jl,jk)
            end do
            lwc_line(:) = xlm1(jl,:,krow) + xim1(jl,:,krow)
            cld_line(:) = aclc(jl,:,krow)
           
            call cloud_mod( amu0_x(jl,krow), cld_line, lwc_line, fac1, albedo(jl,krow), &
                 eff_alb, cld_mult)   
            
         end if
      end if

     ! eth_as_photo-
 
      ! Perform chemical calculations stepwise for each grid box:
      DO jk = 1, nlev

         cloud_mul_d(jl,jk,krow) = cld_mult(jk) ! eth_as_photo
 
         ! Temperature, pressure, molecule number density:                !MSv3
         tw = t(jl,jk)                                              ! K
         pw = prest(jl,jk)                                          ! hPa
         xxd = dens(jl,jk)                                          ! mol/cm**3
 
         ! Stratospheric aerosols:
         sadw = sad_strataer(jl,jk)  ! surface area density  [cm2/cm3]    !MSv3
!!$         xnp  = nd_strataer(jl,jk)   ! number density [1/cm3]             !MSv3
         r_mean= rm_strataer(jl,jk)   ! mean radius [um]
         avol  = av_strataer(jl,jk)   ! volume density [cm3/cm3]
        
         ! Rename photolysis rates and multiply them by delta-e corrections (if
         ! necessary):
         ! eth_as_31052013: multiplication by delta-e corrections eliminated, not necessary any more
         diso(1)  = cld_mult(jk) * tjo2a(jk)
         diso(2)  = cld_mult(jk) * tjo3p(jk)   ! * 0.9_dp
         diso(3)  = cld_mult(jk) * tjo3d(jk)
         diso(4)  = cld_mult(jk) * tjno(jk)
         diso(5)  = cld_mult(jk) * tjno2(jk)   ! * 0.9_dp
 
         diso(6)  = cld_mult(jk) * tjhno3(jk)  ! * 1.2_dp
         diso(7)  = cld_mult(jk) * tjno3a(jk)  ! * 0.9_dp
         diso(8)  = cld_mult(jk) * tjno3b(jk)  ! * 0.9_dp
         diso(9)  = cld_mult(jk) * tjn2o5(jk) ! * xn2o5(jk) ! delta-e-correction eliminated
         diso(10) = cld_mult(jk) * tjn2o(jk) ! * xn2o(jk) ! delta-e-correction eliminated
 
         diso(11) = cld_mult(jk) * tjhno4(jk)
         diso(12) = cld_mult(jk) * tjclno3a(jk)
         diso(13) = cld_mult(jk) * tjh2o2(jk)
         diso(14) = cld_mult(jk) * tjh2o(jk)
         diso(15) = cld_mult(jk) * tjhocl(jk)  ! * 0.9_dp
 
         diso(16) = cld_mult(jk) * tjcl2(jk)    ! * 0.9_dp
         diso(17) = cld_mult(jk) * tjcl2o2(jk)
         diso(18) = cld_mult(jk) * tjbro(jk)    ! * 0.9_dp
         diso(19) = cld_mult(jk) * tjbrno3a(jk) ! * 0.9_dp
         diso(20) = cld_mult(jk) * tjbrno3b(jk) ! * 0.9_dp
 
         diso(21) = cld_mult(jk) * tjbrcl(jk)   ! * 0.9_dp
         diso(22) = cld_mult(jk) * tjhobr(jk)  ! * 0.9_dp
         diso(23) = cld_mult(jk) * tjch2oa(jk) ! * 0.9_dp
         diso(24) = cld_mult(jk) * tjch2ob(jk) ! * 0.9_dp
         diso(25) = cld_mult(jk) * tjco2(jk)
 
         diso(26) = cld_mult(jk) * tjch3o2h(jk)
         diso(27) = cld_mult(jk) * tjch4(jk)
         diso(28) = cld_mult(jk) * tjo2b(jk)
         diso(29) = cld_mult(jk) * tjclno3b(jk)
         diso(30) = cld_mult(jk) * tjhcl(jk)   ! * 1.4_dp
 
         diso(31) = cld_mult(jk) * tjch3o2h(jk)        ! ISO2H
         diso(32) = cld_mult(jk) * 3.7_dp * tjpan(jk)  ! ISON
         diso(33) = cld_mult(jk) * tjmacr(jk)   ! MACR
         diso(34) = cld_mult(jk) * tjpan(jk)     ! MPAN
         diso(35) = cld_mult(jk) * tjch3o2h(jk)  ! MACRO2H

         diso(36) = cld_mult(jk) * tjhac(jk) ! HACET
         diso(37) = cld_mult(jk) * tjmgly(jk)  ! MGLY
         diso(38) = cld_mult(jk) * 0.19_dp * tjch2ob(jk)  ! NALD
         diso(39) = cld_mult(jk) * tjpan(jk) ! PAN
         diso(40) = cld_mult(jk) * tjpaa(jk) ! CH3CO3H
 
         diso(41) = cld_mult(jk) * tjf11(jk)  ! * xf11(jk) ! delta-e-correction eliminated
         diso(42) = cld_mult(jk) * tjf12(jk)  ! * xf12(jk) ! delta-e-correction eliminated
         diso(43) = cld_mult(jk) * tjcbrf3(jk)
         diso(44) = cld_mult(jk) * tjcfc113(jk)
         diso(45) = cld_mult(jk) * tjcfc114(jk)
 
         diso(46) = cld_mult(jk) * tjcfc115(jk)
         diso(47) = cld_mult(jk) * tjccl4(jk)
         diso(48) = cld_mult(jk) * tjch3ccl3(jk)
         diso(49) = cld_mult(jk) * tjhcfc22(jk)
         diso(50) = cld_mult(jk) * tjhcfc141b(jk)
 
         diso(51) = cld_mult(jk) * tjhcfc142b(jk)
         diso(52) = cld_mult(jk) * tjh1211(jk)
         diso(53) = cld_mult(jk) * tjch3br(jk)
         diso(54) = cld_mult(jk) * tjch3cl(jk)
         diso(55) = cld_mult(jk) * tjhcfc21(jk)
 
         diso(56) = cld_mult(jk) * tjhcfc123(jk)
         diso(57) = cld_mult(jk) * tjh2402(jk)
         diso(58) = cld_mult(jk) * tjchbr3(jk)
         diso(59) = cld_mult(jk) * tjch2br2(jk)
 
         ! Initial value for nakop (number to sub-divide model time step length
         ! for chemical calculations):
         nakop = 1
 
         ! Initial value for lexitloop (flag to exit the endless loop below):
         lexitloop = .FALSE.
 
         ! The following endless loop is passed exactly once, if the mass of
         ! Bry, Cly, and NOy remain conserved after chemical calculations.
         ! Otherwise, the timestep for chemical calculations is halved and all
         ! calculations are repeated, i.e. the loop is passed a second time.
         ! This procedure goes on until Bry, Cly, and NOy are conserved or
         ! nakop > nakopmax.
         DO                                                               !MSv3
 
            ! Sub-devide model time step in nakop time steps of length tmww [s]:
            tmww = time_step_len_chem/REAL(nakop,dp)  ! nakop sould be 1, 2, 4,
                                                      ! 8, 16,...
 
            ! Initial values for flags of violation of mass conservation for
            ! Bry, Cly and NOy (if .FALSE. mass remains conserved after chemical
            ! calculations, if .TRUE. it does not):
            kluerrb = .FALSE.  ! Flag of mass conservation of Bry         !MSv3
            kluerrc = .FALSE.  ! Flag of mass conservation of Cly         !MSv3
            kluerrn = .FALSE.  ! Flag of mass conservation of NOy         !MSv3
 
            ! Tranform chemical species from mixing ratio to mol/cm**3:   !MSv3
            gasch( 1) = xt(jl,jk,idt_O3      ) * xxd   !  mol/cm**3
            gasch( 2) = xt(jl,jk,idt_OD      ) * xxd   !  mol/cm**3
            gasch( 3) = xt(jl,jk,idt_O       ) * xxd   !  mol/cm**3
            gasch( 4) = xt(jl,jk,idt_NO      ) * xxd   !  mol/cm**3
            gasch( 5) = xt(jl,jk,idt_HO2     ) * xxd   !  mol/cm**3
            gasch( 6) = xt(jl,jk,idt_CLO     ) * xxd   !  mol/cm**3
            gasch( 7) = xt(jl,jk,idt_NO2     ) * xxd   !  mol/cm**3
            gasch( 8) = xt(jl,jk,idt_OH      ) * xxd   !  mol/cm**3
            gasch( 9) = xt(jl,jk,idt_NO3     ) * xxd   !  mol/cm**3
            gasch(10) = xt(jl,jk,idt_N2O5    ) * xxd   !  mol/cm**3
            gasch(11) = xt(jl,jk,idt_HNO3    ) * xxd   !  mol/cm**3
            gasch(12) = xt(jl,jk,idt_HNO4    ) * xxd   !  mol/cm**3
            gasch(13) = xt(jl,jk,idt_CLNO3   ) * xxd   !  mol/cm**3
            gasch(14) = xt(jl,jk,idt_CL      ) * xxd   !  mol/cm**3
            gasch(15) = xt(jl,jk,idt_N       ) * xxd   !  mol/cm**3
            gasch(16) = xt(jl,jk,idt_H2O2    ) * xxd   !  mol/cm**3
            gasch(17) = xt(jl,jk,idt_H       ) * xxd   !  mol/cm**3
            gasch(18) = xt(jl,jk,idt_HOCL    ) * xxd   !  mol/cm**3
            gasch(19) = xt(jl,jk,idt_CL2     ) * xxd   !  mol/cm**3
            gasch(20) = xt(jl,jk,idt_CL2O2   ) * xxd   !  mol/cm**3
            gasch(21) = xt(jl,jk,idt_HCL     ) * xxd   !  mol/cm**3
            gasch(22) = xt(jl,jk,idt_BR      ) * xxd   !  mol/cm**3
            gasch(23) = xt(jl,jk,idt_CH2O    ) * xxd   !  mol/cm**3
            gasch(24) = xt(jl,jk,idt_BRO     ) * xxd   !  mol/cm**3
            gasch(25) = xt(jl,jk,idt_HBR     ) * xxd   !  mol/cm**3
            gasch(26) = xt(jl,jk,idt_HOBR    ) * xxd   !  mol/cm**3
            gasch(27) = xt(jl,jk,idt_BRNO3   ) * xxd   !  mol/cm**3
            gasch(28) = xt(jl,jk,idt_BRCL    ) * xxd   !  mol/cm**3
            gasch(29) = xt(jl,jk,idt_CH3     ) * xxd   !  mol/cm**3
            gasch(30) = xt(jl,jk,idt_CH3O2   ) * xxd   !  mol/cm**3
            gasch(31) = xt(jl,jk,idt_CH3O    ) * xxd   !  mol/cm**3
            gasch(32) = xt(jl,jk,idt_HCO     ) * xxd   !  mol/cm**3
            gasch(33) = xt(jl,jk,idt_CH3O2H  ) * xxd   !  mol/cm**3
            gasch(34) = xt(jl,jk,idt_C5H8    ) * xxd   !  mol/cm**3
            gasch(35) = xt(jl,jk,idt_ISO2    ) * xxd   !  mol/cm**3
            gasch(36) = xt(jl,jk,idt_ISO2H   ) * xxd   !  mol/cm**3
            gasch(37) = xt(jl,jk,idt_ISON    ) * xxd   !  mol/cm**3
            gasch(38) = xt(jl,jk,idt_MACR    ) * xxd   !  mol/cm**3
            gasch(39) = xt(jl,jk,idt_MACRO2  ) * xxd   !  mol/cm**3
            gasch(40) = xt(jl,jk,idt_MPAN    ) * xxd   !  mol/cm**3
            gasch(41) = xt(jl,jk,idt_MACRO2H ) * xxd   !  mol/cm**3
            gasch(42) = xt(jl,jk,idt_HACET   ) * xxd   !  mol/cm**3
            gasch(43) = xt(jl,jk,idt_MGLY    ) * xxd   !  mol/cm**3
            gasch(44) = xt(jl,jk,idt_NALD    ) * xxd   !  mol/cm**3
            gasch(45) = xt(jl,jk,idt_CH3CO3  ) * xxd   !  mol/cm**3
            gasch(46) = xt(jl,jk,idt_PAN     ) * xxd   !  mol/cm**3
            gasch(47) = xt(jl,jk,idt_CH3CO3H ) * xxd   !  mol/cm**3
            gasch(48) = xt(jl,jk,idt_CH3COOH ) * xxd   !  mol/cm**3
            gasch(49) = xt(jl,jk,idt_HCOOH   ) * xxd   !  mol/cm**3
            gasch(50) = xt(jl,jk,idt_H2O     ) * xxd   !  mol/cm**3
            gasch(51) = xt(jl,jk,idt_N2O     ) * xxd   !  mol/cm**3
            gasch(52) = xt(jl,jk,idt_CH4     ) * xxd   !  mol/cm**3
            gasch(53) = xt(jl,jk,idt_CO      ) * xxd   !  mol/cm**3
            gasch(54) = xt(jl,jk,idt_H2      ) * xxd   !  mol/cm**3
!!$            gasch(74) = xt(jl,jk,idt_odscls  ) * xxd   !  mol/cm**3
!!$            gasch(75) = xt(jl,jk,idt_odscll  ) * xxd   !  mol/cm**3
!!$            gasch(76) = xt(jl,jk,idt_odsbr   ) * xxd   !  mol/cm**3
            gasch(77) = 0.0_dp ! psc1
            gasch(78) = 0.0_dp ! psc2

            ! eth_as_ch4+
            IF (lch4_isotope) THEN
               gasch(120) = xt(jl,jk,idt_12c)  *xxd
               gasch(121) = xt(jl,jk,idt_13c)  *xxd 
            END IF
            ! eth_as_ch4+  
 
!!$            ! Determine individual ODS species from families:             !MSv3
!!$            gasch(55) = ods_fa2mem_chem(jl,jk, 1)*gasch(74)  !MSODS
!!$            gasch(56) = ods_fa2mem_chem(jl,jk, 2)*gasch(75)  !MSODS
!!$            gasch(57) = ods_fa2mem_chem(jl,jk, 3)*gasch(76)  !MSODS
!!$            gasch(58) = ods_fa2mem_chem(jl,jk, 4)*gasch(75)  !MSODS
!!$            gasch(59) = ods_fa2mem_chem(jl,jk, 5)*gasch(75)  !MSODS
!!$            gasch(60) = ods_fa2mem_chem(jl,jk, 6)*gasch(75)  !MSODS
!!$            gasch(61) = ods_fa2mem_chem(jl,jk, 7)*gasch(74)  !MSODS
!!$            gasch(62) = ods_fa2mem_chem(jl,jk, 8)*gasch(74)  !MSODS
!!$            gasch(63) = ods_fa2mem_chem(jl,jk, 9)*gasch(74)  !MSODS
!!$            gasch(64) = ods_fa2mem_chem(jl,jk,10)*gasch(74)  !MSODS
!!$            gasch(65) = ods_fa2mem_chem(jl,jk,11)*gasch(74)  !MSODS
!!$            gasch(66) = ods_fa2mem_chem(jl,jk,12)*gasch(76)  !MSODS
!!$            gasch(67) = ods_fa2mem_chem(jl,jk,13)*gasch(76)  !MSODS
!!$            gasch(68) = ods_fa2mem_chem(jl,jk,14)*gasch(74)  !MSODS
!!$            gasch(69) = ods_fa2mem_chem(jl,jk,15)*gasch(74)  !MSODS
!!$            gasch(70) = ods_fa2mem_chem(jl,jk,16)*gasch(74)  !MSODS
!!$            gasch(71) = ods_fa2mem_chem(jl,jk,17)*gasch(76)  !MSODS
!!$            gasch(72) = ods_fa2mem_chem(jl,jk,18)*gasch(76)  !MSODS
!!$            gasch(73) = ods_fa2mem_chem(jl,jk,19)*gasch(76)  !MSODS
!!$            gasch94_cl= ods_fa2mem_chem(jl,jk,20)*gasch(74)   !H1211_CL
!!$            gasch94_b = gasch(66)
!!$ !!$
            ! ODS species from tracers:                 ! eth_as_odslt
            gasch(55) = xt(jl,jk,idt_f11      ) * xxd   !  mol/cm**3
            gasch(56) = xt(jl,jk,idt_f12      ) * xxd   !  mol/cm**3
            gasch(57) = xt(jl,jk,idt_cbrf3    ) * xxd   !  mol/cm**3 
            gasch(58) = xt(jl,jk,idt_cfc113   ) * xxd   !  mol/cm**3 
            gasch(59) = xt(jl,jk,idt_cfc114   ) * xxd   !  mol/cm**3 
            gasch(60) = xt(jl,jk,idt_cfc115   ) * xxd   !  mol/cm**3 
            gasch(61) = xt(jl,jk,idt_ccl4     ) * xxd   !  mol/cm**3 
            gasch(62) = xt(jl,jk,idt_ch3ccl3  ) * xxd   !  mol/cm**3 
            gasch(63) = xt(jl,jk,idt_hcfc22   ) * xxd   !  mol/cm**3 
            gasch(64) = xt(jl,jk,idt_hcfc141b ) * xxd   !  mol/cm**3 
            gasch(65) = xt(jl,jk,idt_hcfc142b ) * xxd   !  mol/cm**3
            gasch(66) = xt(jl,jk,idt_h1211    ) * xxd   !  mol/cm**3 
            gasch(67) = xt(jl,jk,idt_ch3br    ) * xxd   !  mol/cm**3 
            gasch(68) = xt(jl,jk,idt_ch3cl    ) * xxd   !  mol/cm**3 
            gasch(69) = xt(jl,jk,idt_hcfc21   ) * xxd   !  mol/cm**3 
            gasch(70) = xt(jl,jk,idt_hcfc123  ) * xxd   !  mol/cm**3 
            gasch(71) = xt(jl,jk,idt_h2402    ) * xxd   !  mol/cm**3 
            gasch(72) = xt(jl,jk,idt_chbr3    ) * xxd   !  mol/cm**3 
            gasch(73) = xt(jl,jk,idt_ch2br2   ) * xxd   !  mol/cm**3  

            ! Bry, Cly, NOy before chemical calculations:
!!$            brsumb = (xt(jl,jk,idt_br)+xt(jl,jk,idt_bro)+xt(jl,jk,idt_hbr)+ &
!!$                      xt(jl,jk,idt_hobr)+xt(jl,jk,idt_brno3)+ &
!!$                      xt(jl,jk,idt_brcl))*xxd  !without organic bromines
            brsumb = gasch(22) + gasch(24) + gasch(25)+ gasch(26) + &
                 gasch(27) + gasch(28)
 
!!$            brsumbb = (xt(jl,jk,idt_br)+xt(jl,jk,idt_bro)+xt(jl,jk,idt_hbr)+ &
!!$                       xt(jl,jk,idt_hobr)+xt(jl,jk,idt_brno3)+ &
!!$                       xt(jl,jk,idt_brcl)+xt(jl,jk,idt_odsbr))*xxd
             brsumbb = gasch(22) + gasch(24) + gasch(25) + gasch(26) + &
                  gasch(27) + gasch(28) + gasch(57) + gasch(66)+ &
                  gasch(67) + 2._dp*gasch(71) + 3._dp*gasch(72)+ 2._dp*gasch(73)

!!$            clsumb = (xt(jl,jk,idt_clo)+xt(jl,jk,idt_clno3)+xt(jl,jk,idt_cl)+ &
!!$                      xt(jl,jk,idt_hocl)+2._dp*xt(jl,jk,idt_cl2)+ &
!!$                      2._dp*xt(jl,jk,idt_cl2o2)+xt(jl,jk,idt_hcl)+ &
!!$                      xt(jl,jk,idt_brcl))*xxd  !without inorganic chlorines
            clsumb = gasch( 6) + gasch(13) + gasch(14) + gasch(18) + &
                 2._dp*gasch(19) + 2._dp*gasch(20) + gasch(21) + gasch(28)

!!$            clsumbb = (xt(jl,jk,idt_clo)+xt(jl,jk,idt_clno3)+xt(jl,jk,idt_cl)+ &
!!$                      xt(jl,jk,idt_hocl)+2._dp*xt(jl,jk,idt_cl2)+ &
!!$                      2._dp*xt(jl,jk,idt_cl2o2)+xt(jl,jk,idt_hcl)+ &
!!$                      xt(jl,jk,idt_brcl)+xt(jl,jk,idt_odscls)+ &
!!$                      xt(jl,jk,idt_odscll))*xxd
            clsumbb = gasch( 6) + gasch(13) + gasch(14) + gasch(18) + &
                 2._dp*gasch(19) + 2._dp*gasch(20) + gasch(21) + gasch(28) + &
                 3._dp*gasch(55) + 2._dp*gasch(56) + 3._dp*gasch(58) + &
                 2._dp*gasch(59) + gasch(60) + 4._dp*gasch(61) + 3._dp*gasch(62) + &
                 gasch(63) + 2._dp*gasch(64) + gasch(65) + gasch(66) + &
                 gasch(68) + 2._dp*gasch(69) + 2._dp*gasch(70)
 
            noysumb = (xt(jl,jk,idt_no)+xt(jl,jk,idt_no2)+ &
                       2._dp*xt(jl,jk,idt_n2o5)+xt(jl,jk,idt_hno3)+ &
                       xt(jl,jk,idt_hno4)+xt(jl,jk,idt_no3)+ &
                       xt(jl,jk,idt_clno3)+xt(jl,jk,idt_n)+ &
                       xt(jl,jk,idt_brno3)+2._dp*xt(jl,jk,idt_n2o)+ &
                       xt(jl,jk,idt_ison)+xt(jl,jk,idt_mpan)+ &
                       xt(jl,jk,idt_nald)+xt(jl,jk,idt_pan))*xxd
 
 
            ! Initial values for Cl / Br produced from ODSs and nitrogen oxide
            ! lost by canibalistic reactions:
            clprodods = 0.0_dp                                          !MSCLBI
            brprodods = 0.0_dp                                          !MSCLBI
            noyloss   = 0.0_dp                                          !MSv3
 
            ! Initial values of rates of heterogeneous reactions with H2O:
            ! STS:
            hiso(1) = 0.0_dp                 ! N2O5   + H2O -> 2HNO3
            hiso(2) = 0.0_dp                 ! ClONO2 + H2O -> HOCl + HNO3
            hiso(3) = 0.0_dp                 ! BrONO2 + H2O -> HOBr + HNO3
            ! PSC I:
            hiso(4) = 0.0_dp                 ! N2O5   + H2O -> 2HNO3
            hiso(5) = 0.0_dp                 ! ClONO2 + H2O -> HOCl + HNO3
            ! PSC II:
            hiso(6) = 0.0_dp                 ! N2O5   + H2O -> 2HNO3
            hiso(7) = 0.0_dp                 ! ClONO2 + H2O -> HOCl + HNO3
            hiso(8) = 0.0_dp                 ! BrONO2 + H2O -> HOBr + HNO3
            ! Initial values of rates of heterogeneous reactions with HCl:
            ! STS:
            hiso(9)  = 0.0_dp                ! ClONO2 + HCl -> Cl2  + HNO3
            hiso(10) = 0.0_dp                ! HOCl   + HCl -> Cl2  + H2O
            hiso(11) = 0.0_dp                ! HOBr   + HCl -> BrCl + H2O
            ! PSC I:
            hiso(12) = 0.0_dp                ! ClONO2 + HCl -> Cl2  + HNO3
            hiso(13) = 0.0_dp                ! HOCl   + HCl -> Cl2  + H2O
            ! PSC II:
            hiso(14) = 0.0_dp                ! ClONO2 + HCl -> Cl2  + HNO3
            hiso(15) = 0.0_dp                ! HOCl   + HCl -> Cl2  + H2O
            hiso(16) = 0.0_dp                ! HOBr   + HCl -> BrCl + H2O
 
            ! Calculate heterogeneous reaction rates (hiso), new concentrations
            ! of H2O(g), H2O(s), HNO3(g), and HNO3(s) (if PSC I and/or PSC II
            ! are available) and determine type of liquid aerosols (if
            ! available):
            ! Note: hiso(2), hiso(9), hiso(10), and hiso(11) are NOT calculated
            !       by *het_new*, but are determined later in subroutine
            !       *calc_hetre* (called from *rpart*)
            IF (lhetchem) THEN                                          !MSHET2
!!$               CALL het_new(jl,jk,krow,pw,tw,xxd,sadw,xnp,gasch,hiso)
               CALL het_new(jl,jk,krow,pw,tw,xxd,sadw,r_mean,avol,gasch,hiso)
               
               ! Subtract frozen part of HNO3 and H2O:
               gasch(11) = gasch(11) - gasch(77) ! psc I
               gasch(50) = gasch(50) - gasch(78) ! psc II
 
            END IF
 
            ! Determine concentrations of chemical species at next time step
            ! (model time step length sub-devided in nakop sub-time steps
            ! for chemistry):
            DO it = 1, nakop
 
               ! Determine rates of gas phase reactions and rename reaction
               ! rates of heterogeneous ("hiso") and photolysis ("diso")
               ! reactions:
               CALL fico(tw,xxd,gasch,diso,hiso)                          !MSv3
 
               ! Solve non-linear equation system by Newton-Raphson method:
               CALL newraph(jl,jk,krow,nuit,eeps,tmww)
 
               ! Solve non-linear equations of some long lived species
               ! (H2O, CH4, CO, H2, ODSs) by explicit Eulerian method (instead
               ! of Newton-Raphson):
               CALL cofi(gasch,tmww,clprodods,brprodods,noyloss,jk) !MSCLBI, !MSv3
 
            ENDDO
 
            ! Set negative chemical concentrations to zero:
            WHERE (gasch .LT. 0) gasch = 0.0_dp                           !MSv3
 
            ! Calculate ODSCLS, ODSCLL and ODSBR families (used for advection):
!!$            gasch(74) = 0.0_dp  ! ODSCLS
!!$            gasch(75) = 0.0_dp  ! ODSCLL
!!$            gasch(76) = 0.0_dp  ! ODSBR
!!$ 
!!$            gasch(74) = gasch(74)+gasch(55)*odsnat_chem( 1)  !MSODS
!!$            gasch(75) = gasch(75)+gasch(56)*odsnat_chem( 2)  !MSODS
!!$            gasch(76) = gasch(76)+gasch(57)*odsnat_chem( 3)  !MSODS
!!$            gasch(75) = gasch(75)+gasch(58)*odsnat_chem( 4)  !MSODS
!!$            gasch(75) = gasch(75)+gasch(59)*odsnat_chem( 5)  !MSODS
!!$            gasch(75) = gasch(75)+gasch(60)*odsnat_chem( 6)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(61)*odsnat_chem( 7)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(62)*odsnat_chem( 8)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(63)*odsnat_chem( 9)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(64)*odsnat_chem(10)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(65)*odsnat_chem(11)  !MSODS
!!$            gasch(76) = gasch(76)+gasch(66)*odsnat_chem(12)  !MSODS
!!$            gasch(76) = gasch(76)+gasch(67)*odsnat_chem(13)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(68)*odsnat_chem(14)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(69)*odsnat_chem(15)  !MSODS
!!$            gasch(74) = gasch(74)+gasch(70)*odsnat_chem(16)  !MSODS
!!$            gasch(76) = gasch(76)+gasch(71)*odsnat_chem(17)  !MSODS
!!$            gasch(76) = gasch(76)+gasch(72)*odsnat_chem(18)  !MSODS
!!$            gasch(76) = gasch(76)+gasch(73)*odsnat_chem(19)  !MSODS
!!$ 
!!$            IF (gasch94_b .GE. eps_mr) THEN
!!$               gasch94_cl = gasch94_cl*gasch(66)/gasch94_b
!!$            ENDIF
!!$ 
!!$            gasch(74) = gasch(74)+gasch94_cl*odsnat_chem(20) !H1211_CL
 
            ! Bry, Cly, NOy after chemistry:
            brsuma = &
            gasch(22)+ &
            gasch(24)+ &
            gasch(25)+ &
            gasch(26)+ &
            gasch(27)+ &
            gasch(28)+ &
            -1._dp*brprodods
           !minus bromine produced from ODSBR in the current time step
 
!!$            brsumaa = &
!!$            gasch(22)+ &
!!$            gasch(24)+ &
!!$            gasch(25)+ &
!!$            gasch(26)+ &
!!$            gasch(27)+ &
!!$            gasch(28)+ &
!!$            gasch(76)
!!$
            brsumaa = &
            gasch(22)+ &
            gasch(24)+ &
            gasch(25)+ &
            gasch(26)+ &
            gasch(27)+ &
            gasch(28)+ &
            gasch(57)+ &
            gasch(66)+ &
            gasch(67)+ &
            2._dp*gasch(71)+ &
            3._dp*gasch(72)+ &
            2._dp*gasch(73)
 
            clsuma = &
            gasch( 6)+ &
            gasch(13)+ &
            gasch(14)+ &
            gasch(18)+ &
            2._dp*gasch(19)+ &
            2._dp*gasch(20)+ &
            gasch(21)+ &
            gasch(28)+ &
            -1._dp*clprodods
           !minus chlorine produced from ODS in the current time step
!!$ 
!!$            clsumaa = &
!!$            gasch( 6)+ &
!!$            gasch(13)+ &
!!$            gasch(14)+ &
!!$            gasch(18)+ &
!!$            2._dp*gasch(19)+ &
!!$            2._dp*gasch(20)+ &
!!$            gasch(21)+ &
!!$            gasch(28)+ &
!!$            gasch(74) + gasch(75)

            clsumaa = &
            gasch( 6)+ &
            gasch(13)+ &
            gasch(14)+ &
            gasch(18)+ &
            2._dp*gasch(19)+ &
            2._dp*gasch(20)+ &
            gasch(21)+ &
            gasch(28)+ &
            3._dp*gasch(55)+ &
            2._dp*gasch(56)+ &
            3._dp*gasch(58)+ &
            2._dp*gasch(59)+ &
            gasch(60)+ &
            4._dp*gasch(61)+ &
            3._dp*gasch(62)+ &
            gasch(63)+ &
            2._dp*gasch(64)+ &
            gasch(65)+ &
            gasch(66)+ &
            gasch(68)+ &
            2._dp*gasch(69)+ &
            2._dp*gasch(70) 

            noysuma = &
            gasch( 4)+ &
            gasch( 7)+ &
            gasch( 9)+ &
            2._dp*gasch(10)+ &
            gasch(11)+ &
            gasch(12)+ &
            gasch(13)+ &
            gasch(15)+ &
            gasch(27)+ &
            gasch(37)+ &
            gasch(40)+ &
            gasch(44)+ &
            gasch(46)+ &
            2._dp*gasch(51)+ &
            gasch(77)+noyloss
 
            ! Check if Bry, Cly and NOy remain conserved by chemical
            ! calculations (Newton-Raphson method):
 
            ! Check for Bry (if mass is not conserved, set kluerrb = .TRUE.)
            IF(brsumb .LE. 1.0E4_dp) THEN
               IF (ABS(brsuma-brsumb).GT. 1.E4_dp) kluerrb = .TRUE.
            ELSE
               IF (ABS(brsuma-brsumb) .GT. 0.0001_dp*brsumb) kluerrb = .TRUE.
            ENDIF
            IF(brsumbb .LE. 1.0E4_dp) THEN
               IF (ABS(brsumaa-brsumbb) .GT. 1.E4_dp) kluerrb = .TRUE.
            ELSE
               IF (ABS(brsumaa-brsumbb) .GT. 0.0001*brsumbb) kluerrb = .TRUE.
            ENDIF
 
             ! Check for Cly (if mass is not conserved, set kluerrc = .TRUE.):
            IF(clsumb .LE. 1.0E4_dp) THEN
               IF(ABS(clsuma-clsumb) .GT. 1.E4_dp) kluerrc = .TRUE.
            ELSE
               IF(ABS(clsuma-clsumb) .GT. 0.0001_dp*clsumb) kluerrc = .TRUE.
            ENDIF
            IF(clsumbb .LE. 1.0E4_dp) THEN
               IF(ABS(clsumaa-clsumbb) .GT. 1.E4_dp) kluerrc = .TRUE.
            ELSE
               IF(ABS(clsumaa-clsumbb) .GT. 0.0001_dp*clsumbb) &
                    kluerrc = .TRUE.
            ENDIF
 
            ! Check for NOy (if mass is not conserved, set kluerrn = .TRUE.):
            IF(noysumb .LE. 1.0E4_dp) THEN                                !MSv3
               IF (ABS(noysuma-noysumb) .GT. 1.E4_dp) kluerrn = .TRUE.
            ELSE
               IF (ABS(noysuma-noysumb) .GT. 0.0001_dp*noysumb) kluerrn = .TRUE.
            ENDIF
 
            IF (nakop .LE. nakopmax) THEN                                 !MSv3
 
               ! Checks if nakop <= nakopmax:
 
               IF (kluerrc .OR. kluerrb .OR. kluerrn) THEN
                  ! If Bry, Cly and/or NOy is not conserved after chemical
                  ! calculations, pass this loop another time with halved time
                  ! step:
                  nakop = nakop*2
               ELSE
                  ! Otherwise exit the loop:
                  lexitloop = .TRUE.
               END IF
 
            ELSE
 
               ! Checks if nakop > nakopmax:
               ! If Bry, Cly and/or NOy are still not conserved, use values
               ! of time step before and print a warning.
               ! Exit loop.
 
               ! Bry:
               IF (kluerrb) THEN
                  gasch(22) = xt(jl,jk,idt_BR      ) *xxd
                  gasch(24) = xt(jl,jk,idt_BRO     ) *xxd
                  gasch(25) = xt(jl,jk,idt_HBR     ) *xxd
                  gasch(26) = xt(jl,jk,idt_HOBR    ) *xxd
                  gasch(27) = xt(jl,jk,idt_BRNO3   ) *xxd
                  gasch(28) = xt(jl,jk,idt_BRCL    ) *xxd
!!$                  gasch(76) = xt(jl,jk,idt_odsbr   ) *xxd
                  gasch(57) = xt(jl,jk,idt_cbrf3   ) *xxd
                  gasch(66) = xt(jl,jk,idt_h1211   ) *xxd
                  gasch(67) = xt(jl,jk,idt_ch3br   ) *xxd
                  gasch(71) = xt(jl,jk,idt_h2402   ) *xxd
                  gasch(72) = xt(jl,jk,idt_chbr3   ) *xxd
                  gasch(73) = xt(jl,jk,idt_ch2br2  ) *xxd
 
                  WRITE(*,*)' Problems reported in Bry', jl, krow, jk
               ENDIF
 
               ! Cly:
               IF (kluerrc) THEN
                  gasch( 6) = xt(jl,jk,idt_CLO     ) *xxd
                  gasch(13) = xt(jl,jk,idt_CLNO3   ) *xxd
                  gasch(14) = xt(jl,jk,idt_CL      ) *xxd
                  gasch(18) = xt(jl,jk,idt_HOCL    ) *xxd
                  gasch(19) = xt(jl,jk,idt_CL2     ) *xxd
                  gasch(20) = xt(jl,jk,idt_CL2O2   ) *xxd
                  gasch(21) = xt(jl,jk,idt_HCL     ) *xxd
                  gasch(28) = xt(jl,jk,idt_BRCL    ) *xxd
!!$                  gasch(74) = xt(jl,jk,idt_odscls  ) *xxd
!!$                  gasch(75) = xt(jl,jk,idt_odscll  ) *xxd
                  gasch(55) = xt(jl,jk,idt_f11     ) *xxd
                  gasch(56) = xt(jl,jk,idt_f12     ) *xxd
                  gasch(58) = xt(jl,jk,idt_cfc113  ) *xxd
                  gasch(59) = xt(jl,jk,idt_cfc114  ) *xxd
                  gasch(60) = xt(jl,jk,idt_cfc115  ) *xxd
                  gasch(61) = xt(jl,jk,idt_ccl4    ) *xxd
                  gasch(62) = xt(jl,jk,idt_ch3ccl3 ) *xxd
                  gasch(63) = xt(jl,jk,idt_hcfc22  ) *xxd
                  gasch(64) = xt(jl,jk,idt_hcfc141b) *xxd
                  gasch(65) = xt(jl,jk,idt_hcfc142b) *xxd
                  gasch(66) = xt(jl,jk,idt_h1211   ) *xxd
                  gasch(68) = xt(jl,jk,idt_ch3cl   ) *xxd
                  gasch(69) = xt(jl,jk,idt_hcfc21  ) *xxd
                  gasch(70) = xt(jl,jk,idt_hcfc123 ) *xxd

                  WRITE(*,*)' Problems reported in Cly', jl, krow, jk
               ENDIF
 
               ! NOy:
               IF(kluerrn) THEN
                  gasch( 4) = xt(jl,jk,idt_NO      ) *xxd
                  gasch( 7) = xt(jl,jk,idt_NO2     ) *xxd
                  gasch( 9) = xt(jl,jk,idt_NO3     ) *xxd
                  gasch(10) = xt(jl,jk,idt_N2O5    ) *xxd
                  gasch(11) = xt(jl,jk,idt_HNO3    ) *xxd
                  gasch(12) = xt(jl,jk,idt_HNO4    ) *xxd
                  gasch(13) = xt(jl,jk,idt_CLNO3   ) *xxd
                  gasch(15) = xt(jl,jk,idt_N       ) *xxd
                  gasch(27) = xt(jl,jk,idt_BRNO3   ) *xxd
                  gasch(37) = xt(jl,jk,idt_ISON    ) *xxd
                  gasch(40) = xt(jl,jk,idt_MPAN    ) *xxd
                  gasch(44) = xt(jl,jk,idt_NALD    ) *xxd
                  gasch(46) = xt(jl,jk,idt_PAN     ) *xxd
                  gasch(51) = xt(jl,jk,idt_N2O     ) *xxd
                  gasch(77) = xt(jl,jk,idt_psc1    ) *xxd
 
                  WRITE(*,*)' Problems reported in NOy', jl, krow ,jk
               ENDIF
 
               ! Exit loop:
               lexitloop = .TRUE.
 
            END IF
 
            ! Exit endless loop:
            IF (lexitloop) EXIT
 
         END DO
 
         ! Add frozen part of HNO3 and H2O:
         gasch(11) = gasch(11) + gasch(77)
         gasch(50) = gasch(50) + gasch(78)
 
         ! Tranform chemical species from mol/cm**3 to mixing ratio:      !MSv3
         xt(jl,jk,idt_O3      ) = gasch( 1)/xxd                    ! mol/mol
         xt(jl,jk,idt_OD      ) = gasch( 2)/xxd                    ! mol/mol
         xt(jl,jk,idt_O       ) = gasch( 3)/xxd                    ! mol/mol
         xt(jl,jk,idt_NO      ) = gasch( 4)/xxd                    ! mol/mol
         xt(jl,jk,idt_HO2     ) = gasch( 5)/xxd                    ! mol/mol
         xt(jl,jk,idt_CLO     ) = gasch( 6)/xxd                    ! mol/mol
         xt(jl,jk,idt_NO2     ) = gasch( 7)/xxd                    ! mol/mol
         xt(jl,jk,idt_OH      ) = gasch( 8)/xxd                    ! mol/mol
         xt(jl,jk,idt_NO3     ) = gasch( 9)/xxd                    ! mol/mol
         xt(jl,jk,idt_N2O5    ) = gasch(10)/xxd                    ! mol/mol
         xt(jl,jk,idt_HNO3    ) = gasch(11)/xxd                    ! mol/mol
         xt(jl,jk,idt_HNO4    ) = gasch(12)/xxd                    ! mol/mol
         xt(jl,jk,idt_CLNO3   ) = gasch(13)/xxd                    ! mol/mol
         xt(jl,jk,idt_CL      ) = gasch(14)/xxd                    ! mol/mol
         xt(jl,jk,idt_N       ) = gasch(15)/xxd                    ! mol/mol
         xt(jl,jk,idt_H2O2    ) = gasch(16)/xxd                    ! mol/mol
         xt(jl,jk,idt_H       ) = gasch(17)/xxd                    ! mol/mol
         xt(jl,jk,idt_HOCL    ) = gasch(18)/xxd                    ! mol/mol
         xt(jl,jk,idt_CL2     ) = gasch(19)/xxd                    ! mol/mol
         xt(jl,jk,idt_CL2O2   ) = gasch(20)/xxd                    ! mol/mol
         xt(jl,jk,idt_HCL     ) = gasch(21)/xxd                    ! mol/mol
         xt(jl,jk,idt_BR      ) = gasch(22)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH2O    ) = gasch(23)/xxd                    ! mol/mol
         xt(jl,jk,idt_BRO     ) = gasch(24)/xxd                    ! mol/mol
         xt(jl,jk,idt_HBR     ) = gasch(25)/xxd                    ! mol/mol
         xt(jl,jk,idt_HOBR    ) = gasch(26)/xxd                    ! mol/mol
         xt(jl,jk,idt_BRNO3   ) = gasch(27)/xxd                    ! mol/mol
         xt(jl,jk,idt_BRCL    ) = gasch(28)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH3     ) = gasch(29)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH3O2   ) = gasch(30)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH3O    ) = gasch(31)/xxd                    ! mol/mol
         xt(jl,jk,idt_HCO     ) = gasch(32)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH3O2H  ) = gasch(33)/xxd                    ! mol/mol
         xt(jl,jk,idt_C5H8    ) = gasch(34)/xxd                    ! mol/mol
         xt(jl,jk,idt_ISO2    ) = gasch(35)/xxd                    ! mol/mol
         xt(jl,jk,idt_ISO2H   ) = gasch(36)/xxd                    ! mol/mol
         xt(jl,jk,idt_ISON    ) = gasch(37)/xxd                    ! mol/mol
         xt(jl,jk,idt_MACR    ) = gasch(38)/xxd                    ! mol/mol
         xt(jl,jk,idt_MACRO2  ) = gasch(39)/xxd                    ! mol/mol
         xt(jl,jk,idt_MPAN    ) = gasch(40)/xxd                    ! mol/mol
         xt(jl,jk,idt_MACRO2H ) = gasch(41)/xxd                    ! mol/mol
         xt(jl,jk,idt_HACET   ) = gasch(42)/xxd                    ! mol/mol
         xt(jl,jk,idt_MGLY    ) = gasch(43)/xxd                    ! mol/mol
         xt(jl,jk,idt_NALD    ) = gasch(44)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH3CO3  ) = gasch(45)/xxd                    ! mol/mol
         xt(jl,jk,idt_PAN     ) = gasch(46)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH3CO3H ) = gasch(47)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH3COOH ) = gasch(48)/xxd                    ! mol/mol
         xt(jl,jk,idt_HCOOH   ) = gasch(49)/xxd                    ! mol/mol
         xt(jl,jk,idt_H2O     ) = gasch(50)/xxd                    ! mol/mol
         xt(jl,jk,idt_N2O     ) = gasch(51)/xxd                    ! mol/mol
         xt(jl,jk,idt_CH4     ) = gasch(52)/xxd                    ! mol/mol
         xt(jl,jk,idt_CO      ) = gasch(53)/xxd                    ! mol/mol
         xt(jl,jk,idt_H2      ) = gasch(54)/xxd                    ! mol/mol
!!$         xt(jl,jk,idt_odscls  ) = gasch(74)/xxd                    ! mol/mol
!!$         xt(jl,jk,idt_odscll  ) = gasch(75)/xxd                    ! mol/mol
!!$         xt(jl,jk,idt_odsbr   ) = gasch(76)/xxd                    ! mol/mol
         xt(jl,jk,idt_psc1    ) = gasch(77)/xxd                    ! mol/mol
         xt(jl,jk,idt_psc2    ) = gasch(78)/xxd                    ! mol/mol
         ! ODS species
         xt(jl,jk,idt_f11      ) = gasch(55)/xxd   !  mol/mol
         xt(jl,jk,idt_f12      ) = gasch(56)/xxd   !  mol/mol
         xt(jl,jk,idt_cbrf3    ) = gasch(57)/xxd   !  mol/mol
         xt(jl,jk,idt_cfc113   ) = gasch(58)/xxd   !  mol/mol
         xt(jl,jk,idt_cfc114   ) = gasch(59)/xxd   !  mol/mol
         xt(jl,jk,idt_cfc115   ) = gasch(60)/xxd   !  mol/mol
         xt(jl,jk,idt_ccl4     ) = gasch(61)/xxd   !  mol/mol
         xt(jl,jk,idt_ch3ccl3  ) = gasch(62)/xxd   !  mol/mol
         xt(jl,jk,idt_hcfc22   ) = gasch(63)/xxd   !  mol/mol
         xt(jl,jk,idt_hcfc141b ) = gasch(64)/xxd   !  mol/mol
         xt(jl,jk,idt_hcfc142b ) = gasch(65)/xxd   !  mol/mol
         xt(jl,jk,idt_h1211    ) = gasch(66)/xxd   !  mol/mol 
         xt(jl,jk,idt_ch3br    ) = gasch(67)/xxd   !  mol/mol
         xt(jl,jk,idt_ch3cl    ) = gasch(68)/xxd   !  mol/mol
         xt(jl,jk,idt_hcfc21   ) = gasch(69)/xxd   !  mol/mol
         xt(jl,jk,idt_hcfc123  ) = gasch(70)/xxd   !  mol/mol
         xt(jl,jk,idt_h2402    ) = gasch(71)/xxd   !  mol/mol
         xt(jl,jk,idt_chbr3    ) = gasch(72)/xxd   !  mol/mol
         xt(jl,jk,idt_ch2br2   ) = gasch(73)/xxd   !  mol/mol
 
       ! eth_as_ch4+
        ! chemical methane sinks:
        ch4descl_d(jl,jk,krow)   = -1._dp*gasch(115)  
        ch4desod_1_d(jl,jk,krow) = -1._dp*gasch(116)  
        ch4desod_2_d(jl,jk,krow) = -1._dp*gasch(117)
        ch4desoh_d(jl,jk,krow)   = -1._dp*gasch(118)
        ch4dest01_d(jl,jk,krow)  = -1._dp*gasch(119)
        ! eth_as_ch4-
        
        ! eth_as_ch4+
        IF (lch4_isotope) THEN
           xt(jl,jk,idt_12c)  = gasch(120) /xxd
           xt(jl,jk,idt_13c)  = gasch(121) /xxd
        END IF
        ! eth_as_ch4+ 

        ! eth_as_oh+
        ohprod_rn3_d(jl,jk,krow)   = gasch(122)
        ohdest_rn8_d(jl,jk,krow)   = -1._dp*gasch(123)   
        ohprod_rn11_d(jl,jk,krow)  = gasch(124)
        ohdest_rn12_d(jl,jk,krow)  = -1._dp*gasch(125)
        ohdest_rn13_d(jl,jk,krow)  = -1._dp*gasch(126)
        ohdest_rn21_d(jl,jk,krow)  = -1._dp*gasch(127)
        ohdest_rn24a_d(jl,jk,krow) = -1._dp*gasch(128)
        ohprod_rn27_d(jl,jk,krow)  = gasch(129)
        ohdest_rn31_d(jl,jk,krow)  = -1._dp*gasch(130)
        ohprod_rh1_d(jl,jk,krow)   =  2._dp*gasch(131)
        ohdest_rh2_d(jl,jk,krow)   = -1._dp*gasch(132)
        ohdest_rh3_d(jl,jk,krow)   = -1._dp*gasch(133)
        ohdest_rh4_d(jl,jk,krow)   = -1._dp*gasch(134)
        ohdest_rh13_d(jl,jk,krow)  = -2._dp*gasch(135)
        ohdest_rh11_d(jl,jk,krow)  = -1._dp*gasch(136)
        ohprod_rh5_d(jl,jk,krow)   = gasch(137)
        ohprod_rh6_d(jl,jk,krow)   = gasch(138)
        ohprod_rh8_d(jl,jk,krow)   =  2._dp*gasch(139)
        ohdest_rh12_d(jl,jk,krow)  = -1._dp*gasch(140)
        ohprod_rh12a_d(jl,jk,krow) = gasch(141)
        ohprod_rh9_d(jl,jk,krow)   = gasch(142)
        ohprod_rh14a_d(jl,jk,krow) =  2._dp*gasch(143)
        ohdest_rh15_d(jl,jk,krow)  = -1._dp*gasch(144)
        ohprod_rh16_d(jl,jk,krow)  = gasch(145)
        ohprod_rh18_d(jl,jk,krow)  = gasch(146)
        ohdest_rc5a_d(jl,jk,krow)  = -1._dp*gasch(147)
        ohdest_rc5b_d(jl,jk,krow)  = -1._dp*gasch(148)
        ohprod_rc9_d(jl,jk,krow)   = gasch(149)
        ohdest_rc10_d(jl,jk,krow)  = -1._dp*gasch(150)
        ohprod_rc13a_d(jl,jk,krow) = gasch(151)
        ohdest_rc6_d(jl,jk,krow)   = -1._dp*gasch(152)
        ohprod_rc26_d(jl,jk,krow)  = gasch(153)
        ohdest_rb14_d(jl,jk,krow)  = -1._dp*gasch(154)
        ohprod_rb15_d(jl,jk,krow)  = gasch(155)
        ohprod_rb16_d(jl,jk,krow)  = gasch(156)
        ohprod_rb17_d(jl,jk,krow)  = gasch(157)
        ohprod_rb22_d(jl,jk,krow)  = gasch(158)
        ohprod_rch1_d(jl,jk,krow)  = gasch(159)
        ohdest_rch3_d(jl,jk,krow)  = -1._dp*gasch(160)
        ohdest_rch7_d(jl,jk,krow)  = -1._dp*gasch(161)
        ohprod_rch7a_d(jl,jk,krow) = gasch(162)
        ohdest_rch14_d(jl,jk,krow) = -1._dp*gasch(163)
        ohprod_rch15_d(jl,jk,krow) = gasch(164)
        ohprod_rt06_d(jl,jk,krow)  = gasch(165)
        ohdest_rt07_d(jl,jk,krow)  = -1._dp*gasch(166)
        ohprod_rt08_d(jl,jk,krow)  = gasch(167)
        ohprod_rt11_d(jl,jk,krow)  = gasch(168)
        ohprod_rt12_d(jl,jk,krow)  = gasch(169)
        ohprod_rt17_d(jl,jk,krow)  = gasch(170)
        ohdest_rt18_d(jl,jk,krow)  = -1._dp*gasch(171)
        ohdest_rt19_d(jl,jk,krow)  = -1._dp*gasch(172)
        ohprod_rt20_d(jl,jk,krow)  = gasch(173)
        ohdest_rt26_d(jl,jk,krow)  = -2._dp*gasch(174)
        ohdest_ods10_d(jl,jk,krow) = -1._dp*gasch(175)
        ohdest_ods11_d(jl,jk,krow) = -1._dp*gasch(176)
        ohdest_ods13_d(jl,jk,krow) = -1._dp*gasch(177)
        ohdest_ods14_d(jl,jk,krow) = -1._dp*gasch(178)
        ohdest_ods15_d(jl,jk,krow) = -1._dp*gasch(179)
        ohdest_ods16_d(jl,jk,krow) = -1._dp*gasch(180)
        ohdest_ods18_d(jl,jk,krow) = -1._dp*gasch(181)
        ohdest_ods19_d(jl,jk,krow) = -1._dp*gasch(182)
  
        ho2dest_rn3_d(jl,jk,krow)   = -1._dp*gasch(183)
        ho2prod_rn20_d(jl,jk,krow)  = gasch(184)   
        ho2dest_rn22_d(jl,jk,krow)  = -1._dp*gasch(185)
        ho2prod_rn22a_d(jl,jk,krow) = gasch(186)
        ho2dest_rn27_d(jl,jk,krow)  = -1._dp*gasch(187)
        ho2prod_rh3_d(jl,jk,krow)   = gasch(188) 
        ho2dest_rh4_d(jl,jk,krow)   = -1._dp*gasch(189)
        ho2dest_rh5_d(jl,jk,krow)   = -1._dp*gasch(190)
        ho2dest_rh6_d(jl,jk,krow)   = -1._dp*gasch(191)
        ho2dest_rh7_d(jl,jk,krow)   = -2._dp*gasch(192)
        ho2dest_rh7a_d(jl,jk,krow)  = -2._dp*gasch(193)
        ho2prod_rh12_d(jl,jk,krow)  = gasch(194) 
        ho2prod_rh12a_d(jl,jk,krow) = gasch(195) 
        ho2prod_rh10_d(jl,jk,krow)  = gasch(196)
        ho2dest_rh14_d(jl,jk,krow)  = -1._dp*gasch(197)
        ho2dest_rh14a_d(jl,jk,krow) = -1._dp*gasch(198)
        ho2prod_rc5a_d(jl,jk,krow)  = gasch(199)
        ho2dest_rc11_d(jl,jk,krow)  = -1._dp*gasch(200)
        ho2dest_rc13_d(jl,jk,krow)  = -1._dp*gasch(201)
        ho2dest_rc13a_d(jl,jk,krow) = -1._dp*gasch(202)
        ho2prod_rc14_d(jl,jk,krow)  = gasch(203)
        ho2dest_rb02_d(jl,jk,krow)  = -1._dp*gasch(204)
        ho2prod_rb05_d(jl,jk,krow)  = gasch(205)
        ho2dest_rb13_d(jl,jk,krow)  = -1._dp*gasch(206)
        ho2prod_rch6_d(jl,jk,krow)  = gasch(207)
        ho2prod_rch10_d(jl,jk,krow) = gasch(208)
        ho2dest_rch13_d(jl,jk,krow) = -1._dp*gasch(209)
        ho2prod_rch21_d(jl,jk,krow) = gasch(210)
        ho2prod_rt07_d(jl,jk,krow)  = gasch(211)
        ho2dest_rt08_d(jl,jk,krow)  = -1._dp*gasch(212)
        ho2dest_rt09_d(jl,jk,krow)  = -1._dp*gasch(213)
        ho2prod_rt10_d(jl,jk,krow)  = gasch(214)
        ho2prod_rt19_d(jl,jk,krow)  = gasch(215)

        ! eth_as_oh-   

        !CCMI reactions, LR, 03/06/13
        rxn_1(jl,jk,krow) = gasch(216)
        rxn_2(jl,jk,krow) = gasch(217)
        rxn_3(jl,jk,krow) = gasch(218)
        rxn_4(jl,jk,krow) = gasch(219)
        rxn_5(jl,jk,krow) = gasch(220)

        rxn_6(jl,jk,krow) = gasch(221)
        rxn_7(jl,jk,krow) = gasch(222)
        rxn_8(jl,jk,krow) = gasch(223)
        rxn_9(jl,jk,krow) = gasch(224)
        rxn_10(jl,jk,krow) = gasch(225)

        rxn_11(jl,jk,krow) = gasch(226)
        rxn_12(jl,jk,krow) = gasch(227)
        rxn_13(jl,jk,krow) = gasch(228)
        rxn_14(jl,jk,krow) = gasch(229)
        rxn_15(jl,jk,krow) = gasch(230)

        rxn_16(jl,jk,krow) = gasch(231)
        rxn_17(jl,jk,krow) = gasch(232)
        rxn_18(jl,jk,krow) = gasch(233)
        rxn_19(jl,jk,krow) = gasch(234)
        rxn_20(jl,jk,krow) = gasch(235)
        rxn_21(jl,jk,krow) = gasch(236)
        rxn_22(jl,jk,krow) = gasch(237)
        rxn_23(jl,jk,krow) = gasch(238)
        rxn_24(jl,jk,krow) = gasch(239)
        rxn_25(jl,jk,krow) = gasch(240)

        rxn_26(jl,jk,krow) = gasch(241)
        rxn_27(jl,jk,krow) = gasch(242)
        rxn_28(jl,jk,krow) = gasch(243)
        rxn_29(jl,jk,krow) = gasch(244)
        rxn_30(jl,jk,krow) = gasch(245)

        rxn_31(jl,jk,krow) = gasch(246)
        rxn_32(jl,jk,krow) = gasch(247)
        rxn_33(jl,jk,krow) = gasch(248)
        rxn_34(jl,jk,krow) = gasch(249)

        jo2(jl,jk,krow) = diso(1)
        jo3d(jl,jk,krow) = diso(3)
        jno2(jl,jk,krow) = diso(5)
        jcl2o2(jl,jk,krow) = diso(17)

!!$        jmgly(jl,jk,krow) = diso(37)
!!$        jpan(jl,jk,krow) = diso(39)
!!$        jpaa(jl,jk,krow) = diso(40)
        
        IF (lo3orig) then
!!$           o3losst(jl,jk,krow) = gasch(79)/xxd ! mol/mol/s
           o3losst(jl,jk,krow) = 0._dp
           o3prod(jl,jk,krow) = gasch(80)/xxd  ! mol/mol/s
        endif

      ENDDO
 
   ENDDO
 
 END SUBROUTINE chem

 

 MODULE mo_socol_hetero
 
   ! Description:
 
   ! This module contains the subroutines to calculate the (first order) rate
   ! coefficients of the heterogeneous chemistry.
 
   ! Original code: Eugene Rozanov, PMOD/WRC, Davos,
   !                Martin Schraner, ETH Z?rich,
   !                Christopher Hoyle, PMOD/WRC, Davos
   ! f90 version:   Martin Schraner, ETH Z?rich, May 2009
 
   USE mo_constants,               ONLY: api, argas, g, amw
   USE mo_geoloc,                  ONLY: philon_2d, philat_2d
   USE mo_kind,                    ONLY: dp
   USE mo_socol_constants,         ONLY: atomweight, amh2so4, amhno3, amhcl, &
                                         amhocl, amhbr, amhobr, amnat, &
                                         sigmanat, rhonat, sigmaice, rhoice, &
                                         sigmaliq, psurf_stand
   USE mo_socol_dimensions,        ONLY: nhiso, ngasch, ngazov
   USE mo_socol_grid_calculations, ONLY: hetnat_lowlevind, hetice_lowlevind
   USE mo_socol_namelist,          ONLY: hetnat_uplev, hetnat_north, &
                                         hetnat_south, hetice_uplev, &
                                         hetice_north, hetice_south, &
                                         hetnat_rmode, hetice_ndens
!!$   USE mo_socol_streams,           ONLY: sadsts_d, sadpsc1_d, sadpsc2_d, &
!!$                                         sadrat_d, clacsts_d, clacpsc1_d, &
!!$                                         clacpsc2_d
 
   IMPLICIT NONE
 
   PRIVATE
 
   !--- Module variables:
   INTEGER, PUBLIC :: pscliqtest
   REAL(dp) :: t, aliq, rmean, ah2o, dl, hhocl, hhbr, hhobr, hstarhcl, &
        hstarhocl, hstarhbr, hstarhobr
 
   !--- Subroutines:
   PUBLIC :: het_new, calc_hetre
 
 CONTAINS
 
   SUBROUTINE het_new(jl,jk,krow,pph,temp,denh,sadw,r_mean,avol,gasch,hiso)
 
     ! This subroutine calls the different subroutines to calculate the
     ! (first order) rate coefficients of the heterogeneous chemistry (hiso),
     ! to re-calculate the concentrations of H2O(g), H2O(s), HNO3(g), and
     ! HNO3(s) (if PSC I and/or PSC II are available) and to determine the type
     ! of liquid aerosols (pscliqtest):
     ! The subroutines are called in the following sequence:
     ! 1. *hetice*   (calculate het. reaction rates on ice)
     ! 2. *hetnat*   (calculate het. reaction rates on NAT)
     ! 3. *sadcalc*  (calculate sulphate concentration of liquid aerosols)
     ! 4. *lacomp*   (calculate composition and Henry constant of liq. aerosols)
     ! 5. *hetla*    (calcualte het. reaction rates on liquid aerosols)
     ! Note: hiso(2), hiso(9), hiso(10), and hiso(11) are NOT calculated
     !       by *het_new*, but are determined later in subroutine
     !       *calc_hetre*, which is also part of this module.
 
     ! *het_new* is called from *chem*, src/socol_chemini.f90.
 
     ! Subroutine arguments:
     INTEGER, INTENT(in)     :: jl, jk, krow
     REAL(dp), INTENT(in)    :: pph, temp, denh, sadw, r_mean, avol
     REAL(dp), INTENT(inout) :: gasch(ngasch), hiso(nhiso)
 
     ! Local variables:
     INTEGER, PARAMETER :: iihobr = 1, iihocl = 2, iih2o = 3, iihno3 = 4
     REAL(dp) :: ppa, zh2o, zhno3, zhcl, zhbr, zhocl, zhobr, mh2so4, mhno3, &
          mhcl, mhocl, mhbr, mhobr, ns, ms, mn, ws, wn, wcl, wbr, whocl, &
          whobr, hhcl, rho, vliq, hno3s, h2os, algsigliq, &
          algsigliq2, patm, &
          ak150s, ak151s, ak152s, ak153s, ak154s, ak155s, ak156s, &
          ak157s, ak158s, ak159s, ak160s, ak161s, ak162s, &
          ak150p, ak151p, ak152p, ak153p, ak154p, ak155p, ak156p, &
          ak157p, ak158p, ak159p, ak160p, ak161p, ak162p
     REAL(dp) :: conc(ngazov)
     LOGICAL :: lpsc
 
 
     ! Executable statements:
 
     ! Save temperature and pressure locally:
     t      = temp            ! [K]
     ppa    = pph *100.0_dp   ! [hPa] -> [Pa]
 
     ! Save some chemical species as [mol/mol]:
     zhno3  = gasch(11)/denh ! [molec/cm3] -> [mol/mol]
     zhocl  = gasch(18)/denh ! [molec/cm3] -> [mol/mol]
     zhcl   = gasch(21)/denh ! [molec/cm3] -> [mol/mol]
     zhbr   = gasch(25)/denh ! [molec/cm3] -> [mol/mol]
     zhobr  = gasch(26)/denh ! [molec/cm3] -> [mol/mol]
     zh2o   = gasch(50)/denh ! [molec/cm3] -> [mol/mol]
 
     ! Save some chemical species as [molec/cm3]:
     conc(1) = zhobr*denh     ! HOBr
     conc(2) = zhocl*denh     ! HOCl
     conc(4) = zhno3*denh     ! HNO3
     conc(3) = zh2o *denh     ! H2O
 
     ! Conversion of molecular weights [g/mol] -> [kg/mol]:
     mh2so4 = amh2so4*0.001_dp
     mhno3  = amhno3 *0.001_dp
     mhcl   = amhcl  *0.001_dp
     mhocl  = amhocl *0.001_dp
     mhbr   = amhbr  *0.001_dp
     mhobr  = amhobr *0.001_dp
 
     ! Initialize some variables:
     hno3s      = 0.0_dp      ! HNO3 in the solid phase
     h2os       = 0.0_dp      ! H2O in the solid phase
     hhbr       = 0.0_dp
     hhocl      = 0.0_dp
     hhobr      = 0.0_dp
     hstarhcl   = 0.0_dp
     hstarhocl  = 0.0_dp
     hstarhbr   = 0.0_dp
     hstarhobr  = 0.0_dp
     aliq       = 0.0_dp
     rmean      = 0.0_dp
     ah2o       = 0.0_dp
     dl         = 0.0_dp
     pscliqtest = 0           ! 0 = no liq. aerosols, 1 = binary, 2 = ternary
                              ! liq. aerosols
 
!!$     ! Initialize output streams:
!!$     sadsts_d(jl,jk,krow)   = 0.0_dp
!!$     sadpsc1_d(jl,jk,krow)  = 0.0_dp
!!$     sadpsc2_d(jl,jk,krow)  = 0.0_dp
!!$     sadrat_d(jl,jk,krow)   = 1.0_dp
!!$     clacsts_d(jl,jk,krow)  = 0.0_dp
!!$     clacpsc1_d(jl,jk,krow) = 0.0_dp
!!$     clacpsc2_d(jl,jk,krow) = 0.0_dp
 
     ! Calculate heterogeneous chemistry on ice (if necessary):
     IF (jk .LE. hetice_lowlevind(jl) .AND. pph .GE. hetice_uplev .AND. &
          (philat_2d(jl,krow) .GE. hetice_north .OR. &
          philat_2d(jl,krow) .LE. hetice_south)) THEN
 
        ! Determine lpsc and calculate first order heterogeneous rate
        ! coefficients on ice (if lpsc=.TRUE.):
        CALL hetice
 
        ! Update corresponding values of hiso:
        IF (lpsc) THEN    ! ice available
           hiso(6)  = ak150p
           hiso(7)  = ak151p
           hiso(8)  = ak161p
           hiso(14) = ak152p
           hiso(15) = ak154p
           hiso(16) = ak156p                                             !MSHET2
 
           ! Update H2O(s) and H2O(g):
           gasch(78)= h2os              ! H2O(s) [molec/cm3]
           zh2o     = conc(iih2o)/denh  ! H2O(g) [mol/mol]               !MSHET
        ENDIF
 
     ENDIF
 
     ! Calculate heterogeneous chemistry on ice (if necessary):
     IF (jk .LE. hetnat_lowlevind(jl) .AND. pph .GE. hetnat_uplev .AND. &
          (philat_2d(jl,krow) .GE. hetnat_north .OR. &
          philat_2d(jl,krow) .LE. hetnat_south)) THEN
 
        ! Determine lpsc and calculate first order heterogeneous rate
        ! coefficients on NAT (if lpsc=.TRUE.):
        CALL hetnat
 
        IF (lpsc) THEN    ! NAT available
           ! Update corresponding values of hiso:
           hiso(4)  = ak150p
           hiso(5)  = ak151p
           hiso(12) = ak152p
           hiso(13) = ak154p
 
           ! Update HNO3(s) and HNO3(g):
           gasch(77)= hno3s             ! HNO3(s) [molec/cm3]
           zhno3    = conc(iihno3)/denh ! HNO3(g) [mol/mol]              !MSHET
        ENDIF
 
     ENDIF
 
     ! Calculate heterogeneous chemistry on liquid aerosols (if necessary;
     ! use HNO3(g),H2O(g) calcualted by *hetnat* and *hetice* as input):
     IF (sadw .GT. EPSILON(1._dp) .AND. r_mean .GT. EPSILON(1._dp)) THEN
 
        ! Calculate sulphate concentration ns of liquid aerosols from sadw
        ! (surface area density) and xnp (number density):               !MSHET
        CALL sadcalc                                                     !MSHET
 
        ! Determine psliqtest and calculate composition and Henry constants of
        ! liquid aerosol using ns (if pscliqtest .GE. 1):                !MSHET
        CALL lacomp                                                      !MSHET
 
        IF (pscliqtest .GE. 1) THEN  ! binary or ternary liquid aerosols !MSHET
           ! Use output of *lacomp* to calculate first order heterogeneous
           ! reaction coefficients on binary and ternary liquid aerosols:
           CALL hetla                                           !MSHET, !MSHET2
 
           ! Update corresponding values of hiso:
           hiso(1) = ak150s
           hiso(3) = ak161s
        ENDIF
 
     ENDIF
 
   CONTAINS
 
     SUBROUTINE hetice
 
       ! This subroutine calculates the 1st order heterogeneous reaction
       ! coefficients on ice of the following reactions:
       !
       !      N2O5   + H2O -> 2HNO3             (ak150p/hiso(6))
       !      ClONO2 + H2O -> HOCl     + HNO3   (ak151p/hiso(7))
       !      ClONO2 + HCl -> Cl2      + HNO3   (ak152p/hiso(14))
       !      N2O5   + HCl -> Cl + NO2 + HNO3   (ak153p)
       !      HOCl   + HCl -> Cl2      + H2O    (ak154p/hiso(15))
       !      HOBr   + HBr -> 2Br      + H2O    (ak155p)
       !      HOBr   + HCl -> BrCl     + H2O    (ak156p/hiso(16))
       !      HOCl   + HBr -> BrCl     + H2O    (ak157p)
       !      BrONO2 + HBr -> 2Br      + HNO3   (ak158p)
       !      BrONO2 + HCl -> BrCl     + HNO3   (ak159p)
       !      ClONO2 + HBr -> BrCl     + HNO3   (ak160p)
       !      BrONO2 + H2O -> HOBr     + HNO3   (ak161p/hiso(8))
       !      N2O5   + HBr -> Br + NO2 + HNO3   (ak162p)
 
       ! Local variables:
       REAL(dp) :: algsigice, algsigice2, freqc, freqn, freqh, freqb, freqo, &
            rice, rice2, rs4l, vice, rmodeice, zh2oeq, ztpsc, ptot
 
       ! Local parameters:
 
       ! Gamma values:
       REAL(dp), PARAMETER :: g150 = 0.01_dp, g151 = 0.1_dp, g152 = 0.2_dp, &
            g153 = 0.03_dp, g154 = 0.3_dp, g155 = 0.1_dp, g156 = 0.3_dp, &
            g157 = 0.3_dp,  g158 = 0.0_dp, g159 = 0.0_dp, g160 = 0.0_dp, &
            g161 = 0.3_dp,  g162 = 0.0_dp
 
 
       ! Executable statements:
 
       ztpsc = t
       ztpsc = MIN(ztpsc, 273.0_dp)
 
       ! Equilibrium vapour pressure of water vapour [molec/cm^3]
       ! (*denh/ppa: [Pa] -> [molec/cm^3]):
       zh2oeq = 10.0_dp**(10.431_dp - 2668.70_dp/ztpsc)*psurf_stand/760._dp &
            *denh/ppa
 
       ! Test for ice particles:
       IF(conc(iih2o) .LE. zh2oeq) THEN
          lpsc = .FALSE.
          RETURN      ! Exit subroutine if there are no ice particles
       ENDIF
 
       lpsc = .TRUE.
 
       ! New concentrations of water vapour(g) and ice [molec/cm^3]:
       h2os = conc(iih2o) - zh2oeq
       conc(iih2o) = conc(iih2o) - h2os
 
       ! Algorithm of sigma of log normal distribution:
       algsigice = LOG(sigmaice)
       algsigice2 = algsigice*algsigice                                  !MSHET
 
       ! Pressure in hPa:
       ptot = 0.01_dp*ppa
 
       ! Total ice volume (dimensionless e.g. cm3/cm3 air):
       vice = amw*atomweight*h2os*1000.0_dp/rhoice
 
       ! Mean radius of ice paticles rice [cm] (4.5=9./2.):
       rmodeice = (3.0_dp*vice/(hetice_ndens*4.0_dp*api* &
            EXP(4.5_dp*algsigice2)))**(1.0_dp/3.0_dp)                    !MSHET
       rice = rmodeice*EXP(algsigice2)                                   !MSHET
 
       ! Since ice particles can be very large, gas diffusion limitation
       ! is taken into account. The first order rate constants are
       ! calculated from the equation K=G*api*R**2*CBAR*N/(1+3G.R/(4.L)), where
       ! G=gamma, R=particle radius, CBAR=mean molecular speed, L=mean free
       ! path (Turco et al, JGR, 94, 16493, 1989).
       rice2 = rice*rice*hetice_ndens
       freqc = 4.56E4_dp*SQRT(t/ 97.45_dp)*rice2
       freqn = 4.56E4_dp*SQRT(t/108.0_dp )*rice2
       freqh = 4.56E4_dp*SQRT(t/ 52.5_dp )*rice2
       freqb = 4.56E4_dp*SQRT(t/142.0_dp )*rice2
       freqo = 4.56E4_dp*SQRT(t/ 97.0_dp )*rice2
       rs4l  = 3.3E4_dp*rice*ptot/t
 
       ! 1st order heterogeneous reaction coefficients:
       ak150p = g150*freqn/(1.0_dp + g150*rs4l)
       ak151p = g151*freqc/(1.0_dp + g151*rs4l)
       ak152p = g152*freqc/(1.0_dp + g152*rs4l)
       !ak153p = g153*freqn/(1.0_dp + g153*rs4l)       ! currently not used
       ak154p = g154*freqh/(1.0_dp + g154*rs4l)
       !ak155p = g155*freqo/(1.0_dp + g155*rs4l)       ! currently not used
       ak156p = g156*freqo/(1.0_dp + g156*rs4l)
       !ak157p = g157*freqh/(1.0_dp + g157*rs4l)       ! currently not used
       !ak158p = g158*freqb/(1.0_dp + g158*rs4l)       ! currently not used
       !ak159p = g159*freqb/(1.0_dp + g159*rs4l)       ! currently not used
       !ak160p = g160*freqc/(1.0_dp + g160*rs4l)       ! currently not used
       ak161p = g161*freqb/(1.0_dp + g161*rs4l)
       !ak162p = g162*freqn/(1.0_dp + g162*rs4l)       ! currently not used
 
!!$       ! Save some output:
!!$       sadpsc2_d(jl,jk,krow) = 4.E08_dp*api*rice2 ! SAD [um^2/cm^3]
!!$       clacpsc2_d(jl,jk,krow) = ak151p + ak152p   ! ClONO2 activation [s^-1]
 
     END SUBROUTINE hetice
 
 
     SUBROUTINE hetnat
 
       ! This subroutine calculates the 1st order heterogeneous reaction
       ! coefficients on ice of the following reactions:
       !
       !      N2O5   + H2O -> 2HNO3             (ak150p/hiso(4))
       !      ClONO2 + H2O -> HOCl     + HNO3   (ak151p/hiso(5))
       !      ClONO2 + HCl -> Cl2      + HNO3   (ak152p/hiso(12))
       !      N2O5   + HCl -> Cl + NO2 + HNO3   (ak153p)
       !      HOCl   + HCl -> Cl2      + H2O    (ak154p/hiso(13))
       !      HOBr   + HBr -> 2Br      + H2O    (ak155p)
       !      HOBr   + HCl -> BrCl     + H2O    (ak156p)
       !      HOCl   + HBr -> BrCl     + H2O    (ak157p)
       !      BrONO2 + HBr -> 2Br      + HNO3   (ak158p)
       !      BrONO2 + HCl -> BrCl     + HNO3   (ak159p)
       !      ClONO2 + HBr -> BrCl     + HNO3   (ak160p)
       !      BrONO2 + H2O -> HOBr     + HNO3   (ak161p)
       !      N2O5   + HBr -> Br + NO2 + HNO3   (ak162p)
 
       ! Local variables:
       REAL(dp) :: algsignat, algsignat2, expalgsignat45, nnat, nnat_diff, &
            vnat_diff, hno3s_diff, rmodenat3, rnat, rnat2, vnat, weightrhonat, &
            sice, freqc, freqn, freqh, freqb, freqo, rs4l, zh2ot, zhno3eq, &
            zbt, zmt, ztpsc, delt, ptot, pw, tice, g151, g152
 
       ! Local parameters:
 
       ! Upper boundary for number density of NAT:
       REAL(dp), PARAMETER :: nnat0 = 5.0E-04_dp
 
       ! Gamma values:
       REAL(dp), PARAMETER :: g150 = 3.0E-04_dp, g153 = 0.003_dp, &
                 g154 = 0.1_dp, g155 = 0.1_dp, g156 = 0.1_dp, g157 = 0.1_dp,   &
                 g158 = 0.0_dp, g159 = 0.0_dp, g160 = 0.0_dp, g161 = 0.006_dp, &
                 g162 = 0.0_dp
 
 
       ! Executable statements:
 
       ztpsc = t
 
       ! Equilibrium vapour pressure of HNO3 [molec/cm^3]:
       zh2ot = 0.0075_dp*ppa*conc(iih2o)/denh                   ! H2O in torr
       zh2ot = MAX(zh2ot, 1.0E-12_dp)
       zmt = -2.7836_dp - 0.00088_dp*ztpsc
       zbt = 38.9855_dp - 11397.0_dp/ztpsc + 0.009179_dp*ztpsc
       ! (133.3_dp*denh/ppa: [Torr] -> [molec/cm^3]):
       zhno3eq = 10.0_dp**(zmt*LOG10(zh2ot) + zbt)*133.3_dp*denh/ppa
 
       ! Test for NAT particles:
       IF(conc(iihno3) .LE. zhno3eq) THEN
          lpsc = .FALSE.
          RETURN      ! Exit subroutine if there are no NAT particles
       ENDIF
 
       lpsc = .TRUE.
 
       algsignat = LOG(sigmanat)
       algsignat2 = algsignat*algsignat
       expalgsignat45 = EXP(4.5_dp*algsignat2)                           !MSHET
       rmodenat3 = hetnat_rmode*hetnat_rmode*hetnat_rmode
       weightrhonat = amnat*atomweight*1000.0_dp/rhonat
 
       ! Check if nnat0 or a smaller value has to be used as number density:
       hno3s_diff = conc(iihno3)-zhno3eq                                 !MSHET
       vnat_diff = hno3s_diff*weightrhonat                               !MSHET
       nnat_diff = (vnat_diff*3.0_dp)/(4.0_dp*api*rmodenat3*expalgsignat45)
       nnat = MIN(nnat0, nnat_diff)
 
       ! Total NAT volume  (dimensionless e.g. cm3/cm3 air):             !MSHET
       vnat = 4.0_dp/3.0_dp*api*rmodenat3*nnat*expalgsignat45            !MSHET
 
       ! Calculate corresponding mixing ratio of NAT [mol/cm3]:          !MSHET
       hno3s = vnat/weightrhonat
 
       ! Mean radius of NAT:
       rnat = hetnat_rmode*EXP(algsignat2)                               !MSHET
 
       ! New concentrations of HNO3(g) [molec/cm^3]:
       conc(iihno3) = conc(iihno3) - hno3s
 
       ! Pressure in hPa:
       ptot = 0.01_dp*ppa
 
       ! Calculate gamma of reactions ClONO2 + H2O -> HOCl + HNO3 (g151) and
       ! ClONO2 + HCl -> 2Cl + HNO3 (g152) from data of Hanson and Ravishankara.
       ! For ClONO2+HCl their fit is used;  for ClONO2+H2O, gamma has been
       ! fitted as a function of the saturation ratio with respect to ice:
 
       ! H2O partial pressure [atm]:
       pw = zh2o*ppa/psurf_stand
 
       ! Ice frost point temperature:
       tice = 2668.70_dp/(10.431_dp - 0.434294_dp*(LOG(pw*760.0_dp)))
       delt = t - tice
 
       ! Ice saturation:
       sice = pw/10**(10.431_dp - 2668.70_dp/t)*760.0_dp
       sice = MIN(sice,5.0_dp)
 
       ! Gamma values of reactions 151 and 152:
       g151 = EXP(-9.03_dp + 2.81_dp*sice)
       g152 = 1.0_dp/(4.348_dp + 1.0_dp/(0.7022_dp*EXP(-0.518_dp*delt)))
 
       ! Since NAT particles can be very large, gas diffusion limitation
       ! is taken into account. The first order rate constants are
       ! calculated from the equation K=G*api*R**2*CBAR*N/(1+3G.R/(4.L)), where
       ! G=gamma, R=particle radius, CBAR=mean molecular speed, L=mean free
       ! path (Turco et al, JGR, 94, 16493, 1989).
       rnat2 = rnat*rnat*nnat                                           !MSHET
       freqc = 4.56E4_dp*SQRT(t/ 97.45_dp)*rnat2
       freqn = 4.56E4_dp*SQRT(t/108.0_dp )*rnat2
       freqh = 4.56E4_dp*SQRT(t/ 52.5_dp )*rnat2
       freqb = 4.56E4_dp*SQRT(t/142.0_dp )*rnat2
       freqo = 4.56E4_dp*SQRT(t/ 97.0_dp )*rnat2
       rs4l  = 3.3E4_dp*rnat*ptot/t
 
       ! 1st order heterogeneous reaction coefficients:
       ak150p = g150*freqn/(1.0_dp + g150*rs4l)
       ak151p = g151*freqc/(1.0_dp + g151*rs4l)
       ak152p = g152*freqc/(1.0_dp + g152*rs4l)
       !ak153p = g153*freqn/(1.0_dp + g153*rs4l)       ! currently not used
       ak154p = g154*freqh/(1.0_dp + g154*rs4l)
       !ak155p = g155*freqo/(1.0_dp + g155*rs4l)       ! currently not used
       !ak156p = g156*freqo/(1.0_dp + g156*rs4l)       ! currently not used
       !ak157p = g157*freqh/(1.0_dp + g157*rs4l)       ! currently not used
       !ak158p = g158*freqb/(1.0_dp + g158*rs4l)       ! currently not used
       !ak159p = g159*freqb/(1.0_dp + g159*rs4l)       ! currently not used
       !ak160p = g160*freqc/(1.0_dp + g160*rs4l)       ! currently not used
       !ak161p = g161*freqb/(1.0_dp + g161*rs4l)       ! currently not used
       !ak162p = g162*freqn/(1.0_dp + g162*rs4l)       ! currently not used
 
!!$       ! Save some output:
!!$       sadpsc1_d(jl,jk,krow) = 4.E08_dp*api*rnat2 ! SAD [um^2/cm^3]
!!$       clacpsc1_d(jl,jk,krow) = ak151p + ak152p   ! ClONO2 activation [s^-1]
 
     END SUBROUTINE hetnat
 
 
     SUBROUTINE sadcalc                                                  !MSHET
 
       ! This subroutine calculates the amount of H2SO4 (ns) [mol/m3] from
       ! sulphate aerosol surface area density (sadw) [um/cm3] and number
       ! density of liquid aerosols (xnp) [1/cm3].
 
       ! Original code: Chris Hoyle, PMOD/WRC, Davos, November 2004
 
       REAL(dp) :: rmode, avolume             ! Mode radius, volume of aerosols
       REAL(dp), PARAMETER :: ws_est = 0.7_dp ! Assumed weight fraction of H2SO4
       REAL(dp), PARAMETER :: wn_est = 0.0_dp ! Assumed weight fraction of HNO3
 
 
       ! Executable statements:
 
       ! Logarithm of sigma of lognormal distribution:
       algsigliq = LOG(sigmaliq)
       algsigliq2 = algsigliq*algsigliq
 
       ! Calculate the mode radius of the aerosol droplets (rmode) [cm].
       ! Assuming a lognormal size distribution for the aerosol droplets, the
       ! surface area is given by
       ! sadw= xnp*(4.*pi*rmode**2)*exp**((4./2.)*(alog(sigmaliq))**2).
       ! Re-arranging for rmode (*1.E-04: conversion [um]->[cm]):
!!$       rmode = 1.E-04_dp*SQRT(sadw/(4._dp*api*xnp*EXP(2._dp*algsigliq2)))
 
       ! Calculate the volume of the aerosol:
!!$       avolume = 4._dp/3._dp*api*xnp*rmode*rmode*rmode*EXP(4.5_dp*algsigliq2)
       avolume = avol*1.e-12_dp  ![um3/cm3 -> cm3/cm3]

       ! Calculate the density rho for estimated H2SO4- and HNO3-weight-fractions
       ! ws_est and wn_est (SAGE2 measures background droplets (with no
       ! HNO3-uptake), so let's set ws_est=0.7, wn_est=0.0):             !MSHET
       CALL density_liq(ws_est,wn_est,t,rho)                             !MSHET
 
       ! Calculate the amount of H2SO4 (ns) [mol/m3] of the aerosol
       ! (*1.E06: [mol/cm3]->[mol/m3]):
       ns = 1.E06_dp*ws_est*avolume*rho/amh2so4  ![mol/cm3]              !MSHET
 
     END SUBROUTINE sadcalc
 
 
     SUBROUTINE lacomp
 
       ! This subroutine calculates the composition of aqueous
       ! HNO3/H2SO4/HCl/HOCl stratospheric aerosols
       ! (Carslaw, Luo, Peter - Geophys. Res. Lett., 1995).
       !
       ! HNO3/H2SO4 composition based on thermodynamic model of system
       ! HCl/HBr/HNO3/H2SO4/H2O (Carslaw et al, J. Phys. Chem., 1995 -
       ! A thermodynamic modelof the system HCl-HNO3-H2SO4-H2O,
       ! including solubilities of HBr, from <200K to 328K.)
       !
       ! HCl  solubility parameterisation from Luo et al., GRL, 1995
       ! HOCl solubility parameterisation from Huthwelker et al., JAS., 1995
       !
       ! The model is valid for 2E-5 mb < pw < 2E-3 mb
       ! (pw is the water partial pressure).
       !
       ! The upper temperature limit is 240 K
       ! The lower temperature limit is 185 K or TICE-3 K, (whichever is  higher)
       !
       ! HNO3 solubilities are calculated to a maximum temperature of 215K.
       ! Solubilities are calculated on a molality basis
       ! (moles of solute per kg of water).
       !
       ! The solubilities of HCl and HOCl are assumed not to affect
       ! the uptake of HNO3 or H2O.
       ! This introduces ONLY small errors at the very lowest temperatures,
       ! but for a full calculation WHERE interactions
       ! between all species in solution are considered, USE
       ! the model of Carslaw et al., given above.
       !
       ! This subroutine has been adapted from original code of Carslaw.
       !
       ! Arguments:
       ! 1)  t      Temperature [K]           (input)
       ! 2)  ppa    Pressure    [Pa]          (input)
       ! 3)  zh2o   H2O   volume mixing ratio (input)
       ! 4)  zhno3  HNO3  volume mixing ratio (input)
       ! 5)  ZHCL   HCl   volume mixing ratio (input)
       ! 6)  zhbr   HBr   volume mixing ratio (input)
       ! 7)  zhocl  HOCl  volume mixing ratio (input)
       ! 8)  zhobr  HOBr  volume mixing ratio (input)
       ! 9)  ns     H2SO4 concentration [mol/m3] (input)                !MSHET
       ! 10) ms
       ! 11) mn
       ! 12) ws     composition of aerosol (wt fraction H2SO4)
       ! 13) wn
       ! 14) wcl
       ! 15) wbr
       ! 16) whocl
       ! 17) whobr
       ! 18) rho    Liquid solution density [g/cm3]
       ! 19) vliq   specific volume of aerosol [cm^3 cm^-3]
       ! 20) hhcl   (mol/kg/atm) Henry's law coefficient
       ! 21) hhbr   (mol/kg/atm) Modified Henry's law coefficient
       ! 22) hhocl  (mol/kg/atm) Henry's law coefficient
       ! 23) hhobr  (mol/kg/atm) Henry's law coefficient
       ! 24) pscliqtest: 0=no liq. aerosols, 1=binary, 2=ternary liq. aerosols
       !
       ! hhocl, hhcl etc = Effective Henry's law constants in HNO3-H2SO4-H2O
       !                   solution [mol/kg water/atm]
       ! hnb, hsb = Effective henry's law constant of HNO3 in HNO3-H2O and
       !            H2SO4-H2O solutions [mol/kg water/atm]
       ! ms, mn, mcl, mhocl = Concentrations of H2SO4, HNO3, HCL and HOCL
       !                   [mol/kg of water]
       ! msb, mnb = Binary solution concentrations [mol/kg of water]
       ! ns = Total moles sulphate/m3 air (assumed to be totally in liquid
       !      phase)
       ! ptot = Total atmospheric pressure [mb]
       ! pn0, phcl0 etc = Initial HNO3, HCl etc partial pressure [atm]
       ! qhno3, qhcl etc = Total volume mixing ratios of hno3, hcl etc
       ! pw = Water partial pressure [atm]
       ! tice = Ice frost point temperature [K]
       ! vliq = Specific aerosol volume (dimensionless)
       ! ws, wn etc = Weight fractions of H2SO4, HNO3 etc in the aerosol
 
       ! Local variables:
       REAL(dp) :: tlocal, tt, phcl0, phocl0, phbr0, phobr0, ph2so4, &
            muh2so4, zh2so4, pn0, hsb, sqrtwn, sqrtws, msb, mnb, mcl, mbr, &
            mnb2, mnbmsb, wnws, ws2, wn2, pwlog, tlog, tice, pw, pr, tr, xsb, &
            hnb, xnb, ust, tr2, tr3, pr2, pr3, dg, a, b, c, d, e, phi, wts, &
            algm, awice, vice, vwater
 
       ! Local parameters:
 
       ! Gas constant/101325 Pa [m^3 atm mol^-1 K^-1]:
       REAL(dp), PARAMETER :: r = argas/psurf_stand
 
       ! Limits of H2O [atm]:
       REAL(dp), PARAMETER :: pwl=2.0E-8_dp, pwu=2.0E-6_dp
 
       REAL(dp), PARAMETER :: &
            qn1 = 14.5734_dp   , qn2 = 0.0615994_dp, qn3 = -1.14895_dp , &
            qn4 =  0.691693_dp , qn5 = -0.098863_dp, qn6 = 0.0051579_dp, &
            qn7 = 0.123472_dp  , qn8 = -0.115574_dp, qn9 = 0.0110113_dp, &
            qn10 = 0.0097914_dp, &
            qs1 = 14.4700_dp   , qs2 = 0.0638795_dp, qs3 = -3.29597_dp , &
            qs4 =  1.778224_dp , qs5 = -0.223244_dp, qs6 = 0.0086486_dp, &
            qs7 = 0.536695_dp  , qs8 = -0.335164_dp, qs9 = 0.0265153_dp, &
            qs10 = 0.0157550_dp, &
            kn1 = -39.136_dp   , kn2 =  6358.4_dp  , kn3 =  83.29_dp   , &
            kn4 = -17650.0_dp  , kn5 = 198.53_dp   , kn6 = -11948_dp   , &
            kn7 = -28.469_dp   , &
            ks1 = -21.661_dp   , ks2 =  2724.2_dp  , ks3 =  51.81_dp   , &
            ks4 = -15732.0_dp  , ks5 = 47.004_dp   , ks6 = -6969.0_dp  , &
            ks7 = -4.6183_dp, &
            alg36 = 3.59624_dp , alg80 = 4.39344_dp, alg760 = 6.63332_dp, &
            xpi = 3.14159_dp, t0 = 273.16_dp, vwater185 = 0.000279926814_dp, &
            us298 = 1.0_dp/298.150_dp
 
 
       ! Executable statements:
 
       ! Pressure [atm]:
       patm = ppa/psurf_stand
 
       ! Set tlocal (local temperature) to 185 K if t < 185 K (for stability of
       ! the parameterization) and adapt water vapour (for same water activity)
       ! (Beiping Luo, 7 Dec 2005):                                      !MSHET
       IF (t .GE. 185.0_dp) THEN
          tlocal = t
          ! H2O partial pressure [atm]:
          pw = MAX(EPSILON(1.E-12_dp), MIN(pwu, zh2o*patm))              !MSHET
       ELSE
          tlocal = 185.0_dp
          dg = -210368.0_dp - 131.438_dp*t + 3.32373E06_dp/t + 41729.1_dp*LOG(t)
          awice = EXP(-dg/(8.31441_dp*t))
          vice = 5.75185606E10_dp*EXP(-20.947031_dp*t0/t - &
               3.56654_dp*LOG(t0/t) - 2.01889049_dp/t0*t)
          vwater = vice/awice    ! saturation pressure of water vapour
          ! H2O partial pressure [atm]:                                  !MSHET
          pw = MAX(EPSILON(1.E-12_dp), MIN(pwu, zh2o*patm*vwater185/vwater))
       ENDIF
 
       ! Logarithm of water vapour partial pressure:
       pwlog = LOG(pw)
 
       ! Ice point temperature [K]:
       tice = 2668.70_dp/(10.4310_dp - 0.434294_dp*(pwlog + alg760))
 
       ! Determine pscliqtest (0=no liq aerosol, 1=binary, 2=ternary aerosol):
       IF (t .GE. (tice - 3.0_dp) .AND. t .LE. 240.0_dp .AND. &
            pw .GE. pwl) THEN  ! (where parameterisation is valid)
          ! New criterion for ternary aerosols (Beiping Luo, 2 Sep 2005):
          IF (t .GE. (tice + 6.0_dp)) THEN                          !MSHET
             pscliqtest = 1
          ELSE
             pscliqtest = 2
          ENDIF
       ELSE
          pscliqtest = 0
          RETURN  ! Exit subroutine
       ENDIF
 
       ! Initialize some variables with zero:
       ws    = 0.0_dp
       wn    = 0.0_dp
       wcl   = 0.0_dp
       wbr   = 0.0_dp
       whocl = 0.0_dp
       whobr = 0.0_dp
       ms    = 0.0_dp
       mn    = 0.0_dp
       vliq  = 0.0_dp
       hhcl  = 0.0_dp
 
       ! Logarithm and inverse of (local) temperature:
       tlog = LOG(tlocal)
       ust = 1.0_dp/tlocal
 
       ! Partial pressure of some chemical species [atm]:
       pn0    = zhno3*patm
       phcl0  = zhcl *patm
       phocl0 = zhocl*patm
       phbr0  = zhbr *patm
       phobr0 = zhobr*patm
 
       ! H2SO4 volume mixing ratio [mol/mol]
       ! (Chris Hoyle changed to calculate zh2so4 (16/11/04)):           !MSHET
       zh2so4 = argas*tlocal*ns/ppa
 
       ! Partial pressure of H2SO4 [atm]:
       tt = r*tlocal*ns
 
       pr = pwlog + 18.4_dp
       tr = 1.0E4_dp*ust - 43.4782608_dp
       tr2 = tr *tr
       tr3 = tr2*tr
       pr2 = pr *pr
       pr3 = pr2*pr
 
       ! I. Calculate the pure binary H2SO4/H2O solution composition:
 
       ! Preparations:
       a = ks1 + ks2*ust
       b = ks3 + ks4*ust
       xsb = (-a - SQRT(a*a - 4.0_dp*b*(ks5 + ks6*ust + ks7*tlog - pwlog)))/ &
            (2.0_dp*b)
       msb = 55.51_dp*xsb/(1.0_dp-xsb)
 
       ! Concentration of H2SO4 and HNO3 in the binary solution [mol/kg]:
       ms = msb
       mn = 0.0_dp
 
       ! Weight fraction of H2SO4 and HNO3 in binary solution [kg/kg]:
       ws = msb*mh2so4/(1.0_dp + msb*mh2so4)
       wn = 0.0_dp
 
       ! Calculate H2SO4 equilibrium vapour pressure (ph2so4) [atm] from the
       ! expression for the H2SO4 vapour pressure over H2SO4/H2O solutions
       ! from paper of Ayers et al., GRL, 7, 433-436, 1980.
       ! The Giaugue expression for chemical potential has been fitted
       ! to a simple expression (see the paper). It is likely to be a bit
       ! inaccurate for use at low strat conditions, but should be reasonable:
 
       ! Preparations:
       wts = MAX(41.0_dp, ws*100.0_dp)       ! Weight fraction of H2SO4 [%]
       muh2so4 = 4.184_dp*(1.514E4_dp - 286.0_dp*(wts - 40.0_dp) + 1.080_dp* &
            (wts - 40.0_dp)**2.0_dp - 3941.0_dp/(wts - 40.0_dp)**0.1_dp)
 
       ! H2SO4 equilibrium vapour pressure (ph2so4) [atm]:
       ph2so4 = EXP(-10156.0_dp*ust + 16.2590_dp - muh2so4/(8.314_dp*tlocal))
 
       ! Switch off the aerosols if the partial pressure of H2O4 < saturation
       ! pressure of H2SO4:
       IF (zh2so4*ppa/psurf_stand .LT. ph2so4) THEN
          pscliqtest = 0
          RETURN   ! Exit subroutine
       ENDIF
 
       ! If the uptake of HNO3 is possible:
       IF (pscliqtest .EQ. 2) THEN
          ! II. Calculate the ternary HNO3/H2SO4/H2O solution composition:
 
          ! Effective Henry constant of H2SO4 [mol kg^-1 atm^-1]:
          hsb = EXP(qs1 + qs2*tr2 + (qs3 + qs4*tr + qs5*tr2 + qs6*tr3)*pr + &
               (qs7 + qs8*tr + qs9*tr2)*pr2 + qs10*tr*pr3)
 
          ! Mole fraction of HNO3 in the binary solution:
          a = kn1 + kn2*ust
          b = kn3 + kn4*ust
          xnb = (-a - SQRT(a*a - 4.0_dp*b*(kn5 + kn6*ust + kn7*tlog - &
                pwlog)))/(2.0_dp*b)
 
          ! Concentration of HNO3 [mol/kg]:
          mnb=55.51_dp*xnb/(1.0_dp - xnb)
 
          ! Effective Henry constant of HNO3 [mol kg^-1 atm^-1]:
          hnb = EXP(qn1 + qn2*tr2 + (qn3 + qn4*tr + qn5*tr2 + qn6*tr3)*pr + &
               (qn7 + qn8*tr + qn9*tr2)*pr2 + qn10*tr*pr3)
 
          mnbmsb = mnb*msb
          mnb2 = mnb*mnb
 
          a = (tt*hnb*mnb2 - tt*hsb*mnbmsb - 2.0_dp*mnb*mnbmsb + mnbmsb*msb + &
               hnb*mnbmsb*pn0 - hsb*msb*msb*pn0)/ (mnb2 - mnbmsb)
          b = msb*(-2.0_dp*tt*hnb*mnb + tt*hsb*msb + mnbmsb - hnb*msb*pn0)/ &
                (mnb-msb)
          c = (tt*hnb*mnbmsb*msb)/(mnb - msb)
          d = a*a - 3.0_dp*b
          e = -2.0_dp*a*a*a + 9.0_dp*a*b - 27.0_dp*c
 
          phi= ATAN(SQRT(4.0_dp*d*d*d - e*e)/e)
          IF (phi .LT. 0) phi = phi + xpi
 
          ! Concentration of H2SO4 and HNO3 in the ternary solution [mol/kg]:
          ms = -(1.0_dp/3.0_dp)*(a + 2.0_dp*SQRT(d)*COS((xpi + phi)/3.0_dp))
          mn = mnb*(1.0_dp - ms/msb)
 
          ! Weight fraction of H2SO4 and HNO3 in the ternary solution [kg/kg]:
          ws = ms*mh2so4/(1.0_dp + ms*mh2so4 + mn*mhno3)
          wn = mn*mhno3/(1.0_dp + ms*mh2so4 + mn*mhno3)
       ENDIF
 
       ! Switch off the aerosols if ws is negative:
       IF (ws .LE. 0.0_dp) THEN
          pscliqtest = 0
          RETURN
       ENDIF
 
       ! Calculate the liquid solution density [g/cm3]:
       CALL density_liq(ws,wn,t,rho)
 
       ! Total specific liquid aerosol volume [cm3/cm3]:
       vliq = ns*98.076E-6_dp/(ws*rho)
 
       ! Preparations for calculation of the Henry coefficients:
       sqrtwn = SQRT(wn)
       sqrtws = SQRT(ws)
       wnws = wn*ws
       ws2 = ws*ws
       wn2 = wn*wn
       algm = LOG(1000.0_dp + 98.076_dp*ms + 63.012_dp*mn)
 
       ! Solubility of HCl (H* [mol/kg/atm] adapted from Luo et al., GRL, 1995;
       ! calculated concentrations assume that HCl is trace component of the
       ! aerosol):
       IF (phcl0 .GT. 0.0_dp) THEN
          ! Solubility of HCl in aerosol [mol/kg/atm]:
          hhcl = EXP(-(21.0_dp + 46.610_dp*wn + 4.0690_dp*ws - &
               4.8370_dp*sqrtwn + 2.1860_dp*sqrtws - 63.00_dp*wn2 - &
               40.170_dp*wnws - 1.5710_dp*ws2) - &
               ust*(-7437.0_dp - 8327.80_dp*wn + 1300.90_dp*ws + &
               1087.20_dp*sqrtwn - 242.710_dp*sqrtws + 18749.0_dp*wn2 + &
               18500.0_dp*wnws + 5632.0_dp*ws2) - LOG(wn+0.610_dp*ws) - &
               alg36 + algm)*1.013E3_dp
 
          ! Concentration of HCl in aerosol [mol/kg]:
          mcl = phcl0/(r*tlocal*ns/ms + 1.0_dp/hhcl)
 
          ! Weight fraction of HCl in aerosol [kg/kg]:
          wcl=mcl*mhcl/(1._dp + mh2so4*ms + mhno3*mn)
       ENDIF
 
       ! Solubility of HBr (H* [mol/kg/atm] adapted from Luo et al., GRL, 1995;
       ! Calculated concentrations assume that HBr is trace component of the
       ! aerosol):
       IF (phbr0 .GT. 0.0_dp) THEN
          ! Solubility of HBr in aerosol [mol/kg/atm]:
          hhbr = EXP(-(17.83_dp + 1.02_dp*wn - 1.08_dp*ws + 3.9_dp*sqrtwn + &
               4.38_dp*sqrtws - 8.87_dp*wn2 - 17.0_dp*wnws +3.73_dp*ws2) - &
               ust*(-8220.50_dp - 362.76_dp*wn + 658.93_dp*ws - &
               914.0_dp*sqrtwn - 955.3_dp*sqrtws + 9976.6_dp*wn2 + &
               19778.5_dp*wnws + 7680.0_dp*ws2) - LOG(wn+0.410_dp*ws) - &
               alg80 + algm)*1.013E3_dp
 
          ! Concentration of HBr in aerosol [mol/kg]:
          mbr = phbr0/(r*tlocal*ns/ms + 1.0_dp/hhbr)
 
          ! Weight fraction of HBr in aerosol [kg/kg]:
          wbr = mbr*mhbr/(1._dp + mh2so4*ms + mhno3*mn)
       ENDIF
 
       ! Solubility of HOCl (H* (mol/kg/atm) from Huthwelker et al., JAS, 1995;
       ! as an approximation, assume H* depends upon total molality;
       ! solubility of HOCl is low enough to ignore gas phase removal):
       IF (phocl0 .GT. 0.0_dp) THEN
          ! Solubility of HOCl in aerosol [mol/kg/atm]:
          hhocl = EXP(6.49460_dp - (-0.041070_dp + 54.56_dp*ust)*(ms+mn) - &
               5862.0_dp*(us298 - ust))
 
          ! Concentration of HOCl in aerosol [mol/kg]:
          mhocl = phocl0/(r*tlocal*ns/ms + 1.0_dp/hhocl)
 
          ! Weight fraction of HCl in aerosol [kg/kg]:
          whocl = mhocl*mhocl/(1._dp + mh2so4*ms + mhno3*mn)
       ENDIF
 
       ! The solubility of HOBr (solubility of HOBr is not known for all
       ! stratospheric conditions; limited data (Hanson&Ravishankara 1995 al)
       ! 210 K, 60 wt% H2SO4 indicate that hhobr=approx 18*hhocl. for HOBr, an
       ! effective Henry's law constant =18*hhocl is used; solubility of HOBr
       ! is low enough to ignore gas phase removal):
       IF (phobr0 .GT. 0.0_dp .AND. hhocl .GT. 0.0_dp) THEN
          ! Solubility of HOBr in aerosol [mol/kg/atm]:
          hhobr = 18.0_dp*hhocl
 
          ! Concentration of HOBr in aerosol [mol/kg]:
          mhobr = phobr0/(r*tlocal*ns/ms + 1.0_dp/hhobr)
 
          ! Weight fraction of HOBr in aerosol [kg/kg]:
          whobr = mhobr*mhobr/(1._dp + mh2so4*ms + mhno3*mn)
       ENDIF
 
     END SUBROUTINE lacomp
 
 
     SUBROUTINE hetla
 
       ! This subroutine first calculates the 1st order heterogeneous reaction
       ! coefficients on liquid aerosols of the following reactions:
       !
       ! N2O5   + H2O -> 2HNO3                 (ak150s/hiso(1))
       ! BrONO2 + H2O -> HOBr + HNO3           (ak161s/hiso(3))
       !
       ! Gamma values are indicated by variables with prefix 'G', for
       ! example g154 is the gamma value of HOCl due to reaction with
       ! HCl in the droplets.
       !
       ! From the gamma values, first order rate constants are calculated.
       ! These have the form 'akxxxS' and have units s^-1 or cm^3 mol^-1 s^-1.
       ! For example, the loss of ClONO2 + HCl due to the heterogeneous reaction
       ! ClONO2+HCl -> Cl2+HNO3 is d(ClONO2)/dt (units molecule cm-3 s-1) =
       ! ak152S*[ClONO2]*[HCl], where [ClONO2] and [HCl] are the gas phase
       ! amounts of these species in units molecule cm^-3.
       !
       ! Note: The reaction coefficients of the other heterogeneous reactions
       !       on liquid aerosols are calculated by the subroutine *calc_hetre*,
       !       which is called by *rpart* (update for every iteration step of
       !       Newton-Raphson scheme).                                   !MSHET
 
       ! Local variables:
       REAL(dp) :: pw, rmode, chclliq, choclliq, chbrliq, chobrliq, factor, &
            grxn, g161
 
       ! Local parameters:
       REAL(dp), PARAMETER :: g150 = 0.1_dp
 
 
       ! Executable statements:
 
       ! Calculate the liquid phase diffusion coefficient dl of HOCl in a
       ! H2SO4/HNO3 solution [cm2/s]. This is also used for HBr and HOBr.
       CALL diffusion_hocl(dl,ws,wn,ms,mn,t,rho)
 
       ! H2O partial pressure [atm]:
       pw = zh2o*patm
 
       ! Mode radius of liquid aerosols:
!!$       rmode = (3.0_dp*vliq/(4.0_dp*api*xnp)* &
!!$            EXP(-4.5_dp*(algsigliq**2)))**(1.0_dp/3.0_dp)                !MSHET
 
       ! Mean radius of liquid aerosols [cm]:
!!$       rmean = rmode*EXP(0.5_dp*(algsigliq**2))
       rmean = r_mean*1e-4_dp ! [um->cm]

       ! Total area of liquid aerosols [cm^2 cm^-3 air]:
!!$       aliq = 4.0_dp*api*xnp*(rmode**2)*EXP(2.0_dp*(algsigliq**2))       !MSHET
       aliq = sadw * 1e-8_dp ![um2/cm3->cm2/cm3]

!!$       ! [cm^2 cm^-3] -> [um^2/cm^3] (for output only):
!!$       sadsts_d(jl,jk,krow) = 1.E08_dp*aliq
!!$ 
!!$       ! Ratio STS : SSA (=prescribed as boundary condition; for output only):
!!$       sadrat_d(jl,jk,krow) = sadsts_d(jl,jk,krow)/sadw
 
       ! Liquid phase concentrations [mol dm^-3 atm^-1]:
       chclliq  = wcl  *rho/mhcl
       choclliq = whocl*rho/mhocl
       chbrliq  = wbr  *rho/mhbr
       chobrliq = whobr*rho/mhobr
 
       ! Convert effective Henry's law constants to mol/dm^3/atm:
       factor = rho/(1.0_dp + ms*mh2so4 + mn*mhno3)
       hhcl  = hhcl *factor
       hhocl = hhocl*factor
       hhobr = hhobr*factor
       hhbr  = hhbr *factor
 
       ! Calculate H*:                                                  !MSHET2
       !IF (gasch(18) .GE. 1) THEN
       !   hstarhocl = choclliq/gasch(18)
       !ELSE
       !   hstarhocl = 0.0_dp
       !ENDIF
       IF (gasch(21) .GE. 1) THEN
          hstarhcl = chclliq/gasch(21)
       ELSE
          hstarhcl = 0.0_dp
       ENDIF
       !IF (gasch(25) .GE. 1) THEN
       !   hstarhbr = chbrliq/gasch(25)
       !ELSE
        !   hstarhbr = 0.0_dp
       !ENDIF
       !IF (gasch(26) .GE. 1) THEN
       !   hstarhobr = chobrliq/gasch(26)
       !ELSE
       !   hstarhobr = 0.0_dp
       !ENDIF
 
       ! H2O activity:
       ah2o = 0.01_dp*psurf_stand*pw/(10.0_dp**(9.217_dp - 2190.0_dp/ &
            (t - 12.70_dp)))
 
       ! N2O5 + H2O -> 2HNO3  (ak150s/hiso(1)):
       ak150s = 0.25_dp*g150*1400.1_dp*SQRT(t)*aliq
 
       ! BrONO2 + H2O -> HOBr + HNO3  (ak161s/hiso(3))
       ! (from Hanson et al., JGR, 9063-9069, 1996):
       grxn = 211.0_dp*(ah2o**1.37_dp)
       g161 = 0.84_dp*grxn/(grxn + 0.84_dp)    ! Reactive uptake coeff. gamma
       g161 = MAX(0.1_dp, g161)
       g161 = MIN(0.8_dp, g161)
       ak161s = 0.25_dp*g161*1221.4_dp*SQRT(t)*aliq
 
     END SUBROUTINE hetla
 
   END SUBROUTINE het_new
 
 
   SUBROUTINE calc_hetre(jl,jk,krow,hcl,shs2,hs1,hs2,hs3)
 
     ! This subroutine first calculates the 1st order heterogeneous reaction
     ! coefficients on liquid aerosols of the following reactions:
 
     !      ClONO2 + H2O -> HOCl + HNO3   (ak151s/shs2)
     !      ClONO2 + HCl -> Cl2  + HNO3   (ak152s/hs1)
     !      HOCl   + HCl -> Cl2  + H2O    (ak154s/hs2)
     !      HOBr   + HBr -> 2Br  + H2O    (ak155s)
     !      HOBr   + HCl -> BrCl + H2O    (ak156s/hs3)
     !      HOCl   + HBr -> BrCl + H2O    (ak157s)
     !
     ! Gamma values are indicated by variables with prefix 'G', for
     ! example g154 is the gamma value of HOCl due to reaction with
     ! HCl in the droplets.
     !
     ! From the gamma values, first order rate constants are calculated.
     ! These have the form 'akxxxS' and have units s^-1 or cm^3 mol^-1 s^-1.
     ! For example, the loss of ClONO2 + HCl due to the heterogeneous reaction
     ! ClONO2+HCl -> Cl2+HNO3 is d(ClONO2)/dt (units molecule cm-3 s-1) =
     ! ak152S*[ClONO2]*[HCl], where [ClONO2] and [HCl] are the gas phase
     ! amounts of these species in units molecule cm^-3.
 
     ! t, aliq, rmean, ah2o, dl, hhocl, hhbr, hhobr, hstarhcl, hstarhocl,
     ! hstarhbr, hstarhobr: determined in *helta* (module variables)
 
     ! hcl, hocl, hbr, hobr: concentrations from current iteration step of
     ! Newton-Raphson (*rpart*)
     !
     ! Note: The reaction coefficients of the heterogeneous reactions
     !           N2O5   + H2O -> 2HNO3           (ak150s/hiso(1))
     !           BrONO2 + H2O -> HOBr + HNO3     (ak161s/hiso(3))
     !       on liquid aerosols are calculated by the subroutine *hetla*,
     !       which is called by *het_new*.
 
     ! *calc_hetre* is called from *rpart*, modules/mo_socol_chem.f90
     ! (update of reaction coefficients for every iteration step of
     ! Newton-Raphson).
 
     ! Subroutine arguments:
     INTEGER, INTENT(in) :: jl, jk, krow
     REAL(dp), INTENT(in) :: hcl !, hocl, hbr, hobr
     REAL(dp), INTENT(out) :: shs2, hs1, hs2, hs3
 
     ! Local variables:
     REAL(dp) :: chocl, chobr, adivl, adivl2, cbar, p, q, qq, f, fq, g0, &
          gcalc, ge, gs, chclliq, choclliq, chbrliq, chobrliq, kf154, kf156, &
          kf155, kf157, g151, g152, g154, g155, g156, g157
 
     ! Local parameters:
 
     ! Second order solution rate constants:
     REAL(dp), PARAMETER :: ks154 = 1.E05_dp, ks155 = 1.E07_dp, &
          ks156 = 1.E05_dp, ks157 = 1.E06_dp
 
     REAL(dp), PARAMETER :: ksur = 576.0_dp, alpha = 0.30_dp, rho1 = 2.E03_dp
 
 
     ! Executable statements:
 
     chclliq = hstarhcl *hcl
     !choclliq = hstarhocl*hocl
     !chbrliq  = hstarhbr *hbr
     !chobrliq = hstarhobr*hobr
 
     ! 1st order rate constants for reaction in liquid [dm^3 mol^-1 s^-1]:
     kf154 = ks154*chclliq
     !kf155 = ks155*chobrliq
     kf156 = ks156*chclliq
     !kf157 = ks157*choclliq
 
     ! ClONO2 + HCl -> Cl2  + HNO3 (ak151s/shs2) and
     ! ClONO2 + H2O -> HOCl + HNO3 (ak152s/hs1):
 
     ! Parameterization directly taken from Hanson and Ravishankara, J. Phys.
     ! Chem., 98, 5728, 1994, Except for the function F - see following
     ! comments - and the HCl solubility, which is calculated according to
     ! Luo et al.
     ! The function F: The form of F used by Hanson and Ravishankara can
     ! 'explode' under certain conditions. It has been replaced here by
     ! a stable function that is accurate within about 4%. This is also
     ! the case for other reactions that follow.
 
     g0 = 1.18E-04_dp + 9.1E-03_dp*ah2o + 0.5_dp*ah2o*ah2o
     gs = ah2o*ksur*chclliq
     p = rho1*chclliq/ah2o
     gcalc = g0*SQRT(1.0_dp + p)
     adivl = rmean/(1.4E-06_dp*SQRT(1.0_dp/ah2o))
     adivl2 = adivl*adivl
     f = (adivl + 0.312_dp*adivl2)/(3.0_dp + adivl + 0.312_dp*adivl2)
     ge = 1.0_dp/(1.0_dp/(gs+f*gcalc) + 1.0_dp/alpha)
     g152 = ge*(gs + f*gcalc*p/(1.0_dp + p))/(gs + f*gcalc)
     g151 = ge - g152
     hs1 = 0.25_dp*g152*1474.0_dp*SQRT(t)*aliq    ! ClONO2 + H2O -> HOCl + HNO3
     shs2 = 0.25_dp*g151*1474.0_dp*SQRT(t)*aliq   ! ClONO2 + HCl -> Cl2 + HNO3
!!$ 
!!$     clacsts_d(jl,jk,krow) = hs1 +  shs2  ! ClONO2 activation (for output)
 
     ! HOCl + HCl -> Cl2 + H2O (ak154s/hs2):
     IF(hhocl .GT. 0.0_dp .AND. kf154 .GT. 0.0_dp .AND. dl .GT. 0.0_dp) THEN
        q = rmean*SQRT(kf154/dl)
        qq = q + 0.312_dp*q*q
        fq = qq/(3.0_dp + qq)
        cbar = SQRT(8._dp*argas/(api*amhocl*1.0E-07)*t)
        g154 = fq/(fq + cbar/(4.0_dp*hhocl*0.082_dp*t*SQRT(kf154*dl)))
        hs2 = 0.25_dp*g154*cbar*aliq
     ELSE
        hs2 = 0.0_dp
     ENDIF
 
     ! HOBr + HCl -> BrCl + H2O (ak156s/hs3):
     IF(hhobr .GT. 0.0_dp .AND. kf156 .GT. 0.0_dp .AND. dl .GT. 0.0_dp) THEN
        q = rmean*SQRT(kf156/dl)
        qq = q + 0.312_dp*q*q
        fq = qq /(3.0_dp + qq)
        cbar = SQRT(8._dp*argas/(api*amhobr*1.0E-07)*t)
        g156 = fq/(fq + cbar/(4.0_dp*hhobr*0.082_dp*t*SQRT(kf156*dl)))
        hs3 = 0.25_dp*g156*cbar*aliq
     ELSE
        hs3 = 0.0_dp
     ENDIF
 
 
     !! HOBr + HBr -> 2Br  + H2O (ak155s):              ! currently not used
     !
     !! K155 not multiplied by fhbr as integration is done only with gas
     !! phase concentrations.
     !
     !IF (hhbr .GT. 0.0_dp .AND. kf155 .GT. 0.0_dp) THEN
     !   ! Total molecular concentrations [molecule cm^-3]:
     !   chobr = MAX(conc(iihobr), denh*1.0E-15_dp)
     !
     !   q = rmean*SQRT(kf155/dl)
     !   qq = q + 0.312_dp*q*q
     !   fq = qq/(3.0_dp + qq)
     !   cbar = SQRT(8._dp*argas/(api*amhbr*1.0E-07)*t)
     !   g155 = fq/(fq + cbar/(4.0_dp*hhbr*0.082_dp*t*SQRT(kf155*dl)))
     !   ak155S = 0.25_dp*g155*cbar*aliq/chobr
     !ELSE
     !   ak155S = 0.0_dp
     !ENDIF
     !
     !
     !! HOCl + HBr -> BrCl + H2O (ak157s):              ! currently not used
     !
     !! K155 not multiplied by fhbr as integration is done only with gas
     !! phase concentrations.
     !
     !IF (hhbr .GT. 0.0_dp .AND. kf157 .GT. 0.0_dp) THEN
     !   ! Total molecular concentrations [molecule cm^-3]:
     !   chocl = MAX(conc(iihocl), denh*1.0E-13_dp)
     !
     !   q = rmean*SQRT(kf157/dl)
     !   qq = q + 0.312_dp*q*q
     !   fq = qq/(3.0_dp + qq)
     !   cbar = SQRT(8._dp*argas/(api*amhbr*1.0E-07)*t)
     !   g157 = fq/(fq + cbar/(4.0_dp*hhbr*0.082_dp*t*SQRT(kf157*dl)))
     !   ak157S = 0.25_dp*g157*cbar*aliq/chocl
     !ELSE
     !   ak157S = 0.0_dp
     !ENDIF
 
   END SUBROUTINE calc_hetre
 
 
   SUBROUTINE density_liq (ws,wn,te,rho)
 
     ! This subroutine calculates the density of a ternary solution rho [g/cm3]
     ! for given H2SO4- and HNO3 weight fractions (ws and wn). The following
     ! parameterization is fitted to 0.05<ws+wn<0.70 wt fraction, but
     ! extrapolates well. 185 < T [K].
 
     ! *density_liq* is called from *sadcalc* and *hetla*.
 
     ! Subroutine arguments:
     REAL(dp), INTENT(in) :: ws, wn, te
     REAL(dp), INTENT(out) :: rho
 
     ! Local variables:
     REAL(dp) :: te2, te3, w, w2, wh, v1, vs, vn, vmcal
     REAL(dp), PARAMETER :: &
          x1 = 2.393284E-02    , x2=-4.359335E-05  , x3=7.961181E-08  , &
          x4 = 0.0             , x5 =-0.198716351  , x6=1.39564574E-03, &
          x7 =-2.020633E-06    , x8 = 0.51684706   , x9=-3.0539E-03   , &
          x10= 4.505475E-06    , x11=-0.30119511   , x12=1.840408E-03 , &
          x13=-2.7221253742E-06, x14=-0.11331674116, x15=8.47763E-04  , &
          x16=-1.22336185E-06  , x17= 0.3455282    , x18=-2.2111E-03  , &
          x19= 3.503768245E-06 , x20=-0.2315332    , x21=1.60074E-03  , &
          x22=-2.5827835E-06
 
 
     ! Executable statements:
 
     te2 = te*te
     te3 = te*te2
     w = ws + wn
     w2 = w*w
     wh = 1.0_dp-w
     v1 = x1  + x2*te + x3 *te2 + x4*te3
     vs = x5  + x6*te + x7 *te2 + (x8  + x9 *te + x10*te2)*w + &
          (x11 + x12*te + x13*te2)*w2
     vn = x14 + x15*te + x16*te2 + (x17 + x18*te + x19*te2)*w + &
          (x20 + x21*te + x22*te2)*w2
 
     vmcal = (v1*wh)/18.0160_dp + (vs*ws)/98.080_dp + (vn*wn)/63.0160_dp
     rho = 0.001_dp/vmcal
 
   END SUBROUTINE density_liq
 
 
   SUBROUTINE diffusion_hocl(dl,ws,wn,ms,mn,te,rho)
 
     ! This subroutine calculates the liquid phase diffusion coefficient dl
     ! (diffusivity) of HOCl using the Houghton cubic cell model
     ! (see Luo et al, GRL 1994) with cell dimension 3.65 angstroms.
     !
     ! Note: This is in good agreement (+-10% WITH a composition dependent
     !       cell dimension as given in Huthwelker et al (J. AT. CHEM., 1995)
     !
     ! Parameters:
     !     c:    wt % of H2SO4
     !     visc: Viscosity as mole ratio of H2SO4 and HNO3 (SI units)
     !     te:   Temperature [K]
     !     dl:   Diffusivity of HOCl [cm^2 s^-1]
 
     ! *diffusion_hocl* is called from *hetla*.
 
     ! Subroutine arguments:
     REAL(dp), INTENT(in) :: ws, wn, ms, mn, te, rho
     REAL(dp), INTENT(out) :: dl
 
     ! Local variables:
     REAL(dp) :: a, b, c, c2, c3, xb, xn, t0, visc, visn, viss
 
 
     ! Executable statements:
 
     IF(ms .GT. 0.0_dp .OR. mn .GT. 0.0_dp) THEN
        c = ws*100.0_dp
        c2 = c*c
        c3 = c2*c
        a = EXP(-7.7221330_dp + c*3.773159E-02_dp)
        xb = 623.80820_dp + 5.2216060_dp*c - 8.085769E-02_dp*c2 + &
             2.1769575E-04_dp*c3
        t0 = 154.34660_dp - 0.95216940_dp*c - &
             2.6749929E-03_dp*c**2+1.984055E-04_dp*c3
        xn = 0.51860400_dp
 ! AS, ETH, 15 March 2011
 ! model may crash, if te is equal or below t0
 ! -> set minimum value for te: te_min = t0 + 8 K
 !!$       viss = 1.E-03_dp*a*te**xn*EXP(xb/(te-t0))
        viss = 1.E-03_dp*a*te**xn*EXP(xb/max(8._dp,te-t0))
        a = 0.026560_dp + 0.0019710_dp*mn + 0.00023760_dp*mn*mn
        xb = 735.70_dp
        t0 = 92.890_dp + 0.68480_dp*mn
        xn = -0.012750_dp*mn
 ! AS, ETH, 15 March 2011
 ! model may crash, if te is equal or below t0
 ! -> set minimum value for te: te_min = t0 + 8 K
 !!$       visn = 1.E-03_dp*a*te**xn*EXP(xb/(te-t0))
        visn = 1.E-03_dp*a*te**xn*EXP(xb/max(8._dp,te-t0))
        visc = (viss*ms + visn*mn)/(ms + mn)
        ! dl=3.37E-14_dp*te*rho/visc
        ! !!! correct bug - increase dl by 1E03 ???
        dl = 1.0E03_dp*3.37E-14_dp*te*rho/visc
     ELSE
        dl = 0.0_dp
     ENDIF
 
   END SUBROUTINE diffusion_hocl
 
 END MODULE mo_socol_hetero
 
 
 
 
 

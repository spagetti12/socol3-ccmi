MODULE mo_socol_constants

  ! Constants used for SOCOL not already contained in *mo_constnats*

  ! M. Schraner, ETH Zürich, January 2009

  USE mo_constants,          ONLY: avo, api
  USE mo_kind,               ONLY: dp

  REAL(dp), PARAMETER :: api180  = api/180._dp ! Conversion degrees -> radians

  REAL(dp), PARAMETER :: mro2    = 0.21_dp    ! Mixing ratio of oxygen [mol/mol]

  REAL(dp), PARAMETER :: amo2    = 31.998_dp  ! Molec. weight of ogygen [g/mol]
  REAL(dp), PARAMETER :: amco    = 28.01_dp   ! Molec. weight of CO     [g/mol]
  REAL(dp), PARAMETER :: amno    = 30.01_dp   ! Molec. weight of NO     [g/mol]
  REAL(dp), PARAMETER :: amno2   = 44.01_dp   ! Molec. weight of NOx    [g/mol]
  REAL(dp), PARAMETER :: amh2so4 = 98.076_dp  ! Molec. weight of H2SO4  [g/mol]
  REAL(dp), PARAMETER :: amhno3  = 63.012_dp  ! Molec. weight of HNO3   [g/mol]
  REAL(dp), PARAMETER :: amhcl   = 36.461_dp  ! Molec. weight of HCl    [g/mol]
  REAL(dp), PARAMETER :: amhocl  = 52.46_dp   ! Molec. weight of HOCl   [g/mol]
  REAL(dp), PARAMETER :: amhbr   = 80.918_dp  ! Molec. weight of HBr    [g/mol]
  REAL(dp), PARAMETER :: amhobr  = 96.91_dp   ! Molec. weight of HOBr   [g/mol]
  REAL(dp), PARAMETER :: amnat   = 117.0_dp   ! Molec. weight of NAT    [g/mol]
  REAL(dp), PARAMETER :: amch4   = 16.04_dp   ! Molec. weight of CH4    [g/mol]  ! eth_as_ch4
  REAL(dp), PARAMETER :: amc5h8  = 68.12_dp   ! Molec. weight of C5H8   [g/mol]  ! eth_as_tropchem
  REAL(dp), PARAMETER :: amch2o  = 30.03_dp   ! Molec. weight of CH2O   [g/mol]
  REAL(dp), PARAMETER :: amhcooh = 46.03_dp   ! Molec. weight of HCOOH  [g/mol]
  REAL(dp), PARAMETER :: amch3cooh = 60.05_dp ! Molec. weight of CH3COOH [g/mol]

  REAL(dp), PARAMETER :: rmodenat= 5.0E-04_dp ! Mode radius of NAT particles 
                                              ! [cm]
  REAL(dp), PARAMETER :: sigmanat= 1.0_dp
  REAL(dp), PARAMETER :: rhonat  = 1.62_dp    ! (old: 1.35) Density of NAT 
                                              ! [g/cm3]
  REAL(dp), PARAMETER :: nice    = 1.0E-01_dp ! Number of ice particles [1/cm3]
  REAL(dp), PARAMETER :: sigmaice= 1.0_dp
  REAL(dp), PARAMETER :: rhoice  = 0.928_dp   ! Density of ice [g/cm3]
  REAL(dp), PARAMETER :: sigmaliq= 1.8_dp     ! Sigma of liquid aerosols

  REAL(dp), PARAMETER :: nofracnox_bcond = 0.9_dp  ! Fraction of NO in NOx 
                                              ! used for the boundary conditions
                                              ! (surface, aircraft, lightning)
  REAL(dp), PARAMETER :: no2fracnox_bcond = 0.1_dp ! Fraction of NO2 in NOx 
                                              ! used for the boundary conditions
                                              ! (surface, aircraft, lightning)
 
  REAL(dp), PARAMETER :: avoinv = 1.E-03_dp/avo  ! Inverse Avogadro constant
  REAL(dp), PARAMETER :: atomweight = 1.66054402E-27 ! Atomar mass unit [kg]

  REAL(dp), PARAMETER :: eps_mr = EPSILON(1.0_dp)*1.E-12 ! "Epsilon" for mixing 
                                              ! ratio

  REAL(dp), PARAMETER :: cp = 1007.0_dp       ! Heat capacity of air for 
                                              ! constant pressure [J/kg/K]

  REAL(dp), PARAMETER :: psurf_stand = 101325.0_dp   ! Srandard surface pressure
                                                     ! [Pa]
  REAL(dp), PARAMETER :: dobson = 2.6871E20_dp   ! molec/m^2 of 1 D.U.
 
END MODULE mo_socol_constants
  

MODULE mo_socol_deposition

  ! Contains deposition velocities.
  !
  ! M. Schraner, ETH Zürich, February 2009

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PUBLIC

  ! Module variables:

  ! Deposition velocities over land, sea and ice [m/s] based on Hauglustaine
  ! et al., 2004.
  ! (Deposition over ice is not used in the current model version.)
  REAL(dp), DIMENSION(3), PARAMETER :: & 
       wcodep   = (/ 3.E-04_dp, 0._dp, 0._dp /), &
       wnodep   = (/ 1.6E-04_dp, 3.E-05_dp, 2.E-05_dp /), &
       wno2dep  = (/ 1.E-03_dp, 2.E-04_dp, 1.E-04_dp /), &
       wo3dep   = (/ 4.E-03_dp, 7.E-04_dp, 7.E-04_dp /), &
       whno3dep = (/ 4.E-02_dp, 1.E-02_dp, 5.E-03_dp /), &
       wh2o2dep = (/ 5.E-03_dp, 1.E-02_dp, 3.2E-03_dp /), &
       wh2dep   = (/ 4.5E-04_dp, 0._dp, 0._dp /) ! eth_as_20100527
                                                 ! assume constant ratio wh2dep:wcodep = 1.5
                                                 ! Hauglustaine et Ehhalt, JGR, 2002.

END MODULE mo_socol_deposition

  

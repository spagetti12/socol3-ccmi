MODULE mo_socol_interpo

  USE mo_kind,  ONLY: dp

  IMPLICIT NONE

  PUBLIC
  
  ! Weighting factors and month indices for interpolation in time for SOCOL
  ! boundary conditions. Values are calculated by subroutine 
  ! *time_weights_socol*, called from *time_sets_socol*
  ! (both in mo_socol_time_control.f90).
  !
  ! Legend: "_current": interpolation for fields used at every time step
  !         "_chem"   : interpolation for "chemistry_date", i.e. fields used by
  !                     chemistry module
  !         "_rad"    : interpolation for "radiation_date", i.e. fields used by
  !                     radiation module

  ! Martin Schraner, ETH Zurich, December 2008

  ! Weighting factors:
  REAL (dp) :: wgt1_current, wgt2_current, wgt1_chem, wgt2_chem, wgt1_rad, &
               wgt2_rad

  ! Monthly indices for data of type (preceeding month, current month, 
  ! following month):
  INTEGER :: m3w1_current, m3w2_current, m3w1_chem, m3w2_chem, m3w1_rad, &
             m3w2_rad

  ! Monthly indices for data of type (0,...,13):
  INTEGER :: yw1_current, yw2_current, yw1_chem, yw2_chem, yw1_rad, yw2_rad

END MODULE mo_socol_interpo

SUBROUTINE socol_col_o3 (krow, xt)

  ! Description:

  ! Calculates total ozone [D.U.] and NO2 column [molec/m^2].
  !
  ! *col_o3* is called from *mezon*.
  !
  ! Eugene Rozanov, PMOD/WRC Davos, original code
  ! Martin Schraner, ETH Zurich, June 2009: Modifications for SOCOLvs3.0

  USE mo_constants,               ONLY: avo, amd, g
  USE mo_kind,                    ONLY: dp
  USE mo_socol_constants,         ONLY: dobson
  USE mo_socol_dimensions,        ONLY: nbdim, nproma, nlev, ntrac
  USE mo_socol_gcmfields,         ONLY: presb
  USE mo_socol_streams,           ONLY: totoz_d !!, totno2_d
  USE mo_socol_tracers,           ONLY: idt_o3, idt_no2

  IMPLICIT NONE

  ! Subroutine arguments:
  INTEGER, INTENT(in) :: krow  ! Index of current grid point block
  REAL(dp), INTENT(in) :: xt(nbdim,nlev,ntrac)  ! Chemical species

  ! Local variables:
  INTEGER :: jk
  REAL(dp) :: conv, dpress(nbdim), zo3(nbdim), zno2(nbdim)


  ! Executable statements:

  zo3(1:nproma) = 0.0_dp
  zno2(1:nproma) = 0.0_dp

  conv = avo*1000.0_dp/(g*amd)

  DO jk = 1, nlev
     ! Pressure difference of two vertical model boundaries [Pa]:
     dpress(1:nproma) = 100._dp*(presb(1:nproma,jk+1) - presb(1:nproma,jk))

     zo3(1:nproma) = zo3(1:nproma) + xt(1:nproma,jk,idt_o3)*dpress(1:nproma)
     zno2(1:nproma) = zno2(1:nproma) + xt(1:nproma,jk,idt_no2)*dpress(1:nproma)
  END DO

  ! Transformation to D.U. (ozone) and molec/m^2 (NO2):
  totoz_d(1:nproma,krow) = zo3(1:nproma)*conv/dobson
!!$  totno2_d(1:nproma,krow) = zno2(1:nproma)*conv

END SUBROUTINE socol_col_o3
       

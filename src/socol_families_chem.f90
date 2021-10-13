SUBROUTINE socol_families_chem (kproma, kbdim, klev, ktrac, krow, pxtm1, pxtte)

  ! Description:
 
  ! Calculates chemical families (used for output only).
  !
  ! *socol_families_chem* is called from *call_chem2*, src/call_submodels.f90.   
  !
  ! Martin Schraner, ETH Zurich, July 2009

  USE mo_kind,                    ONLY: dp
!!$  USE mo_socol_streams,           ONLY: clox_d, cly_d, ccly_d, bry_d, cbry_d, &
!!$                                        nox_d, noy_d
  USE mo_socol_tracers,           ONLY: idt_n,     idt_no,     idt_no2,     &
                                        idt_hno3,  idt_no3,    idt_n2o5,    &
                                        idt_hno4,  idt_clno3,  idt_clo,     &
                                        idt_cl,    idt_hocl,   idt_odscls,  &
                                        idt_odscll,idt_hcl,    idt_cl2,     &
                                        idt_cl2o2, idt_br,     idt_bro,     &
                                        idt_hbr,   idt_hobr,   idt_brno3,   &
                                        idt_brcl,  idt_odsbr
  USE mo_time_control,            ONLY: delta_time

  IMPLICIT NONE

  ! Subroutine arguments:
  INTEGER, INTENT(in) :: krow, kproma, kbdim, klev, ktrac
  REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)


  ! Excecutable statements:

!!$  clox_d(1:kproma,:,krow) = pxtm1(1:kproma,:,idt_cl) + &
!!$       pxtm1(1:kproma,:,idt_clo) + pxtm1(1:kproma,:,idt_hocl) + &
!!$       2._dp* pxtm1(1:kproma,:,idt_cl2) + &
!!$       2._dp* pxtm1(1:kproma,:,idt_cl2o2) + pxtm1(1:kproma,:,idt_brcl) + &
!!$       delta_time * (pxtte(1:kproma,:,idt_cl) + pxtte(1:kproma,:,idt_clo) + &
!!$       pxtte(1:kproma,:,idt_hocl) + 2._dp* pxtte(1:kproma,:,idt_cl2) + &
!!$       2._dp* pxtte(1:kproma,:,idt_cl2o2) + pxtte(1:kproma,:,idt_brcl))
!!$  
!!$  cly_d(1:kproma,:,krow)  = clox_d(1:kproma,:,krow) + &
!!$       pxtm1(1:kproma,:,idt_hcl) + pxtm1(1:kproma,:,idt_clno3) + &
!!$       delta_time * (pxtte(1:kproma,:,idt_hcl) + pxtte(1:kproma,:,idt_clno3))
!!$  
!!$  ccly_d(1:kproma,:,krow)  = cly_d(1:kproma,:,krow) + &
!!$       pxtm1(1:kproma,:,idt_odscls) + pxtm1(1:kproma,:,idt_odscll) + &
!!$       delta_time * (pxtte(1:kproma,:,idt_odscls) + &
!!$       pxtte(1:kproma,:,idt_odscll))
!!$  
!!$  bry_d(1:kproma,:,krow) = pxtm1(1:kproma,:,idt_br) + &
!!$       pxtm1(1:kproma,:,idt_bro) + pxtm1(1:kproma,:,idt_hbr) + &
!!$       pxtm1(1:kproma,:,idt_hobr) + pxtm1(1:kproma,:,idt_brno3) + &
!!$       pxtm1(1:kproma,:,idt_brcl) + &
!!$       delta_time * (pxtte(1:kproma,:,idt_br) + pxtte(1:kproma,:,idt_bro) + &
!!$       pxtte(1:kproma,:,idt_hbr) + pxtte(1:kproma,:,idt_hobr) + &
!!$       pxtte(1:kproma,:,idt_brno3) + pxtte(1:kproma,:,idt_brcl))
!!$  
!!$  cbry_d(1:kproma,:,krow)  = bry_d(1:kproma,:,krow) + &
!!$       pxtm1(1:kproma,:,idt_odsbr) + &
!!$       delta_time * pxtte(1:kproma,:,idt_odsbr)
!!$  
!!$  nox_d(1:kproma,:,krow)  = pxtm1(1:kproma,:,idt_n) + &
!!$       pxtm1(1:kproma,:,idt_no) + pxtm1(1:kproma,:,idt_no2) + &
!!$       pxtm1(1:kproma,:,idt_no3) + 2._dp* pxtm1(1:kproma,:,idt_n2o5) + &
!!$       pxtm1(1:kproma,:,idt_hno4) + &
!!$       delta_time * (pxtte(1:kproma,:,idt_n) + pxtte(1:kproma,:,idt_no) + &
!!$       pxtte(1:kproma,:,idt_no2) + pxtte(1:kproma,:,idt_no3) + &
!!$       2._dp* pxtte(1:kproma,:,idt_n2o5) + pxtte(1:kproma,:,idt_hno4))
!!$  
!!$  noy_d(1:kproma,:,krow)  = nox_d(1:kproma,:,krow) + &
!!$       pxtm1(1:kproma,:,idt_hno3) + pxtm1(1:kproma,:,idt_clno3) + &
!!$       pxtm1(1:kproma,:,idt_brno3) + &
!!$       delta_time * (pxtte(1:kproma,:,idt_hno3) + &
!!$       pxtte(1:kproma,:,idt_clno3) + 2._dp* pxtte(1:kproma,:,idt_brno3))
!!$       ! NAT is included in HNO3 (see chem/em_chemini)

END SUBROUTINE socol_families_chem

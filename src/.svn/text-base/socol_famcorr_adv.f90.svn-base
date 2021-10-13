SUBROUTINE prepare_socol_famcorr_adv

  ! Description:
  
  ! Calculates Cly, Bry and NOy from their members before calling the 
  ! advection scheme. After the advection scheme, the transported families are 
  ! used to scale their transported family members, see
  ! see *apply_socol_famcorr_adv*.

  ! *prepare_socol_famcorr_adv* is called from *scan1*.

  ! Martin Schraner, ETH Zurich, December 2009

  USE mo_kind,                    ONLY: dp
  USE mo_memory_g1a,              ONLY: xtm1
  USE mo_socol_tracers,           ONLY: idt_n,     idt_no,     idt_no2,     &
                                        idt_hno3,  idt_no3,    idt_n2o5,    &
                                        idt_hno4,  idt_clno3,  idt_clo,     &
                                        idt_cl,    idt_hocl,   idt_odscls,  &
                                        idt_odscll,idt_hcl,    idt_cl2,     &
                                        idt_cl2o2, idt_br,     idt_bro,     &
                                        idt_hbr,   idt_hobr,   idt_brno3,   &
                                        idt_brcl,  idt_odsbr,               &
                                        idt_cly_adv, idt_bry_adv, idt_noy_adv

  IMPLICIT NONE


  ! Executable statements:

  ! Cly:
  xtm1(:,:,idt_cly_adv,:) = xtm1(:,:,idt_cl,:) + xtm1(:,:,idt_clo,:) + &
       xtm1(:,:,idt_hocl,:) + 2._dp*xtm1(:,:,idt_cl2,:) + &
       2._dp*xtm1(:,:,idt_cl2o2,:) + xtm1(:,:,idt_clno3,:) + &
       xtm1(:,:,idt_hcl,:) + xtm1(:,:,idt_brcl,:)

  ! Bry:
  xtm1(:,:,idt_bry_adv,:) = xtm1(:,:,idt_br,:) + xtm1(:,:,idt_bro,:) + &
       xtm1(:,:,idt_hbr,:) + xtm1(:,:,idt_hobr,:) + xtm1(:,:,idt_brno3,:) + &
       xtm1(:,:,idt_brcl,:)

  ! NOy:
  xtm1(:,:,idt_noy_adv,:) = xtm1(:,:,idt_n,:) + xtm1(:,:,idt_no,:) + &
       xtm1(:,:,idt_no2,:) + xtm1(:,:,idt_no3,:) + &
       2._dp*xtm1(:,:,idt_n2o5,:) + xtm1(:,:,idt_hno3,:) + &
       xtm1(:,:,idt_hno4,:) + xtm1(:,:,idt_clno3,:) + xtm1(:,:,idt_brno3,:)


END SUBROUTINE prepare_socol_famcorr_adv


SUBROUTINE apply_socol_famcorr_adv

  ! Description:
  
  ! Compares advected Cly, Bry, and NOy with the sum of their advected family
  ! members. This ratio is used to correct the individual family members.

  ! *prepare_socol_famcorr_adv* is called from *scan1*.

  ! Martin Schraner, ETH Zurich, December 2009

  USE mo_decomposition,           ONLY: ldc=>local_decomposition
  USE mo_kind,                    ONLY: dp
  USE mo_memory_g1a,              ONLY: xtm1
  USE mo_scan_buffer,             ONLY: xtte
  USE mo_socol_constants,         ONLY: eps_mr
  USE mo_socol_namelist,          ONLY: lfamily_correct_adv
  USE mo_socol_streams,           ONLY: famfixcl_d, famfixbr_d, famfixn_d
  USE mo_socol_tracers,           ONLY: idt_n,     idt_no,     idt_no2,     &
                                        idt_hno3,  idt_no3,    idt_n2o5,    &
                                        idt_hno4,  idt_clno3,  idt_clo,     &
                                        idt_cl,    idt_hocl,   idt_hcl,     &
                                        idt_cl2,   idt_cl2o2,  idt_br,      &
                                        idt_bro,   idt_hbr,    idt_hobr,    &
                                        idt_brno3, idt_brcl,   idt_odsbr,   &
                                        idt_cly_adv, idt_bry_adv, idt_noy_adv
  USE mo_time_control,            ONLY: time_step_len

  IMPLICIT NONE

  ! Local variables:
  REAL(dp) :: time_step_len_inv
  REAL(dp), DIMENSION(ldc%nproma,ldc%nlev,ldc%ngpblks) :: cly_adv, &
       cly_mem_adv, fcl0, fcl1, bry_adv, bry_mem_adv, fbr0, fbr1, noy_adv, &
       noy_mem_adv, fn0, fn1, brcl_adv_corr, brcl_adv_uncorr, clno3_adv_corr, &
       clno3_adv_uncorr, brno3_adv_corr, brno3_adv_uncorr


  ! Executable statements:

  IF (lfamily_correct_adv) THEN

     time_step_len_inv = 1./time_step_len


     ! 1. Compare advected family with sum of advected members at t+dt:

     ! Advected families at t+dt:
     cly_adv(:,:,:) = &
          xtm1(:,:,idt_cly_adv,:) + time_step_len*xtte(:,:,idt_cly_adv,:)
     bry_adv(:,:,:) = &
          xtm1(:,:,idt_bry_adv,:) + time_step_len*xtte(:,:,idt_bry_adv,:)
     noy_adv(:,:,:) = &
          xtm1(:,:,idt_noy_adv,:) + time_step_len*xtte(:,:,idt_noy_adv,:)

     ! Sum of advected members at t+dt:
     cly_mem_adv(:,:,:) = xtm1(:,:,idt_cly_adv,:) + &
          time_step_len * (xtte(:,:,idt_cl,:) + xtte(:,:,idt_clo,:) + &
          xtte(:,:,idt_hocl,:) + 2._dp*xtte(:,:,idt_cl2,:) + &
          2._dp*xtte(:,:,idt_cl2o2,:) + xtte(:,:,idt_clno3,:) + &
          xtte(:,:,idt_hcl,:) + xtte(:,:,idt_brcl,:))
     bry_mem_adv(:,:,:) = xtm1(:,:,idt_bry_adv,:) + &
          time_step_len * (xtte(:,:,idt_br,:) + xtte(:,:,idt_bro,:) + &
          xtte(:,:,idt_hbr,:) + xtte(:,:,idt_hobr,:) + xtte(:,:,idt_brno3,:) + &
          xtte(:,:,idt_brcl,:))
     noy_mem_adv(:,:,:) = xtm1(:,:,idt_noy_adv,:) + &
          time_step_len * (xtte(:,:,idt_n,:) + xtte(:,:,idt_no,:) + &
          xtte(:,:,idt_no2,:) + xtte(:,:,idt_no3,:) + &
          2._dp*xtte(:,:,idt_n2o5,:) + xtte(:,:,idt_hno3,:) + &
          xtte(:,:,idt_hno4,:) + xtte(:,:,idt_clno3,:) + xtte(:,:,idt_brno3,:))

     ! Calculate ratios:
     WHERE (cly_mem_adv .GE. eps_mr)
        fcl1 = cly_adv/cly_mem_adv
     ELSEWHERE
        fcl1(:,:,:) = 1._dp
     END WHERE
     fcl0(:,:,:) = time_step_len_inv*(fcl1(:,:,:) - 1._dp) 

     WHERE (bry_mem_adv .GE. eps_mr)
        fbr1 = bry_adv/bry_mem_adv
     ELSEWHERE
        fbr1 = 1._dp
     END WHERE
     fbr0 = time_step_len_inv*(fbr1 - 1._dp) 

     WHERE (noy_mem_adv .GE. eps_mr)
        fn1 = noy_adv/noy_mem_adv
     ELSEWHERE
        fn1 = 1._dp
     END WHERE
     fn0 = time_step_len_inv*(fn1 - 1._dp) 


     ! 2. Corrections

     ! 2.1 Correct BRCl, ClNO3, and BrNO3 (included in two families):
     
     ! BrCl (Bry and Cly family)
     ! Correct as Bry, because BrCl has a larger contribution there:
     brcl_adv_uncorr(:,:,:) = xtm1(:,:,idt_brcl,:) + &
          time_step_len*xtte(:,:,idt_brcl,:)                   ! uncorrected
     xtte(:,:,idt_brcl,:) = &
          fbr0*xtm1(:,:,idt_brcl,:) + fbr1*xtte(:,:,idt_brcl,:)
     brcl_adv_corr(:,:,:) = xtm1(:,:,idt_brcl,:) + &
          time_step_len*xtte(:,:,idt_brcl,:)                   ! corrected

     ! ClNO3 (Cly and NOy family)
     ! Correct as Cly, because ClNO3 has a larger contribution there:
     clno3_adv_uncorr(:,:,:) = xtm1(:,:,idt_clno3,:) + &
          time_step_len*xtte(:,:,idt_clno3,:)                  ! uncorrected
     xtte(:,:,idt_clno3,:) = &
          fcl0*xtm1(:,:,idt_clno3,:) + fcl1*xtte(:,:,idt_clno3,:)
     clno3_adv_corr(:,:,:) = xtm1(:,:,idt_clno3,:) + &
          time_step_len*xtte(:,:,idt_clno3,:)                  ! corrected

     ! BrNO3 (Bry and NOy family)
     ! Correct as Bry, because BrNO3 has a larger contribution there:
     brno3_adv_uncorr(:,:,:) = xtm1(:,:,idt_brno3,:) + &
          time_step_len*xtte(:,:,idt_brno3,:)                  ! uncorrected
     xtte(:,:,idt_brno3,:) = &
          fbr0*xtm1(:,:,idt_brno3,:) + fbr1*xtte(:,:,idt_brno3,:)
     brno3_adv_corr(:,:,:) = xtm1(:,:,idt_brno3,:) + &
          time_step_len*xtte(:,:,idt_brno3,:)                  ! corrected


     ! 2.2 Correct Cly members (except BrCl and ClNO3):

     cly_adv = cly_adv - brcl_adv_corr - clno3_adv_corr
     cly_mem_adv = cly_mem_adv - brcl_adv_uncorr - clno3_adv_uncorr
           
     ! Re-calculate ratios:
     WHERE (cly_adv .GE. eps_mr .AND. cly_mem_adv .GE. eps_mr)
        fcl1 = cly_adv/cly_mem_adv
     ELSEWHERE
        fcl1 = 1._dp
     END WHERE
     fcl0 = time_step_len_inv*(fcl1 - 1._dp) 

     xtte(:,:,idt_cl,:)    = &
          fcl0*xtm1(:,:,idt_cl,:)    + fcl1*xtte(:,:,idt_cl,:)     
     xtte(:,:,idt_clo,:)   = &
          fcl0*xtm1(:,:,idt_clo,:)   + fcl1*xtte(:,:,idt_clo,:)
     xtte(:,:,idt_hocl,:)  = &
          fcl0*xtm1(:,:,idt_hocl,:)  + fcl1*xtte(:,:,idt_hocl,:)
     xtte(:,:,idt_cl2,:)   = &
          fcl0*xtm1(:,:,idt_cl2,:)   + fcl1*xtte(:,:,idt_cl2,:)     
     xtte(:,:,idt_cl2o2,:) = &
          fcl0*xtm1(:,:,idt_cl2o2,:) + fcl1*xtte(:,:,idt_cl2o2,:)
     xtte(:,:,idt_hcl,:)   = &
          fcl0*xtm1(:,:,idt_hcl,:)   + fcl1*xtte(:,:,idt_hcl,:)


     ! 2.3 Correct Bry members (except BrCl and BrNO3):

     bry_adv = bry_adv - brcl_adv_corr - brno3_adv_corr
     bry_mem_adv = bry_mem_adv - brcl_adv_uncorr - brno3_adv_uncorr
           
     ! Re-calculate ratios:
     WHERE (bry_adv .GE. eps_mr .AND. bry_mem_adv .GE. eps_mr)
        fbr1 = bry_adv/bry_mem_adv
     ELSEWHERE
        fbr1 = 1._dp
     END WHERE
     fbr0 = time_step_len_inv*(fbr1 - 1._dp) 

     xtte(:,:,idt_br,:)    = &
          fbr0*xtm1(:,:,idt_br,:)    + fbr1*xtte(:,:,idt_br,:)     
     xtte(:,:,idt_bro,:)   = &
          fbr0*xtm1(:,:,idt_bro,:)   + fbr1*xtte(:,:,idt_bro,:)
     xtte(:,:,idt_hobr,:)  = &
          fbr0*xtm1(:,:,idt_hobr,:)  + fbr1*xtte(:,:,idt_hobr,:)
     xtte(:,:,idt_hbr,:)   = &
          fbr0*xtm1(:,:,idt_hbr,:)   + fbr1*xtte(:,:,idt_hbr,:)


     ! 2.4 Correct NOy members (except ClNO3 and BrNO3):

     noy_adv = noy_adv - clno3_adv_corr - brno3_adv_corr
     noy_mem_adv = noy_mem_adv - clno3_adv_uncorr - brno3_adv_uncorr
           
     ! Re-calculate ratios:
     WHERE (noy_adv .GE. eps_mr .AND. noy_mem_adv .GE. eps_mr)
        fn1 = noy_adv/noy_mem_adv
     ELSEWHERE
        fn1 = 1._dp
     END WHERE
     fn0 = time_step_len_inv*(fn1 - 1._dp) 

     xtte(:,:,idt_n,:)     = &
          fn0*xtm1(:,:,idt_n,:)      + fn1*xtte(:,:,idt_n,:)     
     xtte(:,:,idt_no,:)    = &
          fn0*xtm1(:,:,idt_no,:)     + fn1*xtte(:,:,idt_no,:)
     xtte(:,:,idt_no2,:)   = &
          fn0*xtm1(:,:,idt_no2,:)    + fn1*xtte(:,:,idt_no2,:)
     xtte(:,:,idt_no3,:)   = &
          fn0*xtm1(:,:,idt_no3,:)    + fn1*xtte(:,:,idt_no3,:)
     xtte(:,:,idt_n2o5,:)  = &
          fn0*xtm1(:,:,idt_n2o5,:)   + fn1*xtte(:,:,idt_n2o5,:)
     xtte(:,:,idt_hno3,:)  = &
          fn0*xtm1(:,:,idt_hno3,:)   + fn1*xtte(:,:,idt_hno3,:)  
     xtte(:,:,idt_hno4,:)  = &
          fn0*xtm1(:,:,idt_hno4,:)   + fn1*xtte(:,:,idt_hno4,:)


     ! 3. Save correction factors as streams:

     famfixcl_d = fcl1
     famfixbr_d = fbr1
     famfixn_d  = fn1

  END IF

END SUBROUTINE apply_socol_famcorr_adv

!*    namelist *QBOCTL* to control options of the QBO parametrisation    !MSQBO
!                                                                       
!      name        type      purpose                                       
!      ----        ----      -------                                       
!                                                                       
!     *lqbo*       LOGICAL   *true for QBO assimilation
!     *nudg0*      INTEGER   *basic amplitude of the nudging field, [1/s] 

NAMELIST /qboctl/ lqbo, qbonudg0
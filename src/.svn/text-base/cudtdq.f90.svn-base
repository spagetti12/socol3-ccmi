SUBROUTINE cudtdq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, ktrac, krow,  & ! eth_as_lpj
                  pxten,                             & ! eth_as_20090731
                  paphp1,   pten,     ptte,     pqte,                  &
                  pxtte,    pxtec,    pmfuxt,   pmfdxt,                &
                  pmfus,    pmfds,    pmfuq,    pmfdq,                 &
                  pmful,    pdmfup,   pdmfdp,   plude,                 &
                  pdpmel,   prfl,     psfl,                            &
                  pcpen,    pqtec,    pqude,                           &
                  prsfc,    pssfc,    paprc,    paprs,                 &
                  ! eth_as_scav
                  xn3depcv, cvdprec,  kconbot)
!
!
!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDTDQ* IS CALLED FROM *CUMASTR*
!
USE mo_kind,         ONLY : dp
USE mo_constants,    ONLY : alv, als, alf, tmelt, g
USE mo_control,        ONLY : lsocol
USE mo_tracer,       ONLY : trlist
USE mo_time_control, ONLY : delta_time, time_step_len ! eth_as_20090731
USE mo_socol_scav,    ONLY : nspec_scav, idx_spec ! eth_as_scav
USE mo_socol_tracers, ONLY : trac_chemspec, n_trac_chemspec
USE mo_socol_namelist,ONLY : lscav ! eth_as_scav
! eth_as_lpj+
USE mo_lpj_streams,  ONLY : tot_precip
! eth_as_lpj-

!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, ktrac, krow ! eth_as_lpj
!
LOGICAL  llo1
!
REAL(dp) :: ptte(kbdim,klev),        pqte(kbdim,klev),                 &
            pten(kbdim,klev),        paphp1(kbdim,klevp1),             &
            paprc(kbdim),            paprs(kbdim),                     &
            prsfc(kbdim),            pssfc(kbdim)
REAL(dp) :: pmfus(kbdim,klev),       pmfds(kbdim,klev),                &
            pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                &
            pmful(kbdim,klev),       plude(kbdim,klev),                &
            pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),               &
            pqtec(kbdim,klev),       pqude(kbdim,klev),                &
            pxtec(kbdim,klev),       prfl(kbdim)
REAL(dp) :: pdpmel(kbdim,klev),      psfl(kbdim)
REAL(dp) :: pcpen(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
REAL(dp) :: zalpha(kbdim,klev,ktrac),zgamma(kbdim,ktrac),     & ! eth_as_20090731
            ztend(kbdim,klev,ktrac), ztenm(kbdim,klev,ktrac), & ! eth_as_20090731
            pxten(kbdim,klev,ktrac)                             ! eth_as_20090731
REAL(dp) :: zeps1, ztmst ! eth_as_20090731
!
REAL(dp) :: zmelt(kbdim)
REAL(dp) :: zsheat(kbdim)
REAL(dp) :: pxtte(kbdim,klev,ktrac), pmfuxt(kbdim,klev,ktrac),         &
            pmfdxt(kbdim,klev,ktrac)
!
REAL(dp) :: zrcpm ! reciprocal value of specific heat of moist air
INTEGER  :: jl, jk, jt
REAL(dp) :: zdiagt, zalv, zdtdt, zdqdt, zdxtdt
!
! eth_as_scav+
INTEGER  :: ji, kconbot(kbdim)
REAL(dp) :: ZDXDEP, zdep, zdep1, zdep2
REAL(dp) :: xn3depcv(kbdim,klev,nspec_scav), cvdprec(kbdim,klev)
!!$REAL(dp) :: xn3depcv(kbdim,klev,3), cvdprec(kbdim,klev)
! eth_as_scav-
!
cvdprec = 0._dp ! eth_as_scav
!
!----------------------------------------------------------------------
!
!*    1.0          SPECIFY PARAMETERS
!                  ------------------
!
100 CONTINUE
  zeps1=(1.-1.E-10)    ! eth_as_20090731
  zdiagt=delta_time
  ztmst=time_step_len  ! eth_as_20090731
!
!
!----------------------------------------------------------------------
!
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------
!
200 CONTINUE
  DO 210 jl=1,kproma
     zmelt(jl)=0._dp
     zsheat(jl)=0._dp
210 END DO
!
! eth_as_20090731+
  DO 2100 jt=1,ktrac
     DO 2101 jl=1,kproma

        ! eth_as_scav+
        ZDEP1=0._dp
        ZDEP2=0._dp
!C
!C     H2O2
!C
!!$        IF (JT .LE. n_trac_chemspec) THEN
!!$           IF (trac_chemspec(jt) .EQ. idx_spec(3)) THEN
!!$              ZDEP1=XN3DEPCV(JL,1,3)
!!$              ZDEP2=XN3DEPCV(JL,KLEV,3)
!!$           END IF
!C
!C     HCL
!C
!!$           IF (trac_chemspec(jt) .EQ. idx_spec(6)) THEN
!!$              ZDEP1=XN3DEPCV(JL,1,2)
!!$              ZDEP2=XN3DEPCV(JL,KLEV,2)
!!$           ENDIF
!C
!C     HNO3
!C
!!$           IF (trac_chemspec(jt) .EQ. idx_spec(4)) THEN
!!$              ZDEP1=XN3DEPCV(JL,1,1)
!!$              ZDEP2=XN3DEPCV(JL,KLEV,1)
!!$           ENDIF
!!$        END IF

        IF (lscav) THEN
          IF (jt .LE. n_trac_chemspec) THEN
             DO ji = 1, nspec_scav
              
                IF (trac_chemspec(jt) .EQ. idx_spec(ji)) THEN
                   zdep1 = XN3DEPCV(JL,1,ji)
                   zdep2 = XN3DEPCV(JL,KLEV,ji)               
                END  IF
             END DO
          END IF
        END IF

        ! eth_as_scav-

        ztend(jl,1,jt)=(pmfuxt(jl,1,jt)+pmfdxt(jl,1,jt) + zdep1)             &
                      *ztmst*g/(paphp1(jl,2)-paphp1(jl,1))
        ztenm(jl,1,jt)=max(ztend(jl,1,jt),pxten(jl,1,jt)*zeps1)
        ztend(jl,klev,jt)=(pmfuxt(jl,klev,jt)+pmfdxt(jl,klev,jt) + zdep2)    &
                         *ztmst*g/(paphp1(jl,klev+1)-paphp1(jl,klev))
        ztenm(jl,klev,jt)=max(ztend(jl,klev,jt)                      &
                             ,pxten(jl,klev,jt)*zeps1)
        zgamma(jl,jt)=999.
2101    ENDDO
2100 ENDDO
!
  DO 2110 jt=1,ktrac
     DO 2111 jk=2,klev-1
        DO 2112 jl=1,kproma

           ! eth_as_scav+
           zdep = 0._dp

!C
!C     H2O2
!C
!!$           IF (jt .LE. n_trac_chemspec) THEN
!!$              IF (trac_chemspec(jt) .EQ. idx_spec(3)) THEN
!!$                 ZDEP=XN3DEPCV(JL,JK,3)
!!$              ENDIF
!C
!C     HCL
!C
!!$              IF (trac_chemspec(jt) .EQ. idx_spec(6)) THEN
!!$                 ZDEP=XN3DEPCV(JL,JK,2)
!!$              ENDIF
!C
!C     HNO3
!C
!!$              IF (trac_chemspec(jt) .EQ. idx_spec(4)) THEN
!!$                 ZDEP=XN3DEPCV(JL,JK,1)
!!$              ENDIF
!!$           END IF

           IF (lscav) THEN
             IF (jt .LE. n_trac_chemspec) THEN
                DO ji = 1, nspec_scav
  
                   IF (trac_chemspec(jt) .EQ. idx_spec(ji)) THEN
                      zdep = XN3DEPCV(JL,JK,ji)
                   END  IF
                END DO
             END IF
           END IF            
  
           ! eth_as_scav-

           ztend(jl,jk,jt)=(pmfuxt(jl,jk,jt)-pmfuxt(jl,jk+1,jt)      &
                          +pmfdxt(jl,jk,jt)-pmfdxt(jl,jk+1,jt) + zdep)      &
                          *ztmst*g/(paphp1(jl,jk+1)-paphp1(jl,jk))
           ztenm(jl,jk,jt)=max(ztend(jl,jk,jt),pxten(jl,jk,jt)*zeps1)
2112       ENDDO
2111    ENDDO
2110 ENDDO
!
  DO 2120 jt=1,ktrac
     DO 2121 jk=1,klev
        DO 2122 jl=1,kproma
           IF (pxten(jl,jk,jt).le.0..and.ztenm(jl,jk,jt).eq.0.) THEN
              zalpha(jl,jk,jt)=1.
           ELSE
              zalpha(jl,jk,jt)=(pxten(jl,jk,jt)*zeps1)/(ztenm(jl,jk,jt))
           ENDIF
           zgamma(jl,jt)=min(zgamma(jl,jt),zalpha(jl,jk,jt))
           zgamma(jl,jt)=min(zgamma(jl,jt),1.)
           zgamma(jl,jt)=max(zgamma(jl,jt),0.)
2122       ENDDO 
2121    ENDDO  
2120 ENDDO
! eth_as_20090731-

  DO 250 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 220 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*         &
                                  (pmfus(jl,jk+1)-pmfus(jl,jk)+        &
                                   pmfds(jl,jk+1)-pmfds(jl,jk)-        &
                                   alf*pdpmel(jl,jk)-                  &
                                   zalv*(pmful(jl,jk+1)-pmful(jl,jk)-  &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk))))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zdqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                                  (pmfuq(jl,jk+1)-pmfuq(jl,jk)+        &
                                   pmfdq(jl,jk+1)-pmfdq(jl,jk)+        &
                                   pmful(jl,jk+1)-pmful(jl,jk)-        &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
220     END DO
!
        IF (trlist% anyconv /= 0) THEN
           DO 2204 jt=1,ktrac
              IF (trlist% ti(jt)% nconv == 1) THEN
                DO 2202 jl=1,kproma
                   IF(ldcum(jl)) THEN
                     zdxtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                                 *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) &
                                  +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt)) &
                                  *zgamma(jl,jt)
                      ! eth_as_scav+
                      IF (lscav) THEN
                        IF (jt .LE. n_trac_chemspec) THEN
                           DO ji = 1, nspec_scav
                              IF (trac_chemspec(jt) .EQ. idx_spec(ji)) THEN
                                 ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
                                      XN3DEPCV(JL,JK,ji)                           
                                 ZDXTDT=ZDXTDT+ZDXDEP
                                 EXIT
                              ENDIF
                           END DO
!!$                         IF (trac_chemspec(jt) .EQ. idx_spec(4)) THEN
!!$                            ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
!!$                                 XN3DEPCV(JL,JK,1)                           
!!$                            ZDXTDT=ZDXTDT+ZDXDEP
!!$                         ENDIF
!!$                         IF (trac_chemspec(jt) .EQ. idx_spec(6)) THEN
!!$                            ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
!!$                                 XN3DEPCV(JL,JK,2)                           
!!$                            ZDXTDT=ZDXTDT+ZDXDEP
!!$                         ENDIF                         
!!$                         IF (trac_chemspec(jt) .EQ. idx_spec(3)) THEN
!!$                            ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
!!$                                 XN3DEPCV(JL,JK,3)                           
!!$                            ZDXTDT=ZDXTDT+ZDXDEP
!!$                         ENDIF
                        END IF
                      END IF
                     ! eth_as_scav-
                     pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
                   ENDIF
2202            END DO
              ENDIF
2204       END DO
        ENDIF
!
!
     ELSE
        DO 230 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*        &
                     (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-     &
                                 zalv*(pmful(jl,jk)+pdmfup(jl,jk)      &
                                +pdmfdp(jl,jk)+plude(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zdqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                        (pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+       &
                        (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
230     END DO
!
        IF (trlist% anyconv /= 0) THEN
           DO 2304 jt=1,ktrac
              IF (trlist% ti(jt)% nconv == 1) THEN
                DO 2302 jl=1,kproma
                   IF(ldcum(jl)) THEN
                      zdxtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                             *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
                      ! eth_as_scav+
                      IF (lscav) THEN
                        IF (jt .LE. n_trac_chemspec) THEN
!!$                         IF (trac_chemspec(jt) .EQ. idx_spec(4)) THEN
!!$                            ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
!!$                                 XN3DEPCV(JL,JK,1)                           
!!$                            ZDXTDT=ZDXTDT+ZDXDEP
!!$                         ENDIF
!!$                         IF (trac_chemspec(jt) .EQ. idx_spec(6)) THEN
!!$                            ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
!!$                                 XN3DEPCV(JL,JK,2)                           
!!$                            ZDXTDT=ZDXTDT+ZDXDEP
!!$                         ENDIF                         
!!$                         IF (trac_chemspec(jt) .EQ. idx_spec(3)) THEN
!!$                            ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
!!$                                 XN3DEPCV(JL,JK,3)                           
!!$                            ZDXTDT=ZDXTDT+ZDXDEP
!!$                         ENDIF
                           DO ji = 1, nspec_scav
                              IF (trac_chemspec(jt) .EQ. idx_spec(ji)) THEN
                                 ZDXDEP=-1.*(G/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK)))* & 
                                      XN3DEPCV(JL,JK,ji)                           
                                 ZDXTDT=ZDXTDT+ZDXDEP
                                 EXIT
                   ENDIF
                           END DO
                        END IF
                      END IF
                     ! eth_as_scav-
                      pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt &
                           * zgamma(jl,jt)
                   ENDIF
2302            END DO
              END IF
2304       END DO
        ENDIF
!
     END IF

     DO 2305 jl = 1, kproma                                  
        ! CVDPREC(JL,JK)=PDMFUP(JL,JK)                    
        CVDPREC(JL,JK)=PDMFUP(JL,JK)+PDMFDP(JL,JK)       
2305 END DO
250 END DO
!
!
!---------------------------------------------------------------------
!
!      3.          UPDATE SURFACE FIELDS
!                  ---------------------
!
300 CONTINUE
  DO 310 jl=1,kproma
     prsfc(jl)=prfl(jl)
     pssfc(jl)=psfl(jl)
     paprc(jl)=paprc(jl)+zdiagt*(prfl(jl)+psfl(jl))
     paprs(jl)=paprs(jl)+zdiagt*psfl(jl)
     ! eth_as_lpj+
     if (lsocol) &
          tot_precip(jl,krow) = tot_precip(jl,krow) + zdiagt* &
          (prfl(jl)+psfl(jl)) * 86400._dp
     ! eth_as_lpj-
310 END DO
!
  RETURN
END SUBROUTINE cudtdq

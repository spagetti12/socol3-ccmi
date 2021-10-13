SUBROUTINE cuflx(    kproma, kbdim, klev, klevp1, krow,                     &
           pqen,     pqsen,    ptenh,    pqenh,                        &
           ktrac,                                                      &
           zpmfun,   &         ! eth_as_20090731
           pxtenh,   pmfuxt,   pmfdxt,                                 &
           paphp1,   pgeoh,                                            &
           kcbot,    kctop,    kdtop,                                  &
           ktype,    lddraf,   ldcum,                                  &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pmful,    plude,                        &
           pdmfup,   pdmfdp,   prfl,     prain,                        &
           pcpcu,                                                      &
           pten,     psfl,     pdpmel,   ktopm2,                       &
           pxten,    &         ! eth_as_20090731
           xn3depcv, cvdprec,  kconbot,                                &
           pxtu,     plu,      totliq ) ! eth_as_scav
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          EXTERNALS
!          ---------
!          NONE
!
USE mo_kind,         ONLY: dp
USE mo_constants,     ONLY: g, alf, cpd, tmelt, vtmpc2, rd
USE mo_physc2,       ONLY: cevapcu
USE mo_time_control, ONLY: time_step_len
USE mo_socol_scav,    ONLY: nspec_scav, idx_spec, zdg, dht, khcp
USE mo_socol_namelist,ONLY: lscav, lchem
USE mo_socol_streams, ONLY: cu_uvelo, xupdr
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac, krow
INTEGER, INTENT (OUT):: ktopm2
!
REAL(dp):: pqen(kbdim,klev),        pqsen(kbdim,klev),                 &
           ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           paphp1(kbdim,klevp1),    pgeoh(kbdim,klev)
!
REAL(dp):: pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           pmfus(kbdim,klev),       pmfds(kbdim,klev),                 &
           pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                 &
           pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                &
           pmful(kbdim,klev),       plude(kbdim,klev),                 &
           prfl(kbdim),             prain(kbdim)
REAL(dp):: pten(kbdim,klev),        pdpmel(kbdim,klev),                &
           psfl(kbdim)
REAL(dp):: pcpcu(kbdim,klev)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           kdtop(kbdim),            ktype(kbdim)
LOGICAL :: ldp(kbdim,klev) ! eth_as_20090731
LOGICAL :: lddraf(kbdim),           ldcum(kbdim)
REAL(dp):: pxtenh(kbdim,klev,ktrac),                                   &
           pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
REAL(dp):: zmass(kbdim,klev),  zpmfun(kbdim,klev) ! eth_as_20090731
REAL(dp):: pxten(kbdim,klev,ktrac) ! eth_as_20090731
!
INTEGER :: jl, jk, jt, ikb
REAL(dp):: zcons1, zcons2, zcucov, ztmelp2, zzp, zfac, zsnmlt          &
         , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum, zdpevap
REAL(dp):: zpsubcl(kbdim)
!
! eth_as_scav+
INTEGER  :: kconbot(kbdim)
REAL(dp) :: fice, scavco, scavco2, tcelc, totwflx, zrhoa, henper, perfracl
REAL(dp) :: totliq(kbdim,klev)
REAL(dp) :: xn3depcv(kbdim,klev,nspec_scav), cvdprec(kbdim,klev) 
!!$REAL(dp) :: xn3depcv(kbdim,klev,3), cvdprec(kbdim,klev) 
REAL(dp) :: pxtu(kbdim,klev,ktrac), plu(kbdim,klev) 
! eth_as_scav-

!
!*             SPECIFY CONSTANTS
!
  zcons1=cpd/(alf*g*time_step_len)
  zcons2=1._dp/(g*time_step_len)
  zcucov=0.05_dp
  ztmelp2=tmelt+2._dp
!
! eth_as_scav+
  kconbot = -1
  xn3depcv = 0._dp
  cvdprec = 0._dp 
! eth_as_scav-
!
!
!
!*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
!                  ---------------------------------
!
100 CONTINUE
!  itop=klev
  DO 110 jl=1,kproma
!     itop=MIN(itop,kctop(jl))
     IF(.NOT.ldcum(jl).OR.kdtop(jl).LT.kctop(jl)) lddraf(jl)=.FALSE.
     IF(.NOT.ldcum(jl)) ktype(jl)=0
110 END DO
  ktopm2=1
! eth_as_20090731+
  DO 1100 jk=1,klev
     DO 1101 jl=1,kproma
        zmass(jl,jk)=0._dp
        ldp(jl,jk)=.FALSE.
1101    ENDDO
1100 ENDDO
!
  DO 1110 jt=1,ktrac
     DO 1111 jk=ktopm2,klev
        DO 1112 jl=1,kproma
           IF (ldcum(jl).and.jk.ge.kctop(jl)-1) THEN
             zmass(jl,jk)=zpmfun(jl,jk)
             IF (lddraf(jl).and.jk.ge.kdtop(jl)) THEN
               zmass(jl,jk)=zpmfun(jl,jk)+pmfd(jl,jk)
             ELSE
               pmfdxt(jl,jk,jt)=0._dp
             ENDIF
           ELSE
             pmfuxt(jl,jk,jt)=0._dp
             pmfdxt(jl,jk,jt)=0._dp
           ENDIF
1112       ENDDO
1111    ENDDO
1110 ENDDO
!
  DO 1120 jk=ktopm2,klev
     DO 1121 jl=1,kproma
           IF (zmass(jl,jk).GE.0.) ldp(jl,jk)=.TRUE.

1121    ENDDO
1120 ENDDO
! eth_as_20090731-
  DO 120 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 115 jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
           pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*                      &
                               (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
           pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
           IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
              pmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)*                   &
                               (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
           ELSE
              pmfd(jl,jk)=0._dp
              pmfds(jl,jk)=0._dp
              pmfdq(jl,jk)=0._dp
              pdmfdp(jl,jk-1)=0._dp
           END IF
        END IF
115  END DO
!
! eth_as_scav+
     IF (lscav) THEN
! original code from ECHAM4.L39(DLR)/CHEM
     DO 1151 jl=1,kproma                                               
        IF (ldcum(jl)) KCONBOT(JL)=KCBOT(JL)                           
        !                                                               
        ! In-cloud scavenging calculated in CUFLX; budget in CUDTDQ           
        !                                                                      
        IF (LDCUM(JL).AND.JK.GE.KCTOP(JL)-1.AND.JK.LE.KCBOT(JL)) THEN  
           ! determine fraction ice                                               
           TCELC=PTENH(JL,JK)-273.15_dp                                    
           FICE=0._dp                                                      
           IF (TCELC.LT.-30._dp) FICE=1._dp                                    
           IF (TCELC.GE.-30._dp .AND. TCELC.LE.0._dp) &                         
                FICE=-1._dp*TCELC/30._dp                                           
           TOTWFLX=PLU(JL,JK)*PMFU(JL,JK)+PLUDE(JL,JK)+PDMFUP(JL,JK)     
           SCAVCO=0._dp                                                   
           IF (TOTWFLX.GT.0._dp) SCAVCO=(1._dp-FICE)*PDMFUP(JL,JK)/TOTWFLX     
           SCAVCO=MIN(SCAVCO,1._dp)                                         
           SCAVCO=MAX(SCAVCO,0._dp)

!C     HNO3
!C
!!$           XN3DEPCV(JL,JK,1)=SCAVCO*PXTU(JL,JK,idx_spec(4))*PMFU(JL,JK)
!!$           PXTU(JL,JK,idx_spec(4))=(1.-SCAVCO)*PXTU(JL,JK,idx_spec(4))
!!$           PMFUXT(JL,JK,idx_spec(4))=PXTU(JL,JK,idx_spec(4))*PMFU(JL,JK)
!C
!C     HCL
!C
!!$           XN3DEPCV(JL,JK,2)=SCAVCO*PXTU(JL,JK,idx_spec(6))*PMFU(JL,JK)
!!$           PXTU(JL,JK,idx_spec(6))=(1.-SCAVCO)*PXTU(JL,JK,idx_spec(6))
!!$           PMFUXT(JL,JK,idx_spec(6))=PXTU(JL,JK,idx_spec(6))*PMFU(JL,JK)
!C
!C     H2O2
!C
!!$           SCAVCO2=0.
!!$           IF (TOTWFLX.GT.0..AND.PMFU(JL,JK).GT.0.) THEN
!!$              ZRHOA=PAPHP1(JL,JK)/(PTENH(JL,JK)*RD)
!!$              TOTLIQ=TOTWFLX/PMFU(JL,JK)*ZRHOA
!!$              HENPER=5.49E-7*EXP(-6620.*(1./PTENH(JL,JK)-1./298.))
!!$              PERFRACL=1.E-3/(HENPER/TOTLIQ+1.E-3)
!!$              SCAVCO2=SCAVCO*PERFRACL
!!$           ENDIF
!!$           XN3DEPCV(JL,JK,3)=SCAVCO2*PXTU(JL,JK,idx_spec(3))*PMFU(JL,JK)
!!$           PXTU(JL,JK,idx_spec(3))=(1.-SCAVCO2)*PXTU(JL,JK,idx_spec(3))
!!$           PMFUXT(JL,JK,idx_spec(3))=PXTU(JL,JK,idx_spec(3))*PMFU(JL,JK)
!!$
!!$
           DO jt = 1, nspec_scav

              SCAVCO2=0._dp                                                  
              IF (TOTWFLX.GT.0._dp.AND.PMFU(JL,JK).GT.0._dp) THEN                 
                 ZRHOA=PAPHP1(JL,JK)/(PTENH(JL,JK)*RD)                        
                 TOTLIQ(JL,JK)=TOTWFLX/PMFU(JL,JK)*ZRHOA                             
                 HENPER=(4.088E-2_dp / khcp(jt))*EXP(-dht(jt)*(1._dp/PTENH(JL,JK)-1._dp/298._dp))         
                 PERFRACL=1.E-3_dp/(HENPER/TOTLIQ(JL,JK)+1.E-3_dp)
                         
                 SCAVCO2=SCAVCO*PERFRACL
              ELSE
                 SCAVCO2 = SCAVCO
              END IF

              XN3DEPCV(JL,JK,jt)=SCAVCO2*PXTU(JL,JK,idx_spec(jt))*PMFU(JL,JK) 
              PXTU(JL,JK,idx_spec(jt))=(1._dp-SCAVCO2)*PXTU(JL,JK,idx_spec(jt))  
              PMFUXT(JL,JK,idx_spec(jt))=PXTU(JL,JK,idx_spec(jt))*PMFU(JL,JK)           
              
           END DO                
        ENDIF
1151 END DO
     END IF
! eth_as_scav-
!
     DO 1154 jt=1,ktrac
        DO 1152 jl=1,kproma
           IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
! eth_as_20090731+
!              pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)                        &
!                                     -pmfu(jl,jk)*pxtenh(jl,jk,jt)
!              IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
!                 pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)                     &
!                                     -pmfd(jl,jk)*pxtenh(jl,jk,jt)
!              ELSE
!                 pmfdxt(jl,jk,jt)=0._dp
!              ENDIF
!           ELSE
!              pmfuxt(jl,jk,jt)=0._dp
!              pmfdxt(jl,jk,jt)=0._dp
!           ENDIF
!
              IF (ldp(jl,jk)) THEN
                 pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)                     &
                                     -zmass(jl,jk)*pxten(jl,jk-1,jt)
              ELSE
                 pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)                     &
                                     -zmass(jl,jk)*pxten(jl,jk,jt)
              ENDIF
           ELSE
              pmfuxt(jl,jk,jt)=0._dp
              pmfdxt(jl,jk,jt)=0._dp
           ENDIF
! eth_as_20090731-
1152    END DO
1154 END DO
!
120 END DO
  DO 130 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 125 jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           ikb=kcbot(jl)
           zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/                     &
                         (paphp1(jl,klevp1)-paphp1(jl,ikb)))
           zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
           pmfu(jl,jk)=pmfu(jl,ikb)*zzp
           pmfus(jl,jk)=pmfus(jl,ikb)*zzp
           pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
           pmful(jl,jk)=pmful(jl,ikb)*zzp
        END IF
125  END DO
!
! eth_as_20090731+
!      DO 1254 jt=1,ktrac
! !DIR$ IVDEP
! !OCL NOVREC
!         DO 1252 jl=1,kproma
!            IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
!               ikb=kcbot(jl)
!               zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/                   &
!                             (paphp1(jl,klevp1)-paphp1(jl,ikb))
!               zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
!               pmfuxt(jl,jk,jt)=pmfuxt(jl,ikb,jt)*zzp
!            ENDIF
! 1252    END DO
! 1254 END DO
! eth_as_20090731
!
130 END DO
!
!
!*    2.            CALCULATE RAIN/SNOW FALL RATES
!*                  CALCULATE MELTING OF SNOW
!*                  CALCULATE EVAPORATION OF PRECIP
!                   -------------------------------
!
  ! eth_as_03062013+
  ! Calculate 'updraft' =  mass flux / density 
  ! xupdr 
  !
  if (lchem) then
     do jl=1,kproma
        if (kctop(jl).ne.0) then
           do jk=kctop(jl),kcbot(jl)
              zrhoa=paphp1(jl,jk)/(ptenh(jl,jk)*rd)
              xupdr(jl,jk,krow)=pmfu(jl,jk)
              cu_uvelo(jl,jk,krow)=pmfu(jl,jk)/zrhoa
           enddo
        endif
     enddo
  end if
  !  eth_as_03062013-

!!$200 CONTINUE

  DO 210 jl=1,kproma
     prfl(jl)=0._dp
     psfl(jl)=0._dp
     prain(jl)=0._dp
210 END DO
  DO 220 jk=ktopm2,klev
     DO 215 jl=1,kproma
        IF(ldcum(jl)) THEN
           prain(jl)=prain(jl)+pdmfup(jl,jk)
           IF(pten(jl,jk).GT.tmelt) THEN
              prfl(jl)=prfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
              IF(psfl(jl).GT.0._dp.AND.pten(jl,jk).GT.ztmelp2) THEN
                 zfac=zcons1*(1._dp+vtmpc2*pqen(jl,jk))                &
                             *(paphp1(jl,jk+1)-paphp1(jl,jk))
                 zsnmlt=MIN(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
                 pdpmel(jl,jk)=zsnmlt
                 psfl(jl)=psfl(jl)-zsnmlt
                 prfl(jl)=prfl(jl)+zsnmlt
              END IF
           ELSE
              psfl(jl)=psfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
           END IF
        END IF
215  END DO
220 END DO
  DO 230 jl=1,kproma
     prfl(jl)=MAX(prfl(jl),0._dp)
     psfl(jl)=MAX(psfl(jl),0._dp)
     zpsubcl(jl)=prfl(jl)+psfl(jl)
230 END DO
  DO 240 jk=ktopm2,klev
     DO 235 jl=1,kproma
       IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_dp)  &
                                                                  THEN
           zrfl=zpsubcl(jl)
           zrnew=(MAX(0._dp,SQRT(zrfl/zcucov)-                         &
                        cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
                        MAX(0._dp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
           zrmin=zrfl-zcucov*MAX(0._dp,0.8_dp*pqsen(jl,jk)-pqen(jl,jk))&
                        *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
           zrnew=MAX(zrnew,zrmin)
           zrfln=MAX(zrnew,0._dp)
           zdrfl=MIN(0._dp,zrfln-zrfl)
           pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
           zpsubcl(jl)=zrfln
       END IF
235  END DO
240 END DO
  DO 250 jl=1,kproma
     zrsum=prfl(jl)+psfl(jl)
     zdpevap=zpsubcl(jl)-zrsum
     prfl(jl)=prfl(jl)+zdpevap*prfl(jl)*(1._dp/MAX(1.e-20_dp,zrsum))
     psfl(jl)=psfl(jl)+zdpevap*psfl(jl)*(1._dp/MAX(1.e-20_dp,zrsum))
250 END DO
!
  RETURN
END SUBROUTINE cuflx

SUBROUTINE dis(jl,krow,o3)

  ! Description:
  !
  ! This subroutine calculates photodissosiation rates (PR). The scheme is 
  ! based on precalculated PR0i = Fi(XO2,XO3), which is saved in a data file. 
  ! Calculation of actual PRi is performed by bi-linear interpolation of the
  ! PR0i data onto actual values of total oxygen (XO2) and ozone (XO3) above 
  ! the treated layer. Data needed: Concentration of O3, oxygen mixing ratio, 
  ! cosine of zenith angle.
  !
  ! *dis* is called from *chem*, src/socol_chemini.f90.
  ! 
  ! Eugene Rozanov, PMOD/WRC Davos, original code
  ! Martin Schraner, ETH Zurich, April 2009: Modifications for SOCOLvs3.0

  USE mo_constants,               ONLY: a
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_geoloc,                  ONLY: philon_2d, philat_2d
  USE mo_kind,                    ONLY: dp
  USE mo_socol_constants,         ONLY: mro2
  USE mo_socol_dimensions,        ONLY: nlev, nlevp1
  USE mo_socol_grid_calculations, ONLY: cosz, dens, zetb
  USE mo_socol_namelist,          ONLY: lsphericdependphot
  USE mo_socol_sun,               ONLY: do2a,       do3p,       do3d,       &
                                        dno,        dno2,       dhno3,      &
                                        dno3a,      dno3b,      dn2o5,      &
                                        dn2o,       dhno4,      dclno3a,    &
                                        dh2o2,      df11,       df12,       &
                                        dhocl,      dh2o,       dco2,       &
                                        dcl2,       dcl2o2,     dch2oa,     &
                                        dch2ob,     dch3o2h,    dbro,       &
                                        dbrno3a,    dbrno3b,    dbrcl,      &
                                        dhobr,      dcbrf3,     do2b,       &
                                        dch4,       dclno3b,    dhcl,       &
                                        dcfc113,    dcfc114,    dcfc115,    &
                                        dccl4,      dch3ccl3,   dhcfc22,    &
                                        dhcfc141b,  dhcfc142b,  dh1211,     &
                                        dch3br,     dch3cl,     dhcfc21,    &
                                        dhcfc123,   dh2402,     dchbr3,     &
                                        dch2br2,                            &
                                        dpan, dmacr, dhac, dmgly, dpaa , &
                                        tjo2a,      tjo3p,      tjo3d,      &
                                        tjno,       tjno2,      tjhno3,     &
                                        tjno3a,     tjno3b,     tjn2o5,     &
                                        tjn2o,      tjhno4,     tjclno3a,   &
                                        tjh2o2,     tjf11,      tjf12,      &
                                        tjhocl,     tjh2o,      tjco2,      &
                                        tjcl2,      tjcl2o2,    tjch2oa,    &
                                        tjch2ob,    tjch3o2h,   tjbro,      &
                                        tjbrno3a,   tjbrno3b,   tjbrcl,     &
                                        tjhobr,     tjcbrf3,    tjo2b,      &
                                        tjch4,      tjclno3b,   tjhcl,      &
                                        tjcfc113,   tjcfc114,   tjcfc115,   &
                                        tjccl4,     tjch3ccl3,  tjhcfc22,   &
                                        tjhcfc141b, tjhcfc142b, tjh1211,    &
                                        tjch3br,    tjch3cl,    tjhcfc21,   &
                                        tjhcfc123,  tjh2402,    tjchbr3,    &
                                        tjch2br2,                           &
                                        tjpan, tjmacr, tjhac, tjmgly, tjpaa, &
                                        xo2, xo3, &
                                        n_photolytab_o2, n_photolytab_o3, &
                                        sun_nopho

  IMPLICIT NONE

  ! Subroutine arguments:
  INTEGER, INTENT(in)  :: jl, krow
  REAL(dp), INTENT(in) :: o3(nlev)   ! Ozone [mol/mol]

  ! Local variables:
  INTEGER :: io2, io3, jk, jkk, jkp1, jkh, no2, no3, nlevsunny
  REAL(dp) :: re, zenit, sinza, sumo2, sumo3, h, hs, x, z1, z2, d1, d2, s1, &
       s2, rxo2, rxo3
  REAL(dp) :: vclo2(nlevp1), vclo3(nlevp1), clo2(nlevp1), clo3(nlevp1), &
       z(nlevp1+1)

  ! Photolysis rates at boundaries of model layers:
  REAL(dp) :: TJO2aL(nlevp1), TJO3DL(nlevp1), TJO3PL(nlevp1), TJNOL(nlevp1), &
       TJNO2L(nlevp1), TJHNO3L(nlevp1), TJNO3AL(nlevp1), TJNO3BL(nlevp1), &
       TJN2O5L(nlevp1), TJN2OL(nlevp1), TJHNO4L(nlevp1), TJCLNO3aL(nlevp1), &
       TJH2O2L(nlevp1), TJF11L(nlevp1), TJF12L(nlevp1), TJHOCLL(nlevp1), & 
       TJH2OL(nlevp1),TJCO2L(nlevp1), TJCL2L(nlevp1),TJCL2O2L(nlevp1), & 
       TJCH2OAL(nlevp1), TJCH2OBL(nlevp1), TJCH3O2HL(nlevp1), TJBROL(nlevp1), & 
       TJBRNO3aL(nlevp1), TJBRNO3bL(nlevp1),TJBRCLL(nlevp1), TJHOBRL(nlevp1), &
       TJCBRF3L(nlevp1), TJCH4L(nlevp1), TJHCLL(nlevp1), TJCLNO3bL(nlevp1), & 
       TJO2bL(nlevp1), TJCFC113L(nlevp1), TJCFC114L(nlevp1), & 
       TJCFC115L(nlevp1), TJCCL4L(nlevp1), TJCH3CCL3L(nlevp1), & 
       TJHCFC22L(nlevp1), TJHCFC141BL(nlevp1), TJHCFC142BL(nlevp1), & 
       TJH1211L(nlevp1), TJCH3BRL(nlevp1), TJCH3CLL(nlevp1), & 
       TJHCFC21L(nlevp1), TJHCFC123L(nlevp1), TJH2402L(nlevp1), & 
       TJCHBR3L(nlevp1), TJCH2BR2L(nlevp1), &
       tjpanl(nlevp1), tjmacrl(nlevp1), tjhacl(nlevp1), tjmglyl(nlevp1), tjpaal(nlevp1) 
 

  ! Executables statements:

  nlevsunny = nlev   ! Number of vertical levels with sun              !MSv3

  ! Initialize photolysis rates with zero:
  TJO2a      (:) = 0.0_dp
  TJO3P      (:) = 0.0_dp
  TJO3D      (:) = 0.0_dp
  TJNO       (:) = 0.0_dp
  TJNO2      (:) = 0.0_dp
  TJHNO3     (:) = 0.0_dp

  TJNO3A     (:) = 0.0_dp
  TJNO3B     (:) = 0.0_dp
  TJN2O5     (:) = 0.0_dp
  TJN2O      (:) = 0.0_dp
  TJHNO4     (:) = 0.0_dp
  TJCLNO3a   (:) = 0.0_dp

  TJH2O2     (:) = 0.0_dp
  TJF11      (:) = 0.0_dp
  TJF12      (:) = 0.0_dp
  TJHOCL     (:) = 0.0_dp
  TJH2O      (:) = 0.0_dp
  TJCO2      (:) = 0.0_dp

  TJCL2      (:) = 0.0_dp
  TJCL2O2    (:) = 0.0_dp  
  TJCH2OA    (:) = 0.0_dp
  TJCH2OB    (:) = 0.0_dp
  TJCH3O2H   (:) = 0.0_dp 
  TJBRO      (:) = 0.0_dp

  TJBRNO3a   (:) = 0.0_dp
  TJBRNO3b   (:) = 0.0_dp
  TJBRCL     (:) = 0.0_dp
  TJHOBR     (:) = 0.0_dp
  TJCBRF3    (:) = 0.0_dp 
  TJO2b      (:) = 0.0_dp

  TJCH4      (:) = 0.0_dp
  TJCLNO3b   (:) = 0.0_dp
  TJHCL      (:) = 0.0_dp 
  TJCFC113   (:) = 0.0_dp                                               !MSODS
  TJCFC114   (:) = 0.0_dp                                               !MSODS
  TJCFC115   (:) = 0.0_dp                                               !MSODS

  TJCCL4     (:) = 0.0_dp                                               !MSODS
  TJCH3CCL3  (:) = 0.0_dp                                               !MSODS
  TJHCFC22   (:) = 0.0_dp                                               !MSODS
  TJHCFC141B (:) = 0.0_dp                                               !MSODS
  TJHCFC142B (:) = 0.0_dp                                               !MSODS
  TJH1211    (:) = 0.0_dp                                               !MSODS

  TJCH3BR    (:) = 0.0_dp                                               !MSODS
  TJCH3CL    (:) = 0.0_dp                                               !MSODS
  TJHCFC21   (:) = 0.0_dp                                               !MSODS
  TJHCFC123  (:) = 0.0_dp                                               !MSODS
  TJH2402    (:) = 0.0_dp                                               !MSODS
  TJCHBR3    (:) = 0.0_dp                                               !MSODS
  TJCH2BR2   (:) = 0.0_dp                                               !MSODS

  tjpan(:) = 0.0_dp
  tjmacr(:) = 0.0_dp
  tjhac(:) = 0.0_dp
  tjmgly(:) = 0.0_dp
  tjpaa(:) = 0.0_dp

  TJO2al     (:) = 0.0_dp
  TJO3Pl     (:) = 0.0_dp
  TJO3Dl     (:) = 0.0_dp
  TJNOl      (:) = 0.0_dp
  TJNO2l     (:) = 0.0_dp
  TJHNO3l    (:) = 0.0_dp

  TJNO3Al    (:) = 0.0_dp
  TJNO3Bl    (:) = 0.0_dp
  TJN2O5l    (:) = 0.0_dp
  TJN2Ol     (:) = 0.0_dp
  TJHNO4l    (:) = 0.0_dp
  TJCLNO3al  (:) = 0.0_dp

  TJH2O2l    (:) = 0.0_dp
  TJF11l     (:) = 0.0_dp
  TJF12l     (:) = 0.0_dp
  TJHOCLl    (:) = 0.0_dp
  TJH2Ol     (:) = 0.0_dp
  TJCO2l     (:) = 0.0_dp

  TJCL2l     (:) = 0.0_dp
  TJCL2O2l   (:) = 0.0_dp
  TJCH2OAL   (:) = 0.0_dp
  TJCH2OBL   (:) = 0.0_dp
  TJCH3O2HL  (:) = 0.0_dp
  TJBROL     (:) = 0.0_dp

  TJBRNO3aL  (:) = 0.0_dp
  TJBRNO3bL  (:) = 0.0_dp
  TJBRCLL    (:) = 0.0_dp
  TJHOBRL    (:) = 0.0_dp
  TJCBRF3L   (:) = 0.0_dp
  TJO2bl     (:) = 0.0_dp

  TJCH4l     (:) = 0.0_dp
  TJCLNO3bl  (:) = 0.0_dp
  TJHCLl     (:) = 0.0_dp
  TJCFC113L  (:) = 0.0_dp                                               !MSODS
  TJCFC114L  (:) = 0.0_dp                                               !MSODS
  TJCFC115L  (:) = 0.0_dp                                               !MSODS

  TJCCL4L    (:) = 0.0_dp                                               !MSODS
  TJCH3CCL3L (:) = 0.0_dp                                               !MSODS
  TJHCFC22L  (:) = 0.0_dp                                               !MSODS
  TJHCFC141BL(:) = 0.0_dp                                               !MSODS
  TJHCFC142BL(:) = 0.0_dp                                               !MSODS
  TJH1211L   (:) = 0.0_dp                                               !MSODS

  TJCH3BRL   (:) = 0.0_dp                                               !MSODS
  TJCH3CLL   (:) = 0.0_dp                                               !MSODS
  TJHCFC21L  (:) = 0.0_dp                                               !MSODS
  TJHCFC123L (:) = 0.0_dp                                               !MSODS
  TJH2402L   (:) = 0.0_dp                                               !MSODS
  TJCHBR3L   (:) = 0.0_dp                                               !MSODS
  TJCH2BR2L  (:) = 0.0_dp                                               !MSODS

  tjpanl(:) = 0.0_dp
  tjmacrl(:) = 0.0_dp
  tjhacl(:) = 0.0_dp
  tjmglyl(:) = 0.0_dp
  tjpaal(:) = 0.0_dp

  ! Solar zenith angle:
  IF (lsphericdependphot .AND. cosz(jl) .GT. 0.0_dp) THEN
     zenit = (1224.0_dp * cosz(jl)* cosz(jl) + 1.0_dp)**0.5_dp / 35.0_dp
  ELSE
     zenit = cosz(jl)
  ENDIF
  sinza = SIN(ACOS(zenit))   ! sin of zenith angle

  ! Radius of the earth [km]:
  re = a*0.001_dp                
  
  ! Distance of model layers at model box boundaries from center of the earth
  ! [km]:
  z(:) = zetb(jl,:)+re

  ! Calculate vertical column of oxygen (vclo2) and ozone (vclo3) [molec/cm^2]
  ! at boundaries of model layers:

  vclo2(1) = 2.4E19_dp                                 ! O2 Column above 85 km
  vclo3(1) = 6.0E13_dp                                 ! O3 Column above 85 km

  DO jk = 2, nlevp1   
     vclo2(jk) = mro2 * &                                       ! [mol/mol]
          dens(jl,jk-1) * &                                     ! [molec/cm^3] 
          (z(jk)-z(jk+1))*1.0E5_dp                              ! [molec/cm^2]
     
     vclo3(jk) = o3(jk-1) * &                                   ! [mol/mol] 
          dens(jl,jk-1) * &                                     ! [molec/cm^3]
          (z(jk)-z(jk+1))*1.0E5_dp                              ! [molec/cm^2]
  ENDDO
  
  DO jk = 1, nlevp1

     ! Calculate column of oxygen (clo2) and ozone (clo3) [molec/cm^2] along the
     ! solar beam at boundaries of model layers:

     jkp1 = jk + 1
     h  = z(jkp1)*sinza
     hs = h*h
     
     ! 0 < zenith angle < 70:
     IF (zenit >= 0.3420201_dp) THEN
        IF(jk .EQ. 1) THEN
           clo2(1) = vclo2(1)/zenit
           clo3(1) = vclo3(1)/zenit
        ELSE
           clo2(jk) = clo2(jk-1) + vclo2(jk)/zenit
           clo3(jk) = clo3(jk-1) + vclo3(jk)/zenit
        ENDIF
     ENDIF

     ! 70 < zenith angle < 90:
     IF (zenit >= 1.6E-3_dp .AND. zenit < 0.3420201_dp) THEN
        sumo2 = 0.0_dp
        sumo3 = 0.0_dp

        DO jkk = 1, jkp1-1
           s1 = SQRT(z(jkk+1)*z(jkk+1)-hs)
           s2 = SQRT(z(jkk)  *z(jkk)-hs)               
           x = (z(jkk)-z(jkk+1))/(s2-s1)                ! New cos(zenith angle)
           sumo2 = sumo2 + vclo2(jkk)/x
           sumo3 = sumo3 + vclo3(jkk)/x           
        ENDDO

        clo2(jk)=sumo2
        clo3(jk)=sumo3        
     ENDIF
     
     ! zenith angle ~90:      
     IF (ABS(zenit) < 1.6E-3_dp ) THEN
        sumo2 = 0.0_dp
        sumo3 = 0.0_dp

        DO jkk = 1, jkp1-1
           x = (z(jkk)-z(jkk+1))/SQRT(z(jkk)*z(jkk)-z(jkk+1)*z(jkk+1))
           sumo2 = sumo2 + vclo2(jkk)/x          
           sumo3 = sumo3 + vclo3(jkk)/x 
        ENDDO     
   
        clo2(jk)=sumo2
        clo3(jk)=sumo3        
     ENDIF

     ! 90 < zenith angle < 98:
     IF (zenit <= -1.6E-3_dp .AND. zenit >= -0.139731_dp) THEN

        IF (h .LE. re) THEN   ! No solar radiation at model layer -> exit loop 
           clo2(jk) = clo2(jk-1)   ! for calculation of tjno(jk) (see below)
           clo3(jk) = clo3(jk-1)   ! for calculation of tjno(jk) (see below)
           nlevsunny = jk-1        ! number of levels with sun           !MSv3
           EXIT
        ENDIF
        
        ! Determine lowermost layer passed by the incoming solar beam
        ! (for zenith angle > 90!) -> vertical layer index jkh:
        DO jkk = 1, nlev+2           
           IF(z(jkk) .LE. h) EXIT ! Exit sub-loop
        ENDDO
        jkh = jkk

        ! O2/O3 column passed from sun to layer jkh:
        sumo2 = 0.0_dp
        sumo3 = 0.0_dp

        DO jkk = 1, jkh-1
           z2 = z(jkk)
           z1 = z(jkk+1)           
           s1 = (z1*z1-hs)
           s2 = (z2*z2-hs)        
           IF(s1 > 0.0_dp) THEN              
              x = SQRT(s2)-SQRT(s1)
           ELSE
              x = SQRT(s2)              
           ENDIF           
           x = (z2-z1)/x                    ! New cos(zenith angle)
           sumo2 = sumo2 + vclo2(jkk)/x
           sumo3 = sumo3 + vclo3(jkk)/x
        ENDDO
        
        ! O2/O3 column passed from layer jkh to layer jk:
        DO jkk = jkp1, jkh-1
           z2 = z(jkk)
           z1 = z(jkk+1)
           s1 = z1*z1-hs
           s2 = z2*z2-hs                    
           IF(s1 > 0.0_dp) THEN
              x = SQRT(s2)-SQRT(s1)
           ELSE
              x = SQRT(s2)
           ENDIF           
           x = (z2-z1)/x                    ! New cos(zenith angle)
           sumo2 = sumo2 + vclo2(jkk)/x
           sumo3 = sumo3 + vclo3(jkk)/x
        ENDDO
        
        clo2(jk)=sumo2
        clo3(jk)=sumo3
     ENDIF

     ! zenith angle > 98:
     IF(zenit < -0.139731_dp) THEN
        ! No solar radiation for any model layer -> exit subroutine
        RETURN
     ENDIF   
     
     ! Determine photolysis rates according to the O2/O3 column of the grid box
     ! by bilinear interpolation of the photolysis rates look up table data:
     ! (Probably this way is not the best one, isn't it?)

     ! Initial values for look up indices:
     no2 = 0
     no3 = 0
     
     ! Determine appropriate lookup indices for O3 (=no2) and O3 (=no3):
     DO io2 = 1, n_photolytab_o2-1
        IF(clo2(jk) .GE. xo2(io2) .AND. clo2(jk) .LE. xo2(io2+1)) THEN
           ! appropriate lookup index found
           no2 = io2
           EXIT
        ENDIF
     ENDDO
     DO io3 = 1, n_photolytab_o3-1
        IF(clo3(jk) .GE. xo3(io3) .AND. clo3(jk) .LE. xo3(io3+1)) THEN
           ! appropriate lookup index found
           no3 = io3
           EXIT
        ENDIF
     ENDDO
     
     ! Error message if no appropriate index is found:
     IF (no2 .EQ. 0) THEN
        WRITE(99,*) 'jk = ',jk
        CALL flush(99)
        WRITE(99,*) 'vclo2 = ',vclo2(jk)
        CALL flush(99)
        WRITE(99,*) 'clo2 = ',clo2(jk)
        CALL flush(99)
        WRITE(99,*) 'ZEN = ',zenit
        CALL flush(99)
        WRITE(99,*) &
             'Lon = ', philon_2d(jl,krow), 'Lat = ', philat_2d(jl,krow)
        CALL flush(99)
        WRITE(99,*) 'xo2 = ', xo2
        CALL flush(99)
        
        WRITE(message_text,*) 'jk = ',jk
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'vclo2 = ',vclo2(jk)
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'clo2 = ',clo2(jk)
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'ZEN = ',zenit
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) &
             'Lon = ', philon_2d(jl,krow), 'Lat = ', philat_2d(jl,krow)
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'xo2 = ', xo2
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'Out of table in dis xo2 '
        CALL message('',TRIM(message_text))
        CALL finish('dis','Run terminated')
     ENDIF     
     IF (no3 .EQ. 0) THEN
        WRITE(message_text,*) 'jk = ',jk
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'vclo3 = ',vclo3(jk)
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'clo3 = ',clo3(jk)
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'ZEN = ',zenit
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) &
             'Lon = ', philon_2d(jl,krow), 'Lat = ', philat_2d(jl,krow)
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'xo3 = ', xo3
        CALL message('',TRIM(message_text))
        WRITE(message_text,*) 'Out of table in dis xo3 '
        CALL message('',TRIM(message_text))
        CALL finish('dis','Run terminated')
     ENDIF

     ! Start bilinear interpolation:
     rxo2 = ((xo2(no2+1)-xo2(no2))/(clo2(jk)-xo2(no2)))    
     rxo3 = ((xo3(no3+1)-xo3(no3))/(clo3(jk)-xo3(no3)))
     
     d1 = DO2a(no2,no3) + (DO2a(no2+1,no3)-DO2a(no2,no3))/rxo2
     d2 = DO2a(no2,no3+1) + (DO2a(no2+1,no3+1)-DO2a(no2,no3+1))/rxo2
     TJO2aL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DO3D(no2,no3) + (DO3D(no2+1,no3)-DO3D(no2,no3))/rxo2
     d2 = DO3D(no2,no3+1) + (DO3D(no2+1,no3+1)-DO3D(no2,no3+1))/rxo2
     TJO3DL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DO3P(no2,no3) + (DO3P(no2+1,no3)-DO3P(no2,no3))/rxo2
     d2 = DO3P(no2,no3+1) + (DO3P(no2+1,no3+1)-DO3P(no2,no3+1))/rxo2
     TJO3PL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DNO(no2,no3) + (DNO(no2+1,no3)-DNO(no2,no3))/rxo2
     d2 = DNO(no2,no3+1) + (DNO(no2+1,no3+1)-DNO(no2,no3+1))/rxo2
     TJNOL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DNO2(no2,no3) + (Dno2(no2+1,no3)-DNO2(no2,no3))/rxo2
     d2 = DNO2(no2,no3+1) + (DNO2(no2+1,no3+1)-DNO2(no2,no3+1))/rxo2
     TJNO2L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DHNO3(no2,no3) + (DHNO3(no2+1,no3)-DHNO3(no2,no3))/rxo2
     d2 = DHNO3(no2,no3+1) + (DHNO3(no2+1,no3+1)-DHNO3(no2,no3+1))/rxo2
     TJHNO3L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DNO3a(no2,no3) + (DNO3a(no2+1,no3)-DNO3a(no2,no3))/rxo2
     d2 = DNO3a(no2,no3+1) + (DNO3a(no2+1,no3+1)-DNO3a(no2,no3+1))/rxo2
     TJNO3aL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DNO3b(no2,no3) + (DNO3b(no2+1,no3)-DNO3b(no2,no3))/rxo2
     d2 = DNO3b(no2,no3+1) + (DNO3b(no2+1,no3+1)-DNO3b(no2,no3+1))/rxo2
     TJNO3bL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DN2O5(no2,no3) + (DN2O5(no2+1,no3)-DN2O5(no2,no3))/rxo2
     d2 = DN2O5(no2,no3+1) + (DN2O5(no2+1,no3+1)-DN2O5(no2,no3+1))/rxo2
     TJN2O5L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DN2O(no2,no3) + (DN2O(no2+1,no3)-DN2O(no2,no3))/rxo2
     d2 = DN2O(no2,no3+1) + (DN2O(no2+1,no3+1)-DN2O(no2,no3+1))/rxo2
     TJN2OL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DHNO4(no2,no3) + (DHNO4(no2+1,no3)-DHNO4(no2,no3))/rxo2
     d2 = DHNO4(no2,no3+1) + (DHNO4(no2+1,no3+1)-DHNO4(no2,no3+1))/rxo2
     TJHNO4L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCLNO3a(no2,no3) + (DCLNO3a(no2+1,no3)-DCLNO3a(no2,no3))/rxo2
     d2 = DCLNO3a(no2,no3+1) + (DCLNO3a(no2+1,no3+1)-DCLNO3a(no2,no3+1))/rxo2
     TJCLNO3aL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DH2O2(no2,no3) + (DH2O2(no2+1,no3)-DH2O2(no2,no3))/rxo2
     d2 = DH2O2(no2,no3+1) + (DH2O2(no2+1,no3+1)-DH2O2(no2,no3+1))/rxo2
     TJH2O2L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DF11(no2,no3) + (DF11(no2+1,no3)-DF11(no2,no3))/rxo2
     d2 = DF11(no2,no3+1) + (DF11(no2+1,no3+1)-DF11(no2,no3+1))/rxo2
     TJF11L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DF12(no2,no3) + (DF12(no2+1,no3)-DF12(no2,no3))/rxo2
     d2 = DF12(no2,no3+1) + (DF12(no2+1,no3+1)-DF12(no2,no3+1))/rxo2
     TJF12L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DHOCL(no2,no3) + (DHOCL(no2+1,no3)-DHOCL(no2,no3))/rxo2
     d2 = DHOCL(no2,no3+1) + (DHOCL(no2+1,no3+1)-DHOCL(no2,no3+1))/rxo2
     TJHOCLL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DH2O(no2,no3) + (DH2O(no2+1,no3)-DH2O(no2,no3))/rxo2
     d2 = DH2O(no2,no3+1) + (DH2O(no2+1,no3+1)-DH2O(no2,no3+1))/rxo2
     TJH2OL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCO2(no2,no3) + (DCO2(no2+1,no3)-DCO2(no2,no3))/rxo2
     d2 = DCO2(no2,no3+1) + (DCO2(no2+1,no3+1)-DCO2(no2,no3+1))/rxo2
     TJCO2L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCL2(no2,no3) + (DCL2(no2+1,no3)-DCL2(no2,no3))/rxo2
     d2 = DCL2(no2,no3+1) + (DCL2(no2+1,no3+1)-DCL2(no2,no3+1))/rxo2
     TJCL2L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCL2O2(no2,no3) + (DCL2O2(no2+1,no3)-DCL2O2(no2,no3))/rxo2
     d2 = DCL2O2(no2,no3+1) + (DCL2O2(no2+1,no3+1)-DCL2O2(no2,no3+1))/rxo2
     TJCL2O2L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCH2OA(no2,no3) + (DCH2OA(no2+1,no3)-DCH2OA(no2,no3))/rxo2
     d2 = DCH2OA(no2,no3+1) + (DCH2OA(no2+1,no3+1)-DCH2OA(no2,no3+1))/rxo2
     TJCH2OAL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCH2OB(no2,no3) + (DCH2OB(no2+1,no3)-DCH2OB(no2,no3))/rxo2
     d2 = DCH2OB(no2,no3+1) + (DCH2OB(no2+1,no3+1)-DCH2OB(no2,no3+1))/rxo2
     TJCH2OBL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCH3O2H(no2,no3) + (DCH3O2H(no2+1,no3)-DCH3O2H(no2,no3))/rxo2
     d2=DCH3O2H(no2,no3+1)+(DCH3O2H(no2+1,no3+1)-DCH3O2H(no2,no3+1))/rxo2
     TJCH3O2HL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DBRO(no2,no3) + (DBRO(no2+1,no3)-DBRO(no2,no3))/rxo2
     d2 = DBRO(no2,no3+1) + (DBRO(no2+1,no3+1)-DBRO(no2,no3+1))/rxo2
     TJBROL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DBRNO3a(no2,no3) + (DBRNO3a(no2+1,no3)-DBRNO3a(no2,no3))/rxo2
     d2 = DBRNO3a(no2,no3+1) + (DBRNO3a(no2+1,no3+1)-DBRNO3a(no2,no3+1))/rxo2
     TJBRNO3aL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DBRNO3b(no2,no3) + (DBRNO3b(no2+1,no3)-DBRNO3b(no2,no3))/rxo2
     d2 = DBRNO3b(no2,no3+1) + (DBRNO3b(no2+1,no3+1)-DBRNO3b(no2,no3+1))/rxo2
     TJBRNO3bL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DBRCL(no2,no3) + (DBRCL(no2+1,no3)-DBRCL(no2,no3))/rxo2
     d2 = DBRCL(no2,no3+1) + (DBRCL(no2+1,no3+1)-DBRCL(no2,no3+1))/rxo2
     TJBRCLL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DHOBR(no2,no3) + (DHOBR(no2+1,no3)-DHOBR(no2,no3))/rxo2
     d2 = DHOBR(no2,no3+1) + (DHOBR(no2+1,no3+1)-DHOBR(no2,no3+1))/rxo2
     TJHOBRL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCBRF3(no2,no3) + (DCBRF3(no2+1,no3)-DCBRF3(no2,no3))/rxo2
     d2 = DCBRF3(no2,no3+1) + (DCBRF3(no2+1,no3+1)-DCBRF3(no2,no3+1))/rxo2
     TJCBRF3L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCLNO3b(no2,no3) + (DCLNO3b(no2+1,no3)-DCLNO3b(no2,no3))/rxo2
     d2 = DCLNO3b(no2,no3+1) + (DCLNO3b(no2+1,no3+1)-DCLNO3b(no2,no3+1))/rxo2
     TJCLNO3bL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DO2b(no2,no3) + (DO2b(no2+1,no3)-DO2b(no2,no3))/rxo2
     d2 = DO2b(no2,no3+1) + (DO2b(no2+1,no3+1)-DO2b(no2,no3+1))/rxo2
     TJO2bL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCH4(no2,no3) + (DCH4(no2+1,no3)-DCH4(no2,no3))/rxo2
     d2 = DCH4(no2,no3+1) + (DCH4(no2+1,no3+1)-DCH4(no2,no3+1))/rxo2
     TJCH4L(jk) = d1 + (d2 - d1)/rxo3

     d1 = DHCL(no2,no3) + (DHCL(no2+1,no3)-DHCL(no2,no3))/rxo2
     d2 = DHCL(no2,no3+1) + (DHCL(no2+1,no3+1)-DHCL(no2,no3+1))/rxo2
     TJHCLL(jk) = d1 + (d2 - d1)/rxo3

     d1 = DCFC113(no2,no3) + (DCFC113(no2+1,no3)-DCFC113(no2,no3))/rxo2 !MSODS
     d2 = DCFC113(no2,no3+1) + (DCFC113(no2+1,no3+1)-DCFC113(no2,no3+1))/rxo2
     TJCFC113L(jk) = d1 + (d2 - d1)/rxo3                                !MSODS

     d1 = DCFC114(no2,no3) + (DCFC114(no2+1,no3)-DCFC114(no2,no3))/rxo2 !MSODS
     d2 = DCFC114(no2,no3+1) + (DCFC114(no2+1,no3+1)-DCFC114(no2,no3+1))/rxo2
     TJCFC114L(jk) = d1 + (d2 - d1)/rxo3                                !MSODS

     d1 = DCFC115(no2,no3) + (DCFC115(no2+1,no3)-DCFC115(no2,no3))/rxo2 !MSODS
     d2 = DCFC115(no2,no3+1) + (DCFC115(no2+1,no3+1)-DCFC115(no2,no3+1))/rxo2
     TJCFC115L(jk) = d1 + (d2 - d1)/rxo3                                !MSODS

     d1 = DCCL4(no2,no3) + (DCCL4(no2+1,no3)-DCCL4(no2,no3))/rxo2       !MSODS
     d2 = DCCL4(no2,no3+1) + (DCCL4(no2+1,no3+1)-DCCL4(no2,no3+1))/rxo2
     TJCCL4L(jk) = d1 + (d2 - d1)/rxo3                                  !MSODS

     d1 = DCH3CCL3(no2,no3) + (DCH3CCL3(no2+1,no3)-DCH3CCL3(no2,no3))/rxo2
     d2 = DCH3CCL3(no2,no3+1) + &
          (DCH3CCL3(no2+1,no3+1)-DCH3CCL3(no2,no3+1))/rxo2              !MSODS
     TJCH3CCL3L(jk) = d1 + (d2 - d1)/rxo3                               !MSODS

     d1 = DHCFC22(no2,no3) + (DHCFC22(no2+1,no3)-DHCFC22(no2,no3))/rxo2 !MSODS
     d2 = DHCFC22(no2,no3+1) + (DHCFC22(no2+1,no3+1)-DHCFC22(no2,no3+1))/rxo2
     TJHCFC22L(jk) = d1 + (d2 - d1)/rxo3                                !MSODS

     d1 = DHCFC141B(no2,no3) + (DHCFC141B(no2+1,no3)-DHCFC141B(no2,no3))/rxo2
     d2 = DHCFC141B(no2,no3+1) + &
          (DHCFC141B(no2+1,no3+1)-DHCFC141B(no2,no3+1))/rxo2            !MSODS
     TJHCFC141BL(jk) = d1 + (d2 - d1)/rxo3                              !MSODS

     d1 = DHCFC142B(no2,no3) + (DHCFC142B(no2+1,no3)-DHCFC142B(no2,no3))/rxo2
     d2 = DHCFC142B(no2,no3+1) + &
          (DHCFC142B(no2+1,no3+1)-DHCFC142B(no2,no3+1))/rxo2            !MSODS
     TJHCFC142BL(jk) = d1 + (d2 - d1)/rxo3                              !MSODS

     d1 = DH1211(no2,no3) + (DH1211(no2+1,no3)-DH1211(no2,no3))/rxo2    !MSODS
     d2 = DH1211(no2,no3+1) + (DH1211(no2+1,no3+1)-DH1211(no2,no3+1))/rxo2
     TJH1211L(jk) = d1 + (d2 - d1)/rxo3                                 !MSODS

     d1 = DCH3BR(no2,no3) + (DCH3BR(no2+1,no3)-DCH3BR(no2,no3))/rxo2    !MSODS
     d2 = DCH3BR(no2,no3+1) + (DCH3BR(no2+1,no3+1)-DCH3BR(no2,no3+1))/rxo2
     TJCH3BRL(jk) = d1 + (d2 - d1)/rxo3                                 !MSODS

     d1 = DCH3CL(no2,no3) + (DCH3CL(no2+1,no3)-DCH3CL(no2,no3))/rxo2    !MSODS
     d2 = DCH3CL(no2,no3+1) + (DCH3CL(no2+1,no3+1)-DCH3CL(no2,no3+1))/rxo2
     TJCH3CLL(jk) = d1 + (d2 - d1)/rxo3                                 !MSODS

     d1 = DHCFC21(no2,no3) + (DHCFC21(no2+1,no3)-DHCFC21(no2,no3))/rxo2 !MSODS
     d2 = DHCFC21(no2,no3+1) + (DHCFC21(no2+1,no3+1)-DHCFC21(no2,no3+1))/rxo2
     TJHCFC21L(jk) = d1 + (d2 - d1)/rxo3                                !MSODS

     d1 = DHCFC123(no2,no3) + (DHCFC123(no2+1,no3)-DHCFC123(no2,no3))/rxo2
     d2 = DHCFC123(no2,no3+1) + &
          (DHCFC123(no2+1,no3+1)-DHCFC123(no2,no3+1))/rxo2              !MSODS
     TJHCFC123L(jk) = d1 + (d2 - d1)/rxo3                               !MSODS

     d1 = DH2402(no2,no3) + (DH2402(no2+1,no3)-DH2402(no2,no3))/rxo2    !MSODS
     d2 = DH2402(no2,no3+1) + (DH2402(no2+1,no3+1)-DH2402(no2,no3+1))/rxo2
     TJH2402L(jk) = d1 + (d2 - d1)/rxo3                                 !MSODS

     d1 = DCHBR3(no2,no3) + (DCHBR3(no2+1,no3)-DCHBR3(no2,no3))/rxo2    !MSODS
     d2 = DCHBR3(no2,no3+1) + (DCHBR3(no2+1,no3+1)-DCHBR3(no2,no3+1))/rxo2
     TJCHBR3L(jk) = d1 + (d2 - d1)/rxo3                                 !MSODS

     d1 = DCH2BR2(no2,no3) + (DCH2BR2(no2+1,no3)-DCH2BR2(no2,no3))/rxo2 !MSODS
     d2 = DCH2BR2(no2,no3+1) + (DCH2BR2(no2+1,no3+1)-DCH2BR2(no2,no3+1))/rxo2
     TJCH2BR2L(jk) = d1 + (d2 - d1)/rxo3                                !MSODS

     d1 = dpan(no2,no3) + (dpan(no2+1,no3)-dpan(no2,no3))/rxo2 
     d2 = dpan(no2,no3+1) + (dpan(no2+1,no3+1)-dpan(no2,no3+1))/rxo2
     tjpanl(jk) = d1 + (d2 - d1)/rxo3 
     
     d1 = dmacr(no2,no3) + (dmacr(no2+1,no3)-dmacr(no2,no3))/rxo2 
     d2 = dmacr(no2,no3+1) + (dmacr(no2+1,no3+1)-dmacr(no2,no3+1))/rxo2
     tjmacrl(jk) = d1 + (d2 - d1)/rxo3 
     
     d1 = dhac(no2,no3) + (dhac(no2+1,no3)-dhac(no2,no3))/rxo2 
     d2 = dhac(no2,no3+1) + (dhac(no2+1,no3+1)-dhac(no2,no3+1))/rxo2
     tjhacl(jk) = d1 + (d2 - d1)/rxo3 
     
     d1 = dmgly(no2,no3) + (dmgly(no2+1,no3)-dmgly(no2,no3))/rxo2 
     d2 = dmgly(no2,no3+1) + (dmgly(no2+1,no3+1)-dmgly(no2,no3+1))/rxo2
     tjmglyl(jk) = d1 + (d2 - d1)/rxo3 
     
     d1 = dpaa(no2,no3) + (dpaa(no2+1,no3)-dpaa(no2,no3))/rxo2 
     d2 = dpaa(no2,no3+1) + (dpaa(no2+1,no3+1)-dpaa(no2,no3+1))/rxo2
     tjpaal(jk) = d1 + (d2 - d1)/rxo3 

  ENDDO

  ! Calculate photolysis rates at centre of model layers:
  DO jk = 1, nlevsunny                                                  !MSv3

     TJNO(jk) = sun_nopho * & 
          EXP(-1.0E-8_dp*(0.5_dp*(clo2(jk)+clo2(jk+1)))**0.38_dp) * & 
          EXP(-5.E-19_dp*0.5_dp*(clo3(jk)+clo3(jk+1)))           !ROZSOL-MSSOL
     TJO2a(jk)     = 0.5_dp * (TJO2aL(jk)+TJO2aL(jk+1))
     TJO3D(jk)     = 0.5_dp * (TJO3DL(jk)+TJO3DL(jk+1))
     TJO3P(jk)     = 0.5_dp * (TJO3PL(jk)+TJO3PL(jk+1))
     TJNO2(jk)     = 0.5_dp * (TJNO2L(jk)+TJNO2L(jk+1))
     TJHNO3(jk)    = 0.5_dp * (TJHNO3L(jk)+TJHNO3L(jk+1))

     TJNO3a(jk)    = 0.5_dp * (TJNO3aL(jk)+TJNO3aL(jk+1))
     TJNO3b(jk)    = 0.5_dp * (TJNO3bL(jk)+TJNO3bL(jk+1))
     TJN2O5(jk)    = 0.5_dp * (TJN2O5L(jk)+TJN2O5L(jk+1))
     TJN2O(jk)     = 0.5_dp * (TJN2OL(jk)+TJN2OL(jk+1))
     TJHNO4(jk)    = 0.5_dp * (TJHNO4L(jk)+TJHNO4L(jk+1))
     TJCLNO3a(jk)  = 0.5_dp * (TJCLNO3aL(jk)+TJCLNO3aL(jk+1))

     TJH2O2(jk)    = 0.5_dp * (TJH2O2L(jk)+TJH2O2L(jk+1))
     TJF11(jk)     = 0.5_dp * (TJF11L(jk)+TJF11L(jk+1))
     TJF12(jk)     = 0.5_dp * (TJF12L(jk)+TJF12L(jk+1))
     TJHOCL(jk)    = 0.5_dp * (TJHOCLL(jk)+TJHOCLL(jk+1))
     TJH2O(jk)     = 0.5_dp * (TJH2OL(jk)+TJH2OL(jk+1))
     TJCO2(jk)     = 0.5_dp * (TJCO2L(jk)+TJCO2L(jk+1))

     TJCL2(jk)     = 0.5_dp * (TJCL2L(jk)+TJCL2L(jk+1))
     TJCL2O2(jk)   = 0.5_dp * (TJCL2O2L(jk)+TJCL2O2L(jk+1))
     TJCH2OA(jk)   = 0.5_dp * (TJCH2OAL(jk)+TJCH2OAL(jk+1))
     TJCH2OB(jk)   = 0.5_dp * (TJCH2OBL(jk)+TJCH2OBL(jk+1))
     TJCH3O2H(jk)  = 0.5_dp * (TJCH3O2HL(jk)+TJCH3O2HL(jk+1))
     TJBRO(jk)     = 0.5_dp * (TJBROL(jk)+TJBROL(jk+1))

     TJBRNO3a(jk)  = 0.5_dp * (TJBRNO3aL(jk)+TJBRNO3aL(jk+1))
     TJBRNO3b(jk)  = 0.5_dp * (TJBRNO3bL(jk)+TJBRNO3bL(jk+1))
     TJBRCL(jk)    = 0.5_dp * (TJBRCLL(jk)+TJBRCLL(jk+1))
     TJHOBR(jk)    = 0.5_dp * (TJHOBRL(jk)+TJHOBRL(jk+1))
     TJCBRF3(jk)   = 0.5_dp * (TJCBRF3L(jk)+TJCBRF3L(jk+1))
     TJO2b(jk)     = 0.5_dp * (TJO2bL(jk)+TJO2bL(jk+1))

     TJCLNO3b(jk)  = 0.5_dp * (TJCLNO3bL(jk)+TJCLNO3bL(jk+1))
     TJCH4(jk)     = 0.5_dp * (TJCH4L(jk)+TJCH4L(jk+1))
     TJHCL(jk)     = 0.5_dp * (TJHCLL(jk)+TJHCLL(jk+1))
     TJCFC113(jk)  = 0.5_dp * (TJCFC113L(jk)+TJCFC113L(jk+1))           !MSODS
     TJCFC114(jk)  = 0.5_dp * (TJCFC114L(jk)+TJCFC114L(jk+1))           !MSODS
     TJCFC115(jk)  = 0.5_dp * (TJCFC115L(jk)+TJCFC115L(jk+1))           !MSODS

     TJCCL4(jk)    = 0.5_dp * (TJCCL4L(jk)+TJCCL4L(jk+1))               !MSODS
     TJCH3CCL3(jk) = 0.5_dp * (TJCH3CCL3L(jk)+TJCH3CCL3L(jk+1))         !MSODS
     TJHCFC22(jk)  = 0.5_dp * (TJHCFC22L(jk)+TJHCFC22L(jk+1))           !MSODS
     TJHCFC141B(jk)= 0.5_dp * (TJHCFC141BL(jk)+TJHCFC141BL(jk+1))       !MSODS
     TJHCFC142B(jk)= 0.5_dp * (TJHCFC142BL(jk)+TJHCFC142BL(jk+1))       !MSODS
     TJH1211(jk)   = 0.5_dp * (TJH1211L(jk)+TJH1211L(jk+1))             !MSODS

     TJCH3BR(jk)   = 0.5_dp * (TJCH3BRL(jk)+TJCH3BRL(jk+1))             !MSODS
     TJCH3CL(jk)   = 0.5_dp * (TJCH3CLL(jk)+TJCH3CLL(jk+1))             !MSODS
     TJHCFC21(jk)  = 0.5_dp * (TJHCFC21L(jk)+TJHCFC21L(jk+1))           !MSODS
     TJHCFC123(jk) = 0.5_dp * (TJHCFC123L(jk)+TJHCFC123L(jk+1))         !MSODS
     TJH2402(jk)   = 0.5_dp * (TJH2402L(jk)+TJH2402L(jk+1))             !MSODS
     TJCHBR3(jk)   = 0.5_dp * (TJCHBR3L(jk)+TJCHBR3L(jk+1))             !MSODS
     TJCH2BR2(jk)  = 0.5_dp * (TJCH2BR2L(jk)+TJCH2BR2L(jk+1))           !MSODS

     tjpan(jk) = 0.5_dp * (tjpanl(jk) + tjpanl(jk+1))
     tjmacr(jk) = 0.5_dp * (tjmacrl(jk) + tjmacrl(jk+1))
     tjhac(jk) = 0.5_dp * (tjhacl(jk) + tjhacl(jk+1))
     tjmgly(jk) = 0.5_dp * (tjmglyl(jk) + tjmglyl(jk+1))
     tjpaa(jk)  = 0.5_dp * (tjpaal(jk) + tjpaal(jk+1))

  ENDDO

END SUBROUTINE dis

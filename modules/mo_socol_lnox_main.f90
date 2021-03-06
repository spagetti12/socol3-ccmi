MODULE mo_socol_lnox_main
!
!
!         Calculation of lightning produced NOx
!          
!         Lightnings are calculated using parameterization based on convective cloud top height. (Price & Rind, 1992)
!
!       Authors :
!          Timofei Sukhodolov PMOD / ETH IAC, 2013
!          Andrea Stenke PMOD / ETH IAC, 2013
 
  USE mo_kind,                      ONLY: dp
  USE mo_mpi,                       ONLY: p_pe, p_io
  USE mo_control,                   ONLY: nlev

  IMPLICIT NONE
  
  REAL(dp), ALLOCATABLE :: lightn_fac(:,:,:)          ! lightning factors

  LOGICAL, PUBLIC :: llnox_calc

  PUBLIC :: socol_lnox
  PUBLIC :: lnox_initialize
  PUBLIC :: lnox_close

  REAL(dp) :: PLT(16,4)
  data plt/1.0, 2.1, 3.9, 5.8, 7.7, 9.3, 10.5, 11.0, & !lnox vertical profiles
     11.0, 10.4, 9.2, 7.5, 5.5, 3.4, 1.5, 0.2, &
               2.4, 5.0, 7.4, 9.3, 10.6, 11.4, 11.5, 11.0, &
     9.9, 8.3, 6.3, 4.2, 2.2, 0.5, 0.0, 0.0, &
               0.2, 0.5, 0.6, 1.4, 2.7, 4.0, 5.0, 6.2, 8.6, &
     10.3, 11.6, 12.4, 12.7, 12.4, 7.6, 3.8, &
               0.6, 1.5, 2.9, 4.3, 5.4, 6.7, 7.7, 8.5, 9.6, &
     10.2, 10.5, 10.2, 8.2, 6.5, 4.5, 2.7/

CONTAINS

!==================================================================================
  subroutine socol_lnox(kproma,  ptenh,   klev,   kbdim,    klevp1,              &
       papp1,  paphp1,  kcbot,  kctop,    krow,      pten)
    !
    ! calculation of lightning flashes and NOx produced in these flashes
    ! 
    ! subroutine is called from cumastr.f90
    !
    USE mo_constants,     ONLY: api, rd, g
    USE mo_memory_g3b,    ONLY: slm                          !surface type
    USE mo_geoloc,        ONLY: gboxarea_2d, philat_2d       !grid box area, latitude
    USE mo_socol_streams,     ONLY: lfp_d, lnox_d
    
    
    IMPLICIT NONE
    
    INTEGER,  INTENT(IN) :: kproma, krow, kbdim,  klev, klevp1
    
    REAL(dp) :: ptenh(kbdim,klev),       paphp1(kbdim,klevp1),                  &  !temperature, pressure
         papp1(kbdim,klev),                                                              &
         grheight(kbdim, klev),   pten(kbdim,klev)                          !gridbox height, temperature
    
    INTEGER  :: kcbot(kbdim),            kctop(kbdim),                          &  !cloud top level, bottom level
         cu_freeze(kbdim)                                                   !0deg level
    
    
    INTEGER  :: jl, jk, jt, level

    REAL(dp) :: lfp, lnox
    
    REAL(dp) :: zrhoa
    
    REAL :: zcth                 ! cloud top height [m]
    REAL :: zcbh                 ! cloud bottom height [m]
    REAL :: zcdh                 ! 0degC height [m]
    REAL :: zcdhkm               ! 0degC height [km]
    REAL :: zcthkm               ! cloud top height [km]
    REAL :: pg                   ! cloud to ground lightning fraction
    REAL :: ffp                  ! flash frequency in air column ff [s^-1] Price

    REAL(dp) :: d2r
    
    REAL(dp), PARAMETER :: nox_per_flash = 2900. ! nox emission per flash in the tropics

    d2r = api/180._dp

    !  Calculate layers depth
    DO jl = 1, kproma
       DO jk = klev, 1, -1
          zrhoa=papp1(jl,jk)/(pten(jl,jk)*rd)
          grheight(jl,jk) =  (paphp1(jl,jk+1)-paphp1(jl,jk))/(zrhoa*g) ! m
       ENDDO
    ENDDO

    do jl=1,kproma
       !  Determine 0C cloud level for calculating fraction
       !  of InterCloud lightning
       cu_freeze(jl)=kctop(jl)
       do jk=kctop(jl),kcbot(jl)
          if (pten(jl,jk).le.273.15_dp) cu_freeze(jl)=jk
       enddo
    enddo

    vector_loop: DO jl=1,kproma

       lfp =  0._dp
       lnox =  0._dp

      ! if convective activity, calculate NOx prod. from lightning
      convect: IF ((kctop(jl) /= 0) .AND. (kctop(jl) < kcbot(jl)) .AND. (pten(jl,39) > 233._dp)) THEN
     
          ! calculate
          ! - cloud top height     : zcth [m],
          ! - cloud bottom height  : zcbh [m]
          ! - cloudheight below 0C : zcdh [m]
          zcth=0._dp
          zcbh=0._dp
          zcdh=0._dp
          
          DO level=klev,kctop(jl),-1
             zcth=zcth+grheight(jl,level)
          END DO
          DO level=klev,kcbot(jl),-1
             zcbh=zcbh+grheight(jl,level)
          END DO
          DO level=kctop(jl),cu_freeze(jl)
             zcdh=zcdh+grheight(jl,level)
          END DO
          
          !******************************************
          ! Calculate flash frequencies in air column

          !******************************************

          ffp = 0._dp

          !convert zcth from m to km
          zcthkm=zcth*1.e-3_dp ! m -> km

          if (paphp1(jl,kctop(jl)-1) < 40000.) then ! only, if cloud top above 400 hPa
               IF (slm(jl,krow) >= 0.1) THEN
                  ! land
                  ffp = zcthkm**4.9_dp * 3.44e-5_dp
               ELSE
                  ! sea and ice
                  ffp = zcthkm**1.73_dp * 6.4e-4_dp
               END IF
          endif
  
          !************************************************
          ! Calculate cloud to ground lightning fraction pg
          ! (in case IC and CG are distinguished)
          !   (lightning = flash frequency)
          !************************************************

             ! IC + CG
             IF (zcdh < 5500._dp) THEN
                ! only CG if 0degC level is below 5.5 km
                pg=1._dp
             ELSE
                IF (zcdh < 14000._dp) THEN
                   ! 0degC level is below 14 km
                   zcdhkm=zcdh*1.e-3_dp
                   pg=0.021_dp*zcdhkm-0.648_dp
                   pg=pg*zcdhkm+7.49_dp
                   pg=pg*zcdhkm-36.54_dp
                   pg=pg*zcdhkm+63.09_dp
                   pg=1._dp/pg
                ELSE
                   ! constant fraction of 2% CG in clouds where
                   ! 0degC level is above 14 km
                   pg=0.02_dp
                END IF
             END IF

             lfp =  ffp * cos(d2r*philat_2d(jl,krow))/cos(d2r*30.)*lightn_fac(jl,krow,1)

             !***********************************************
             !calculation of NOx 
             !flashes/min => g[NOx]/m**2/s
             !***********************************************
             lnox=lfp/gboxarea_2d(jl,krow)/60._dp
             
             if ((philat_2d(jl,krow).ge.-20).and.(philat_2d(jl,krow).le.20)) lnox=lnox*nox_per_flash
             
             if (((philat_2d(jl,krow).gt.-40).and.(philat_2d(jl,krow).lt.-20)).or.((philat_2d(jl,krow).gt.20).and.(philat_2d(jl,krow).lt.40))) &
                  lnox=lnox*nox_per_flash*(1+0.05*(abs(philat_2d(jl,krow))-20))
             
             if ((philat_2d(jl,krow).ge.40).or.(philat_2d(jl,krow).le.-40)) lnox=lnox*nox_per_flash*2.
             
             llnox_calc = .TRUE.
             
          END IF convect ! convective activity

          if (lfp.ne.lfp.or.lfp.lt.1.e-10_dp.or.lfp.gt.1.e10_dp) lfp=0._dp
          if (lnox.ne.lnox.or.lnox.lt.1.e-10_dp.or.lnox.gt.1.e15_dp) lnox=0._dp
          
          lfp_d(jl,krow) = lfp/60._dp ! 1/min -> 1/s
          lnox_d(jl,krow) = lnox
          
       END DO vector_loop
    
  end subroutine socol_lnox
! ***********************************************************************

subroutine lnox_initialize
!
! initialization of lightning parameterization
! called from stepon

  USE mo_decomposition, ONLY:  local_decomposition
  USE mo_time_control,  ONLY: current_date, get_date_components, next_date
  USE mo_socol_readfile,     ONLY: socol_read_netcdf

  IMPLICIT NONE
  INTEGER :: icurrentyear, inextyear, icurrentmonth,  inextmonth, inextday, icurrentday
  CHARACTER(25), PARAMETER :: lightn = 'lightn_fac'
  CHARACTER(31) :: varname
  CHARACTER(50) :: varname_longname

  IF (p_pe == p_io) write(*,*)'lnox initialization'
 
  CALL get_date_components (next_date,  month=inextmonth, year=inextyear,  &
       day=inextday                                               )

  varname = 'LFAC'
  varname_longname = 'Lightning factors'
  CALL socol_read_netcdf(lightn, varname, 'LONLAT', mo=inextmonth,&
  varname_longname=varname_longname, data3d=lightn_fac)

!  write(*,*)sum(lightn_fac(:,:,0)),sum(lightn_fac(:,:,1)),sum(lightn_fac(:,:,2))
!  stop 'TIMO'

  llnox_calc = .FALSE.

  end subroutine lnox_initialize

! ***********************************************************************
  subroutine lnox_close
    !
    ! called from stepon
    
    DEALLOCATE(lightn_fac)
  
  end subroutine lnox_close
! ***********************************************************************
  SUBROUTINE LNT1(FIN,ZIN,NZIN,FOUT,ZOUT,NZOUT)
    !
    !interpolation of vertical lnox profiles
    !called from socol_boundary_conditions
    INTEGER :: I, II, NZIN, NZOUT
    REAL(dp) :: FIN(NZIN),ZIN(NZIN),FOUT(NZOUT),ZOUT(NZOUT)
    REAL(dp) :: ZXP, ZXM, ZNZ
    DO 1 I=1,NZOUT
       II=1
       DO WHILE(ZIN(II)+0.1.LT.ZOUT(I))
          II=II+1
          IF(II.GT.NZIN) GOTO 1
       ENDDO
       IF(II.EQ.1)II=II+1
       ZXP=ZIN(II)-ZOUT(I)
       ZXM=ZOUT(I)-ZIN(II-1)
       ZNZ=1./(ZIN(II)-ZIN(II-1))
       FOUT(I)=(FIN(II-1)*ZXP+FIN(II)*ZXM)*ZNZ
1      CONTINUE
       
     END SUBROUTINE LNT1
     ! ***********************************************************************
     
   END MODULE mo_socol_lnox_main
   ! ***********************************************************************

MODULE mo_socol_qbo

  ! Description:

  ! Organises the QBO nudging.
  !
  ! Original code: Marco Giorgetta, MPI for Meteorology, Hamburg, ????
  ! Adaptions for SOCOLvs3.0: Martin Schraner, ETH Zurich, February 2009

  USE mo_exception,               ONLY: message, message_text
  USE mo_geoloc,                  ONLY: philat_2d
  USE mo_kind,                    ONLY: dp
  USE mo_socol_interpo,           ONLY: wgt1_current, wgt2_current, &
                                        m3w1_current, m3w2_current
  USE mo_socol_namelist
  USE mo_socol_grid_calculations, ONLY: presf_stand_socol, &
                                        socol_interpolate_index
  USE mo_socol_readfile,          ONLY: socol_read_netcdf
  USE mo_time_control,            ONLY: time_step_len

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  ! Pressure grid of data set:
  REAL(dp), ALLOCATABLE :: qbolev(:)   ! pressure grid of QBO data [hPa]
  INTEGER               :: nqbolev     ! extent of the pressure dimension
  INTEGER :: k0_qbonudg, k1_qbonudg    ! index of first/last level of ECHAM5
                                       ! grid in range of <qbolev>
  INTEGER, ALLOCATABLE  :: qboindex(:) ! index of next following level of 
                                       ! <qbolev> with respect to a given 
                                       ! ECHAM level

  ! Level depending half width of QBO [m/s]:
  REAL(dp), ALLOCATABLE :: qbohw_qbolev(:)   ! data pressure grid
  REAL(dp), ALLOCATABLE :: qbohw(:,:)        ! model pressure grid

  ! Damping constant for full nudging [1/s]:
  REAL(dp) :: qbodamp0

  ! Level depending damping constant of QBO [-]:
  REAL(dp), ALLOCATABLE :: qbodampv_qbolev(:)   ! data pressure grid
  REAL(dp), ALLOCATABLE :: qbodampv(:,:)        ! model pressure grid

  ! Monthly mean zonal wind observations at Equator [m/s]:
  ! (a) Preceeding, current and following month, pressure grid of data set:
  REAL(dp), ALLOCATABLE :: qboeq_qbolev_m3(:,:)
  ! (b) Interpolated to current time step, pressure grid of ECHAM5:
  REAL(dp), ALLOCATABLE :: qboeq(:,:)

  ! Latitudinal limits of QBO nudging: 
  ! Full nudging for |latitude| <= qbolat1 [degrees]
  ! No QBO nudging for |latitude| > qbolat2 [degrees]
  REAL(dp), PARAMETER ::  qbolat1=10._dp, qbolat2=20._dp

  PUBLIC :: read_socol_qbo, interpolate_socol_qbo, qbonudg_socol, &
       cleanup_socol_qbo

CONTAINS

  SUBROUTINE read_socol_qbo(yr,mo)

    ! Reads QBO observations for current, preceeding and following month.

    ! *read_socol_qbo* is called from *read_socol_bcond_m*,
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo

    ! Local variables:
    CHARACTER(3), PARAMETER  :: fn0 = 'qbo'
    CHARACTER(7)  :: fn 
    CHARACTER(31) :: varname
    CHARACTER(50) :: varname_longname

    ! Executable statements:

    IF (.NOT. lqbonudg) THEN
       WRITE (message_text,*) 'QBO nudging switched off'
       CALL message('',TRIM(message_text))

    ELSE
       WRITE (message_text,*) 'QBO nudging switched on'
       CALL message('',TRIM(message_text))

       ! Damping constant for full nudging [1/s]:
       qbodamp0 = 1./tauqbonudg/24._dp/3600._dp

       ! Full filename:
       WRITE (fn, '("qbo", i4)') yr
       
       ! Monthly mean zonal wind observations at Equator [m/s] ('QBO'):
       varname = 'QBOEQ'
       varname_longname = 'zonal wind observations at Equator'
       CALL socol_read_netcdf(fn0, varname, 'LEV', llevnoecham=.TRUE., yr=yr, &
            mo=mo, varname_longname=varname_longname, data2d=qboeq_qbolev_m3, &
            nlevdimnoecham=nqbolev, levdimnoecham=qbolev)
       
       ! Level depending half width of QBO [m/s]:
       varname = 'QBOHW'   
       varname_longname = 'level depending half width of QBO'
       CALL socol_read_netcdf(fn, varname, 'LEV', llevnoecham=.TRUE., &
            varname_longname=varname_longname, data1d=qbohw_qbolev)
       
       ! Level depending damping constant of QBO [-]:
       varname = 'QBODAMPV'   
       varname_longname = 'level depending damping constant of QBO'
       CALL socol_read_netcdf(fn, varname, 'LEV', llevnoecham=.TRUE., &
            varname_longname=varname_longname, data1d=qbodampv_qbolev)

       ! Find first and last indices of ECHAM grid in range of QBO grid,
       ! and find for each level of the ECHAM grid the index of the next
       ! following level of the QBO grid:
       CALL socol_interpolate_index(qbolev, presf_stand_socol, &
            k0_qbonudg, k1_qbonudg, qboindex)

    ENDIF

  END SUBROUTINE read_socol_qbo


  SUBROUTINE interpolate_socol_qbo(krow, kproma, kbdim, klev, papm1, papp1)

    ! Interpolates monthly monthly QBO observations to the current time step.
    ! (The interpolatation to the model pressure grid is performed

    ! *interpolate_socol_qbo* is called from *qbonudg_socol*.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim, klev
    REAL(dp), INTENT(in) :: papm1(kbdim,klev), papp1(kbdim,klev)

    ! Local variables:
    INTEGER  :: jl, jk
    INTEGER  :: iind0(k0_qbonudg:k1_qbonudg), iind1(k0_qbonudg:k1_qbonudg)
    REAL(dp) :: pfp(k0_qbonudg:k1_qbonudg), w0(k0_qbonudg:k1_qbonudg), &
                w1(k0_qbonudg:k1_qbonudg), dwinv(k0_qbonudg:k1_qbonudg)
    REAL(dp) :: qboeq_qbolev(nqbolev)

    ! Executable statements:

    ! 1. Allocate memory:
    IF (.NOT. ALLOCATED(qboeq)) &
         ALLOCATE(qboeq(kbdim,k0_qbonudg:k1_qbonudg))
    IF (.NOT. ALLOCATED(qbohw)) &
         ALLOCATE(qbohw(kbdim,k0_qbonudg:k1_qbonudg))
    IF (.NOT. ALLOCATED(qbodampv)) &
         ALLOCATE(qbodampv(kbdim,k0_qbonudg:k1_qbonudg))

    ! 2. Initialize fields with 0:
    qboeq(1:kproma,:) = 0._dp
    qbohw(1:kproma,:) = 0._dp
    qbodampv(1:kproma,:) = 0._dp
           
    ! 3. Interpolation to current time step:
    qboeq_qbolev(:) = &
         wgt1_current*qboeq_qbolev_m3(:,m3w1_current) + &
         wgt2_current*qboeq_qbolev_m3(:,m3w2_current)

    ! 4. Interpolation to current model pressure grid:

    ! Level indices for vertical interpolation:
    iind0(:) = qboindex(k0_qbonudg:k1_qbonudg)-1
    iind1(:) = qboindex(k0_qbonudg:k1_qbonudg)

    DO jl = 1, kproma

       ! Nudging for gridboxes equatorwards of qbolat2 (otherwise no nudging):
       IF (ABS(philat_2d(jl,krow)) .LE. qbolat2) THEN

          ! Pressure at full levels at time t [hPa]:
          pfp(:) = 0.005*(papm1(jl,k0_qbonudg:k1_qbonudg)+ &
               papp1(jl,k0_qbonudg:k1_qbonudg))

          ! Weights:
          w0(:) = qbolev(iind1(:))-pfp(:)
          w1(:) = pfp(:)-qbolev(iind0(:))
          dwinv(:) = 1./(qbolev(iind1(:))-qbolev(iind0(:)))

          ! Interpolation:
          qboeq(jl,:) = (w0(:)*qboeq_qbolev(iind0(:))+ &
               w1(:)*qboeq_qbolev(iind1(:))) * dwinv(:)
          qbohw(jl,:) = (w0(:)*qbohw_qbolev(iind0(:))+ &
               w1(:)*qbohw_qbolev(iind1(:))) * dwinv(:)
          qbodampv(jl,:) = (w0(:)*qbodampv_qbolev(iind0(:)) + &
               w1(:)*qbodampv_qbolev(iind1(:))) * dwinv(:)
       ENDIF
    ENDDO

  END SUBROUTINE interpolate_socol_qbo


  SUBROUTINE qbonudg_socol (krow, kproma, kbdim, klev, pum1, pvom, papm1, papp1)

    ! Computes tendency du/dt due to the QBO nudging.
    ! *qbonudg_socol* is called from *physc*.

    ! Subroutine arguments:
    INTEGER, INTENT(in)  :: krow, kproma, kbdim, klev
    REAL(dp), INTENT(in) :: pum1(kbdim,klev)    ! zonal wind
    REAL(dp), INTENT(inout) :: pvom(kbdim,klev) ! tendency of zonal wind
    REAL(dp), INTENT(in) :: papm1(kbdim,klev), papp1(kbdim,klev)  ! pressure

    ! Local variables:
    INTEGER  :: jl
    REAL(dp) :: latabs, latabsq, qbodamph, log5
    REAL(dp) :: qbo(k0_qbonudg:k1_qbonudg), qbodamp(k0_qbonudg:k1_qbonudg), &
                qbohwq(k0_qbonudg:k1_qbonudg), &
                dudt_qbonudg(k0_qbonudg:k1_qbonudg)

    ! Executable statements:

    ! Interpolate prescribed QBO observations to current time step:
    IF (lqbonudg) CALL interpolate_socol_qbo(krow, kproma, kbdim, klev, &
         papm1, papp1)

    log5 = LOG(0.5_dp)
    
    DO jl = 1, kproma
       latabs = ABS(philat_2d(jl,krow))

       ! Nudging for gridboxes equatorwards of qbolat2 (otherwise no nudging):
       IF (latabs .LE. qbolat2) THEN 
          latabsq = latabs*latabs
          qbohwq(:) = qbohw(jl,:)*qbohw(jl,:)

          ! Compute undamped, latitude and pressure depending QBO:
          qbo(:) = qboeq(jl,:)*EXP(log5*latabsq/qbohwq(:))

          ! Compute latitude depending damping constant:
          qbodamph = MIN((qbolat2-latabs)/(qbolat2-qbolat1), 1._dp)

          ! Compute latitude and pressure depending nudging strength 
          ! (damping factor):
          qbodamp(:) = qbodamp0*qbodampv(jl,:)*qbodamph

          ! Compute tendency du/dt due to QBO nudging:     
          dudt_qbonudg(:) = (qbo(:)-(pum1(jl,k0_qbonudg:k1_qbonudg)+ &
               time_step_len*pvom(jl,k0_qbonudg:k1_qbonudg))) / &
               (1._dp+time_step_len*qbodamp(:)) * qbodamp(:)

          ! Add dudt_qbonudg to pvom:
          pvom(jl,k0_qbonudg:k1_qbonudg) = pvom(jl,k0_qbonudg:k1_qbonudg) + &
               dudt_qbonudg(:)
       ENDIF
    ENDDO

  END SUBROUTINE qbonudg_socol


  SUBROUTINE cleanup_socol_qbo

    ! Deallocates module variables.
    
    ! *cleanup_socol_qbo* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(qbolev)) DEALLOCATE(qbolev)
    IF (ALLOCATED(qboindex)) DEALLOCATE(qboindex)
    IF (ALLOCATED(qbohw_qbolev)) DEALLOCATE(qbohw_qbolev)
    IF (ALLOCATED(qbohw)) DEALLOCATE(qbohw)
    IF (ALLOCATED(qbodampv_qbolev)) DEALLOCATE(qbodampv_qbolev)
    IF (ALLOCATED(qbodampv)) DEALLOCATE(qbodampv)
    IF (ALLOCATED(qboeq_qbolev_m3)) DEALLOCATE(qboeq_qbolev_m3)
    IF (ALLOCATED(qboeq)) DEALLOCATE(qboeq)

  END SUBROUTINE cleanup_socol_qbo

END MODULE mo_socol_qbo

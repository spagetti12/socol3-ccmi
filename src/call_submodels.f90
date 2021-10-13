!
! This file holds the calls to the entries of the submodels 
! attached to ECHAM.
!
! Change this file to attach submodels to ECHAM.
!
!==============================================================================
  SUBROUTINE call_init_submodels
  !
  ! This routine is called from 'initialize'
  !
  ! new modules must be introduced from here (call new_submodel)
  !
    USE mo_control,            ONLY: lsocol                           !MSSOCOL
    USE mo_socol_namelist,     ONLY: init_socol, lchem, lo3orig                !MSSOCOL
    USE mo_socol_time_control, ONLY: init_chem_events                 !MSSOCOL
    USE mo_socol_o3orig,       ONLY: o3orig_initialize                !O3 transport

    IMPLICIT NONE

    IF (lsocol) THEN                                                  !MSSOCOL
       IF (lchem) THEN                                                !MSSOCOL
          ! Read namelist:                                            !MSSOCOL
          CALL init_socol                                             !MSSOCOL
          IF(lo3orig) CALL o3orig_initialize                          !O3 transport
          
          ! Initialize events to call chemistry module:
          CALL init_chem_events                                       !MSSOCOL

       ENDIF                                                          !MSSOCOL
    ENDIF                                                             !MSSOCOL

  END SUBROUTINE call_init_submodels 
!==============================================================================
  SUBROUTINE call_request_tracer
  !
  ! This routine is called from 'initrac', module 'mo_tracer'
  ! 
  ! tracer request routines are called from here
  !
  ! Template for a tracer request routine:
  !
  !    SUBROUTINE request_tracer
  !      USE mo_tracer, ONLY: new_tracer
  !      IMPLICIT NONE
  !      CALL new_tracer('tracername1' ,'modulename', ...)
  !      CALL new_tracer('tracername2' ,'modulename', ...)
  !    END SUBROUTINE request_tracer
  !
    USE mo_control,         ONLY: lsocol                              !MSSOCOL
    USE mo_socol_namelist,  ONLY: lchem, lch4_isotope, &              !MSSOCOL ! eth_as_ch4
         lscav, lo3orig, &                                           ! eth_as_scav
         lsynth                                                       !CCMI tracers
    USE mo_socol_tracers,   ONLY: request_socol_tracers               !MSSOCOL
    USE mo_socol_o3orig,    ONLY: o3orig_new_tracer               !O3 transport
    USE mo_socol_isotope,   ONLY: request_socol_isotope_tracers       ! eth_as_ch4
    USE mo_sub_nml,         ONLY: request_tracer_nml, set_tracer_nml
    USE mo_socol_scav,      ONLY: init_scav                           ! eth_as_scav
    USE mo_socol_synth_tracers, ONLY: request_synth_trac              !CCMI tracers

    IMPLICIT NONE

    ! Call subroutine to request chemical species (tracers) for MEZON:!MSSOCOL
    IF (lsocol) THEN                                                  !MSSOCOL
       IF (lchem) then                                                !MSSOCOL

          CALL request_socol_tracers                                  !MSSOCOL

          IF (lo3orig)  CALL o3orig_new_tracer                        !O3 transport
          IF (lch4_isotope) CALL request_socol_isotope_tracers        ! eth_as_ch4
          IF (lscav) CALL init_scav                                   ! eth_as_scav
          IF (lsynth) CALL request_synth_trac                         !CCMI tracers

       END IF                                                         !MSSOCOL
    ENDIF                                                             !MSSOCOL

    !
    ! routine to request tracers via namelist group /NEW_TRACER/
    !
    CALL request_tracer_nml
    CALL set_tracer_nml
  END SUBROUTINE call_request_tracer
!==============================================================================
  SUBROUTINE call_init_submodel_memory
  !
  ! This routine is called from subroutine 'init_memory', module
  ! 'mo_memory_streams'. Routines are called from here to allocate memory
  ! and to define output streams.
  !
    USE mo_control,         ONLY: lsocol                              !MSSOCOL
    USE mo_socol_gcmfields, ONLY: allocate_socol_gcmfields            !MSSOCOL
    USE mo_socol_namelist,  ONLY: lchem, lscav, lo3orig, lnudg_2h     !MSSOCOL
    USE mo_socol_no_gcm_condens_pscs, ONLY: allocate_no_gcm_condens_pscs  !MSSO 
    USE mo_socol_streams,   ONLY: construct_stream_chem_m, &
                                  construct_stream_chem2, &              !MSSOCOL
                                  construct_stream_chem2d,&
                                  construct_stream_nudg_2h
    USE mo_socol_sun,       ONLY: allocate_socol_photolysis           !MSSOCOL
    USE mo_socol_ch4_streams, ONLY: construct_stream_ch4              ! eth_as_ch4
    USE mo_socol_ch4,       ONLY: allocate_ch4_clim                   ! eth_as_ch4 
    USE mo_socol_scav,      ONLY: construct_stream_scav               ! eth_as_scav
    USE mo_socol_photo,     ONLY: construct_stream_photo              ! eth_as_photo
    USE mo_lpj_streams,     ONLY: construct_stream_lpj                ! eth_as_lpj
    USE mo_socol_o3orig,    ONLY: o3orig_init_memory                  !O3 transport

    IMPLICIT NONE

    IF (lsocol) THEN                                                  !MSSOCOL
       CALL construct_stream_lpj      ! eth_as_lpj
       IF (lchem) THEN                                                !MSSOCOL
          CALL construct_stream_chem_m   ! monthly mean output of MEZON!MSSOCOL
          CALL construct_stream_chem2    ! 12 h output of MEZON       !MSSOCOL
          CALL construct_stream_chem2d
          IF (lnudg_2h) CALL construct_stream_nudg_2h    
          CALL allocate_socol_photolysis ! photolysis rates           !MSSOCOL
          CALL allocate_socol_gcmfields  ! fields from GCM used in MEZON!MSSOCOL
          CALL allocate_no_gcm_condens_pscs ! q, xite, xlte before condensation
          CALL construct_stream_ch4      ! eth_as_ch4
          IF (lscav) CALL construct_stream_scav     ! eth_as_scav
          CALL construct_stream_photo    ! eth_as_photo
          CALL allocate_ch4_clim         ! eth_as_ch4
          IF (lo3orig) CALL o3orig_init_memory        !O3 transport
       ENDIF                                                          !MSSOCOL
    ENDIF                                                             !MSSOCOL

  END SUBROUTINE call_init_submodel_memory
!------------------------------------------------------------------------------
  SUBROUTINE call_free_submodel_memory
  !
  ! Routines are called from here to deallocate memory
  !

    USE mo_control,                 ONLY: lsocol                      !MSSOCOL
    USE mo_socol_gcmfields,         ONLY: cleanup_socol_gcmfields     !MSSOCOL
    USE mo_socol_grid_calculations, ONLY: cleanup_socol_grid_calc     !MSSOCOL
    USE mo_socol_namelist,          ONLY: lchem, lscav, lo3orig,&
                                          lsynth, lnudg_2h            !MSSOCOL
    USE mo_socol_no_gcm_condens_pscs, ONLY: cleanup_no_gcm_condens_pscs!MSSOCOL
    USE mo_socol_streams,           ONLY: destruct_stream_chem_m, &
                                          destruct_stream_chem2, &       !MSSOCOL
                                          destruct_stream_chem2d,&
                                          destruct_stream_nudg_2h
    USE mo_socol_sun,               ONLY: cleanup_socol_photolysis    !MSSOCOL
    USE mo_socol_tracers,           ONLY: destruct_socol_tracers      !MSSOCOL
    USE mo_socol_ch4_streams,       ONLY: destruct_stream_ch4         ! eth_as_ch4  
    USE mo_socol_ch4,               ONLY: cleanup_ch4_clim            ! eth_as_ch4 
    USE mo_socol_scav,              ONLY: destruct_stream_scav        ! eth_as_scav  
    USE mo_socol_photo,             ONLY: destruct_stream_photo       ! eth_as_photo  
    USE mo_lpj_streams,             ONLY: destruct_stream_lpj         ! eth_as_lpj 
    USE mo_socol_o3orig,            ONLY: o3orig_free_memory          !O3 transport
    USE mo_socol_synth_tracers,     ONLY: cleanup_synth_trac          !CCMI tracers

    IMPLICIT NONE

    !  External subroutine:                                           !MSSOCOL
    EXTERNAL :: cleanup_socol_bcond                                   !MSSOCOL
                                                                      !MSSOCOL
    IF (lsocol) THEN                                                  !MSSOCOL
       CALL cleanup_socol_grid_calc     ! grid related quantities     !MSSOCOL
       CALL cleanup_socol_bcond         ! boundary conditions         !MSSOCOL
       CALL destruct_stream_lpj      ! eth_as_lpj
       IF (lchem) THEN                                                !MSSOCOL
          CALL destruct_stream_chem_m   ! monthly mean output of MEZON!MSSOCOL
          CALL destruct_stream_chem2    ! 12 h output of MEZON        !MSSOCOL
          CALL destruct_stream_chem2d 
          IF (lnudg_2h) CALL destruct_stream_nudg_2h                
          CALL destruct_socol_tracers                                 !MSSOCOL
          CALL cleanup_socol_photolysis ! photolysis rates            !MSSOCOL
          CALL cleanup_socol_gcmfields  ! fields from GCM used in MEZON!MSSOCOL
          CALL cleanup_no_gcm_condens_pscs ! q, xite, xlte before condens!MSSOC
          CALL destruct_stream_ch4      ! eth_as_ch4
          IF (lscav) CALL destruct_stream_scav     ! eth_as_scav
          CALL destruct_stream_photo    ! eth_as_photo
          CALL cleanup_ch4_clim         ! eth_as_ch4
          IF(lo3orig) CALL o3orig_free_memory !O3transport
          IF (lsynth) CALL cleanup_synth_trac !CCMI tracers
       ENDIF                                                          !MSSOCOL
    ENDIF                                                             !MSSOCOL

  END SUBROUTINE call_free_submodel_memory
!==============================================================================
  SUBROUTINE call_chem1 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, &
                         ptte, ptm1, pxtm1, pxtte, pqm1, pqte, ktrpwmo, krow, kglat)  !MSSOCOL
  !SUBROUTINE call_chem1 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, &
  !                       ptte, ptm1, pxtm1, pxtte)                   !MSSOCOL

    !
    ! This subroutine calls the different parts of the chemistry
    ! module.
    !
    ! It is called by xtdriver (module mo_tracer) which is called by 'physc'
    !
    USE mo_control,   ONLY: lsocol                                    !MSSOCOL
    USE mo_kind,      ONLY : dp
    USE mo_sub_echam, ONLY : radionucl_sink

    IMPLICIT NONE

    ! External subroutines:                                           !MSSOCOL
    EXTERNAL :: mezon                                                 !MSSOCOL

    ! Scalar arguments
    INTEGER, INTENT(in) :: kproma, kbdim, klev, klevp1, ktrac, krow, kglat   !MSSOCOL

    ! Array arguments
    INTEGER, INTENT(in) :: ktrpwmo(kbdim)                             !MSSOCOL
    REAL(dp), INTENT(in) :: paphp1(kbdim,klevp1), papp1(kbdim,klev), &
         ptte(kbdim,klev), ptm1(kbdim,klev), pqm1(kbdim,klev), &
         pqte(kbdim,klev)                                             !MSSOCOL
    REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac), &
         pxtte(kbdim,klev,ktrac)                                      !MSSOCOL

    !! Scalar arguments                                               !MSSOCOL
    !INTEGER :: kproma, kbdim, klev, klevp1, ktrac
    !
    !! Array arguments
    !REAL(dp):: paphp1(kbdim,klevp1), papp1(kbdim,klev), ptte(kbdim,klev),  &
    !           ptm1(kbdim,klev)
    !REAL(dp):: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)


    ! Call MEZON:     !MSSOCOL


    IF (lsocol) THEN                                                  !MSSOCOL
       CALL mezon(kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, ptte, &
            ptm1, pxtm1, pxtte, pqm1, pqte, ktrpwmo, krow, kglat)            !MSSOCOL
    ENDIF                                                             !MSSOCOL

  

    !
    ! commonly used processes calculated by ECHAM:
    !
    CALL radionucl_sink (kproma, kbdim, klev, pxtm1, pxtte)

  

  END SUBROUTINE call_chem1
!------------------------------------------------------------------------------
  SUBROUTINE call_chem2 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, &
                         ptte, ptm1, pxtm1, pxtte, krow, kglat, &              !MSSOCOL
                         ptvm1, &                 ! eth_as_scav
                         zmratep, zfprec, zfevap, zaclc, zmlwc, & ! eth_as_scav
                         xn3depcv, cvdprec, kconbot, totliq) ! eth_as_scav

!  SUBROUTINE call_chem2 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, &
!                         ptte, ptm1, pxtm1, pxtte)
    !
    ! similar to call_chem1, but called from other place within physc
    USE mo_control,        ONLY : lsocol                              !MSSOCOL
    USE mo_kind,           ONLY : dp
    USE mo_sub_echam,      ONLY : radionucl_sink
    USE mo_socol_namelist, ONLY : lchem, lscav                        !MSSOCOL
    USE mo_socol_tracers,  ONLY : pos_socol_tracers                   !MSSOCOL
    USE mo_socol_scav,     ONLY : nspec_scav, socol_scav              ! eth_as_scav   

    IMPLICIT NONE

    ! Scalar arguments
    INTEGER :: kproma, kbdim, klev, klevp1, ktrac, krow, kglat               !MSSOCOL
    ! eth_as_scav+
    INTEGER  :: kconbot(kproma)
    ! eth_as_scav-

    ! Array arguments
    REAL(dp):: paphp1(kbdim,klevp1), papp1(kbdim,klev), ptte(kbdim,klev),  &
               ptm1(kbdim,klev)
    REAL(dp):: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)
    ! eth_as_scav+
    REAL(dp) :: ptvm1(kbdim,klev)
    REAL(dp) :: xn3depcv(kbdim,klev,nspec_scav), cvdprec(kbdim,klev) 
    REAL(dp) :: zmratep(kbdim,klev), zfprec(kbdim,klev), zfevap(kbdim,klev)
    REAL(dp) :: zmlwc(kbdim,klev), zaclc(kbdim,klev) 
    REAL(dp) :: totliq(kbdim,klev)
    ! eth_as_scav-


    IF (lsocol) THEN                                                  !MSSOCOL
       IF (lchem) THEN
          ! eth_as_scav
          IF (lscav) CALL socol_scav(krow, kproma, kbdim, klev, ktrac, papp1, paphp1, ptvm1, pxtm1,&
              pxtte, ptm1, zaclc, &
              zmratep, zfprec, zfevap, zmlwc, xn3depcv, cvdprec, kconbot, totliq)
          !
          CALL pos_socol_tracers(kproma, kbdim, klev, ktrac, &
               krow, pxtm1, pxtte)
       END IF
    ENDIF                                                             !MSSOCOL

  END SUBROUTINE call_chem2
!==============================================================================
  SUBROUTINE call_diagn (kproma, kbdim, klev, ktrac, krow, pxtm1, pxtte)
  !SUBROUTINE call_diagn                                              !MSSOCOL
  !
  ! This routine is called by xtdiagn (mo_tracer) called at the end of physc
  !  physc. Diagnostics at the end of a time step may be called from here
  !
    USE mo_control,        ONLY : lsocol                              !MSSOCOL
    USE mo_kind,           ONLY : dp
    USE mo_socol_namelist, ONLY : lchem, lscav                        !MSSOCOL
    USE mo_socol_streams,  ONLY : accumulate_stream_chem_m,  accumulate_stream_chem2           !MSSOCOL
    USE mo_socol_ch4_streams,  ONLY : accumulate_stream_ch4           ! eth_as_ch4
    USE mo_socol_scav,     ONLY : accumulate_stream_scav              ! eth_as_scav
    USE mo_socol_photo,    ONLY : accumulate_stream_photo             ! eth_as_photo
    USE mo_lpj_streams,    ONLY : accumulate_stream_lpj               ! eth_as_lpj


    IMPLICIT NONE

    ! Scalar arguments
    INTEGER, INTENT(in) :: kproma, kbdim, klev, ktrac, krow           !MSSOCOL

    ! Array arguments                                                 !MSSOCOL
    REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)


    IF (lsocol) THEN                                                  !MSSOCOL
       CALL accumulate_stream_lpj (kproma, krow)                   ! eth_as_lpj
       IF (lchem) THEN                                                !MSSOCOL

!!$          ! Calculate chemical families at time t for output:         !MSSOCOL
!!$          CALL socol_families_chem (kproma, kbdim, klev, ktrac, krow, &
!!$               pxtm1, pxtte)                                          !MSSOCOL

          ! Accumulate values for monthly mean output:                !MSSOCOL
          CALL accumulate_stream_chem_m (kproma, kbdim, klev, ktrac, krow, &
               pxtm1, pxtte)                                          !MSSOCOL

          CALL accumulate_stream_chem2(kproma, kbdim, klev, ktrac, krow)

          CALL accumulate_stream_ch4 (kproma, krow)                   ! eth_as_ch4

          IF (lscav) CALL accumulate_stream_scav (kproma, krow)       ! eth_as_scav

          CALL accumulate_stream_photo (kproma, krow)                 ! eth_as_photo

       ENDIF
    ENDIF                                                             !MSSOCOL

 END SUBROUTINE call_diagn
!==============================================================================
  SUBROUTINE call_init_tracers
  !  
  ! This routine is called from xtini to allow initialization of tracer
  ! variables besides initialization with constant values or from the
  ! rerun file.
  !
  ! Base the decision whether to initialize on the following conditions:
  !
  ! trlist% ti(jt)% init   == 0       (not initialized so far)
  ! trlist% ti(jt)% ninit: >= INITIAL (initialisation requested)
  !
  ! Set trlist%ti(jt) to a value =/0 afterwards.
  !
    USE mo_control,        ONLY : lsocol                              !MSSOCOL
    USE mo_socol_namelist, ONLY : lchem, lch4_isotope, &              !MSSOCOL ! eth_as_ch4
                                  lscav, &                            ! eth_as_scav
                                  lsynth                              !CCMI tracers
    USE mo_socol_streams,  ONLY : init_stream_chem_m, &               !MSSOCOL
                                  init_stream_chem2, &                   !MSSOCOL
                                  init_stream_chem2d
    USE mo_socol_tracers,  ONLY : init_socol_tracers                  !MSSOCOL
    USE mo_socol_ch4_streams,  ONLY : init_stream_ch4                 ! eth_as_ch4
    USE mo_socol_scav,     ONLY : init_stream_scav                    ! eth_as_scav
    USE mo_socol_isotope,  ONLY : init_socol_isotope_tracers          ! eth_as_ch4
    USE mo_socol_photo,    ONLY : init_stream_photo                   ! eth_as_photo
    USE mo_lpj_streams,    ONLY : init_stream_lpj                     ! eth_as_lpj
    USE mo_socol_synth_tracers,     ONLY: init_synth_trac             !CCMI tracers
    IMPLICIT NONE

    IF (lsocol) THEN                                                  !MSSOCOL
       CALL init_stream_lpj                                        ! eth_as_lpj
       IF (lchem) THEN                                                !MSSOCOL
          CALL init_socol_tracers                                     !MSSOCOL
          CALL init_stream_chem_m                                     !MSSOCOL
          CALL init_stream_chem2                                      !MSSOCOL
          CALL init_stream_chem2d
          CALL init_stream_ch4                                        ! eth_as_ch4
          CALL init_stream_photo                                      ! eth_as_photo
          IF (lscav) CALL init_stream_scav                            ! eth_as_scav

          IF (lch4_isotope) CALL init_socol_isotope_tracers           ! eth_as_ch4
          IF (lsynth) CALL init_synth_trac                            !CCMI tracers

       ENDIF                                                          !MSSOCOL
    ENDIF                                                             !MSSOCOL
      
  END SUBROUTINE call_init_tracers
!==============================================================================
  SUBROUTINE call_read_bcond
  !
  ! This routine is called from 'stepon'. Routines which read forcing
  ! data may be called from here.
  !
    USE mo_control,         ONLY: lsocol                              !MSSOCOL
    USE mo_socol_namelist,  ONLY: linit_socol_bcond                   !MSSOCOL
    USE mo_time_control,    ONLY: current_date, get_date_components, &
                                  next_date                           !MSSOCOL

    IMPLICIT NONE

    !INTEGER :: imonth
    INTEGER :: icurrentyear, inextyear, icurrentmonth,  inextmonth, inextday, icurrentday
    LOGICAL :: lnewmonth, lnewyear, lnewday

    EXTERNAL read_socol_bcond_y, read_socol_bcond_m, &                !MSSOCOL
         read_socol_bcond_d, interpolate_socol_bcond_global                           !MSSOCOL

    !CALL get_date_components (current_date, month=imonth)            !MSSOCOL
    CALL get_date_components (current_date, month=icurrentmonth, &
         year=icurrentyear, day=icurrentday)                                           !MSSOCOL
    CALL get_date_components (next_date,  month=inextmonth, &
         year=inextyear, day=inextday)                                              !MSSOCOL

    lnewyear = icurrentyear .NE. inextyear                            !MSSOCOL
    lnewmonth = icurrentmonth .NE. inextmonth                         !MSSOCOL
    lnewday = icurrentday .NE. inextday

    ! Read boundary conditions for SOCOL:                             !MSSOCOL

    IF (lsocol) THEN                                                  !MSSOCOL
       ! a) Initialization: boundary conditions to be read every year:
       IF (linit_socol_bcond .OR. lnewyear) &
            CALL read_socol_bcond_y(inextyear)                        !MSSOCOL
       
       ! b) Initialization: boundary conditions to be read every month:
       IF (linit_socol_bcond .OR. lnewmonth) &
            CALL read_socol_bcond_m(inextyear,inextmonth, inextday)             !MSSOCOL

       ! c) Initialization: boundary conditions to be read every day:
       IF (linit_socol_bcond .OR. lnewday) THEN
          CALL read_socol_bcond_d(inextyear,inextmonth, inextday)
       ENDIF

       linit_socol_bcond = .FALSE.                                    !MSSOCOL

       ! c) Interpolation to the current time step (for 1d-fields):   !MSSOCOL
       CALL interpolate_socol_bcond_global                            !MSSOCOL

       ! d) Interpolation to the current time step (for 3d-fields):   !MSSOCOL
       ! Call from *physc*.                                           !MSSOCOL
    ENDIF                                                             !MSSOCOL

  END SUBROUTINE call_read_bcond
!==============================================================================
  SUBROUTINE call_chem_bcond (kproma, kbdim, klev, pxtte)

    !
    ! This is the interface routine between ECHAM's standard tracer treatment
    ! and the chemistry routines.
    !
    ! This routine is called from vdiff
    !
    ! M. Schultz and Hans-Stefan Bauer, MPI Hamburg, November 2000
    !

    USE mo_kind,            ONLY: dp
    USE mo_tracer,          ONLY: ntrac
    IMPLICIT NONE

    ! Arguments

    INTEGER, INTENT(in)    :: kproma, kbdim, klev
    REAL(dp),INTENT(inout) :: pxtte(kbdim,klev,ntrac)

  END SUBROUTINE call_chem_bcond
!==============================================================================

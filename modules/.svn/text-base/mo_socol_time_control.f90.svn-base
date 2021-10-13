MODULE mo_socol_time_control

  ! Time control for SOCOL. This module is similar to *mo_time_control.f90*.

  ! Martin Schraner, ETH Zurich, December 2008

  USE mo_exception,       ONLY: finish, message
  USE mo_kind,            ONLY: dp
  USE mo_socol_interpo,   ONLY: wgt1_current, wgt2_current, &
                                wgt1_chem, wgt2_chem, wgt1_rad, wgt2_rad, &
                                m3w1_current, m3w2_current, &
                                m3w1_chem, m3w2_chem, m3w1_rad, m3w2_rad, &
                                yw1_current, yw2_current, &
                                yw1_chem, yw2_chem, yw1_rad, yw2_rad
  USE mo_socol_namelist,  ONLY: delta_time_chem, lchem
  USE mo_time_conversion, ONLY: time_days, add_date, time_native, month_len
  USE mo_time_event,      ONLY: time_event, event_state, io_time_event, &
                                TRIG_FIRST, TIME_INC_SECONDS
  USE mo_time_control,    ONLY: delta_time, echam_ev_init, echam_ev_print, &
                                ev_trigrad, l_trigrad, ldebugev, &
                                current_date, next_date, write_date, lstart, &
                                NDAYLEN, radiation_date, get_date_components, &
                                get_year_day, get_year_len, tc_convert, &
                                tc_set, tc_get

  IMPLICIT NONE

  PRIVATE

  TYPE(time_days), PUBLIC, SAVE  :: chemistry_date      ! date corresponding to
                                                        ! call MEZON
  TYPE(time_event), PUBLIC, SAVE :: ev_trigchem
  TYPE(io_time_event), PUBLIC, SAVE :: trigchem
  LOGICAL, PUBLIC                :: l_trigchem = .TRUE.

  INTEGER, PUBLIC  :: chemistry_date_yr      ! year of chemistry date
  INTEGER, PUBLIC  :: chemistry_date_mo      ! month of chemistry date
  INTEGER, PUBLIC  :: chemistry_date_dy      ! days of chemistry date
  INTEGER, PUBLIC  :: chemistry_date_hr      ! hours of chemistry date
  INTEGER, PUBLIC  :: chemistry_date_mn      ! minutes of chemistry date
  INTEGER, PUBLIC  :: chemistry_date_se      ! seconds of chemistry date
  REAL(dp), PUBLIC :: chemistry_date_dyofyr  ! chemistry_date as number of days
                                             ! of the current year
  REAL(dp), PUBLIC :: chemistry_date_daypart ! hours, minutes, seconds as part 
                                             ! of the current day
  REAL(dp), PUBLIC :: ndays_currentyear      ! Number of days of the current 
                                             ! year (including part of current
                                             ! day)
  INTEGER, PUBLIC  :: m3w1, m3w2             ! weight indices for interpolation
                                             ! of SOCOL boundary condition data
                                             ! set of current, preceding and 
                                             ! following month
  REAL(dp), PUBLIC :: time_step_len_chem     ! integration time length of 
                                             ! chemistry step

  PUBLIC :: init_chem_events, time_set_socol

CONTAINS
    
  SUBROUTINE init_chem_events

    ! Initializes event to call chemistry module.

    ! *init_chem_events* is called from *call_init_submodels*, 
    ! src/call_submodels.f90.

    CHARACTER(len=256)          :: m_text
    CHARACTER(len=*), PARAMETER :: EV_TLEV_PRES = 'present' ! check event 
                                                       ! with present date
    
    ! Chemistry module has to be called at least every model time step:
    IF (delta_time_chem .LT. delta_time) delta_time_chem = delta_time

    ! Set time event:
    trigchem = io_time_event (delta_time_chem,TIME_INC_SECONDS,TRIG_FIRST,0) 

    ! Calculate intitial state of event dependent on the initial date 
    ! of a run (at rerun the event triggers will be recalculated from 
    ! beginning):
    CALL echam_ev_init &
         (ev_trigchem,trigchem,'chemistry computation', EV_TLEV_PRES)

    ! Print settings in debugging mode (as in *construct_events*,
    ! modules/mo_time_control.f90):
    IF (ldebugev) THEN
       CALL message('',''); CALL message('','&SOCOLCTL')
       WRITE(m_text,*) 'TRIGCHEM=    ',trigchem,   ','; CALL message('',m_text)
    ENDIF

    ! Set integration time length of chemistry step:
    time_step_len_chem = 2.0_dp*delta_time_chem 
       
  END SUBROUTINE init_chem_events


  SUBROUTINE time_set_socol

    ! Sets time events for calling MEZON. Besides sets time weights for
    ! SOCOL boundary conditions and saves them in *mo_socol_interpo*.

    ! *time_set_socol* is called from *stepon*, immediatly after call
    ! of *time_set*.

    INTEGER :: yr, mo, dy, hr, mn, se, ichemhlen
    TYPE (time_days)   :: date_time_days
    TYPE (time_native) :: date_mon, date_time_native

    IF (delta_time_chem .GT. NINT(delta_time)) THEN ! otherwise keep 
                                                    ! l_trigchem = .TRUE.
                                                    ! and call MEZON at every
                                                    ! timestep
       l_trigchem = event_state(ev_trigchem, current_date) .OR. lstart
    ENDIF

    IF (.NOT. lchem) l_trigchem = .FALSE.

    ! Print settings in debugging mode:
    IF (ldebugev) THEN
       CALL message('time_set','----------------------------------------')
       CALL write_date(current_date,'Current date: ')
       CALL message('','Events triggered with current date (2) ...')
       IF (l_trigchem) CALL echam_ev_print(ev_trigchem)
    ENDIF

    ! Set chemistry date:
    IF (lchem) chemistry_date = current_date
    IF (l_trigchem) THEN
       ichemhlen = INT(0.5_dp*event_state(ev_trigchem,.TRUE.))
       CALL add_date(0,ichemhlen,chemistry_date)
       CALL write_date(chemistry_date,'Chemistry calculated for : ')

       ! Get year, month, hour, minute and seconds of chemistry_date:
       CALL get_date_components(chemistry_date, year=chemistry_date_yr, &
            month=chemistry_date_mo, day=chemistry_date_dy, &
            hour=chemistry_date_hr, minute=chemistry_date_mn, &
            second=chemistry_date_se)

       ! Calculate hours, minutes and seconds as part of the current day:
       chemistry_date_daypart = REAL(chemistry_date_hr,dp)/24._dp + &
            REAL(chemistry_date_mn,dp)/1440._dp + &
            REAL(chemistry_date_se,dp)/86400._dp

       ! Caluclate chemistry_date as number of days of the current year:
       chemistry_date_dyofyr = get_year_day(chemistry_date) 

       ! Number of days of current year:
       ndays_currentyear = get_year_len(chemistry_date_yr)
    ENDIF

    ! Calculate time weights:
    CALL time_weights_socol(current_date, wgt1_current, &
         wgt2_current, m3w1_current, m3w2_current, yw1_current, yw2_current)
    IF (lchem) CALL time_weights_socol(chemistry_date, wgt1_chem, &
         wgt2_chem, m3w1_chem, m3w2_chem, yw1_chem, yw2_chem)
    IF (l_trigrad) CALL time_weights_socol(radiation_date, wgt1_rad, &
         wgt2_rad, m3w1_rad, m3w2_rad, yw1_rad, yw2_rad)

    CONTAINS

      SUBROUTINE time_weights_socol(weighting_date, wgt1_loc, wgt2_loc, &
           m3w1_loc, m3w2_loc, yw1_loc, yw2_loc)

        ! Calculates weighting factors for boundary conditions and saves them
        ! in *mo_socol_interpo*.

        ! *time_weights_socol* is called from *time_set_socol*.

        ! Subroutine arguments:
        TYPE(time_days), INTENT(in) :: weighting_date
        INTEGER, INTENT(out)  :: m3w1_loc, m3w2_loc, yw1_loc, yw2_loc
        REAL(dp), INTENT(out) :: wgt1_loc, wgt2_loc
        
        ! Local variables:
        TYPE (time_native) :: date_monm1, date_mon, date_monp1
        INTEGER   :: yr, mo, dy, hr, mn, se
        INTEGER   :: days, seconds, isec
        INTEGER   :: imp1, imm1, imp1cl, imm1cl, imlenm1, imlen, imlenp1
        REAL (dp) :: zsec, zdayl
        REAL (dp) :: zmohlf, zmohlfp1, zmohlfm1
        REAL (dp) :: zdh, zdhp1, zdhm1
        
        ! Set calendar related parameters:
        CALL TC_convert(weighting_date,date_mon)
        CALL TC_get(date_mon, yr, mo, dy, hr, mn, se)

        ! Month index for boundary condition data of type (0,...,13):
        imp1 = mo+1
        imm1 = mo-1

        ! Month index for cyclic climatological data (1,...,12)
        imp1cl = mo+1
        imm1cl = mo-1
        IF (imp1cl > 12) imp1cl= 1
        IF (imm1cl <  1) imm1cl=12

        ! Determine length of months and position within current month:
        CALL TC_set(yr, imm1cl, 1, 0, 0, 0, date_monm1)
        CALL TC_set(yr, imp1cl, 1, 0, 0, 0, date_monp1)
        imlenm1 = month_len(date_monm1)
        imlen   = month_len(date_mon)
        imlenp1 = month_len(date_monp1)
      
        zdayl    = REAL(NDAYLEN,dp)
        zmohlfm1 = imlenm1*zdayl*0.5_dp
        zmohlf   = imlen  *zdayl*0.5_dp
        zmohlfp1 = imlenp1*zdayl*0.5_dp
        
        ! Weighting factors for first/second half of month:
        yw1_loc  = mo
        m3w1_loc = 1
        
        ! Seconds in the present month:
        CALL TC_get (weighting_date, days, seconds)
        isec = (dy-1)*NDAYLEN + seconds
        zsec = REAL(isec,dp)
        
        IF(zsec <= zmohlf) THEN                     ! first part of month
           wgt1_loc   = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
           wgt2_loc   = 1._dp-wgt1_loc
           yw2_loc    = imm1
        ELSE                                        ! second part of month
           wgt2_loc   = (zsec-zmohlf)/(zmohlf+zmohlfp1)
           wgt1_loc   = 1._dp-wgt2_loc
           yw2_loc    = imp1
        ENDIF
        
        m3w2_loc = yw2_loc-mo+1
        
      END SUBROUTINE time_weights_socol
   
  END SUBROUTINE time_set_socol

END MODULE mo_socol_time_control
    


MODULE mo_lpj_streams

  ! *mo_lpj_streams* contains subroutines to define streams of 
  ! SOCOL output suitable as LPJ input
  !
  ! Andrea Stenke, ETH Zurich, March 2010

  USE mo_kind,          ONLY: dp      
  USE mo_linked_list,   ONLY: HYBRID, NETCDF
  USE mo_memory_base,   ONLY: new_stream, add_stream_element,  &
                              default_stream_setting, add_stream_reference, &
                              delete_stream, t_stream
  USE mo_time_control,  ONLY: delta_time, lstart
  USE mo_time_event,    ONLY: io_time_event, TIME_INC_HOURS, TRIG_FIRST

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_stream_lpj   ! construct stream
  PUBLIC :: destruct_stream_lpj    ! destruct stream
  PUBLIC :: init_stream_lpj        ! initialize stream
  PUBLIC :: accumulate_stream_lpj  ! accumulate stream elements

  TYPE (t_stream), PUBLIC, POINTER :: lpj     !the lpj stream

  REAL(dp), PUBLIC, POINTER :: tsurf_cel(:,:)
  REAL(dp), PUBLIC, POINTER :: tot_precip(:,:)
  REAL(dp), PUBLIC, POINTER :: sundur(:,:)
  REAL(dp), PUBLIC, POINTER :: sundur_d(:,:)
  REAL(dp), PUBLIC, POINTER :: sec_day(:,:)

!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
  SUBROUTINE construct_stream_lpj

    ! Allocates output streams

    ! *construct_stream_lpj* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (lpj ,'lpj', filetype=NETCDF, interval=io_time_event(24,TIME_INC_HOURS,TRIG_FIRST,12) )     
  
  ! add entries for geopotential and log surface pressure (used by the
  ! afterburner):
  CALL add_stream_reference (lpj, 'geosp'   ,'g3b'   ,lpost=.FALSE.)
  CALL add_stream_reference (lpj, 'lsp'     ,'sp'    ,lpost=.FALSE.)

  CALL default_stream_setting (lpj, lrerun    = .TRUE. , &
                                    leveltype = HYBRID , &
                                    table     = 199,     &
                                    laccu     = .TRUE.)

  CALL add_stream_element (lpj, 'tsurf', tsurf_cel, lpost=.FALSE., &
                           longname='surface temperature', &
                           units='celsius', code=100)
  CALL add_stream_element (lpj, 'tot_precip', tot_precip, lpost=.FALSE., &
                           longname='total precipitation', &
                           units='mm', code=101)
  CALL add_stream_element (lpj, 'sundur', sundur, lpost=.FALSE., &
                           laccu = .FALSE., &
                           longname='sunshine duration', &
                           units='%', code=102)
  CALL add_stream_element (lpj, 'sundur_d', sundur_d, lpost=.FALSE., &
                           laccu = .FALSE., &
                           longname='sunshine duration', &
                           units='%', code=103)
  CALL add_stream_element (lpj, 'sec_day', sec_day, lpost=.FALSE., &
                           laccu = .FALSE., &
                           longname='seconds per day', &
                           units='s', code=104)


END SUBROUTINE construct_stream_lpj

!---------------------------------------------------------------------------------

SUBROUTINE destruct_stream_lpj

  ! Deallocates memory.

  ! *destruct_stream_* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(lpj)

END SUBROUTINE destruct_stream_lpj

!---------------------------------------------------------------------------------

SUBROUTINE init_stream_lpj

  ! Initializes streams with zero.

  ! *init_stream_* is called from *call_init_tracers*,
  ! src/call_submodels.f90.

  IF (lstart) THEN     ! use restart fields otherwise
     tsurf_cel(:,:)  = 0._dp
     tot_precip(:,:) = 0._dp
     sundur(:,:)     = 0._dp
     sundur_d(:,:)   = 0._dp
     sec_day(:,:)    = 0._dp
  ENDIF
 
END SUBROUTINE init_stream_lpj
!---------------------------------------------------------------------------------

SUBROUTINE accumulate_stream_lpj(kproma,krow)

  ! *accumulate_stream* is called from *call_diagn*, 
  ! src/call_submodels

  USE mo_time_control, ONLY: current_date, get_date_components

  ! Scalar arguments
  INTEGER, INTENT(in) :: kproma, krow

  INTEGER :: year, month, day, hour, minute, second

  call get_date_components(current_date, year, month, &
            day, hour, minute, second)

  IF (hour == 0 .and. minute == 0 .and. second == 0) THEN
     sundur(1:kproma,krow) = sundur_d(1:kproma,krow)/sec_day(1:kproma,krow) * 100._dp
     sundur_d(1:kproma,krow) = 0._dp
     sec_day(1:kproma,krow)  = 0._dp
  END IF

 
END SUBROUTINE accumulate_stream_lpj


END MODULE mo_lpj_streams

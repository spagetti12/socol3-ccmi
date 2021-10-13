PROGRAM timesocol2hiphop

  ! Umformung von "%Y%m%d" zu "month since 15-Dec-SYEAR"  und Pa -> mb
  ! für NetCdf-File
  ! Input: ncname, syear

  USE netcdf

  IMPLICIT NONE

  INTEGER :: i, nt, nlev, m, y, month, syear, time_int, oldsyear, oldsmonth
  INTEGER :: ncid, status, varid, len, timeid, timevarid, year, months, &
       levid, levvarid
  REAL :: d
  REAL, DIMENSION(:), ALLOCATABLE :: time, lev
  CHARACTER*80 ncname
  CHARACTER*30 oldtimeunits
  CHARACTER*23 newtimeunits
  CHARACTER*11 newtimeorigin
  CHARACTER*4  syear_st, syearm1_st
  LOGICAL :: lecham5, llev

  llev=.TRUE.

  CALL GETARG(1,ncname)
  CALL GETARG(2,syear_st)
  READ(syear_st,'(i4.4)') syear

  ! open netcdf-file:   

  status=NF90_OPEN(TRIM(ncname), NF90_WRITE, ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)


  ! get information about dimensions:
  
  status=NF90_INQ_DIMID(ncid, 'time', timeid)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQUIRE_DIMENSION(ncid, timeid, len = nt)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQ_DIMID(ncid, 'lev', levid)
  IF (status /= NF90_NOERR) llev=.FALSE.
  !IF (status /= NF90_NOERR) CALL handle_err(status)
  IF (llev) THEN
     status=NF90_INQUIRE_DIMENSION(ncid, levid, len = nlev)
     IF (status /= NF90_NOERR) CALL handle_err(status)
  ENDIF

  ALLOCATE(time(nt))
  IF (llev) ALLOCATE(lev(nlev))

  status=NF90_INQ_VARID(ncid, 'time', timevarid)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_GET_VAR(ncid, timevarid, time) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_GET_ATT(ncid, timevarid, 'units', oldtimeunits)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  IF (llev) THEN
     status=NF90_INQ_VARID(ncid, 'lev', levvarid)
     IF (status /= NF90_NOERR) CALL handle_err(status)
     status=NF90_GET_VAR(ncid, levvarid, lev) 
     IF (status /= NF90_NOERR) CALL handle_err(status)
  ENDIF
  IF (oldtimeunits(1:10) .EQ. 'days since' .OR. &
       oldtimeunits(1:9) .EQ. 'day since') THEN
     lecham5=.TRUE.
  ELSE
     lecham5=.FALSE.
  ENDIF

  ! neue Zeit:

  IF (lecham5) THEN  !ECHAM5

     IF (oldtimeunits(1:10) .EQ. 'days since') THEN 
        READ(oldtimeunits(12:15),'(i4.4)') oldsyear
        READ(oldtimeunits(17:18),'(i2.2)') oldsmonth
     ELSE
        READ(oldtimeunits(11:14),'(i4.4)') oldsyear
        READ(oldtimeunits(16:17),'(i2.2)') oldsmonth
     ENDIF

     DO i=1, nt
        time_int=NINT((time(i))/30.)
        time(i) = REAL(time_int)+REAL(oldsmonth)-1.+(oldsyear-syear)*12.
     ENDDO

  ELSE  !ECHAM4

     DO i=1, nt
        y=time(i)/10000   !Jahr
        m=time(i)/100-100*y   !Monat
        d=time(i)-10000.*y-100.*m    !Tag
        
        time(i)=(y-1)*12. + m-1. + (d-1.)/30.
     ENDDO

  END IF

  WRITE(syearm1_st,'(i4)') syear-1

  newtimeunits(1:19)= 'month since 15-DEC-'
  newtimeunits(20:23)=syearm1_st
  newtimeorigin(1:7)='15-DEC-'
  newtimeorigin(8:11)=syearm1_st


  ! neue Levels:
  IF (llev) lev(:)=0.01*lev(:)

  ! Oeffne Def-Mode und definiere Zeit- und Level-Attribute neu:
  status=NF90_REDEF(ncid) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, timevarid, 'long_name', 'model time')
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_DEL_ATT(ncid, timevarid, 'units') 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, timevarid, 'units', newtimeunits)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, timevarid, 'time_origin', newtimeorigin) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  IF (llev) THEN
     status=NF90_DEL_ATT(ncid, levvarid, 'units')
     IF (status /= NF90_NOERR) CALL handle_err(status)
     status=NF90_PUT_ATT(ncid, levvarid, 'units', 'mb')
     IF (status /= NF90_NOERR) CALL handle_err(status)
     status=NF90_PUT_ATT(ncid, levvarid, 'positive', 'down')
     IF (status /= NF90_NOERR) CALL handle_err(status)
  ENDIF
  status=NF90_ENDDEF(ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Oeffne Data-Mode und definiere die Zeit und Levels neu:
  status=NF90_PUT_VAR(ncid, timevarid, time) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  IF (llev) THEN
     status=NF90_PUT_VAR(ncid, levvarid, lev) 
     IF (status /= NF90_NOERR) CALL handle_err(status)
  ENDIF
  status = NF90_CLOSE(ncid)
  IF (status /= NF90_NOERR) CALL HANDLE_ERR(status)

  CONTAINS

    SUBROUTINE HANDLE_ERR(status)
  
      INTEGER, INTENT (IN) :: status

      IF (status /= NF90_NOERR) THEN
         PRINT *, TRIM(NF90_STRERROR(status))
         STOP 'Stopped'
      END IF

    END SUBROUTINE HANDLE_ERR

END PROGRAM timesocol2hiphop

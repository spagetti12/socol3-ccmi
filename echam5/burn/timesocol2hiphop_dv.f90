PROGRAM timesocol2hiphop_dv

  ! Umformung von "%Y%m%d" zu "month since 15-Dec-SYEAR"  und Pa -> mb
  ! für NetCdf-File
  ! Input: ncname, syear

  USE netcdf

  IMPLICIT NONE

  INTEGER :: i, nt, nlev, m, y, month, syear, smonth
  INTEGER :: ncid, status, varid, len, timeid, timevarid, year, months, &
       levid, levvarid
  REAL :: d
  REAL, DIMENSION(:), ALLOCATABLE :: time, lev
  CHARACTER*80 ncname
  CHARACTER*30 oldtimeunits
  CHARACTER*27 newtimeunits
  CHARACTER*17 newtimeorigin
  CHARACTER*4  syear_st, syearm1_st
  CHARACTER*3  smonth_st
  LOGICAL :: lecham5

  CALL GETARG(1,ncname) 

  ! open netcdf-file:   

  status=NF90_OPEN(ncname, NF90_WRITE, ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)


  ! get information about dimensions:
  
  status=NF90_INQ_DIMID(ncid, 'time', timeid)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQUIRE_DIMENSION(ncid, timeid, len = nt)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQ_DIMID(ncid, 'lev', levid)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQUIRE_DIMENSION(ncid, levid, len = nlev)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ALLOCATE(time(nt))
  ALLOCATE(lev(nlev))

  status=NF90_INQ_VARID(ncid, 'time', timevarid)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_GET_VAR(ncid, timevarid, time) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_GET_ATT(ncid, timevarid, 'units', oldtimeunits)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQ_VARID(ncid, 'lev', levvarid)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_GET_VAR(ncid, levvarid, lev) 
  IF (status /= NF90_NOERR) CALL handle_err(status)

  IF (oldtimeunits(1:10) .EQ. 'days since') THEN
     lecham5=.TRUE.
  ELSE
     lecham5=.FALSE.
  ENDIF


  ! neue Zeit:

  IF (lecham5) THEN  !ECHAM5

     
     READ(oldtimeunits(17:18),'(i2.2)') smonth  !smonth
     
     SELECT CASE (smonth)
        CASE (1)
           smonth_st='JAN'
        CASE (2)
           smonth_st='FEB'
        CASE (3)
           smonth_st='MAR'
        CASE (4)
           smonth_st='APR'
        CASE (5)
           smonth_st='MAY'
        CASE (6)
           smonth_st='JUN'
        CASE (7)
           smonth_st='JUL'
        CASE (8)
           smonth_st='AUG'
        CASE (9)
           smonth_st='SEP'
        CASE (10)
           smonth_st='OCT'
        CASE (11)
           smonth_st='NOV'
        CASE (12)
           smonth_st='DEC'
     END SELECT

     newtimeorigin(1:2) = oldtimeunits(20:21)    !sday
     newtimeorigin(3:3) = '-'
     newtimeorigin(4:6) = smonth_st
     newtimeorigin(7:7) = '-'
     newtimeorigin(8:11) = oldtimeunits(12:15)   !syear
     newtimeorigin(12:12) = ' '
     newtimeorigin(13:17) = oldtimeunits(23:27)  !shour

  ELSE

     CALL GETARG(2,syear_st)

     DO i=1, nt
        y=time(i)/10000   !Jahr
        m=time(i)/100-100*y   !Monat
        d=time(i)-10000.*y-100.*m-0.01041666633    !Tag inkl. Bruchteil (mit Korrektur)
        
        time(i)=(y-1)*360. + (m-1.)*12. + d-1.
     ENDDO

     newtimeorigin(1:7)  = '01-JAN-'
     newtimeorigin(8:11) = syear_st
     newtimeorigin(12:17)= ' 00:00'
  ENDIF

  newtimeunits(1:10)= 'day since '
  newtimeunits(11:27)= newtimeorigin(1:17)


  ! neue Levels:

  lev(:)=0.01*lev(:)

  ! Oeffne Def-Mode und definiere Zeit- und Level-Attribute neu:
  status=NF90_REDEF(ncid) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, timevarid, 'long_name', 'model time')
  IF (status /= NF90_NOERR) CALL handle_err(status)
  IF (NOT(lecham5)) THEN
     status=NF90_PUT_ATT(ncid, timevarid, 'calendar', '360_day')
     IF (status /= NF90_NOERR) CALL handle_err(status)
  ENDIF
  status=NF90_DEL_ATT(ncid, timevarid, 'units') 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, timevarid, 'units', newtimeunits)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, timevarid, 'time_origin', newtimeorigin) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_DEL_ATT(ncid, levvarid, 'units')
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, levvarid, 'units', 'mb')
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_ATT(ncid, levvarid, 'positive', 'down')
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_ENDDEF(ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Oeffne Data-Mode und definiere die Zeit und Levels neu:

  status=NF90_PUT_VAR(ncid, timevarid, time) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_PUT_VAR(ncid, levvarid, lev) 
  IF (status /= NF90_NOERR) CALL handle_err(status)

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


END PROGRAM timesocol2hiphop_dv

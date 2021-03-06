PROGRAM echamnetcdf_afterburner

  ! Creates netcdf file with transformed pressure and time grid:
  ! (a) Pressure grid: interpolation from sigma model levels (depending
  !     on surface pressure) to pressure grid defined in a namelist
  ! (b) Time units: given as "month since start year"
  !
  ! Input: 
  !
  ! $1 : Input file name (netcdf file to be transformed)
  ! $2 : Output file name (newly created netcdf file)
  ! $3 : File name of namelist containing pressure levels of output file
  ! $4 : Start year for Netcdf files (uniform for all output files)
  !
  ! M. Schraner, ETH Zurich, 10.3.2009

  USE netcdf

  IMPLICIT NONE

  INTEGER :: i, ji, jl, jk, jt, nvar, nlon, nlat, nlev_in, nlev_out, ntime, &
       m, y, month, syear, syear_in, time_int, year, months, ioerror, &
       status, ncid_in, ncid_out, lonid_in, lonid_out, latid_in, &
       latid_out, mlevid_in, levid_out, timeid_in, timeid_out, &
       nvartot, varnatts, vardimids(5), varndims, vartype, &
       indfirst, indlast, ind, iind
  REAL, ALLOCATABLE :: lon(:), lat(:), hyam(:), hybm(:), &
       aps_in(:,:,:), lev_in(:,:,:,:), levstand_in(:), lev_out(:), &
       time_in(:), time_out(:), var3(:,:,:), &
       var4_in(:,:,:,:), var4_out(:,:,:,:), &
       w0(:,:,:,:), w1(:,:,:,:), dwinv(:,:,:,:), &
       log_lev_in0(:,:,:,:), log_lev_in1(:,:,:,:), log_lev_out(:), &
       log_lev_out3d(:,:,:,:)
  LOGICAL, ALLOCATABLE :: ll(:,:,:,:)
  INTEGER, ALLOCATABLE :: varid_out(:), ind0(:), ind1(:), iind0(:,:,:,:), &
       iind1(:,:,:,:)
  CHARACTER(100) :: ncname_in, ncname_out, namelistfile, varname_in, varname_out
  CHARACTER(23) :: timeunits_in, timeunits_out
  CHARACTER(11) :: timeorigin_out
  CHARACTER(4) :: syear_st, syearm1_st
  LOGICAL :: laps=.FALSE.
  REAL, PARAMETER :: ps_stand = 101325  ! Standard surface pressure [Pa]
  INTEGER, PARAMETER :: lun=10

  INTEGER :: levels(99), levelscount(99)
  NAMELIST /afterburner_levels/ levels

  ! External subroutines:
  EXTERNAL interpolate_index

  ! Executable statements:

  CALL GETARG(1, ncname_in)
  CALL GETARG(2, ncname_out)
  CALL GETARG(3, namelistfile)
  CALL GETARG(4, syear_st)

  READ(syear_st,'(i4.4)') syear

  ! 1. Open / create file:

  ! 1.1 Open input file:
  status=NF90_OPEN(TRIM(ncname_in), NF90_NOWRITE, ncid_in)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! 1.2 Create output file:
  status=NF90_CREATE(TRIM(ncname_out), NF90_CLOBBER, ncid_out)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! 2. Get / set dimensions:

  ! 2.1 Get dimension lengths of input file:
  status=NF90_INQ_DIMID(ncid_in, 'lon', lonid_in)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQUIRE_DIMENSION(ncid_in, lonid_in, len = nlon)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQ_DIMID(ncid_in, 'lat', latid_in)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQUIRE_DIMENSION(ncid_in, latid_in, len = nlat)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQ_DIMID(ncid_in, 'mlev', mlevid_in)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQUIRE_DIMENSION(ncid_in, mlevid_in, len = nlev_in)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQ_DIMID(ncid_in, 'time', timeid_in)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_INQUIRE_DIMENSION(ncid_in, timeid_in, len = ntime)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Allocate arrays:
  ALLOCATE(lon(nlon))
  ALLOCATE(lat(nlat))
  ALLOCATE(hyam(nlev_in))
  ALLOCATE(hybm(nlev_in))
  ALLOCATE(levstand_in(nlev_in))
  ALLOCATE(lev_in(nlon,nlat,nlev_in,ntime))
  ALLOCATE(time_in(ntime))
  ALLOCATE(time_out(ntime))
  ALLOCATE(aps_in(nlon,nlat,ntime))
  ALLOCATE(var3(nlon,nlat,ntime))
  ALLOCATE(var4_in(nlon,nlat,nlev_in,ntime))

  ! 2.2 Define dimensions of output file:

  ! 2.2.1 Read pressure levels from namelist:
  levels(:) = 0
  levelscount(:) = 0
  ioerror = 0
  OPEN (lun, ERR=110, FILE=TRIM(namelistfile), STATUS='old', &
       FORM='formatted', ACTION='read')
  GOTO 111
110 ioerror = 1
111 CONTINUE
  IF (ioerror == 1) THEN
     PRINT*, 'Could not open', TRIM(namelistfile)
     STOP 'Stopped'
  ENDIF  
  READ(lun, NML=afterburner_levels)
  CLOSE(lun)

  ! 2.2.2 Determine number of levels:
  WHERE (levels .GT. 0) levelscount=1
  nlev_out=SUM(levelscount)

  ! 2.2.3 Allocate arrays:
  ALLOCATE(lev_out(nlev_out))
  ALLOCATE(var4_out(nlon,nlat,nlev_out,ntime))
  ALLOCATE(ind0(nlev_out))
  ALLOCATE(ind1(nlev_out))
  ALLOCATE(iind0(nlon,nlat,nlev_out,ntime))
  ALLOCATE(iind1(nlon,nlat,nlev_out,ntime))
  ALLOCATE(log_lev_in0(nlon,nlat,nlev_out,ntime))
  ALLOCATE(log_lev_in1(nlon,nlat,nlev_out,ntime))
  ALLOCATE(log_lev_out(nlev_out))
  ALLOCATE(log_lev_out3d(nlon,nlat,nlev_out,ntime))
  ALLOCATE(w0(nlon,nlat,nlev_out,ntime))
  ALLOCATE(w1(nlon,nlat,nlev_out,ntime))
  ALLOCATE(dwinv(nlon,nlat,nlev_out,ntime))
  ALLOCATE(ll(nlon,nlat,nlev_out,ntime))

  ! 2.2.4 Set out pressure grid:
  lev_out(:)=0.01*REAL(levels(1:nlev_out))  !hPa

  ! 2.2.5 Create dimensions of output file:
  status = NF90_DEF_DIM(ncid_out, 'lon', nlon, lonid_out)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_DIM(ncid_out, 'lat', nlat, latid_out)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_DIM(ncid_out, 'lev', nlev_out, levid_out)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_DIM(ncid_out, 'time', NF90_UNLIMITED, timeid_out)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! 3. Define variables in output file and get / set attributes:

  ! 3.1.1 Inquire number of variables (input file):
  status=NF90_INQUIRE(ncid_in, nvariables=nvartot)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ALLOCATE(varid_out(nvartot))


  write(*,*)nvartot
  DO i=1, nvartot

     ! 3.2.1 Get variable name, type, number of dimesions, dimension IDs, and
     ! number of variable attributes (input file):
     status=NF90_INQUIRE_VARIABLE(ncid_in, i, varname_in, vartype, varndims, &
          vardimids, varnatts)
     IF (status /= NF90_NOERR) CALL handle_err(status)

     ! Remove _m in variable name of output file (if present):
     ind = INDEX(TRIM(varname_in), '_m')
     IF (ind .GE. 1) THEN
        varname_out=varname_in(1:ind-1)
     ELSE
        varname_out=TRIM(varname_in)
     ENDIF

     ! 3.3.1 Define variables and set attributes (output file):

     IF (varndims .EQ. 2) THEN
        PRINT *, '2 dimensional arrays not allowed'
        STOP 'Stopped'
     ENDIF
     IF (varndims .EQ. 1) THEN
        SELECT CASE (TRIM(varname_in))

        CASE ('lon')
           ! Define variable:           
           status=NF90_DEF_VAR(ncid_out, TRIM(varname_out), NF90_FLOAT, &
                (/ lonid_out /), varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Set attributes:
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "long_name", "longitude")
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "units", "degrees east")
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "standard_name", "longitude")
           IF (status /= NF90_NOERR) CALL handle_err(status)

        CASE ('lat')
           ! Define variable:
           status=NF90_DEF_VAR(ncid_out, TRIM(varname_out), NF90_FLOAT, &
                (/ latid_out /), varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Set attributes:
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "long_name", "latitude")
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "units", "degrees north")
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "standard_name", "latitude")
           IF (status /= NF90_NOERR) CALL handle_err(status)

        CASE ('mlev')
           ! Define variable:
           status=NF90_DEF_VAR(ncid_out, 'lev', NF90_FLOAT, &
                (/ levid_out /), varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Set attributes:
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "long_name", "pressure")
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "units", "mb")
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "positive", "down")
           IF (status /= NF90_NOERR) CALL handle_err(status)

        CASE ('time')
           ! Get time units:
           status=NF90_GET_ATT(ncid_in, i, 'units', timeunits_in)
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Transform time:
           CALL transform_timeunits
           ! Define variable:
           status=NF90_DEF_VAR(ncid_out, TRIM(varname_out), NF90_FLOAT, &
                (/ timeid_out /), varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Set attributes:
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "long_name", "model time")
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "units", timeunits_out)
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_PUT_ATT(ncid_out, varid_out(i), "time_origin", &
                timeorigin_out)
           IF (status /= NF90_NOERR) CALL handle_err(status)
        END SELECT

     ELSE

        IF (varndims .EQ. 3) THEN
           ! Define variable:
           status=NF90_DEF_VAR(ncid_out, TRIM(varname_out), NF90_FLOAT, &
                (/ lonid_out, latid_out, timeid_out /), varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
        ELSE  ! varndims = 4
           ! Define variable:
           status=NF90_DEF_VAR(ncid_out, TRIM(varname_out), NF90_FLOAT, &
                (/ lonid_out, latid_out, levid_out, timeid_out /), varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
        ENDIF



        ! Copy attributes:
        status=NF90_COPY_ATT(ncid_in, i, "long_name", ncid_out, varid_out(i))
        IF (status /= NF90_NOERR) CALL handle_err(status)
!        status=NF90_COPY_ATT(ncid_in, i, "units", ncid_out, varid_out(i))
!        IF (status /= NF90_NOERR) CALL handle_err(status)
        IF (TRIM(varname_in) .NE. 'aps') THEN
           status=NF90_COPY_ATT(ncid_in, i, "code", ncid_out, varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_COPY_ATT(ncid_in, i, "table", ncid_out, varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
!           status=NF90_COPY_ATT(ncid_in, i, "axis", ncid_out, varid_out(i))
!           IF (status /= NF90_NOERR) CALL handle_err(status)
           status=NF90_COPY_ATT(ncid_in, i, "grid_type", ncid_out, varid_out(i))
           IF (status /= NF90_NOERR) CALL handle_err(status)
        ENDIF

write(*,*)varname_in
!pause 'tim'
     ENDIF

  ENDDO
!pause 'tim1'
  ! 3.4 Copy global attributes:
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "Conventions", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "title", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "source", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "echam_version", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "institution", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "advection", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "physics", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "radiation", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "date_time", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "operating_system", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "user_name", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "host_name", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status=NF90_COPY_ATT(ncid_in, NF90_GLOBAL, "truncation", ncid_out, NF90_GLOBAL)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! 3.5.2 End definition modus of output file:
  status=NF90_ENDDEF(ncid_out)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! 4. Get / set variable values:

  DO i=1, nvartot

     ! 4.1.1 Get variable name, type, number of dimesions, dimension IDs, and
     ! number of variable attributes (input file):
     status=NF90_INQUIRE_VARIABLE(ncid_in, i, varname_in, vartype, varndims, &
          vardimids, varnatts)
     IF (status /= NF90_NOERR) CALL handle_err(status)

     ! Remove _m in variable name of output file (if present):
     ind = INDEX(TRIM(varname_in), '_m')
     IF (ind .GE. 1) THEN
        varname_out=varname_in(1:ind-1)
     ELSE
        varname_out=TRIM(varname_in)
     ENDIF

     ! 4.2.1 Get variables values from input file and put values to output file
     !       (transformed to output grid, if necessary):

     IF (varndims .EQ. 1) THEN
        SELECT CASE (TRIM(varname_in))

        CASE ('lon')
           ! Get values from input file:
           status=NF90_GET_VAR(ncid_in, i, lon)
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Put values to output file:
           status=NF90_PUT_VAR(ncid_out, varid_out(i), lon)
           IF (status /= NF90_NOERR) CALL handle_err(status)

        CASE ('lat')
           ! Get values from input file:
           status=NF90_GET_VAR(ncid_in, i, lat)
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Put values to output file:
           status=NF90_PUT_VAR(ncid_out, varid_out(i), lat)
           IF (status /= NF90_NOERR) CALL handle_err(status)

        CASE ('mlev')
           ! Put values to output file:
           status=NF90_PUT_VAR(ncid_out, varid_out(i), lev_out)
           IF (status /= NF90_NOERR) CALL handle_err(status)

        CASE ('time')
           ! Get values from input file:
           status=NF90_GET_VAR(ncid_in, i, time_in)
           IF (status /= NF90_NOERR) CALL handle_err(status)
           ! Transform time:
           DO jt=1, ntime
              y=INT(time_in(jt)/365.25)  ! ~number of full years
              time_int=NINT((time_in(jt)-y*365.25)/30.)+12*y ! [days->months]
              time_in(jt) = REAL(time_int)
           ENDDO
           time_out(:)=time_in(:)+12.*(syear_in-syear)
           ! Put values to output file:
           status=NF90_PUT_VAR(ncid_out, varid_out(i), time_out)
           IF (status /= NF90_NOERR) CALL handle_err(status)

        CASE ('hyam')
           ! Get values from input file:
           status=NF90_GET_VAR(ncid_in, i, hyam)

        CASE ('hybm')
           ! Get values from input file:
           status=NF90_GET_VAR(ncid_in, i, hybm)

        END SELECT

     ELSEIF (varndims .EQ. 3) THEN
        ! Get values from input file:
        status=NF90_GET_VAR(ncid_in, i, var3)
        IF (status /= NF90_NOERR) CALL handle_err(status)
        ! Put values to output file:
        status=NF90_PUT_VAR(ncid_out, varid_out(i), var3)
        IF (status /= NF90_NOERR) CALL handle_err(status)
        ! Call *prepare_tranform_lev* if TRIM(varname_in)=='aps':
        IF (TRIM(varname_in) .EQ. 'aps') THEN
           laps = .TRUE.
           aps_in = var3
           CALL prepare_transform_lev
        ENDIF

     ELSE  ! varndims = 4
        IF (.NOT. laps) THEN
           PRINT *, 'Surface pressure field "aps" has to be read before'
           STOP 'Stopped'
        ENDIF
        ! Get values from input file:
        status=NF90_GET_VAR(ncid_in, i, var4_in)
        IF (status /= NF90_NOERR) CALL handle_err(status)
        ! Interpolate var4_in to pressure grid of output grid:
        CALL transform_lev
        ! Put values to output file:
        status=NF90_PUT_VAR(ncid_out, varid_out(i), var4_out)
        IF (status /= NF90_NOERR) CALL handle_err(status)
     ENDIF
  ENDDO

  ! 5. Close files:
  
  status = NF90_CLOSE(ncid_in)
  IF (status /= NF90_NOERR) CALL HANDLE_ERR(status)
  status = NF90_CLOSE(ncid_out)
  IF (status /= NF90_NOERR) CALL HANDLE_ERR(status)

  ! 6. Deallocate memory:

  DEALLOCATE(lon)
  DEALLOCATE(lat)
  DEALLOCATE(hyam)
  DEALLOCATE(hybm)
  DEALLOCATE(aps_in)
  DEALLOCATE(lev_in)
  DEALLOCATE(levstand_in)
  DEALLOCATE(lev_out)
  DEALLOCATE(time_in)
  DEALLOCATE(time_out)
  DEALLOCATE(iind0)
  DEALLOCATE(iind1)
  DEALLOCATE(var3)
  DEALLOCATE(var4_in)
  DEALLOCATE(var4_out)
  DEALLOCATE(w0)
  DEALLOCATE(w1)
  DEALLOCATE(dwinv)
  DEALLOCATE(log_lev_in0)
  DEALLOCATE(log_lev_in1)
  DEALLOCATE(log_lev_out)
  DEALLOCATE(log_lev_out3d)
  DEALLOCATE(ll)
  DEALLOCATE(varid_out)
  DEALLOCATE(ind0)
  DEALLOCATE(ind1)

CONTAINS

  SUBROUTINE transform_timeunits

    ! Transforms time record to output units "months since <syear>".

    IF (timeunits_in(1:10) .EQ. 'days since') THEN 
       READ(timeunits_in(12:15),'(i4.4)') syear_in
    ELSEIF (timeunits_in(1:9) .EQ. 'day since') THEN
       READ(timeunits_in(11:14),'(i4.4)') syear_in    
    ELSE
       PRINT *, 'Invalid time units'
       STOP 'Stopped'
    ENDIF

    ! Set out time units:
    WRITE(syearm1_st,'(i4)') syear-1
    timeunits_out(1:19)= 'month since 15-DEC-'
    timeunits_out(20:23)=syearm1_st
    timeorigin_out(1:7)='15-DEC-'
    timeorigin_out(8:11)=syearm1_st

  END SUBROUTINE transform_timeunits

  SUBROUTINE prepare_transform_lev

    ! Preparates level transformation: Calculates interpolation indices
    ! and interpolation weights.

    ! Pressure grid [hPa] for standard surface pressure:
    levstand_in(:)=0.01*(hyam(:)+hybm(:)*ps_stand)

    ! Calculate pressure interpolation indices:
    CALL interpolate_index(nlev_in, nlev_out, levstand_in, lev_out, &
         indfirst, indlast, ind1)
    ! iind1: index of next following level of input gird with respect
    !        to a given levels of output grid
    ind0(:)=ind1(:)-1
    WHERE (ind0 .GT. nlev_in) ind0=nlev_in
    WHERE (ind1 .GT. nlev_in) ind1=nlev_in
    WHERE (ind0 .LT. 1) ind0=1
    WHERE (ind1 .LT. 1) ind1=1

    ! Spread pressure levels onto horizontal grid:
    FORALL (jl=1:nlon, jk=1:nlat, jt=1:ntime)
       iind0(jl,jk,:,jt)=ind0(:)
       iind1(jl,jk,:,jt)=ind1(:)
    END FORALL

    ! Pressure grid [hPa] for surface pressure 'aps':
    FORALL (jk= 1:nlev_in) &
         lev_in(:,:,jk,:)=0.01*(hyam(jk)+hybm(jk)*aps_in(:,:,:))

    ! Correct iind0, iind1 upwards / downwards if necessary:
    DO ji=1, nlon
       DO jk=1, nlat
          DO jl=1, nlev_out
             DO jt=1, ntime
                DO
                   iind=iind1(ji,jk,jl,jt)
                   IF (lev_in(ji,jk,iind,jt) .LT. lev_out(jl) &
                        .AND. iind .LT. nlev_in) THEN
                      iind0(ji,jk,jl,jt)=iind0(ji,jk,jl,jt)+1
                      iind1(ji,jk,jl,jt)=iind1(ji,jk,jl,jt)+1
                   ELSE
                      IF (lev_in(ji,jk,iind,jt) .LT. lev_out(jl)) &
                           iind0(ji,jk,jl,jt)=nlev_in
                      EXIT
                   ENDIF
                ENDDO

                DO
                   iind=iind0(ji,jk,jl,jt)
                   IF (lev_out(jl) .LT. lev_in(ji,jk,iind,jt) &
                        .AND. iind .GT. 1) THEN
                      iind0(ji,jk,jl,jt)=iind0(ji,jk,jl,jt)-1
                      iind1(ji,jk,jl,jt)=iind1(ji,jk,jl,jt)-1
                   ELSE
                      IF (lev_out(jl) .LT. lev_in(ji,jk,iind,jt)) &
                           iind1(ji,jk,jl,jt)=1
                      EXIT
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
       ENDDO
    ENDDO

    WHERE (iind0 .GT. nlev_in) iind0=nlev_in
    WHERE (iind1 .GT. nlev_in) iind1=nlev_in
    WHERE (iind0 .LT. 1) iind0=1
    WHERE (iind1 .LT. 1) iind1=1

    ! Weights:

    ! Negative logarithm of pressure grid for interpolation
    ! <=> interpolation on height grid
    DO ji=1, nlon
       DO jk=1, nlat
          DO jt=1, ntime
            ind0(:)=iind0(ji,jk,:,jt)
            ind1(:)=iind1(ji,jk,:,jt)
            log_lev_in0(ji,jk,:,jt)=-LOG(lev_in(ji,jk,ind0(:),jt))
            log_lev_in1(ji,jk,:,jt)=-LOG(lev_in(ji,jk,ind1(:),jt))
          ENDDO
       ENDDO
    ENDDO
    log_lev_out(:)=-LOG(lev_out(:))

    ! Spread pressure interpolation indices onto horizontal grid:
    FORALL (jl=1:nlon, jk=1:nlat, jt=1:ntime) &
       log_lev_out3d(jl,jk,:,jt)=log_lev_out(:)

    ll(:,:,:,:)=iind0(:,:,:,:) .LT. iind1(:,:,:,:)
    w0=MERGE(log_lev_in1-log_lev_out3d, 1., ll)
    w1=MERGE(log_lev_out3d-log_lev_in0, 0., ll)
    WHERE (ll)   ! MERGE does not work here (division by 0 for ll=.FALSE.)
       dwinv=1./(log_lev_in1-log_lev_in0)
    ELSEWHERE
       dwinv=1.
    END WHERE

  END SUBROUTINE prepare_transform_lev

  SUBROUTINE transform_lev

    ! Pressure levels: interpolation from input to output grid.

    DO ji=1, nlon
       DO jk=1, nlat
          DO jt=1, ntime
             ind0(:)=iind0(ji,jk,:,jt)
             ind1(:)=iind1(ji,jk,:,jt)
             var4_out(ji,jk,:,jt) = &
                  (w0(ji,jk,:,jt)*var4_in(ji,jk,ind0(:),jt)+ &
                  w1(ji,jk,:,jt)*var4_in(ji,jk,ind1(:),jt)) &
                  * dwinv(ji,jk,:,jt)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE transform_lev

  SUBROUTINE HANDLE_ERR(status)
    
    INTEGER, INTENT(IN) :: status
    
    IF (status /= NF90_NOERR) THEN
       PRINT *, TRIM(NF90_STRERROR(status))
       STOP
    END IF
    
  END SUBROUTINE HANDLE_ERR
  
END PROGRAM echamnetcdf_afterburner

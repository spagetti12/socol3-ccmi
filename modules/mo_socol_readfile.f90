MODULE mo_socol_readfile

  ! Description:

  ! This module contains subroutines for reading data from files.

  USE mo_control,                 ONLY: nlon, ngl, nlev
  USE mo_decomposition,           ONLY: lc => local_decomposition, &
                                        global_decomposition
  USE mo_io
  USE mo_mpi,                     ONLY: p_io, p_parallel, p_parallel_io, &
                                        p_bcast
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_kind,                    ONLY: dp
  USE mo_netcdf,                  ONLY: FILE_INFO, io_get_vara_double
  USE mo_socol_grid_calculations, ONLY: check_dim_socol, &
                                        calculate_pres_stand_socol
  USE mo_transpose,               ONLY: scatter_gp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: socol_ascii_yearmo_tab_y, socol_read_netcdf, socol_ascii_ap, socol_ascii_ionpair, &
              socol_read_ionization_netcdf_1d, socol_read_ionization_netcdf_3d, &
              socol_ascii_geomag

CONTAINS
 
  FUNCTION socol_ascii_yearmo_tab_y(nval, yr, filename)

    ! Returns monthly data values for year <yr> plus December of preceeding /
    ! January of following year from the ascii file <filename>, where the
    ! data in the file are grouped as follows:
    ! 
    ! YYYYMM value(1) value(2) ... value(nval)
    ! YYYYMM value(1) value(2) ... value(nval)
    !   :       :      :          :
    
    ! M. Schraner, ETH Z?rich, January 2009
    
    IMPLICIT NONE
    
    ! Function arguments:
    INTEGER, INTENT(in) :: nval, yr
    CHARACTER(*), INTENT(in) :: filename
    
    ! Function output:
    REAL(dp), DIMENSION(nval,0:13) :: socol_ascii_yearmo_tab_y
    
    ! Local variables:
    INTEGER :: ioerror, yearmo, yearmo0, m
    INTEGER, PARAMETER :: lun = 20
    CHARACTER(4) :: yr_st
    CHARACTER(2) :: mo_st


    ! Executable statements: 
    
    ioerror = 0
    
    IF (p_parallel_io) THEN
       
       ! Open file:
       OPEN (lun, ERR=110, FILE=filename, STATUS='old', &
            FORM='formatted', ACTION='read')
       GOTO 111
110    ioerror = 1
111    CONTINUE
       
       ! Read file:
       IF (ioerror .EQ. 0) THEN
          yearmo0 = 100*(yr-1)+12
          DO
             READ (lun,*, END=120, ERR=120) &
                  yearmo, socol_ascii_yearmo_tab_y(1:nval,0)
             IF (yearmo .EQ. yearmo0) THEN
                EXIT
             ENDIF
          ENDDO
120       CONTINUE
          IF (yearmo .NE. yearmo0) ioerror=2
       ENDIF
       IF (ioerror .EQ. 0) THEN
          DO m = 1, 13
             IF (m .EQ. 13) THEN
                yearmo0 = 100*(yr+1)+1
             ELSE
                yearmo0 = 100*yr+m
             ENDIF
             READ (lun,*, END=121, ERR=121) &
                  yearmo, socol_ascii_yearmo_tab_y(1:nval,m)
121          CONTINUE
             IF (yearmo .NE. yearmo0) THEN
                ioerror = 2
                EXIT
             ENDIF
          ENDDO
       ENDIF
       
       ! Close file:
       CLOSE (lun)
       
    ENDIF
    
    ! Terminate if ioerror > 0:
    IF (p_parallel) CALL p_bcast (ioerror, p_io)
    SELECT CASE (ioerror)
    CASE(1)
       WRITE(message_text,*) 'Could not open ', TRIM(filename)
       CALL message('',TRIM(message_text))
       CALL finish('socol_ascii_yearmo_tab_y','Run terminated')
    CASE(2)
       WRITE(yr_st, '(i4)') yearmo0/100
       IF (yearmo0-yearmo0/100 .LE. 9) THEN
          WRITE(mo_st, '(i1)') yearmo0-yearmo0/100
       ELSE
          WRITE(mo_st, '(i2)') yearmo0-yearmo0/100
       ENDIF
       WRITE(message_text,*) &
            'Year ', yr_st, ', month ', TRIM(mo_st), ' not found in file <', &
            TRIM(filename),'>'
       CALL message('',TRIM(message_text))
       CALL finish ('socol_ascii_yearmo_tab_y', 'Run terminated.')
    END SELECT
    
    ! Distribute data:
    IF (p_parallel) CALL p_bcast (socol_ascii_yearmo_tab_y, p_io)    
    
  END FUNCTION socol_ascii_yearmo_tab_y


  SUBROUTINE socol_read_netcdf(filename0, varname, dimnames, extradimname, &
       extradim2name, llevnoecham, yr, mo, mo14, varname_longname, data6d, &
       data5d, data4d, data3d, data2d, data1d, lowlevind, uplevind, &
       nextradim, extradim, nextradim2, extradim2, nlevdimnoecham, &
       levdimnoecham)

    ! Reads data array from netcdf file depending on the following settings.
    !
    ! INPUT:
    !
    ! <filename0>   : (STRING): Name of the data file. If <yr> is present, the 
    !                 filename must have the format filename0YYYY with 
    !                 YYYY=<yr>, i.e. <filename0> is the filename without the 
    !                 year YYYY. Otherwise filename=<filename0>
    ! <varname>     : (STRING): The name of the variable (array) to be read
    ! <dimnames>    : (STRING): Allowed values are:
    !                 - 'lonlatp', 'LONLATP', 'lonlatlev', or 'LONLATLEV', 
    !                   if data=data(lon,lat,lev) (according to the dimension
    !                   names in the netcdf file, i.e. if LONLATLEV, the 
    !                   corresponding dimension names of the netcdf files are 
    !                   'LON', 'LAT', and 'LEV')
    !                 - 'lonlat' or 'LONLAT', if data=data(lon,lat)
    !                 - 'latp', 'LATP', 'latlev', or 'LATLEV', 
    !                   if data=data(lat,lev)
    !                 - 'p', 'P', 'lev', or 'LEV', if data=data(lev)
    !                 - '-', if there is not a geographical data set, i.e. 
    !                   only additional dimension(s) are present and neihter 
    !                   'lon' nor 'lat' nor 'lev' are present
    ! <extradimname>: (STRING, OPTIONAL): Dimension name of a (possible) 
    !                 additional dimension (e.g. wave length)
    ! <extradim2name>: (STRING, OPTIONAL): Dimension name of a (possible) 
    !                 second additional dimension (e.g. wave length)
    ! <llevnoecham> : (LOGICAL, OPTIONAL): If true, the pressure grid of the
    !                 data set is not identical to the one of ECHAM5. The 
    !                 pressure grid of the data set is given as output by the
    !                 the vector <levdimnoecham> and it's extent
    !                 <nlevdimnoecham>. Default: false.
    ! <yr>          : (INTEGER, OPTIONAL): The current year (for yearly varying
    !                 data sets) 
    ! <mo>          : (INTEGER, OPTIONAL): The current month (for monthly 
    !                 varying data sets). In this case, the output array has
    !                 the additional dimension <m>, with m=0, 1, 2, where
    !                 - m=1: current month
    !                 - m=0: preceeding month
    !                 - m=2: following month
    ! <mo14>        : (LOGICAL, OPTIONAL): Applicable instead of <mo>: If
    !                 <mo14>=.TRUE., the dimension length of m is 14 containing
    !                 all months of a year plus December of preceeding year/
    !                 January of following year are written in the output array.
    !
    ! <varaname_longname>: (STING, OPTIONAL): Long name of the variable (used
    !                 for messages). If omitted, <varname> is used.
    !                 If set to '-', no message is printed.
    ! 
    ! OUTPUT:
    !
    ! data6d        : (REAL(:,:,:,:,:,:), OPTIONAL) if the output array has 
    !                                              rank 6
    ! OR:
    ! data5d        : (REAL(:,:,:,:,:), OPTIONAL) if the output array has rank 5
    ! OR:
    ! data4d        : (REAL(:,:,:,:), OPTIONAL)   if the output array has rank 4
    ! OR:
    ! data3d        : (REAL(:,:,:), OPTIONAL)     if the output array has rank 3
    ! OR:
    ! data2d        : (REAL(:,:), OPTIONAL)       if the output array has rank 2
    ! OR:
    ! data1d        : (REAL(:), OPTIONAL)         if the output array has rank 1
    ! <lowlevind>, <uplevind> : (INTEGER, OPTIONAL): If the data array
    !                 depends on the pressure dimension and 
    !                 <llevnoecham>=.false., <lowlevind> / <uplevind> indicate 
    !                 the lower / upper bounds of the vertical dimension of the
    !                 output array such that the vertical indices of the model
    !                 output array agree with the ones of ECHAM5.
    ! <nextradim>   : (INTEGER, OPTIONAL): Extent of dimension <extradimname>
    ! <extradim>    : (REAL(:), OPTIONAL): Grid of dimension <extradimname>
    ! <nextradim2>  : (INTEGER, OPTIONAL): Extent of dimension <extradim2name>
    ! <extradim2>   : (REAL(:), OPTIONAL): Grid of dimension <extradim2name>
    ! <nlevdimnoecham> : (INTEGER, OPTIONAL): Extent of dimension 
    !                 <levdimnoecham>.
    ! <levdimnoecham> : (REAL(:), OPTIONAL): Grid of pressure dimension of the 
    !                 data set (in the case that <llevnoecham>=.true., i.e.
    !                 the pressure grid of the data set does not agree with 
    !                 the one of ECHAM5).

    ! *socol_read_netcdf* is called from several boundary condition 
    ! subroutines.
    
    ! M. Schraner, ETH Z?rich, February 2009
    
    IMPLICIT NONE
    
    ! Subroutine arguments:
    CHARACTER(*), INTENT(in)           :: filename0, varname, dimnames
    CHARACTER(*), OPTIONAL, INTENT(in) :: extradimname, extradim2name
    CHARACTER(*), OPTIONAL, INTENT(in) :: varname_longname
    INTEGER, OPTIONAL, INTENT(in)      :: yr, mo
    LOGICAL, OPTIONAL, INTENT(in)      :: mo14
    LOGICAL, OPTIONAL, INTENT(in)      :: llevnoecham
    REAL(dp), ALLOCATABLE, OPTIONAL, INTENT(out) :: data6d(:,:,:,:,:,:), &
         data5d(:,:,:,:,:), data4d(:,:,:,:), data3d(:,:,:), data2d(:,:), &
         data1d(:), extradim(:), extradim2(:), levdimnoecham(:)
    INTEGER, OPTIONAL, INTENT(out)     :: lowlevind, uplevind, nextradim, &
                                          nextradim2, nlevdimnoecham
    
    ! Local variables:
    ! m0, m1, f0, f1, mof0, mof1, moa0, mo1: 
    ! 0: lower boundary, 1: upper boundary
    INTEGER       :: m0, m1               ! array dimension <m>
    INTEGER       :: f0, f1               ! loop over files
    INTEGER       :: mof0(0:2), mof1(0:2) ! month index in netdf file
    INTEGER       :: moa0(0:2), moa1(0:2) ! corresponding month index in array
    INTEGER       :: f, i, i0, ii, m, io_var_id, rank, rankf, nlonloc, &
                     nlatloc, lowlevindloc, uplevindloc, nextradimloc, &
                     nextradim2loc, nlevdimnoechamloc
    INTEGER       :: start(6), count(6)
    REAL(dp), ALLOCATABLE :: zin(:,:,:,:,:,:)
    REAL(dp), ALLOCATABLE, TARGET :: zin_tr(:,:,:,:,:,:)
    REAL(dp), POINTER :: data_gl(:,:,:,:)
    REAL(dp), ALLOCATABLE :: data_loc(:,:,:,:,:,:)
    REAL(dp), ALLOCATABLE :: data_loc4d(:,:,:,:)
    REAL(dp), ALLOCATABLE :: extradimloc(:), extradimbefore(:), &
         extradim2loc(:), extradim2before(:), levdimnoechamloc(:), &
         levdimnoechambefore(:)
    TYPE (FILE_INFO) :: socolnc
    LOGICAL       :: lex, mo14loc
    LOGICAL       :: llondim, llatdim, llevdim, lextradim, lextradim2, &
         llevnoechamloc
    CHARACTER(3)  :: londimname, latdimname, levdimname
    CHARACTER(32) :: fn(0:2)
    CHARACTER(40) :: vn
    CHARACTER(4)  :: yr0_st, yr1_st, yr2_st
    CHARACTER(2)  :: mo0_st, mo1_st, mo2_st

    ! Executable statements:

    IF (PRESENT(mo14)) THEN
       mo14loc=.TRUE. 
    ELSE
       mo14loc=.FALSE.
    ENDIF

    IF (PRESENT(llevnoecham)) THEN
       llevnoechamloc=llevnoecham
    ELSE
       llevnoechamloc=.FALSE.
    ENDIF

    IF (PRESENT(extradimname)) THEN
       lextradim=.TRUE.
    ELSE
       lextradim=.FALSE.
       nextradimloc=1
    ENDIF

    IF (PRESENT(extradim2name)) THEN
       lextradim2=.TRUE.
       IF (.NOT. lextradim) THEN
          WRITE (message_text,*) 'if extradim2name is present, also extradimname has to be present'
          CALL message('',TRIM(message_text))
       ENDIF
    ELSE
       lextradim2=.FALSE.
       nextradim2loc=1
    ENDIF

    ! Initial value for open/close status of socolnc:
    socolnc%opened = .FALSE.

    ! Set array rank according to dimnames:
    IF (TRIM(dimnames) .EQ. 'lonlatp' .OR. TRIM(dimnames) .EQ. 'LONLATP' .OR. &
         TRIM(dimnames) .EQ. 'lonlatlev' .OR. TRIM(dimnames) .EQ. 'LONLATLEV') &
         rank=4

    IF (TRIM(dimnames) .EQ. 'lonlat' .OR. TRIM(dimnames) .EQ. 'LONLAT') rank=3
    
    IF (TRIM(dimnames) .EQ. 'latp' .OR. TRIM(dimnames) .EQ. 'LATP' .OR. &
         TRIM(dimnames) .EQ. 'latlev' .OR. TRIM(dimnames) .EQ. 'LATLEV') rank=3

    IF (TRIM(dimnames) .EQ. 'p' .OR. TRIM(dimnames) .EQ. 'P' .OR. &
         TRIM(dimnames) .EQ. 'lev' .OR. TRIM(dimnames) .EQ. 'LEV') rank=2

    IF (TRIM(dimnames) .EQ. '-') rank=1

    ! Enhance rank by 1 if extradim is present:
    IF (lextradim) rank=rank+1

    ! Enhance rank by 1 if extradim2 is present:
    IF (lextradim2) rank=rank+1

    ! Set dimension names and dimension switches according to dimnames:
    IF (TRIM(dimnames) .EQ. 'lonlatp' .OR. TRIM(dimnames) .EQ. 'lonlatlev' &
         .OR. TRIM(dimnames) .EQ. 'lonlat') THEN
       londimname='lon'
       llondim=.TRUE.
    ELSEIF (TRIM(dimnames) .EQ. 'LONLATP' .OR. TRIM(dimnames) .EQ. &
         'LONLATLEV' .OR. TRIM(dimnames) .EQ. 'LONLAT') THEN
       londimname='LON'
       llondim=.TRUE.
    ELSE
       llondim=.FALSE.
    ENDIF

    ! Set local number of longitudes:
    IF (llondim) THEN
       nlonloc = nlon
    ELSE
       nlonloc = 1
    ENDIF

    IF (TRIM(dimnames) .EQ. 'lonlatp' .OR. TRIM(dimnames) .EQ. 'lonlatlev' &
         .OR. TRIM(dimnames) .EQ. 'latp' .OR. TRIM(dimnames) .EQ. 'latlev' &
         .OR. TRIM(dimnames) .EQ. 'lonlat') THEN
       latdimname='lat'
       llatdim=.TRUE.
    ELSEIF (TRIM(dimnames) .EQ. 'LONLATP' .OR. TRIM(dimnames) .EQ. 'LONLATLEV' &
         .OR. TRIM(dimnames) .EQ. 'LATP' .OR. TRIM(dimnames) .EQ. 'LATLEV' &
         .OR. TRIM(dimnames) .EQ. 'LONLAT') THEN
       latdimname='LAT'
       llatdim=.TRUE.
    ELSE
       llatdim=.FALSE.
    ENDIF

    ! Set local number of latitudes:
    IF (llatdim) THEN
       nlatloc = ngl
    ELSE
       nlatloc = 1
    ENDIF

    IF (TRIM(dimnames) .EQ. 'lonlatp' .OR. TRIM(dimnames) .EQ. 'latp' .OR. &
         TRIM(dimnames) .EQ. 'p' .OR. TRIM(dimnames) .EQ. 'LONLATP' .OR. &
         TRIM(dimnames) .EQ. 'LATP' .OR. TRIM(dimnames) .EQ. 'P' .OR. &
         TRIM(dimnames) .EQ. 'lonlatlev' .OR. TRIM(dimnames) .EQ. 'latlev' .OR. &
         TRIM(dimnames) .EQ. 'lev' .OR. TRIM(dimnames) .EQ. 'LONLATLEV' .OR. &
         TRIM(dimnames) .EQ. 'LATLEV' .OR. TRIM(dimnames) .EQ. 'LEV') THEN

       IF (TRIM(dimnames) .EQ. 'lonlatp' .OR. TRIM(dimnames) .EQ. 'latp' .OR. &
            TRIM(dimnames) .EQ. 'p') levdimname='p'
       IF (TRIM(dimnames) .EQ. 'LONLATP' .OR. TRIM(dimnames) .EQ. 'LATP' .OR. &
            TRIM(dimnames) .EQ. 'P') levdimname='P'
       IF (TRIM(dimnames) .EQ. 'lonlatlev' .OR. TRIM(dimnames) .EQ. 'latlev' &
            .OR. TRIM(dimnames) .EQ. 'lev') levdimname='lev'
       IF (TRIM(dimnames) .EQ. 'LONLATLEV' .OR. TRIM(dimnames) .EQ. 'LATLEV' &
            .OR. TRIM(dimnames) .EQ. 'LEV') levdimname='LEV'

       llevdim=.TRUE.
    ELSE
       llevdim=.FALSE.
    ENDIF

    ! No levels:
    IF (.NOT. llevdim) THEN
       lowlevindloc=1
       uplevindloc=1
    ENDIF
    
    ! Set dimension lengths for months:
    IF (mo14loc) THEN
       m0=0
       m1=13
       f0=0
       f1=2
       mof0(0:2)=(/ 12, 1, 1 /)
       mof1(0:2)=(/ 12, 12, 1 /)
       moa0(0:2)=(/ 0, 1, 13 /) 
       moa1(0:2)=(/ 0, 12, 13 /)
    ELSEIF (PRESENT(mo)) THEN
       m0=0
       m1=2
       IF (mo .EQ. 1) THEN
          f0=0
          mof0(0:1)=(/ 12, 1 /)
          mof1(0)=12
          moa0(0:1)=(/ 0, 1 /)
          moa1(0)=0
       ELSE 
          f0=1
          mof0(1)=mo-1
          moa0(1)=0
       ENDIF
       IF (mo .EQ. 12) THEN
          f1=2
          mof0(2)=1
          mof1(1:2)=(/ 12, 1 /)
          moa0(2)=2
          moa1(1:2)=(/ 1, 2 /)
       ELSE 
          f1=1
          mof1(1)=mo+1
          moa1(1)=2
       ENDIF         
    ELSE      
       m0=1
       m1=1
       f0=1
       f1=1
       mof0(1)=1
       mof1(1)=1
       moa0(1)=1
       moa1(1)=1
    ENDIF

    ! Determine start(1:rank-1) (used for reading netcdf) and rankf:
    IF (p_parallel_io) THEN
       start(1:rank-1) = 1
       IF (m0 .EQ. m1) THEN
          rankf=rank-1
       ELSE
          rankf=rank
       ENDIF
    ENDIF

    ! Determine filenames:
    IF (PRESENT(yr)) THEN
       WRITE(fn(0), '(a,i4)') TRIM(filename0), yr-1
       WRITE(fn(1), '(a,i4)') TRIM(filename0), yr
       WRITE(fn(2), '(a,i4)') TRIM(filename0), yr+1
    ELSE
       fn(1) = TRIM(filename0)
       fn(0) = fn(1)
       fn(2) = fn(1)
    ENDIF

    ! Loop over files:
    DO f = f0, f1
       IF (p_parallel_io) INQUIRE (FILE=TRIM(fn(f)), EXIST=lex)
          
       ! Terminate if file does not exist:
       IF (p_parallel) CALL p_bcast (lex, p_io)
       IF (.NOT. lex) THEN
          WRITE (message_text,*) 'Could not open file <', TRIM(fn(f)), '>'
          CALL message('',TRIM(message_text))
          CALL finish ('socol_read_netcdf', 'Run terminated.')
       ENDIF      

       ! Open file:
       IF (p_parallel_io) CALL IO_open (TRIM(fn(f)), socolnc, IO_READ)

       ! Check if data grid agrees with ECHAM grid, determine uplevind /
       ! lowlevind (for pressure grids if <llevdimnoecham>=.false.), 
       ! nlevdimnoechamloc / levdimnoechamloc (for pressure grids if
       ! <llevdimnoecham>=.true.), nextradimloc / extradimloc (if dimension
       ! <extradimname> is present), nextradim2loc / extradim2loc (if dimension
       ! <extradimname> is present) and set lowlevindloc / uplevindloc (for
       ! pressure grids):
       IF (llondim) CALL check_dim_netcdf(londimname)
       IF (llatdim) CALL check_dim_netcdf(latdimname)
       IF (llevdim) CALL check_dim_netcdf(levdimname)
       IF (lextradim) CALL check_dim_netcdf(extradimname)
       IF (lextradim2) CALL check_dim_netcdf(extradim2name)

       ! Allocate memory for global and local arrays:
       IF (.NOT. llondim .AND. .NOT. llatdim) THEN 
          IF (p_parallel_io .AND. .NOT. ALLOCATED(zin)) &
               ALLOCATE(zin(1,1,lowlevindloc:uplevindloc,nextradimloc, &
               nextradim2loc,m0:m1))
          IF (.NOT. ALLOCATED(data_loc)) &
               ALLOCATE(data_loc(1,lowlevindloc:uplevindloc,nextradimloc, &
               nextradim2loc,1,m0:m1))
          IF (.NOT. ALLOCATED(data_loc4d)) &
               ALLOCATE(data_loc4d(lowlevindloc:uplevindloc,nextradimloc, &
               nextradim2loc,m0:m1))
       ELSE
          IF (p_parallel_io .AND. .NOT. ALLOCATED(zin)) &
               ALLOCATE(zin(nlon,ngl,lowlevindloc:uplevindloc,nextradimloc, &
                     nextradim2loc,m0:m1))
          IF (.NOT. ALLOCATED(data_loc)) &
               ALLOCATE(data_loc(lc%nproma,lowlevindloc:uplevindloc, &
               nextradimloc,nextradim2loc,lc%ngpblks,m0:m1))
       ENDIF
       IF (p_parallel_io .AND. .NOT. ALLOCATED(zin_tr)) &
            ALLOCATE(zin_tr(SIZE(zin,1),LBOUND(zin,3):UBOUND(zin,3), &
            SIZE(zin,4),SIZE(zin,5),SIZE(zin,2),LBOUND(zin,6):UBOUND(zin,6)))

       IF (p_parallel_io) THEN
          ! Determine count(1:rank-1):
          ! Skip longitudinal dimension for 'latlev'/'LATLEV', '-':
          IF (TRIM(dimnames) .EQ. 'latlev' .OR. TRIM(dimnames) .EQ. 'LATLEV') &
               THEN
             i0 = 2
          ELSE
             i0 = 1
          ENDIF
          ii = 1
          DO i = i0, 5             
             IF (SIZE(zin,i) .GE. 2) THEN
                count(ii)=SIZE(zin,i)
                ii=ii+1
             ENDIF
          ENDDO

          ! Read data:
          start(rank)=mof0(f)
          count(rank)=mof1(f)-mof0(f)+1
          CALL IO_inq_varid (socolnc%file_id, TRIM(varname), io_var_id)
          CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
               count(1:rankf), zin(1:nlonloc,1:nlatloc,:,:,:,moa0(f):moa1(f)))

          ! Scatter array over longitudinal grid, if data set depends on
          ! latitudes, but does not depend on longitudes:
          IF (.NOT. llondim .AND. llatdim) THEN
             DO i = 1, nlon   !FORALL (i = 1:nlon)
                zin(i,:,:,:,:,:) = zin(1,:,:,:,:,:)
             END DO
          ENDIF

          ! Close file:
          CALL IO_close(socolnc)
       ENDIF
    ENDDO
    
    IF (llondim .OR. llatdim) THEN

       IF (p_parallel_io) THEN
          ! Permute (lat) and (lev, extradim, extradim2) (order required for 
          ! *scatter_gp*):
          DO i = LBOUND(zin,2), UBOUND(zin,2)   !FORALL (i = LBOUND(zin,2):UBOUND(zin,2))
             zin_tr(:,:,:,:,i,:)=zin(:,i,:,:,:,:)
          END DO
       ENDIF
       
       ! Scatter fields over processors :
       DO i = 1, nextradimloc
          DO ii = 1, nextradim2loc
             DO m = m0, m1
                NULLIFY(data_gl)
                IF (p_parallel_io) data_gl => zin_tr(:,:,i,ii,:,m:m)
                CALL scatter_gp &
                     (data_gl, data_loc(:,:,i,ii,:,m:m), global_decomposition)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       IF (p_parallel_io) data_loc4d(:,:,:,:) = zin(1,1,:,:,:,:)

       ! Distribute data:
       IF (p_parallel) CALL p_bcast (data_loc4d, p_io)
       data_loc(1,:,:,:,1,:) = data_loc4d(:,:,:,:)
    ENDIF

    ! Allocate data to output array:
    IF (.NOT. llondim .AND. llatdim) rank=rank+1

    IF (lextradim) THEN
       IF (lextradim2) THEN
          SELECT CASE(rank)
          CASE(6)   
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data5d)) DEALLOCATE(data5d)
                ALLOCATE(data5d(SIZE(data_loc,1), &
                     LBOUND(data_loc,2):UBOUND(data_loc,2),SIZE(data_loc,3), &
                     SIZE(data_loc,4),SIZE(data_loc,5)))
                data5d(:,:,:,:,:) = data_loc(:,:,:,:,:,1)
             ELSE
                IF (ALLOCATED(data6d)) DEALLOCATE(data6d)
                ALLOCATE(data6d(SIZE(data_loc,1), &
                     LBOUND(data_loc,2):UBOUND(data_loc,2),SIZE(data_loc,3), &
                     SIZE(data_loc,4),SIZE(data_loc,5), &
                     LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data6d(:,:,:,:,:,:) = data_loc(:,:,:,:,:,:)
             ENDIF
          CASE(5)  !lev is missing
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data4d)) DEALLOCATE(data4d)
                ALLOCATE(data4d(SIZE(data_loc,1),SIZE(data_loc,3), &
                     SIZE(data_loc,4),SIZE(data_loc,5)))
                data4d(:,:,:,:) = data_loc(:,1,:,:,:,1)
             ELSE
                IF (ALLOCATED(data5d)) DEALLOCATE(data5d)
                ALLOCATE(data5d(SIZE(data_loc,1),SIZE(data_loc,3), &
                     SIZE(data_loc,4),SIZE(data_loc,5), &
                     LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data5d(:,:,:,:,:) = data_loc(:,1,:,:,:,:)
             ENDIF
          CASE(4)  !lon, lat are missing
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data3d)) DEALLOCATE(data3d)
                ALLOCATE(data3d(LBOUND(data_loc,2):UBOUND(data_loc,2), &
                     SIZE(data_loc,3),SIZE(data_loc,4)))
                data3d(:,:,:) = data_loc(1,:,:,:,1,1)
             ELSE
                IF (ALLOCATED(data4d)) DEALLOCATE(data4d)
                ALLOCATE(data4d(LBOUND(data_loc,2):UBOUND(data_loc,2), &
                     SIZE(data_loc,3),SIZE(data_loc,4), &
                     LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data4d(:,:,:,:) = data_loc(1,:,:,:,1,:)
             ENDIF
          CASE(3)  !no geographical dimensions
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data2d)) DEALLOCATE(data2d)
                ALLOCATE(data2d(SIZE(data_loc,3),SIZE(data_loc,4)))
                data2d(:,:) = data_loc(1,1,:,:,1,1)
             ELSE
                IF (ALLOCATED(data3d)) DEALLOCATE(data3d)
                ALLOCATE(data3d(SIZE(data_loc,3),SIZE(data_loc,4), &
                     LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data3d(:,:,:) = data_loc(1,1,:,:,1,:)
             ENDIF
          END SELECT
       ELSE
          SELECT CASE(rank)
          CASE(5)   
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data4d)) DEALLOCATE(data4d)
                ALLOCATE(data4d(SIZE(data_loc,1), &
                     LBOUND(data_loc,2):UBOUND(data_loc,2),SIZE(data_loc,3), &
                     SIZE(data_loc,5)))
                data4d(:,:,:,:) = data_loc(:,:,:,1,:,1)
             ELSE
                IF (ALLOCATED(data5d)) DEALLOCATE(data5d)
                ALLOCATE(data5d(SIZE(data_loc,1), &
                     LBOUND(data_loc,2):UBOUND(data_loc,2),SIZE(data_loc,3), &
                     SIZE(data_loc,5),LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data5d(:,:,:,:,:) = data_loc(:,:,:,1,:,:)
             ENDIF
          CASE(4)  !lev is missing
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data3d)) DEALLOCATE(data3d)
                ALLOCATE(data3d(SIZE(data_loc,1),SIZE(data_loc,3), &
                     SIZE(data_loc,5)))
                data3d(:,:,:) = data_loc(:,1,:,1,:,1)
             ELSE
                IF (ALLOCATED(data4d)) DEALLOCATE(data4d)
                ALLOCATE(data4d(SIZE(data_loc,1),SIZE(data_loc,3), &
                     SIZE(data_loc,5),LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data4d(:,:,:,:) = data_loc(:,1,:,1,:,:)
             ENDIF
          CASE(3)  !lon, lat are missing
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data2d)) DEALLOCATE(data2d)
                ALLOCATE(data2d(LBOUND(data_loc,2):UBOUND(data_loc,2), &
                     SIZE(data_loc,3)))
                data2d(:,:) = data_loc(1,:,:,1,1,1)
             ELSE
                IF (ALLOCATED(data3d)) DEALLOCATE(data3d)
                ALLOCATE(data3d(LBOUND(data_loc,2):UBOUND(data_loc,2), &
                     SIZE(data_loc,3),LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data3d(:,:,:) = data_loc(1,:,:,1,1,:)
             ENDIF
          CASE(2)  !no geographical dimensions
             IF (m0 .EQ. m1) THEN
                IF (ALLOCATED(data1d)) DEALLOCATE(data1d)
                ALLOCATE(data1d(SIZE(data_loc,3)))
                data1d(:) = data_loc(1,1,:,1,1,1)
             ELSE
                IF (ALLOCATED(data2d)) DEALLOCATE(data2d)
                ALLOCATE(data2d(SIZE(data_loc,3), &
                     LBOUND(data_loc,6):UBOUND(data_loc,6)))
                data2d(:,:) = data_loc(1,1,:,1,1,:)
             ENDIF
          END SELECT
       ENDIF
    ELSE
       SELECT CASE(rank)
       CASE(4)   
          IF (m0 .EQ. m1) THEN
             IF (ALLOCATED(data3d)) DEALLOCATE(data3d)
             ALLOCATE(data3d(SIZE(data_loc,1), &
                  LBOUND(data_loc,2):UBOUND(data_loc,2),SIZE(data_loc,5)))
             data3d(:,:,:) = data_loc(:,:,1,1,:,1)
          ELSE
             IF (ALLOCATED(data4d)) DEALLOCATE(data4d)
             ALLOCATE(data4d(SIZE(data_loc,1), &
                  LBOUND(data_loc,2):UBOUND(data_loc,2),SIZE(data_loc,5), &
                  LBOUND(data_loc,6):UBOUND(data_loc,6)))
             data4d(:,:,:,:) = data_loc(:,:,1,1,:,:)
          ENDIF
       CASE(3)  !lev is missing
          IF (m0 .EQ. m1) THEN
             IF (ALLOCATED(data2d)) DEALLOCATE(data2d)
             ALLOCATE(data2d(SIZE(data_loc,1),SIZE(data_loc,5)))
             data2d(:,:) = data_loc(:,1,1,1,:,1)
          ELSE
             IF (ALLOCATED(data3d)) DEALLOCATE(data3d)
             ALLOCATE(data3d(SIZE(data_loc,1),SIZE(data_loc,5), &
                  LBOUND(data_loc,6):UBOUND(data_loc,6)))
             data3d(:,:,:) = data_loc(:,1,1,1,:,:)
          ENDIF
       CASE(2)  !lon, lat are missing
          IF (m0 .EQ. m1) THEN
             IF (ALLOCATED(data1d)) DEALLOCATE(data1d)
             ALLOCATE(data1d(LBOUND(data_loc,2):UBOUND(data_loc,2)))
             data1d(:) = data_loc(1,:,1,1,1,1)
          ELSE
             IF (ALLOCATED(data2d)) DEALLOCATE(data2d)
             ALLOCATE(data2d(LBOUND(data_loc,2):UBOUND(data_loc,2), &
                  LBOUND(data_loc,6):UBOUND(data_loc,6)))
             data2d(:,:) = data_loc(1,:,1,1,1,:)
          ENDIF
       END SELECT
    ENDIF

    ! Set <levdimnoecham>, <nlevdimnoecham> if necessary:
    IF (PRESENT(nlevdimnoecham)) nlevdimnoecham=nlevdimnoechamloc
    IF (PRESENT(levdimnoecham)) THEN
       ALLOCATE(levdimnoecham(nlevdimnoecham))
       IF (p_parallel_io) levdimnoecham=levdimnoechamloc
       IF (p_parallel) CALL p_bcast (levdimnoecham, p_io)
    ENDIF

    ! Set <extradim>, <nextradim> if necessary:
    IF (PRESENT(nextradim)) nextradim=nextradimloc
    IF (PRESENT(extradim)) THEN
       ALLOCATE(extradim(nextradim))
       IF (p_parallel_io) extradim=extradimloc
       IF (p_parallel) CALL p_bcast (extradim, p_io)
    ENDIF

    ! Set <extradim2>, <nextradim2> if necessary:
    IF (PRESENT(nextradim2)) nextradim2=nextradim2loc
    IF (PRESENT(extradim2)) THEN
       ALLOCATE(extradim2(nextradim2))
       IF (p_parallel_io) extradim2=extradim2loc
       IF (p_parallel) CALL p_bcast (extradim2, p_io)
    ENDIF

    ! Deallocate memory:
    IF (p_parallel_io) THEN
       DEALLOCATE (zin)
       DEALLOCATE (zin_tr)
    ENDIF
    DEALLOCATE (data_loc)
    IF (ALLOCATED(data_loc4d)) DEALLOCATE (data_loc4d)
    IF (ALLOCATED(extradimloc)) DEALLOCATE(extradimloc)
    IF (ALLOCATED(extradimbefore)) DEALLOCATE(extradimbefore)
    IF (ALLOCATED(extradim2loc)) DEALLOCATE(extradim2loc)
    IF (ALLOCATED(extradim2before)) DEALLOCATE(extradim2before)
    IF (ALLOCATED(levdimnoechamloc)) DEALLOCATE(levdimnoechamloc)
    IF (ALLOCATED(levdimnoechambefore)) DEALLOCATE(levdimnoechambefore)

    ! Print message:
    IF (p_parallel_io) THEN

       IF (PRESENT(varname_longname)) THEN
          vn=TRIM(varname_longname)
       ELSE
          vn=TRIM(varname)
       ENDIF

          IF (vn .NE. '-') THEN
          
             IF (mo14loc) THEN
             IF (PRESENT(yr)) THEN
                WRITE(yr0_st, '(i4)') yr-1
                WRITE(yr1_st, '(i4)') yr
                WRITE(yr2_st, '(i4)') yr+1
                WRITE (message_text,*) &
                     'Reading ', TRIM(vn),' from files ', TRIM(fn(0)), ', ', &
                     TRIM(fn(1)), ', ', TRIM(fn(2))
                CALL message('',TRIM(message_text))
                WRITE (message_text,*) ' for year ', yr0_st, ', month 12', &
                     ', for year ', yr1_st, ', months 1-12, and for year ', &
                     yr2_st, ', month 1'
                CALL message('',TRIM(message_text))
             ELSE
                WRITE (message_text,*) &
                     'Reading ', TRIM(vn),' from file ', TRIM(fn(1)) , &
                     ' for year ', yr1_st
                CALL message('',TRIM(message_text))
             ENDIF
   
          ELSEIF (PRESENT(mo)) THEN
             IF (mo-1 .LE. 9) THEN
                WRITE(mo0_st, '(i1)') mo-1
             ELSE
                WRITE(mo0_st, '(i2)') mo-1
             ENDIF
             IF (mo .LE. 9) THEN
                WRITE(mo1_st, '(i1)') mo
             ELSE
                WRITE(mo1_st, '(i2)') mo
             ENDIF
             IF (mo+1 .LE. 9) THEN
                WRITE(mo2_st, '(i1)') mo+1
             ELSE
                WRITE(mo2_st, '(i2)') mo+1
             ENDIF
             IF (PRESENT(yr)) THEN
                WRITE(yr0_st, '(i4)') yr-1
                WRITE(yr1_st, '(i4)') yr
                WRITE(yr2_st, '(i4)') yr+1
                IF (mo .NE. 1 .AND. mo .NE. 12) THEN
                   WRITE (message_text,*) &
                        'Reading ', TRIM(vn),' from file ', TRIM(fn(1))
                   CALL message('',TRIM(message_text))
                   WRITE (message_text,*) ' for year ', yr1_st, ', months ', &
                        TRIM(mo0_st), ', ', TRIM(mo1_st), ', ', TRIM(mo2_st)
                ELSE
                   IF (mo .EQ. 1) THEN
                      WRITE (message_text,*) &
                           'Reading ', TRIM(vn),' from files ', TRIM(fn(0)), &
                           ', ', TRIM(fn(1))
                      CALL message('',TRIM(message_text))
                      WRITE (message_text,*) ' for year ', yr0_st, &
                           ', month 12, and for year ', yr1_st, ', months ', &
                           TRIM(mo1_st), ', ', TRIM(mo2_st)
                   ELSE
                      WRITE (message_text,*) &
                           'Reading ', TRIM(vn),' from files ', TRIM(fn(1)), &
                           ', ', TRIM(fn(2))
                      CALL message('',TRIM(message_text))
                      WRITE (message_text,*) ' for year ', yr1_st, &
                           ', months ', TRIM(mo0_st), ', ', TRIM(mo1_st), &
                           ', and for year ', yr2_st, ', month 1'
                   ENDIF
                ENDIF
                CALL message('',TRIM(message_text))
             ELSE
                WRITE (message_text,*) &
                     'Reading yearly constant ', TRIM(vn),' from file ', &
                     TRIM(fn(1))
                CALL message('',TRIM(message_text))
                IF (mo .NE. 1 .AND. mo .NE. 12) THEN
                   WRITE (message_text,*) ' for months ', TRIM(mo0_st), ', ', &
                        TRIM(mo1_st), ', ', TRIM(mo2_st)
                ELSE
                   IF (mo .EQ. 1) THEN
                      WRITE (message_text,*) ' for months ', TRIM(mo1_st), &
                           ', ', TRIM(mo2_st), ', 12'
                   ELSE
                      WRITE (message_text,*) ' for months 1, ', TRIM(mo0_st), &
                           ', ', TRIM(mo1_st)
                   ENDIF
                ENDIF
                CALL message('',TRIM(message_text))
             ENDIF
   
          ELSE
             WRITE (message_text,*) &
                  'Reading ', TRIM(vn),' from file ', TRIM(fn(1))      
             CALL message('',TRIM(message_text))
          ENDIF

       ENDIF

    ENDIF

  CONTAINS

    SUBROUTINE check_dim_netcdf(dimname)
      
      ! Checks for the dimension <dimname> if the data grid agrees with the one 
      ! of ECHAM5. Besides <lowlevind> / <uplevind>, <lowlevindloc> / 
      ! <uplevindloc> (in case of a pressure grid) and <extradim>, <nextradim>, 
      ! <nextradimloc> (in case of dimension <extradim>), <extradim2>, 
      ! <nextradim2>, <nextradim2loc> (in case of dimension <extradim2>) are 
      ! set. 

      ! Subroutine arguments:    
      CHARACTER(*), INTENT(in) :: dimname
  
      ! Local variables:
      INTEGER               :: io_var_id, ndim_f, dimerror
      REAL(dp), ALLOCATABLE :: dim_f(:)
      INTEGER, ALLOCATABLE  :: dimerrorv(:)
      LOGICAL, ALLOCATABLE  :: ll(:)
      CHARACTER(10)         :: dimname_e
      CHARACTER(12)         :: long_dimname
      INTEGER, SAVE         :: lowlevindbefore, uplevindbefore, &
           nextradimbefore, nextradim2before

      ! Executable statements:
      
      ! Calculate pressure at model levels for standard surface pressure:
      IF (llevdim) CALL calculate_pres_stand_socol

      IF (p_parallel_io) THEN
         ! Read data grid of dimension <dimname>:
         CALL IO_inq_dimid  (socolnc%file_id, dimname, io_var_id)
         CALL IO_inq_dimlen (socolnc%file_id, io_var_id, ndim_f)

         ! Allocate memory:
         ALLOCATE(dim_f(ndim_f))

         CALL IO_inq_varid (socolnc%file_id, dimname, io_var_id)
         CALL IO_get_var_double (socolnc%file_id, io_var_id, dim_f)
         
         ! Check if grids agree:
         dimname_e = dimname
         IF (TRIM(dimname) .EQ. 'lon' .OR. TRIM(dimname) .EQ. 'LON') THEN
            dimname_e='lon'
         ELSEIF (TRIM(dimname) .EQ. 'lat' .OR. TRIM(dimname) .EQ. 'LAT') THEN
            dimname_e='lat' 
         ELSEIF (TRIM(dimname) .EQ. 'lev' .OR. TRIM(dimname) .EQ. 'LEV' .OR. &
              TRIM(dimname) .EQ. 'p' .OR. TRIM(dimname) .EQ. 'P') THEN
            IF (.NOT. llevnoechamloc) THEN
               dimname_e='lev'
            ELSE
               dimname_e='levnoecham'
            ENDIF
         ELSEIF (TRIM(dimname) .EQ. TRIM(extradimname)) THEN
            dimname_e='extradim'
         ELSEIF (TRIM(dimname) .EQ. TRIM(extradim2name)) THEN
            dimname_e='extradim2'
         ENDIF

         ! Check dimensions and determine lowlevindloc and uplevindloc for 
         ! vertical grids:
         IF (dimname_e .NE. 'extradim' .AND. dimname_e .NE. 'extradim2' .AND. &
              dimname_e .NE. 'levnoecham') THEN
            CALL check_dim_socol(dimname_e, dim_f, dimerror, &
                 lowlevind=lowlevindloc, uplevind=uplevindloc)
         ELSE
            ! Initialization value for dimerror if no call to *check_dim_socol*:
            dimerror=0 
         ENDIF

         IF (dimname_e .EQ. 'lev') THEN
            ! Set <lowlevind> / <uplevind>:
            IF (f .EQ. f0) THEN
               IF (PRESENT(lowlevind)) lowlevind=lowlevindloc
               IF (PRESENT(uplevind)) uplevind=uplevindloc
            ENDIF
            IF (f .GT. f0) THEN
               IF (lowlevindbefore .NE. lowlevindloc .OR. &
                    uplevindbefore .NE. uplevindloc) dimerror=-100
            ENDIF
            lowlevindbefore=lowlevindloc
            uplevindbefore=uplevindloc
         ELSEIF (dimname_e .EQ. 'levnoecham') THEN
            ! Set nlevdimloc and levdimloc:
            nlevdimnoechamloc=ndim_f
            lowlevindloc=1
            uplevindloc=nlevdimnoechamloc
            IF (.NOT. ALLOCATED(levdimnoechamloc)) &
                 ALLOCATE(levdimnoechamloc(nlevdimnoechamloc))
            IF (.NOT. ALLOCATED(levdimnoechambefore)) &
                 ALLOCATE(levdimnoechambefore(nlevdimnoechamloc))
            levdimnoechamloc=dim_f
            IF (f .GT. f0) THEN
               IF (uplevindbefore .NE. nlevdimnoechamloc) dimerror=-100
               IF (dimerror .NE. -100) THEN
                  ALLOCATE(ll(nlevdimnoechamloc))
                  ALLOCATE(dimerrorv(nlevdimnoechamloc))
                  ll(:) = (levdimnoechamloc(:) .NE. levdimnoechambefore(:))
                  dimerrorv(:) = MERGE(-100,0,ll(:))
                  dimerror = SUM(dimerrorv)
               ENDIF
            ENDIF
            uplevindbefore=nlevdimnoechamloc
            levdimnoechambefore=levdimnoechamloc
         ENDIF

         ! Set extradimloc and nextradimloc:
         IF (dimname_e .EQ. 'extradim') THEN
            nextradimloc=ndim_f
            IF (.NOT. ALLOCATED(extradimloc)) &
                 ALLOCATE(extradimloc(nextradimloc))
            IF (.NOT. ALLOCATED(extradimbefore)) &
                 ALLOCATE(extradimbefore(nextradimloc))
            extradimloc=dim_f
            IF (f .GT. f0) THEN
               IF (nextradimbefore .NE. nextradimloc) dimerror=-1
               IF (dimerror .NE. -1) THEN
                  ALLOCATE(ll(nextradimloc))
                  ALLOCATE(dimerrorv(nextradimloc))
                  ll(:) = (extradimloc(:) .NE. extradimbefore(:))
                  dimerrorv(:) = MERGE(-1,0,ll(:))
                  dimerror = SUM(dimerrorv)
               ENDIF
            ENDIF
            nextradimbefore=nextradimloc
            extradimbefore=extradimloc
         ENDIF

         ! Set extradim2loc and nextradim2loc:
         IF (dimname_e .EQ. 'extradim2') THEN
            nextradim2loc=ndim_f
            IF (.NOT. ALLOCATED(extradim2loc)) &
                 ALLOCATE(extradim2loc(nextradim2loc))
            IF (.NOT. ALLOCATED(extradim2before)) &
                 ALLOCATE(extradim2before(nextradim2loc))
            extradim2loc=dim_f
            IF (f .GT. f0) THEN
               IF (nextradim2before .NE. nextradim2loc) dimerror=-1
               IF (dimerror .NE. -1) THEN
                  ALLOCATE(ll(nextradim2loc))
                  ALLOCATE(dimerrorv(nextradim2loc))
                  ll(:) = (extradim2loc(:) .NE. extradim2before(:))
                  dimerrorv(:) = MERGE(-1,0,ll(:))
                  dimerror = SUM(dimerrorv)
               ENDIF
            ENDIF
            nextradim2before=nextradim2loc
            extradim2before=extradim2loc
         ENDIF

         ! Deallocate memory:
         DEALLOCATE (dim_f)
         IF (ALLOCATED(ll)) DEALLOCATE(ll)
         IF (ALLOCATED(dimerrorv)) DEALLOCATE(dimerrorv)
      ENDIF

      ! Scatter dimname_e, dimerror, lowlevindloc, uplevindloc, 
      ! nlevdimnoechamloc, nextradimloc and nextradim2loc over processors:
      IF (p_parallel) THEN
         CALL p_bcast (dimname_e, p_io)
         CALL p_bcast (dimerror, p_io)
         IF (f .EQ. f0) THEN
            IF (dimname_e .EQ. 'lev' .OR. dimname_e .EQ. 'levnoecham') THEN
               CALL p_bcast (lowlevindloc, p_io)
               CALL p_bcast (uplevindloc, p_io)
               CALL p_bcast (nlevdimnoechamloc, p_io)
               IF (PRESENT(lowlevind)) CALL p_bcast (lowlevind, p_io)
               IF (PRESENT(uplevind)) CALL p_bcast (uplevind, p_io)
            ENDIF
            IF (dimname_e .EQ. 'extradim') CALL p_bcast (nextradimloc, p_io)
            IF (dimname_e .EQ. 'extradim2') CALL p_bcast (nextradim2loc, p_io)
         ENDIF
      ENDIF

      ! Terminate if dimerror .NE. 0:
      IF (dimerror .GE. 1) THEN
         IF (dimerror .EQ. 9999) THEN
            WRITE(message_text,*) 'Invalid dimension'
         ELSE
            SELECT CASE (dimname_e)
            CASE ('lon')
               long_dimname = 'Longitudinal'
            CASE ('lat')
               long_dimname = 'Latitudinal'
            CASE ('lev')
               long_dimname = 'Pressure'
            ENDSELECT
            
            WRITE(message_text,*) TRIM(long_dimname), ' grid of file <', &
                 TRIM(fn(f)),'> does not agree with the one of ECHAM5'
         ENDIF

         CALL message('',TRIM(message_text))
         CALL finish('socol_read_netcdf','Run terminated')
      ENDIF

      IF (dimerror .LE. -100) THEN
         WRITE(message_text,*) 'Pressure grid of file <', TRIM(fn(f)), &
              '> is not identical to the one of file <', TRIM(fn(f-1)), '>'
         CALL message('',TRIM(message_text))
         CALL finish('socol_read_netcdf','Run terminated') 
      ENDIF

      IF (dimerror .LE. -1 .AND. dimerror .GT. -100) THEN
         WRITE(message_text,*) 'Pressure grid of file <', TRIM(fn(f)), &
              '> is not identical to the one of file <', TRIM(fn(f-1)), '>'
         CALL message('',TRIM(message_text))
         CALL finish('socol_read_netcdf','Run terminated') 
      ENDIF

    END SUBROUTINE check_dim_netcdf

  END SUBROUTINE socol_read_netcdf

! ###########################################################################################
  FUNCTION socol_ascii_ionpair(yr, filename)

    ! Returns daily data values for year <yr> from the
    ! ascii file <filename>, where the data in the file
    ! are grouped as follows:
    !
    ! Pressh  Day1  Day2  Day3  Day4  Day5  Day6
    ! Press   Day1 ...
    ! Press
    ! Pressh  Day7   ..
    !  ..      ..    ..
    !
    ! Author: Julien G?rard Anet ETH Zuerich, 09/2010
    ! Template from Martin Schraner

    IMPLICIT NONE

    ! Function arguments:
    INTEGER, INTENT(in) :: yr
    CHARACTER(*), INTENT(in) :: filename

    ! Function output:
    REAL(dp), DIMENSION(59,59,367) :: socol_ascii_ionpair ! Dimension s...(:,1,1) is for pressi() !!!

    ! Local variables:
    INTEGER :: ioerror, i, j, Y, leapyear
    INTEGER, PARAMETER :: lun = 20
    REAL :: pressi(58), ionpair_y(58,366)

    Y = yr

IF (p_parallel_io == .TRUE.) THEN
    ! Executable statements:

    ioerror = 0

     ! Open file:
       OPEN (lun, ERR=110, FILE=filename, STATUS='old', &
            FORM='formatted', ACTION='read')
       GOTO 111
110    ioerror = 1
111    CONTINUE

    IF(MOD(Y,100).NE.0.AND.MOD(Y,4).EQ.0) THEN
       leapyear=1
    ELSEIF(MOD(Y,400).EQ.0) THEN
       leapyear=1
    ELSE
       leapyear=0
    ENDIF

       ! Read file:
     IF (leapyear==1) THEN
       IF (ioerror .EQ. 0) THEN
         READ(lun,'(1a)')
         READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
             DO j=1,58
                 READ(lun,'(4x,e11.7,4x,366(e12.5,4x))') pressi(59-j),ionpair_y(59-j,1:366)
             ENDDO
        ENDIF
      ELSEIF (leapyear==0) THEN
        IF (ioerror .EQ. 0) THEN
         READ(lun,'(1a)')
         READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
             DO j=1,58
                 READ(lun,'(4x,e11.7,4x,365(e12.5,4x))') pressi(59-j),ionpair_y(59-j,1:365)
             ENDDO
        ENDIF
       ENDIF

120      CONTINUE

 DO i=2,59
       socol_ascii_ionpair(i,1,1) = pressi(i-1)
    ENDDO
    DO i=2,59
       DO j=2,367
           socol_ascii_ionpair(1,i,j) = ionpair_y(i-1,j-1)
        ENDDO
    ENDDO
       ! Close file:
    CLOSE (lun)

ELSE  ! Not Parallel!

    ! Executable statements:

    ioerror = 0

     ! Open file:
       OPEN (lun, ERR=510, FILE=filename, STATUS='old', &
            FORM='formatted', ACTION='read')
       GOTO 511
510    ioerror = 1
511    CONTINUE

    IF(MOD(Y,100).NE.0.AND.MOD(Y,4).EQ.0) THEN
       leapyear=1
    ELSEIF(MOD(Y,400).EQ.0) THEN
       leapyear=1
    ELSE
       leapyear=0
    ENDIF

       ! Read file:
     IF (leapyear==1) THEN
       IF (ioerror .EQ. 0) THEN
         READ(lun,'(1a)')
         READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
             DO j=1,58
                 READ(lun,'(4x,e11.7,4x,366(e12.5,4x))') pressi(59-j),ionpair_y(59-j,1:366)
             ENDDO
        ENDIF
      ELSEIF (leapyear==0) THEN
        IF (ioerror .EQ. 0) THEN
         READ(lun,'(1a)')
         READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
	 READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
         READ(lun,'(1a)')
             DO j=1,58
                 READ(lun,'(4x,e11.7,4x,365(e12.5,4x))') pressi(59-j),ionpair_y(59-j,1:365)
             ENDDO
        ENDIF
       ENDIF

130      CONTINUE

 DO i=2,59
       socol_ascii_ionpair(i,1,1) = pressi(i-1)
    ENDDO
    DO i=2,59
       DO j=2,367
           socol_ascii_ionpair(1,i,j) = ionpair_y(i-1,j-1)
        ENDDO
    ENDDO
       ! Close file:
       CLOSE (lun)

ENDIF

    ! Terminate if ioerror > 0:
    IF (p_parallel) THEN
       CALL p_bcast (ioerror, p_io)
    ELSE
       CALL p_bcast (ioerror,p_io)
    ENDIF
    SELECT CASE (ioerror)
    CASE(1)
       WRITE(message_text,*) 'Could not open ', TRIM(filename)
       CALL message('',TRIM(message_text))
       CALL finish('socol_ascii_ionpair','Run terminated')
    END SELECT

    ! Distribute data:

    IF (p_parallel) THEN
       CALL p_bcast (socol_ascii_ionpair, p_io)
    ELSE
       CALL p_bcast (socol_ascii_ionpair, p_io)
    ENDIF


RETURN

  END FUNCTION socol_ascii_ionpair
  ! ###########################################################################################
  FUNCTION socol_ascii_geomag(yr, filename)

    ! Returns daily data values for year <yr> from the
    ! ascii file <filename>, where the data in the file
    ! are grouped as follows:
    ! YYYY MAGLAT MAGLON MAGSTRENTH
    ! ...
    !
    ! Author: Julien G?rard Anet ETH Zuerich, 09/2010
    ! Template from Martin Schraner

    IMPLICIT NONE

    ! Function arguments:
    INTEGER, INTENT(in) :: yr
    CHARACTER(*), INTENT(in) :: filename

    ! Function output:
    REAL(dp) :: socol_ascii_geomag(3)

    ! Local variables:
    INTEGER :: ioerror, i, j, year0, year
    INTEGER, PARAMETER :: lun = 87
    REAL :: maglat(1), maglon(1), magstrength(1)

    ioerror=0

    IF (p_parallel_io) THEN


       ! Open file:
       OPEN (lun, ERR=810, FILE=filename, STATUS='unknown')
       GOTO 811
810    ioerror = 1
811    CONTINUE

       ! Read file:
       IF (ioerror .EQ. 0) THEN
          year0 = yr
          DO
             READ (lun,*, END=820, ERR=820) year, maglat, maglon, magstrength
             IF (year .EQ. year0) THEN
                EXIT
             ENDIF
          ENDDO
820       CONTINUE
          IF (year .NE. year0) ioerror=2
       ENDIF

       ! Close file:
       CLOSE (lun)

    ELSE ! NOT Parallel_IO


       ! Open file:
       OPEN (lun, ERR=890, FILE=filename, STATUS='unknown')
       GOTO 812
890    ioerror = 1
812    CONTINUE

       ! Read file:
       IF (ioerror .EQ. 0) THEN
          year0 = yr
          DO
             READ (lun,*, END=840, ERR=840) &
                  year, maglat, maglon, magstrength
             IF (year .EQ. year0) THEN
                EXIT
             ENDIF
          ENDDO
840       CONTINUE
          IF (year .NE. year0) ioerror=2
       ENDIF

       ! Close file:
       CLOSE (lun)

    ENDIF

   socol_ascii_geomag(1) = maglat(1)
   socol_ascii_geomag(2) = maglon(1)
   socol_ascii_geomag(3) = magstrength(1)

    ! Terminate if ioerror > 0:
    IF (p_parallel) THEN
       CALL p_bcast (ioerror, p_io)
    ELSE
       CALL p_bcast (ioerror,p_io)
    ENDIF
    SELECT CASE (ioerror)
    CASE(1)
       WRITE(message_text,*) 'Could not open ', TRIM(filename)
       CALL message('',TRIM(message_text))
       CALL finish('socol_ascii_geomag','Run terminated')
    END SELECT

    ! Distribute data:

    IF (p_parallel) THEN
       CALL p_bcast (socol_ascii_geomag, p_io)
    ELSE
       CALL p_bcast (socol_ascii_geomag, p_io)
    ENDIF


RETURN

  END FUNCTION socol_ascii_geomag
! ###########################################################################################
 FUNCTION socol_ascii_ap(yr,filename)

    ! Returns yearly data values for year <yrs> to year <yre>
    ! from the ascii file <filename>, where the data in the
    ! file are grouped as follows:
    !
    ! Apyears
    ! Apyears+1
    ! ...
    ! Apyeare
    !
    ! Author: Julien G?rard Anet ETH Zuerich, 09/2010
    ! Template from Martin Schraner

     IMPLICIT NONE

    INTEGER, INTENT(in) :: yr
    CHARACTER(*), INTENT(in) :: filename

    ! Function output:
    REAL(dp), DIMENSION(1) :: socol_ascii_ap

    ! Local variables:
    INTEGER :: ioerror, year, year0
    INTEGER, PARAMETER :: lun = 20
    CHARACTER(4) :: yr_st

    ! Executable statements:

    ioerror = 0

    IF (p_parallel_io) THEN

       ! Open file:
       OPEN (lun, ERR=110, FILE=filename, STATUS='unknown')
       GOTO 111
110    ioerror = 1
111    CONTINUE

       ! Read file:
       IF (ioerror .EQ. 0) THEN
          year0 = yr
          DO
             READ (lun,*, END=120, ERR=120) &
                  year, socol_ascii_ap

             IF (year .EQ. year0) THEN
                EXIT
             ENDIF
          ENDDO
120       CONTINUE
          IF (year .NE. year0) ioerror=2
       ENDIF

       ! Close file:
       CLOSE (lun)

    ELSE


       ! Open file:
       OPEN (lun, ERR=190, FILE=filename, STATUS='unknown')
       GOTO 112
190    ioerror = 1
112    CONTINUE

       ! Read file:
       IF (ioerror .EQ. 0) THEN
          year0 = yr
          DO
             READ (lun,*, END=140, ERR=140) &
                  year, socol_ascii_ap

             IF (year .EQ. year0) THEN
                EXIT
             ENDIF
          ENDDO
140       CONTINUE
          IF (year .NE. year0) ioerror=2
       ENDIF

       ! Close file:
       CLOSE (lun)

    ENDIF


    ! Terminate if ioerror > 0:
    IF (p_parallel) CALL p_bcast (ioerror, p_io)
    SELECT CASE (ioerror)
    CASE(1)
       WRITE(message_text,*) 'Could not open ', TRIM(filename)
       CALL message('',TRIM(message_text))
       CALL finish('socol_ascii_ap','Run terminated')
    CASE(2)
       WRITE(yr_st, '(i4)') year0
       WRITE(message_text,*) &
            'Year ', yr_st, ' not found in file <', &
            TRIM(filename),'>'
       CALL message('',TRIM(message_text))
       CALL finish ('socol_ascii_ap', 'Run terminated.')
    END SELECT

    ! Distribute data:
    IF (p_parallel) THEN
       CALL p_bcast (socol_ascii_ap, p_io)
    ELSE
       CALL p_bcast (socol_ascii_ap, p_io)
    ENDIF


  END FUNCTION socol_ascii_ap
! ###########################################################################################
 FUNCTION socol_read_ionization_netcdf_1d(filename1, arrayd1, a, p, g, iru0, pt)

 IMPLICIT NONE

    ! Function input arguments:
    CHARACTER(*), INTENT(in)           :: filename1
    INTEGER, INTENT(in)      :: arrayd1
    LOGICAL, OPTIONAL, INTENT(in) :: a, p, g, iru0, pt

    ! Function output
    REAL(dp) :: socol_read_ionization_netcdf_1d(arrayd1)

    ! Local variables:
    TYPE (FILE_INFO) :: socolnc
    CHARACTER(10	) :: xdimname, ydimname, zdimname
    CHARACTER(32) :: fn(0:2)
    CHARACTER(40) :: vn2,vn3,vn4,vn5
    INTEGER :: i,j,k,rankf,rank,f,lowlevindloc,uplevindloc,start(6),count(6),moyr
    LOGICAL :: lex

    ! Initial value for open/close status of socolnc:
    socolnc%opened = .FALSE.
    vn2='phitime'
    vn3='phi'
    vn4='alt'
    vn5='geomagnetic_cutoff'

    ! Set rank for ir_u0
    rank=4

    ! Determine start(1:rank-1) (used for reading netcdf) and rankf:
    IF (p_parallel) THEN
       IF (p_parallel_io) THEN
          start(1:rank-1) = 1
          rankf=rank
       ENDIF
    ELSE
       start(1:rank-1) = 1
       rankf=rank
    ENDIF

    ! Determine filenames:
       fn(1) = TRIM(filename1)
       fn(0) = fn(1)
       fn(2) = fn(1)

    ! Level heights
    lowlevindloc=1
    ! Loop over files:
    f = 1
    IF (p_parallel) THEN
       IF (p_parallel_io) INQUIRE (FILE=TRIM(fn(f)), EXIST=lex)

       ! Terminate if file does not exist:
       IF (p_parallel) CALL p_bcast (lex, p_io)
       IF (.NOT. lex) THEN
          WRITE (message_text,*) 'Could not open file <', TRIM(fn(f)), '>'
          CALL message('',TRIM(message_text))
          CALL finish ('socol_read_netcdf', 'Run terminated.')
       ENDIF

       ! Open file:
       IF (p_parallel_io) CALL IO_open (TRIM(fn(f)), socolnc, IO_READ)

       IF (p_parallel_io) THEN
          ! Read data:
          start(rank)=1
          rankf=1
          count(1) = 6000
          count(2) = 0
          count(3) = 0
          IF (PRESENT(pt)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn2), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF
          count(1) = 31
          IF (PRESENT(p)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn3), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF
          count(1) = 131
          IF (PRESENT(a)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn4), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF
          count(1) = 77
          IF (PRESENT(g)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn5), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF

          ! Close file:
          CALL IO_close(socolnc)
       ENDIF
    ELSE ! Not parallel

       INQUIRE (FILE=TRIM(fn(f)), EXIST=lex)

       ! Terminate if file does not exist:
       CALL p_bcast (lex, p_io)
       IF (.NOT. lex) THEN
          WRITE (message_text,*) 'Could not open file <', TRIM(fn(f)), '>'
          CALL message('',TRIM(message_text))
          CALL finish ('socol_read_netcdf', 'Run terminated.')
       ENDIF

       ! Open file:
      CALL IO_open (TRIM(fn(f)), socolnc, IO_READ)


          ! Read data:
          start(rank)=1
          rankf=1
          count(1) = arrayd1
          count(2) = 0
          count(3) = 0
          IF (PRESENT(pt)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn2), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF
          count(1) = 31
          IF (PRESENT(p)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn3), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF
          count(1) = 131
          IF (PRESENT(a)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn4), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF
          count(1) = 77
          IF (PRESENT(g)) THEN
             CALL IO_inq_varid (socolnc%file_id, TRIM(vn5), io_var_id)
             CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:rankf), &
                  count(1:rankf), socol_read_ionization_netcdf_1d(:))
          ENDIF

          ! Close file:
          CALL IO_close(socolnc)
ENDIF


    ! Distribute data:
    IF (p_parallel) THEN
       CALL p_bcast (socol_read_ionization_netcdf_1d, p_io)
    ELSE
       CALL p_bcast (socol_read_ionization_netcdf_1d, p_io)
    ENDIF

 END FUNCTION
! ###########################################################################################
FUNCTION socol_read_ionization_netcdf_3d(filename1, arrayd1, arrayd2, arrayd3)

 IMPLICIT NONE

    ! Function input arguments:
    CHARACTER(*), INTENT(in)           :: filename1
    INTEGER, INTENT(in)      :: arrayd1, arrayd2, arrayd3

    ! Function output
    REAL(dp), DIMENSION(arrayd1,arrayd2,arrayd3) :: socol_read_ionization_netcdf_3d

    ! Local variables:
    TYPE (FILE_INFO) :: socolnc
    CHARACTER(10) :: xdimname, ydimname, zdimname
    CHARACTER(32) :: fn(0:2)
    CHARACTER(40) :: vn1
    INTEGER :: i,j,k,rankf,rank,f,lowlevindloc,uplevindloc,start(6),count(6),moyr
    LOGICAL :: lex

    ! Initial value for open/close status of socolnc:
    socolnc%opened = .FALSE.
    vn1='ionization_rate'

    ! Set rank for ir_u0
    rank=4

    ! Determine start(1:rank-1) (used for reading netcdf) and rankf:
     IF (p_parallel) THEN
       IF (p_parallel_io) THEN
          start(1:rank-1) = 1
          rankf=rank
       ENDIF
     ELSE
       start(1:rank-1) = 1
       rankf=rank
     ENDIF

    ! Determine filenames:
       fn(1) = TRIM(filename1)
       fn(0) = fn(1)
       fn(2) = fn(1)

    ! Level heights
    lowlevindloc=1
    ! Loop over files:
    f = 1
IF (p_parallel) THEN
    IF (p_parallel_io) INQUIRE (FILE=TRIM(fn(f)), EXIST=lex)

       ! Terminate if file does not exist:
       IF (p_parallel) CALL p_bcast (lex, p_io)
       IF (.NOT. lex) THEN
          WRITE (message_text,*) 'Could not open file <', TRIM(fn(f)), '>'
          CALL message('',TRIM(message_text))
          CALL finish ('socol_read_netcdf', 'Run terminated.')
       ENDIF

       ! Open file:
       IF (p_parallel_io) CALL IO_open (TRIM(fn(f)), socolnc, IO_READ)

       IF (p_parallel_io) THEN
          ! Read data:
          start(rank)=1
          count(1) = 77
          count(2) = 131
          count(3) = 31

          CALL IO_inq_varid (socolnc%file_id, TRIM(vn1), io_var_id)
          CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:3), &
                  count(1:3), socol_read_ionization_netcdf_3d(:,:,:))

          ! Close file:
          CALL IO_close(socolnc)
       ENDIF
ELSE ! Not parallel
   INQUIRE (FILE=TRIM(fn(f)), EXIST=lex)

       ! Terminate if file does not exist:
       CALL p_bcast (lex, p_io)
       IF (.NOT. lex) THEN
          WRITE (message_text,*) 'Could not open file <', TRIM(fn(f)), '>'
          CALL message('',TRIM(message_text))
          CALL finish ('socol_read_netcdf', 'Run terminated.')
       ENDIF

       ! Open file:
       CALL IO_open (TRIM(fn(f)), socolnc, IO_READ)

          ! Read data:
          start(rank)=1
          count(1) = 77
          count(2) = 131
          count(3) = 31

          CALL IO_inq_varid (socolnc%file_id, TRIM(vn1), io_var_id)
          CALL IO_get_vara_double (socolnc%file_id, io_var_id, start(1:3), &
                  count(1:3), socol_read_ionization_netcdf_3d(:,:,:))

          ! Close file:
          CALL IO_close(socolnc)
ENDIF
    ! Distribute data:
    IF (p_parallel) THEN
       CALL p_bcast (socol_read_ionization_netcdf_3d, p_io)
    ELSE
       CALL p_bcast (socol_read_ionization_netcdf_3d, p_io)
    ENDIF

 END FUNCTION

END MODULE mo_socol_readfile

MODULE mo_ohclim

  !- Description:
  !
  !  This module provides a 3 dimensional monthly mean OH climatology.
  !  Module is based on mo_o3clim
  !
  !  This module contains:
  !
  !  A) internal variables
  !  B) the subroutine su_ohclim to read the OH file and
  !     initialize the OH climatology
  !  C) the function ohclim to get a longitude height ozone section
  !  D) the subroutine cleanup_ohclim to deallocate memory
  !
  !- Author:
  !
  !  A. Stenke, ETHZ, Feb 2010
  !
  !=======================================================================

  USE mo_kind,       ONLY: dp
  USE mo_interpo,    ONLY: wgt1,wgt2,nmw1,nmw2,nmw1cl,nmw2cl
  USE mo_radiation,  ONLY: nmonth

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: su_ohclim, ohclim, pre_ohclim, cleanup_ohclim
  PUBLIC :: zohc_x

  !=======================================================================
  !
  ! A)
  !
  ! unit number of ozone file
  INTEGER, PARAMETER :: nioh=22

  ! pressure level grid of ozone climatology in Pa
  INTEGER            :: noh          ! number of ozone full levels
  REAL(dp), ALLOCATABLE  :: pfoh(:)      ! full levels 
  REAL(dp), ALLOCATABLE  :: phoh(:)      ! half levels

  ! ozone climatology in mass mixing ratio in layers, latitudes and months
  REAL(dp), ALLOCATABLE  :: ohcli(:,:,:,:) ! mass mixing ratio at full levels
  REAL(dp), ALLOCATABLE, TARGET  :: zohc_x(:,:,:)

  !=======================================================================

CONTAINS

  !=======================================================================
  !
  ! B)

  SUBROUTINE su_ohclim

    !- Description:
    !
    !  su_ohclim provides a monthly mean
    !  climatology of OH volume mixing ratio (mol/mol).
    !
    !  The NetCDF file accessed by unit nioh must contain ozone
    !  volume mixing ratios in ppmv at pressure levels.
    !
    !  The pressure levels are ordered from top to surface.
    !  At each level the data are given for ngl latitudes
    !  and all months.
    !
    !  The arrays are extended in time such that:
    !  month 0 corresponds to month 12 and month 13 to month 1
    !  in order to facilitate the interpolation in time.
    !
    !  Reading the NetCDF file is based on readozone in mo_midatm.f90
    !  (version ECHAM5.0.03) of U. Schulzweida
    !
    !- Author:
    !
    !  A. Stenke, ETHZ, Feb 2010
    !

    USE mo_control,       ONLY: ngl, nlon
    USE mo_doctor,        ONLY: nout, nerr
    USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
         &                      io_var_id
    USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
         &                      io_inq_varid, io_get_var_double, &
         &                      io_get_vara_double, FILE_INFO
    USE mo_decomposition, ONLY: lc => local_decomposition, gl_dc => global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_filename,      ONLY: NETCDF  

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:,:), zin2(:,:,:,:)
    REAL(dp), POINTER :: gl_oh(:,:,:,:)

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER               :: jk ,i     ! loop index

    INTEGER :: IO_file_id
    INTEGER        :: iy
    TYPE (FILE_INFO) :: fohclim

    ! Read ozone file
    ! ===============

    IF (p_pe==p_io) THEN

      fohclim%opened = .FALSE.

       WRITE(nout,'(/,A,I2)') ' Read OH climatology from unit ', nioh

       fohclim%format = NETCDF
       CALL io_open_unit (nioh, fohclim, io_read)
       io_file_id = fohclim%file_id

      ! Check resolution
      CALL io_inq_dimid  (io_file_id, 'lat', io_var_id)
      CALL io_inq_dimlen (io_file_id, io_var_id, io_ngl)
      CALL io_inq_dimid  (io_file_id, 'lon', io_var_id)
      CALL io_inq_dimlen (io_file_id, io_var_id, io_nlon)
      CALL io_inq_dimid  (io_file_id, 'lev_p', io_var_id)
      CALL io_inq_dimlen (io_file_id, io_var_id, noh)

      IF (io_ngl/=ngl) THEN
        WRITE(nerr,*) 'su_ohclim: unexpected resolution ',io_nlon,io_ngl
        WRITE(nerr,*) 'expected number of latitudes = ',ngl
        WRITE(nerr,*) 'number of latitudes of OH climatology = ',io_ngl
        CALL finish ('su_ohclim','unexpected resolution')
      END IF
       
      IF (io_nlon == 1) THEN
        WRITE(nerr,*) 'OH data has zonal mean climatology'
        CALL finish ('su_ohclim','OH zonal mean')
      END IF
       
    END IF

    CALL p_bcast (noh, p_io)

    ! ALLOCATE memory
    IF ( .NOT. ALLOCATED(pfoh))  ALLOCATE(pfoh(noh))
    IF ( .NOT. ALLOCATED(phoh))  ALLOCATE(phoh(noh+1))

    !     Allocate memory for oh per PE

    IF ( .NOT. ALLOCATED(ohcli)) ALLOCATE(ohcli(lc%nproma,noh,lc%ngpblks,0:13))

    IF (p_pe == p_io) THEN

       ALLOCATE (zin(nlon,noh,ngl,1:12))
       ALLOCATE (zin2(nlon,noh,ngl,0:13))

       ! read full pressure levels of OH climatology
       CALL io_inq_varid      (io_file_id, 'lev_p', io_var_id)
       CALL io_get_var_double (io_file_id, io_var_id, pfoh(:))

       WRITE(nout,*) 'number of pressure levels of OH climatology: ',noh
       WRITE(nout,*) 'pressure levels in [hPa]:'
       WRITE(nout,'(38F8.2)') pfoh(:)
       
       ! hPa -> Pa
       pfoh(:) = pfoh(:) * 100.
       
       ! read OH volume mixing ratio at full pressure levels in mol/mol
       ! read 12 month
       CALL io_inq_varid (io_file_id, 'OHMR', io_var_id)
       CALL io_get_var_double(io_file_id,io_var_id,zin)
       
       CALL io_close (fohclim)
      
       WRITE(nout,*) 'OH climatology read successfully'

       DO i = 0, 13
          IF ( i == 0) THEN
             zin2(:,:,:,i) = zin(:,:,:,12)
          ELSEIF (i == 13 ) THEN
             zin2(:,:,:,i) = zin(:,:,:,1)
          ELSE
             zin2(:,:,:,i) = zin(:,:,:,i)
          END IF
       END DO
      
    END IF
   
    NULLIFY (gl_oh)

    DO i = 0, 13
       IF (p_pe == p_io) gl_oh => zin2(:,:,:,i:i)
       CALL scatter_gp(gl_oh, ohcli(:,:,:,i:i), gl_dc)
    END DO
    
    IF (p_pe == p_io) THEN
       DEALLOCATE (zin)
       DEALLOCATE (zin2)
    END IF

    CALL p_bcast (pfoh, p_io)


    ! define half levels of ozone pressure grid
    ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
    ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
    phoh(1)=60._dp
    phoh(2:noh)=(pfoh(1:noh-1)+pfoh(2:noh))/2._dp
    phoh(noh+1)=125000._dp
    
  END SUBROUTINE su_ohclim


  !=======================================================================
  !
  ! C)

  FUNCTION ohclim(krow,kproma,kbdim,klev,pph,ppf)

    !- Description:
    !
    !  The OH monthly mean climatology is interpolated
    !  to the actual time of the integration and integrated from p=0 to
    !  the surface p=ps.
    !  The time interpolated profile of OH is interpolated to the
    !  model full levels and integrated again from p=60. to p=ps.
    !  Finally the OH profile on the model levels is normalized such
    !  that the integrated amount on the model grid is identical to
    !  that on the grid of the climatology.
    !
    !- Author:
    !
    !  A. Stenke, ETHZ, Feb 2010

    ! INPUT
    ! -----

    INTEGER, INTENT(in)                      :: krow   ! local latitude index
    INTEGER, INTENT(in)                      :: kproma ! number of local longitudes
    INTEGER, INTENT(in)                      :: kbdim  ! first dimension of 2-d arrays
    INTEGER, INTENT(in)                      :: klev   ! number of levels
    REAL(dp), INTENT(in),DIMENSION(kbdim,klev)   :: ppf  ! full level pressure
    REAL(dp), INTENT(in),DIMENSION(kbdim,klev+1) :: pph  ! half level pressure

    ! OUTPUT
    ! ------
    REAL(dp), DIMENSION(kproma,klev)              :: ohclim ! OH in mol/mol


    ! LOCAL
    ! -----

    ! pressure integrated ozone at half levels of ozone grid,
    ! and integral at surface
    REAL(dp), DIMENSION(kproma)                 :: zohintc

    ! time interpolated ozone at full levels of model grid,
    ! pressure integrated ozone at half levels of model,
    ! and integral at surface
    REAL(dp), DIMENSION(kproma,klev)            :: zohm !kk
    REAL(dp), DIMENSION(kproma)                 :: zohintm

    REAL(dp) :: zdp1,zdp2 ! pressure weights in linear interpolation
    INTEGER  :: jl,jk,jkk ! loop indices
    INTEGER,DIMENSION(kproma)  :: jk1,jkn   ! first and last model level in interpolation

    INTEGER,DIMENSION(kproma)            :: kwork
    LOGICAL,DIMENSION(kproma)            :: kk_flag


!kk NEC version                                      Klaus Ketelsen
!kk   1. Remove EXIT from loop for vectorizing
!kk   2. DO jl=1,kproma as inner loop to get sufficiant vector length

!kk Now some loops have higher operation cout, but they are vectorizing
!kk The logic is the same as in the original version using EXIT in loops


    ! interpolate ozone profile to model grid
    ! ---------------------------------------
    ! set OH concentration at levels above the uppermost level of
    ! the OH climatology to the value in the uppermost level of 
    ! the OH climatology

    jk1(:)     = 1
    kk_flag(:) = .TRUE.
    DO jk = 1,klev
       DO jl=1,kproma
          IF (ppf(jl,jk)<=pfoh(1) .AND. kk_flag(jl)) THEN 
             zohm(jl,jk)=zohc_x(jl,1,krow)
             jk1(jl)=jk+1
          ELSE
             kk_flag(jl) = .FALSE.
          END IF
       END DO
    END DO

    ! set OH concentration at levels below the lowermost level of
    ! the OH climatology to the value in the lowermost level of 
    ! the OH climatology
    jkn(:)=klev
    kk_flag(:) = .TRUE.
    DO jk = klev,1,-1
       DO jl=1,kproma
          IF (ppf(jl,jk)>=pfoh(noh).AND. kk_flag(jl)) THEN
             zohm(jl,jk)=zohc_x(jl,noh,krow)
             jkn(jl)=jk-1
          ELSE
             kk_flag(jl) = .FALSE.
          END IF
       END DO
    ENDDO

    DO jk=1,klev
       kk_flag(:) = .TRUE.
       kwork(:)   = 1
       DO jkk = 1,noh
          DO jl=1,kproma
             IF(jk >= jk1(jl) .AND. jk <= jkn(jl))  THEN
                IF (ppf(jl,jk) <= pfoh(jkk) .AND. jkk >= kwork(jl) &
                                            .AND. kk_flag(jl)) THEN
                   kwork(jl)   = jkk
                   kk_flag(jl) = .FALSE.
                END IF
             END IF
          END DO
       END DO

       DO jl=1,kproma
          IF(jk >= jk1(jl) .AND. jk <= jkn(jl))  THEN
                jkk = kwork(jl)
                ! model level is in interval ]pfoh(jkk-1),pfoh(jkk)]
                ! -> make interpolation
                zdp1=pfoh(jkk)-ppf(jl,jk)
                zdp2=ppf(jl,jk)-pfoh(jkk-1)
                zohm(jl,jk)=(zdp1*zohc_x(jl,jkk-1,krow) &
                      +zdp2*zohc_x(jl,jkk,krow))/ &
                     &      (zdp1+zdp2)
          END IF
       END DO
    END DO

       ! integrate ozone profile on grid of climatology
       ! from top to surface
       ! ----------------------------------------------
    zohintc=0._dp
    kk_flag(:) = .TRUE.
    jk1(:)     = 2
    DO jk=2,noh+1
          ! integrate layers of climatology above surface
       DO jl=1,kproma
          IF (phoh(jk)<=pph(jl,klev) .AND. kk_flag(jl) ) THEN
              zohintc(jl)=zohintc(jl)+ &
               &  zohc_x(jl,jk-1,krow)*(phoh(jk)-phoh(jk-1))
              jk1(jl) = jk+1
          ELSE
              kk_flag(jl) = .FALSE.
          END IF
       END DO
    END DO
       ! integrate layer of climatology that is intersected
       ! by the surface from upper boundary to surface
    DO jl=1,kproma
       zohintc(jl)=zohintc(jl)+ &
            &  zohc_x(jl,jk1(jl)-1,krow)*(pph(jl,klev)-phoh(jk1(jl)-1))
    END DO

       ! integrate ozone profile on grid of model
       ! from top to surface
       ! ----------------------------------------
    zohintm=0._dp
    DO jk=2,klev+1
       DO jl=1,kproma
         zohintm(jl)=zohintm(jl) + zohm(jl,jk-1)*(pph(jl,jk)-pph(jl,jk-1))
       END DO
    END DO

       ! normalize interpolated ozone profile such that the
       ! ozone integral computed on the model grid is equal
       ! to that integrated on the grid of the climatology
       ! --------------------------------------------------
    DO jk=1,klev
       DO jl=1,kproma
         ohclim(jl,jk)=zohm(jl,jk)/zohintm(jl) * zohintc(jl)
         if (ohclim(jl,jk) .lt. 0._dp) ohclim(jl,jk) = 0._dp
       END DO
    END DO


  END FUNCTION ohclim

  !=======================================================================

  SUBROUTINE pre_ohclim

    USE mo_decomposition,  ONLY: ldc => local_decomposition

    INTEGER  :: jl, jk, jrow
    INTEGER  :: ngpblks, nproma

    !  Executable statements

    ngpblks = ldc% ngpblks

    IF ( .NOT. ALLOCATED(zohc_x)) ALLOCATE(zohc_x(ldc%nproma,noh,ldc%ngpblks))

    DO jrow = 1, ngpblks

      IF ( jrow == ldc% ngpblks ) THEN
        nproma = ldc% npromz
      ELSE
        nproma = ldc% nproma
      END IF

      DO jk=1,noh
        DO jl=1,nproma
          zohc_x(jl,jk,jrow) = wgt1*ohcli(jl,jk,jrow,nmw1)+wgt2*ohcli(jl,jk,jrow,nmw2)
        END DO
      END DO

    END DO

  END SUBROUTINE pre_ohclim
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_ohclim
    !----------------------------
    ! deallocate module variables
    !----------------------------
    IF (ALLOCATED(zohc_x))    DEALLOCATE (zohc_x)
    IF (ALLOCATED(ohcli))     DEALLOCATE (ohcli)
    IF (ALLOCATED(pfoh))      DEALLOCATE (pfoh)
    IF (ALLOCATED(phoh))      DEALLOCATE (phoh)

  END SUBROUTINE cleanup_ohclim
!------------------------------------------------------------------------------
END MODULE mo_ohclim


      MODULE mo_socol_ionization

!+
!
!       Summary :
!          Calculation of SPE/EEP/GCR flux from the modeltop
!          downwards. Leads to NO/N/OH "creation"
!       -> EEP flux Parameterization based on Baumgaertner (2008)
!          uses Ap averaged for May-July for 1976-2005
!       -> SPE from Jackman (M Calisto for SOCOLv2)
!       -> GCR from Usoskin (M Calisto for SOCOLv2)
!       -> 17.5.2013:: If SPE is not true, GCR did not produce any NO
!          This bug has been corrected now
!
!         *mo_socol_eep_spe* is called from *socol_mezon*
!
!          Also includes an interpolation subroutine and a
!          geomagnetic transposer, to get the geomagnetic
!          latitude coordinate on the earth.
!
!       Authors :
!          Marco Calisto ETH/IAC, 2009
!          Julien Anet ETH/IAC 2010/2011
!          Ivo Suter, ETH/IAC 2013

    USE mo_socol_readfile,            ONLY: socol_ascii_ap, socol_ascii_ionpair, &
                                             socol_read_ionization_netcdf_1d, socol_read_netcdf, &
                                             socol_read_ionization_netcdf_3d, socol_ascii_geomag
    USE mo_kind,                      ONLY: dp
    USE mo_socol_namelist,            ONLY: spe, eep, gcr

IMPLICIT NONE

      ! This is the lower N/S latitudinal boarder below with no ionization is possible for SPEs
      REAL, PARAMETER :: lat_border_spe=60.

      ! Declaration of the needed arrays to read in the ionization data for EPP
      REAL(dp), PUBLIC :: r_ion(58,367), pressi_loc(58), ap(1), &
                           gmc(77), phi_time(6000), phi(31), &
                           alt(131), ir_u0(77,131,31), phi_mod(1)
      REAL(dp), ALLOCATABLE, PUBLIC, DIMENSION(:,:,:)           :: ir_u2_loc, spen_loc
      REAL(dp), ALLOCATABLE, PUBLIC, DIMENSION(:,:)             :: mlat, ir_u1,  ap_flux_tot

      ! Target to repartition the ionization information globally (by scatter_gp)
      REAL(dp), ALLOCATABLE, PUBLIC, DIMENSION(:,:,:), TARGET  :: ir_u2, spen
      REAL(dp), ALLOCATABLE, PUBLIC, DIMENSION(:,:), TARGET    :: ap_flux_tot_loc

      ! Variables to save the information about the state of the Earths magnetic field
      REAL, PUBLIC ::  mplat, mplon, magstrength

      ! Declaration of subroutines
      PUBLIC :: read_socol_ionization_data,  socol_ionization_repartition_month, &
                 socol_ionization_repartition_day, socol_gcr_spe_ionization

CONTAINS

! ###########################################################################################
    SUBROUTINE read_socol_ionization_data(yr,mo,dy,newmonth)

    ! read_socol_ionization_data is called from *socol_read_bcond_m and _d*
    ! Needs to be called every day for SPE. For GCR / EEP every month.
    ! Calls mo_socol_readfile in order to read in the needed data from either
    ! NetCDF (ionization_info) or ASCII (ap-index, geomagetic data etc...)

    USE mo_time_control,              ONLY: dt_start, dt_stop, next_date, get_year_day, get_date_components
    USE mo_exception,                 ONLY: message_text, message
    USE mo_socol_grid_calculations,   ONLY: calculate_pres_stand_socol
    USE mo_mpi,                       ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe

    IMPLICIT NONE

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr,mo,dy
    LOGICAL, INTENT(in) :: newmonth

    ! Local variables:
    INTEGER, PARAMETER :: nsw_data=6
    INTEGER :: iyr, imo, day, ihr, imn, ise
    REAL(dp) :: ionization_info(59,59,367),geomag_data(3)
    CHARACTER(20) :: fn
    INTEGER :: i,j,timerange,timeap
    CHARACTER(LEN=100) :: year

 ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PARALLEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF (p_parallel) THEN

       IF (newmonth) THEN

          IF (spe) THEN

             ! Allocate variables
             alt(:) = 0.0_dp
             phi(:) = 0.0_dp
             gmc(:) = 0.0_dp
             phi_time(:)=0.0_dp
             ir_u0(:,:,:)=0.0_dp

             ! Writes Integer2String
             WRITE(year, '(i8)') yr
             year = TRIM(adjustl(year))

             ! Read in the ionization rate from the ASCII tables
             ionization_info(:,:,:) = socol_ascii_ionpair(yr,''//TRIM(year)//'_ionrate-nitrate.txt')
             WRITE (message_text,*) &
                  'Reading file ',TRIM(year),'_ionrate-nitrate.txt'
             CALL message('',TRIM(message_text))

             ! Split up the variable, future name of the array: SPEN(:,:)
             DO i=2,59
                pressi_loc(i-1) = ionization_info(i,1,1)
             ENDDO
             DO i=2,59
                DO j=2,367
                   r_ion(i-1,j-1) = ionization_info(1,i,j)
                ENDDO
             ENDDO

          ENDIF ! No SPE calculation
          
       ENDIF ! Not present day

       ! Read in geomagnetic latitude/longitude/strength of field
       geomag_data(:) = socol_ascii_geomag(yr,'geomag_data_1600_2100')
       WRITE (message_text,*) &
            'Reading file geomag_data_1600_2100'
       CALL message('',TRIM(message_text))

       ! Split up the variable
       mplat = geomag_data(1)
       mplon = geomag_data(2)
       magstrength = geomag_data(3)
       
       WRITE (message_text,*) 'MPLAT, MPLON and Strength for year ', yr ,' is: ', mplat, mplon, magstrength
       CALL message ('',TRIM(message_text))

       ! Get DOY
       CALL get_date_components(next_date,iyr,imo,day,ihr,imn,ise)
       day = get_year_day(next_date)
       
       IF (.NOT.(newmonth)) THEN

          IF (spe) THEN
             CALL socol_ionization_repartition_day(yr,mo,day,mplat,mplon,pressi_loc,r_ion)
          ELSEIF (eep) THEN ! We do only repartition for EEP
             CALL socol_ionization_repartition_day(yr,mo,day,mplat,mplon)
          ENDIF
          
       ENDIF ! Day present
    
       IF (newmonth) THEN
          
          ! Read in the auroral index from the table
          IF (eep) THEN
             ap(:)=socol_ascii_ap(yr,'monthly_mean_ap')
             WRITE (message_text,*) &
                  'Read file monthly_mean_ap for year ', yr, '. Ap index is ', ap
             CALL message('',TRIM(message_text))
          ENDIF ! EEP on
          IF (gcr) THEN
             
             ! Read NetCDF File (GCR informations, originally from EAWAG, F. Steinhilber)
             fn = 'gcr_data'
             ! Ionization rates
             ir_u0(:,:,:) = socol_read_ionization_netcdf_3d(TRIM(fn), 77, 131, 31)
             ! Altitude level information of the ionization rates
             alt(:) = socol_read_ionization_netcdf_1d(TRIM(fn), 131, a=.TRUE.)
             ! Range of the possible solar modulation potential
             phi(:) = socol_read_ionization_netcdf_1d(TRIM(fn), 31, p=.TRUE.)
             ! Range of the geomagnetic activity
             gmc(:) = socol_read_ionization_netcdf_1d(TRIM(fn), 77, g=.TRUE.)
             ! Range of the number of months available in the NetCDF file
             phi_time(:) =  socol_read_ionization_netcdf_1d(TRIM(fn), 6000, pt=.TRUE.)
             
             ! Correcting the Phi-LIS to Phi-US
             DO i = 1,6000
                phi_time(i) = 1._dp/1.04 * (phi_time(i)+73._dp)
                IF (phi_time(i) .LT. 0 ) THEN
                   phi_time(i) = 0
                ENDIF
             ENDDO
             
             ! Repartition of the ionization information over the globe
             CALL socol_ionization_repartition_month(yr,mo,phi,gmc,ir_u0,mplat,mplon,magstrength)
             
          ENDIF ! GCR on
          
       ENDIF ! newmonth
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOT PARALLEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ELSE ! Not Parallel

    WRITE(*,*) 'Init of EPP is not called parallely. Is this wanted?'

    IF (newmonth) THEN
       IF (spe) THEN

          ! Allocate variables
          alt(:) = 0.0_dp
          phi(:) = 0.0_dp
          gmc(:) = 0.0_dp
          phi_time(:)=0.0_dp
          ir_u0(:,:,:)=0.0_dp
          
          ! Writes Integer2String
          WRITE(year, '(i8)') yr
          year = TRIM(adjustl(year))
          
          ! Read in the ionization rate from the tables
          ionization_info(:,:,:) = socol_ascii_ionpair(yr,''//TRIM(year)//'_ionrate-nitrate.txt')
          WRITE(*,*) 'Exited readin'
          WRITE (message_text,*) &
               'Reading file ',TRIM(year),'_ionrate-nitrate.txt'
          CALL message('',TRIM(message_text))
          
          ! Split up the variable, future name of the array is SPEN(:,:)
          DO i=2,59
             pressi_loc(i-1) = ionization_info(i,1,1)
          ENDDO
          DO i=2,59
             DO j=2,367
                r_ion(i-1,j-1) = ionization_info(1,i,j)
             ENDDO
          ENDDO
       ENDIF ! No SPE calculation
    ENDIF ! Not present day
    
    ! Read in geomagnetic latitude/longitude/strength of field
    geomag_data(:) = socol_ascii_geomag(yr,'geomag_data_1600_2100')
    WRITE (message_text,*) &
         'Reading file geomag_data_1600_2100'
    CALL message('',TRIM(message_text))
    
    ! Split up the variable into magnetic lat/lon & strength
    mplat = geomag_data(1)
    mplon = geomag_data(2)
    magstrength = geomag_data(3)

    WRITE (message_text,*) 'MPLAT, MPLON and Strength for year ', yr ,' is: ', mplat, mplon, magstrength
    CALL message ('',TRIM(message_text))
    
    
    IF (.NOT.(newmonth)) THEN
       
       IF (spe) THEN
          CALL socol_ionization_repartition_day(yr,mo,day,mplat,mplon,pressi_loc,r_ion)
       ELSE ! We do only repartition for EEP
          CALL socol_ionization_repartition_day(yr,mo,day,mplat,mplon)
       ENDIF
    ENDIF ! Day present
    
    IF (newmonth) THEN
       
       ! Read in the auroral index from the table
       IF (eep) THEN
          ! Read in the auroral index from the ASCII table
          ap(:)=socol_ascii_ap(yr,'monthly_mean_ap')
          WRITE (message_text,*) &
               'Reading file monthly_mean_ap for year ', yr, '. Ap index is ', ap
          CALL message('',TRIM(message_text))
       ENDIF ! EEP on
       
       IF (gcr) THEN
          ! Read NetCDF File (GCR informations, originally from EAWAG, F. Steinhilber)
          fn = 'gcr_data'
          ! Ionization rates
          ir_u0(:,:,:) = socol_read_ionization_netcdf_3d(TRIM(fn), 77, 131, 31)
          ! Altitude level information of the ionization rates
          alt(:) = socol_read_ionization_netcdf_1d(TRIM(fn), 131, a=.TRUE.)
          ! Range of the possible solar modulation potential
          phi(:) = socol_read_ionization_netcdf_1d(TRIM(fn), 31, p=.TRUE.)
          ! Range of the geomagnetic activity
          gmc(:) = socol_read_ionization_netcdf_1d(TRIM(fn), 77, g=.TRUE.)
          ! Range of the number of months available in the NetCDF file
          phi_time(:) =  socol_read_ionization_netcdf_1d(TRIM(fn), 6000, pt=.TRUE.)
          
          ! Correcting the Phi-LIS to Phi-US
          DO i = 1,6000
             phi_time(i) = 1._dp/1.04 * (phi_time(i)+73._dp)
             IF (phi_time(i) .LT. 0 ) THEN
                phi_time(i) = 0
             ENDIF
          ENDDO

          ! Calling the repartition routine which interpolated the data over the globe
          CALL socol_ionization_repartition_month(yr,mo,phi,gmc,ir_u0,mplat,mplon,magstrength)
          
       ENDIF ! GCR on
    ENDIF ! No newmonth
 ENDIF ! Not parallel


    END SUBROUTINE read_socol_ionization_data
! ###########################################################################################
    SUBROUTINE socol_ionization_repartition_month(yr,mo,phi,gmc,ir_u0,mplat,mplon,magstrength)


!          Description :
!          This routine interpolates the GCR ionization rates
!
!      1.  Over the model altitude levels
!      2.  Over the globe (lat & lon)
!      3.  With respect of the geomagnetic latitude

       ! Use of Pi
       USE mo_constants,                 ONLY: api

       ! Needed for the density information of the atmosphere
       USE mo_socol_grid_calculations,   ONLY: dens, zetb

       ! How many lon/lat/lev do we have? As well, izcalc gives the information about
       ! the starting year of the simulation / first year of the GCR NetCDF file
       ! Should your GCR file start e.g. in year 1000, please change IZCALC_START in
       ! the namelist of SOCOL appropriately!
       USE mo_control,                   ONLY: nlon, ngl, nlev, izcalc_start, nproma

       ! Decomposition of the Latitude / Longitude over the number of cores
       USE mo_decomposition,             ONLY: lc => local_decomposition, dcg => global_decomposition

       ! Which doy do we have?
       USE mo_socol_time_control,        ONLY: ndays_currentyear

       ! Scatter routine to repartition the data over the globe
       USE mo_transpose,                 ONLY: scatter_gp

       ! Reading the pressure information of every model level
       USE mo_socol_gcmfields,           ONLY: prest ! (nbdim,nlev)

       ! MPI parallelization routines
       USE mo_mpi,                       ONLY: p_parallel, p_parallel_io, p_io, p_pe


       IMPLICIT NONE
!*********************************************************************
       ! Variables needed for GCRs
       REAL(dp), INTENT(in)  :: phi(31), gmc(77), ir_u0(77,131,31)
       REAL, INTENT(in)      :: mplat, mplon, magstrength

      ! Size of phi, gmc and alt dimensions of your NetCDF GCR file (CHANGE appropriately)
      INTEGER, PARAMETER :: n_phi=31, n_gmc=77, n_alt=131

      ! output of ionization rate for the model for any time step
      REAL(dp) :: ir_socol(nlon,ngl,nlev)

      ! GMC array of geomagnetic latitudes
      REAL(dp) :: gmc_mod(nlon,ngl)

      ! Pointer & Target for scattering
      REAL(dp), POINTER :: ir_u2_point(:,:,:)
      REAL(dp), TARGET :: ir_u2_regrid(nlon,n_alt,ngl)

      ! Pressure height arrays temporary needed
      REAL(dp), TARGET :: prest_talc(nproma,nlev),prestglob(nlon,nlev,ngl)
      REAL(dp), ALLOCATABLE :: prest_gl(:,:,:)
! ********************************************************************
      ! Loop indexes
      INTEGER   :: i,j,k

      ! General integers for time calculation or else
      INTEGER :: yr,mo
      INTEGER :: timeidx, output, deltat

! ********************************************************************

      ! Get magnetic latitudes and longitudes for GCR
      IF (.NOT. ALLOCATED(ir_u1)) &
       ALLOCATE(ir_u1(n_gmc,n_alt))
      IF (.NOT. ALLOCATED(ir_u2)) &
         ALLOCATE(ir_u2(nlon,ngl,n_alt))
      IF (.NOT. ALLOCATED(prest_gl)) &
         ALLOCATE(prest_gl(nlon,nlev,ngl))

      ! We put in mplat and mplon we get back the magn. lat (mlat)
      CALL geo2mag(mplat,mplon)

      ! Zero-setting of the GCR ionization rates arrays
      ir_u1(:,:) = 0.0_dp
      ir_u2(:,:,:) = 0.0_dp

      ! Current GM rigidity in Cb * m^2 * s^-1 (Resistance)
      DO i = 1, nlon
         DO j = 1, ngl
            gmc_mod(i,j) = 1.9*(magstrength/1e22)*(cos(mlat(i,j)/180.*api))**4
         ENDDO
      ENDDO

      ! Interpolate the GCR ionization rates over the "phi-time"
      ! (REMEMBER: IZCALC_START is set to first year of your GCR NetCDF-File!)
      timeidx = (yr-izcalc_start)*12+mo
      phi_mod = 1/1.04 * (phi_time(timeidx)+73)

      IF (p_pe == p_io) THEN
         WRITE(*,*) 'Phi for timeindex ',timeidx, ' is ',phi_mod
      ENDIF
      ! Interpolation of ir_u0(phi 'time') in our model time steps -> ir_u1(modeltimestep)
      DO i = 1,n_gmc
         DO j = 1,n_alt
            CALL line(n_phi,phi,ir_u0(i,j,:),1,phi_mod,ir_u1(i,j))
         ENDDO
      ENDDO

      ! Interpolate the GCR ionization rates over our *!!*horizontal*!!* model grid
      DO k = 1, n_alt
         DO i = 1, nlon
            DO j = 1, ngl
               ir_u2(i,j,k) = 0.
               CALL line(n_gmc,gmc,ir_u1(:,k),1,gmc_mod(i,j),ir_u2(i,j,k))
            ENDDO
         ENDDO
      ENDDO

      ! Switch dimensions of the GCR ionization array
      DO i = 1,nlon
          DO j = 1,n_alt
             DO k = 1,ngl
                ir_u2_regrid(i,j,k) = 0.0_dp
                ir_u2_regrid(i,j,k) = ir_u2(i,k,j)
             ENDDO
          ENDDO
       ENDDO

   ! Prepare for scattering of GCR fields
   IF (p_parallel .AND. .NOT. ALLOCATED(ir_u2_loc)) &
      ALLOCATE(ir_u2_loc(lc%nproma,n_alt,lc%ngpblks))
   IF (.NOT. p_parallel .AND. .NOT. ALLOCATED(ir_u2_loc)) &
      ALLOCATE(ir_u2_loc(lc%nproma,n_alt,lc%ngpblks))

   DO i = 1,lc%nproma
      DO j = 1,n_alt
         DO k = 1,lc%ngpblks
            ir_u2_loc(i,j,k) = 0.0_dp
         ENDDO
      ENDDO
   ENDDO
   NULLIFY(ir_u2_point)

     	 IF (p_pe == p_io) ir_u2_point => ir_u2_regrid(:,:,:)
         CALL scatter_gp(ir_u2_point,ir_u2_loc(:,:,:),dcg)
   ! Clean up

      IF (ALLOCATED(ir_u1)) &
         DEALLOCATE(ir_u1)
      IF (ALLOCATED(ir_u2)) &
         DEALLOCATE(ir_u2)
      IF (ALLOCATED(prest_gl)) &
         DEALLOCATE(prest_gl)


    END SUBROUTINE socol_ionization_repartition_month
! ###########################################################################################
    SUBROUTINE socol_ionization_repartition_day(yr,mo,day,mplat,mplon,pressi_loc,r_ion)

    ! Interpolate EEP/SPE ionization rates to model pressure levels.
    ! Called by read_socol_ionization_data

    ! How many lon/lat/lev do we have?
    USE mo_control,                      ONLY: nlon, ngl, nlev, nproma

    ! Pressure levels in SOCOL
    USE mo_socol_grid_calculations,      ONLY: presf_stand_socol

    ! Scattering routines to scatter arrays over the globe
    USE mo_transpose,                    ONLY: scatter_gp
    USE mo_decomposition,                ONLY: lc => local_decomposition, global_decomposition

    ! Constant Pi
    USE mo_constants,                    ONLY: api

    ! Get calendar information
    USE mo_time_base,                    ONLY: Set_JulianDay, Set_Ly360Day
    USE mo_time_control,                 ONLY: get_clock,current_date, get_date_components, Get_Ly360YearDay, ly360_date
    USE mo_socol_time_control,           ONLY: chemistry_date_dyofyr

    ! MPI Parallelization routines
    USE mo_mpi,                          ONLY: p_parallel, p_parallel_io, p_io, p_pe, p_bcast

    IMPLICIT NONE
!*********************************************************************
    ! Variables used for calculations of ionization by SPE

    ! Input data which was read in before. Optional as SPE can be disabled
    REAL(dp), INTENT(in), OPTIONAL :: pressi_loc(58), r_ion(58,367)

    ! Lower N/S latitudinal boarder below which no SPE can ionize the atmosphere
    REAL, PARAMETER :: lat_border_spe=60.

    ! Temporary arrays to regrid the SPE data
    REAL(dp), TARGET  :: spen_regrid(nlon,nlev,ngl)
    REAL(dp), POINTER :: spen_point(:,:,:)
    REAL(dp)           :: spe_loc(nlev)

    ! Resized SPE data arrays
    REAL(dp) :: pressi(58), r_ion_rs(58,366)
!********************************************************************
   ! Variables for calculation for EEP flux by Baumg??rtner
    REAL, PARAMETER :: lat_border_eep=55.,c=0.23
    REAL(dp), TARGET :: ap_flux_tot_tar(nlon,ngl)
    REAL(dp), POINTER :: ap_flux_tot_point(:,:)

    ! Date & loop integers
    INTEGER, INTENT(in) :: mo,day,yr
    REAL, INTENT(in) :: mplat, mplon
    INTEGER :: m,dy,y,se,hr,mn,i,j,k,sh_flag,nday
    REAL(dp) :: zclock

    ! Ap-flux computation
    REAL(dp) :: ap_flux(nlon,ngl), max_cos(nlon,ngl)

    ! Tool for calendar recognition
    TYPE(ly360_date) :: ly360_date_now
!********************************************************************

   ! Allocate SPE arrays
   IF (spe) THEN
       IF (.NOT. ALLOCATED(spen)) &
          ALLOCATE(spen(nlon,ngl,nlev))
       IF (.NOT. ALLOCATED(spen_loc)) &
          ALLOCATE(spen_loc(lc%nproma,nlev,lc%ngpblks))
   ENDIF ! SPE True

   IF (eep) THEN
   ! Allocate AP-index arrays
   IF (.NOT. ALLOCATED(ap_flux_tot)) &
       ALLOCATE(ap_flux_tot(nlon,ngl))
   IF (.NOT. ALLOCATED(ap_flux_tot_loc)) &
       ALLOCATE(ap_flux_tot_loc(lc%nproma,lc%ngpblks))
   ENDIF ! EEP True

   ! Getting the geaomagnetic latitude
   CALL geo2mag(mplat,mplon)

   IF (spe) THEN
      DO i=1,58
         pressi(i) = pressi_loc(i)
         r_ion_rs(i,day) = r_ion(i,day)
      ENDDO

      ! Nullify arrays
      spe_loc(:) = 0.0_dp
      spen(:,:,:) = 0.0_dp

      ! Interpolate SPE ionization rates over our model pressure grid (vertical)
      CALL line(58,pressi,r_ion_rs(:,day),nlev,presf_stand_socol(:),spe_loc)

      ! Write ionization data only if we are above 60 N/S
      DO i = 1,nlon
         DO j = 1,ngl
            spen(i,j,:) = 0.0_dp
            IF (abs(mlat(i,j)) >= lat_border_spe) THEN ! Inside polar cup
               spen(i,j,:) = spe_loc(:)
            ENDIF
         ENDDO
      ENDDO

      DO i = 1,nlon
         DO j = 1,nlev
            DO k = 1,ngl
               spen_regrid(i,j,k) = 0.0_dp
               spen_regrid(i,j,k) = spen(i,k,j)
            ENDDO
         ENDDO
      ENDDO
   ENDIF ! SPE True

   IF (eep) THEN
      ! Flag to see if we are on the southern hemisphere or not
      sh_flag = ngl/2

      DO i = 1,nlon
         DO j = 1,ngl
            ! Nullify AP-index arrays
            ap_flux(i,j) = 0.0_dp
            ap_flux_tot(i,j) = 0.0_dp
            ap_flux_tot_tar(i,j) = 0.0_dp
            max_cos(i,j) = 0.0_dp
            ! Are we in the polar cup?
            IF(abs(mlat(i,j)) >= lat_border_eep) THEN
               ap_flux(i,j)=ap(1)**2.5*c*2.2*10**5    ! Baumgaertner 2009, Equation 4, Atmos. Chem. Phys., 9, 2729???2740, 2009
               ! Southern hemisphere?
               IF (j > sh_flag) THEN
                  max_cos(i,j)=max(0.1,cos(api/182.625*(chemistry_date_dyofyr-172.625)))!southern hemisphere because of 172.625
               ELSE
                  max_cos(i,j)=max(0.1,cos(api/182.625*(chemistry_date_dyofyr-355.25))) !northern hemisphere because of 355.25
               ENDIF

            ! For those interested: flux given in mol(NO)*cm-2*s-1
            ap_flux_tot(i,j)=ap_flux(i,j)*max_cos(i,j)

            ENDIF ! Polar cup
         ENDDO ! ngl
      ENDDO ! nlon
   ENDIF ! EEP True

   IF (spe) THEN
      NULLIFY(spen_point)
      IF(p_pe == p_io)  spen_point => spen_regrid(:,:,:)
      CALL scatter_gp(spen_point,spen_loc(:,:,:),global_decomposition)
   ENDIF ! SPE True

   IF (eep) THEN
      ap_flux_tot_tar = ap_flux_tot

      NULLIFY(ap_flux_tot_point)

      IF(p_pe == p_io)  ap_flux_tot_point => ap_flux_tot_tar(:,:)
      CALL scatter_gp(ap_flux_tot_point,ap_flux_tot_loc(:,:),global_decomposition)

   ENDIF ! EEP
   ! Clean up

    IF (ALLOCATED(spen)) &
         DEALLOCATE(spen)


    END SUBROUTINE socol_ionization_repartition_day
! ###########################################################################################
    SUBROUTINE socol_gcr_spe_ionization(krow,klev,kbdim,ktrac,prest,kproma,pxtte,pxtm1)


    ! Finally, after all those repartitions and interpolations, we are able to compute the net
    ! NO, N and OH flux out of the EPP

    ! Gather model level information, incl. pressure and density on model levels
    USE mo_socol_grid_calculations, ONLY: calculate_pres_stand_socol, zlevb, zetb, dens, presf_stand_socol

    ! Get the tracer arrays, which we will feed with the additional flux
    USE mo_socol_tracers,           ONLY: idt_no, idt_n, idt_oh
    USE mo_mpi,                     ONLY: p_pe, p_io

!*********************************************************************
    ! Input indexes
    INTEGER, INTENT(in) :: krow,klev,kbdim,kproma,ktrac
    REAL(dp), INTENT(in), OPTIONAL, DIMENSION(kproma,klev) :: prest

    ! This is going in, will be updated, and is going out
    REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)

    ! Delta-arrays to memorize the NO/N/OH production
    REAL(dp)   :: deltNO(krow,kproma,klev), deltN(krow,kproma,klev), deltOH(krow,kproma,klev)

    ! Indexes which "limit" the OH production induced by LEEP / GCR
    REAL(dp)   :: RXO2, RXO3, D1,D2,D3

    ! Helper arrays for NO production by LEE. Dz is for finding out the two upperst levels
    REAL(dp)   :: ap_NO(krow,kproma,klev), dz(1)

    ! Needed to know the altitude number of the GCR array, as well: Local ionization rate array
    INTEGER, PARAMETER :: n_alt = 131
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ir_socol_loc

    ! Loop indexes
    INTEGER :: i, j, k, iz, iz1, ii1, m

    ! Flux arrays of OH production by GCR/SPE
    REAL(dp) :: f_oh_gcr,f_oh_spe

    ! Data arrays, which tell the model how many HOxes is produced per ion pair
    REAL(dp) ::  t_hox(19,6), z_hox(19), i_hox(6)
    DATA t_hox(:,1) /2.,2.,2.,2.,2.,2.,2.,2.,2.,1.99,1.99,1.99,1.98,&
&                       1.98,1.96,1.93,1.84,1.6,.8/
    DATA t_hox(:,2) /2.,2.,2.,2.,2.,2.,2.,2.,2.,1.99,1.99,1.99,1.98,&
&                     1.98,1.96,1.93,1.84,1.6,.8/
    DATA t_hox (:,3) /2.,2.,2.,2.,2.,1.99, 1.99, 1.99, 1.99, 1.98,&
&                   1.98, 1.98, 1.96, 1.94, 1.9, 1.84, 1.72, 1.4,.6/
    DATA t_hox (:,4) /2.,2.,1.99,1.99,1.99,1.99,1.99,1.98,1.98,1.97,&
&                     1.96,1.94,1.91,1.87,1.8,1.73,1.6,1.2,0.4/
    DATA t_hox(:,5)/2.,1.99,1.99,1.99,1.99,1.98,1.98,1.97,1.95,1.94,&
&                     1.9,1.87,1.82,1.77,1.7,1.6,1.4,.95,.15/
    DATA t_hox (:,6) /2.00,1.99,1.99,1.99,1.98,1.98,1.95,1.93,1.89,&
&                     1.85,1.81,1.77,1.72,1.64,1.5,1.3,.93,.4,0./
    DATA z_hox /0.0,40.0,42.5,45.0,47.5,50.0,52.5,55.0,57.5,60.0,62.5,&
&                      65.0,67.5,70.0,72.5,75.0,77.5,80.0,100/
    DATA i_hox /0.,1e1,1e2,1e3,1e4,1e5/ ! Ionization levels
!*********************************************************************

    IF (gcr) THEN ! eth_as_03062013
       ! Allocate local GCR array
       IF (.NOT. ALLOCATED(ir_socol_loc)) &
            ALLOCATE(ir_socol_loc(kproma,klev))
    ENDIF ! GCR True

    ! Fetch model levels in pressure data
    CALL calculate_pres_stand_socol

    IF (gcr) THEN
      DO i = 1,kproma
         DO k = 1,klev
            ir_socol_loc(i,k)=0.0_dp
         ENDDO
         ! Interpolate one last time vertically over all pressure levels
         CALL line(n_alt,alt,ir_u2_loc(i,:,krow),klev,presf_stand_socol(:),ir_socol_loc(i,:))
         DO m = 1,klev
            ! Should there be one value over one billion, something has gone seriously wrong!
            IF (ir_socol_loc(i,m) >= 10E6) THEN
               WRITE(*,*) 'SEVERE ERROR, EXITING DUE TO HIGH GCR IONIZATION RATE'
               WRITE(*,*) 'ir_u2_loc at krow was ', ir_u2_loc(i,:,krow), krow
               WRITE(*,*) 'ir_socol_loc at krow and nlev/nproma was ', ir_socol_loc(i,m), krow, i, m
               WRITE(*,*) 'presf_stand_socol(m) is ', presf_stand_socol(m)
               STOP
            ENDIF ! No error
         ENDDO ! Check over all levels
      ENDDO ! kproma
    ENDIF ! GCR True

    ! Initialize arrays
    deltNO(:,:,:) = 0.0_dp
    deltOH(:,:,:) = 0.0_dp
    deltN(:,:,:) = 0.0_dp

     DO i = 1,kproma
       DO k = 1,klev
         IF (gcr) THEN
            ir_socol_loc(i,k)=ir_socol_loc(i,k)*dens(i,k)*29.*1.661*1e-24 !convert /g/s
            ! Caution, there might be an error
            IF (ir_socol_loc(i,k) .GT. 10E6) THEN
               WRITE(*,*) 'WARNING: Ionization too high, kproma klev is', i, k, ir_socol_loc(i,k), dens(i,k)
            ENDIF
            DO iz = 1,19
               IF(z_hox(iz) < zlevb(i,k) ) THEN
                  iz1 = iz
               ENDIF ! ionization of hox
            ENDDO ! iz

         ! HOx production by GCRs
            DO iz = 1,5
               IF(i_hox(iz) <= ir_socol_loc(i,k) ) THEN
                   ii1 = iz
               ENDIF
            ENDDO

            RXO2 = ((z_hox(iz1+1)-z_hox(iz1))/(zlevb(i,k)-z_hox(iz1)))
            RXO3 = ((i_hox(ii1+1)-i_hox(ii1))/(ir_socol_loc(i,k)-i_hox(ii1)))
            D1 = t_hox(iz1,ii1) + (t_hox(iz1+1,ii1)-t_hox(iz1,ii1))/ RXO2
            D2 = t_hox(iz1,ii1+1) + (t_hox(iz1+1,ii1+1)-t_hox(iz1,ii1+1))/RXO2
            D3 = D1 + (D2 - D1)/RXO3
            f_oh_gcr =  D3
          ELSE ! GCR True
            f_oh_gcr = 0.0_dp
          ENDIF ! GCR false

          ! HOx production by SPEs
          IF (spe) THEN

     		  ! For SPE

              DO iz = 1,19
               IF(z_hox(iz) < zlevb(i,k) ) THEN
                  iz1 = iz
               ENDIF ! ionization of hox
              ENDDO ! iz

              DO iz = 1,5
                 IF(i_hox(iz) <= spen_loc(i,k,krow)) THEN
                    ii1 = iz
                 ENDIF
              ENDDO

              RXO2 = ((z_hox(iz1+1)-z_hox(iz1))/(zlevb(i,k)-z_hox(iz1)))
              RXO3 = ((i_hox(ii1+1)-i_hox(ii1))/(spen_loc(i,k,krow)-i_hox(ii1)))
              D1 = t_hox(iz1,ii1) + (t_hox(iz1+1,ii1)-t_hox(iz1,ii1))/ RXO2
              D2 = t_hox(iz1,ii1+1) + (t_hox(iz1+1,ii1+1)-t_hox(iz1,ii1+1))/RXO2
              D3 = D1 + (D2 - D1)/RXO3
              f_oh_spe =  D3                                                   ! GCR -> SPE

              ! Something has gone seriously wrong
              IF (spen_loc(i,k,krow) >= 1E6) THEN
                 WRITE(*,*) 'SPE ionization rate FAILURE                     ', spen_loc(i,k,krow)
                 WRITE(*,*) 'kproma klev krow were                           ', i, k, krow
                 STOP
              ENDIF

              ! Compute additional NO now.
              deltNO(krow,i,k)  =  (0.7*spen_loc(i,k,krow))&
                        /dens(i,k)
              deltN(krow,i,k)   =  (0.55*spen_loc(i,k,krow))&
                         /dens(i,k)
          ELSE ! SPE False
             f_oh_spe = 0.0_dp
          ENDIF !SPE True

          ! Add the additional OH now.
          IF (gcr) THEN
             IF (spe) THEN
                deltNO(krow,i,k)  =  (0.7*(ir_socol_loc(i,k)+spen_loc(i,k,krow))&
                        /dens(i,k))
                deltN(krow,i,k)   =  (0.55*(ir_socol_loc(i,k)+spen_loc(i,k,krow))&
                         /dens(i,k))
             ELSE
                deltNO(krow,i,k)  =  (0.7*ir_socol_loc(i,k))&
                        /dens(i,k)
                deltN(krow,i,k)   =  (0.55*ir_socol_loc(i,k))&
                         /dens(i,k)
             ENDIF
             deltOH(krow,i,k)  =  ((f_oh_gcr*ir_socol_loc(i,k))/dens(i,k))
          ENDIF !GCR True
          IF (spe) THEN
             deltOH(krow,i,k)  =  ((f_oh_spe*spen_loc(i,k,krow))/dens(i,k))
          ENDIF ! SPE True

          ! Update global tracer variables
          pxtte(i,k,idt_no) = deltNO(krow,i,k) + pxtte(i,k,idt_no)
          pxtte(i,k,idt_n) = deltN(krow,i,k) + pxtte(i,k,idt_n)
          pxtte(i,k,idt_oh) = deltOH(krow,i,k) + pxtte(i,k,idt_oh)


       ENDDO ! klev

       ! Add NO production by EEP only in the two uppermost layers
       IF (eep) THEN
          ap_NO(krow,i,1)=0.0_dp
          ap_NO(krow,i,2)=0.0_dp
          dz = (zlevb(i,1)-zlevb(i,2))*1e5
          ap_NO(krow,i,1)  =  (ap_flux_tot_loc(i,krow)*0.5/dens(i,1)/dz(1))
          ! Update global tracer variables
          pxtte(i,1,idt_no) = pxtte(i,1,idt_no) + ap_NO(krow,i,1)

          dz = (zlevb(i,2)-zlevb(i,3))*1e5
          ap_NO(krow,i,2)  =  (ap_flux_tot_loc(i,krow)*0.5/dens(i,2)/dz(1))
          ! Update global tracer variables
          pxtte(i,2,idt_no) = pxtte(i,2,idt_no) + ap_NO(krow,i,2)
       ENDIF ! eep

    ENDDO ! kproma

    END SUBROUTINE socol_gcr_spe_ionization
! ###########################################################################################
    SUBROUTINE geo2mag(mplat,mplon)

      ! Subroutine geo2mag
      ! convert from geographic to geomagnetic coordinates
      ! written by pascal saint-hilaire (Saint-Hilaire@astro.phys.ethz.ch) May 2002
      ! rewritten from IDL to fortran90 by Marco Calisto IAC-ETHZ July 2007
      ! called once a month from em_init_actm.F

     ! Using modules for grid/coord. conversion
      USE mo_constants,                 ONLY: api
      USE mo_control,                   ONLY: ngl, nlon

      IMPLICIT NONE
    

      ! Declaring variables needed for this subroutine

      REAL       :: glon(nlon),glat(ngl) ! geographical lon and lat (deg.)
      REAL       :: mplon,mplat          ! geomagnetic lon and lat (deg.) of the North magnetic Pole
      REAL       :: Dlong,Dlat, rlat, rlon
      REAL       :: galt, x, y, z
      REAL, PARAMETER :: d2r=api/180
      REAL  :: geolong2maglong(3,3), tomaglat(3,3) ! conversion vars
      REAL  :: coord(3), out(3)
      INTEGER   :: i, j
      REAL, DIMENSION(ngl)   :: convglat ! Gaussian latitudes, inequally spaced, 96*48 grid
      REAL, DIMENSION(nlon)    :: convglon

      IF (.NOT. ALLOCATED(mlat)) &
         ALLOCATE(mlat(nlon,ngl))   ! geomagnetic lat (deg.)

      glon = 0.
      glat = 0.

      IF (ngl==48) THEN
      	! Gaussian latitudes data in T31
      	convglat=(/87.1591, 83.4789, 79.7770, 76.0702, 72.3616,&
	&                   68.6520, 64.9419, 61.2316, 57.5210, 53.8103,&
	&                   50.0995, 46.3886, 42.6776, 38.9666,35.2556,31.5445,&
	&                   27.8334, 24.1223,20.4112, 16.7001,12.9890,  9.2779,&
	&                   5.5667,  1.8556,-1.8556,-5.5667, -9.2779,-12.9890,&
	&                   -16.7001,-20.4112,-24.1223,-27.8334,-31.5445,&
	&                   -35.2556,-38.9666,-42.6776,-46.3886,-50.0995,&
	&                   -53.8103,-57.5210,-61.2316,-64.9419,-68.6520,&
	&                  -72.3616,-76.0702,-79.7770,-83.4789,-87.1591/)

	   ! Make longitude mesh in T31 resolution

	  	DO i = 2,nlon
        	 glon(i) = glon(i-1) + 3.75
      	ENDDO

      ELSEIF (ngl==64) THEN
      ! Gaussian latitudes data in T42
      convglat=(/87.863, 85.096, 82.312, 79.525, 76.736, 73.947, 71.157, 68.367,&
	&                    65.577, 62.787, 59.997, 57.206, 54.416, 51.625, 48.835, 46.044,&
	&                    43.254, 40.463, 37.673, 34.882, 32.091, 29.301, 26.510, 23.720,&
	&                    20.929, 18.138, 15.348, 12.557,  9.767,  6.976,  4.185,  1.395,&
	&                    -1.395,  -4.185,  -6.976,  -9.767, -12.557, -15.348, -18.138, -20.929,&
	&                    -23.720, -26.510, -29.301, -32.091, -34.882, -37.673, -40.463, -43.254,&
	&                    -46.044, -48.835, -51.625, -54.416, -57.206, -59.997, -62.787, -65.577,&
	&                    -68.367, -71.157, -73.947, -76.736, -79.525, -82.312, -85.096, -87.863/)


    ! Make longitude mesh in T42 resolution, not needed anymore -> convglon
      glon(1) = 0.
	      DO i = 2,nlon
	         glon(i) = glon(i-1) + 2.8
	      ENDDO
     ENDIF


    glat(:) = convglat(:)

    !convert to radians

      Dlong = mplon*d2r
      Dlat  = mplat*d2r

      DO i = 1, nlon
         DO j = 1, ngl
            mlat(i,j) = 0.
            rlat = glat(j)*d2r
            rlon = glon(i)*d2r
            galt = 1.

            !convert to rectangular coordinates
            !X-axis: defined by the vector going form Earth's center towards
            !the intersection of the equator and greenwich's meridian.
            !Z-axis: axis of the geographic poles
            !Y-axis: defined by Y=Z**X

            coord(1) = galt*cos(rlat)*cos(rlon)
            coord(2) = galt*cos(rlat)*sin(rlon)
            coord(3) = galt*sin(rlat)

            !compute 1st rotation matrix: rotation around plane of the equatro,
            !from greenwich meridian to meridian containing magnetic dipole pole.

            geolong2maglong(1,1) =  cos(Dlong)
            geolong2maglong(1,2) =  sin(Dlong)
            geolong2maglong(2,1) = -sin(Dlong)
            geolong2maglong(2,2) =  cos(Dlong)
            geolong2maglong(3,3) = 1.

            out=matmul(geolong2maglong,coord)

            !second rotation: in the plane of the current meridian from geographic
            !pole to magnetic dipole pole.

            tomaglat(1,1) =  cos(api/2. -Dlat)
            tomaglat(1,3) = -sin(api/2. -Dlat)
            tomaglat(3,1) =  sin(api/2. -Dlat)
            tomaglat(3,3) =  cos(api/2. -Dlat)
            tomaglat(2,2) = 1.

            out=matmul(tomaglat,out)

            !convert back to latitude, longitude and altitude

            mlat(i,j)=atan2(out(3),sqrt(out(1)**2+(out(2)**2)))/d2r

        ENDDO
	 ENDDO

     RETURN

     END SUBROUTINE geo2mag

    SUBROUTINE line(N1,X1,Y1,N2,X2,Y2)


!     *****************************************************
!     *                                                   *
!     *      SUB. LINEAR INTERPOLATION                   *
!     *      Y1,X1,N1 - FUNCTION,ARGUMENTS AND NUMBER OF  *
!     *                 INPUT DATA                        *
!     *      X2,N2    - ARGUMENTS AND NUMBERS OF GRID     *
!     *      Y2       - OUTPUT                            *
!     *****************************************************

     IMPLICIT NONE

     INTEGER :: i,i1,i2,j2
     INTEGER, INTENT(in) :: N1,N2
     REAL(dp), INTENT(in) :: X1(N1), X2(N2), Y1(N1)
     REAL(dp), INTENT(out) :: Y2(N2)
     REAL :: S,X

    DO i=1,N2
       Y2(i) = 0.0_dp
    ENDDO

 IF(N1.LE.0.OR.N2.LE.0.OR.N1.GT.1000000.OR.N2.GT.1000000) GOTO 9
      IF( N1.EQ.1) GOTO 7
      S = SIGN(1.0,X1(N1) - X1(1))
      i2 = 1
      i1 = 2
    1 IF((X2(i2) - X1(1))*S) 2,2,3
    2 Y2(i2) = Y1(1)
      i2 = i2 + 1
      IF(i2.GT.N2) RETURN
      GOTO 1
    3 IF((X2(i2) - X1(i1))*S) 4,4,5
    4 X = (X2(i2) - X1(i1-1))/(X1(i1)-X1(i1-1))
      Y2 (i2) = Y1(i1) *X + Y1(i1 - 1) * (1.0 - X)
      i2 = i2 + 1
      IF(i2.GT.N2) RETURN
      GOTO 3
    5 i1 = i1 + 1
      IF(i1.LE.N1) GOTO 3
      DO 6 j2 = i2,N2
    6 Y2(j2) = Y1(N1)
      RETURN
    7 DO 8 i = 1,N2
    8 Y2(i) = Y1(i)
      RETURN
    9 WRITE (*,*)' ERROR LIN'
      STOP

  END SUBROUTINE line


END MODULE mo_socol_ionization

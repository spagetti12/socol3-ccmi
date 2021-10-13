MODULE mo_socol_strataerosols

  ! Description:

  ! Organises stratospheric aerosols (reading data set; interpolation to
  ! current time step; transformation of extinction coefficients to optical
  ! depth). Stratospheric aerosols are used for radiation (LW: optical depth, 
  ! see *rad_int*; SW: optical depth, single scattering factor, 
  ! asymmetry factor, see *swclr*) and for heterogeneous chemistry (surface
  ! area density, number density, see *em_chemini* and *em_hetero*).
  ! Data source: modified SAGE 2 retrieval (Larry Thomason) for SAD and
  ! effective radius. Extinction coefficients via Mie calculations (Beiping
  ! Luo).
  ! In contrast to ECHAM4, all fields are stored from top to bottom!
  !
  ! M. Schraner, ETH Zürich, November 2008

  USE mo_control,                 ONLY: nlon, ngl, nlev
  USE mo_decomposition,           ONLY: lc => local_decomposition, &
                                        global_decomposition
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_io
  USE mo_kind,                    ONLY: dp
  USE mo_mpi,                     ONLY: p_io, p_parallel, p_parallel_io, &
                                        p_bcast
  USE mo_parrrtm,                 ONLY: jpband       ! number of bands in rrtm
  USE mo_socol_grid_calculations, ONLY: zlevb, zaetr_socol
  USE mo_socol_interpo,           ONLY: wgt1_chem, wgt2_chem, wgt1_rad, &
                                        wgt2_rad, m3w1_chem, m3w2_chem, &
                                        m3w1_rad, m3w2_rad
  USE mo_socol_namelist
  USE mo_socol_readfile,          ONLY: socol_read_netcdf
  USE mo_socol_time_control,      ONLY: l_trigchem
  USE mo_sw,                      ONLY: nsw          ! number of SW bands
  USE mo_time_control,            ONLY: get_time_step, l_trigrad
  USE mo_transpose,               ONLY: scatter_gp    

  IMPLICIT NONE

  PRIVATE

  !--- Module variables:

  ! Vertical dimension of data set: first / last corresponding ECHAM-index:
  INTEGER :: k0_strataer, k1_strataer
  
  ! Preceeding, current and following month:
  REAL(dp), ALLOCATABLE :: ext_lw_m3(:,:,:,:,:)   ! LW ext coeff
  REAL(dp), ALLOCATABLE :: ext_sw_m3(:,:,:,:,:)   ! SW ext coeff
  REAL(dp), ALLOCATABLE :: omega_sw_m3(:,:,:,:,:) ! Single scattering albedo
  REAL(dp), ALLOCATABLE :: g_sw_m3(:,:,:,:,:)     ! Asymmetry factor
  REAL(dp), ALLOCATABLE :: sad_m3(:,:,:,:)      ! Surface area density [cm2/cm3]
!!$  REAL(dp), ALLOCATABLE :: nd_m3(:,:,:,:)       ! Number density [1/cm3]
  REAL(dp), ALLOCATABLE :: rm_m3(:,:,:,:)       ! Mean Radius [um]
  REAL(dp), ALLOCATABLE :: av_m3(:,:,:,:)       ! Volume density [um3/cm3]

  ! Interpolated to current time step:
  REAL(dp), ALLOCATABLE, PUBLIC :: tau_strataer_lw(:,:,:)   ! Optical depth
  REAL(dp), ALLOCATABLE, PUBLIC :: tau_strataer_sw(:,:,:)   ! Optical depth
  REAL(dp), ALLOCATABLE, PUBLIC :: omega_strataer_sw(:,:,:) ! Single scattering
                                                            ! albedo
  REAL(dp), ALLOCATABLE, PUBLIC :: g_strataer_sw(:,:,:)     ! Asymmetry factor
  REAL(dp), ALLOCATABLE, PUBLIC :: sad_strataer(:,:) ! Surface area density
!!$  REAL(dp), ALLOCATABLE, PUBLIC :: nd_strataer(:,:)  ! Number density
  REAL(dp), ALLOCATABLE, PUBLIC :: rm_strataer(:,:)  ! Mean Radius [um]
  REAL(dp), ALLOCATABLE, PUBLIC :: av_strataer(:,:)  ! Volume density [um3/cm3]

  PUBLIC :: read_socol_strataerosols, interpolate_socol_strataerosols, &
       cleanup_socol_strataerosols

  ! Intrinsic functions: 
  INTRINSIC MAX, MIN

CONTAINS

  SUBROUTINE read_socol_strataerosols(yr,mo)

    ! Reads surface area density, number density, extinction coefficients,
    ! single scattering albedo and asymmetry factor of stratospheric aersols 
    ! (data of current year, climatology of background aerosols or zero values,
    ! according to setting in namelist).

    ! *read_socol_strataerosols* is called from *read_socol_bcond_m*,
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER :: yr, mo

    ! Local variables:
    INTEGER :: nsw_file, nlw_file
    CHARACTER(31), PARAMETER :: fnvar = 'strataer', fnbg = 'strataerbg'
    CHARACTER(31) :: varname
    CHARACTER(50) :: varname_longname


    ! Executable statements:

    ! Initial values of k0_strataer, k1_strataer, nlw_file and nsw_file
    ! (only used if no data from files are read):
    k0_strataer = 1
    k1_strataer = nlev
    nsw_file = nsw
    nlw_file = jpband

    ! Extinction coefficients (LW):
    varname = 'EXT_LW'
    varname_longname = 'extinction coefficients (LW)'
    IF (lstrataer_var_rad) THEN
       CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
            varname_longname=varname_longname, extradimname='LW', &
            data5d=ext_lw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nlw_file)
    ELSE IF (lstrataer_bg_rad) THEN
       CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
            varname_longname=varname_longname, extradimname='LW', &
            data5d=ext_lw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nlw_file)
    ENDIF

    ! Extinction coefficients (SW):
    varname = 'EXT_SW'
    varname_longname = 'extinction coefficients (SW)'
    IF (lstrataer_var_rad) THEN
       CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
            varname_longname=varname_longname, extradimname='SW', &
            data5d=ext_sw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nsw_file)
    ELSE IF (lstrataer_bg_rad) THEN
       CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
            varname_longname=varname_longname, extradimname='SW', &
            data5d=ext_sw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nsw_file)
    ENDIF

    ! Single scattering albedo (SW):
    varname = 'OMEGA_SW'
    varname_longname = 'single scattering albedo (SW)'
    IF (lstrataer_var_rad) THEN
       CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
            varname_longname=varname_longname, extradimname='SW', &
            data5d=omega_sw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nsw_file)
    ELSE IF (lstrataer_bg_rad) THEN
       CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
            varname_longname=varname_longname, extradimname='SW', &
            data5d=omega_sw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nsw_file)
    ENDIF

    ! Asymmetry factor (SW):
    varname = 'G_SW'
    varname_longname = 'asymmetry factor (SW)'
    IF (lstrataer_var_rad) THEN
       CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
            varname_longname=varname_longname, extradimname='SW', &
            data5d=g_sw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nsw_file)
    ELSE IF (lstrataer_bg_rad) THEN
       CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
            varname_longname=varname_longname, extradimname='SW', &
            data5d=g_sw_m3, lowlevind=k0_strataer, uplevind=k1_strataer, &
            nextradim=nsw_file)
    ENDIF

    IF (lchem) THEN
       ! Surface area density:
       varname = 'SAD'
       varname_longname = 'surface area density'
       IF (lstrataer_var_chem) THEN
          CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
               varname_longname=varname_longname, data4d=sad_m3, &
               lowlevind=k0_strataer, uplevind=k1_strataer)
       ELSE IF (lstrataer_bg_chem) THEN
          CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
               varname_longname=varname_longname, data4d=sad_m3, &
               lowlevind=k0_strataer, uplevind=k1_strataer)
       ENDIF

!!$       ! Number density:
!!$       varname = 'ND'
!!$       varname_longname = 'number density'
!!$       IF (lstrataer_var_chem) THEN
!!$          CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
!!$               varname_longname=varname_longname, data4d=nd_m3, &
!!$               lowlevind=k0_strataer, uplevind=k1_strataer)
!!$       ELSE IF (lstrataer_bg_chem) THEN
!!$          CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
!!$               varname_longname=varname_longname, data4d=nd_m3, &
!!$               lowlevind=k0_strataer, uplevind=k1_strataer)
!!$       ENDIF
!!$    ENDIF

       ! Mean Radius:
       varname = 'R_MEAN'
       varname_longname = 'Mean radius'
       IF (lstrataer_var_chem) THEN
          CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
               varname_longname=varname_longname, data4d=rm_m3, &
               lowlevind=k0_strataer, uplevind=k1_strataer)
       ELSE IF (lstrataer_bg_chem) THEN
          CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
               varname_longname=varname_longname, data4d=rm_m3, &
               lowlevind=k0_strataer, uplevind=k1_strataer)
       ENDIF
       
       ! Volume density:
       varname = 'VOLUME_DENSITY'
       varname_longname = 'Volume density'
       IF (lstrataer_var_chem) THEN
          CALL socol_read_netcdf(fnvar, varname, 'LATLEV', yr=yr, mo=mo, &
               varname_longname=varname_longname, data4d=av_m3, &
               lowlevind=k0_strataer, uplevind=k1_strataer)
       ELSE IF (lstrataer_bg_chem) THEN
          CALL socol_read_netcdf(fnbg, varname, 'LATLEV', mo=mo, &
               varname_longname=varname_longname, data4d=av_m3, &
               lowlevind=k0_strataer, uplevind=k1_strataer)
       ENDIF
    ENDIF

    ! No stratospheric aerosols:
    IF (.NOT. lstrataer_var_rad .AND. .NOT. lstrataer_bg_rad) THEN
       IF (.NOT. ALLOCATED(ext_lw_m3)) &
            ALLOCATE(ext_lw_m3(lc%nproma,k0_strataer:k1_strataer, &
            jpband,lc%ngpblks,0:2))
       IF (.NOT. ALLOCATED(ext_sw_m3)) &
            ALLOCATE(ext_sw_m3(lc%nproma,k0_strataer:k1_strataer, &
            nsw,lc%ngpblks,0:2))
       IF (.NOT. ALLOCATED(omega_sw_m3)) &
            ALLOCATE(omega_sw_m3(lc%nproma,k0_strataer:k1_strataer, &
            nsw,lc%ngpblks,0:2))
       IF (.NOT. ALLOCATED(g_sw_m3)) &
            ALLOCATE(g_sw_m3(lc%nproma,k0_strataer:k1_strataer, &
            nsw,lc%ngpblks,0:2))

       WRITE (message_text,*) &
            'NO stratospheric aerosols for radiation'
       CALL message('',TRIM(message_text))
       ext_lw_m3(:,:,:,:,:)   = 0._dp
       ext_sw_m3(:,:,:,:,:)   = 0._dp
       omega_sw_m3(:,:,:,:,:) = 0._dp
       g_sw_m3(:,:,:,:,:)     = 0._dp
    ENDIF

    IF (lchem) THEN
       IF (.NOT. lstrataer_var_chem .AND. .NOT. lstrataer_bg_chem) THEN      
          IF (.NOT. ALLOCATED(sad_m3)) &
               ALLOCATE(sad_m3(lc%nproma, k0_strataer:k1_strataer, &
               lc%ngpblks,0:2))
!!$          IF (.NOT. ALLOCATED(nd_m3)) &
!!$               ALLOCATE(nd_m3(lc%nproma, k0_strataer:k1_strataer, &
!!$               lc%ngpblks,0:2))
          IF (.NOT. ALLOCATED(rm_m3)) &
               ALLOCATE(rm_m3(lc%nproma, k0_strataer:k1_strataer, &
               lc%ngpblks,0:2))
          IF (.NOT. ALLOCATED(av_m3)) &
               ALLOCATE(av_m3(lc%nproma, k0_strataer:k1_strataer, &
               lc%ngpblks,0:2))      

          WRITE (message_text,*) &   
               'NO stratospheric aerosols for chemistry' 
          CALL message('',TRIM(message_text))
          sad_m3(:,:,:,:) = 0._dp
!!$          nd_m3(:,:,:,:)  = 0._dp
          rm_m3(:,:,:,:)  = 0._dp
          av_m3(:,:,:,:)  = 0._dp
       ENDIF
    ENDIF

    ! Check number of LW/SW bands:
    IF (nlw_file .NE. jpband) THEN
       WRITE(message_text,*) 'Number of LW bands in data file: ', nlw_file
       CALL message('',TRIM(message_text))
       WRITE(message_text,*) 'Number of LW bands in ECHAM5: ', jpband
       CALL message('',TRIM(message_text))
       CALL finish('read_strataerosols','Run terminated')      
    ENDIF
    IF (nsw_file .NE. nsw) THEN
       WRITE(message_text,*) 'Number of SW bands in data file: ', nsw_file
       CALL message('',TRIM(message_text))
       WRITE(message_text,*) 'Number of SW bands in ECHAM5: ', nsw
       CALL message('',TRIM(message_text))
       CALL finish('read_strataerosols','Run terminated')      
    ENDIF

  END SUBROUTINE read_socol_strataerosols
 

  SUBROUTINE interpolate_socol_strataerosols(krow, kproma, kbdim, klev)

    ! Interpolates monthly stratospheric aerosol data to the current time step
    ! and multiplies data by zaetr_socol (to avoid stratospheric aerosols in
    ! the troposphere). Besides extinction coefficients are multiplied by
    ! thickness of model layers (-> transformation to optical depths).

    ! *interpolate_socol_strataerosols* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.

    ! Subroutine arguments:
    INTEGER, INTENT(in)   :: krow, kproma, kbdim, klev

    ! Local variables:
    INTEGER :: jl, jk
    REAL(dp) :: dz

    IF (l_trigrad) THEN
       IF (.NOT. ALLOCATED(tau_strataer_lw)) &
            ALLOCATE(tau_strataer_lw(kbdim,klev,jpband))
       IF (.NOT. ALLOCATED(tau_strataer_sw)) &
            ALLOCATE(tau_strataer_sw(kbdim,klev,nsw))
       IF (.NOT. ALLOCATED(omega_strataer_sw)) &
            ALLOCATE(omega_strataer_sw(kbdim,klev,nsw))
       IF (.NOT. ALLOCATED(g_strataer_sw)) &
            ALLOCATE(g_strataer_sw(kbdim,klev,nsw))

       ! No stratospheric aerosols outside of data set:
       tau_strataer_lw(1:kproma,1:k0_strataer-1,:) = 0._dp
       tau_strataer_lw(1:kproma,k1_strataer+1:klev,:) = 0._dp
       tau_strataer_sw(1:kproma,1:k0_strataer-1,:) = 0._dp
       tau_strataer_sw(1:kproma,k1_strataer+1:klev,:) = 0._dp
       omega_strataer_sw(1:kproma,1:k0_strataer-1,:) = 0.5_dp
       omega_strataer_sw(1:kproma,k1_strataer+1:klev,:) = 0.5_dp
       g_strataer_sw(1:kproma,1:k0_strataer-1,:) = 0._dp
       g_strataer_sw(1:kproma,k1_strataer+1:klev,:) = 0._dp
       
       tau_strataer_lw(1:kproma,k0_strataer:k1_strataer,:) = wgt1_rad * &
            ext_lw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w1_rad) + &
            wgt2_rad * &
            ext_lw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w2_rad)
       tau_strataer_sw(1:kproma,k0_strataer:k1_strataer,:) = wgt1_rad * &
            ext_sw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w1_rad) + &
            wgt2_rad * &
            ext_sw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w2_rad)
       omega_strataer_sw(1:kproma,k0_strataer:k1_strataer,:) = wgt1_rad * &
            omega_sw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w1_rad) + &
            wgt2_rad * &
            omega_sw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w2_rad)
       g_strataer_sw(1:kproma,k0_strataer:k1_strataer,:) = wgt1_rad * &
            g_sw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w1_rad) + &
            wgt2_rad * &
            g_sw_m3(1:kproma,k0_strataer:k1_strataer,:,krow,m3w2_rad)
       
       ! Multiplication by thickness of the layer [km]
       ! (extinction -> optical depth) and by zaetr_socol:
       DO jl = 1, kproma
          DO jk = k0_strataer, k1_strataer 
             dz = zlevb(jl,jk)-zlevb(jl,jk+1)
             tau_strataer_lw(jl,jk,:) = &
                  dz*zaetr_socol(jl,jk)*tau_strataer_lw(jl,jk,:)
             tau_strataer_sw(jl,jk,:) = &
                  dz*zaetr_socol(jl,jk)*tau_strataer_sw(jl,jk,:)
          ENDDO
       ENDDO
       
       ! Optical depths are at least epsilon:
       tau_strataer_lw(1:kproma,:,:) = &
            MAX(tau_strataer_lw(1:kproma,:,:),EPSILON(1._dp))
       tau_strataer_sw(1:kproma,:,:) = &
            MAX(tau_strataer_sw(1:kproma,:,:),EPSILON(1._dp))
       
       ! Single scattering albedo is at least epsilon and at most 1-epsilon:
       omega_strataer_sw(1:kproma,k0_strataer:k1_strataer,:) = &
            MAX(omega_strataer_sw(1:kproma,k0_strataer:k1_strataer,:), &
            EPSILON(1._dp))
       omega_strataer_sw(1:kproma,k0_strataer:k1_strataer,:) = &
            MIN(omega_strataer_sw(1:kproma,k0_strataer:k1_strataer,:), &
            1._dp-EPSILON(1._dp))
    ENDIF

    IF (l_trigchem) THEN
       IF (.NOT. ALLOCATED(sad_strataer)) ALLOCATE(sad_strataer(kbdim,klev))
!!$       IF (.NOT. ALLOCATED(nd_strataer)) ALLOCATE(nd_strataer(kbdim,klev))
       IF (.NOT. ALLOCATED(rm_strataer)) ALLOCATE(rm_strataer(kbdim,klev))
       IF (.NOT. ALLOCATED(av_strataer)) ALLOCATE(av_strataer(kbdim,klev))

       ! No stratospheric aerosols outside of data set:
       sad_strataer(1:kproma,1:k0_strataer) = 0._dp
       sad_strataer(1:kproma,k1_strataer+1:klev) = 0._dp
!!$       nd_strataer(1:kproma,1:k0_strataer-1) = 0._dp
!!$       nd_strataer(1:kproma,k1_strataer+1:klev) = 0._dp
       rm_strataer(1:kproma,1:k0_strataer-1) = 0._dp
       rm_strataer(1:kproma,k1_strataer+1:klev) = 0._dp
       av_strataer(1:kproma,1:k0_strataer-1) = 0._dp
       av_strataer(1:kproma,k1_strataer+1:klev) = 0._dp

       sad_strataer(1:kproma,k0_strataer:k1_strataer) = wgt1_chem * &
            sad_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w1_chem) + &
            wgt2_chem * sad_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w2_chem)
!!$       nd_strataer(1:kproma,k0_strataer:k1_strataer) = wgt1_chem * &
!!$            nd_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w1_chem) + &
!!$            wgt2_chem * nd_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w2_chem)
       rm_strataer(1:kproma,k0_strataer:k1_strataer) = wgt1_chem * &
            rm_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w1_chem) + &
            wgt2_chem * rm_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w2_chem)
       av_strataer(1:kproma,k0_strataer:k1_strataer) = wgt1_chem * &
            av_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w1_chem) + &
            wgt2_chem * av_m3(1:kproma,k0_strataer:k1_strataer,krow,m3w2_chem)       

       ! Multiplication with zaetr_socol:
       DO jl = 1, kproma
          DO jk = k0_strataer, k1_strataer
             sad_strataer(jl,jk) = zaetr_socol(jl,jk)*sad_strataer(jl,jk)
          ENDDO
       ENDDO

       ! Surface area density is at least epsilon:
       sad_strataer(1:kproma,:) = MAX(sad_strataer(1:kproma,:),EPSILON(1._dp))
    ENDIF

  END SUBROUTINE interpolate_socol_strataerosols


  SUBROUTINE cleanup_socol_strataerosols

    ! Deallocates module variables.
    
    ! *cleanup_socol_strataerosols* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(ext_lw_m3)) DEALLOCATE(ext_lw_m3)
    IF (ALLOCATED(ext_sw_m3)) DEALLOCATE(ext_sw_m3)
    IF (ALLOCATED(omega_sw_m3)) DEALLOCATE(omega_sw_m3)
    IF (ALLOCATED(g_sw_m3)) DEALLOCATE(g_sw_m3)
    IF (ALLOCATED(sad_m3)) DEALLOCATE(sad_m3)
!!$    IF (ALLOCATED(nd_m3)) DEALLOCATE(nd_m3)
    IF (ALLOCATED(rm_m3)) DEALLOCATE(rm_m3)
    IF (ALLOCATED(av_m3)) DEALLOCATE(av_m3)
    IF (ALLOCATED(tau_strataer_lw)) DEALLOCATE(tau_strataer_lw)
    IF (ALLOCATED(tau_strataer_sw)) DEALLOCATE(tau_strataer_sw)
    IF (ALLOCATED(omega_strataer_sw)) DEALLOCATE(omega_strataer_sw)
    IF (ALLOCATED(g_strataer_sw)) DEALLOCATE(g_strataer_sw)
    IF (ALLOCATED(sad_strataer)) DEALLOCATE(sad_strataer)
!!$    IF (ALLOCATED(nd_strataer)) DEALLOCATE(nd_strataer)
    IF (ALLOCATED(rm_strataer)) DEALLOCATE(rm_strataer)
    IF (ALLOCATED(av_strataer)) DEALLOCATE(av_strataer)

  END SUBROUTINE cleanup_socol_strataerosols
      
END MODULE mo_socol_strataerosols

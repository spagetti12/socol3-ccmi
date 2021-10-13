MODULE mo_socol_photo

  ! Description:
  !
  ! Computes effect of cloudiness on photolysis rates
  ! Parameterization: Chang et al., JGR, 92 (D12), 1987.
  !
  ! Andrea Stenke, ETH Zurich, May 2010
  !
  ! Andrea Stenke, ETH Zurich, August 2010
  ! 
  ! Reads in pre-calculated photolysis rates for PAN, CH3CO3H, CH3COCHO
  ! photolysis rates based on MESSy-simulation

  USE mo_kind,           ONLY: dp
  USE mo_socol_interpo,  ONLY: wgt1_chem, wgt2_chem, m3w1_chem, m3w2_chem
  USE mo_linked_list,    ONLY: HYBRID, NETCDF
  USE mo_memory_base,    ONLY: new_stream, add_stream_element,  &
                               default_stream_setting, add_stream_reference, &
                               delete_stream, t_stream
  USE mo_time_control,   ONLY: delta_time, lstart
  USE mo_socol_readfile, ONLY: socol_read_netcdf

  IMPLICIT NONE

  PRIVATE

  TYPE (t_stream), PUBLIC, POINTER :: photo     !the photo stream

  REAL(dp), PUBLIC, POINTER :: cloud_mul(:,:,:)

  REAL(dp), PUBLIC, POINTER :: cloud_mul_d(:,:,:)

  ! eth_as_tropchem+
  REAL(dp), ALLOCATABLE :: pan_jval_m3(:,:,:,:,:),  & 
                           ch3co3h_jval_m3(:,:,:,:,:),  &
                           mgly_jval_m3(:,:,:,:,:)

  REAL(dp), ALLOCATABLE, PUBLIC :: pan_jval(:,:),  & 
                                   ch3co3h_jval(:,:),  &
                                   mgly_jval(:,:)
  ! eth_as_tropchem-

  PUBLIC :: cloud_mod
  PUBLIC :: construct_stream_photo   ! construct stream
  PUBLIC :: destruct_stream_photo    ! destruct stream
  PUBLIC :: init_stream_photo        ! initialize stream
  PUBLIC :: accumulate_stream_photo  ! accumulate stream elements

  ! eth_as_tropchem+
  PUBLIC :: read_socol_photo
  PUBLIC :: interpolate_socol_photo  
  PUBLIC :: cleanup_socol_photo  
  ! eth_as_tropchem-

CONTAINS

      subroutine cloud_mod( zen_angle, clouds, lwc, delp, srf_alb, &
                            eff_alb, cld_mult )
!-----------------------------------------------------------------------
! 	... cloud alteration factors for photorates and albedo
!-----------------------------------------------------------------------

      use mo_control, only : nlev

      implicit none

      real(dp), parameter :: gi = 1._dp/9.80616_dp

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      real(dp), intent(in)    ::  zen_angle         ! zenith angle
      real(dp), intent(in)    ::  srf_alb           ! surface albedo
      real(dp), intent(in)    ::  clouds(nlev)       ! cloud fraction
      real(dp), intent(in)    ::  lwc(nlev)          ! liquid water content (mass mr)
      real(dp), intent(in)    ::  delp(nlev)         ! del press about midpoint in pascals
      real(dp), intent(out)   ::  eff_alb(nlev)      ! effective albedo
      real(dp), intent(out)   ::  cld_mult(nlev)     ! photolysis mult factor
      real(dp)                ::  cld_mult1(nlev)     ! photolysis mult factor

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer :: k, nlevm
      real(dp)    :: coschi
      real(dp)    :: del_lwp(nlev)
      real(dp)    :: del_tau(nlev)
      real(dp)    :: above_tau(nlev)
      real(dp)    :: below_tau(nlev)
      real(dp)    :: above_cld(nlev)
      real(dp)    :: below_cld(nlev)
      real(dp)    :: above_tra(nlev)
      real(dp)    :: below_tra(nlev)
      real(dp)    :: fac1(nlev)
      real(dp)    :: fac2(nlev)
      real(dp)    :: fac3(nlev) ! eth_as

      nlevm = nlev - 1

!---------------------------------------------------------
!	... modify lwc for cloud fraction and form
!	    liquid water path for each layer
!---------------------------------------------------------
      where( clouds(:) /= 0._dp )
         del_lwp(:) = gi * lwc(:) * delp(:) * 1.e3_dp / clouds(:)
      elsewhere
	 del_lwp(:) = 0._dp
      endwhere
!---------------------------------------------------------
!    	... form tau for each model layer
!---------------------------------------------------------
      where( clouds(:) /= 0._dp )
	 del_tau(:) = del_lwp(:) *.155_dp * clouds(:)**1.5_dp
      elsewhere
	 del_tau(:) = 0._dp
      end where
!---------------------------------------------------------
!    	... form integrated tau from top down
!---------------------------------------------------------
      above_tau(1) = 0._dp
      do k = 1,nlevm
	 above_tau(k+1) = del_tau(k) + above_tau(k)
      end do
!---------------------------------------------------------
!    	... form integrated tau from bottom up
!---------------------------------------------------------
      below_tau(nlev) = 0._dp
      do k = nlevm,1,-1
	 below_tau(k) = del_tau(k+1) + below_tau(k+1)
      end do
!---------------------------------------------------------
!	... form vertically averaged cloud cover above and below
!---------------------------------------------------------
      above_cld(1) = 0._dp
      do k = 1,nlevm
	 above_cld(k+1) = clouds(k) * del_tau(k) + above_cld(k)
      end do
      do k = 2,nlev
	 if( above_tau(k) /= 0._dp ) then
	    above_cld(k) = above_cld(k) / above_tau(k)
	 else
	    above_cld(k) = above_cld(k-1)
	 end if
      end do
      below_cld(nlev) = 0._dp
      do k = nlevm,1,-1
	 below_cld(k) = clouds(k+1) * del_tau(k+1) + below_cld(k+1)
      end do
      do k = nlevm,1,-1
	 if( below_tau(k) /= 0._dp ) then
	    below_cld(k) = below_cld(k) / below_tau(k)
	 else
	    below_cld(k) = below_cld(k+1)
	 end if
      end do
!---------------------------------------------------------
!	... modify above_tau and below_tau via jfm
!---------------------------------------------------------
      where( above_cld(2:nlev) /= 0._dp )
	 above_tau(2:nlev) = above_tau(2:nlev) / above_cld(2:nlev)
      end where
      where( below_cld(:nlevm) /= 0._dp )
         below_tau(:nlevm) = below_tau(:nlevm) / below_cld(:nlevm)
      end where
      where( above_tau(2:nlev) < 5._dp )
	    above_cld(2:nlev) = 0._dp
      end where
      where( below_tau(:nlevm) < 5._dp )
	 below_cld(:nlevm) = 0._dp
      end where
!---------------------------------------------------------
!	... form transmission factors
!---------------------------------------------------------
      above_tra(:) = 11.905_dp / (9.524_dp + above_tau(:))
      below_tra(:) = 11.905_dp / (9.524_dp + below_tau(:))
!---------------------------------------------------------
!	... form effective albedo
!---------------------------------------------------------
      where( below_cld(:) /= 0._dp )
	 eff_alb(:) = srf_alb + below_cld(:) * (1._dp - below_tra(:)) &
                                             * (1._dp - srf_alb)
      elsewhere
	 eff_alb(:) = srf_alb
      end where
!++jsr
!      coschi = max( cos( zen_angle ),.5 )
      coschi = max ( zen_angle, .5_dp )
!--jsr
      where( del_lwp(:)*.155_dp < 5._dp )
	 fac1(:) = 0._dp
      elsewhere
	 fac1(:) = 1.4_dp * coschi - 1._dp
      end where
      fac2(:)     = min( 0._dp, 1.6_dp*coschi*above_tra(:) - 1._dp )
      fac3(:)     = min( 0._dp, 1._dp*(1._dp - below_tra(:))*coschi )

      cld_mult(:) = 1._dp + fac1(:) * clouds(:) + fac2(:) * above_cld(:) + fac3(:) * below_cld(:)
      cld_mult(:) = max( .05_dp,cld_mult(:) )

      end subroutine cloud_mod


  !-----------------------------------------------------------------------------------

  SUBROUTINE construct_stream_photo

    ! Allocates output streams

    ! *construct_stream_* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90.

  !--- Create new stream:

  CALL new_stream (photo ,'photo', filetype=NETCDF)      
  
  ! add entries for geopotential and log surface pressure (used by the
  ! afterburner):
  CALL add_stream_reference (photo, 'geosp'   ,'g3b'   ,lpost=.FALSE.)
  CALL add_stream_reference (photo, 'lsp'     ,'sp'    ,lpost=.FALSE.)

  CALL default_stream_setting (photo, lrerun    = .TRUE. , &
                                       leveltype = HYBRID , &
                                       table     = 199,     &
                                       laccu     = .TRUE.)

  CALL add_stream_element (photo, 'Cloud_mul', cloud_mul, lpost=.FALSE., &
                           longname='Photolysis cloud modification factor', &
                           units='-', code=100)

  CALL add_stream_element (photo, 'Cloud_mul_d', cloud_mul_d, lpost=.FALSE., &
                           laccu = .FALSE., longname='Photolysis cloud modification factor', &
                           units='-', code=106)

END SUBROUTINE construct_stream_photo

!---------------------------------------------------------------------------------

SUBROUTINE destruct_stream_photo

  ! Deallocates memory.

  ! *destruct_stream_* is called from *call_free_submodel_memory*,
  ! src/call_submodels.f90.

  CALL delete_stream(photo)

END SUBROUTINE destruct_stream_photo

!---------------------------------------------------------------------------------

SUBROUTINE init_stream_photo

  ! Initializes streams with zero.

  ! *init_stream_* is called from *call_init_tracers*,
  ! src/call_submodels.f90.

  IF (lstart) THEN     ! use restart fields otherwise

     cloud_mul(:,:,:)   = 0._dp
     cloud_mul_d(:,:,:) = 0._dp

  ENDIF
 
END SUBROUTINE init_stream_photo

!---------------------------------------------------------------------------------

SUBROUTINE accumulate_stream_photo (kproma, krow)

  USE mo_socol_namelist, ONLY: lch4_wetland

  ! This subroutine accumulates the current value of a variable at every time 
  ! step, such that finally monthly streams can be calculated.

  ! *accumulate_stream_* is called from *call_diagn*, 
  ! src/call_submodels

  ! Scalar arguments
  INTEGER, INTENT(in) :: kproma, krow


  ! Executable statements:

  cloud_mul(1:kproma,:,krow)  = cloud_mul(1:kproma,:,krow)  + delta_time * &
       cloud_mul_d(1:kproma,:,krow)

END SUBROUTINE accumulate_stream_photo

!---------------------------------------------------------------------------------

  SUBROUTINE read_socol_photo(yr,mo)

    ! Reads precalculated photolysis rates for current, preceeding and following month.

    ! *read_socol_photo* is called from *read_socol_bcond_m*,
    ! src/socol_boundary_conditions.90.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: yr, mo
    INTEGER :: nhr

    ! Local variables:
    CHARACTER(25), PARAMETER :: fn_photo = 'photolysis_rates_messy'

    CHARACTER(31) :: varname
    CHARACTER(50) :: varname_longname

    nhr = 5

    ! Executable statements:

    ! Photolysis rates PAN
    varname = 'J_PAN'
    varname_longname = 'J(PAN)'
    CALL socol_read_netcdf( fn_photo, varname, 'LONLATLEV', mo=mo, &
         varname_longname=varname_longname, extradimname='hour', data5d=pan_jval_m3, nextradim=nhr)

    ! Photolysis rates CH3CO3H = PAA
    varname = 'J_PAA'
    varname_longname = 'J(PAA)'
    CALL socol_read_netcdf( fn_photo, varname, 'LONLATLEV', mo=mo, &
         varname_longname=varname_longname, extradimname='hour', data5d=ch3co3h_jval_m3, nextradim=nhr)

    ! Photolysis rates CH3COCHO = MGLY
    varname = 'J_CH3COCHO'
    varname_longname = 'J(CH3COCHO)'
    CALL socol_read_netcdf( fn_photo, varname, 'LONLATLEV', mo=mo, &
         varname_longname=varname_longname, extradimname='hour', data5d=mgly_jval_m3, nextradim=nhr)

  END SUBROUTINE read_socol_photo

!---------------------------------------------------------------------------------

  SUBROUTINE interpolate_socol_photo(krow, kproma, kbdim, klev)

    ! Interpolates monthly photolysis rates to the current time step.

    ! *interpolate_socol_nmvoc* is called from 
    ! *interpolate_socol_bcond_local*, src/socol_boundary_conditions.

    USE mo_time_control, ONLY : get_date_components, current_date

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: krow, kproma, kbdim, klev
    INTEGER :: yr, mo, dy, hr, mn, se
    INTEGER :: ih1, ih2
  
    REAL(dp) :: wgt1, wgt2
    REAL(dp), DIMENSION(kbdim,klev,0:2) :: pan_jval_d, ch3co3h_jval_d, mgly_jval_d

    pan_jval_d = 0._dp
    ch3co3h_jval_d = 0._dp
    mgly_jval_d = 0._dp

    IF (.NOT. ALLOCATED(pan_jval)) &
         ALLOCATE(pan_jval(kbdim,klev))

    IF (.NOT. ALLOCATED(ch3co3h_jval)) &
         ALLOCATE(ch3co3h_jval(kbdim,klev))

    IF (.NOT. ALLOCATED(mgly_jval)) &
         ALLOCATE(mgly_jval(kbdim,klev))

    CALL get_date_components(current_date, yr, mo, dy, hr, mn, se)

    IF (hr .le. 5) THEN
       ih1 = 1
       ih2 = 2
       wgt1 = (5._dp - hr)/5._dp
       wgt2 = 1._dp - wgt1
    ELSEIF(hr .gt. 5 .and. hr .le. 10) THEN
       ih1 = 2
       ih2 = 3
       wgt1 = (10._dp - hr)/5._dp
       wgt2 = 1._dp - wgt1
    ELSEIF(hr .gt. 10 .and. hr .le. 15) THEN
       ih1 = 3
       ih2 = 4
       wgt1 = (15._dp - hr)/5._dp
       wgt2 = 1._dp - wgt1
    ELSEIF(hr .gt. 15 .and. hr .le. 20) THEN
       ih1 = 4
       ih2 = 5
       wgt1 = (20._dp - hr)/5._dp
       wgt2 = 1._dp - wgt1
    ELSEIF(hr .gt. 20 .and. hr .le. 24) THEN
       ih1 = 5
       ih2 = 1
       wgt1 = (24._dp - hr)/4._dp
       wgt2 = 1._dp - wgt1
    END IF

    ! interpolate daily cycle

    pan_jval_d(1:kproma,:,m3w1_chem) = &
         wgt1*pan_jval_m3(1:kproma,:,ih1,krow,m3w1_chem) + wgt2*pan_jval_m3(1:kproma,:,ih1,krow,m3w1_chem)
    pan_jval_d(1:kproma,:,m3w2_chem) = &
         wgt1*pan_jval_m3(1:kproma,:,ih1,krow,m3w2_chem) + wgt2*pan_jval_m3(1:kproma,:,ih1,krow,m3w2_chem)

    ch3co3h_jval_d(1:kproma,:,m3w1_chem) = &
         wgt1*ch3co3h_jval_m3(1:kproma,:,ih1,krow,m3w1_chem) + wgt2*ch3co3h_jval_m3(1:kproma,:,ih1,krow,m3w1_chem)
    ch3co3h_jval_d(1:kproma,:,m3w2_chem) = &
         wgt1*ch3co3h_jval_m3(1:kproma,:,ih1,krow,m3w2_chem) + wgt2*ch3co3h_jval_m3(1:kproma,:,ih1,krow,m3w2_chem)
 
    mgly_jval_d(1:kproma,:,m3w1_chem) = &
         wgt1*mgly_jval_m3(1:kproma,:,ih1,krow,m3w1_chem) + wgt2*mgly_jval_m3(1:kproma,:,ih1,krow,m3w1_chem)
    mgly_jval_d(1:kproma,:,m3w2_chem) = &
         wgt1*mgly_jval_m3(1:kproma,:,ih1,krow,m3w2_chem) + wgt2*mgly_jval_m3(1:kproma,:,ih1,krow,m3w2_chem)

    pan_jval(1:kproma,:) = &
         wgt1_chem*pan_jval_d(1:kproma,:,m3w1_chem) + &
         wgt2_chem*pan_jval_d(1:kproma,:,m3w2_chem)

    ch3co3h_jval(1:kproma,:) = &
         wgt1_chem*ch3co3h_jval_d(1:kproma,:,m3w1_chem) + &
         wgt2_chem*ch3co3h_jval_d(1:kproma,:,m3w2_chem)

    mgly_jval(1:kproma,:) = &
         wgt1_chem*mgly_jval_d(1:kproma,:,m3w1_chem) + &
         wgt2_chem*mgly_jval_d(1:kproma,:,m3w2_chem)


  END SUBROUTINE interpolate_socol_photo

!---------------------------------------------------------------------------------
  SUBROUTINE cleanup_socol_photo

    ! Deallocates module variables.
    
    ! *cleanup_socol_photo* is called from *cleanup_socol_bcond*,
    ! src/socol_boundary_conditions.f90.

    IF (ALLOCATED(pan_jval_m3))     DEALLOCATE(pan_jval_m3)
    IF (ALLOCATED(ch3co3h_jval_m3)) DEALLOCATE(ch3co3h_jval_m3)
    IF (ALLOCATED(mgly_jval_m3))    DEALLOCATE(mgly_jval_m3)

    IF (ALLOCATED(pan_jval))     DEALLOCATE(pan_jval)
    IF (ALLOCATED(ch3co3h_jval)) DEALLOCATE(ch3co3h_jval)
    IF (ALLOCATED(mgly_jval))    DEALLOCATE(mgly_jval)

  END SUBROUTINE cleanup_socol_photo

END MODULE mo_socol_photo

  

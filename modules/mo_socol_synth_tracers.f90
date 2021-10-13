MODULE mo_socol_synth_tracers
!
!       CCMI synthetic tracers
!       NH_5           CCMI tracer 30°-50°N, 5-day decay 
!       NH_50         CCMI tracer 30°-50°N, 50-day decay 
!       NH_50W     CCMI tracer 30°-50°N, 50-day decay + wet deposition 
!       ST80_25      CCMI tracer above 80hPa, 25-day troposphere decay 
!       CO_25         CCMI tracer anthropogenic CO, 25-day decay 
!       CO_50         CCMI tracer anthropogenic CO, 50-day decay 
!
!       Authors: 
!                Sukhodolov Timofei PMOD/WRC / IAC ETH, 2013
!                Stenke Andrea IAC ETH, 2013

    USE mo_kind,                      ONLY: dp
    USE mo_exception,                 ONLY: finish, message, message_text

IMPLICIT NONE

PRIVATE

   REAL(dp), ALLOCATABLE :: CO_emiss(:,:)          ! CO emissions from CCMI file

   INTEGER, PUBLIC :: idt_NH_5, idt_NH_50, idt_NH_50W, idt_ST80_25, idt_CO_25, idt_CO_50


   PUBLIC    :: request_synth_trac, init_synth_trac, synth_trac_flux, cleanup_synth_trac

 
CONTAINS

! ###########################################################################################

SUBROUTINE request_synth_trac

    ! *request_synth_trac* is called from *call_request_tracer*,
    ! src/call_submodels.f90.

    USE mo_advection,         ONLY: iadvec   ! selected advection scheme
    USE mo_tracer,            ONLY: new_tracer,CONSTANT,INITIAL,ON,OFF, &
                                     RESTART,GAS

   INTEGER                   :: ierr

    CALL new_tracer('NH_5','MEZON', idx=idt_NH_5,     nwrite=ON, ierr=ierr,  &
         longname='CCMI tracer 30-50°, 5-day decay',     units='mol/mol', &
         code=212, table=131, ninit=RESTART+CONSTANT+INITIAL, nrerun=ON, &
         ntran=iadvec, nvdiff=ON, nconv=ON, nphase=GAS, tdecay=4.32e+05_dp)

    CALL new_tracer('NH_50','MEZON', idx=idt_NH_50,     nwrite=ON, ierr=ierr,  &
         longname='CCMI tracer 30-50°, 50-day decay',     units='mol/mol', &
         code=213, table=131, ninit=RESTART+CONSTANT+INITIAL, nrerun=ON, &
         ntran=iadvec, nvdiff=ON, nconv=ON, nphase=GAS, tdecay=4.32e+06_dp)

    CALL new_tracer('NH_50W','MEZON', idx=idt_NH_50W,     nwrite=ON, ierr=ierr,  &
         longname='CCMI tracer 30°-50°N, 50-day decay + wet deposition',     units='mol/mol', &
         code=214, table=131, ninit=RESTART+CONSTANT+INITIAL, nrerun=ON, &
         ntran=iadvec, nvdiff=ON, nconv=ON, nphase=GAS, tdecay=4.32e+06_dp)

    CALL new_tracer('ST80_25','MEZON', idx=idt_ST80_25,     nwrite=ON, ierr=ierr,  &
         longname='CCMI tracer above 80hPa, 25-day troposphere decay', units='mol/mol', &
         code=215, table=131, ninit=RESTART+CONSTANT+INITIAL, nrerun=ON, &
         ntran=iadvec, nvdiff=ON, nconv=ON, nphase=GAS)

    CALL new_tracer('CO_25','MEZON', idx=idt_CO_25,     nwrite=ON, ierr=ierr,  &
         longname='CCMI tracer anthropogenic CO, 25-day decay', units='mol/mol', &
         code=216, table=131, ninit=RESTART+CONSTANT+INITIAL, nrerun=ON, &
         ntran=iadvec, nvdiff=ON, nconv=ON, nphase=GAS, tdecay=2.16e+06_dp)

    CALL new_tracer('CO_50','MEZON', idx=idt_CO_50,     nwrite=ON, ierr=ierr,  &
         longname='CCMI tracer anthropogenic CO, 50-day decay', units='mol/mol', &
         code=217, table=131, ninit=RESTART+CONSTANT+INITIAL, nrerun=ON, &
         ntran=iadvec, nvdiff=ON, nconv=ON, nphase=GAS, tdecay=4.32e+06_dp)

      IF (ierr /= 0 ) THEN
          CALL finish ('request_tracers_synthetic',' failed to define synthetic tracers')
      END IF

END SUBROUTINE request_synth_trac

! ###########################################################################################

SUBROUTINE init_synth_trac

    ! *init_socol_tracers* is called from *call_init_tracers*,
    ! src/call_submodels.f90.

  USE mo_decomposition,      ONLY: local_decomposition
  USE mo_socol_readfile,     ONLY: socol_read_netcdf
  USE mo_mpi,                       ONLY: p_pe, p_io

  IMPLICIT NONE
  CHARACTER(25), PARAMETER :: CO_syn_emiss = 'CO_syn_emiss'
  CHARACTER(31) :: varname

    IF (p_pe == p_io) write(*,*)'synth tracers initialization'

    varname = 'emission_flux'
    CALL socol_read_netcdf(CO_syn_emiss, varname, 'lonlat', data2d=CO_emiss)

END SUBROUTINE init_synth_trac

! ###########################################################################################

SUBROUTINE synth_trac_flux(krow, kproma, kbdim, klev, ktrac, pxtm1, pxtte, papp1, ktrpwmo, airmolec_col)

    ! Introduces CCMI synthetic tracers emissions and VMR.

    ! *synth_trac_flux* is called from *set_socol_bcond_fluxes*,
    ! src/socol_boundary_conditions.f90.

    USE mo_geoloc,                  ONLY: gboxarea_2d, philat_2d       !grid box area, latitude
    USE mo_socol_grid_calculations, ONLY: calculate_zaetr_socol,zaetr_socol
    USE mo_time_control,            ONLY: time_step_len
    USE mo_constants,               ONLY: g
    USE mo_socol_constants,         ONLY: atomweight, amco

!!$    USE mo_socol_grid_calculations, ONLY: zaetr_socol

    ! Subroutine arguments:
    INTEGER, INTENT(in)     :: krow, kproma, kbdim, klev, ktrac, ktrpwmo(kbdim)
    REAL(dp) ::    papp1(kbdim,klev)
    REAL(dp), INTENT(inout)  :: pxtm1(kbdim,klev,ktrac)
    REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)
    REAL(dp) :: airmolec_col(kbdim)

    ! Local variables:
    REAL(dp):: zxtp1(kproma,klev) ! tracer concentration at t+dt
    INTEGER  :: jl, jk
    INTEGER, PARAMETER  :: ST80_tdecay = 2.16e+06_dp

    INTEGER :: ktrpwmo1(kbdim)    

    DO jl=1,kproma
       ktrpwmo1(jl) = ktrpwmo(jl) + 1
    ENDDO

   CALL calculate_zaetr_socol(kproma, kbdim, klev, ktrpwmo1) ! calculation of tropopause level

   DO jl = 1, kproma

     if ((philat_2d(jl,krow).ge.30).and.(philat_2d(jl,krow).le.50)) then
        pxtm1(jl,klev,idt_NH_5)   = 100.e-09_dp
        pxtm1(jl,klev,idt_NH_50)  = 100.e-09_dp
        pxtm1(jl,klev,idt_NH_50W) = 100.e-09_dp

        pxtte(jl,klev,idt_NH_5)    = 0._dp   
        pxtte(jl,klev,idt_NH_50)   = 0._dp
        pxtte(jl,klev,idt_NH_50W)  = 0._dp

     endif

     DO jk = 1, klev
        if (papp1(jl,jk).lt.8000._dp) then
           pxtm1(jl,jk,idt_ST80_25) = 200.e-09_dp
           pxtte(jl,jk,idt_ST80_25) = 0._dp
        end if

       IF (zaetr_socol(jl,jk) .gt. 1._dp) THEN    ! Decay in troposphere only. Done like in *radionucl_sink* mo_sub_echam.f90
         zxtp1(jl,jk) = pxtm1(jl,jk,idt_ST80_25) + pxtte(jl,jk,idt_ST80_25) * time_step_len
         zxtp1(jl,jk) = EXP (-time_step_len/ST80_tdecay) * zxtp1(jl,jk)
         pxtte(jl,jk,idt_ST80_25) = (zxtp1(jl,jk)-pxtm1(jl,jk,idt_ST80_25))* 1.0_dp/time_step_len
       endif
     ENDDO

     pxtte(jl,klev,idt_CO_25) = pxtte(jl,klev,idt_CO_25)+CO_emiss(jl,krow)/airmolec_col(jl)/(atomweight*amco)*time_step_len  !CO units are recalculated from kg/m2/s to mol/mol like in *set_socol_bcond_fluxes* socol_boundary_conditions.f90.
     pxtte(jl,klev,idt_CO_50) = pxtte(jl,klev,idt_CO_50)+CO_emiss(jl,krow)/airmolec_col(jl)/(atomweight*amco)*time_step_len

   END DO

END SUBROUTINE synth_trac_flux

! ###########################################################################################

SUBROUTINE cleanup_synth_trac

    ! Deallocates module variables.

    ! *cleanup_synth_trac* is called from *call_free_submodel_memory*,
    ! src/socol_submodels.f90.

       DEALLOCATE(CO_emiss)

 END SUBROUTINE cleanup_synth_trac

   END MODULE mo_socol_synth_tracers

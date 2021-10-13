MODULE mo_socol_isotope
  
  ! Contains subroutines for diagnostic of CH4 isotopic composition
  ! Andrea Stenke, ETH Zurich, November 2009

  USE mo_advection,         ONLY: iadvec   ! selected advection scheme
  USE mo_kind,              ONLY: dp
  USE mo_memory_gl,         ONLY: xt       ! tracer field array
  USE mo_memory_g1a,        ONLY: xtm1     ! tracer field array (t-1)
  USE mo_semi_impl,         ONLY: eps
  USE mo_time_control,      ONLY: lresume
  USE mo_tracdef,           ONLY: trlist, t_trlist, t_trinfo
  USE mo_tracer,            ONLY: ntrac, & ! number of tracers
                                  new_tracer, & ! request resources for a tracer
                                  GAS, & ! phase indicators
                                  OFF, ON, &
                                  trlist, RESTART, INITIAL
  
  IMPLICIT NONE
  
  !--- Module variables:
  
  ! Vector with indices in order of initialization file:
  INTEGER, PUBLIC :: isotope_trac(2)
  
  ! Tracer indices (C-isotopes) [mass mixing ratio]:
  INTEGER, PUBLIC :: idt_12c,  idt_13c

  ! Specific isotopic signature of each emission type, delta13C [permille]
  ! values are taken from Lassey et al., ACP, 7, 2119-2139, 2007
  ! CH4 sources
  ! natural
  REAL, PARAMETER :: delta13C_wetlands   = -60.0_dp
  REAL, PARAMETER :: delta13C_termites   = -57.0_dp
  REAL, PARAMETER :: delta13C_wildfire   = -25.0_dp 
  REAL, PARAMETER :: delta13C_oceans     = -40.0_dp 
  REAL, PARAMETER :: delta13C_animals    = -62.0_dp
  REAL, PARAMETER :: delta13C_geologic   = -40.0_dp

  REAL, PARAMETER :: delta13C_natural    = -57.4_dp  

  ! anthropogenic
  REAL, PARAMETER :: delta13C_rice       = -64.0_dp
  REAL, PARAMETER :: delta13C_livestock  = -62.0_dp
  REAL, PARAMETER :: delta13C_bb         = -25.0_dp
  REAL, PARAMETER :: delta13C_savburn    = -12.0_dp
  REAL, PARAMETER :: delta13C_waste      = -55.0_dp
  REAL, PARAMETER :: delta13C_animalwaste= -55.0_dp
  REAL, PARAMETER :: delta13C_manure     = -55.0_dp
  REAL, PARAMETER :: delta13C_coal       = -35.0_dp
  REAL, PARAMETER :: delta13C_fossils    = -40.0_dp

  REAL, PARAMETER :: delta13C_anthrop    = -47.0_dp  

  ! kinetic isotope effects of chemical CH4 sinks
  REAL, PARAMETER :: kie_oh  = 1.0039_dp ! Saueressig et al., JGR, 106, 23127-23138, 2001
  REAL, PARAMETER :: kie_o1d = 1.013_dp ! Saueressig et al., JGR, 106, 23127-23138, 2001
  REAL, PARAMETER :: kie_cl  = 1.066_dp ! Crowley et al., Chemical Physics Letters, 303, 268-274, 1999

!!$  REAL, PARAMETER :: frac_cl_12c = kie_cl / (kie_cl + 1._dp) 
!!$  REAL, PARAMETER :: frac_cl_13c = 1._dp / (kie_cl + 1._dp) 
!!$  REAL, PARAMETER :: frac_oh_12c = kie_oh / (kie_oh + 1._dp) 
!!$  REAL, PARAMETER :: frac_oh_13c = 1._dp / (kie_oh + 1._dp) 
!!$  REAL, PARAMETER :: frac_o1d_12c = kie_o1d / (kie_o1d + 1._dp) 
!!$  REAL, PARAMETER :: frac_o1d_13c = 1._dp / (kie_o1d + 1._dp)

  REAL, PARAMETER :: frac_cl_12c = 1._dp / (kie_cl + 1._dp) 
  REAL, PARAMETER :: frac_cl_13c = kie_cl / (kie_cl + 1._dp) 
  REAL, PARAMETER :: frac_oh_12c = 1._dp / (kie_oh + 1._dp) 
  REAL, PARAMETER :: frac_oh_13c = kie_oh / (kie_oh + 1._dp) 
  REAL, PARAMETER :: frac_o1d_12c = 1._dp / (kie_o1d + 1._dp) 
  REAL, PARAMETER :: frac_o1d_13c = kie_o1d / (kie_o1d + 1._dp)

  ! 13C/12C ratio in the stable carbon isotope standard Vienne Peedee belemnite (VPDB)
  ! Lassey et al., ACP, 7, 2119-2139, 2007
  REAL, PARAMETER :: r_ref = 0.0112372 
  
  CONTAINS

    SUBROUTINE request_socol_isotope_tracers

      ! Defines 12C and 13C tracer for CH4 isotopic composition

      ! *request_socol_isotope_tracers* is called from *call_request_tracer*,
      ! src/call_submodels.f90.

      INTEGER :: trac_tran  ! switch for transport in ECHAM5
      INTEGER :: trac_vdiff ! switch for vertical diffusion in ECHAM5
      INTEGER :: trac_conv  ! switch for convection in ECHAM5

      trac_tran=iadvec
      trac_vdiff=ON
      trac_conv=ON

    CALL new_tracer('C12',     'MEZON', idx=idt_12c,     nwrite=ON,  &
         longname='12C-(CH4) volume mixing ratio',     units='mol/mol', &
         code=54, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)
    CALL new_tracer('C13',      'MEZON', idx=idt_13c,      nwrite=ON, &
         longname='13C-(CH4) volume mixing ratio',      units='mol/mol', &
         code=55, table=131, ninit=RESTART+INITIAL, nrerun=ON, &
         ntran=trac_tran, nvdiff=trac_vdiff, nconv=trac_conv, nphase=GAS)

    ! Vector with tracer indices in order of initialisation file:
    isotope_trac(1) = idt_12c
    isotope_trac(2) = idt_13c

  END SUBROUTINE request_socol_isotope_tracers

  !------------------------------------------------------------------------------

  SUBROUTINE init_socol_isotope_tracers

    ! Reads initial file of chemistry module MEZON if no restart file and
    ! initialize tracer tendencies due to chemistry.

    ! *init_socol_isotope_tracers* is called from *call_init_tracers*,
    ! src/call_submodels.f90.

    USE mo_socol_tracers, ONLY: idt_ch4

    IMPLICIT NONE

    INTEGER :: jt

    ! Executable statements:

    DO jt = 1, 2
       IF (trlist% ti(isotope_trac(jt))% init .EQ. 0 .AND. &
            IAND (trlist% ti(isotope_trac(jt))% ninit, INITIAL) /= 0) THEN

          ! Allocate field to xt:
          IF (jt == 1) THEN ! 12C
!!$             xt(:,:,isotope_trac(jt),:) = 0.99 * xt(:,:,idt_ch4,:)
             xt(:,:,isotope_trac(jt),:) = 0._dp
          ELSE ! 13C
!!$             xt(:,:,isotope_trac(jt),:) = 0.01 * xt(:,:,idt_ch4,:)
             xt(:,:,isotope_trac(jt),:) = 0._dp
          END IF
                
          IF (lresume) xtm1(:,:,isotope_trac(jt),:) = (1._dp - eps) &
               * xt(:,:,isotope_trac(jt),:)

          trlist% ti(isotope_trac(jt))% init = INITIAL

       ENDIF
    ENDDO

  END SUBROUTINE init_socol_isotope_tracers

  !------------------------------------------------------------------------------





END MODULE mo_socol_isotope

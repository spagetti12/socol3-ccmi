MODULE mo_socol_gcmfields

  ! Description:

  ! Contains fields calculated by the GCM for the use in MEZON.
  !
  ! Martin Schraner, ETH Zurich, April 2009

  USE mo_kind,             ONLY: dp
  USE mo_control,          ONLY: nlev, nlevp1
  USE mo_decomposition,    ONLY: ldc=>local_decomposition

  IMPLICIT NONE

  PUBLIC

  !--- Module variables:

  REAL(dp), ALLOCATABLE :: t(:,:)            ! Temperature [K]
  REAL(dp), ALLOCATABLE :: p(:,:)            ! Pressure [Pa] at full levels
  REAL(dp), ALLOCATABLE :: ph(:,:)           ! Pressure [Pa] at half levels
  REAL(dp), ALLOCATABLE :: prest(:,:)        ! Pressure [hPa] at full levels
  REAL(dp), ALLOCATABLE :: presb(:,:)        ! Pressure [hPa] at half levels

  LOGICAL, SAVE, PRIVATE :: lnot_used = .TRUE.
 
CONTAINS

  SUBROUTINE allocate_socol_gcmfields
    
    ! Allocates module variables.
    
    ! *allocate_socol_gcmfields* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90

    INTEGER :: nbdim
    
    IF (lnot_used) THEN
       
       nbdim = ldc%nproma
       
       ALLOCATE(t(nbdim,nlev))           ; t(:,:)      = 0.0_dp
       ALLOCATE(p(nbdim,nlev))           ; p(:,:)      = 0.0_dp
       ALLOCATE(ph(nbdim,nlevp1))        ; ph(:,:)     = 0.0_dp
       ALLOCATE(prest(nbdim,nlev))       ; prest(:,:)  = 0.0_dp
       ALLOCATE(presb(nbdim,nlevp1))     ; presb(:,:)  = 0.0_dp
       
       lnot_used = .FALSE.
       
    ENDIF
    
  END SUBROUTINE allocate_socol_gcmfields
  
  
  SUBROUTINE cleanup_socol_gcmfields
    
    ! Deallocates module variables.
    
    ! *cleanup_socol_gcmfields* is called from *call_free_submodel_memory*,
    ! src/socol_submodels.f90.
    
    IF (.NOT. lnot_used) THEN
       
       DEALLOCATE(t)
       DEALLOCATE(p)
       DEALLOCATE(ph)
       DEALLOCATE(prest)
       DEALLOCATE(presb)
       
       lnot_used = .TRUE.
       
    ENDIF
    
  END SUBROUTINE cleanup_socol_gcmfields
  
END MODULE mo_socol_gcmfields

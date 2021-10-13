MODULE mo_socol_no_gcm_condens_pscs

  ! This module is to avoid that the GCM performs condensation in the PSC
  ! regions.
  !
  ! M. Schraner, ETH Zurich, January 2010

  USE mo_control,                 ONLY: nlev
  USE mo_decomposition,           ONLY: ldc => local_decomposition
  USE mo_geoloc,                  ONLY: philat_2d
  USE mo_kind,                    ONLY: dp
  USE mo_socol_grid_calculations, ONLY: hetice_lowlevind
  USE mo_socol_namelist,          ONLY: hetice_uplev, hetice_north, &
                                        hetice_south
                                          
 
  IMPLICIT NONE

  PUBLIC

  ! Module variables:

  REAL(dp), ALLOCATABLE, PRIVATE :: pqte_befcond(:,:), pxlte_befcond(:,:), &
       pxite_befcond(:,:)
  
  LOGICAL, SAVE, PRIVATE :: lnot_used = .TRUE.

CONTAINS

  SUBROUTINE allocate_no_gcm_condens_pscs

    ! Allocate module variables.

    ! *allocate_no_gcm_condens_pscs* is called from *call_init_submodel_memory*,
    ! src/call_submodels.f90
    
    INTEGER :: nbdim
    
    IF (lnot_used) THEN
       
       nbdim = ldc%nproma
       
       ALLOCATE(pqte_befcond(nbdim,nlev))    ; pqte_befcond(:,:)  = 0.0_dp
       ALLOCATE(pxlte_befcond(nbdim,nlev))   ; pxlte_befcond(:,:) = 0.0_dp
       ALLOCATE(pxite_befcond(nbdim,nlev))   ; pxite_befcond(:,:) = 0.0_dp
       
       lnot_used = .FALSE.
       
    ENDIF
    
  END SUBROUTINE allocate_no_gcm_condens_pscs
  
  SUBROUTINE prepare_no_gcm_condens_pscs (kbdim, kproma, klev, pqte, pxlte, &
       pxite)

    ! Save values of qte, xlte and xite befcondore call of GCM condensation 
    ! subroutines.

    ! *prepare_no_gcm_condens_pscs* is called from *physc*.

    ! Subroutine arguments:
    INTEGER, INTENT(in) :: kbdim, kproma, klev
    REAL(dp), INTENT(in) :: pqte(kbdim,klev), pxlte(kbdim,klev), &
         pxite(kbdim,klev)


    ! Executable statements:
    
    pqte_befcond(1:kproma,:) = pqte(1:kproma,:)
    pxlte_befcond(1:kproma,:) = pxlte(1:kproma,:)
    pxite_befcond(1:kproma,:) = pxite(1:kproma,:)   

  END SUBROUTINE prepare_no_gcm_condens_pscs


  SUBROUTINE apply_no_gcm_condens_pscs (krow, kbdim, kproma, klev, papp1, &
       pqte, pxlte, pxite)
 
    ! Avoids condensation calculated by GCM in the region of the lower
    ! polar stratosphere.
    
    ! *avoid_gcm_condensation_psc* is called from *physc*.
 
    ! Subroutine arguments:
    INTEGER, INTENT(in) :: krow, kbdim, kproma, klev
    REAL(dp), INTENT(in) :: papp1(kbdim,klev)
    REAL(dp), INTENT(inout) :: pqte(kbdim,klev), pxlte(kbdim,klev), &
         pxite(kbdim,klev)

    ! Local variables:
    REAL(dp) :: hetice_uplev_pa
    INTEGER  :: jl, jk
 

    ! Executable statements:

    ! [hPa] -> [Pa]:
    hetice_uplev_pa  = 100._dp*hetice_uplev

    ! Keep tendencies calculated before condenstation in PSC region:
    DO jl = 1, kproma 
       IF (philat_2d(jl,krow) .GE. hetice_north .OR. &
            philat_2d(jl,krow) .LE. hetice_south) THEN
          DO jk = 1, nlev 
             IF (papp1(jl,jk) .GE. hetice_uplev_pa) THEN
                IF (jk .GT. hetice_lowlevind(jl)) EXIT
                pqte(jl,jk) = pqte_befcond(jl,jk)           
                pxlte(jl,jk) = pxlte_befcond(jl,jk)
                pxite(jl,jk) = pxite_befcond(jl,jk)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    
  END SUBROUTINE apply_no_gcm_condens_pscs

  SUBROUTINE cleanup_no_gcm_condens_pscs
    
    ! Deallocates module variables.
    
    ! *cleanup_no_gcm_condens_pscs* is called from *call_free_submodel_memory*,
    ! src/socol_submodels.f90.
    
    IF (.NOT. lnot_used) THEN
       
       DEALLOCATE(pqte_befcond)
       DEALLOCATE(pxlte_befcond)
       DEALLOCATE(pxite_befcond)
       
       lnot_used = .TRUE.
       
    ENDIF
    
  END SUBROUTINE cleanup_no_gcm_condens_pscs

END MODULE mo_socol_no_gcm_condens_pscs

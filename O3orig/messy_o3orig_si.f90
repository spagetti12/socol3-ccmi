#include "messy_main_ppd_si.inc"

! **********************************************************************
! INTERFACE TO ECHAM-5 FOR
! OZONE ORIGIN TRACERS, Details see messy_o3orig.f90 
! Version: see 'modver' in 
!
! Authors: Volker Grewe,    DLR,   April 2006        First Draft for QUANTIFY
!          Patrick Joeckel, MPICH, March 2006        Corrections
!          Stefanie Meul, 
!          Sophie Oberlaender, FUB, October 2010     for Messy 1.7 and 
!                                                    comparable to Grewe (2006)
!          Volker Grewe 
!          Patrick J�ckel,    DLR,   October 2010    inclusion in Messy1.11/2.41
!          Volker Grewe,      DLR,   April 2011      Tests for Messy2
!
! References:
!     Grewe, V., The origin of ozone, ACP 6, 1495-1511, 2006. 
!     http://www.atmos-chem-phys.net/6/1495/2006/acp-6-1495-2006.pdf
!
! TO DO
!     Redefine input-regions with all information on boxes and boundaries. 3D.
!     (i_trac_orig should be derived from input field)
!     more Integration schemes ? 
!     More flexibility in defining the regions, should be done via namelist
!           - 3D regions as in Grewe 2006
!           - use of  TP as boundary (as done) or not?
!           - Check whether errors are implemented correctly. pos. error still looks patchy
!           - ...
!     Remove some test-outputfields
! Tested: so far LINUX 1 Processor T21L19, only 
!
! **********************************************************************

! **********************************************************************
MODULE messy_o3orig_si
! **********************************************************************

  ! USE TOOLS
! op_pj_20101028+
!!$! ak_sm_11082010+
!!$  ! USE messy_main_tools_e5,   ONLY: start_message_e5, end_message_e5
!!$  USE messy_main_tools_bi,       ONLY: start_message_e5, end_message_e5
!!$! ak_sm_11082010-
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
! op_pj_20101028-

  USE messy_main_channel,    ONLY: t_chaobj_cpl 
  USE messy_o3orig              ! global data and public 

  IMPLICIT NONE
! PRIVATE

  ! ... CHANNEL OBJECTS: (DIAGNOSED FIELDS)
! Input to module
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3=>NULL()                ! kproma*nlev
  REAL(DP), DIMENSION(:,:),   POINTER :: tp=>NULL()                ! kproma*ngl    Tropopause in Pa
  REAL(DP), DIMENSION(:,:),   POINTER :: o3orig_regions=>NULL()    ! kproma*ngl
  REAL(DP), DIMENSION(:,:),   POINTER :: o3orig_regions2=>NULL()    ! kproma*ngl

! Newly defined in module?
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3prod=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3losst=>NULL(),o3lossd=>NULL()  ! kproma*nlev t=tracer d=diagnosed
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3m1=>NULL()              ! kproma*nlev
  REAL(DP), DIMENSION(:,:,:), POINTER :: co3=>NULL()              ! kproma*nlev

  REAL(DP), DIMENSION(:,:,:), POINTER :: o3_orig=>NULL(),o3_origm1=>NULL()          !nproma*nlev*n_trac_orig


  ! GLOBAL COUPLING SWITCHES
  ! um_ak_20110719+
  !CHARACTER(len=16),dimension(2) :: c_tropop
  TYPE(t_chaobj_cpl), SAVE :: c_tropop
  TYPE(t_chaobj_cpl), SAVE :: c_O3orig
  ! um_ak_20110719-

! tracer_numbers
  INTEGER, DIMENSION(n_trac_orig) :: itrac_o3orig                
! PTR TO OZONE TRACER/CHANNEL OBJECT

  INTEGER               :: idx_o3, idx_o3losst, idx_o3prod   !   Tracer identifier for tracer_gp channel

  ! -------------------------------------------------------------

  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: o3orig_initialize            ! global initialisation of module
  PUBLIC :: o3orig_new_tracer            ! define new tracers
  PUBLIC :: o3orig_init_memory           ! allocate memory and define channels
  PUBLIC :: o3orig_init_coupling         ! INitialize ...
  PUBLIC :: o3orig_physc                 ! integrate one time step
  PUBLIC :: o3orig_free_memory


CONTAINS

! ************************************************************************
! PUBLIC ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ------------------------------------------------------------------------
  SUBROUTINE  o3orig_initialize

    ! OZONE ORIGIN MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    !
    ! Author: Volker Grewe, DLR, March 2006

    ! ECHAM5
    ! op_pj_20100922+ cmodified for MESSy conformity
!!$    USE mo_mpi,                ONLY: p_parallel_io, p_io, p_bcast
!!$    USE mo_filename,           ONLY: find_next_free_unit
!!$    USE mo_exception,          ONLY: finish
    ! um_ak_20110719+ 
    !USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    ! um_ak_20110719-
    USE messy_main_tools,      ONLY: find_next_free_unit
    ! op_pj_20100922-

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'o3orig_initialize'
    INTEGER     :: iou    ! I/O unit
    INTEGER     :: status ! error status
    INTEGER     :: i_trac
   
    CALL start_message_bi(modstr, 'O3ORIG INITIALIZATION', substr)

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL o3orig_read_nml(status, iou)
       ! um_ak_20110719 IF (status /= 0) CALL finish(substr)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(i_integrate, p_io)
    CALL p_bcast(i_trac_orig, p_io)
    do i_trac=1,i_trac_orig
      CALL p_bcast(sn_o3orig(i_trac), p_io)
    enddo
    CALL p_bcast(l_err, p_io)

    if (i_trac_orig.gt.n_trac_orig) then
       ! um_ak_20110719 CALL finish(substr,'Number of regions exceeds limit')
       CALL error_bi('Number of regions exceeds limit', substr)
    endif

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! um_ak_20110719 CALL o3orig_read_nml_si(status, iou)
       CALL o3orig_read_nml_cpl(status, iou)
       ! um_ak_20110719 IF (status /= 0) CALL finish(substr,'Error in Reading Namelist CPL')
       IF (status /= 0) CALL error_bi('Error in Reading Namelist CPL', substr)
    END IF
    ! um_ak_20110719+
    CALL p_bcast(c_tropop%CHA, p_io)
    CALL p_bcast(c_tropop%OBJ, p_io)
    CALL p_bcast(c_o3orig%CHA, p_io)
    CALL p_bcast(c_o3orig%OBJ, p_io)
    ! um_ak_20110719-
    
    CALL end_message_bi(modstr, 'O3ORIG INITIALIZATION', substr)

!  
  END SUBROUTINE o3orig_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE o3orig_new_tracer

    ! ECHAM5/MESSy
! ak_sm_11082010+
!    USE messy_main_mpi_e5,        ONLY: p_parallel_io
! um_ak_20110719+
    !USE  messy_main_mpi_bi      , ONLY: p_parallel_io, finish ! op_pj_20100922
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: error_bi
! um_ak_20110719-
! ak_sm_11082010-
! ak_03082010+     change USE messy_main_tracer_mem_e5 into USE messy_main_tracer_mem_bi; AUG 3, 2010
!                  change USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, AIR, OFF, ON into
!                         USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
!                         USE messy_main_tracer       , ONLY: AIR, OFF, ON
!    USE messy_main_tracer_mem_e5, ONLY: GPTRSTR, AIR, OFF, ON !, LGTRSTR    
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
    USE messy_main_tracer_bi,    ONLY:  tracer_halt ! Error handling
    USE messy_main_tracer       , ONLY: new_tracer,set_tracer, R_MOLARMASS, R_HENRY, R_dryreac_sf,  &
                                                  I_DRYDEP, I_SCAV, I_WETDEP, AIR, OFF, ON !, LGTRSTR
! ak_03082010-
    USE messy_main_constants_mem, ONLY: MO
    ! MESSy
!!$    USE mo_exception,             ONLY: finish ! op_pj_20100922

    IMPLICIT NONE


    INTRINSIC   :: TRIM

    ! LOCAL
    INTEGER :: status,sum_status
    CHARACTER(LEN=*), PARAMETER :: substr = 'o3orig_new_tracer'
    INTEGER                     :: i_trac

!   IF (.NOT.lo3orig) RETURN     ! wo ist das definiert    ???VG

    CALL start_message_bi(modstr, 'TRACER REQUEST', substr)

       sum_status=0
!      write(*,*) 'start to define',i_trac_orig,'tracers for O3ORIG'
       do i_trac=1,i_trac_orig
 
         status=0
!        write(*,*) "i_trac=",i_trac,"Name:",TRIM(sn_o3orig(i_trac)),"#",status,sum_status
         CALL new_tracer(status, GPTRSTR, TRIM(sn_o3orig(i_trac)), modstr    &
               , unit='mol/mol', idx=itrac_o3orig(i_trac)                    &
               , longname='O3orig tracer - '//TRIM(sn_o3orig(i_trac))         &
                 )
!        write(*,*) "i_trac=",i_trac,TRIM(sn_o3orig(i_trac)),status,sum_status
         CALL tracer_halt(substr, status)
         if (sum_status/=0) sum_status=sum_status+1

         status=0
         CALL set_tracer(status, GPTRSTR, itrac_o3orig(i_trac)    &
                            , R_molarmass, r=MO*3._dp)
!        write(*,*) "i_trac=",i_trac,TRIM(sn_o3orig(i_trac)),status,sum_status
         CALL tracer_halt(substr, status)
         if (sum_status/=0) sum_status=sum_status+1

         status=0
         CALL set_tracer(status, GPTRSTR, itrac_o3orig(i_trac)    &
                            , I_Drydep, i=ON)
!        write(*,*) "i_trac=",i_trac,TRIM(sn_o3orig(i_trac)),status,sum_status
         CALL tracer_halt(substr, status)
         if (sum_status/=0) sum_status=sum_status+1
 
         status=0
         CALL set_tracer(status, GPTRSTR, itrac_o3orig(i_trac)    &
                            , I_wetdep,  i=OFF )
!         write(*,*) "i_trac=",i_trac,TRIM(sn_o3orig(i_trac)),status,sum_status
         CALL tracer_halt(substr, status)
         if (sum_status/=0) sum_status=sum_status+1

         status=0
         CALL set_tracer(status, GPTRSTR, itrac_o3orig(i_trac)    &
                            , I_scav, i=OFF  )
!         write(*,*) "i_trac=",i_trac,TRIM(sn_o3orig(i_trac)),status,sum_status
         CALL tracer_halt(substr, status)
         if (sum_status/=0) sum_status=sum_status+1

         status=0
         CALL set_tracer(status, GPTRSTR, itrac_o3orig(i_trac)    &
                            , R_HENRY, r=0.01_dp )
!         write(*,*) "i_trac=",i_trac,TRIM(sn_o3orig(i_trac)),status,sum_status
         CALL tracer_halt(substr, status)
         if (sum_status/=0) sum_status=sum_status+1

         status=0
         CALL set_tracer(status, GPTRSTR, itrac_o3orig(i_trac)    &
                            , R_dryreac_sf,  r= 1._dp     )
!         write(*,*) "i_trac=",i_trac,TRIM(sn_o3orig(i_trac)),status,sum_status
         CALL tracer_halt(substr, status)
         if (sum_status/=0) sum_status=sum_status+1
         
         
         IF (sum_status/=0) CALL error_bi(TRIM(sn_o3orig(i_trac))//' not allocated ',substr)
         IF (p_parallel_io) WRITE(*,*) ' ... New tracer ',TRIM(sn_o3orig(i_trac)),' defined'
       enddo

    CALL end_message_bi(modstr, 'TRACER REQUEST', substr)

  END SUBROUTINE o3orig_new_tracer
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE o3orig_init_memory
! Wird eigentlich nicht gebraucht Loss and Production sowie origin tracers �ber 
! Tracer channel definiert, weas braucht man snst - nichts. VG!
    ! Ozone Origin MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! define O3ORIG specific channel(s) and allocate memory for
    ! global fields
    !
    ! Author:  Volker Grewe, DLR, March 2006

    ! ECHAM5
! op_pj_20101028+
!!$    USE mo_memory_base,   ONLY: t_stream, new_stream         &
!!$                               ,add_stream_element           &
!!$                               ,default_stream_setting       
    USE messy_main_channel_bi, ONLY: channel_halt, GP_3D_MID, GP_2D_HORIZONTAL
    ! MESSy
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute
! op_pj_20101028-
! op_pj_20100922+ modified for MESSy conformity
!!$    USE mo_mpi,           ONLY: p_parallel_io
    USE messy_main_mpi_bi,  ONLY: p_parallel_io
    USE messy_main_data_bi,    ONLY: nproma, nlev
! op_pj_20100922-

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'o3orig_init_memory'
! op_pj_20101028+
!!$    ! DEFINE NEW DIAGNOSTIC STREAM (FOR OUTPUT) ...
!!$    TYPE (t_stream), POINTER            :: stream_ptr
    INTEGER :: status
! op_pj_20101028-

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

! op_pj_20101028+
!!$    ! define new stream for O3ORIG
!!$    CALL new_stream (stream_ptr,modstr)
!!$    CALL default_stream_setting(stream_ptr            &
!!$                                 ,units =''           &
!!$                                 ,lpost =.TRUE.       &
!!$                                 ,lrerun = .false.    &
!!$                                 ,lav = .true.        &
!!$                                 ,lsd = .true.        &
!!$                                )
    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
! op_pj_20101028-

    ! ELEMENTS USED FOR GP AND LG
    IF (p_parallel_io) WRITE(*,*) ' ... '

! op_pj_20101028+
!!$    CALL add_stream_element(stream_ptr      &
!!$         , 'O3Losst'                        &
!!$         , o3losst                          &
!!$         , longname='Ozone loss from tracer' &
!!$         , units = 'mol/mol/s'              &
!!$         )
!!$
!!$    IF (p_parallel_io) WRITE(*,*) ' ... O3Loss from tracer added to stream ',modstr

    CALL new_channel_object(status, modstr, 'O3Losst' &
         , p3 = o3losst )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3Losst'   &
         , 'long_name', c='Ozone loss from tracer' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3Losst'   &
         , 'units', c='mol/mol/s' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) ' ... O3Loss from tracer added to channel ',modstr
! op_pj_20101028-

! op_pj_20101028+
!!$    CALL add_stream_element(stream_ptr      &
!!$         , 'O3Lossd'                        &
!!$         , o3lossd                          &
!!$         , longname='Ozone loss diagnosed  ' &
!!$         , units = 'mol/mol/s'              &
!!$         )
!!$
!!$    IF (p_parallel_io) WRITE(*,*) ' ... O3Loss diagnosed   added to stream ',modstr

    CALL new_channel_object(status, modstr, 'O3Lossd' &
         , p3 = o3lossd )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3Lossd'   &
         , 'long_name', c='Ozone loss diagnosed' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3Lossd'   &
         , 'units', c='mol/mol/s' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) ' ... O3Loss diagnosed   added to channel ',modstr
! op_pj_20101028-

! op_pj_20101028+
!!$    CALL add_stream_element(stream_ptr      &
!!$         , 'O3Prod'                         &
!!$         , o3prod                           &
!!$         , longname='Ozone production     ' &
!!$         , units = 'mol/mol/s'              &
!!$         )
!!$
!!$    IF (p_parallel_io) WRITE(*,*) ' ... O3Prod added to stream ',modstr

    CALL new_channel_object(status, modstr, 'O3Prod' &
         , p3 = o3prod )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3Prod'   &
         , 'long_name', c='Ozone production' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3Prod'   &
         , 'units', c='mol/mol/s' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) ' ... O3Prod added to channel ',modstr
! op_pj_20101028-

! op_pj_20101028+
!!$    CALL add_stream_element(stream_ptr      &
!!$         , 'O3M1'                           &
!!$         , o3m1                             &
!!$         , longname='Pre-Chem ozone       ' &
!!$         , units = 'mol/mol'                &
!!$         )
!!$
!!$    IF (p_parallel_io) WRITE(*,*) ' ... O3M1   added to stream ',modstr

    CALL new_channel_object(status, modstr, 'O3M1' &
         , p3 = o3m1 )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3M1'   &
         , 'long_name', c='Pre-Chem ozone ' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3M1'   &
         , 'units', c='mol/mol' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) ' ... O3M1   added to channel ',modstr
! op_pj_20101028-

! op_pj_20101028+
!!$    CALL add_stream_element(stream_ptr      &
!!$         , 'O3'                             &
!!$         , o3                               &
!!$         , longname='Post-Chem ozone      ' &
!!$         , units = 'mol/mol'                &
!!$         )
!!$
!!$    IF (p_parallel_io) WRITE(*,*) ' ... O3     added to stream ',modstr

    CALL new_channel_object(status, modstr, 'O3' &
         , p3 = o3 )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3'   &
         , 'long_name', c='Post-Chem ozone' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'O3'   &
         , 'units', c='mol/mol' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) ' ... O3     added to channel ',modstr
! op_pj_20101028-

! op_pj_20101028+
!!$    CALL add_stream_element(stream_ptr      &
!!$         , 'cO3'                            &
!!$         , co3                              &
!!$         , longname='Chem ozone tendency  ' &
!!$         , units = 'mol/mol/s'              &
!!$         )
!!$
!!$    IF (p_parallel_io) WRITE(*,*) ' ... cO3    added to stream ',modstr

    CALL new_channel_object(status, modstr, 'cO3' &
         , p3 = co3 )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cO3'   &
         , 'long_name', c='Chem ozone tendency' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cO3'   &
         , 'units', c='mol/mol/s' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) ' ... cO3    added to channel ',modstr
! op_pj_20101028-

! op_pj_20101028+
!!$    CALL add_stream_element(stream_ptr      &
!!$         , 'Reg'                            &
!!$         , o3orig_regions2                  &
!!$         , longname='Regions              ' &
!!$         , units = 'none   '                &
!!$         )
!!$
!!$    IF (p_parallel_io) WRITE(*,*) ' ... Reg    added to stream ',modstr

    CALL new_channel_object(status, modstr, 'Reg' &
         , p2 = o3orig_regions2, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Reg'   &
         , 'long_name', c='Regions' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Reg'   &
         , 'units', c='none' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) ' ... Reg    added to channel ',modstr
! op_pj_20101028-

    ALLOCATE(o3_orig(nproma,nlev,n_trac_orig))
    ALLOCATE(o3_origm1(nproma,nlev,n_trac_orig))

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE o3orig_init_memory
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------

  SUBROUTINE o3orig_init_coupling

    ! ECHAM5
    ! op_pj_20100922+ modified for MESSy conformance
!!$    USE mo_mpi,                   ONLY: p_parallel_io, p_bcast
! op_pj_20101028+
!!$    USE mo_memory_base,           ONLY: t_stream, get_stream, get_stream_element
    USE messy_main_channel_bi,    ONLY: channel_halt
    USE messy_main_channel,       ONLY: get_channel_object
! op_pj_20101028-
!!$    USE mo_exception,             ONLY: finish
    ! um_ak_20110719+
    !USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_bcast, finish
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_bcast
    USE messy_main_blather_bi,    ONLY: error_bi
    ! um_ak_20110719-
    ! op_pj_20100922-
! ak_sm_11082010+ 
!    USE messy_main_tracer_mem_e5, ONLY: GPTRSTR
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
! ak_sm_11082010- 
    USE messy_main_tracer,        ONLY: get_tracer

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'o3orig_init_coupling'
    INTEGER                       :: ierr    ! error flag
!!$    TYPE (t_stream), POINTER      :: tropop, emission  ! op_pj_20101028

    ! PROCESS NAMELIST

       CALL start_message_bi(modstr, 'COUPLING INITIALIZING', substr)

       ierr=0
! Get Ozone
       CALL get_tracer(ierr, GPTRSTR, 'O3',idx=idx_O3)
       ! um_ak_20110719 if (ierr/=0) CALL FINISH(substr,'Tracer O3 not found')
       if (ierr/=0) CALL error_bi('Tracer O3 not found',substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fetched O3 from '


! Get Loss (Lost?)

       CALL get_tracer(ierr, GPTRSTR, 'LossO3',idx=idx_o3losst)
       ! um_ak_20110719if (ierr/=0) CALL FINISH(substr,'Tracer LossO3 not found')
       if (ierr/=0) CALL error_bi('Tracer LossO3 not found',substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fetched LossO3'

! Get Prod

       CALL get_tracer(ierr, GPTRSTR, 'ProdO3',idx=idx_o3prod)
       ! um_ak_20110719if (ierr/=0) CALL FINISH(substr,'Tracer PRODO3 not found')
       if (ierr/=0) CALL error_bi('Tracer PRODO3 not found',substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fetched ProdO3'
         

! Get Tropopause in [Pa]
! op_pj_20101028+
!!$       CALL get_stream(tropop,c_tropop(1),ierr) 
!!$       if (ierr/=0) CALL FINISH(substr,'Stream not found '//c_tropop(1))
!!$
!!$       CALL get_stream_element(tropop,c_tropop(2),tp,ierr)
!!$       if (ierr/=0) CALL FINISH(substr,'Stream element not found '//c_tropop(2))
! um_ak_20110719+
       !CALL get_channel_object(ierr, TRIM(c_tropop(1)), TRIM(c_tropop(2)) &
       !     , p2=tp)
       CALL get_channel_object(ierr, TRIM(c_tropop%CHA), TRIM(c_tropop%OBJ) &
            , p2=tp)
! um_ak_20110719-
       CALL channel_halt(substr, ierr)
! op_pj_20101028-
       IF (p_parallel_io) WRITE(*,*) ' ... fetched Tropop'

! Get region index field
!!$! op_pj_20101028+
!!$       CALL get_stream(emission,'offlem',ierr)
!!$       if (ierr/=0) CALL FINISH(substr,'Stream not found offlem')
!!$
!!$! op_pj_20100922+ more robust with name (see offlem.nml)
!!$       CALL get_stream_element(emission,'RGT0033_orig',o3orig_regions,ierr)
!!$       CALL get_stream_element(emission,'O3ORIG_orig_regions',o3orig_regions,ierr)
!!$! op_pj_20100922-
!!$       if (ierr/=0) CALL FINISH(substr,'Stream element O3ORIG_orig_regions not found o3orig')

! um_ak_20110719+
       !CALL get_channel_object(ierr, 'import_rgt', 'O3ORIG_orig' &
       !     , p2=o3orig_regions)
       CALL get_channel_object(ierr, TRIM(c_o3orig%CHA), TRIM(c_o3orig%OBJ) &
            , p2=o3orig_regions)
! um_ak_20110719-
       CALL channel_halt(substr, ierr)
! op_pj_20101028-
       IF (p_parallel_io) WRITE(*,*) ' ... fetched Regions'

       CALL end_message_bi(modstr, 'COUPLING INITIALIZING', substr)

  END SUBROUTINE o3orig_init_coupling

  !***************************************************************************


  SUBROUTINE o3orig_physc(i_o3orig_ctrl)

    ! ECHAM5
    ! op_pj_20100922+ modified for MESSy conformaty
!!$    USE mo_control,            ONLY: nrow,ngl,nlev,nlon
!!$    USE mo_decomposition,      ONLY: dcl => local_decomposition
!!$    USE mo_exception,          ONLY: finish
    ! um_ak_20110719+
    ! USE messy_main_mpi_bi,          ONLY: finish
    USE messy_main_blather_bi,    ONLY: error_bi, warning_bi
    ! um_ak_20110719-
    ! op_pj_20100922-
!ak_sm_11082010+
!    USE messy_main_data_e5,    ONLY: time_step_len,                &
!                                     jrow, kproma, nproma,         &
!                                     press_3d                         ! Better in init_coupling?
    USE messy_main_data_bi,    ONLY: jrow, kproma, nproma,         &
                                     nlev, nlon, & ! op_pj_20100922
                                     press_3d                         ! Better in init_coupling?
#ifndef MESSYTIMER
    USE messy_main_data_bi,    ONLY: time_step_len
#else
    USE messy_main_timer,      ONLY: time_step_len
#endif
!ak_sm_11082010-

! ak_sm_11082010+
!    USE messy_main_tracer_mem_e5,    ONLY: pxtte => qxtte, pxtm1 => qxtm1 &
!                                             , ntrac_gp, ti_gp
    USE messy_main_tracer_mem_bi,    ONLY: pxtte => qxtte, pxtm1 => qxtm1 &
                                             , ntrac_gp, ti_gp
!ak_sm_11082010-


    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)            :: i_o3orig_ctrl             ! 1: pre mecca; 2: post Mecca

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'o3orig_physc'

    INTEGER     :: status

    REAL(DP), DIMENSION (kproma,nlev)          ::  o3l, o3p

   ! one could try to use znrbc and zmrac from messy_mecca_e5.f90

!   REAL(DP), DIMENSION(kproma,nlev,n_trac_orig)     :: o3origM1  !sm_18082010 sm_02092010

    INTEGER                               :: jt
    INTEGER ::  idt ! um_ak_20110719

    status=0

    SELECT CASE(i_o3orig_ctrl)
    CASE(1)                    ! Pre-MECCA
                               ! Set ozone and tendencies
       ! um_ak_20110719+
       ! o3m1(1:kproma,:,jrow) = pxtm1(1:kproma,:,idx_O3) & 
       !          + pxtte(1:kproma,:,idx_O3) * time_step_len
       idt = idx_O3
       o3m1(1:kproma,_RI_YC_) = pxtm1(1:kproma,_RI_CD_) & 
                 + pxtte(1:kproma,_RI_CD_) * time_step_len
       ! um_ak_20110719-
        o3orig_regions2(1:kproma,jrow) = o3orig_regions(1:kproma,jrow)

    CASE(2)                    ! Post-MECCA
     
       ! um_ak_20110719+
!!$        o3(1:kproma,:,jrow)   = pxtm1(1:kproma,:,idx_O3) & 
!!$                   + pxtte(1:kproma,:,idx_O3) * time_step_len
!!$        co3(:,:,jrow)=(o3(:,:,jrow)-o3m1(:,:,jrow))/time_step_len  ! per second!
!!$        do jt=1,i_trac_orig
!!$           o3_orig(1:kproma,:,jt)  = pxtm1(1:kproma,:,itrac_o3orig(jt))  &
!!$              + pxtte(1:kproma,:,itrac_o3orig(jt)) * time_step_len
!!$        enddo
       idt = idx_O3
        o3(1:kproma,_RI_YC_)   = pxtm1(1:kproma,_RI_CD_) & 
                   + pxtte(1:kproma,_RI_CD_) * time_step_len
        co3(:,_RI_YC_)=(o3(:,_RI_YC_)-o3m1(:,_RI_YC_))/time_step_len  ! per second!
        do jt=1,i_trac_orig
           idt = itrac_o3orig(jt)
           o3_orig(1:kproma,:,jt)  = pxtm1(1:kproma,_RI_CD_)  &
              + pxtte(1:kproma,_RI_CD_) * time_step_len
        enddo
       ! um_ak_20110719-
       
        o3_origM1(:,:,:)=o3_orig(:,:,:)

       ! um_ak_20110719+
!!$        o3losst(1:kproma,:,jrow)   = pxtte(1:kproma,:,idx_O3losst) 
!!$        o3prod(1:kproma,:,jrow)    = pxtte(1:kproma,:,idx_O3prod) 
!!$! Can be via namelist, nyd 
!!$        o3l(1:kproma,:)    = pxtte(1:kproma,:,idx_O3losst)
!!$        o3p(1:kproma,:)    = pxtte(1:kproma,:,idx_O3prod)
        idt = idx_O3losst
        o3losst(1:kproma,_RI_YC_) = pxtte(1:kproma,_RI_CD_)
        o3l(1:kproma,:)           = pxtte(1:kproma,_RI_CD_)

        idt = idx_O3prod
        o3prod(1:kproma,_RI_YC_) = pxtte(1:kproma,_RI_CD_) 
        o3p(1:kproma,:)          = pxtte(1:kproma,_RI_CD_)
       ! um_ak_20110719-

        CALL o3orig_integrate(kproma,nlev,               & ! Input  (Grid)
             time_step_len,                              & ! Input  (Timestep)
             o3orig_regions(1:kproma,jrow),tp(1:kproma,jrow), & ! Input  (Regions; tropopause pressure)
! um_ak_20110719 press_3d(:,:,jrow),                         & ! Input  (Pressure)
             press_3d(1:kproma,_RI_YC_),                 & ! Input  (Pressure)
             o3l(1:kproma,:),o3p(1:kproma,:),            & ! Input  (Chem)
! um_ak_20110719 o3(:,:,jrow),o3m1(:,:,jrow),            & ! Input  (Chem)
             o3(1:kproma,_RI_YC_),o3m1(1:kproma,_RI_YC_),& ! Input  (Chem)
             o3_orig(1:kproma,:,:),                      & ! In/Output (O3Orig)
             status)
! um_ak_20110719+
#ifdef ECHAM5
        IF (status /= 0) CALL error_bi('error in o3orig_integrate', substr)
#endif
#ifdef COSMO
        IF (status /= 0) CALL warning_bi('error in o3orig_integrate', substr)
#endif
!         o3lossd(1:kproma,:,jrow)    = o3l(1:kproma,:)
!         o3prod(1:kproma,:,jrow)     = o3p(1:kproma,:)
         o3lossd(1:kproma,_RI_YC_)    = o3l(1:kproma,:)
         o3prod(1:kproma,_RI_YC_)     = o3p(1:kproma,:)
! um_ak_20110719-

! um_ak_20110719+
!         do jt=1,i_trac_orig
!          pxtte(1:kproma,:,itrac_o3orig(jt))=pxtte(1:kproma,:,itrac_o3orig(jt))  &
!             + (o3_orig(1:kproma,:,jt)-o3_origM1(1:kproma,:,jt)) / time_step_len
         do jt=1,i_trac_orig
            idt = itrac_o3orig(jt)
            pxtte(1:kproma,_RI_CD_)=pxtte(1:kproma,_RI_CD_)  &
             + (o3_orig(1:kproma,:,jt)-o3_origM1(1:kproma,:,jt)) / time_step_len
! um_ak_20110719-
         enddo
    CASE DEFAULT
         status=1
    ! um_ak_20110719 IF (status/=0) CALL FINISH(substr,'incorect value for i_o3orig_ctrl')
    IF (status/=0) THEN
       write (0,*) 'incorrect value for i_o3orig_ctrl ', i_o3orig_ctrl
       CALL error_bi(&
         'incorrect value for i_o3orig_ctrl',substr)
    ENDIF
    END SELECT

  END SUBROUTINE o3orig_physc
! ----------------------------------------------------------------------

! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  !um_ak_20110719 SUBROUTINE o3orig_read_nml_e5(status, iou)
  SUBROUTINE o3orig_read_nml_cpl(status, iou)

    ! O3ORIG MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5/ATTILA
    !
    ! Author: Volker Grewe, DLR, 2006

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
  
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    ! um_ak_20110719 CHARACTER(LEN=*), PARAMETER :: substr = 'o3orig_read_nml_e5'
    CHARACTER(LEN=*), PARAMETER :: substr = 'o3orig_read_nml_cpl'

    NAMELIST /CPL/ c_tropop, c_o3orig

    ! LOCAL
    LOGICAL          :: lex      ! file exists ?
    INTEGER          :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
       ! um_ak_20110719WRITE(*,*) 'C_tropop: ',C_tropop(1),C_tropop(2)
       WRITE(*,*) 'C_tropop: ',TRIM(C_tropop%CHA),' ',TRIM(C_tropop%obj)
!      WRITE(*,*) 'C_B: ',C_B(1),C_B(2)

       WRITE(*,*) 'C_o3orig: ',TRIM(C_o3orig%cha),' ',TRIM(c_o3orig%obj)
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE o3orig_read_nml_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE o3orig_free_memory

    ! O3ORIG MODULE ROUTINE 
    ! Author: Volker Grewe, DLR, 2010
  implicit none

  CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_free_memory'

  DEALLOCATE(o3_orig)
  DEALLOCATE(o3_origm1)

  END SUBROUTINE o3orig_free_memory
! **********************************************************************
END MODULE messy_o3orig_si
! **********************************************************************

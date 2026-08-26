module fks_model_state_module
  use process_dimensions, only: nexternal, validate_process_dimensions
  implicit none
  private

  double precision, pointer, public :: strong_coupling => null()
  double precision, public :: active_flavours = 0d0
  double precision, allocatable, public :: external_masses(:)
  logical, save :: initialized = .false.

  public :: initialize_fks_model_state, validate_fks_model_state

contains

  subroutine initialize_fks_model_state(coupling, flavours, masses)
    implicit none
    double precision, target, intent(inout) :: coupling
    double precision, intent(in) :: flavours, masses(:)

    call validate_process_dimensions()
    if (size(masses) /= nexternal) then
      write (*, '(a)') 'ERROR in fks_model_state_module: invalid model-data shape'
      stop 1
    end if

    strong_coupling => coupling
    active_flavours = flavours
    if (.not. allocated(external_masses)) allocate (external_masses(nexternal))
    external_masses = masses
    initialized = .true.
  end subroutine initialize_fks_model_state

  subroutine validate_fks_model_state()
    implicit none

    if (.not. initialized .or. .not. associated(strong_coupling) .or. &
        .not. allocated(external_masses)) then
      write (*, '(a)') 'ERROR in fks_model_state_module: state is not initialized'
      stop 1
    end if
  end subroutine validate_fks_model_state

end module fks_model_state_module

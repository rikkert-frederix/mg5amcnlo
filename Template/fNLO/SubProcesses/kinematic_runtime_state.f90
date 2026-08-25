module kinematic_runtime_state
  use process_dimensions, only: nexternal, validate_process_dimensions
  implicit none
  private

  logical, allocatable, public :: is_a_j(:)
  logical, allocatable, public :: is_a_lp(:)
  logical, allocatable, public :: is_a_lm(:)
  logical, allocatable, public :: is_a_ph(:)

  double precision, allocatable, public :: fxfx_ren_scales(:)
  double precision, allocatable, public :: fxfx_fac_scale(:)
  integer, public :: nfxfx_ren_scales = 0

  logical, public :: kinematic_state_initialized = .false.

  public :: init_kinematic_state
  public :: sync_kinematic_state
  public :: reset_kinematic_state
  public :: finalize_kinematic_state
  public :: validate_kinematic_state

contains

  subroutine init_kinematic_state()
    implicit none

    call validate_process_dimensions()
    if (kinematic_state_initialized) then
      call validate_kinematic_state()
      return
    end if

    allocate(is_a_j(nexternal), is_a_lp(nexternal))
    allocate(is_a_lm(nexternal), is_a_ph(nexternal))
    allocate(fxfx_ren_scales(0:nexternal), fxfx_fac_scale(2))
    kinematic_state_initialized = .true.
    call reset_kinematic_state()
  end subroutine init_kinematic_state


  subroutine sync_kinematic_state(is_a_j_in, is_a_lp_in, &
       is_a_lm_in, is_a_ph_in, fxfx_ren_scales_in, &
       nfxfx_ren_scales_in, fxfx_fac_scale_in)
    implicit none
    logical, intent(in) :: is_a_j_in(:), is_a_lp_in(:)
    logical, intent(in) :: is_a_lm_in(:), is_a_ph_in(:)
    double precision, intent(in) :: fxfx_ren_scales_in(0:)
    integer, intent(in) :: nfxfx_ren_scales_in
    double precision, intent(in) :: fxfx_fac_scale_in(:)

    call init_kinematic_state()
    if (size(is_a_j_in) /= nexternal .or. &
        size(is_a_lp_in) /= nexternal .or. &
        size(is_a_lm_in) /= nexternal .or. &
        size(is_a_ph_in) /= nexternal) then
      call fail_validation('particle-classification shape mismatch')
    end if
    if (ubound(fxfx_ren_scales_in, 1) /= nexternal) then
      call fail_validation('FxFx renormalization-scale shape mismatch')
    end if
    if (size(fxfx_fac_scale_in) /= 2) then
      call fail_validation('FxFx factorization-scale shape mismatch')
    end if
    if (nfxfx_ren_scales_in < 0 .or. &
        nfxfx_ren_scales_in > nexternal) then
      call fail_validation('invalid number of FxFx scales')
    end if

    is_a_j = is_a_j_in
    is_a_lp = is_a_lp_in
    is_a_lm = is_a_lm_in
    is_a_ph = is_a_ph_in
    fxfx_ren_scales = fxfx_ren_scales_in
    fxfx_fac_scale = fxfx_fac_scale_in
    nfxfx_ren_scales = nfxfx_ren_scales_in
    call validate_kinematic_state()
  end subroutine sync_kinematic_state


  subroutine reset_kinematic_state()
    implicit none

    if (.not. kinematic_state_initialized) then
      call fail_validation('cannot reset uninitialized state')
    end if
    is_a_j = .false.
    is_a_lp = .false.
    is_a_lm = .false.
    is_a_ph = .false.
    fxfx_ren_scales = 0d0
    fxfx_fac_scale = 0d0
    nfxfx_ren_scales = 0
  end subroutine reset_kinematic_state


  subroutine validate_kinematic_state()
    implicit none

    call validate_process_dimensions()
    if (.not. kinematic_state_initialized) then
      call fail_validation('state is not initialized')
    end if
    if (.not. allocated(is_a_j) .or. .not. allocated(is_a_lp) .or. &
        .not. allocated(is_a_lm) .or. .not. allocated(is_a_ph)) then
      call fail_validation('particle-classification storage is absent')
    end if
    if (.not. allocated(fxfx_ren_scales) .or. &
        .not. allocated(fxfx_fac_scale)) then
      call fail_validation('FxFx scale storage is absent')
    end if
    if (size(is_a_j) /= nexternal .or. size(is_a_lp) /= nexternal .or. &
        size(is_a_lm) /= nexternal .or. size(is_a_ph) /= nexternal) then
      call fail_validation('particle-classification storage is invalid')
    end if
    if (lbound(fxfx_ren_scales, 1) /= 0 .or. &
        ubound(fxfx_ren_scales, 1) /= nexternal .or. &
        size(fxfx_fac_scale) /= 2) then
      call fail_validation('FxFx scale storage is invalid')
    end if
    if (nfxfx_ren_scales < 0 .or. nfxfx_ren_scales > nexternal) then
      call fail_validation('stored number of FxFx scales is invalid')
    end if
  end subroutine validate_kinematic_state


  subroutine finalize_kinematic_state()
    implicit none

    if (allocated(is_a_j)) deallocate(is_a_j)
    if (allocated(is_a_lp)) deallocate(is_a_lp)
    if (allocated(is_a_lm)) deallocate(is_a_lm)
    if (allocated(is_a_ph)) deallocate(is_a_ph)
    if (allocated(fxfx_ren_scales)) deallocate(fxfx_ren_scales)
    if (allocated(fxfx_fac_scale)) deallocate(fxfx_fac_scale)
    nfxfx_ren_scales = 0
    kinematic_state_initialized = .false.
  end subroutine finalize_kinematic_state


  subroutine fail_validation(message)
    implicit none
    character(len=*), intent(in) :: message

    write(*, '(a)') 'ERROR in kinematic_runtime_state: ' // trim(message)
    stop 1
  end subroutine fail_validation

end module kinematic_runtime_state

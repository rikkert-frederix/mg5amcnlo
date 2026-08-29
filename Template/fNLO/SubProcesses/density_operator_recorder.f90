module density_operator_recorder
  use multiplicative_density_terms, only: block_distribution_term, &
       density_primitive_descriptor, density_scale_coefficient_count, &
       finalize_block_distribution_term
  use fnlo_process_common, only: soft_counterevent, real_event
  implicit none
  private

  ! The recorder is deliberately scoped to one mapped momentum point.  A
  ! caller must finish (or cancel) the current term before selecting another
  ! event slot.  In particular, real and counterevent operators can never be
  ! accumulated in this buffer and mistaken for one numerical density.
  logical, save :: recording = .false.
  integer, save :: recorded_block = -1
  integer, save :: recorded_event_slot = -1
  integer, save :: recorded_nlo_order = -1
  integer, save :: recorded_primitive_count = 0
  type(density_primitive_descriptor), allocatable, save :: &
       recorded_primitives(:)

  public :: begin_density_operator_recording
  public :: density_operator_is_recording
  public :: recorded_density_operator_count
  public :: recorded_density_operator_nonzero_count
  public :: record_density_operator_primitive
  public :: scale_recorded_density_operator
  public :: finish_density_operator_recording
  public :: cancel_density_operator_recording

contains

  subroutine begin_density_operator_recording(block, event_slot, nlo_order)
    integer, intent(in) :: block, event_slot, nlo_order

    if (recording) then
      call fail_density_operator_recorder( &
           'a same-momentum term is already being recorded')
    end if
    if (block < 0 .or. event_slot < 0) then
      call fail_density_operator_recorder( &
           'a density term has a negative block or event slot')
    end if
    if (nlo_order < 0 .or. nlo_order > 1) then
      call fail_density_operator_recorder( &
           'a density term must be LO or one-NLO-order')
    end if
    if (allocated(recorded_primitives)) deallocate(recorded_primitives)
    recorded_block = block
    recorded_event_slot = event_slot
    recorded_nlo_order = nlo_order
    recorded_primitive_count = 0
    recording = .true.
  end subroutine begin_density_operator_recording


  logical function density_operator_is_recording()
    density_operator_is_recording = recording
  end function density_operator_is_recording


  integer function recorded_density_operator_count()
    if (.not. recording) then
      recorded_density_operator_count = 0
    else
      recorded_density_operator_count = recorded_primitive_count
    end if
  end function recorded_density_operator_count


  integer function recorded_density_operator_nonzero_count()
    integer :: primitive

    recorded_density_operator_nonzero_count = 0
    if (.not. recording) return
    do primitive = 1, recorded_primitive_count
      if (any(recorded_primitives(primitive)%scale_coefficients /= &
              (0d0, 0d0))) then
        recorded_density_operator_nonzero_count = &
             recorded_density_operator_nonzero_count + 1
      end if
    end do
  end function recorded_density_operator_nonzero_count


  subroutine record_density_operator_primitive( &
       insertion_kind, insertion_identifier, insertion_rank, &
       correlation_leg, scale_coefficients, laurent_poles_cancelled)
    integer, intent(in) :: insertion_kind, insertion_identifier
    integer, intent(in) :: insertion_rank, correlation_leg
    complex(kind=8), intent(in) :: scale_coefficients( &
         density_scale_coefficient_count)
    logical, intent(in) :: laurent_poles_cancelled
    type(density_primitive_descriptor), allocatable :: enlarged(:)
    integer :: new_count

    call require_recording()
    new_count = recorded_primitive_count + 1
    allocate(enlarged(new_count))
    if (recorded_primitive_count > 0) then
      enlarged(1:recorded_primitive_count) = recorded_primitives
    end if
    enlarged(new_count)%insertion_kind = insertion_kind
    enlarged(new_count)%insertion_identifier = insertion_identifier
    enlarged(new_count)%insertion_rank = insertion_rank
    enlarged(new_count)%correlation_leg = correlation_leg
    enlarged(new_count)%nlo_order = recorded_nlo_order
    enlarged(new_count)%scale_coefficients = scale_coefficients
    enlarged(new_count)%laurent_poles_cancelled = &
         laurent_poles_cancelled
    call move_alloc(enlarged, recorded_primitives)
    recorded_primitive_count = new_count
  end subroutine record_density_operator_primitive


  subroutine scale_recorded_density_operator( &
       scale_coefficients, first_primitive, last_primitive)
    complex(kind=8), intent(in) :: scale_coefficients( &
         density_scale_coefficient_count)
    integer, intent(in), optional :: first_primitive, last_primitive
    complex(kind=8) :: original(density_scale_coefficient_count)
    complex(kind=8) :: quadratic_rr, quadratic_rf, quadratic_ff
    integer :: primitive, first, last

    call require_recording()
    if (recorded_primitive_count < 1) then
      call fail_density_operator_recorder( &
           'an empty density term cannot be scaled')
    end if
    first = 1
    last = recorded_primitive_count
    if (present(first_primitive)) first = first_primitive
    if (present(last_primitive)) last = last_primitive
    if (first < 1 .or. last > recorded_primitive_count .or. &
        first > last) then
      call fail_density_operator_recorder( &
           'a density-primitive scaling range is invalid')
    end if
    do primitive = first, last
      original = recorded_primitives(primitive)%scale_coefficients
      quadratic_rr = original(2)*scale_coefficients(2)
      quadratic_rf = original(2)*scale_coefficients(3) + &
           original(3)*scale_coefficients(2)
      quadratic_ff = original(3)*scale_coefficients(3)
      if (quadratic_rr /= (0d0, 0d0) .or. &
          quadratic_rf /= (0d0, 0d0) .or. &
          quadratic_ff /= (0d0, 0d0)) then
        call fail_density_operator_recorder( &
             'one block primitive would acquire quadratic scale logs')
      end if
      recorded_primitives(primitive)%scale_coefficients(1) = &
           original(1)*scale_coefficients(1)
      recorded_primitives(primitive)%scale_coefficients(2) = &
           original(1)*scale_coefficients(2) + &
           original(2)*scale_coefficients(1)
      recorded_primitives(primitive)%scale_coefficients(3) = &
           original(1)*scale_coefficients(3) + &
           original(3)*scale_coefficients(1)
    end do
  end subroutine scale_recorded_density_operator


  subroutine finish_density_operator_recording( &
       term, sign, luminosity_configuration)
    type(block_distribution_term), intent(out) :: term
    integer, intent(in) :: sign
    integer, intent(in), optional :: luminosity_configuration

    integer :: primitive, retained
    integer :: radiation_group

    call require_recording()
    if (abs(sign) /= 1) then
      call fail_density_operator_recorder( &
           'a recorded density term has an invalid sign')
    end if
    retained = recorded_density_operator_nonzero_count()
    if (retained < 1) then
      call fail_density_operator_recorder( &
           'a recorded density term contains no nonzero primitives')
    end if

    term%block = recorded_block
    term%event_slot = recorded_event_slot
    if (present(luminosity_configuration)) then
      if (luminosity_configuration < 0) then
        call fail_density_operator_recorder( &
             'a recorded density term has an invalid luminosity owner')
      end if
      term%luminosity_configuration = luminosity_configuration
    end if
    term%sign = sign
    term%nlo_order = recorded_nlo_order
    term%primitive_count = retained
    allocate(term%primitives(retained))
    radiation_group = 1
    if (recorded_event_slot == real_event .or. &
        (recorded_event_slot == soft_counterevent .and. sign < 0)) then
      radiation_group = 2
    else if (recorded_event_slot > soft_counterevent .and. &
             recorded_event_slot < real_event) then
      radiation_group = 3
    end if
    retained = 0
    do primitive = 1, recorded_primitive_count
      if (.not. any(recorded_primitives(primitive)%scale_coefficients /= &
                    (0d0, 0d0))) cycle
      retained = retained + 1
      term%primitives(retained) = recorded_primitives(primitive)
      term%primitives(retained)%radiation_group = radiation_group
    end do
    call finalize_block_distribution_term(term)
    call clear_recording()
  end subroutine finish_density_operator_recording


  subroutine cancel_density_operator_recording()
    call clear_recording()
  end subroutine cancel_density_operator_recording


  subroutine require_recording()
    if (.not. recording) then
      call fail_density_operator_recorder( &
           'no same-momentum density term is being recorded')
    end if
  end subroutine require_recording


  subroutine clear_recording()
    if (allocated(recorded_primitives)) deallocate(recorded_primitives)
    recorded_block = -1
    recorded_event_slot = -1
    recorded_nlo_order = -1
    recorded_primitive_count = 0
    recording = .false.
  end subroutine clear_recording


  subroutine fail_density_operator_recorder(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in density_operator_recorder: '//trim(message)
    stop 1
  end subroutine fail_density_operator_recorder
end module density_operator_recorder

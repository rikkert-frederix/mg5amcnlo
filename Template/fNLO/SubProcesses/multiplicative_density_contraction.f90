module multiplicative_density_contraction
  use process_dimensions, only: nexternal
  use spin_density_matrix_results, only: spin_density_no_insertion
  use multiplicative_density_terms, only: block_nlo_distribution, &
       multiplicative_density_tuple, density_cartesian_tuple_count, &
       decode_density_cartesian_tuple, &
       evaluate_density_primitive_coefficient
  use multiplicative_kinematics, only: realize_factorized_event_tuple
  implicit none
  private

  public :: contract_multiplicative_density_tuple

  interface
    integer function sdm_multiplicative_block_count()
    end function sdm_multiplicative_block_count

    integer function sdm_multiplicative_physical_block(position)
      integer, intent(in) :: position
    end function sdm_multiplicative_physical_block

    integer function sdm_multiplicative_component_position(block)
      integer, intent(in) :: block
    end function sdm_multiplicative_component_position

    subroutine sdm_multiplicative_contraction( &
         event_slots, insertion_kinds, insertion_ids, insertion_ranks, &
         correlation_legs, result, precision_asked, precision_found, &
         return_code)
      integer, intent(in) :: event_slots(*), insertion_kinds(*)
      integer, intent(in) :: insertion_ids(*), insertion_ranks(*)
      integer, intent(in) :: correlation_legs(*)
      complex(kind=8), intent(out) :: result
      double precision, intent(in) :: precision_asked
      double precision, intent(out) :: precision_found
      integer, intent(out) :: return_code
    end subroutine sdm_multiplicative_contraction
  end interface

contains

  subroutine contract_multiplicative_density_tuple( &
       distributions, tuple_index, logarithmic_mu2_r, &
       logarithmic_mu2_f, precision_asked, result, nlo_order, &
       event_slots, precision_found, return_code, coupling_rescaling, &
       already_realized)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    integer, intent(in) :: tuple_index
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: precision_asked
    complex(kind=8), intent(out) :: result
    integer, intent(out) :: nlo_order
    integer, intent(out) :: event_slots(0:nexternal)
    double precision, intent(out) :: precision_found
    integer, intent(out) :: return_code
    double precision, intent(in), optional :: coupling_rescaling(0:, 0:)
    logical, intent(in), optional :: already_realized
    logical :: kinematics_are_realized

    type(multiplicative_density_tuple) :: tuple
    integer :: block_count, distribution, physical_block, position
    integer, allocatable :: component_event_slots(:)
    integer, allocatable :: insertion_kinds(:), insertion_ids(:)
    integer, allocatable :: insertion_ranks(:), correlation_legs(:)
    logical, allocatable :: component_is_owned(:)
    logical :: pass

    if (tuple_index < 1 .or. &
        tuple_index > density_cartesian_tuple_count(distributions)) then
      call fail_density_contraction('a density tuple index is out of range')
    end if
    if (ubound(logarithmic_mu2_r, 1) < nexternal .or. &
        ubound(logarithmic_mu2_f, 1) < nexternal) then
      call fail_density_contraction( &
           'a block-local scale-log vector has the wrong size')
    end if
    if (present(coupling_rescaling)) then
      if (ubound(coupling_rescaling, 1) < nexternal .or. &
          ubound(coupling_rescaling, 2) < 1) then
        call fail_density_contraction( &
             'a block-local coupling-rescaling table has the wrong size')
      end if
    end if
    call decode_density_cartesian_tuple(distributions, tuple_index, tuple)
    event_slots = tuple%event_slots
    nlo_order = tuple%nlo_order

    ! This realizes only the mutually compatible boosted block caches.  The
    ! visible event is deliberately not assembled until after all matrix
    ! elements in this tuple have been evaluated.
    kinematics_are_realized = .false.
    if (present(already_realized)) then
      kinematics_are_realized = already_realized
    end if
    pass = .true.
    if (.not. kinematics_are_realized) then
      call realize_factorized_event_tuple(event_slots, pass)
    end if
    if (.not. pass) then
      result = (0d0, 0d0)
      precision_found = 0d0
      return_code = 0
      return
    end if

    block_count = sdm_multiplicative_block_count()
    if (block_count /= size(distributions)) then
      call fail_density_contraction( &
           'the density distributions do not cover every matrix block')
    end if
    if (allocated(component_event_slots)) deallocate(component_event_slots)
    if (allocated(insertion_kinds)) deallocate(insertion_kinds)
    if (allocated(insertion_ids)) deallocate(insertion_ids)
    if (allocated(insertion_ranks)) deallocate(insertion_ranks)
    if (allocated(correlation_legs)) deallocate(correlation_legs)
    if (allocated(component_is_owned)) deallocate(component_is_owned)
    allocate(component_event_slots(block_count))
    allocate(insertion_kinds(block_count))
    allocate(insertion_ids(block_count))
    allocate(insertion_ranks(block_count))
    allocate(correlation_legs(block_count))
    allocate(component_is_owned(block_count))
    component_is_owned = .false.
    insertion_kinds = spin_density_no_insertion
    insertion_ids = 0
    insertion_ranks = 0
    correlation_legs = 0
    do position = 1, block_count
      physical_block = sdm_multiplicative_physical_block(position)
      if (physical_block < 0 .or. physical_block > nexternal) then
        call fail_density_contraction( &
             'a generated component owns an invalid physical block')
      end if
      component_event_slots(position) = event_slots(physical_block)
    end do
    do distribution = 1, size(distributions)
      physical_block = distributions(distribution)%block
      position = sdm_multiplicative_component_position(physical_block)
      if (position < 1 .or. position > block_count .or. &
          component_is_owned(position)) then
        call fail_density_contraction( &
             'the density-to-component map is not one-to-one')
      end if
      component_is_owned(position) = .true.
    end do
    if (.not. all(component_is_owned)) then
      call fail_density_contraction( &
           'a matrix component has no block distribution')
    end if

    result = (0d0, 0d0)
    precision_found = 0d0
    return_code = 0
    call enumerate_primitives(1, (1d0, 0d0))
    result = tuple%sign*result

  contains

    recursive subroutine enumerate_primitives( &
         distribution_index, coefficient)
      integer, intent(in) :: distribution_index
      complex(kind=8), intent(in) :: coefficient
      integer :: term_index, primitive_index, component_position
      integer :: selected_block, local_return_code
      double precision :: local_precision, term_coupling_rescaling
      complex(kind=8) :: local_result, local_coefficient

      if (distribution_index > size(distributions)) then
        call sdm_multiplicative_contraction( &
             component_event_slots, insertion_kinds, insertion_ids, &
             insertion_ranks, correlation_legs, local_result, &
             precision_asked, local_precision, local_return_code)
        result = result + coefficient*local_result
        precision_found = max(precision_found, local_precision)
        return_code = max(return_code, local_return_code)
        return
      end if

      selected_block = distributions(distribution_index)%block
      component_position = &
           sdm_multiplicative_component_position(selected_block)
      term_index = tuple%term_indices(distribution_index)
      term_coupling_rescaling = 1d0
      if (present(coupling_rescaling)) then
        term_coupling_rescaling = coupling_rescaling( &
             selected_block, distributions(distribution_index)% &
             terms(term_index)%nlo_order)
      end if
      do primitive_index = 1, &
           distributions(distribution_index)%terms(term_index)% &
           primitive_count
        associate(primitive => distributions(distribution_index)% &
                  terms(term_index)%primitives(primitive_index))
          insertion_kinds(component_position) = primitive%insertion_kind
          insertion_ids(component_position) = &
               primitive%insertion_identifier
          insertion_ranks(component_position) = primitive%insertion_rank
          correlation_legs(component_position) = primitive%correlation_leg
          local_coefficient = evaluate_density_primitive_coefficient( &
               primitive, logarithmic_mu2_r(selected_block), &
               logarithmic_mu2_f(selected_block))
          call enumerate_primitives( &
               distribution_index + 1, coefficient*local_coefficient* &
               term_coupling_rescaling)
        end associate
      end do
    end subroutine enumerate_primitives
  end subroutine contract_multiplicative_density_tuple


  subroutine fail_density_contraction(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_density_contraction: '//trim(message)
    stop 1
  end subroutine fail_density_contraction
end module multiplicative_density_contraction

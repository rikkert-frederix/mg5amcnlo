module multiplicative_density_contraction
  use process_dimensions, only: nexternal
  use spin_density_matrix_results, only: spin_density_no_insertion, &
       spin_density_virtual_insertion
  use multiplicative_density_terms, only: block_nlo_distribution, &
       multiplicative_density_tuple, density_scale_coefficient_count
  use multiplicative_kinematics, only: realize_factorized_event_tuple
  use multiplicative_generated_metadata, only: &
       multiplicative_block_count, multiplicative_physical_blocks, &
       multiplicative_component_position
  implicit none
  private

  ! One tuple's matrix-valued coefficient basis. Provider calls populate the
  ! generated backing store once; scale points then form scalar coefficients
  ! and contract the already-loaded matrices.
  type, public :: multiplicative_density_basis
    logical :: prepared = .false.
    logical :: virtual_primitives_loaded = .false.
    integer :: block_count = 0
    integer :: maximum_primitives = 0
    integer :: nlo_order = 0
    integer :: sign = 1
    integer, allocatable :: event_slots(:)
    integer, allocatable :: component_event_slots(:)
    integer, allocatable :: primitive_counts(:)
    integer, allocatable :: component_nlo_orders(:)
    integer, allocatable :: insertion_kinds(:, :)
    integer, allocatable :: insertion_ids(:, :)
    integer, allocatable :: insertion_ranks(:, :)
    integer, allocatable :: correlation_legs(:, :)
    complex(kind=8), allocatable :: scale_coefficients(:, :, :)
    complex(kind=8), allocatable :: coefficients(:, :)
    integer :: scale_monomial_count = 0
    integer, allocatable :: scale_mode_counts(:)
    integer, allocatable :: scale_modes(:, :)
    integer, allocatable :: scale_monomial_modes(:, :)
    complex(kind=8), allocatable :: full_scale_polynomial(:)
    complex(kind=8), allocatable :: novirtual_scale_polynomial(:)
    logical :: full_scale_polynomial_prepared = .false.
    logical :: novirtual_scale_polynomial_prepared = .false.
    logical, allocatable :: component_is_owned(:)
    double precision :: precision_found = 0d0
    integer :: return_code = 0
  end type multiplicative_density_basis

  public :: prepare_multiplicative_density_basis
  public :: evaluate_multiplicative_density_basis
  public :: prepare_multiplicative_scale_polynomial
  public :: evaluate_multiplicative_scale_polynomial

  interface
    subroutine sdm_multiplicative_prepare_basis( &
         event_slots, maximum_primitives, primitive_counts, &
         insertion_kinds, insertion_ids, insertion_ranks, &
         correlation_legs, include_virtual, precision_asked, &
         precision_found, return_code)
      integer, intent(in) :: event_slots(*), maximum_primitives
      integer, intent(in) :: primitive_counts(*)
      integer, intent(in) :: insertion_kinds(maximum_primitives, *)
      integer, intent(in) :: insertion_ids(maximum_primitives, *)
      integer, intent(in) :: insertion_ranks(maximum_primitives, *)
      integer, intent(in) :: correlation_legs(maximum_primitives, *)
      logical, intent(in) :: include_virtual
      double precision, intent(in) :: precision_asked
      double precision, intent(out) :: precision_found
      integer, intent(out) :: return_code
    end subroutine sdm_multiplicative_prepare_basis

    subroutine sdm_multiplicative_evaluate_basis( &
         maximum_primitives, primitive_counts, coefficients, result)
      integer, intent(in) :: maximum_primitives
      integer, intent(in) :: primitive_counts(*)
      complex(kind=8), intent(in) :: coefficients(maximum_primitives, *)
      complex(kind=8), intent(out) :: result
    end subroutine sdm_multiplicative_evaluate_basis
  end interface

contains

  subroutine prepare_multiplicative_density_basis( &
       distributions, tuple, precision_asked, basis, already_realized, &
       include_virtual)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(multiplicative_density_tuple), intent(in) :: tuple
    double precision, intent(in) :: precision_asked
    type(multiplicative_density_basis), intent(inout) :: basis
    logical, intent(in), optional :: already_realized
    logical, intent(in), optional :: include_virtual
    logical :: kinematics_are_realized, pass
    logical :: load_virtual_primitives
    integer :: block_count, distribution, physical_block, position
    integer :: primitive_index, term_index, maximum_primitives

    call validate_tuple_shape(distributions, tuple)
    basis%prepared = .false.
    load_virtual_primitives = .true.
    if (present(include_virtual)) load_virtual_primitives = include_virtual
    basis%virtual_primitives_loaded = load_virtual_primitives
    if (.not. allocated(basis%event_slots)) then
      allocate(basis%event_slots(0:nexternal))
    else if (lbound(basis%event_slots, 1) /= 0 .or. &
             ubound(basis%event_slots, 1) /= nexternal) then
      deallocate(basis%event_slots)
      allocate(basis%event_slots(0:nexternal))
    end if
    basis%event_slots = tuple%event_slots
    basis%nlo_order = tuple%nlo_order
    basis%sign = tuple%sign
    basis%precision_found = 0d0
    basis%return_code = 0

    ! Matrix elements consume only boosted block-local momenta. Visible-event
    ! assembly remains a later, independent operation for cuts and analysis.
    kinematics_are_realized = .false.
    if (present(already_realized)) &
         kinematics_are_realized = already_realized
    pass = .true.
    if (.not. kinematics_are_realized) then
      call realize_factorized_event_tuple(basis%event_slots, pass)
    end if
    if (.not. pass) return

    block_count = multiplicative_block_count
    if (block_count /= size(distributions)) then
      call fail_density_contraction( &
           'the density distributions do not cover every matrix block')
    end if
    maximum_primitives = 0
    do distribution = 1, size(distributions)
      term_index = tuple%term_indices(distribution)
      if (term_index < 1 .or. &
          term_index > distributions(distribution)%term_count) then
        call fail_density_contraction( &
             'a density tuple selects an invalid block term')
      end if
      maximum_primitives = max(maximum_primitives, &
           distributions(distribution)%terms(term_index)%primitive_count)
    end do
    if (maximum_primitives < 1) then
      call fail_density_contraction('a density tuple has no primitives')
    end if
    call ensure_basis_storage(basis, block_count, maximum_primitives)
    basis%component_is_owned = .false.
    basis%insertion_kinds = spin_density_no_insertion
    basis%insertion_ids = 0
    basis%insertion_ranks = 0
    basis%correlation_legs = 0
    basis%primitive_counts = 0
    basis%component_nlo_orders = 0
    basis%scale_coefficients = (0d0, 0d0)
    basis%coefficients = (0d0, 0d0)
    basis%full_scale_polynomial_prepared = .false.
    basis%novirtual_scale_polynomial_prepared = .false.

    do position = 1, block_count
      physical_block = multiplicative_physical_blocks(position)
      if (physical_block < 0 .or. physical_block > nexternal) then
        call fail_density_contraction( &
             'a generated component owns an invalid physical block')
      end if
      basis%component_event_slots(position) = &
           basis%event_slots(physical_block)
    end do
    do distribution = 1, size(distributions)
      physical_block = distributions(distribution)%block
      position = multiplicative_component_position(physical_block)
      if (position < 1 .or. position > block_count .or. &
          basis%component_is_owned(position)) then
        call fail_density_contraction( &
             'the density-to-component map is not one-to-one')
      end if
      basis%component_is_owned(position) = .true.
      term_index = tuple%term_indices(distribution)
      basis%primitive_counts(position) = distributions(distribution)% &
           terms(term_index)%primitive_count
      basis%component_nlo_orders(position) = distributions(distribution)% &
           terms(term_index)%nlo_order
      do primitive_index = 1, basis%primitive_counts(position)
        associate(primitive => distributions(distribution)% &
                  terms(term_index)%primitives(primitive_index))
          basis%insertion_kinds(primitive_index, position) = &
               primitive%insertion_kind
          basis%insertion_ids(primitive_index, position) = &
               primitive%insertion_identifier
          basis%insertion_ranks(primitive_index, position) = &
               primitive%insertion_rank
          basis%correlation_legs(primitive_index, position) = &
               primitive%correlation_leg
          basis%scale_coefficients(:, primitive_index, position) = &
               primitive%scale_coefficients
        end associate
      end do
    end do
    if (.not. all(basis%component_is_owned)) then
      call fail_density_contraction( &
           'a matrix component has no block distribution')
    end if
    call initialize_scale_polynomial_layout(basis)

    call sdm_multiplicative_prepare_basis( &
         basis%component_event_slots, basis%maximum_primitives, &
         basis%primitive_counts, basis%insertion_kinds, &
         basis%insertion_ids, basis%insertion_ranks, &
         basis%correlation_legs, load_virtual_primitives, precision_asked, &
         basis%precision_found, basis%return_code)
    basis%prepared = .true.
  end subroutine prepare_multiplicative_density_basis


  subroutine evaluate_multiplicative_density_basis( &
       basis, logarithmic_mu2_r, logarithmic_mu2_f, result, &
       coupling_rescaling, include_virtual)
    type(multiplicative_density_basis), intent(inout) :: basis
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    complex(kind=8), intent(out) :: result
    double precision, intent(in), optional :: coupling_rescaling(0:, 0:)
    logical, intent(in), optional :: include_virtual
    double precision :: coupling_factor
    integer :: position, physical_block, primitive
    logical :: use_virtual_primitives

    if (.not. basis%prepared) then
      call fail_density_contraction('an unprepared density basis was used')
    end if
    use_virtual_primitives = .true.
    if (present(include_virtual)) use_virtual_primitives = include_virtual
    if (use_virtual_primitives .and. &
        .not. basis%virtual_primitives_loaded) then
      do position = 1, basis%block_count
        do primitive = 1, basis%primitive_counts(position)
          if (basis%insertion_kinds(primitive, position) == &
              spin_density_virtual_insertion) then
            call fail_density_contraction( &
                 'a virtual primitive was not loaded in this basis')
          end if
        end do
      end do
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

    basis%coefficients = (0d0, 0d0)
    do position = 1, basis%block_count
      physical_block = multiplicative_physical_blocks(position)
      coupling_factor = 1d0
      if (present(coupling_rescaling)) then
        coupling_factor = coupling_rescaling( &
             physical_block, basis%component_nlo_orders(position))
      end if
      do primitive = 1, basis%primitive_counts(position)
        if (.not. use_virtual_primitives .and. &
            basis%insertion_kinds(primitive, position) == &
            spin_density_virtual_insertion) cycle
        basis%coefficients(primitive, position) = coupling_factor*( &
             basis%scale_coefficients(1, primitive, position) + &
             basis%scale_coefficients(2, primitive, position)* &
             logarithmic_mu2_r(physical_block) + &
             basis%scale_coefficients(3, primitive, position)* &
             logarithmic_mu2_f(physical_block))
      end do
    end do
    call sdm_multiplicative_evaluate_basis( &
         basis%maximum_primitives, basis%primitive_counts, &
         basis%coefficients, result)
    result = basis%sign*result
  end subroutine evaluate_multiplicative_density_basis


  subroutine prepare_multiplicative_scale_polynomial( &
       basis, include_virtual)
    type(multiplicative_density_basis), intent(inout) :: basis
    logical, intent(in), optional :: include_virtual
    logical :: use_virtual_primitives
    integer :: monomial, position, primitive, mode
    complex(kind=8) :: coefficient

    if (.not. basis%prepared) then
      call fail_density_contraction( &
           'an unprepared density basis was used for scale algebra')
    end if
    use_virtual_primitives = .true.
    if (present(include_virtual)) use_virtual_primitives = include_virtual
    if (use_virtual_primitives .and. &
        basis%full_scale_polynomial_prepared) return
    if (.not. use_virtual_primitives .and. &
        basis%novirtual_scale_polynomial_prepared) return
    if (use_virtual_primitives .and. &
        .not. basis%virtual_primitives_loaded) then
      do position = 1, basis%block_count
        do primitive = 1, basis%primitive_counts(position)
          if (basis%insertion_kinds(primitive, position) == &
              spin_density_virtual_insertion) then
            call fail_density_contraction( &
                 'a virtual scale polynomial was requested without its matrices')
          end if
        end do
      end do
    end if

    if (basis%scale_monomial_count == 0) then
      if (use_virtual_primitives) then
        basis%full_scale_polynomial_prepared = .true.
      else
        basis%novirtual_scale_polynomial_prepared = .true.
      end if
      return
    end if
    basis%coefficients = (0d0, 0d0)
    do monomial = 1, basis%scale_monomial_count
      do position = 1, basis%block_count
        mode = basis%scale_monomial_modes(position, monomial)
        do primitive = 1, basis%primitive_counts(position)
          if (.not. use_virtual_primitives .and. &
              basis%insertion_kinds(primitive, position) == &
              spin_density_virtual_insertion) then
            coefficient = (0d0, 0d0)
          else
            coefficient = basis%scale_coefficients( &
                 mode, primitive, position)
          end if
          basis%coefficients(primitive, position) = coefficient
        end do
      end do
      if (use_virtual_primitives) then
        call sdm_multiplicative_evaluate_basis( &
             basis%maximum_primitives, basis%primitive_counts, &
             basis%coefficients, basis%full_scale_polynomial(monomial))
        basis%full_scale_polynomial(monomial) = basis%sign* &
             basis%full_scale_polynomial(monomial)
      else
        call sdm_multiplicative_evaluate_basis( &
             basis%maximum_primitives, basis%primitive_counts, &
             basis%coefficients, basis%novirtual_scale_polynomial(monomial))
        basis%novirtual_scale_polynomial(monomial) = basis%sign* &
             basis%novirtual_scale_polynomial(monomial)
      end if
    end do
    if (use_virtual_primitives) then
      basis%full_scale_polynomial_prepared = .true.
    else
      basis%novirtual_scale_polynomial_prepared = .true.
    end if
  end subroutine prepare_multiplicative_scale_polynomial


  subroutine evaluate_multiplicative_scale_polynomial( &
       basis, logarithmic_mu2_r, logarithmic_mu2_f, result, &
       coupling_rescaling, include_virtual)
    type(multiplicative_density_basis), intent(inout) :: basis
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    complex(kind=8), intent(out) :: result
    double precision, intent(in), optional :: coupling_rescaling(0:, 0:)
    logical, intent(in), optional :: include_virtual
    logical :: use_virtual_primitives
    integer :: monomial, position, physical_block, mode
    double precision :: coupling_factor, monomial_factor

    use_virtual_primitives = .true.
    if (present(include_virtual)) use_virtual_primitives = include_virtual
    if (ubound(logarithmic_mu2_r, 1) < nexternal .or. &
        ubound(logarithmic_mu2_f, 1) < nexternal) then
      call fail_density_contraction( &
           'a scale-polynomial log vector has the wrong size')
    end if
    if (present(coupling_rescaling)) then
      if (ubound(coupling_rescaling, 1) < nexternal .or. &
          ubound(coupling_rescaling, 2) < 1) then
        call fail_density_contraction( &
             'a scale-polynomial coupling table has the wrong size')
      end if
    end if
    call prepare_multiplicative_scale_polynomial( &
         basis, use_virtual_primitives)

    coupling_factor = 1d0
    do position = 1, basis%block_count
      physical_block = multiplicative_physical_blocks(position)
      if (present(coupling_rescaling)) then
        coupling_factor = coupling_factor*coupling_rescaling( &
             physical_block, basis%component_nlo_orders(position))
      end if
    end do
    result = (0d0, 0d0)
    do monomial = 1, basis%scale_monomial_count
      monomial_factor = 1d0
      do position = 1, basis%block_count
        physical_block = multiplicative_physical_blocks(position)
        mode = basis%scale_monomial_modes(position, monomial)
        select case (mode)
        case (1)
        case (2)
          monomial_factor = monomial_factor* &
               logarithmic_mu2_r(physical_block)
        case (3)
          monomial_factor = monomial_factor* &
               logarithmic_mu2_f(physical_block)
        case default
          call fail_density_contraction( &
               'a scale polynomial contains an invalid logarithm mode')
        end select
      end do
      if (use_virtual_primitives) then
        result = result + monomial_factor* &
             basis%full_scale_polynomial(monomial)
      else
        result = result + monomial_factor* &
             basis%novirtual_scale_polynomial(monomial)
      end if
    end do
    result = coupling_factor*result
  end subroutine evaluate_multiplicative_scale_polynomial


  subroutine initialize_scale_polynomial_layout(basis)
    type(multiplicative_density_basis), intent(inout) :: basis
    integer :: position, primitive, mode, count, monomial, remainder

    if (.not. allocated(basis%scale_mode_counts)) &
         allocate(basis%scale_mode_counts(basis%block_count))
    if (.not. allocated(basis%scale_modes)) &
         allocate(basis%scale_modes(density_scale_coefficient_count, &
                                    basis%block_count))
    basis%scale_mode_counts = 0
    basis%scale_modes = 0
    basis%scale_monomial_count = 1
    do position = 1, basis%block_count
      count = 0
      do mode = 1, density_scale_coefficient_count
        do primitive = 1, basis%primitive_counts(position)
          if (basis%scale_coefficients(mode, primitive, position) /= &
              (0d0, 0d0)) then
            count = count + 1
            basis%scale_modes(count, position) = mode
            exit
          end if
        end do
      end do
      basis%scale_mode_counts(position) = count
      if (count == 0) then
        basis%scale_monomial_count = 0
      else if (basis%scale_monomial_count > 0) then
        basis%scale_monomial_count = basis%scale_monomial_count*count
      end if
    end do
    if (allocated(basis%scale_monomial_modes)) then
      if (size(basis%scale_monomial_modes, 1) /= basis%block_count .or. &
          size(basis%scale_monomial_modes, 2) /= &
          basis%scale_monomial_count) &
           deallocate(basis%scale_monomial_modes)
    end if
    if (.not. allocated(basis%scale_monomial_modes)) &
         allocate(basis%scale_monomial_modes( &
              basis%block_count, basis%scale_monomial_count))
    do monomial = 1, basis%scale_monomial_count
      remainder = monomial - 1
      do position = 1, basis%block_count
        count = basis%scale_mode_counts(position)
        basis%scale_monomial_modes(position, monomial) = &
             basis%scale_modes(mod(remainder, count) + 1, position)
        remainder = remainder/count
      end do
    end do
    call ensure_complex_vector( &
         basis%full_scale_polynomial, basis%scale_monomial_count)
    call ensure_complex_vector( &
         basis%novirtual_scale_polynomial, basis%scale_monomial_count)
    basis%full_scale_polynomial = (0d0, 0d0)
    basis%novirtual_scale_polynomial = (0d0, 0d0)
  end subroutine initialize_scale_polynomial_layout


  subroutine ensure_complex_vector(values, count)
    complex(kind=8), allocatable, intent(inout) :: values(:)
    integer, intent(in) :: count

    if (allocated(values)) then
      if (size(values) /= count) deallocate(values)
    end if
    if (.not. allocated(values)) allocate(values(count))
  end subroutine ensure_complex_vector


  subroutine ensure_basis_storage(basis, block_count, maximum_primitives)
    type(multiplicative_density_basis), intent(inout) :: basis
    integer, intent(in) :: block_count, maximum_primitives
    logical :: wrong_shape

    wrong_shape = basis%block_count /= block_count .or. &
         basis%maximum_primitives /= maximum_primitives
    if (wrong_shape) then
      if (allocated(basis%component_event_slots)) &
           deallocate(basis%component_event_slots)
      if (allocated(basis%primitive_counts)) &
           deallocate(basis%primitive_counts)
      if (allocated(basis%component_nlo_orders)) &
           deallocate(basis%component_nlo_orders)
      if (allocated(basis%insertion_kinds)) &
           deallocate(basis%insertion_kinds)
      if (allocated(basis%insertion_ids)) deallocate(basis%insertion_ids)
      if (allocated(basis%insertion_ranks)) &
           deallocate(basis%insertion_ranks)
      if (allocated(basis%correlation_legs)) &
           deallocate(basis%correlation_legs)
      if (allocated(basis%scale_coefficients)) &
           deallocate(basis%scale_coefficients)
      if (allocated(basis%coefficients)) deallocate(basis%coefficients)
      if (allocated(basis%scale_mode_counts)) &
           deallocate(basis%scale_mode_counts)
      if (allocated(basis%scale_modes)) deallocate(basis%scale_modes)
      if (allocated(basis%scale_monomial_modes)) &
           deallocate(basis%scale_monomial_modes)
      if (allocated(basis%full_scale_polynomial)) &
           deallocate(basis%full_scale_polynomial)
      if (allocated(basis%novirtual_scale_polynomial)) &
           deallocate(basis%novirtual_scale_polynomial)
      if (allocated(basis%component_is_owned)) &
           deallocate(basis%component_is_owned)
    end if
    if (.not. allocated(basis%component_event_slots)) then
      allocate(basis%component_event_slots(block_count))
      allocate(basis%primitive_counts(block_count))
      allocate(basis%component_nlo_orders(block_count))
      allocate(basis%insertion_kinds(maximum_primitives, block_count))
      allocate(basis%insertion_ids(maximum_primitives, block_count))
      allocate(basis%insertion_ranks(maximum_primitives, block_count))
      allocate(basis%correlation_legs(maximum_primitives, block_count))
      allocate(basis%scale_coefficients( &
           density_scale_coefficient_count, maximum_primitives, block_count))
      allocate(basis%coefficients(maximum_primitives, block_count))
      allocate(basis%component_is_owned(block_count))
    end if
    basis%block_count = block_count
    basis%maximum_primitives = maximum_primitives
  end subroutine ensure_basis_storage


  subroutine validate_tuple_shape(distributions, tuple)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(multiplicative_density_tuple), intent(in) :: tuple

    if (tuple%distribution_count /= size(distributions) .or. &
        .not. allocated(tuple%term_indices) .or. &
        size(tuple%term_indices) /= size(distributions) .or. &
        .not. allocated(tuple%event_slots) .or. &
        lbound(tuple%event_slots, 1) /= 0 .or. &
        ubound(tuple%event_slots, 1) < nexternal) then
      call fail_density_contraction('a density tuple has the wrong shape')
    end if
  end subroutine validate_tuple_shape


  subroutine fail_density_contraction(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_density_contraction: '//trim(message)
    stop 1
  end subroutine fail_density_contraction
end module multiplicative_density_contraction

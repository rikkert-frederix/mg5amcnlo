module spin_density_fks_matrices
  use process_dimensions, only: amp_split_size, amp_split_orders, &
       nsplitorders, qcd_pos
  implicit none
  private

  logical, save :: initialized = .false.
  logical, save :: collection_enabled = .false.
  integer, save :: maximum_open_size = 0
  integer, save :: active_open_size = 0
  integer, save :: active_component_position = 0
  integer, save :: born_amp_position = 0
  integer, save :: nlo_amp_position = 0
  logical, save :: born_matrix_is_available = .false.
  logical, save :: real_matrix_is_available = .false.
  logical, save :: color_matrix_is_available = .false.
  logical, save :: virtual_matrix_is_available = .false.
  complex(kind=8), allocatable, save :: born_density(:, :, :)
  complex(kind=8), allocatable, save :: real_density(:, :, :)
  complex(kind=8), allocatable, save :: color_density(:, :)
  complex(kind=8), allocatable, save :: virtual_density(:, :, :)
  complex(kind=8), allocatable, save :: reduced_density(:, :)
  complex(kind=8), allocatable, save :: degenerate_density(:, :, :)
  complex(kind=8), allocatable, save :: integrated_density(:, :, :)
  logical, save :: reduced_matrix_is_available = .false.
  logical, save :: degenerate_matrix_is_available = .false.
  logical, save :: integrated_matrix_is_available = .false.

  public :: reset_spin_density_fks_matrices
  public :: set_spin_density_fks_collection
  public :: spin_density_fks_collection_enabled
  public :: load_spin_density_born_matrix
  public :: load_spin_density_real_matrix
  public :: load_spin_density_color_matrix
  public :: load_spin_density_virtual_matrix
  public :: set_spin_density_virtual_matrix
  public :: reset_spin_density_virtual_matrix
  public :: reduce_spin_density_virtual_matrix
  public :: spin_density_active_open_size
  public :: spin_density_active_component_position
  public :: spin_density_born_amp_position
  public :: spin_density_nlo_amp_position
  public :: spin_density_born_matrix_available
  public :: spin_density_real_matrix_available
  public :: spin_density_color_matrix_available
  public :: spin_density_virtual_matrix_available
  public :: get_spin_density_born_matrix
  public :: get_spin_density_real_matrix
  public :: get_spin_density_color_matrix
  public :: get_spin_density_virtual_matrix
  public :: reset_spin_density_reduced_matrix
  public :: set_spin_density_reduced_from_real
  public :: set_spin_density_reduced_from_born
  public :: add_spin_density_reduced_color
  public :: get_spin_density_reduced_matrix
  public :: reset_spin_density_degenerate_matrix
  public :: set_spin_density_degenerate_from_born
  public :: get_spin_density_degenerate_matrix
  public :: spin_density_reduced_matrix_available
  public :: spin_density_degenerate_matrix_available
  public :: reset_spin_density_integrated_matrix
  public :: add_spin_density_integrated_born
  public :: add_spin_density_integrated_color
  public :: get_spin_density_integrated_matrix
  public :: spin_density_integrated_matrix_available

  interface
    integer function sdm_branch_max_open_size()
    end function sdm_branch_max_open_size

    integer function sdm_contribution_component_position(contribution)
      integer, intent(in) :: contribution
    end function sdm_contribution_component_position

    integer function sdm_fks_correlation_leg(configuration)
      integer, intent(in) :: configuration
    end function sdm_fks_correlation_leg

    integer function sdm_local_correlation_leg(contribution, global_leg)
      integer, intent(in) :: contribution, global_leg
    end function sdm_local_correlation_leg
  end interface

contains

  subroutine set_spin_density_fks_collection(enabled)
    logical, intent(in) :: enabled
    collection_enabled = enabled
    if (enabled) call ensure_spin_density_fks_matrices()
  end subroutine set_spin_density_fks_collection


  logical function spin_density_fks_collection_enabled()
    spin_density_fks_collection_enabled = collection_enabled
  end function spin_density_fks_collection_enabled

  subroutine ensure_spin_density_fks_matrices()
    if (initialized) return
    maximum_open_size = sdm_branch_max_open_size()
    if (maximum_open_size < 1) then
      call fail_spin_density_fks('the generated open size is invalid')
    end if
    allocate(born_density(2, maximum_open_size, maximum_open_size))
    allocate(real_density(2, maximum_open_size, maximum_open_size))
    allocate(color_density(maximum_open_size, maximum_open_size))
    allocate(virtual_density(3, maximum_open_size, maximum_open_size))
    allocate(reduced_density(maximum_open_size, maximum_open_size))
    allocate(degenerate_density(3, maximum_open_size, maximum_open_size))
    allocate(integrated_density(3, maximum_open_size, maximum_open_size))
    call identify_density_amp_positions()
    initialized = .true.
    call reset_spin_density_fks_matrices()
  end subroutine ensure_spin_density_fks_matrices


  subroutine reset_spin_density_fks_matrices()
    if (.not. initialized) then
      call ensure_spin_density_fks_matrices()
      return
    end if
    born_density = (0d0, 0d0)
    real_density = (0d0, 0d0)
    color_density = (0d0, 0d0)
    virtual_density = (0d0, 0d0)
    reduced_density = (0d0, 0d0)
    degenerate_density = (0d0, 0d0)
    integrated_density = (0d0, 0d0)
    active_open_size = 0
    active_component_position = 0
    born_matrix_is_available = .false.
    real_matrix_is_available = .false.
    color_matrix_is_available = .false.
    virtual_matrix_is_available = .false.
    reduced_matrix_is_available = .false.
    degenerate_matrix_is_available = .false.
    integrated_matrix_is_available = .false.
  end subroutine reset_spin_density_fks_matrices


  subroutine load_spin_density_born_matrix( &
       contribution, configuration, event_slot)
    integer, intent(in) :: contribution, configuration, event_slot
    integer :: global_leg, local_leg

    call ensure_spin_density_fks_matrices()
    born_density = (0d0, 0d0)
    born_matrix_is_available = .false.
    real_matrix_is_available = .false.
    global_leg = sdm_fks_correlation_leg(configuration)
    local_leg = sdm_local_correlation_leg(contribution, global_leg)
    call sdm_born_block_density( &
         contribution, event_slot, local_leg, active_open_size, &
         born_density)
    active_component_position = &
         sdm_contribution_component_position(contribution)
    call validate_loaded_density('Born')
    born_matrix_is_available = .true.
  end subroutine load_spin_density_born_matrix


  subroutine load_spin_density_real_matrix(configuration, event_slot, &
                                           contribution)
    integer, intent(in) :: configuration, event_slot, contribution

    call ensure_spin_density_fks_matrices()
    real_density = (0d0, 0d0)
    real_matrix_is_available = .false.
    call sdm_real_block_density( &
         configuration, event_slot, active_open_size, real_density)
    active_component_position = &
         sdm_contribution_component_position(contribution)
    call validate_loaded_density('real')
    real_matrix_is_available = .true.
  end subroutine load_spin_density_real_matrix


  subroutine load_spin_density_color_matrix( &
       contribution, first, second, event_slot)
    integer, intent(in) :: contribution, first, second, event_slot

    call ensure_spin_density_fks_matrices()
    color_density = (0d0, 0d0)
    color_matrix_is_available = .false.
    call sdm_color_block_density_pair( &
         contribution, first, second, event_slot, active_open_size, &
         color_density)
    active_component_position = &
         sdm_contribution_component_position(contribution)
    call validate_loaded_density('color')
    color_matrix_is_available = .true.
  end subroutine load_spin_density_color_matrix


  subroutine load_spin_density_virtual_matrix( &
       contribution, event_slot, precision_asked, precision, return_code)
    integer, intent(in) :: contribution, event_slot
    double precision, intent(in) :: precision_asked
    double precision, intent(out) :: precision
    integer, intent(out) :: return_code

    call ensure_spin_density_fks_matrices()
    virtual_density = (0d0, 0d0)
    virtual_matrix_is_available = .false.
    call sdm_virtual_block_density( &
         contribution, event_slot, precision_asked, active_open_size, &
         virtual_density, precision, return_code)
    active_component_position = &
         sdm_contribution_component_position(contribution)
    call validate_loaded_density('virtual')
    virtual_matrix_is_available = .true.
  end subroutine load_spin_density_virtual_matrix


  subroutine set_spin_density_virtual_matrix(contribution, density)
    integer, intent(in) :: contribution
    complex(kind=8), intent(in) :: density(:, :, :)

    call ensure_spin_density_fks_matrices()
    if (size(density, 1) /= 3 .or. &
        size(density, 2) /= size(density, 3) .or. &
        size(density, 2) < 1 .or. &
        size(density, 2) > maximum_open_size) then
      call fail_spin_density_fks( &
           'a supplied virtual density has the wrong shape')
    end if
    virtual_density = (0d0, 0d0)
    active_open_size = size(density, 2)
    active_component_position = &
         sdm_contribution_component_position(contribution)
    virtual_density(:, 1:active_open_size, 1:active_open_size) = density
    call validate_loaded_density('virtual')
    virtual_matrix_is_available = .true.
  end subroutine set_spin_density_virtual_matrix


  subroutine reset_spin_density_virtual_matrix()
    call ensure_spin_density_fks_matrices()
    virtual_density = (0d0, 0d0)
    virtual_matrix_is_available = .false.
  end subroutine reset_spin_density_virtual_matrix


  subroutine reduce_spin_density_virtual_matrix( &
       averaged_virtual, sampling_fraction)
    double precision, intent(in) :: averaged_virtual, sampling_fraction

    call ensure_spin_density_fks_matrices()
    if (.not. virtual_matrix_is_available) then
      call fail_spin_density_fks('no virtual density is available to reduce')
    end if
    if (.not. born_matrix_is_available) then
      call fail_spin_density_fks( &
           'no Born density is available to reduce a virtual density')
    end if
    if (sampling_fraction <= 0d0) then
      call fail_spin_density_fks( &
           'the virtual sampling fraction is not positive')
    end if
    virtual_density(1, 1:active_open_size, 1:active_open_size) = ( &
         virtual_density(1, 1:active_open_size, 1:active_open_size) - &
         averaged_virtual* &
         born_density(1, 1:active_open_size, 1:active_open_size))/ &
         sampling_fraction
  end subroutine reduce_spin_density_virtual_matrix


  integer function spin_density_active_open_size()
    call ensure_spin_density_fks_matrices()
    spin_density_active_open_size = active_open_size
  end function spin_density_active_open_size


  integer function spin_density_active_component_position()
    call ensure_spin_density_fks_matrices()
    spin_density_active_component_position = active_component_position
  end function spin_density_active_component_position


  integer function spin_density_born_amp_position()
    call ensure_spin_density_fks_matrices()
    spin_density_born_amp_position = born_amp_position
  end function spin_density_born_amp_position


  integer function spin_density_nlo_amp_position()
    call ensure_spin_density_fks_matrices()
    spin_density_nlo_amp_position = nlo_amp_position
  end function spin_density_nlo_amp_position


  logical function spin_density_born_matrix_available()
    spin_density_born_matrix_available = born_matrix_is_available
  end function spin_density_born_matrix_available


  logical function spin_density_real_matrix_available()
    spin_density_real_matrix_available = real_matrix_is_available
  end function spin_density_real_matrix_available


  logical function spin_density_color_matrix_available()
    spin_density_color_matrix_available = color_matrix_is_available
  end function spin_density_color_matrix_available


  logical function spin_density_virtual_matrix_available()
    spin_density_virtual_matrix_available = virtual_matrix_is_available
  end function spin_density_virtual_matrix_available


  logical function spin_density_reduced_matrix_available()
    spin_density_reduced_matrix_available = reduced_matrix_is_available
  end function spin_density_reduced_matrix_available


  logical function spin_density_degenerate_matrix_available()
    spin_density_degenerate_matrix_available = degenerate_matrix_is_available
  end function spin_density_degenerate_matrix_available


  logical function spin_density_integrated_matrix_available()
    spin_density_integrated_matrix_available = integrated_matrix_is_available
  end function spin_density_integrated_matrix_available


  subroutine reset_spin_density_reduced_matrix()
    call ensure_spin_density_fks_matrices()
    reduced_density = (0d0, 0d0)
    reduced_matrix_is_available = .false.
  end subroutine reset_spin_density_reduced_matrix


  subroutine set_spin_density_reduced_from_real(multiplier)
    double precision, intent(in) :: multiplier

    call ensure_spin_density_fks_matrices()
    if (.not. real_matrix_is_available) then
      call fail_spin_density_fks('no real density is available for reduction')
    end if
    reduced_density = (0d0, 0d0)
    reduced_density(1:active_open_size, 1:active_open_size) = &
         multiplier*real_density(1, 1:active_open_size, 1:active_open_size)
    reduced_matrix_is_available = .true.
  end subroutine set_spin_density_reduced_from_real


  subroutine set_spin_density_reduced_from_born(ordinary, correlated)
    complex(kind=8), intent(in) :: ordinary, correlated

    call ensure_spin_density_fks_matrices()
    if (.not. born_matrix_is_available) then
      call fail_spin_density_fks('no Born density is available for reduction')
    end if
    reduced_density = (0d0, 0d0)
    reduced_density(1:active_open_size, 1:active_open_size) = &
         ordinary*born_density(1, 1:active_open_size, 1:active_open_size) + &
         correlated*born_density(2, 1:active_open_size, 1:active_open_size)
    reduced_matrix_is_available = .true.
  end subroutine set_spin_density_reduced_from_born


  subroutine add_spin_density_reduced_color(multiplier)
    complex(kind=8), intent(in) :: multiplier

    call ensure_spin_density_fks_matrices()
    if (.not. color_matrix_is_available) then
      call fail_spin_density_fks('no color density is available for reduction')
    end if
    if (.not. reduced_matrix_is_available) then
      reduced_density = (0d0, 0d0)
    end if
    reduced_density(1:active_open_size, 1:active_open_size) = &
         reduced_density(1:active_open_size, 1:active_open_size) + &
         multiplier*color_density(1:active_open_size, 1:active_open_size)
    reduced_matrix_is_available = .true.
  end subroutine add_spin_density_reduced_color


  subroutine get_spin_density_reduced_matrix(density)
    complex(kind=8), intent(out) :: density(:, :)

    call ensure_spin_density_fks_matrices()
    if (.not. reduced_matrix_is_available) then
      call fail_spin_density_fks('no reduced density matrix is available')
    end if
    if (size(density, 1) /= active_open_size .or. &
        size(density, 2) /= active_open_size) then
      call fail_spin_density_fks('a reduced-density target has wrong shape')
    end if
    density = reduced_density(1:active_open_size, 1:active_open_size)
  end subroutine get_spin_density_reduced_matrix


  subroutine set_spin_density_degenerate_from_born(multipliers)
    complex(kind=8), intent(in) :: multipliers(3)
    integer :: coefficient

    call ensure_spin_density_fks_matrices()
    if (.not. born_matrix_is_available) then
      call fail_spin_density_fks( &
           'no Born density is available for a degenerate remainder')
    end if
    degenerate_density = (0d0, 0d0)
    do coefficient = 1, 3
      degenerate_density(coefficient, 1:active_open_size, &
                         1:active_open_size) = multipliers(coefficient)* &
           born_density(1, 1:active_open_size, 1:active_open_size)
    end do
    degenerate_matrix_is_available = .true.
  end subroutine set_spin_density_degenerate_from_born


  subroutine reset_spin_density_degenerate_matrix()
    call ensure_spin_density_fks_matrices()
    degenerate_density = (0d0, 0d0)
    degenerate_matrix_is_available = .false.
  end subroutine reset_spin_density_degenerate_matrix


  subroutine get_spin_density_degenerate_matrix(density)
    complex(kind=8), intent(out) :: density(:, :, :)

    call ensure_spin_density_fks_matrices()
    if (.not. degenerate_matrix_is_available) then
      call fail_spin_density_fks('no degenerate density matrix is available')
    end if
    if (size(density, 1) /= 3 .or. &
        size(density, 2) /= active_open_size .or. &
        size(density, 3) /= active_open_size) then
      call fail_spin_density_fks( &
           'a degenerate-density target has wrong shape')
    end if
    density = degenerate_density(:, 1:active_open_size, 1:active_open_size)
  end subroutine get_spin_density_degenerate_matrix


  subroutine reset_spin_density_integrated_matrix()
    call ensure_spin_density_fks_matrices()
    integrated_density = (0d0, 0d0)
    integrated_matrix_is_available = .false.
  end subroutine reset_spin_density_integrated_matrix


  subroutine add_spin_density_integrated_born(coefficients)
    complex(kind=8), intent(in) :: coefficients(3)
    integer :: coefficient

    call ensure_spin_density_fks_matrices()
    if (.not. born_matrix_is_available) then
      call fail_spin_density_fks( &
           'no Born density is available for an integrated term')
    end if
    do coefficient = 1, 3
      integrated_density(coefficient, 1:active_open_size, &
                         1:active_open_size) = &
           integrated_density(coefficient, 1:active_open_size, &
                              1:active_open_size) + &
           coefficients(coefficient)* &
           born_density(1, 1:active_open_size, 1:active_open_size)
    end do
    integrated_matrix_is_available = .true.
  end subroutine add_spin_density_integrated_born


  subroutine add_spin_density_integrated_color(coefficients)
    complex(kind=8), intent(in) :: coefficients(3)
    integer :: coefficient

    call ensure_spin_density_fks_matrices()
    if (.not. color_matrix_is_available) then
      call fail_spin_density_fks( &
           'no color density is available for an integrated term')
    end if
    do coefficient = 1, 3
      integrated_density(coefficient, 1:active_open_size, &
                         1:active_open_size) = &
           integrated_density(coefficient, 1:active_open_size, &
                              1:active_open_size) + &
           coefficients(coefficient)* &
           color_density(1:active_open_size, 1:active_open_size)
    end do
    integrated_matrix_is_available = .true.
  end subroutine add_spin_density_integrated_color


  subroutine get_spin_density_integrated_matrix(density)
    complex(kind=8), intent(out) :: density(:, :, :)

    call ensure_spin_density_fks_matrices()
    if (.not. integrated_matrix_is_available) then
      call fail_spin_density_fks('no integrated density matrix is available')
    end if
    if (size(density, 1) /= 3 .or. &
        size(density, 2) /= active_open_size .or. &
        size(density, 3) /= active_open_size) then
      call fail_spin_density_fks( &
           'an integrated-density target has wrong shape')
    end if
    density = integrated_density(:, 1:active_open_size, 1:active_open_size)
  end subroutine get_spin_density_integrated_matrix


  subroutine get_spin_density_born_matrix(density)
    complex(kind=8), intent(out) :: density(:, :, :)
    call copy_ranked_density(born_density, density, 2)
  end subroutine get_spin_density_born_matrix


  subroutine get_spin_density_real_matrix(density)
    complex(kind=8), intent(out) :: density(:, :, :)
    call copy_ranked_density(real_density, density, 2)
  end subroutine get_spin_density_real_matrix


  subroutine get_spin_density_color_matrix(density)
    complex(kind=8), intent(out) :: density(:, :)

    call ensure_spin_density_fks_matrices()
    if (active_open_size < 1) then
      call fail_spin_density_fks('no color density is loaded')
    end if
    if (size(density, 1) /= active_open_size .or. &
        size(density, 2) /= active_open_size) then
      call fail_spin_density_fks('a color-density target has the wrong shape')
    end if
    density = color_density(1:active_open_size, 1:active_open_size)
  end subroutine get_spin_density_color_matrix


  subroutine get_spin_density_virtual_matrix(density)
    complex(kind=8), intent(out) :: density(:, :, :)
    call copy_ranked_density(virtual_density, density, 3)
  end subroutine get_spin_density_virtual_matrix


  subroutine copy_ranked_density(source, target, rank_count)
    complex(kind=8), intent(in) :: source(:, :, :)
    complex(kind=8), intent(out) :: target(:, :, :)
    integer, intent(in) :: rank_count

    call ensure_spin_density_fks_matrices()
    if (active_open_size < 1) then
      call fail_spin_density_fks('no block density is loaded')
    end if
    if (size(target, 1) /= rank_count .or. &
        size(target, 2) /= active_open_size .or. &
        size(target, 3) /= active_open_size) then
      call fail_spin_density_fks('a block-density target has the wrong shape')
    end if
    target = source(1:rank_count, 1:active_open_size, 1:active_open_size)
  end subroutine copy_ranked_density


  subroutine validate_loaded_density(label)
    character(len=*), intent(in) :: label

    if (active_open_size < 1 .or. &
        active_open_size > maximum_open_size) then
      call fail_spin_density_fks( &
           trim(label)//' block density has an invalid open size')
    end if
    if (active_component_position < 1) then
      call fail_spin_density_fks( &
           trim(label)//' block density has an invalid component')
    end if
  end subroutine validate_loaded_density


  subroutine identify_density_amp_positions()
    integer :: candidate, position
    integer :: target_orders(nsplitorders)

    born_amp_position = 0
    nlo_amp_position = 0
    ! fNLO bundles carry one LO order and its QCD NLO successor.  Select the
    ! pair by their two-unit squared-order separation, independent of array
    ! ordering.  If only one slot exists it is both positions (LO-only use).
    do candidate = 1, amp_split_size
      target_orders = amp_split_orders(candidate, :)
      target_orders(qcd_pos) = target_orders(qcd_pos) + 2
      do position = 1, amp_split_size
        if (.not. all(amp_split_orders(position, :) == target_orders)) cycle
        if (born_amp_position /= 0) then
          call fail_spin_density_fks( &
               'more than one LO/NLO amplitude-order pair is present')
        end if
        born_amp_position = candidate
        nlo_amp_position = position
      end do
    end do
    if (born_amp_position == 0) then
      if (amp_split_size == 1) then
        born_amp_position = 1
        nlo_amp_position = 1
      else
        call fail_spin_density_fks( &
             'the unique LO/NLO amplitude-order pair is absent')
      end if
    end if
  end subroutine identify_density_amp_positions


  subroutine fail_spin_density_fks(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in spin_density_fks_matrices: '//trim(message)
    stop 1
  end subroutine fail_spin_density_fks

end module spin_density_fks_matrices

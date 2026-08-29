module multiplicative_point_workspace
  use process_dimensions, only: nexternal
  use multiplicative_process_plan, only: multiplicative_process_plan_type
  use multiplicative_phase_space, only: multiplicative_phase_space_assembly
  use multiplicative_density_terms, only: block_nlo_distribution, &
       multiplicative_density_tuple, density_tuple_schedule
  use multiplicative_density_contraction, only: multiplicative_density_basis
  use multiplicative_runtime, only: multiplicative_event_evaluation
  use multiplicative_reweighter, only: multiplicative_partonic_reweight
  use multiplicative_lambda_validation, only: &
       multiplicative_lambda_accumulator, &
       initialize_multiplicative_lambda_accumulator
  implicit none
  private

  type, public :: multiplicative_point_workspace_type
    type(multiplicative_phase_space_assembly) :: assembly
    type(block_nlo_distribution), allocatable :: distributions(:)
    type(multiplicative_density_tuple) :: tuple
    type(density_tuple_schedule) :: tuple_schedule
    type(multiplicative_density_basis) :: density_basis
    type(multiplicative_event_evaluation) :: evaluation
    type(multiplicative_partonic_reweight) :: variation_reweight
    type(multiplicative_lambda_accumulator) :: lambda_validation
    integer, allocatable :: sampled_fks(:)
    integer, allocatable :: sampled_integer(:)
    integer, allocatable :: sampled_dimension(:)
    double precision, allocatable :: sampled_volume(:)
    double precision, allocatable :: radiation_grid_weights(:)
    double precision, allocatable :: radiation_grid_groups(:, :)
    double precision, allocatable :: linear_radiation_groups(:, :)
    logical, allocatable :: component_owned(:)
    integer, allocatable :: factor_indices(:)
    double precision, allocatable :: plotted_weight(:)
    integer, allocatable :: validation_block_orders(:)
    integer, allocatable :: channel_event_slots(:)
    integer, allocatable :: previous_event_slots(:)
    integer, allocatable :: decay_block_factor_indices(:)
  end type multiplicative_point_workspace_type

  type(multiplicative_point_workspace_type), target, save :: point_workspace

  public :: acquire_multiplicative_point_workspace

contains

  subroutine acquire_multiplicative_point_workspace( &
       plan, plot_weight_count, decay_factor_count, workspace)
    type(multiplicative_process_plan_type), intent(in) :: plan
    integer, intent(in) :: plot_weight_count, decay_factor_count
    type(multiplicative_point_workspace_type), pointer, intent(out) :: workspace

    if (.not. plan%initialized) then
      call fail_point_workspace('an uninitialized process plan was supplied')
    end if
    if (plot_weight_count < 1) then
      call fail_point_workspace('the plot-weight capacity is not positive')
    end if
    if (decay_factor_count < 0 .or. &
        decay_factor_count > plan%decay_scale_factor_count) then
      call fail_point_workspace('the decay-factor capacity is invalid')
    end if
    call ensure_integer_vector( &
         point_workspace%sampled_fks, plan%contribution_count)
    call ensure_integer_vector( &
         point_workspace%sampled_integer, plan%contribution_count)
    call ensure_integer_vector( &
         point_workspace%sampled_dimension, plan%contribution_count)
    call ensure_real_vector( &
         point_workspace%sampled_volume, plan%contribution_count)
    call ensure_real_vector( &
         point_workspace%radiation_grid_weights, plan%contribution_count)
    call ensure_real_matrix( &
         point_workspace%radiation_grid_groups, &
         plan%contribution_count, 3)
    call ensure_real_matrix( &
         point_workspace%linear_radiation_groups, &
         plan%contribution_count, 3)
    call ensure_logical_vector( &
         point_workspace%component_owned, plan%component_count)
    call ensure_integer_vector( &
         point_workspace%validation_block_orders, plan%component_count)
    call ensure_integer_vector( &
         point_workspace%factor_indices, decay_factor_count)
    call ensure_real_vector( &
         point_workspace%plotted_weight, plot_weight_count)
    call ensure_block_vector(point_workspace%channel_event_slots)
    call ensure_block_vector(point_workspace%previous_event_slots)
    call ensure_block_vector(point_workspace%decay_block_factor_indices)
    if (allocated(point_workspace%distributions)) then
      if (size(point_workspace%distributions) /= plan%component_count) &
           deallocate(point_workspace%distributions)
    end if
    if (.not. allocated(point_workspace%distributions)) &
         allocate(point_workspace%distributions(plan%component_count))

    point_workspace%sampled_fks = 0
    point_workspace%sampled_integer = 0
    point_workspace%sampled_dimension = 0
    point_workspace%sampled_volume = 0d0
    point_workspace%radiation_grid_weights = 0d0
    point_workspace%radiation_grid_groups = 0d0
    point_workspace%linear_radiation_groups = 0d0
    point_workspace%component_owned = .false.
    point_workspace%factor_indices = 1
    point_workspace%plotted_weight = 0d0
    point_workspace%validation_block_orders = 0
    point_workspace%channel_event_slots = 0
    point_workspace%previous_event_slots = 0
    point_workspace%decay_block_factor_indices = 1
    call initialize_multiplicative_lambda_accumulator( &
         point_workspace%lambda_validation, plan%component_count)
    workspace => point_workspace
  end subroutine acquire_multiplicative_point_workspace


  subroutine ensure_integer_vector(values, count)
    integer, allocatable, intent(inout) :: values(:)
    integer, intent(in) :: count

    if (allocated(values)) then
      if (size(values) /= count) deallocate(values)
    end if
    if (.not. allocated(values)) allocate(values(count))
  end subroutine ensure_integer_vector


  subroutine ensure_block_vector(values)
    integer, allocatable, intent(inout) :: values(:)

    if (allocated(values)) then
      if (lbound(values, 1) /= 0 .or. &
          ubound(values, 1) /= nexternal) deallocate(values)
    end if
    if (.not. allocated(values)) allocate(values(0:nexternal))
  end subroutine ensure_block_vector


  subroutine ensure_logical_vector(values, count)
    logical, allocatable, intent(inout) :: values(:)
    integer, intent(in) :: count

    if (allocated(values)) then
      if (size(values) /= count) deallocate(values)
    end if
    if (.not. allocated(values)) allocate(values(count))
  end subroutine ensure_logical_vector


  subroutine ensure_real_vector(values, count)
    double precision, allocatable, intent(inout) :: values(:)
    integer, intent(in) :: count

    if (allocated(values)) then
      if (size(values) /= count) deallocate(values)
    end if
    if (.not. allocated(values)) allocate(values(count))
  end subroutine ensure_real_vector


  subroutine ensure_real_matrix(values, rows, columns)
    double precision, allocatable, intent(inout) :: values(:, :)
    integer, intent(in) :: rows, columns

    if (allocated(values)) then
      if (size(values, 1) /= rows .or. &
          size(values, 2) /= columns) deallocate(values)
    end if
    if (.not. allocated(values)) allocate(values(rows, columns))
  end subroutine ensure_real_matrix


  subroutine fail_point_workspace(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in multiplicative_point_workspace: '//trim(message)
    stop 1
  end subroutine fail_point_workspace

end module multiplicative_point_workspace

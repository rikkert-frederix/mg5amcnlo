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
    double precision, allocatable :: family_single_weights(:)
    integer, allocatable :: channel_event_slots(:)
    integer, allocatable :: previous_event_slots(:)
    integer, allocatable :: decay_block_factor_indices(:)
    ! Accepted exact event families are staged until every partonic density
    ! and scale polynomial has been evaluated.  PDF members can then be
    ! initialized outside the family loop and each member is reused for the
    ! complete batch before analyses are called.
    integer :: family_count = 0
    integer :: family_capacity = 0
    integer :: family_event_capacity = 0
    integer :: family_plot_weight_count = 0
    double precision, allocatable :: family_momenta(:, :, :)
    double precision, allocatable :: family_y_to_lab(:)
    double precision, allocatable :: family_bjorken_x(:, :)
    double precision, allocatable :: family_mu2_f(:)
    double precision, allocatable :: family_pdf_partonic_factor(:)
    double precision, allocatable :: family_pdf_luminosity(:)
    double precision, allocatable :: family_scale_mu2_f(:, :)
    double precision, allocatable :: family_scale_luminosity(:, :)
    double precision, allocatable :: family_plot_weights(:, :)
    integer, allocatable :: family_pdgs(:, :)
    integer, allocatable :: family_origin_blocks(:, :)
    integer, allocatable :: family_luminosity_configuration(:)
    integer, allocatable :: family_production_event_slot(:)
    integer, allocatable :: family_luminosity_owner(:)
  end type multiplicative_point_workspace_type

  type(multiplicative_point_workspace_type), target, save :: point_workspace

  public :: acquire_multiplicative_point_workspace
  public :: prepare_multiplicative_family_batch
  public :: stage_multiplicative_family

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
    call ensure_real_vector( &
         point_workspace%family_single_weights, plan%component_count)
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
    point_workspace%family_single_weights = 0d0
    point_workspace%channel_event_slots = 0
    point_workspace%previous_event_slots = 0
    point_workspace%decay_block_factor_indices = 1
    call initialize_multiplicative_lambda_accumulator( &
         point_workspace%lambda_validation, plan%component_count)
    workspace => point_workspace
  end subroutine acquire_multiplicative_point_workspace


  subroutine prepare_multiplicative_family_batch( &
       workspace, family_capacity, event_capacity, plot_weight_count)
    type(multiplicative_point_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: family_capacity, event_capacity
    integer, intent(in) :: plot_weight_count
    logical :: resize

    if (family_capacity < 1 .or. event_capacity < 1 .or. &
        plot_weight_count < 1) then
      call fail_point_workspace('an event-family batch has invalid capacity')
    end if
    resize = workspace%family_capacity /= family_capacity .or. &
         workspace%family_event_capacity /= event_capacity .or. &
         workspace%family_plot_weight_count /= plot_weight_count
    if (resize) then
      if (allocated(workspace%family_momenta)) &
           deallocate(workspace%family_momenta)
      if (allocated(workspace%family_y_to_lab)) &
           deallocate(workspace%family_y_to_lab)
      if (allocated(workspace%family_bjorken_x)) &
           deallocate(workspace%family_bjorken_x)
      if (allocated(workspace%family_mu2_f)) &
           deallocate(workspace%family_mu2_f)
      if (allocated(workspace%family_pdf_partonic_factor)) &
           deallocate(workspace%family_pdf_partonic_factor)
      if (allocated(workspace%family_pdf_luminosity)) &
           deallocate(workspace%family_pdf_luminosity)
      if (allocated(workspace%family_scale_mu2_f)) &
           deallocate(workspace%family_scale_mu2_f)
      if (allocated(workspace%family_scale_luminosity)) &
           deallocate(workspace%family_scale_luminosity)
      if (allocated(workspace%family_plot_weights)) &
           deallocate(workspace%family_plot_weights)
      if (allocated(workspace%family_pdgs)) &
           deallocate(workspace%family_pdgs)
      if (allocated(workspace%family_origin_blocks)) &
           deallocate(workspace%family_origin_blocks)
      if (allocated(workspace%family_luminosity_configuration)) &
           deallocate(workspace%family_luminosity_configuration)
      if (allocated(workspace%family_production_event_slot)) &
           deallocate(workspace%family_production_event_slot)
      if (allocated(workspace%family_luminosity_owner)) &
           deallocate(workspace%family_luminosity_owner)
      allocate(workspace%family_momenta( &
           0:3, event_capacity, family_capacity))
      allocate(workspace%family_y_to_lab(family_capacity))
      allocate(workspace%family_bjorken_x(2, family_capacity))
      allocate(workspace%family_mu2_f(family_capacity))
      allocate(workspace%family_pdf_partonic_factor(family_capacity))
      allocate(workspace%family_pdf_luminosity(family_capacity))
      allocate(workspace%family_scale_mu2_f( &
           plot_weight_count, family_capacity))
      allocate(workspace%family_scale_luminosity( &
           plot_weight_count, family_capacity))
      allocate(workspace%family_plot_weights( &
           plot_weight_count, family_capacity))
      allocate(workspace%family_pdgs(event_capacity, family_capacity))
      allocate(workspace%family_origin_blocks( &
           event_capacity, family_capacity))
      allocate(workspace%family_luminosity_configuration(family_capacity))
      allocate(workspace%family_production_event_slot(family_capacity))
      allocate(workspace%family_luminosity_owner(family_capacity))
      workspace%family_capacity = family_capacity
      workspace%family_event_capacity = event_capacity
      workspace%family_plot_weight_count = plot_weight_count
    end if
    workspace%family_count = 0
  end subroutine prepare_multiplicative_family_batch


  subroutine stage_multiplicative_family( &
       workspace, evaluation, production_mu2_f, &
       luminosity_configuration, pdf_partonic_factor, family_index)
    type(multiplicative_point_workspace_type), intent(inout) :: workspace
    type(multiplicative_event_evaluation), intent(in) :: evaluation
    double precision, intent(in) :: production_mu2_f, pdf_partonic_factor
    integer, intent(in) :: luminosity_configuration
    integer, intent(out) :: family_index
    integer :: previous_family

    if (.not. evaluation%available) then
      call fail_point_workspace('an unavailable event family was staged')
    end if
    if (workspace%family_capacity < 1 .or. &
        .not. allocated(workspace%family_momenta)) then
      call fail_point_workspace('the event-family batch is not prepared')
    end if
    if (.not. allocated(evaluation%momenta) .or. &
        .not. allocated(evaluation%pdgs) .or. &
        .not. allocated(evaluation%origin_blocks) .or. &
        .not. allocated(evaluation%event_slots)) then
      call fail_point_workspace('an event family is not allocated')
    end if
    if (lbound(evaluation%momenta, 1) /= 0 .or. &
        ubound(evaluation%momenta, 1) /= 3 .or. &
        size(evaluation%momenta, 2) /= workspace%family_event_capacity .or. &
        size(evaluation%pdgs) /= workspace%family_event_capacity .or. &
        size(evaluation%origin_blocks) /= workspace%family_event_capacity .or. &
        lbound(evaluation%event_slots, 1) /= 0 .or. &
        ubound(evaluation%event_slots, 1) /= nexternal) then
      call fail_point_workspace('an event family has the wrong shape')
    end if
    if (production_mu2_f <= 0d0 .or. luminosity_configuration < 1) then
      call fail_point_workspace('an event-family luminosity key is invalid')
    end if

    workspace%family_count = workspace%family_count + 1
    family_index = workspace%family_count
    if (family_index > workspace%family_capacity) then
      call fail_point_workspace('the event-family batch overflowed')
    end if
    workspace%family_momenta(:, :, family_index) = evaluation%momenta
    workspace%family_y_to_lab(family_index) = evaluation%y_to_lab
    workspace%family_bjorken_x(:, family_index) = evaluation%bjorken_x
    workspace%family_mu2_f(family_index) = production_mu2_f
    workspace%family_pdf_partonic_factor(family_index) = &
         pdf_partonic_factor
    workspace%family_pdf_luminosity(family_index) = 0d0
    workspace%family_scale_mu2_f(:, family_index) = 0d0
    workspace%family_scale_luminosity(:, family_index) = 0d0
    workspace%family_plot_weights(:, family_index) = 0d0
    workspace%family_pdgs(:, family_index) = evaluation%pdgs
    workspace%family_origin_blocks(:, family_index) = &
         evaluation%origin_blocks
    workspace%family_luminosity_configuration(family_index) = &
         luminosity_configuration
    workspace%family_production_event_slot(family_index) = &
         evaluation%event_slots(0)
    workspace%family_luminosity_owner(family_index) = family_index

    ! Store a stable owner for the exact production-luminosity key.  Decay
    ! families commonly differ in visible momenta while sharing all incoming
    ! data, so their scale and PDF variations can reuse one provider result.
    do previous_family = family_index - 1, 1, -1
      if (workspace%family_production_event_slot(previous_family) /= &
          workspace%family_production_event_slot(family_index)) cycle
      if (workspace%family_luminosity_configuration(previous_family) /= &
          luminosity_configuration) cycle
      if (workspace%family_mu2_f(previous_family) /= production_mu2_f) cycle
      if (any(workspace%family_bjorken_x(:, previous_family) /= &
          evaluation%bjorken_x)) cycle
      workspace%family_luminosity_owner(family_index) = &
           workspace%family_luminosity_owner(previous_family)
      exit
    end do
  end subroutine stage_multiplicative_family


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

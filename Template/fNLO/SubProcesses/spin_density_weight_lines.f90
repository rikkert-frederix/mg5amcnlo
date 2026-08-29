module spin_density_weight_lines
  use spin_density_matrix_results, only: spin_density_bornlike_branch, &
       spin_density_real_branch
  use multiplicative_nlo_decay, only: multiplicative_nlo_workspace, &
       add_multiplicative_block_density
  use decay_chain_parameters, only: decay_renormalization_scale, &
       decay_scale_species_index, decay_multiplicative_width_rescaling, &
       multiplicative_nlo_enabled
  use alfas_functions_module, only: alphas
  implicit none
  private

  logical, save :: initialized = .false.
  integer, save :: maximum_open_size = 0
  integer, save :: line_capacity = 0
  integer, save :: weight_capacity = 0
  logical, allocatable, save :: line_present(:)
  logical, allocatable, save :: weight_present(:, :)
  integer, allocatable, save :: line_component(:)
  integer, allocatable, save :: line_branch(:)
  integer, allocatable, save :: line_open_size(:)
  integer, allocatable, save :: line_qcd_power(:)
  integer, allocatable, save :: line_scale_pdg(:)
  logical, allocatable, save :: line_is_production(:)
  complex(kind=8), allocatable, save :: line_coefficients(:, :, :, :)
  complex(kind=8), allocatable, save :: evaluated_density(:, :, :, :)

  public :: clear_spin_density_weight_lines
  public :: record_spin_density_weight_line
  public :: copy_spin_density_weight_line
  public :: spin_density_weight_line_present
  public :: spin_density_weight_line_is_production
  public :: evaluate_spin_density_weight_line
  public :: spin_density_weight_line_multiplier
  public :: aggregate_spin_density_weight_lines

  interface
    integer function sdm_branch_max_open_size()
    end function sdm_branch_max_open_size
  end interface

contains

  subroutine initialize_spin_density_weight_lines()
    if (initialized) return
    maximum_open_size = sdm_branch_max_open_size()
    if (maximum_open_size < 1) then
      call fail_density_lines('the generated open size is invalid')
    end if
    initialized = .true.
  end subroutine initialize_spin_density_weight_lines


  subroutine clear_spin_density_weight_lines()
    ! The number and shape of generated lines are process-fixed.  Retain the
    ! backing arrays between integrand points and clear only their presence
    ! masks; each recorded line overwrites its active coefficient/density
    ! slices before aggregation.
    if (allocated(line_present)) line_present = .false.
    if (allocated(weight_present)) weight_present = .false.
  end subroutine clear_spin_density_weight_lines


  subroutine record_spin_density_weight_line( &
       line, component, branch, coefficients, qcd_power, is_production, &
       scale_pdg)
    integer, intent(in) :: line, component, branch, qcd_power, scale_pdg
    logical, intent(in) :: is_production
    complex(kind=8), intent(in) :: coefficients(:, :, :)
    integer :: open_size

    call initialize_spin_density_weight_lines()
    if (line < 1 .or. component < 1 .or. qcd_power < 0) then
      call fail_density_lines('a density weight-line identity is invalid')
    end if
    if (.not. is_production .and. scale_pdg == 0 .and. qcd_power > 0) then
      call fail_density_lines( &
           'a QCD decay density line has no scale species')
    end if
    if (branch /= spin_density_bornlike_branch .and. &
        branch /= spin_density_real_branch) then
      call fail_density_lines('a density weight-line branch is invalid')
    end if
    if (size(coefficients, 1) /= 3 .or. &
        size(coefficients, 2) /= size(coefficients, 3)) then
      call fail_density_lines('density weight coefficients have wrong shape')
    end if
    open_size = size(coefficients, 2)
    if (open_size < 1 .or. open_size > maximum_open_size) then
      call fail_density_lines('a density weight-line open size is invalid')
    end if
    call ensure_density_line_capacity(line, max(1, weight_capacity))
    if (line_present(line)) then
      call fail_density_lines('a density weight line was recorded twice')
    end if
    line_present(line) = .true.
    line_component(line) = component
    line_branch(line) = branch
    line_open_size(line) = open_size
    line_qcd_power(line) = qcd_power
    line_scale_pdg(line) = abs(scale_pdg)
    line_is_production(line) = is_production
    line_coefficients(:, :, :, line) = (0d0, 0d0)
    line_coefficients(:, 1:open_size, 1:open_size, line) = coefficients
    weight_present(:, line) = .false.
    evaluated_density(:, :, :, line) = (0d0, 0d0)
  end subroutine record_spin_density_weight_line


  subroutine copy_spin_density_weight_line(source_line, target_line)
    integer, intent(in) :: source_line, target_line

    call initialize_spin_density_weight_lines()
    if (source_line < 1 .or. source_line > line_capacity .or. &
        .not. line_present(source_line)) then
      call fail_density_lines('a source density weight line is absent')
    end if
    call ensure_density_line_capacity(target_line, max(1, weight_capacity))
    if (line_present(target_line)) then
      call fail_density_lines('a target density weight line already exists')
    end if
    line_present(target_line) = .true.
    line_component(target_line) = line_component(source_line)
    line_branch(target_line) = line_branch(source_line)
    line_open_size(target_line) = line_open_size(source_line)
    line_qcd_power(target_line) = line_qcd_power(source_line)
    line_scale_pdg(target_line) = line_scale_pdg(source_line)
    line_is_production(target_line) = line_is_production(source_line)
    line_coefficients(:, :, :, target_line) = &
         line_coefficients(:, :, :, source_line)
    weight_present(:, target_line) = weight_present(:, source_line)
    evaluated_density(:, :, :, target_line) = &
         evaluated_density(:, :, :, source_line)
  end subroutine copy_spin_density_weight_line


  logical function spin_density_weight_line_present(line)
    integer, intent(in) :: line

    spin_density_weight_line_present = initialized .and. line >= 1 .and. &
         line <= line_capacity
    if (spin_density_weight_line_present) then
      spin_density_weight_line_present = line_present(line)
    end if
  end function spin_density_weight_line_present


  logical function spin_density_weight_line_is_production(line)
    integer, intent(in) :: line

    if (.not. spin_density_weight_line_present(line)) then
      call fail_density_lines('a queried density weight line is absent')
    end if
    spin_density_weight_line_is_production = line_is_production(line)
  end function spin_density_weight_line_is_production


  double precision function spin_density_weight_line_multiplier( &
       line, production_g, factor_indices, parton_luminosity)
    integer, intent(in) :: line
    double precision, intent(in) :: production_g, parton_luminosity
    integer, intent(in) :: factor_indices(:)
    integer :: factor_index, species_index
    double precision :: local_g, scale
    double precision, parameter :: pi = 3.14159265358979323846d0

    if (.not. spin_density_weight_line_present(line)) then
      call fail_density_lines('a multiplied density weight line is absent')
    end if
    if (line_is_production(line)) then
      spin_density_weight_line_multiplier = parton_luminosity* &
           production_g**line_qcd_power(line)
      if (multiplicative_nlo_enabled()) then
        spin_density_weight_line_multiplier = &
             spin_density_weight_line_multiplier* &
             decay_multiplicative_width_rescaling(factor_indices)
      end if
      return
    end if
    local_g = 1d0
    if (line_qcd_power(line) > 0) then
      factor_index = 1
      if (size(factor_indices) > 0) then
        species_index = decay_scale_species_index(line_scale_pdg(line))
        if (species_index < 1 .or. &
            species_index > size(factor_indices)) then
          call fail_density_lines( &
               'a decay density line has no scale-factor index')
        end if
        factor_index = factor_indices(species_index)
      end if
      scale = decay_renormalization_scale( &
           line_scale_pdg(line), factor_index)
      local_g = sqrt(4d0*pi*alphas(scale))
    end if
    ! PDFs, incoming flux and the global integration weight belong solely to
    ! the production block and must not be repeated by a decay density.
    spin_density_weight_line_multiplier = &
         local_g**line_qcd_power(line)
  end function spin_density_weight_line_multiplier


  subroutine evaluate_spin_density_weight_line( &
       line, weight, logarithmic_mu2_r, mu2_f, reference_scale2, &
       multiplier)
    integer, intent(in) :: line, weight
    double precision, intent(in) :: logarithmic_mu2_r, mu2_f
    double precision, intent(in) :: reference_scale2, multiplier
    integer :: open_size
    complex(kind=8) :: logarithmic_density(maximum_open_size, &
                                           maximum_open_size)

    call initialize_spin_density_weight_lines()
    if (line < 1 .or. line > line_capacity .or. &
        .not. line_present(line)) then
      call fail_density_lines('an evaluated density weight line is absent')
    end if
    if (weight < 1 .or. reference_scale2 <= 0d0) then
      call fail_density_lines('an evaluated density weight is invalid')
    end if
    call ensure_density_line_capacity(line, weight)
    open_size = line_open_size(line)
    logarithmic_density = (0d0, 0d0)
    logarithmic_density(1:open_size, 1:open_size) = &
         line_coefficients(1, 1:open_size, 1:open_size, line) + &
         line_coefficients(2, 1:open_size, 1:open_size, line)* &
         log(logarithmic_mu2_r/reference_scale2) + &
         line_coefficients(3, 1:open_size, 1:open_size, line)* &
         log(mu2_f/reference_scale2)
    evaluated_density(weight, :, :, line) = &
         multiplier*logarithmic_density
    weight_present(weight, line) = .true.
  end subroutine evaluate_spin_density_weight_line


  subroutine aggregate_spin_density_weight_lines(workspace)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer :: line, open_size

    if (.not. initialized .or. line_capacity == 0) return
    do line = 1, line_capacity
      if (.not. line_present(line)) cycle
      if (line_component(line) > workspace%component_count) then
        call fail_density_lines( &
             'a density line component is absent from the workspace')
      end if
      if (workspace%weight_count > weight_capacity .or. &
          .not. all(weight_present(1:workspace%weight_count, line))) then
        call fail_density_lines( &
             'a density line has incomplete reweighting information')
      end if
      open_size = line_open_size(line)
      if (open_size /= &
          workspace%component_open_sizes(line_component(line))) then
        call fail_density_lines( &
             'a density line and component have different open sizes')
      end if
      call add_multiplicative_block_density( &
           workspace, line_component(line), line_branch(line), &
           evaluated_density(1:workspace%weight_count, &
                             1:open_size, 1:open_size, line))
    end do
  end subroutine aggregate_spin_density_weight_lines


  subroutine ensure_density_line_capacity(required_lines, required_weights)
    integer, intent(in) :: required_lines, required_weights
    integer :: old_lines, old_weights, new_lines, new_weights
    logical, allocatable :: old_line_present(:), old_weight_present(:, :)
    integer, allocatable :: old_component(:), old_branch(:), old_open(:)
    integer, allocatable :: old_qcd_power(:), old_scale_pdg(:)
    logical, allocatable :: old_is_production(:)
    complex(kind=8), allocatable :: old_coefficients(:, :, :, :)
    complex(kind=8), allocatable :: old_evaluated(:, :, :, :)

    call initialize_spin_density_weight_lines()
    if (required_lines <= line_capacity .and. &
        required_weights <= weight_capacity) return
    old_lines = line_capacity
    old_weights = weight_capacity
    new_lines = max(required_lines, max(1, 2*old_lines))
    new_weights = max(required_weights, max(1, 2*old_weights))
    if (old_lines > 0) then
      allocate(old_line_present(old_lines))
      allocate(old_weight_present(old_weights, old_lines))
      allocate(old_component(old_lines), old_branch(old_lines), &
               old_open(old_lines))
      allocate(old_qcd_power(old_lines), old_scale_pdg(old_lines))
      allocate(old_is_production(old_lines))
      allocate(old_coefficients(3, maximum_open_size, maximum_open_size, &
                                old_lines))
      allocate(old_evaluated(old_weights, maximum_open_size, &
                             maximum_open_size, old_lines))
      old_line_present = line_present
      old_weight_present = weight_present
      old_component = line_component
      old_branch = line_branch
      old_open = line_open_size
      old_qcd_power = line_qcd_power
      old_scale_pdg = line_scale_pdg
      old_is_production = line_is_production
      old_coefficients = line_coefficients
      old_evaluated = evaluated_density
      deallocate(line_present, weight_present, line_component, line_branch, &
                 line_open_size, line_qcd_power, line_scale_pdg, &
                 line_is_production, line_coefficients, evaluated_density)
    end if
    allocate(line_present(new_lines))
    allocate(weight_present(new_weights, new_lines))
    allocate(line_component(new_lines), line_branch(new_lines), &
             line_open_size(new_lines))
    allocate(line_qcd_power(new_lines), line_scale_pdg(new_lines))
    allocate(line_is_production(new_lines))
    allocate(line_coefficients(3, maximum_open_size, maximum_open_size, &
                               new_lines))
    allocate(evaluated_density(new_weights, maximum_open_size, &
                               maximum_open_size, new_lines))
    line_present = .false.
    weight_present = .false.
    line_component = 0
    line_branch = -1
    line_open_size = 0
    line_qcd_power = 0
    line_scale_pdg = 0
    line_is_production = .false.
    line_coefficients = (0d0, 0d0)
    evaluated_density = (0d0, 0d0)
    if (old_lines > 0) then
      line_present(1:old_lines) = old_line_present
      weight_present(1:old_weights, 1:old_lines) = old_weight_present
      line_component(1:old_lines) = old_component
      line_branch(1:old_lines) = old_branch
      line_open_size(1:old_lines) = old_open
      line_qcd_power(1:old_lines) = old_qcd_power
      line_scale_pdg(1:old_lines) = old_scale_pdg
      line_is_production(1:old_lines) = old_is_production
      line_coefficients(:, :, :, 1:old_lines) = old_coefficients
      evaluated_density(1:old_weights, :, :, 1:old_lines) = old_evaluated
      deallocate(old_line_present, old_weight_present, old_component, &
                 old_branch, old_open, old_qcd_power, old_scale_pdg, &
                 old_is_production, old_coefficients, old_evaluated)
    end if
    line_capacity = new_lines
    weight_capacity = new_weights
  end subroutine ensure_density_line_capacity


  subroutine fail_density_lines(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in spin_density_weight_lines: '//trim(message)
    stop 1
  end subroutine fail_density_lines

end module spin_density_weight_lines

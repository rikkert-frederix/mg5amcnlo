module check_poles_module
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
       validate_process_dimensions
  use fks_metadata, only: fks_i_d, pdg_type_d, validate_fks_metadata
  use run_state, only: lpp, ebeam
  use mint_module, only: iconfig, ichan, iconfigs
  use FKSParams, only: paramFileName, IRPoleCheckThreshold, &
       FKSParamReader
  implicit none
  private

  double precision, parameter :: pi = 3.1415926535897932385d0
  double precision, parameter :: zero = 0d0
  double precision, parameter :: rambo_accuracy = 1d-14
  integer, parameter :: rambo_iteration_limit = 10
  integer, parameter :: rambo_max_particles = 100

  double precision, allocatable :: generated_masses(:)
  logical :: generated_data_initialized = .false.

  double precision, allocatable :: virtual_weights(:, :)
  double precision, allocatable :: pole_accuracies(:)
  logical, allocatable :: kept_orders(:)
  logical :: pole_work_initialized = .false.

  double precision, allocatable :: rambo_log_weights(:)
  integer :: rambo_warnings(5) = 0
  double precision :: rambo_two_pi = 0d0
  double precision :: rambo_log_pi = 0d0
  logical :: rambo_initialized = .false.

  logical :: calculatedborn
  common /ccalculatedborn/ calculatedborn

  integer :: nfksprocess
  common /c_nfksprocess/ nfksprocess

  double precision :: qes2
  common /coupl_es/ qes2

  logical :: force_polecheck, polecheck_passed
  common /to_polecheck/ force_polecheck, polecheck_passed

  integer :: ret_code_ml
  common /to_ret_code/ ret_code_ml

  public :: run_check_poles
  public :: initialize_check_poles_data
  public :: finalize_check_poles_data
  public :: rambo_impl
  public :: rans_impl

  interface
    subroutine init_process_dimensions_bridge()
    end subroutine init_process_dimensions_bridge

    subroutine init_born_dimensions_bridge()
    end subroutine init_born_dimensions_bridge

    subroutine init_fks_metadata_bridge()
    end subroutine init_fks_metadata_bridge

    subroutine init_check_poles_data_bridge()
    end subroutine init_check_poles_data_bridge

    subroutine check_poles_set_model_scale(scale_value)
      double precision, intent(in) :: scale_value
    end subroutine check_poles_set_model_scale

    subroutine setrun_model_strong_coupling(value)
      double precision, intent(out) :: value
    end subroutine setrun_model_strong_coupling

    subroutine get_nsqso_loop(number_of_orders)
      integer, intent(out) :: number_of_orders
    end subroutine get_nsqso_loop

    subroutine get_answer_dimension(answer_dimension)
      integer, intent(out) :: answer_dimension
    end subroutine get_answer_dimension

    subroutine setrun()
    end subroutine setrun

    subroutine setpara(parameter_card)
      character(len=*), intent(in) :: parameter_card
    end subroutine setpara

    subroutine setcuts()
    end subroutine setcuts

    subroutine sync_cuts_bridge_state()
    end subroutine sync_cuts_bridge_state

    subroutine printout()
    end subroutine printout

    subroutine run_printout()
    end subroutine run_printout

    subroutine fks_inc_chooser()
    end subroutine fks_inc_chooser

    subroutine leshouche_inc_chooser()
    end subroutine leshouche_inc_chooser

    subroutine setfksfactor()
    end subroutine setfksfactor

    subroutine force_stability_check(enabled)
      logical, intent(in) :: enabled
    end subroutine force_stability_check

    subroutine collier_compute_uv_poles(enabled)
      logical, intent(in) :: enabled
    end subroutine collier_compute_uv_poles

    subroutine collier_compute_ir_poles(enabled)
      logical, intent(in) :: enabled
    end subroutine collier_compute_ir_poles

    subroutine update_as_param()
    end subroutine update_as_param

    subroutine sborn(momentum, born_weight)
      double precision, intent(in) :: momentum(0:3, *)
      double precision, intent(out) :: born_weight
    end subroutine sborn

    subroutine binothlha(momentum, born_weight, virtual_weight)
      double precision, intent(in) :: momentum(0:3, *)
      double precision, intent(inout) :: born_weight
      double precision, intent(out) :: virtual_weight
    end subroutine binothlha

    double precision function ran2()
    end function ran2
  end interface

contains

  subroutine initialize_check_poles_data(masses)
    implicit none
    double precision, intent(in) :: masses(:)

    call validate_process_dimensions()
    if (size(masses) /= nexternal) then
      call fail_check_poles('generated mass array has the wrong size')
    end if

    if (generated_data_initialized) then
      if (any(abs(generated_masses - masses) > 0d0)) then
        call fail_check_poles('generated masses changed after initialization')
      end if
      return
    end if

    allocate(generated_masses(nexternal))
    generated_masses = masses
    generated_data_initialized = .true.
  end subroutine initialize_check_poles_data


  subroutine initialize_pole_work()
    implicit none
    integer :: number_of_orders, answer_dimension

    if (pole_work_initialized) return
    call get_nsqso_loop(number_of_orders)
    call get_answer_dimension(answer_dimension)
    allocate(virtual_weights(0:3, 0:answer_dimension))
    allocate(pole_accuracies(0:number_of_orders))
    allocate(kept_orders(number_of_orders))
    pole_work_initialized = .true.
  end subroutine initialize_pole_work


  subroutine finalize_check_poles_data()
    implicit none

    if (allocated(generated_masses)) deallocate(generated_masses)
    if (allocated(virtual_weights)) deallocate(virtual_weights)
    if (allocated(pole_accuracies)) deallocate(pole_accuracies)
    if (allocated(kept_orders)) deallocate(kept_orders)
    if (allocated(rambo_log_weights)) deallocate(rambo_log_weights)
    generated_data_initialized = .false.
    pole_work_initialized = .false.
    rambo_initialized = .false.
    rambo_warnings = 0
    rambo_two_pi = 0d0
    rambo_log_pi = 0d0
  end subroutine finalize_check_poles_data


  subroutine run_check_poles(p_born)
    implicit none
    double precision, intent(inout) :: p_born(0:, :)
    double precision, allocatable :: momentum(:, :), rambo_masses(:)
    double precision, allocatable :: rambo_momenta(:, :)
    double precision :: tolerance, tolerance_default
    double precision :: renormalization_scale, energy
    double precision :: born_weight, virtual_weight, total_mass
    double precision :: strong_coupling
    integer :: npoints, npoints_checked, nfail
    integer :: particle, component

    call init_process_dimensions_bridge()
    call init_born_dimensions_bridge()
    call init_fks_metadata_bridge()
    call validate_process_dimensions(require_born=.true.)
    call validate_fks_metadata()

    force_polecheck = .true.
    call initialize_pole_work()

    call setrun()
    call setpara('param_card.dat')
    call setcuts()
    call sync_cuts_bridge_state()
    call printout()
    call run_printout()
    call init_check_poles_data_bridge()
    if (.not. generated_data_initialized) then
      call fail_check_poles('generated masses are not initialized')
    end if

    call FKSParamReader(paramFileName, .true., .false.)
    tolerance_default = IRPoleCheckThreshold
    iconfig = 1
    ichan = 1
    iconfigs(1) = iconfig

    total_mass = sum(generated_masses)
    if (abs(lpp(1)) == 1 .and. abs(lpp(2)) == 1) then
      energy = max((ebeam(1) + ebeam(2)) / 20d0, 2d0 * total_mass)
    else if (lpp(1) == 2 .and. lpp(2) == 2) then
      energy = max((ebeam(1) + ebeam(2)) / 200d0, &
           2d0 * total_mass)
    else
      energy = ebeam(1) + ebeam(2)
    end if
    renormalization_scale = energy / 2d0

    write (*, *) ' Insert the number of points to test'
    read (*, *) npoints
    write (*, *) 'Insert the relative tolerance'
    write (*, *) ' A negative number will mean use the default one: ', &
         tolerance_default
    read (*, *) tolerance
    if (tolerance <= zero) then
      tolerance = tolerance_default
    else
      IRPoleCheckThreshold = tolerance
    end if

    call check_poles_set_model_scale(renormalization_scale)
    qes2 = renormalization_scale**2

    allocate(rambo_masses(rambo_max_particles))
    allocate(rambo_momenta(0:3, rambo_max_particles))
    allocate(momentum(0:3, nexternal))
    do particle = nincoming + 1, nexternal - 1
      rambo_masses(particle - nincoming) = generated_masses(particle)
    end do

    iconfig = 1
    ichan = 1
    iconfigs(1) = iconfig
    do nfksprocess = 1, fks_configs
      call fks_inc_chooser()
      particle = fks_i_d(nfksprocess)
      if (abs(pdg_type_d(nfksprocess, particle)) == 21 .or. &
          pdg_type_d(nfksprocess, particle) == 22) exit
    end do
    if (nfksprocess > fks_configs) return

    call fks_inc_chooser()
    call leshouche_inc_chooser()
    call setfksfactor()

    nfail = 0
    npoints_checked = 0
    call force_stability_check(.true.)
    call collier_compute_uv_poles(.true.)
    call collier_compute_ir_poles(.true.)

200 continue
    calculatedborn = .false.

    if (nincoming == 1) then
      call rambo_impl(0, nexternal - nincoming - 1, &
           generated_masses(1), rambo_masses, rambo_momenta)
      p_born(0, 1) = generated_masses(1)
      p_born(1, 1) = 0d0
      p_born(2, 1) = 0d0
      p_born(3, 1) = 0d0
    else if (nincoming == 2) then
      if (nexternal - nincoming - 1 == 1) then
        p_born(0, 1) = generated_masses(3) / 2d0
        p_born(1, 1) = 0d0
        p_born(2, 1) = 0d0
        p_born(3, 1) = generated_masses(3) / 2d0
        if (generated_masses(1) > 0d0) then
          p_born(3, 1) = sqrt(generated_masses(3)**2 / 4d0 - &
               generated_masses(1)**2)
        end if
        p_born(0, 2) = generated_masses(3) / 2d0
        p_born(1, 2) = 0d0
        p_born(2, 2) = 0d0
        p_born(3, 2) = -generated_masses(3) / 2d0
        if (generated_masses(2) > 0d0) then
          p_born(3, 2) = -sqrt(generated_masses(3)**2 / 4d0 - &
               generated_masses(1)**2)
        end if
        rambo_momenta(0, 1) = generated_masses(3)
        rambo_momenta(1, 1) = 0d0
        rambo_momenta(2, 1) = 0d0
        rambo_momenta(3, 1) = 0d0
      else
        call rambo_impl(0, nexternal - nincoming - 1, energy, &
             rambo_masses, rambo_momenta)
        p_born(0, 1) = energy / 2d0
        p_born(1, 1) = 0d0
        p_born(2, 1) = 0d0
        p_born(3, 1) = energy / 2d0
        if (generated_masses(1) > 0d0) then
          p_born(3, 1) = sqrt(energy**2 / 4d0 - &
               generated_masses(1)**2)
        end if
        p_born(0, 2) = energy / 2d0
        p_born(1, 2) = 0d0
        p_born(2, 2) = 0d0
        p_born(3, 2) = -energy / 2d0
        if (generated_masses(2) > 0d0) then
          p_born(3, 2) = -sqrt(energy**2 / 4d0 - &
               generated_masses(1)**2)
        end if
      end if
    else
      write (*, *) 'INVALID NUMBER OF INCOMING PARTICLES', nincoming
      stop
    end if

    do component = 0, 3
      do particle = nincoming + 1, nexternal - 1
        p_born(component, particle) = &
             rambo_momenta(component, particle - nincoming)
      end do
    end do

    call update_as_param()
    call sborn(p_born, born_weight)
    call binothlha(p_born, born_weight, virtual_weight)
    if (npoints_checked == 0) then
      if (mod(ret_code_ml, 100) / 10 == 3 .or. &
          mod(ret_code_ml, 100) / 10 == 4) then
        write (*, *) 'INITIALIZATION POINT.'
        write (*, *) 'RESULTS FROM INITIALIZATION POINTS WILL NOT ' // &
             'BE USED FOR STATISTICS'
        go to 200
      end if
    end if

    call setrun_model_strong_coupling(strong_coupling)
    write (*, *) 'MU_R    = ', renormalization_scale
    write (*, *) 'ALPHA_S = ', strong_coupling**2 / (4d0 * pi)
    npoints_checked = npoints_checked + 1

    do component = 0, 3
      do particle = 1, nexternal - 1
        momentum(component, particle) = p_born(component, particle)
      end do
      momentum(component, nexternal) = 0d0
    end do

    if (tolerance < 0d0) then
      write (*, *) 'PASSED', tolerance
    else
      if (polecheck_passed) then
        write (*, *) 'PASSED', tolerance
      else
        write (*, *) 'FAILED', tolerance
        nfail = nfail + 1
      end if
    end if
    write (*, *)

    if (npoints_checked < npoints) go to 200

    write (*, *) 'NUMBER OF POINTS PASSING THE CHECK', npoints - nfail
    write (*, *) 'NUMBER OF POINTS FAILING THE CHECK', nfail
    write (*, *) 'TOLERANCE ', tolerance

  end subroutine run_check_poles


  subroutine initialize_rambo()
    implicit none
    integer :: particle

    if (rambo_initialized) return
    allocate(rambo_log_weights(rambo_max_particles))
    rambo_two_pi = 8d0 * atan(1d0)
    rambo_log_pi = log(rambo_two_pi / 4d0)
    rambo_log_weights(2) = rambo_log_pi
    do particle = 3, rambo_max_particles
      rambo_log_weights(particle) = &
           rambo_log_weights(particle - 1) + rambo_log_pi - &
           2d0 * log(dble(particle - 2))
    end do
    do particle = 3, rambo_max_particles
      rambo_log_weights(particle) = rambo_log_weights(particle) - &
           log(dble(particle - 1))
    end do
    rambo_warnings = 0
    rambo_initialized = .true.
  end subroutine initialize_rambo


  subroutine rambo_impl(lflag, number_of_particles, total_energy, &
       masses, momenta)
    implicit none
    integer, intent(in) :: lflag, number_of_particles
    double precision, intent(in) :: total_energy
    double precision, intent(in) :: masses(rambo_max_particles)
    double precision, intent(inout) :: &
         momenta(0:3, rambo_max_particles)
    double precision, allocatable :: q(:, :), momentum_squared(:)
    double precision, allocatable :: mass_squared(:), energies(:)
    double precision, allocatable :: momentum_norm(:)
    double precision :: random_one, random_two, random_three, random_four
    double precision :: cosine_theta, sine_theta, azimuth
    double precision :: total_vector(4), boost(3), boost_dot_q
    double precision :: invariant_mass, gamma, factor_a, scale_factor
    double precision :: total_mass, massless_weight, mass_weight
    double precision :: weight_product, weight_sum, event_weight
    double precision :: coefficient_one, coefficient_two
    double precision :: coefficient_three, trial_function
    double precision :: trial_derivative, scale_squared
    double precision :: maximum_scale, target_accuracy
    integer :: nonzero_masses, iteration, particle, component

    call initialize_rambo()
    if (number_of_particles <= 1 .or. &
        number_of_particles >= rambo_max_particles + 1) then
      print 1001, number_of_particles
      stop
    end if

    allocate(q(4, number_of_particles))
    allocate(momentum_squared(number_of_particles))
    allocate(mass_squared(number_of_particles))
    allocate(energies(number_of_particles))
    allocate(momentum_norm(number_of_particles))

    total_mass = 0d0
    nonzero_masses = 0
    do particle = 1, number_of_particles
      if (abs(masses(particle)) > 0d0) then
        nonzero_masses = nonzero_masses + 1
      end if
      total_mass = total_mass + abs(masses(particle))
    end do
    if (total_mass > total_energy) then
      print 1002, total_mass, total_energy
      stop
    end if

    if (lflag == 1) then
      massless_weight = exp((2d0 * number_of_particles - 4d0) * &
           log(total_energy) + rambo_log_weights(number_of_particles))
      do particle = 1, number_of_particles
        momentum_norm(particle) = sqrt(momenta(1, particle)**2 + &
             momenta(2, particle)**2 + momenta(3, particle)**2)
      end do
      coefficient_one = 0d0
      coefficient_three = 0d0
      coefficient_two = 1d0
      do particle = 1, number_of_particles
        coefficient_one = coefficient_one + &
             momentum_norm(particle) / total_energy
        coefficient_two = coefficient_two * &
             momentum_norm(particle) / momenta(0, particle)
        coefficient_three = coefficient_three + &
             momentum_norm(particle)**2 / momenta(0, particle) / &
             total_energy
      end do
      mass_weight = coefficient_one**(2 * number_of_particles - 3) * &
           coefficient_two / coefficient_three
      event_weight = 1d0 / massless_weight / mass_weight
      deallocate(q, momentum_squared, mass_squared, energies, momentum_norm)
      return
    end if

    do particle = 1, number_of_particles
      call rans_impl(random_one)
      call rans_impl(random_two)
      call rans_impl(random_three)
      call rans_impl(random_four)
      cosine_theta = 2d0 * random_one - 1d0
      sine_theta = sqrt(1d0 - cosine_theta * cosine_theta)
      azimuth = rambo_two_pi * random_two
      q(4, particle) = -log(random_three * random_four)
      q(3, particle) = q(4, particle) * cosine_theta
      q(2, particle) = q(4, particle) * sine_theta * cos(azimuth)
      q(1, particle) = q(4, particle) * sine_theta * sin(azimuth)
    end do

    total_vector = 0d0
    do particle = 1, number_of_particles
      do component = 1, 4
        total_vector(component) = total_vector(component) + &
             q(component, particle)
      end do
    end do
    invariant_mass = sqrt(total_vector(4)**2 - total_vector(3)**2 - &
         total_vector(2)**2 - total_vector(1)**2)
    do component = 1, 3
      boost(component) = -total_vector(component) / invariant_mass
    end do
    gamma = total_vector(4) / invariant_mass
    factor_a = 1d0 / (1d0 + gamma)
    scale_factor = total_energy / invariant_mass

    do particle = 1, number_of_particles
      boost_dot_q = boost(1) * q(1, particle) + &
           boost(2) * q(2, particle) + boost(3) * q(3, particle)
      do component = 1, 3
        momenta(component, particle) = scale_factor * &
             (q(component, particle) + boost(component) * &
             (q(4, particle) + factor_a * boost_dot_q))
      end do
      momenta(0, particle) = scale_factor * &
           (gamma * q(4, particle) + boost_dot_q)
    end do

    massless_weight = rambo_log_pi
    if (number_of_particles /= 2) then
      massless_weight = (2d0 * number_of_particles - 4d0) * &
           log(total_energy) + rambo_log_weights(number_of_particles)
    end if
    if (massless_weight < -180d0) then
      if (rambo_warnings(1) <= 5) print 1004, massless_weight
      rambo_warnings(1) = rambo_warnings(1) + 1
    end if
    if (massless_weight > 174d0) then
      if (rambo_warnings(2) <= 5) print 1005, massless_weight
      rambo_warnings(2) = rambo_warnings(2) + 1
    end if
    if (nonzero_masses == 0) then
      event_weight = 1d0 / exp(massless_weight)
      deallocate(q, momentum_squared, mass_squared, energies, momentum_norm)
      return
    end if

    maximum_scale = sqrt(1d0 - (total_mass / total_energy)**2)
    do particle = 1, number_of_particles
      mass_squared(particle) = masses(particle)**2
      momentum_squared(particle) = momenta(0, particle)**2
    end do
    iteration = 0
    scale_factor = maximum_scale
    target_accuracy = total_energy * rambo_accuracy

302 continue
    trial_function = -total_energy
    trial_derivative = 0d0
    scale_squared = scale_factor * scale_factor
    do particle = 1, number_of_particles
      energies(particle) = sqrt(mass_squared(particle) + &
           scale_squared * momentum_squared(particle))
      trial_function = trial_function + energies(particle)
      trial_derivative = trial_derivative + &
           momentum_squared(particle) / energies(particle)
    end do
    if (abs(trial_function) > target_accuracy) then
      iteration = iteration + 1
      if (iteration <= rambo_iteration_limit) then
        scale_factor = scale_factor - trial_function / &
             (scale_factor * trial_derivative)
        go to 302
      end if
      print 1006, rambo_iteration_limit
    end if

    do particle = 1, number_of_particles
      momentum_norm(particle) = scale_factor * momenta(0, particle)
      do component = 1, 3
        momenta(component, particle) = &
             scale_factor * momenta(component, particle)
      end do
      momenta(0, particle) = energies(particle)
    end do

    weight_product = 1d0
    weight_sum = 0d0
    do particle = 1, number_of_particles
      weight_product = weight_product * &
           momentum_norm(particle) / energies(particle)
      weight_sum = weight_sum + &
           momentum_norm(particle)**2 / energies(particle)
    end do
    mass_weight = (2d0 * number_of_particles - 3d0) * &
         log(scale_factor) + &
         log(weight_product / weight_sum * total_energy)
    massless_weight = massless_weight + mass_weight
    if (massless_weight < -180d0) then
      if (rambo_warnings(3) <= 5) print 1004, massless_weight
      rambo_warnings(3) = rambo_warnings(3) + 1
    end if
    if (massless_weight > 174d0) then
      if (rambo_warnings(4) <= 5) print 1005, massless_weight
      rambo_warnings(4) = rambo_warnings(4) + 1
    end if
    event_weight = 1d0 / exp(massless_weight)
    deallocate(q, momentum_squared, mass_squared, energies, momentum_norm)
    return

1001 format(' RAMBO FAILS: # OF PARTICLES =', I5, ' IS NOT ALLOWED')
1002 format(' RAMBO FAILS: TOTAL MASS =', D15.6, ' IS NOT', &
         ' SMALLER THAN TOTAL ENERGY =', D15.6)
1004 format(' RAMBO WARNS: WEIGHT = EXP(', F20.9, ') MAY UNDERFLOW')
1005 format(' RAMBO WARNS: WEIGHT = EXP(', F20.9, ') MAY  OVERFLOW')
1006 format(' RAMBO WARNS:', I3, ' ITERATIONS DID NOT GIVE THE', &
         ' DESIRED ACCURACY =', D15.6)
  end subroutine rambo_impl


  subroutine rans_impl(random_value)
    implicit none
    double precision, intent(out) :: random_value

    random_value = ran2()
  end subroutine rans_impl


  subroutine fail_check_poles(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, *) 'ERROR in check_poles: ', trim(message)
    stop 1
  end subroutine fail_check_poles

end module check_poles_module

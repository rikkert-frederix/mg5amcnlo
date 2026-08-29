module multiplicative_event_materialization
  use multiplicative_nlo_decay, only: multiplicative_nlo_workspace
  use spin_density_matrix_results, only: spin_density_real_branch
  use factorized_phase_space, only: fetch_factorized_block_momenta, &
       fetch_factorized_embedded_momenta
  use factorized_block_kinematics, only: &
       boost_factorized_momentum_from_rest, factorized_minkowski_square
  implicit none
  private

  public :: multiplicative_leaf_particle_capacity
  public :: materialize_multiplicative_leaf

  interface
    integer function sdm_leaf_max_particle_count()
    end function sdm_leaf_max_particle_count

    integer function sdm_leaf_max_block_size()
    end function sdm_leaf_max_block_size

    subroutine sdm_leaf_block_layout( &
         position, branch, configuration, particle_count, statuses, &
         target_kinds, target_ids, pdgs)
      integer, intent(in) :: position, branch, configuration
      integer, intent(out) :: particle_count
      integer, intent(out) :: statuses(*), target_kinds(*)
      integer, intent(out) :: target_ids(*), pdgs(*)
    end subroutine sdm_leaf_block_layout
  end interface

contains

  integer function multiplicative_leaf_particle_capacity()
    multiplicative_leaf_particle_capacity = &
         sdm_leaf_max_particle_count()
    if (multiplicative_leaf_particle_capacity < 1) then
      call fail_materialization( &
           'the generated leaf particle capacity is invalid')
    end if
  end function multiplicative_leaf_particle_capacity


  subroutine materialize_multiplicative_leaf( &
       workspace, event_slot, momenta, statuses, pdgs, &
       particle_from_decay, particle_count)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    integer, intent(in) :: event_slot
    double precision, intent(out) :: momenta(0:, :)
    integer, intent(out) :: statuses(:), pdgs(:)
    logical, intent(out) :: particle_from_decay(:)
    integer, intent(out) :: particle_count
    integer :: production_position, position
    double precision :: unused_parent(0:3)

    if (size(momenta, 1) < 4 .or. &
        size(momenta, 2) < multiplicative_leaf_particle_capacity() .or. &
        size(statuses) < size(momenta, 2) .or. &
        size(pdgs) < size(momenta, 2) .or. &
        size(particle_from_decay) < size(momenta, 2)) then
      call fail_materialization('leaf scratch storage is too small')
    end if
    if (workspace%component_count < 1 .or. &
        .not. allocated(workspace%component_ids) .or. &
        .not. allocated(workspace%branch_by_component) .or. &
        .not. allocated(workspace%real_configuration_by_component)) then
      call fail_materialization('the multiplicative workspace is invalid')
    end if

    production_position = 0
    do position = 1, workspace%component_count
      if (workspace%component_ids(position) == 0) then
        production_position = position
        exit
      end if
    end do
    if (production_position == 0) then
      call fail_materialization('the production component is absent')
    end if

    momenta = 0d0
    statuses = 0
    pdgs = 0
    particle_from_decay = .false.
    particle_count = 0
    unused_parent = 0d0
    call append_component(production_position, unused_parent, .false.)

  contains

    recursive subroutine append_component( &
         component_position, new_parent, has_parent)
      integer, intent(in) :: component_position
      double precision, intent(in) :: new_parent(0:3)
      logical, intent(in) :: has_parent
      double precision, allocatable :: block_momenta(:, :)
      integer, allocatable :: block_statuses(:), target_kinds(:)
      integer, allocatable :: target_ids(:), block_pdgs(:)
      double precision :: old_parent(0:3), inverse_parent(0:3)
      double precision :: rest_momentum(0:3), event_momentum(0:3)
      double precision :: parent_mass2, parent_mass, tolerance
      integer :: block_size, block_count, branch, configuration
      integer :: leg, initial_count, child_position
      logical :: available

      block_size = sdm_leaf_max_block_size()
      if (block_size < 1) then
        call fail_materialization('the generated block size is invalid')
      end if
      allocate(block_momenta(0:3, block_size))
      allocate(block_statuses(block_size), target_kinds(block_size))
      allocate(target_ids(block_size), block_pdgs(block_size))
      branch = workspace%branch_by_component(component_position)
      configuration = 0
      if (branch == spin_density_real_branch) then
        configuration = workspace%real_configuration_by_component( &
             component_position)
        if (configuration < 1) then
          call fail_materialization( &
               'a selected real branch has no FKS configuration')
        end if
      end if
      call sdm_leaf_block_layout( &
           component_position, branch, configuration, block_count, &
           block_statuses, target_kinds, target_ids, block_pdgs)
      if (block_count < 1 .or. block_count > block_size) then
        call fail_materialization('a generated block count is invalid')
      end if
      if (branch == spin_density_real_branch) then
        call fetch_factorized_embedded_momenta( &
             event_slot, workspace%component_ids(component_position), &
             block_count, block_momenta(:, 1:block_count), available)
      else
        ! Counterevents may retain a real-context zero-momentum leg in the
        ! embedding cache.  The reduced B branch is deliberately represented
        ! by its canonical n-body matrix-element block: IR-safe observables
        ! cannot distinguish the projected real layout from this one.
        call fetch_factorized_block_momenta( &
             event_slot, workspace%component_ids(component_position), &
             block_count, block_momenta(:, 1:block_count), available)
      end if
      if (.not. available) then
        write (*, '(a,i0,a,i0,a,i0)') &
             'ERROR in multiplicative_event_materialization: block ', &
             workspace%component_ids(component_position), &
             ' selected branch ', branch, ' requires ', block_count
        stop 1
      end if

      if (.not. has_parent) then
        do leg = 1, block_count
          if (target_kinds(leg) == 1) then
            if (block_statuses(leg) < 0) then
              call fail_materialization( &
                   'an incoming production leg targets a decay node')
            end if
            child_position = component_position_for_id(target_ids(leg))
            call append_component( &
                 child_position, block_momenta(:, leg), .true.)
          else
            call append_particle(block_momenta(:, leg), &
                 block_statuses(leg), block_pdgs(leg), .false.)
          end if
        end do
        deallocate(block_momenta, block_statuses, target_kinds, &
                   target_ids, block_pdgs)
        return
      end if

      initial_count = 0
      old_parent = 0d0
      do leg = 1, block_count
        if (block_statuses(leg) >= 0) cycle
        initial_count = initial_count + 1
        old_parent = block_momenta(:, leg)
      end do
      if (initial_count /= 1) then
        call fail_materialization( &
             'a decay block does not have exactly one parent')
      end if
      parent_mass2 = factorized_minkowski_square(old_parent)
      if (parent_mass2 <= 0d0) then
        call fail_materialization('a decay parent is not time-like')
      end if
      parent_mass = sqrt(parent_mass2)
      tolerance = 1d-8*max(1d0, parent_mass2)
      if (abs(factorized_minkowski_square(new_parent) - parent_mass2) > &
          tolerance) then
        call fail_materialization( &
             'combined branches disagree on a resonance mass')
      end if
      inverse_parent = old_parent
      inverse_parent(1:3) = -inverse_parent(1:3)

      do leg = 1, block_count
        if (block_statuses(leg) < 0) cycle
        call boost_factorized_momentum_from_rest( &
             block_momenta(:, leg), inverse_parent, parent_mass, &
             rest_momentum)
        call boost_factorized_momentum_from_rest( &
             rest_momentum, new_parent, parent_mass, event_momentum)
        if (target_kinds(leg) == 1) then
          child_position = component_position_for_id(target_ids(leg))
          call append_component(child_position, event_momentum, .true.)
        else
          call append_particle(event_momentum, block_statuses(leg), &
               block_pdgs(leg), .true.)
        end if
      end do
      deallocate(block_momenta, block_statuses, target_kinds, &
                 target_ids, block_pdgs)
    end subroutine append_component


    integer function component_position_for_id(component_id)
      integer, intent(in) :: component_id
      integer :: candidate

      component_position_for_id = 0
      do candidate = 1, workspace%component_count
        if (workspace%component_ids(candidate) == component_id) then
          component_position_for_id = candidate
          return
        end if
      end do
      call fail_materialization('a leaf targets an unknown decay node')
    end function component_position_for_id


    subroutine append_particle(momentum, status, pdg, from_decay)
      double precision, intent(in) :: momentum(0:3)
      integer, intent(in) :: status, pdg
      logical, intent(in) :: from_decay

      if (status /= -1 .and. status /= 1) then
        call fail_materialization('a materialized particle status is invalid')
      end if
      if (particle_count >= size(momenta, 2)) then
        call fail_materialization('a materialized leaf exceeds its capacity')
      end if
      particle_count = particle_count + 1
      momenta(:, particle_count) = momentum
      statuses(particle_count) = status
      pdgs(particle_count) = pdg
      particle_from_decay(particle_count) = from_decay
    end subroutine append_particle

  end subroutine materialize_multiplicative_leaf


  subroutine fail_materialization(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_event_materialization: '//trim(message)
    stop 1
  end subroutine fail_materialization

end module multiplicative_event_materialization

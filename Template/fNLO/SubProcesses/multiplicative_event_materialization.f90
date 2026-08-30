module multiplicative_event_materialization
  use multiplicative_nlo_decay, only: multiplicative_nlo_workspace
  use spin_density_matrix_results, only: spin_density_real_branch
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
       workspace, momenta, statuses, pdgs, &
       particle_from_decay, particle_count)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    double precision, intent(out) :: momenta(0:, :)
    integer, intent(out) :: statuses(:), pdgs(:)
    logical, intent(out) :: particle_from_decay(:)
    integer, intent(out) :: particle_count
    integer :: production_position, position
    double precision :: unused_parent(0:3)
    double precision :: scratch_momenta( &
         0:3, size(momenta, 2), workspace%component_count)
    integer :: scratch_statuses( &
         size(momenta, 2), workspace%component_count)
    integer :: scratch_target_kinds( &
         size(momenta, 2), workspace%component_count)
    integer :: scratch_target_ids( &
         size(momenta, 2), workspace%component_count)
    integer :: scratch_pdgs( &
         size(momenta, 2), workspace%component_count)

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
      if (block_size > size(momenta, 2)) then
        call fail_materialization('the block scratch storage is too small')
      end if
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
           scratch_statuses(1, component_position), &
           scratch_target_kinds(1, component_position), &
           scratch_target_ids(1, component_position), &
           scratch_pdgs(1, component_position))
      if (block_count < 1 .or. block_count > block_size) then
        call fail_materialization('a generated block count is invalid')
      end if
      associate(snapshot => workspace%snapshots( &
           branch, component_position))
        if (branch == spin_density_real_branch) then
          available = snapshot%has_embedded .and. &
               snapshot%embedded_count == block_count
          if (available) scratch_momenta(:, 1:block_count, &
               component_position) = &
               snapshot%embedded_momenta(:, 1:block_count)
        else
          ! Counterevents may retain a real-context zero-momentum leg in the
          ! embedding cache.  The reduced B branch is deliberately represented
          ! by its canonical n-body matrix-element block: IR-safe observables
          ! cannot distinguish the projected real layout from this one.
          available = snapshot%has_block .and. &
               snapshot%block_count == block_count
          if (available) scratch_momenta(:, 1:block_count, &
               component_position) = &
               snapshot%block_momenta(:, 1:block_count)
        end if
      end associate
      if (.not. available) then
        write (*, '(a,i0,a,i0,a,i0)') &
             'ERROR in multiplicative_event_materialization: block ', &
             workspace%component_ids(component_position), &
             ' selected branch ', branch, ' requires ', block_count
        stop 1
      end if

      if (.not. has_parent) then
        do leg = 1, block_count
          if (scratch_target_kinds(leg, component_position) == 1) then
            if (scratch_statuses(leg, component_position) < 0) then
              call fail_materialization( &
                   'an incoming production leg targets a decay node')
            end if
            child_position = component_position_for_id( &
                 scratch_target_ids(leg, component_position))
            call append_component( &
                 child_position, &
                 scratch_momenta(:, leg, component_position), .true.)
          else
            call append_particle( &
                 scratch_momenta(:, leg, component_position), &
                 scratch_statuses(leg, component_position), &
                 scratch_pdgs(leg, component_position), .false.)
          end if
        end do
        return
      end if

      initial_count = 0
      old_parent = 0d0
      do leg = 1, block_count
        if (scratch_statuses(leg, component_position) >= 0) cycle
        initial_count = initial_count + 1
        old_parent = scratch_momenta(:, leg, component_position)
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
        write (*, '(a,i0,a,i0,a,es24.16,a,es24.16)') &
             'ERROR in multiplicative_event_materialization: component ', &
             workspace%component_ids(component_position), ' branch ', &
             branch, ' has block mass^2 ', parent_mass2, &
             ' but parent mass^2 ', &
             factorized_minkowski_square(new_parent)
        call fail_materialization( &
             'combined branches disagree on a resonance mass')
      end if
      inverse_parent = old_parent
      inverse_parent(1:3) = -inverse_parent(1:3)

      do leg = 1, block_count
        if (scratch_statuses(leg, component_position) < 0) cycle
        call boost_factorized_momentum_from_rest( &
             scratch_momenta(:, leg, component_position), &
             inverse_parent, parent_mass, &
             rest_momentum)
        call boost_factorized_momentum_from_rest( &
             rest_momentum, new_parent, parent_mass, event_momentum)
        if (scratch_target_kinds(leg, component_position) == 1) then
          child_position = component_position_for_id( &
               scratch_target_ids(leg, component_position))
          call append_component(child_position, event_momentum, .true.)
        else
          call append_particle(event_momentum, &
               scratch_statuses(leg, component_position), &
               scratch_pdgs(leg, component_position), .true.)
        end if
      end do
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

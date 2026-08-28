module multiplicative_kinematics
  use process_dimensions, only: nexternal, validate_process_dimensions
  use factorized_block_kinematics, only: &
       boost_from_rest => boost_factorized_momentum_from_rest, &
       minkowski_square => factorized_minkowski_square
  use kin_functions_module, only: transverse_energy => et_impl
  use factorized_phase_space, only: factorized_no_target, &
       factorized_visible_target, factorized_block_target, &
       factorized_local_count, fetch_factorized_local_momenta, &
       fetch_factorized_local_layout, store_factorized_block_momenta, &
       fetch_factorized_block_momenta, &
       store_factorized_embedded_momenta
  implicit none
  private

  public :: realize_factorized_event_tuple
  public :: materialize_factorized_event_tuple
  public :: factorized_production_scale_sums
  public :: factorized_visible_scale_sums

contains

  subroutine realize_factorized_event_tuple(event_slots, pass)
    integer, intent(in) :: event_slots(0:)
    logical, intent(out) :: pass
    double precision :: root_momenta(0:3, nexternal)
    integer :: root_count, leg
    integer :: pdgs(nexternal), target_kinds(nexternal)
    integer :: target_ids(nexternal)
    logical :: particle_is_final(nexternal), available
    logical :: active(0:nexternal), realized(0:nexternal)

    call validate_process_dimensions()
    call validate_slot_tuple(event_slots)
    active = .false.
    realized = .false.
    root_count = factorized_local_count(event_slots(0), 0)
    if (root_count < 1) then
      call fail_multiplicative_kinematics( &
           'the selected production-local configuration is unavailable')
    end if
    call fetch_factorized_local_momenta( &
         event_slots(0), 0, root_count, root_momenta(:, 1:root_count), &
         available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'the selected production-local momenta are unavailable')
    end if
    call fetch_factorized_local_layout( &
         event_slots(0), 0, root_count, pdgs(1:root_count), &
         particle_is_final(1:root_count), target_kinds(1:root_count), &
         target_ids(1:root_count), available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'the selected production-local layout is unavailable')
    end if

    call store_factorized_block_momenta( &
         event_slots(0), 0, root_count, root_momenta)
    call store_factorized_embedded_momenta( &
         event_slots(0), 0, root_count, root_momenta)
    active(0) = .true.
    realized(0) = .true.
    pass = .true.
    do leg = 1, root_count
      if (.not. particle_is_final(leg)) cycle
      if (target_kinds(leg) == factorized_block_target) then
        call realize_decay_block( &
             event_slots, target_ids(leg), root_momenta(:, leg), &
             active, realized, pass)
        if (.not. pass) return
      end if
    end do
    active(0) = .false.
  end subroutine realize_factorized_event_tuple


  recursive subroutine realize_decay_block( &
       event_slots, block, parent_momentum, active, realized, pass)
    integer, intent(in) :: event_slots(0:), block
    double precision, intent(in) :: parent_momentum(0:3)
    logical, intent(inout) :: active(0:), realized(0:)
    logical, intent(out) :: pass
    double precision :: local_momenta(0:3, nexternal)
    double precision :: boosted_momenta(0:3, nexternal)
    integer :: count, leg, incoming_leg, incoming_count, event_slot
    integer :: pdgs(nexternal), target_kinds(nexternal)
    integer :: target_ids(nexternal)
    logical :: particle_is_final(nexternal), available
    double precision :: parent_rest(0:3), parent_mass, tolerance

    if (block < 1 .or. block > nexternal) then
      call fail_multiplicative_kinematics( &
           'a decay target has an invalid block index')
    end if
    if (active(block)) then
      call fail_multiplicative_kinematics( &
           'the factorized block layout contains a cycle')
    end if
    if (realized(block)) then
      call fail_multiplicative_kinematics( &
           'one physical decay block has two selected parents')
    end if
    active(block) = .true.
    event_slot = event_slots(block)
    count = factorized_local_count(event_slot, block)
    if (count < 2) then
      call fail_multiplicative_kinematics( &
           'a selected decay-local configuration is unavailable')
    end if
    call fetch_factorized_local_momenta( &
         event_slot, block, count, local_momenta(:, 1:count), available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'the selected decay-local momenta are unavailable')
    end if
    call fetch_factorized_local_layout( &
         event_slot, block, count, pdgs(1:count), &
         particle_is_final(1:count), target_kinds(1:count), &
         target_ids(1:count), available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'the selected decay-local layout is unavailable')
    end if

    incoming_count = 0
    incoming_leg = 0
    do leg = 1, count
      if (particle_is_final(leg)) cycle
      incoming_count = incoming_count + 1
      incoming_leg = leg
    end do
    if (incoming_count /= 1) then
      call fail_multiplicative_kinematics( &
           'a decay-local block does not have exactly one parent')
    end if
    parent_rest = local_momenta(:, incoming_leg)
    parent_mass = parent_rest(0)
    tolerance = 1d-10*max(1d0, abs(parent_mass))
    if (parent_mass <= 0d0 .or. &
        maxval(abs(parent_rest(1:3))) > tolerance .or. &
        abs(minkowski_square(parent_rest) - parent_mass**2) > &
        tolerance*max(1d0, parent_mass)) then
      call fail_multiplicative_kinematics( &
           'a decay-local parent is not a physical rest-frame momentum')
    end if
    if (abs(minkowski_square(parent_momentum) - parent_mass**2) > &
        1d-8*max(1d0, parent_mass**2)) then
      call fail_multiplicative_kinematics( &
           'a selected parent and child block have different masses')
    end if

    boosted_momenta = 0d0
    do leg = 1, count
      if (particle_is_final(leg)) then
        call boost_from_rest( &
             local_momenta(:, leg), parent_momentum, parent_mass, &
             boosted_momenta(:, leg))
      else
        boosted_momenta(:, leg) = parent_momentum
      end if
    end do
    call validate_decay_conservation( &
         count, particle_is_final, boosted_momenta, parent_momentum)
    call store_factorized_block_momenta( &
         event_slot, block, count, boosted_momenta)
    call store_factorized_embedded_momenta( &
         event_slot, block, count, boosted_momenta)
    realized(block) = .true.

    pass = .true.
    do leg = 1, count
      if (.not. particle_is_final(leg)) cycle
      if (target_kinds(leg) == factorized_block_target) then
        call realize_decay_block( &
             event_slots, target_ids(leg), boosted_momenta(:, leg), &
             active, realized, pass)
        if (.not. pass) return
      end if
    end do
    active(block) = .false.
  end subroutine realize_decay_block


  subroutine materialize_factorized_event_tuple( &
       event_slots, visible_capacity, visible_count, momenta, pdgs, pass, &
       origin_blocks)
    integer, intent(in) :: event_slots(0:), visible_capacity
    integer, intent(out) :: visible_count
    double precision, intent(out) :: momenta(0:3, visible_capacity)
    integer, intent(out) :: pdgs(visible_capacity)
    logical, intent(out) :: pass
    integer, intent(out), optional :: origin_blocks(visible_capacity)
    logical :: occupied(visible_capacity), active(0:nexternal)
    integer :: materialized_origins(visible_capacity)
    integer :: target

    call validate_slot_tuple(event_slots)
    if (visible_capacity < 1) then
      call fail_multiplicative_kinematics( &
           'the visible-event capacity is empty')
    end if
    momenta = 0d0
    pdgs = 0
    materialized_origins = -1
    occupied = .false.
    active = .false.
    call materialize_block( &
         event_slots, 0, visible_capacity, momenta, pdgs, occupied, &
         materialized_origins, active, pass)
    if (.not. pass) return
    visible_count = 0
    do target = 1, visible_capacity
      if (occupied(target)) visible_count = target
    end do
    if (present(origin_blocks)) origin_blocks = materialized_origins
  end subroutine materialize_factorized_event_tuple


  subroutine factorized_production_scale_sums( &
       event_slots, sum_transverse_energy, sum_transverse_mass)
    integer, intent(in) :: event_slots(0:)
    double precision, intent(out) :: sum_transverse_energy
    double precision, intent(out) :: sum_transverse_mass

    call validate_slot_tuple(event_slots)
    sum_transverse_energy = 0d0
    sum_transverse_mass = 0d0
    call accumulate_block_scale_sums( &
         event_slots, 0, .false., sum_transverse_energy, &
         sum_transverse_mass)
  end subroutine factorized_production_scale_sums


  subroutine factorized_visible_scale_sums( &
       event_slots, sum_transverse_energy, sum_transverse_mass)
    integer, intent(in) :: event_slots(0:)
    double precision, intent(out) :: sum_transverse_energy
    double precision, intent(out) :: sum_transverse_mass

    call validate_slot_tuple(event_slots)
    sum_transverse_energy = 0d0
    sum_transverse_mass = 0d0
    call accumulate_block_scale_sums( &
         event_slots, 0, .true., sum_transverse_energy, &
         sum_transverse_mass)
  end subroutine factorized_visible_scale_sums


  recursive subroutine accumulate_block_scale_sums( &
       event_slots, block, expand_decays, sum_transverse_energy, &
       sum_transverse_mass)
    integer, intent(in) :: event_slots(0:), block
    logical, intent(in) :: expand_decays
    double precision, intent(inout) :: sum_transverse_energy
    double precision, intent(inout) :: sum_transverse_mass
    double precision :: block_momenta(0:3, nexternal)
    integer :: count, leg, event_slot
    integer :: local_pdgs(nexternal), target_kinds(nexternal)
    integer :: target_ids(nexternal)
    logical :: particle_is_final(nexternal), available

    event_slot = event_slots(block)
    count = factorized_local_count(event_slot, block)
    if (count < 1) then
      call fail_multiplicative_kinematics( &
           'a scale-sum block has no selected local configuration')
    end if
    call fetch_factorized_block_momenta( &
         event_slot, block, count, block_momenta(:, 1:count), available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'scale sums require a tuple-realized block')
    end if
    call fetch_factorized_local_layout( &
         event_slot, block, count, local_pdgs(1:count), &
         particle_is_final(1:count), target_kinds(1:count), &
         target_ids(1:count), available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'a scale-sum block has no local layout')
    end if

    do leg = 1, count
      if (.not. particle_is_final(leg)) cycle
      if (expand_decays .and. &
          target_kinds(leg) == factorized_block_target) then
        call accumulate_block_scale_sums( &
             event_slots, target_ids(leg), .true., &
             sum_transverse_energy, sum_transverse_mass)
      else
        sum_transverse_energy = sum_transverse_energy + &
             transverse_energy(block_momenta(:, leg))
        sum_transverse_mass = sum_transverse_mass + sqrt(max(0d0, &
             (block_momenta(0, leg) + block_momenta(3, leg))* &
             (block_momenta(0, leg) - block_momenta(3, leg))))
      end if
    end do
  end subroutine accumulate_block_scale_sums


  recursive subroutine materialize_block( &
       event_slots, block, visible_capacity, momenta, pdgs, occupied, &
       origin_blocks, active, pass)
    integer, intent(in) :: event_slots(0:), block, visible_capacity
    double precision, intent(inout) :: momenta(0:3, visible_capacity)
    integer, intent(inout) :: pdgs(visible_capacity)
    logical, intent(inout) :: occupied(visible_capacity)
    integer, intent(inout) :: origin_blocks(visible_capacity)
    logical, intent(inout) :: active(0:)
    logical, intent(out) :: pass
    double precision :: block_momenta(0:3, nexternal)
    integer :: count, leg, target, event_slot
    integer :: local_pdgs(nexternal), target_kinds(nexternal)
    integer :: target_ids(nexternal)
    logical :: particle_is_final(nexternal), available

    if (block < 0 .or. block > nexternal) then
      call fail_multiplicative_kinematics( &
           'a materialized block index is invalid')
    end if
    if (active(block)) then
      call fail_multiplicative_kinematics( &
           'the materialized block layout contains a cycle')
    end if
    active(block) = .true.
    event_slot = event_slots(block)
    count = factorized_local_count(event_slot, block)
    if (count < 1) then
      call fail_multiplicative_kinematics( &
           'a materialized block has no selected local configuration')
    end if
    call fetch_factorized_block_momenta( &
         event_slot, block, count, block_momenta(:, 1:count), available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'a materialized block has not been tuple-realized')
    end if
    call fetch_factorized_local_layout( &
         event_slot, block, count, local_pdgs(1:count), &
         particle_is_final(1:count), target_kinds(1:count), &
         target_ids(1:count), available)
    if (.not. available) then
      call fail_multiplicative_kinematics( &
           'a materialized block has no local layout')
    end if

    pass = .true.
    do leg = 1, count
      if (.not. particle_is_final(leg)) then
        ! The two incoming production legs are part of the observable event.
        ! A decay parent is not: it is replaced recursively by its children.
        ! Production layouts historically leave incoming targets unset, so
        ! their canonical leg number is the unambiguous visible target.
        if (block /= 0) cycle
        target = leg
        if (target < 1 .or. target > visible_capacity) then
          call fail_multiplicative_kinematics( &
               'an incoming target exceeds the event capacity')
        end if
        if (occupied(target)) then
          call fail_multiplicative_kinematics( &
               'an incoming target is already occupied')
        end if
        momenta(:, target) = block_momenta(:, leg)
        pdgs(target) = local_pdgs(leg)
        origin_blocks(target) = 0
        occupied(target) = .true.
        cycle
      end if
      target = target_ids(leg)
      select case (target_kinds(leg))
      case (factorized_visible_target)
        if (target < 1 .or. target > visible_capacity) then
          call fail_multiplicative_kinematics( &
               'a visible target exceeds the event capacity')
        end if
        if (occupied(target)) then
          write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') &
               'ERROR: target ', target, ' from block ', block, &
               ', local leg ', leg, ', PDG ', local_pdgs(leg), &
               ' is already occupied by PDG ', pdgs(target)
          call fail_multiplicative_kinematics( &
               'two block-local particles share one visible target')
        end if
        momenta(:, target) = block_momenta(:, leg)
        pdgs(target) = local_pdgs(leg)
        origin_blocks(target) = block
        occupied(target) = .true.
      case (factorized_block_target)
        call materialize_block( &
             event_slots, target, visible_capacity, momenta, pdgs, &
             occupied, origin_blocks, active, pass)
        if (.not. pass) return
      case (factorized_no_target)
        call fail_multiplicative_kinematics( &
             'a final block-local particle has no materialization target')
      case default
        call fail_multiplicative_kinematics( &
             'a final block-local particle has an unknown target')
      end select
    end do
    active(block) = .false.
  end subroutine materialize_block


  subroutine validate_decay_conservation( &
       count, particle_is_final, momenta, parent)
    integer, intent(in) :: count
    logical, intent(in) :: particle_is_final(count)
    double precision, intent(in) :: momenta(0:3, nexternal)
    double precision, intent(in) :: parent(0:3)
    double precision :: total(0:3), scale
    integer :: leg

    total = 0d0
    do leg = 1, count
      if (particle_is_final(leg)) total = total + momenta(:, leg)
    end do
    scale = max(1d0, maxval(abs(parent)))
    if (maxval(abs(total - parent)) > 1d-8*scale) then
      call fail_multiplicative_kinematics( &
           'a tuple-realized decay block violates momentum conservation')
    end if
  end subroutine validate_decay_conservation


  subroutine validate_slot_tuple(event_slots)
    integer, intent(in) :: event_slots(0:)

    if (ubound(event_slots, 1) < nexternal) then
      call fail_multiplicative_kinematics( &
           'an event-slot tuple has the wrong size')
    end if
  end subroutine validate_slot_tuple


  subroutine fail_multiplicative_kinematics(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_kinematics: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_kinematics
end module multiplicative_kinematics

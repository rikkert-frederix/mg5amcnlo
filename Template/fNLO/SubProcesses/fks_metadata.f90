module fks_metadata
  use process_dimensions, only: nexternal, fks_configs, nsplitorders, &
       validate_process_dimensions
  implicit none
  private

  logical, save :: metadata_initialized = .false.
  integer, allocatable, save :: fks_i_values(:), fks_j_values(:)
  integer, allocatable, save :: extra_cnt_values(:)
  integer, allocatable, save :: isplitorder_born_values(:)
  integer, allocatable, save :: isplitorder_cnt_values(:)
  integer, allocatable, save :: fks_j_from_i_values(:, :, :)
  integer, allocatable, save :: particle_type_values(:, :)
  integer, allocatable, save :: pdg_type_values(:, :)
  logical, allocatable, save :: split_type_values(:, :)
  logical, allocatable, save :: need_color_links_values(:)

  public :: initialize_fks_metadata
  public :: validate_fks_metadata
  public :: fks_i_d, fks_j_d
  public :: extra_cnt_d, isplitorder_born_d, isplitorder_cnt_d
  public :: fks_j_from_i_d
  public :: particle_type_d, pdg_type_d, split_type_d
  public :: need_color_links_d

contains

  subroutine initialize_fks_metadata(fks_i_in, fks_j_in, extra_cnt_in, &
       isplitorder_born_in, isplitorder_cnt_in, fks_j_from_i_in, &
       particle_type_in, pdg_type_in, split_type_in, need_color_links_in)
    integer, intent(in) :: fks_i_in(:), fks_j_in(:)
    integer, intent(in) :: extra_cnt_in(:)
    integer, intent(in) :: isplitorder_born_in(:)
    integer, intent(in) :: isplitorder_cnt_in(:)
    integer, intent(in) :: fks_j_from_i_in(:, :, 0:)
    integer, intent(in) :: particle_type_in(:, :), pdg_type_in(:, :)
    logical, intent(in) :: split_type_in(:, :)
    logical, intent(in) :: need_color_links_in(:)
    integer :: configuration, emitter, count

    call validate_process_dimensions()
    call validate_input_shapes(fks_i_in, fks_j_in, extra_cnt_in, &
         isplitorder_born_in, isplitorder_cnt_in, fks_j_from_i_in, &
         particle_type_in, pdg_type_in, split_type_in, &
         need_color_links_in)
    call validate_input_values(fks_i_in, fks_j_in, extra_cnt_in, &
         isplitorder_born_in, isplitorder_cnt_in, fks_j_from_i_in, &
         split_type_in)

    if (metadata_initialized) then
      if (.not. same_metadata(fks_i_in, fks_j_in, extra_cnt_in, &
           isplitorder_born_in, isplitorder_cnt_in, fks_j_from_i_in, &
           particle_type_in, pdg_type_in, split_type_in, &
           need_color_links_in)) then
        call fail_validation('metadata were reinitialized with ' // &
             'different values')
      end if
      call validate_fks_metadata()
      return
    end if

    allocate(fks_i_values(fks_configs), fks_j_values(fks_configs))
    allocate(extra_cnt_values(fks_configs))
    allocate(isplitorder_born_values(fks_configs))
    allocate(isplitorder_cnt_values(fks_configs))
    allocate(fks_j_from_i_values(fks_configs, nexternal, 0:nexternal))
    allocate(particle_type_values(fks_configs, nexternal))
    allocate(pdg_type_values(fks_configs, nexternal))
    allocate(split_type_values(fks_configs, nsplitorders))
    allocate(need_color_links_values(fks_configs))

    fks_i_values = fks_i_in
    fks_j_values = fks_j_in
    extra_cnt_values = extra_cnt_in
    isplitorder_born_values = isplitorder_born_in
    isplitorder_cnt_values = isplitorder_cnt_in
    particle_type_values = particle_type_in
    pdg_type_values = pdg_type_in
    split_type_values = split_type_in
    need_color_links_values = need_color_links_in

    ! fks_info.inc DATA-initializes only the emitter row of this table.
    ! Canonical zeroes make all unspecified rows portable and retain the
    ! behavior relied upon by the chooser.
    fks_j_from_i_values = 0
    do configuration = 1, fks_configs
      emitter = fks_i_in(configuration)
      count = fks_j_from_i_in(configuration, emitter, 0)
      fks_j_from_i_values(configuration, emitter, 0:count) = &
           fks_j_from_i_in(configuration, emitter, 0:count)
    end do

    metadata_initialized = .true.
    call validate_fks_metadata()
  end subroutine initialize_fks_metadata


  subroutine validate_fks_metadata()
    integer :: configuration, particle, position, count
    logical :: found_fks_partner

    call validate_process_dimensions()
    if (.not. metadata_initialized) then
      call fail_validation('metadata are not initialized')
    end if
    if (.not. storage_shapes_are_valid()) then
      call fail_validation('allocated metadata have inconsistent shapes')
    end if
    if (any(fks_i_values < 1) .or. any(fks_i_values > nexternal)) then
      call fail_validation('FKS_I_D contains an invalid particle index')
    end if
    if (any(fks_j_values < 1) .or. any(fks_j_values > nexternal)) then
      call fail_validation('FKS_J_D contains an invalid particle index')
    end if
    if (any(extra_cnt_values < 0)) then
      call fail_validation('EXTRA_CNT_D contains a negative identifier')
    end if
    if (any(isplitorder_born_values < 0) .or. &
        any(isplitorder_born_values > nsplitorders) .or. &
        any(isplitorder_cnt_values < 0) .or. &
        any(isplitorder_cnt_values > nsplitorders)) then
      call fail_validation('an extra-counterterm split order is invalid')
    end if

    do configuration = 1, fks_configs
      if (.not. any(split_type_values(configuration, :))) then
        call fail_validation('an FKS configuration has no splitting type')
      end if
      do particle = 1, nexternal
        count = fks_j_from_i_values(configuration, particle, 0)
        if (count < 0 .or. count > nexternal) then
          call fail_validation('FKS_J_FROM_I_D has an invalid count')
        end if
        do position = 1, count
          if (fks_j_from_i_values(configuration, particle, position) < 1 &
              .or. fks_j_from_i_values(configuration, particle, &
              position) > nexternal) then
            call fail_validation('FKS_J_FROM_I_D has an invalid partner')
          end if
        end do
      end do

      particle = fks_i_values(configuration)
      count = fks_j_from_i_values(configuration, particle, 0)
      found_fks_partner = .false.
      do position = 1, count
        if (fks_j_from_i_values(configuration, particle, position) == &
            fks_j_values(configuration)) found_fks_partner = .true.
      end do
      if (.not. found_fks_partner) then
        call fail_validation('FKS_J_D is absent from FKS_J_FROM_I_D')
      end if
    end do
  end subroutine validate_fks_metadata






  integer function fks_i_d(configuration)
    integer, intent(in) :: configuration
    call check_configuration(configuration)
    fks_i_d = fks_i_values(configuration)
  end function fks_i_d


  integer function fks_j_d(configuration)
    integer, intent(in) :: configuration
    call check_configuration(configuration)
    fks_j_d = fks_j_values(configuration)
  end function fks_j_d


  integer function extra_cnt_d(configuration)
    integer, intent(in) :: configuration
    call check_configuration(configuration)
    extra_cnt_d = extra_cnt_values(configuration)
  end function extra_cnt_d


  integer function isplitorder_born_d(configuration)
    integer, intent(in) :: configuration
    call check_configuration(configuration)
    isplitorder_born_d = isplitorder_born_values(configuration)
  end function isplitorder_born_d


  integer function isplitorder_cnt_d(configuration)
    integer, intent(in) :: configuration
    call check_configuration(configuration)
    isplitorder_cnt_d = isplitorder_cnt_values(configuration)
  end function isplitorder_cnt_d


  integer function fks_j_from_i_d(configuration, particle, position)
    integer, intent(in) :: configuration, particle, position
    call check_configuration(configuration)
    call check_particle(particle)
    if (position < 0 .or. position > nexternal) then
      call fail_validation('FKS_J_FROM_I_D position is out of range')
    end if
    fks_j_from_i_d = fks_j_from_i_values(configuration, particle, position)
  end function fks_j_from_i_d


  integer function particle_type_d(configuration, particle)
    integer, intent(in) :: configuration, particle
    call check_configuration(configuration)
    call check_particle(particle)
    particle_type_d = particle_type_values(configuration, particle)
  end function particle_type_d


  integer function pdg_type_d(configuration, particle)
    integer, intent(in) :: configuration, particle
    call check_configuration(configuration)
    call check_particle(particle)
    pdg_type_d = pdg_type_values(configuration, particle)
  end function pdg_type_d


  logical function split_type_d(configuration, split_order)
    integer, intent(in) :: configuration, split_order
    call check_configuration(configuration)
    if (split_order < 1 .or. split_order > nsplitorders) then
      call fail_validation('SPLIT_TYPE_D order is out of range')
    end if
    split_type_d = split_type_values(configuration, split_order)
  end function split_type_d


  logical function need_color_links_d(configuration)
    integer, intent(in) :: configuration
    call check_configuration(configuration)
    need_color_links_d = need_color_links_values(configuration)
  end function need_color_links_d


  subroutine validate_input_shapes(fks_i_in, fks_j_in, extra_cnt_in, &
       isplitorder_born_in, isplitorder_cnt_in, fks_j_from_i_in, &
       particle_type_in, pdg_type_in, split_type_in, &
       need_color_links_in)
    integer, intent(in) :: fks_i_in(:), fks_j_in(:), extra_cnt_in(:)
    integer, intent(in) :: isplitorder_born_in(:), isplitorder_cnt_in(:)
    integer, intent(in) :: fks_j_from_i_in(:, :, 0:)
    integer, intent(in) :: particle_type_in(:, :), pdg_type_in(:, :)
    logical, intent(in) :: split_type_in(:, :)
    logical, intent(in) :: need_color_links_in(:)

    if (size(fks_i_in) /= fks_configs .or. &
        size(fks_j_in) /= fks_configs .or. &
        size(extra_cnt_in) /= fks_configs .or. &
        size(isplitorder_born_in) /= fks_configs .or. &
        size(isplitorder_cnt_in) /= fks_configs .or. &
        size(need_color_links_in) /= fks_configs) then
      call fail_validation('a configuration table has the wrong extent')
    end if
    if (size(fks_j_from_i_in, 1) /= fks_configs .or. &
        size(fks_j_from_i_in, 2) /= nexternal .or. &
        size(fks_j_from_i_in, 3) /= nexternal + 1 .or. &
        lbound(fks_j_from_i_in, 3) /= 0 .or. &
        ubound(fks_j_from_i_in, 3) /= nexternal) then
      call fail_validation('FKS_J_FROM_I_D has the wrong shape')
    end if
    if (size(particle_type_in, 1) /= fks_configs .or. &
        size(particle_type_in, 2) /= nexternal .or. &
        size(pdg_type_in, 1) /= fks_configs .or. &
        size(pdg_type_in, 2) /= nexternal) then
      call fail_validation('a particle metadata table has the wrong shape')
    end if
    if (size(split_type_in, 1) /= fks_configs .or. &
        size(split_type_in, 2) /= nsplitorders) then
      call fail_validation('SPLIT_TYPE_D has the wrong shape')
    end if
  end subroutine validate_input_shapes


  subroutine validate_input_values(fks_i_in, fks_j_in, extra_cnt_in, &
       isplitorder_born_in, isplitorder_cnt_in, fks_j_from_i_in, &
       split_type_in)
    integer, intent(in) :: fks_i_in(:), fks_j_in(:), extra_cnt_in(:)
    integer, intent(in) :: isplitorder_born_in(:), isplitorder_cnt_in(:)
    integer, intent(in) :: fks_j_from_i_in(:, :, 0:)
    logical, intent(in) :: split_type_in(:, :)
    integer :: configuration, emitter, position, count
    logical :: found_fks_partner

    if (any(fks_i_in < 1) .or. any(fks_i_in > nexternal) .or. &
        any(fks_j_in < 1) .or. any(fks_j_in > nexternal)) then
      call fail_validation('an FKS particle index is out of range')
    end if
    if (any(extra_cnt_in < 0)) then
      call fail_validation('EXTRA_CNT_D contains a negative identifier')
    end if
    if (any(isplitorder_born_in < 0) .or. &
        any(isplitorder_born_in > nsplitorders) .or. &
        any(isplitorder_cnt_in < 0) .or. &
        any(isplitorder_cnt_in > nsplitorders)) then
      call fail_validation('an extra-counterterm split order is invalid')
    end if

    do configuration = 1, fks_configs
      if (.not. any(split_type_in(configuration, :))) then
        call fail_validation('an FKS configuration has no splitting type')
      end if
      emitter = fks_i_in(configuration)
      count = fks_j_from_i_in(configuration, emitter, 0)
      if (count < 1 .or. count > nexternal) then
        call fail_validation('FKS_J_FROM_I_D has an invalid emitter count')
      end if
      found_fks_partner = .false.
      do position = 1, count
        if (fks_j_from_i_in(configuration, emitter, position) < 1 .or. &
            fks_j_from_i_in(configuration, emitter, position) > &
            nexternal) then
          call fail_validation('FKS_J_FROM_I_D has an invalid partner')
        end if
        if (fks_j_from_i_in(configuration, emitter, position) == &
            fks_j_in(configuration)) found_fks_partner = .true.
      end do
      if (.not. found_fks_partner) then
        call fail_validation('FKS_J_D is absent from FKS_J_FROM_I_D')
      end if
    end do
  end subroutine validate_input_values


  logical function storage_shapes_are_valid()
    storage_shapes_are_valid = .false.
    if (.not. allocated(fks_i_values)) return
    if (.not. allocated(fks_j_values)) return
    if (.not. allocated(extra_cnt_values)) return
    if (.not. allocated(isplitorder_born_values)) return
    if (.not. allocated(isplitorder_cnt_values)) return
    if (.not. allocated(fks_j_from_i_values)) return
    if (.not. allocated(particle_type_values)) return
    if (.not. allocated(pdg_type_values)) return
    if (.not. allocated(split_type_values)) return
    if (.not. allocated(need_color_links_values)) return
    if (size(fks_i_values) /= fks_configs) return
    if (size(fks_j_values) /= fks_configs) return
    if (size(extra_cnt_values) /= fks_configs) return
    if (size(isplitorder_born_values) /= fks_configs) return
    if (size(isplitorder_cnt_values) /= fks_configs) return
    if (size(fks_j_from_i_values, 1) /= fks_configs) return
    if (size(fks_j_from_i_values, 2) /= nexternal) return
    if (lbound(fks_j_from_i_values, 3) /= 0) return
    if (ubound(fks_j_from_i_values, 3) /= nexternal) return
    if (size(particle_type_values, 1) /= fks_configs) return
    if (size(particle_type_values, 2) /= nexternal) return
    if (size(pdg_type_values, 1) /= fks_configs) return
    if (size(pdg_type_values, 2) /= nexternal) return
    if (size(split_type_values, 1) /= fks_configs) return
    if (size(split_type_values, 2) /= nsplitorders) return
    if (size(need_color_links_values) /= fks_configs) return
    storage_shapes_are_valid = .true.
  end function storage_shapes_are_valid


  logical function same_metadata(fks_i_in, fks_j_in, extra_cnt_in, &
       isplitorder_born_in, isplitorder_cnt_in, fks_j_from_i_in, &
       particle_type_in, pdg_type_in, split_type_in, &
       need_color_links_in)
    integer, intent(in) :: fks_i_in(:), fks_j_in(:), extra_cnt_in(:)
    integer, intent(in) :: isplitorder_born_in(:), isplitorder_cnt_in(:)
    integer, intent(in) :: fks_j_from_i_in(:, :, 0:)
    integer, intent(in) :: particle_type_in(:, :), pdg_type_in(:, :)
    logical, intent(in) :: split_type_in(:, :)
    logical, intent(in) :: need_color_links_in(:)
    integer :: configuration, emitter, count

    same_metadata = .false.
    if (.not. storage_shapes_are_valid()) return
    if (.not. all(fks_i_values == fks_i_in)) return
    if (.not. all(fks_j_values == fks_j_in)) return
    if (.not. all(extra_cnt_values == extra_cnt_in)) return
    if (.not. all(isplitorder_born_values == isplitorder_born_in)) return
    if (.not. all(isplitorder_cnt_values == isplitorder_cnt_in)) return
    if (.not. all(particle_type_values == particle_type_in)) return
    if (.not. all(pdg_type_values == pdg_type_in)) return
    if (.not. all(split_type_values .eqv. split_type_in)) return
    if (.not. all(need_color_links_values .eqv. need_color_links_in)) return
    do configuration = 1, fks_configs
      emitter = fks_i_in(configuration)
      count = fks_j_from_i_in(configuration, emitter, 0)
      if (.not. all(fks_j_from_i_values(configuration, emitter, &
           0:count) == fks_j_from_i_in(configuration, emitter, &
           0:count))) return
    end do
    same_metadata = .true.
  end function same_metadata


  subroutine check_configuration(configuration)
    integer, intent(in) :: configuration
    if (.not. metadata_initialized) then
      call fail_validation('metadata are not initialized')
    end if
    if (configuration < 1 .or. configuration > fks_configs) then
      call fail_validation('FKS configuration is out of range')
    end if
  end subroutine check_configuration


  subroutine check_particle(particle)
    integer, intent(in) :: particle
    if (particle < 1 .or. particle > nexternal) then
      call fail_validation('particle index is out of range')
    end if
  end subroutine check_particle


  subroutine fail_validation(message)
    character(len=*), intent(in) :: message
    write (*, *) 'ERROR in fks_metadata: ', message
    stop 1
  end subroutine fail_validation

end module fks_metadata

module decay_chain_parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, node_pdg
  use nlo_decay_metadata, only: has_nlo_decay, corrected_parent_pdg
  implicit none
  private

  logical, save :: initialized = .false.
  double precision, save :: dummy_width_ratio_value = 0d0
  logical, save :: use_decayed_production_momenta_value = .false.
  integer, allocatable, save :: width_pdgs(:)
  double precision, allocatable, save :: width_values(:)
  logical, allocatable, save :: width_is_nlo(:)
  integer, allocatable, save :: scale_pdgs(:)
  double precision, allocatable, save :: scale_values(:)

  public :: initialize_decay_chain_parameters
  public :: decay_dummy_width_ratio, decay_physical_width
  public :: decay_renormalization_scale
  public :: use_decayed_production_ren_scale_momenta

contains

  subroutine initialize_decay_chain_parameters()
    logical :: exists, end_seen, format_seen, ratio_seen, momentum_mode_seen
    integer :: unit_number, ios, card_format, width_count, width_index
    integer :: scale_count, scale_index
    integer :: pdg, node, previous
    double precision :: value
    character(len=512) :: line
    character(len=32) :: keyword, momentum_mode

    if (initialized) return
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) then
      call fail_parameters('decay metadata are absent')
    end if

    inquire(file='decay_card.dat', exist=exists)
    if (.not. exists) call fail_parameters('decay_card.dat is absent')
    open(newunit=unit_number, file='decay_card.dat', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) call fail_parameters('cannot open decay_card.dat')

    width_count = 0
    scale_count = 0
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_parameters('cannot read decay_card.dat')
      if (skip_line(line)) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_parameters('malformed decay-card record')
      if (trim(keyword) == 'DECAY_WIDTH' .or. &
          trim(keyword) == 'NLO_DECAY_WIDTH') then
        width_count = width_count + 1
      end if
      if (trim(keyword) == 'DECAY_REN_SCALE') scale_count = scale_count + 1
    end do
    if (width_count < 1) then
      call fail_parameters('no physical-width records are present')
    end if
    if (scale_count < 1) then
      call fail_parameters('no DECAY_REN_SCALE records are present')
    end if
    allocate(width_pdgs(width_count))
    allocate(width_values(width_count))
    allocate(width_is_nlo(width_count))
    allocate(scale_pdgs(scale_count))
    allocate(scale_values(scale_count))
    width_pdgs = 0
    width_values = 0d0
    width_is_nlo = .false.
    scale_pdgs = 0
    scale_values = 0d0

    rewind(unit_number)
    card_format = 0
    width_index = 0
    scale_index = 0
    dummy_width_ratio_value = 0d0
    use_decayed_production_momenta_value = .false.
    end_seen = .false.
    format_seen = .false.
    ratio_seen = .false.
    momentum_mode_seen = .false.
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_parameters('cannot read decay-card body')
      if (skip_line(line)) cycle
      if (end_seen) call fail_parameters('record found after END')
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_parameters('malformed decay-card keyword')
      select case (trim(keyword))
      case ('FORMAT')
        if (format_seen) call fail_parameters('duplicate FORMAT record')
        read(line, *, iostat=ios) keyword, card_format
        format_seen = .true.
      case ('DUMMY_WIDTH_RATIO')
        if (ratio_seen) then
          call fail_parameters('duplicate DUMMY_WIDTH_RATIO record')
        end if
        read(line, *, iostat=ios) keyword, dummy_width_ratio_value
        ratio_seen = .true.
      case ('PRODUCTION_REN_SCALE_MOMENTA')
        if (momentum_mode_seen) then
          call fail_parameters(&
               'duplicate PRODUCTION_REN_SCALE_MOMENTA record')
        end if
        read(line, *, iostat=ios) keyword, momentum_mode
        if (ios == 0) then
          select case (trim(momentum_mode))
          case ('CORE')
            use_decayed_production_momenta_value = .false.
          case ('DECAYED')
            use_decayed_production_momenta_value = .true.
          case default
            call fail_parameters(&
                 'production scale momenta must be CORE or DECAYED')
          end select
        end if
        momentum_mode_seen = .true.
      case ('DECAY_WIDTH', 'NLO_DECAY_WIDTH')
        read(line, *, iostat=ios) keyword, pdg, value
        if (ios == 0) then
          pdg = abs(pdg)
          if (pdg == 0) call fail_parameters('a width has zero PDG code')
          width_index = width_index + 1
          do previous = 1, width_index - 1
            if (width_pdgs(previous) == pdg) then
              call fail_parameters('duplicate physical-width record')
            end if
          end do
          width_pdgs(width_index) = pdg
          width_values(width_index) = value
          width_is_nlo(width_index) = trim(keyword) == 'NLO_DECAY_WIDTH'
        end if
      case ('DECAY_REN_SCALE')
        read(line, *, iostat=ios) keyword, pdg, value
        if (ios == 0) then
          pdg = abs(pdg)
          if (pdg == 0) then
            call fail_parameters('a decay scale has zero PDG code')
          end if
          scale_index = scale_index + 1
          do previous = 1, scale_index - 1
            if (scale_pdgs(previous) == pdg) then
              call fail_parameters('duplicate DECAY_REN_SCALE record')
            end if
          end do
          scale_pdgs(scale_index) = pdg
          scale_values(scale_index) = value
        end if
      case ('END')
        end_seen = .true.
      case default
        call fail_parameters('unknown keyword '//trim(keyword))
      end select
      if (ios /= 0) call fail_parameters('malformed decay-card record')
    end do
    close(unit_number)

    if (.not. end_seen) call fail_parameters('END record is absent')
    if (card_format /= 3) call fail_parameters('FORMAT 3 is required')
    if (.not. ratio_seen .or. dummy_width_ratio_value <= 0d0 .or. &
        .not. ieee_is_finite(dummy_width_ratio_value)) then
      call fail_parameters('the dummy-width ratio must be finite and positive')
    end if
    if (.not. momentum_mode_seen) then
      call fail_parameters('PRODUCTION_REN_SCALE_MOMENTA record is absent')
    end if
    do width_index = 1, size(width_values)
      if (width_values(width_index) <= 0d0 .or. &
          .not. ieee_is_finite(width_values(width_index))) then
        call fail_parameters('physical widths must be finite and positive')
      end if
    end do
    do scale_index = 1, size(scale_values)
      if (scale_values(scale_index) <= 0d0 .or. &
          .not. ieee_is_finite(scale_values(scale_index))) then
        call fail_parameters('decay scales must be finite and positive')
      end if
    end do
    if (has_decay_chains()) then
      do node = 1, decay_node_count()
        width_index = find_pdg(node_pdg(node), width_pdgs)
        if (width_index == 0) then
          call fail_parameters('a decay node has no physical width')
        end if
        if (width_is_nlo(width_index)) then
          call fail_parameters('an LO decay node has an NLO width record')
        end if
        if (find_pdg(node_pdg(node), scale_pdgs) == 0) then
          call fail_parameters('a decay node has no renormalisation scale')
        end if
      end do
    else if (has_nlo_decay()) then
      width_index = find_pdg(corrected_parent_pdg(), width_pdgs)
      if (width_index == 0) then
        call fail_parameters('the corrected decay has no physical width')
      end if
      if (.not. width_is_nlo(width_index)) then
        call fail_parameters(&
             'the corrected decay requires an NLO_DECAY_WIDTH record')
      end if
      if (find_pdg(corrected_parent_pdg(), scale_pdgs) == 0) then
        call fail_parameters('the corrected decay has no renormalisation scale')
      end if
    end if
    initialized = .true.
  end subroutine initialize_decay_chain_parameters


  double precision function decay_dummy_width_ratio()
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_dummy_width_ratio = dummy_width_ratio_value
  end function decay_dummy_width_ratio


  double precision function decay_physical_width(pdg)
    integer, intent(in) :: pdg
    integer :: width_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    width_index = find_pdg(pdg, width_pdgs)
    if (width_index == 0) then
      call fail_parameters('requested particle has no physical width')
    end if
    decay_physical_width = width_values(width_index)
  end function decay_physical_width


  double precision function decay_renormalization_scale(pdg)
    integer, intent(in) :: pdg
    integer :: scale_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    scale_index = find_pdg(pdg, scale_pdgs)
    if (scale_index == 0) then
      call fail_parameters('requested particle has no decay scale')
    end if
    decay_renormalization_scale = scale_values(scale_index)
  end function decay_renormalization_scale


  logical function use_decayed_production_ren_scale_momenta()
    if (.not. initialized) call initialize_decay_chain_parameters()
    use_decayed_production_ren_scale_momenta = &
         use_decayed_production_momenta_value
  end function use_decayed_production_ren_scale_momenta


  integer function find_pdg(pdg, pdgs)
    integer, intent(in) :: pdg, pdgs(:)
    integer :: index

    find_pdg = 0
    do index = 1, size(pdgs)
      if (pdgs(index) == abs(pdg)) then
        find_pdg = index
        return
      end if
    end do
  end function find_pdg


  logical function skip_line(line)
    character(len=*), intent(in) :: line
    character(len=len(line)) :: adjusted

    adjusted = adjustl(line)
    skip_line = len_trim(adjusted) == 0
    if (.not. skip_line) then
      skip_line = adjusted(1:1) == '#' .or. adjusted(1:1) == '!'
    end if
  end function skip_line


  subroutine fail_parameters(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in decay_chain_parameters: '//trim(message)
    stop 1
  end subroutine fail_parameters

end module decay_chain_parameters

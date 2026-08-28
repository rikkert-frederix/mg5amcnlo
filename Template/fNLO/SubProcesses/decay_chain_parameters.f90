module decay_chain_parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, node_pdg
  use nlo_decay_metadata, only: has_nlo_decay, corrected_parent_pdg, &
       nlo_decay_node_count, nlo_decay_node_pdg
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       bundle_species_is_nlo
  implicit none
  private

  logical, save :: initialized = .false.
  double precision, save :: dummy_width_ratio_value = 0d0
  logical, save :: use_decayed_production_momenta_value = .false.
  integer, save :: number_of_width_species = 0
  integer, allocatable, save :: width_pdgs(:)
  double precision, allocatable, save :: lo_width_values(:)
  double precision, allocatable, save :: nlo_width_values(:)
  logical, allocatable, save :: has_lo_width(:)
  logical, allocatable, save :: has_nlo_width(:)
  integer, allocatable, save :: scale_pdgs(:)
  double precision, allocatable, save :: scale_values(:)

  public :: initialize_decay_chain_parameters
  public :: decay_dummy_width_ratio, decay_physical_width
  public :: decay_lo_width, decay_nlo_width
  public :: decay_width_expansion_coefficient
  public :: decay_renormalization_scale
  public :: use_decayed_production_ren_scale_momenta

contains

  subroutine initialize_decay_chain_parameters()
    logical :: exists, end_seen, format_seen, ratio_seen, momentum_mode_seen
    logical :: legacy_width_seen, explicit_lo_width_seen
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
          trim(keyword) == 'LO_DECAY_WIDTH' .or. &
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
    allocate(lo_width_values(width_count))
    allocate(nlo_width_values(width_count))
    allocate(has_lo_width(width_count))
    allocate(has_nlo_width(width_count))
    allocate(scale_pdgs(scale_count))
    allocate(scale_values(scale_count))
    width_pdgs = 0
    lo_width_values = 0d0
    nlo_width_values = 0d0
    has_lo_width = .false.
    has_nlo_width = .false.
    scale_pdgs = 0
    scale_values = 0d0

    rewind(unit_number)
    card_format = 0
    number_of_width_species = 0
    scale_index = 0
    dummy_width_ratio_value = 0d0
    use_decayed_production_momenta_value = .false.
    end_seen = .false.
    format_seen = .false.
    ratio_seen = .false.
    momentum_mode_seen = .false.
    legacy_width_seen = .false.
    explicit_lo_width_seen = .false.
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
      case ('DECAY_WIDTH', 'LO_DECAY_WIDTH', 'NLO_DECAY_WIDTH')
        read(line, *, iostat=ios) keyword, pdg, value
        if (ios == 0) then
          if (trim(keyword) == 'DECAY_WIDTH') legacy_width_seen = .true.
          if (trim(keyword) == 'LO_DECAY_WIDTH') then
            explicit_lo_width_seen = .true.
          end if
          pdg = abs(pdg)
          if (pdg == 0) call fail_parameters('a width has zero PDG code')
          width_index = find_pdg(pdg, width_pdgs)
          if (width_index == 0) then
            number_of_width_species = number_of_width_species + 1
            width_index = number_of_width_species
            width_pdgs(width_index) = pdg
          end if
          if (trim(keyword) == 'NLO_DECAY_WIDTH') then
            if (has_nlo_width(width_index)) then
              call fail_parameters('duplicate NLO physical-width record')
            end if
            nlo_width_values(width_index) = value
            has_nlo_width(width_index) = .true.
          else
            if (has_lo_width(width_index)) then
              call fail_parameters('duplicate LO physical-width record')
            end if
            lo_width_values(width_index) = value
            has_lo_width(width_index) = .true.
          end if
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
    if (card_format /= 3 .and. card_format /= 4) then
      call fail_parameters('FORMAT 3 or FORMAT 4 is required')
    end if
    if (card_format == 3 .and. explicit_lo_width_seen) then
      call fail_parameters(&
           'FORMAT 3 requires legacy DECAY_WIDTH records')
    end if
    if (card_format == 4 .and. legacy_width_seen) then
      call fail_parameters(&
           'FORMAT 4 requires explicit LO_DECAY_WIDTH records')
    end if
    if (has_nlo_contribution_bundle() .and. card_format /= 4) then
      call fail_parameters(&
           'a full NLO contribution bundle requires FORMAT 4 with both LO and NLO widths')
    end if
    if (.not. ratio_seen .or. dummy_width_ratio_value <= 0d0 .or. &
        .not. ieee_is_finite(dummy_width_ratio_value)) then
      call fail_parameters('the dummy-width ratio must be finite and positive')
    end if
    if (.not. momentum_mode_seen) then
      call fail_parameters('PRODUCTION_REN_SCALE_MOMENTA record is absent')
    end if
    do width_index = 1, number_of_width_species
      if (has_lo_width(width_index)) then
        if (lo_width_values(width_index) <= 0d0 .or. &
            .not. ieee_is_finite(lo_width_values(width_index))) then
          call fail_parameters('LO physical widths must be finite and positive')
        end if
      end if
      if (has_nlo_width(width_index)) then
        if (nlo_width_values(width_index) <= 0d0 .or. &
            .not. ieee_is_finite(nlo_width_values(width_index))) then
          call fail_parameters('NLO physical widths must be finite and positive')
        end if
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
        if (.not. has_lo_width(width_index)) then
          call fail_parameters('a decay node has no LO physical width')
        end if
        if (has_nlo_contribution_bundle() .and. &
            bundle_species_is_nlo(node_pdg(node))) then
          if (.not. has_nlo_width(width_index)) then
            call fail_parameters(&
                 'a corrected decay species has no NLO physical width')
          end if
        else if (has_nlo_width(width_index)) then
          call fail_parameters('an uncorrected decay node has an NLO width record')
        end if
        if (find_pdg(node_pdg(node), scale_pdgs) == 0) then
          call fail_parameters('a decay node has no renormalisation scale')
        end if
      end do
    else if (has_nlo_decay()) then
      do node = 1, nlo_decay_node_count()
        pdg = nlo_decay_node_pdg(node)
        width_index = find_pdg(pdg, width_pdgs)
        if (width_index == 0) then
          call fail_parameters('an NLO-decay topology node has no width')
        end if
        if ((has_nlo_contribution_bundle() .and. &
             bundle_species_is_nlo(pdg)) .or. &
            (.not. has_nlo_contribution_bundle() .and. &
             abs(pdg) == abs(corrected_parent_pdg()))) then
          if (.not. has_nlo_width(width_index)) then
            call fail_parameters(&
                 'the corrected species requires an NLO_DECAY_WIDTH record')
          end if
          if (has_nlo_contribution_bundle() .and. &
              .not. has_lo_width(width_index)) then
            call fail_parameters(&
                 'the corrected species requires an LO_DECAY_WIDTH record')
          end if
        else if (has_nlo_width(width_index)) then
          call fail_parameters('an uncorrected node has an NLO width record')
        else if (.not. has_lo_width(width_index)) then
          call fail_parameters('an uncorrected node has no LO width record')
        end if
        if (find_pdg(pdg, scale_pdgs) == 0) then
          call fail_parameters('an NLO-decay node has no renormalisation scale')
        end if
      end do
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
    if (has_nlo_contribution_bundle()) then
      decay_physical_width = decay_lo_width(pdg)
    else if (has_nlo_width(width_index) .and. &
             .not. has_lo_width(width_index)) then
      ! Compatibility with FORMAT 3 standalone NLO-decay cards.
      decay_physical_width = nlo_width_values(width_index)
    else if (has_nlo_decay() .and. &
             abs(pdg) == abs(corrected_parent_pdg()) .and. &
             has_nlo_width(width_index)) then
      decay_physical_width = nlo_width_values(width_index)
    else
      decay_physical_width = decay_lo_width(pdg)
    end if
  end function decay_physical_width


  double precision function decay_lo_width(pdg)
    integer, intent(in) :: pdg
    integer :: width_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    width_index = find_pdg(pdg, width_pdgs)
    if (width_index == 0 .or. .not. has_lo_width(width_index)) then
      call fail_parameters('requested particle has no LO physical width')
    end if
    decay_lo_width = lo_width_values(width_index)
  end function decay_lo_width


  double precision function decay_nlo_width(pdg)
    integer, intent(in) :: pdg
    integer :: width_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    width_index = find_pdg(pdg, width_pdgs)
    if (width_index == 0 .or. .not. has_nlo_width(width_index)) then
      call fail_parameters('requested particle has no NLO physical width')
    end if
    decay_nlo_width = nlo_width_values(width_index)
  end function decay_nlo_width


  double precision function decay_width_expansion_coefficient()
    integer :: node, pdg

    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_width_expansion_coefficient = 0d0
    if (.not. has_nlo_contribution_bundle() .or. &
        .not. has_decay_chains()) return
    ! Expand every physical 1/Gamma denominator whose species is corrected.
    ! This is deliberately a decay-node count, not an NLO-contribution count:
    ! an identical resonance assigned a QCD-inert numerator decay still has
    ! the same QCD-corrected total width in its NWA denominator.
    do node = 1, decay_node_count()
      pdg = node_pdg(node)
      if (.not. bundle_species_is_nlo(pdg)) cycle
      decay_width_expansion_coefficient = &
           decay_width_expansion_coefficient - &
           (decay_nlo_width(pdg) - decay_lo_width(pdg))/ &
           decay_lo_width(pdg)
    end do
  end function decay_width_expansion_coefficient


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

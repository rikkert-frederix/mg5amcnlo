module decay_chain_parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, &
       node_pdg, node_qcd_order
  use nlo_decay_metadata, only: has_nlo_decay, corrected_parent_pdg, &
       nlo_decay_node_count, nlo_decay_node_pdg, nlo_decay_node_qcd_order
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       bundle_species_is_nlo
  implicit none
  private

  logical, save :: initialized = .false.
  integer, parameter, public :: decay_scale_none = 0
  integer, parameter, public :: decay_scale_correlated = 1
  integer, parameter, public :: decay_scale_independent = 2
  integer, parameter, public :: nlo_combination_additive = 0
  integer, parameter, public :: nlo_combination_multiplicative = 1
  integer, save :: nlo_combination_mode_value = nlo_combination_additive
  double precision, save :: multiplicative_virtual_fraction_value = 1d0
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
  integer, save :: scale_variation_mode_value = decay_scale_none
  integer, save :: number_of_scale_factors = 1
  integer, save :: number_of_scale_species = 0
  double precision, allocatable, save :: scale_factor_values(:)
  integer, allocatable, save :: scale_species_pdgs(:)
  integer, save :: number_of_lo_width_variations = 0
  integer, save :: number_of_nlo_width_variations = 0
  integer, allocatable, save :: lo_variation_pdgs(:)
  integer, allocatable, save :: nlo_variation_pdgs(:)
  double precision, allocatable, save :: lo_variation_factors(:)
  double precision, allocatable, save :: nlo_variation_factors(:)
  double precision, allocatable, save :: lo_variation_values(:)
  double precision, allocatable, save :: nlo_variation_values(:)

  public :: initialize_decay_chain_parameters
  public :: decay_dummy_width_ratio, decay_physical_width
  public :: decay_lo_width, decay_nlo_width
  public :: decay_width_expansion_coefficient
  public :: decay_width_denominator_rescaling
  public :: decay_renormalization_scale
  public :: use_decayed_production_ren_scale_momenta
  public :: decay_scale_variation_mode, decay_scale_variation_enabled
  public :: decay_scale_factor_count, decay_scale_factor
  public :: decay_scale_species_count, decay_scale_species
  public :: decay_scale_species_index
  public :: nlo_combination_mode, uses_multiplicative_nlo_combination
  public :: multiplicative_virtual_sampling_fraction

contains

  subroutine initialize_decay_chain_parameters()
    logical :: exists, end_seen, format_seen, ratio_seen, momentum_mode_seen
    logical :: legacy_width_seen, explicit_lo_width_seen
    logical :: variation_mode_seen, scale_factors_seen
    logical :: combination_mode_seen, virtual_sampling_seen
    integer :: unit_number, ios, card_format, width_count, width_index
    integer :: scale_count, scale_index, factor_count, factor_index
    integer :: lo_variation_count, nlo_variation_count
    integer :: pdg, node, previous, declared_count
    double precision :: value, factor
    character(len=512) :: line
    character(len=32) :: keyword, momentum_mode, variation_mode
    character(len=32) :: combination_mode

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
    factor_count = 0
    lo_variation_count = 0
    nlo_variation_count = 0
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
      if (trim(keyword) == 'DECAY_SCALE_FACTORS') then
        if (factor_count /= 0) then
          call fail_parameters('duplicate DECAY_SCALE_FACTORS record')
        end if
        read(line, *, iostat=ios) keyword, factor_count
        if (ios /= 0 .or. factor_count < 1) then
          call fail_parameters('malformed DECAY_SCALE_FACTORS record')
        end if
      end if
      if (trim(keyword) == 'LO_DECAY_WIDTH_VARIATION') then
        lo_variation_count = lo_variation_count + 1
      end if
      if (trim(keyword) == 'NLO_DECAY_WIDTH_VARIATION') then
        nlo_variation_count = nlo_variation_count + 1
      end if
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
    allocate(scale_factor_values(max(1, factor_count)))
    allocate(lo_variation_pdgs(max(1, lo_variation_count)))
    allocate(nlo_variation_pdgs(max(1, nlo_variation_count)))
    allocate(lo_variation_factors(max(1, lo_variation_count)))
    allocate(nlo_variation_factors(max(1, nlo_variation_count)))
    allocate(lo_variation_values(max(1, lo_variation_count)))
    allocate(nlo_variation_values(max(1, nlo_variation_count)))
    width_pdgs = 0
    lo_width_values = 0d0
    nlo_width_values = 0d0
    has_lo_width = .false.
    has_nlo_width = .false.
    scale_pdgs = 0
    scale_values = 0d0
    scale_factor_values = 0d0
    scale_factor_values(1) = 1d0
    lo_variation_pdgs = 0
    nlo_variation_pdgs = 0
    lo_variation_factors = 0d0
    nlo_variation_factors = 0d0
    lo_variation_values = 0d0
    nlo_variation_values = 0d0

    rewind(unit_number)
    card_format = 0
    number_of_width_species = 0
    scale_index = 0
    number_of_scale_factors = 1
    number_of_lo_width_variations = 0
    number_of_nlo_width_variations = 0
    scale_variation_mode_value = decay_scale_none
    nlo_combination_mode_value = nlo_combination_additive
    multiplicative_virtual_fraction_value = 1d0
    dummy_width_ratio_value = 0d0
    use_decayed_production_momenta_value = .false.
    end_seen = .false.
    format_seen = .false.
    ratio_seen = .false.
    momentum_mode_seen = .false.
    legacy_width_seen = .false.
    explicit_lo_width_seen = .false.
    variation_mode_seen = .false.
    scale_factors_seen = .false.
    combination_mode_seen = .false.
    virtual_sampling_seen = .false.
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
      case ('NLO_COMBINATION_MODE')
        if (combination_mode_seen) then
          call fail_parameters('duplicate NLO_COMBINATION_MODE record')
        end if
        read(line, *, iostat=ios) keyword, combination_mode
        if (ios == 0) then
          select case (trim(combination_mode))
          case ('ADDITIVE')
            nlo_combination_mode_value = nlo_combination_additive
          case ('MULTIPLICATIVE')
            nlo_combination_mode_value = nlo_combination_multiplicative
          case default
            call fail_parameters(&
                 'NLO combination mode must be ADDITIVE or MULTIPLICATIVE')
          end select
        end if
        combination_mode_seen = .true.
      case ('MULTIPLICATIVE_VIRTUAL_FRACTION')
        if (virtual_sampling_seen) then
          call fail_parameters( &
               'duplicate MULTIPLICATIVE_VIRTUAL_FRACTION record')
        end if
        read(line, *, iostat=ios) keyword, &
             multiplicative_virtual_fraction_value
        virtual_sampling_seen = .true.
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
      case ('DECAY_SCALE_VARIATION_MODE')
        if (variation_mode_seen) then
          call fail_parameters(&
               'duplicate DECAY_SCALE_VARIATION_MODE record')
        end if
        read(line, *, iostat=ios) keyword, variation_mode
        if (ios == 0) then
          select case (trim(variation_mode))
          case ('NONE')
            scale_variation_mode_value = decay_scale_none
          case ('CORRELATED')
            scale_variation_mode_value = decay_scale_correlated
          case ('INDEPENDENT')
            scale_variation_mode_value = decay_scale_independent
          case default
            call fail_parameters(&
                 'decay scale mode must be NONE, CORRELATED or INDEPENDENT')
          end select
        end if
        variation_mode_seen = .true.
      case ('DECAY_SCALE_FACTORS')
        if (scale_factors_seen) then
          call fail_parameters('duplicate DECAY_SCALE_FACTORS record')
        end if
        read(line, *, iostat=ios) keyword, declared_count
        if (ios /= 0 .or. declared_count /= factor_count) then
          call fail_parameters('malformed DECAY_SCALE_FACTORS record')
        end if
        read(line, *, iostat=ios) keyword, declared_count, &
             (scale_factor_values(factor_index), &
              factor_index=1, declared_count)
        if (ios == 0) number_of_scale_factors = declared_count
        scale_factors_seen = .true.
      case ('LO_DECAY_WIDTH_VARIATION')
        read(line, *, iostat=ios) keyword, pdg, factor, value
        if (ios == 0) then
          number_of_lo_width_variations = &
               number_of_lo_width_variations + 1
          lo_variation_pdgs(number_of_lo_width_variations) = abs(pdg)
          lo_variation_factors(number_of_lo_width_variations) = factor
          lo_variation_values(number_of_lo_width_variations) = value
        end if
      case ('NLO_DECAY_WIDTH_VARIATION')
        read(line, *, iostat=ios) keyword, pdg, factor, value
        if (ios == 0) then
          number_of_nlo_width_variations = &
               number_of_nlo_width_variations + 1
          nlo_variation_pdgs(number_of_nlo_width_variations) = abs(pdg)
          nlo_variation_factors(number_of_nlo_width_variations) = factor
          nlo_variation_values(number_of_nlo_width_variations) = value
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
    if (card_format /= 3 .and. card_format /= 4 .and. card_format /= 5) then
      call fail_parameters('FORMAT 3, FORMAT 4 or FORMAT 5 is required')
    end if
    if (card_format == 3 .and. explicit_lo_width_seen) then
      call fail_parameters(&
           'FORMAT 3 requires legacy DECAY_WIDTH records')
    end if
    if ((card_format == 4 .or. card_format == 5) .and. &
        legacy_width_seen) then
      call fail_parameters(&
           'FORMAT 4/5 requires explicit LO_DECAY_WIDTH records')
    end if
    if (has_nlo_contribution_bundle() .and. &
        card_format /= 4 .and. card_format /= 5) then
      call fail_parameters(&
           'a full NLO contribution bundle requires FORMAT 4/5 with both LO and NLO widths')
    end if
    if (.not. ratio_seen .or. dummy_width_ratio_value <= 0d0 .or. &
        .not. ieee_is_finite(dummy_width_ratio_value)) then
      call fail_parameters('the dummy-width ratio must be finite and positive')
    end if
    if (.not. momentum_mode_seen) then
      call fail_parameters('PRODUCTION_REN_SCALE_MOMENTA record is absent')
    end if
    if (nlo_combination_mode_value == nlo_combination_multiplicative .and. &
        .not. has_nlo_contribution_bundle()) then
      call fail_parameters(&
           'MULTIPLICATIVE mode requires a full NLO contribution bundle')
    end if
    if (.not. ieee_is_finite(multiplicative_virtual_fraction_value) .or. &
        multiplicative_virtual_fraction_value <= 0d0 .or. &
        multiplicative_virtual_fraction_value > 1d0) then
      call fail_parameters( &
           'MULTIPLICATIVE_VIRTUAL_FRACTION must be in (0,1]')
    end if
    if (multiplicative_virtual_fraction_value < 1d0 .and. &
        nlo_combination_mode_value /= nlo_combination_multiplicative) then
      call fail_parameters( &
           'stochastic virtual sampling requires MULTIPLICATIVE mode')
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
    call initialize_scale_species()
    call validate_scale_variations(card_format, variation_mode_seen, &
                                   scale_factors_seen)
    initialized = .true.
  end subroutine initialize_decay_chain_parameters


  double precision function decay_dummy_width_ratio()
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_dummy_width_ratio = dummy_width_ratio_value
  end function decay_dummy_width_ratio


  double precision function decay_physical_width(pdg, use_nlo_width)
    integer, intent(in) :: pdg
    logical, intent(in), optional :: use_nlo_width
    integer :: width_index
    logical :: select_nlo_width

    if (.not. initialized) call initialize_decay_chain_parameters()
    width_index = find_pdg(pdg, width_pdgs)
    if (width_index == 0) then
      call fail_parameters('requested particle has no physical width')
    end if
    if (has_nlo_contribution_bundle()) then
      if (uses_multiplicative_nlo_combination() .and. &
          bundle_species_is_nlo(pdg)) then
        decay_physical_width = decay_nlo_width(pdg)
      else
        decay_physical_width = decay_lo_width(pdg)
      end if
    else if (has_nlo_width(width_index) .and. &
             .not. has_lo_width(width_index)) then
      ! Compatibility with FORMAT 3 standalone NLO-decay cards.
      decay_physical_width = nlo_width_values(width_index)
    else
      select_nlo_width = .false.
      if (present(use_nlo_width)) select_nlo_width = use_nlo_width
      if (select_nlo_width .and. has_nlo_width(width_index)) then
        decay_physical_width = nlo_width_values(width_index)
      else
        decay_physical_width = decay_lo_width(pdg)
      end if
    end if
  end function decay_physical_width


  double precision function decay_lo_width(pdg, factor_index)
    integer, intent(in) :: pdg
    integer, intent(in), optional :: factor_index
    integer :: width_index, variation_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    width_index = find_pdg(pdg, width_pdgs)
    if (width_index == 0 .or. .not. has_lo_width(width_index)) then
      call fail_parameters('requested particle has no LO physical width')
    end if
    if (.not. present(factor_index) .or. factor_index == 1) then
      decay_lo_width = lo_width_values(width_index)
      return
    end if
    call validate_factor_index(factor_index)
    variation_index = find_width_variation(&
         pdg, scale_factor_values(factor_index), &
         number_of_lo_width_variations, lo_variation_pdgs, &
         lo_variation_factors)
    if (variation_index == 0) then
      call fail_parameters('requested LO width variation is absent')
    end if
    decay_lo_width = lo_variation_values(variation_index)
  end function decay_lo_width


  double precision function decay_nlo_width(pdg, factor_index)
    integer, intent(in) :: pdg
    integer, intent(in), optional :: factor_index
    integer :: width_index, variation_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    width_index = find_pdg(pdg, width_pdgs)
    if (width_index == 0 .or. .not. has_nlo_width(width_index)) then
      call fail_parameters('requested particle has no NLO physical width')
    end if
    if (.not. present(factor_index) .or. factor_index == 1) then
      decay_nlo_width = nlo_width_values(width_index)
      return
    end if
    call validate_factor_index(factor_index)
    variation_index = find_width_variation(&
         pdg, scale_factor_values(factor_index), &
         number_of_nlo_width_variations, nlo_variation_pdgs, &
         nlo_variation_factors)
    if (variation_index == 0) then
      call fail_parameters('requested NLO width variation is absent')
    end if
    decay_nlo_width = nlo_variation_values(variation_index)
  end function decay_nlo_width


  double precision function decay_width_expansion_coefficient(&
       factor_indices)
    integer, intent(in), optional :: factor_indices(:)
    integer :: node, pdg, factor_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_width_expansion_coefficient = 0d0
    if (.not. has_nlo_contribution_bundle() .or. &
        .not. has_decay_chains()) return
    if (uses_multiplicative_nlo_combination()) return
    ! Expand every physical 1/Gamma denominator whose species is corrected.
    ! This is deliberately a decay-node count, not an NLO-contribution count:
    ! an identical resonance assigned a QCD-inert numerator decay still has
    ! the same QCD-corrected total width in its NWA denominator.
    do node = 1, decay_node_count()
      pdg = node_pdg(node)
      if (.not. bundle_species_is_nlo(pdg)) cycle
      factor_index = selected_factor_index(pdg, factor_indices)
      decay_width_expansion_coefficient = &
           decay_width_expansion_coefficient - &
           (decay_nlo_width(pdg, factor_index) - &
            decay_lo_width(pdg, factor_index))/ &
           decay_lo_width(pdg, factor_index)
    end do
  end function decay_width_expansion_coefficient


  double precision function decay_width_denominator_rescaling(&
       factor_indices)
    integer, intent(in), optional :: factor_indices(:)
    integer :: node, pdg, factor_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_width_denominator_rescaling = 1d0
    if (.not. has_decay_chains()) return
    do node = 1, decay_node_count()
      pdg = node_pdg(node)
      factor_index = selected_factor_index(pdg, factor_indices)
      if (factor_index == 1) cycle
      if (uses_multiplicative_nlo_combination() .and. &
          bundle_species_is_nlo(pdg)) then
        decay_width_denominator_rescaling = &
             decay_width_denominator_rescaling*decay_nlo_width(pdg)/ &
             decay_nlo_width(pdg, factor_index)
      else
        decay_width_denominator_rescaling = &
             decay_width_denominator_rescaling*decay_lo_width(pdg)/ &
             decay_lo_width(pdg, factor_index)
      end if
    end do
  end function decay_width_denominator_rescaling


  integer function nlo_combination_mode()
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) then
      nlo_combination_mode = nlo_combination_additive
      return
    end if
    if (.not. initialized) call initialize_decay_chain_parameters()
    nlo_combination_mode = nlo_combination_mode_value
  end function nlo_combination_mode


  logical function uses_multiplicative_nlo_combination()
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) then
      uses_multiplicative_nlo_combination = .false.
      return
    end if
    if (.not. initialized) call initialize_decay_chain_parameters()
    uses_multiplicative_nlo_combination = &
         nlo_combination_mode_value == nlo_combination_multiplicative
  end function uses_multiplicative_nlo_combination


  double precision function multiplicative_virtual_sampling_fraction()
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) then
      multiplicative_virtual_sampling_fraction = 1d0
      return
    end if
    if (.not. initialized) call initialize_decay_chain_parameters()
    multiplicative_virtual_sampling_fraction = &
         multiplicative_virtual_fraction_value
  end function multiplicative_virtual_sampling_fraction


  double precision function decay_renormalization_scale(pdg, factor_index)
    integer, intent(in) :: pdg
    integer, intent(in), optional :: factor_index
    integer :: scale_index, selected_index

    if (.not. initialized) call initialize_decay_chain_parameters()
    scale_index = find_pdg(pdg, scale_pdgs)
    if (scale_index == 0) then
      call fail_parameters('requested particle has no decay scale')
    end if
    selected_index = 1
    if (present(factor_index)) selected_index = factor_index
    call validate_factor_index(selected_index)
    decay_renormalization_scale = scale_values(scale_index)* &
         scale_factor_values(selected_index)
  end function decay_renormalization_scale


  integer function decay_scale_variation_mode()
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_scale_variation_mode = scale_variation_mode_value
  end function decay_scale_variation_mode


  logical function decay_scale_variation_enabled()
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_scale_variation_enabled = &
         scale_variation_mode_value /= decay_scale_none .and. &
         number_of_scale_factors > 1 .and. number_of_scale_species > 0
  end function decay_scale_variation_enabled


  integer function decay_scale_factor_count()
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_scale_factor_count = number_of_scale_factors
  end function decay_scale_factor_count


  double precision function decay_scale_factor(index)
    integer, intent(in) :: index
    if (.not. initialized) call initialize_decay_chain_parameters()
    call validate_factor_index(index)
    decay_scale_factor = scale_factor_values(index)
  end function decay_scale_factor


  integer function decay_scale_species_count()
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_scale_species_count = number_of_scale_species
  end function decay_scale_species_count


  integer function decay_scale_species(index)
    integer, intent(in) :: index
    if (.not. initialized) call initialize_decay_chain_parameters()
    if (index < 1 .or. index > number_of_scale_species) then
      call fail_parameters('decay scale species index is out of range')
    end if
    decay_scale_species = scale_species_pdgs(index)
  end function decay_scale_species


  integer function decay_scale_species_index(pdg)
    integer, intent(in) :: pdg
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_scale_species_index = find_pdg(pdg, scale_species_pdgs)
  end function decay_scale_species_index


  logical function use_decayed_production_ren_scale_momenta()
    if (.not. initialized) call initialize_decay_chain_parameters()
    use_decayed_production_ren_scale_momenta = &
         use_decayed_production_momenta_value
  end function use_decayed_production_ren_scale_momenta


  subroutine initialize_scale_species()
    integer :: candidate_count, node, pdg, index, previous
    integer, allocatable :: candidates(:)
    logical :: scale_dependent

    if (has_decay_chains()) then
      candidate_count = decay_node_count()
    else
      candidate_count = nlo_decay_node_count()
    end if
    allocate(candidates(max(1, candidate_count)))
    candidates = 0
    number_of_scale_species = 0
    do node = 1, candidate_count
      if (has_decay_chains()) then
        pdg = abs(node_pdg(node))
        scale_dependent = node_qcd_order(node) > 0
        if (has_nlo_contribution_bundle()) then
          scale_dependent = scale_dependent .or. bundle_species_is_nlo(pdg)
        end if
      else
        pdg = abs(nlo_decay_node_pdg(node))
        scale_dependent = nlo_decay_node_qcd_order(node) > 0 .or. &
             pdg == abs(corrected_parent_pdg())
      end if
      if (.not. scale_dependent) cycle
      if (find_pdg(pdg, candidates) /= 0) cycle
      number_of_scale_species = number_of_scale_species + 1
      candidates(number_of_scale_species) = pdg
    end do
    allocate(scale_species_pdgs(max(1, number_of_scale_species)))
    scale_species_pdgs = 0
    if (number_of_scale_species > 0) then
      scale_species_pdgs(1:number_of_scale_species) = &
           candidates(1:number_of_scale_species)
      do index = 2, number_of_scale_species
        pdg = scale_species_pdgs(index)
        previous = index - 1
        do while (previous >= 1 .and. &
                  scale_species_pdgs(previous) > pdg)
          scale_species_pdgs(previous + 1) = &
               scale_species_pdgs(previous)
          previous = previous - 1
        end do
        scale_species_pdgs(previous + 1) = pdg
      end do
    end if
    deallocate(candidates)
  end subroutine initialize_scale_species


  subroutine validate_scale_variations(card_format, variation_mode_seen, &
                                       scale_factors_seen)
    integer, intent(in) :: card_format
    logical, intent(in) :: variation_mode_seen, scale_factors_seen
    integer :: index, factor_index, pdg

    if (card_format /= 5) then
      if (variation_mode_seen .or. scale_factors_seen .or. &
          number_of_lo_width_variations /= 0 .or. &
          number_of_nlo_width_variations /= 0) then
        call fail_parameters(&
             'decay-scale variation records require FORMAT 5')
      end if
      scale_variation_mode_value = decay_scale_none
      number_of_scale_factors = 1
      scale_factor_values(1) = 1d0
      return
    end if
    if (.not. has_nlo_contribution_bundle()) then
      call fail_parameters(&
           'FORMAT 5 scale variations require a full NLO bundle')
    end if
    if (.not. variation_mode_seen .or. .not. scale_factors_seen) then
      call fail_parameters(&
           'FORMAT 5 requires decay scale mode and factors')
    end if
    if (scale_variation_mode_value == decay_scale_none) then
      call fail_parameters('FORMAT 5 requires a non-NONE scale mode')
    end if
    if (number_of_scale_species < 1 .or. number_of_scale_factors < 2) then
      call fail_parameters(&
           'decay-scale variation requires species and noncentral factors')
    end if
    do factor_index = 1, number_of_scale_factors
      if (scale_factor_values(factor_index) <= 0d0 .or. &
          .not. ieee_is_finite(scale_factor_values(factor_index))) then
        call fail_parameters('decay scale factors must be finite and positive')
      end if
      do index = 1, factor_index - 1
        if (same_factor(scale_factor_values(index), &
                        scale_factor_values(factor_index))) then
          call fail_parameters('duplicate decay scale factor')
        end if
      end do
    end do
    if (.not. same_factor(scale_factor_values(1), 1d0)) then
      call fail_parameters('the first decay scale factor must be one')
    end if

    call validate_variation_records(&
         number_of_lo_width_variations, lo_variation_pdgs, &
         lo_variation_factors, lo_variation_values, .false.)
    call validate_variation_records(&
         number_of_nlo_width_variations, nlo_variation_pdgs, &
         nlo_variation_factors, nlo_variation_values, .true.)
    do index = 1, number_of_scale_species
      pdg = scale_species_pdgs(index)
      do factor_index = 2, number_of_scale_factors
        if (find_width_variation(&
             pdg, scale_factor_values(factor_index), &
             number_of_lo_width_variations, lo_variation_pdgs, &
             lo_variation_factors) == 0) then
          call fail_parameters(&
               'a varied decay species has no LO width at every factor')
        end if
        if (bundle_species_is_nlo(pdg) .and. &
            find_width_variation(&
             pdg, scale_factor_values(factor_index), &
             number_of_nlo_width_variations, nlo_variation_pdgs, &
             nlo_variation_factors) == 0) then
          call fail_parameters(&
               'a corrected decay species has no NLO width at every factor')
        end if
      end do
    end do
  end subroutine validate_scale_variations


  subroutine validate_variation_records(count, pdgs, factors, values, &
                                        require_nlo_species)
    integer, intent(in) :: count, pdgs(:)
    double precision, intent(in) :: factors(:), values(:)
    logical, intent(in) :: require_nlo_species
    integer :: index, previous, factor_index

    do index = 1, count
      if (find_pdg(pdgs(index), scale_species_pdgs) == 0) then
        call fail_parameters(&
             'a width variation refers to a non-varied decay species')
      end if
      if (require_nlo_species .and. &
          .not. bundle_species_is_nlo(pdgs(index))) then
        call fail_parameters(&
             'an NLO width variation refers to an uncorrected species')
      end if
      factor_index = find_factor_index(factors(index))
      if (factor_index < 2) then
        call fail_parameters(&
             'a width variation has an unknown or central scale factor')
      end if
      if (values(index) <= 0d0 .or. &
          .not. ieee_is_finite(values(index))) then
        call fail_parameters('varied widths must be finite and positive')
      end if
      do previous = 1, index - 1
        if (pdgs(previous) == pdgs(index) .and. &
            same_factor(factors(previous), factors(index))) then
          call fail_parameters('duplicate physical-width variation')
        end if
      end do
    end do
  end subroutine validate_variation_records


  integer function selected_factor_index(pdg, factor_indices)
    integer, intent(in) :: pdg
    integer, intent(in), optional :: factor_indices(:)
    integer :: species_index, index

    selected_factor_index = 1
    species_index = find_pdg(pdg, scale_species_pdgs)
    if (species_index == 0 .or. .not. present(factor_indices)) return
    if (size(factor_indices) == 0) return
    if (size(factor_indices) /= number_of_scale_species) then
      call fail_parameters('decay factor-index array has the wrong size')
    end if
    do index = 1, size(factor_indices)
      call validate_factor_index(factor_indices(index))
    end do
    selected_factor_index = factor_indices(species_index)
  end function selected_factor_index


  subroutine validate_factor_index(index)
    integer, intent(in) :: index
    if (index < 1 .or. index > number_of_scale_factors) then
      call fail_parameters('decay scale factor index is out of range')
    end if
  end subroutine validate_factor_index


  integer function find_factor_index(factor)
    double precision, intent(in) :: factor
    integer :: index
    find_factor_index = 0
    do index = 1, number_of_scale_factors
      if (same_factor(factor, scale_factor_values(index))) then
        find_factor_index = index
        return
      end if
    end do
  end function find_factor_index


  integer function find_width_variation(pdg, factor, count, pdgs, factors)
    integer, intent(in) :: pdg, count, pdgs(:)
    double precision, intent(in) :: factor, factors(:)
    integer :: index
    find_width_variation = 0
    do index = 1, count
      if (pdgs(index) == abs(pdg) .and. &
          same_factor(factors(index), factor)) then
        find_width_variation = index
        return
      end if
    end do
  end function find_width_variation


  logical function same_factor(first, second)
    double precision, intent(in) :: first, second
    same_factor = abs(first - second) <= &
         1d-12*max(1d0, abs(first), abs(second))
  end function same_factor


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

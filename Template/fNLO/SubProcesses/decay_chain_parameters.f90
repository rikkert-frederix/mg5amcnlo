module decay_chain_parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, &
       node_pdg, node_qcd_order
  use nlo_decay_metadata, only: has_nlo_decay, corrected_parent_pdg, &
       nlo_decay_node_count, nlo_decay_node_pdg, nlo_decay_node_qcd_order, &
       nlo_decay_corrected_node
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       bundle_species_is_nlo, contribution_is_nlo_decay, &
       contribution_parent_pdg, active_nlo_contribution
  use alfas_functions_module, only: alphas
  implicit none
  private

  logical, save :: initialized = .false.
  integer, parameter, public :: decay_scale_none = 0
  integer, parameter, public :: decay_scale_correlated = 1
  integer, parameter, public :: decay_scale_independent = 2
  integer, parameter, public :: nlo_decay_additive = 0
  integer, parameter, public :: nlo_decay_multiplicative = 1
  logical, save :: use_decayed_production_momenta_value = .false.
  logical, save :: production_w_system_value = .false.
  logical, save :: production_current_sampling_value = .false.
  double precision, save :: production_sampling_mass_value = 0d0
  double precision, save :: production_sampling_width_value = 0d0
  integer, save :: number_of_width_species = 0
  integer, allocatable, save :: width_pdgs(:)
  double precision, allocatable, save :: lo_width_values(:)
  double precision, allocatable, save :: nlo_width_values(:)
  logical, allocatable, save :: has_lo_width(:)
  logical, allocatable, save :: has_nlo_width(:)
  integer, allocatable, save :: scale_pdgs(:)
  double precision, allocatable, save :: scale_values(:)
  integer, allocatable, save :: dynamic_scale_choices(:)
  logical, allocatable, save :: automatic_widths(:), decay_orders(:)
  logical, save :: production_nlo = .true., decays_nlo = .true.
  logical, save :: lo_run = .false.
  integer, save :: scale_variation_mode_value = decay_scale_none
  ! Widths and reference-scale prescriptions remain absolute-PDG inputs.
  ! Only the reweighting axes optionally distinguish particle/antiparticle.
  logical, save :: signed_scale_axes = .false.
  integer, save :: nlo_combination_value = nlo_decay_additive
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
  public :: decay_physical_width
  public :: decay_lo_width, decay_nlo_width
  public :: decay_width_expansion_coefficient
  public :: decay_width_denominator_rescaling
  public :: decay_multiplicative_width_rescaling
  public :: nlo_decay_combination_mode, multiplicative_nlo_enabled
  public :: decay_renormalization_scale
  public :: use_decayed_production_ren_scale_momenta
  public :: production_w_system_scale
  public :: production_current_proposal
  public :: decay_scale_variation_mode, decay_scale_variation_enabled
  public :: decay_scale_factor_count, decay_scale_factor
  public :: decay_scale_species_count, decay_scale_species
  public :: decay_scale_species_index
  public :: decay_dynamical_scale_choice, decay_automatic_width
  public :: decay_species_nlo_enabled, nlo_correction_enabled
  public :: set_decay_run_order, decay_node_width_rescaling

contains

  subroutine initialize_decay_chain_parameters()
    logical :: exists, momentum_mode_seen, production_grouping_seen
    logical :: variation_mode_seen, scale_factors_seen, scale_grouping_seen
    logical :: combination_mode_seen, production_order_seen, decay_order_seen
    logical :: sampling_mode_seen, sampling_mass_seen, sampling_width_seen
    integer :: unit_number, ios, width_count, width_index
    integer :: scale_count, scale_index, factor_count, factor_index
    integer :: lo_variation_count, nlo_variation_count
    integer :: pdg, node, previous, declared_count
    double precision :: value, factor
    character(len=2048) :: line
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
      call normalize_decay_card_record(line)
      if (skip_line(line)) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_parameters('malformed decay-card record')
      if (trim(keyword) == 'LO_DECAY_WIDTH' .or. &
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
    allocate(dynamic_scale_choices(scale_count), automatic_widths(scale_count))
    allocate(decay_orders(scale_count))
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
    dynamic_scale_choices = 0
    automatic_widths = .false.
    decay_orders = .true.
    scale_factor_values = 0d0
    scale_factor_values(1) = 1d0
    lo_variation_pdgs = 0
    nlo_variation_pdgs = 0
    lo_variation_factors = 0d0
    nlo_variation_factors = 0d0
    lo_variation_values = 0d0
    nlo_variation_values = 0d0

    rewind(unit_number)
    number_of_width_species = 0
    scale_index = 0
    number_of_scale_factors = 1
    number_of_lo_width_variations = 0
    number_of_nlo_width_variations = 0
    scale_variation_mode_value = decay_scale_none
    signed_scale_axes = .false.
    nlo_combination_value = nlo_decay_additive
    use_decayed_production_momenta_value = .false.
    production_w_system_value = .false.
    production_grouping_seen = .false.
    production_current_sampling_value = .false.
    production_sampling_mass_value = 0d0
    production_sampling_width_value = 0d0
    sampling_mode_seen = .false.
    sampling_mass_seen = .false.
    sampling_width_seen = .false.
    momentum_mode_seen = .false.
    variation_mode_seen = .false.
    scale_grouping_seen = .false.
    scale_factors_seen = .false.
    combination_mode_seen = .false.
    production_order_seen = .false.
    decay_order_seen = .false.
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_parameters('cannot read decay-card body')
      call normalize_decay_card_record(line)
      if (skip_line(line)) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_parameters('malformed decay-card keyword')
      select case (trim(keyword))
      case ('PRODUCTION_PHASE_SPACE_SAMPLING')
        if (sampling_mode_seen) call fail_parameters('duplicate PRODUCTION_PHASE_SPACE_SAMPLING record')
        read(line, *, iostat=ios) keyword, momentum_mode
        if (ios == 0) then
          select case (trim(momentum_mode))
          case ('FLAT')
            production_current_sampling_value = .false.
          case ('W_CURRENT')
            production_current_sampling_value = .true.
          case default
            call fail_parameters('production phase-space sampling must be FLAT or W_CURRENT')
          end select
        end if
        sampling_mode_seen = .true.
      case ('PRODUCTION_SAMPLING_MASS')
        if (sampling_mass_seen) call fail_parameters('duplicate PRODUCTION_SAMPLING_MASS record')
        read(line, *, iostat=ios) keyword, production_sampling_mass_value
        sampling_mass_seen = .true.
      case ('PRODUCTION_SAMPLING_WIDTH')
        if (sampling_width_seen) call fail_parameters('duplicate PRODUCTION_SAMPLING_WIDTH record')
        read(line, *, iostat=ios) keyword, production_sampling_width_value
        sampling_width_seen = .true.
      case ('PRODUCTION_SCALE_GROUPING')
        if (production_grouping_seen) &
             call fail_parameters('duplicate PRODUCTION_SCALE_GROUPING record')
        read(line, *, iostat=ios) keyword, momentum_mode
        if (ios == 0) then
          select case (trim(momentum_mode))
          case ('NONE')
            production_w_system_value = .false.
          case ('W_SYSTEM')
            production_w_system_value = .true.
          case default
            call fail_parameters('production scale grouping must be NONE or W_SYSTEM')
          end select
        end if
        production_grouping_seen = .true.
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
      case ('LO_DECAY_WIDTH', 'NLO_DECAY_WIDTH')
        read(line, *, iostat=ios) keyword, pdg, value
        if (ios == 0) then
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
      case ('DECAY_SCALE_GROUPING')
        if (scale_grouping_seen) then
          call fail_parameters('duplicate DECAY_SCALE_GROUPING record')
        end if
        read(line, *, iostat=ios) keyword, variation_mode
        if (ios == 0) then
          select case (trim(variation_mode))
          case ('SPECIES')
            signed_scale_axes = .false.
          case ('SIGNED_PDG')
            signed_scale_axes = .true.
          case default
            call fail_parameters('decay scale grouping must be SPECIES or SIGNED_PDG')
          end select
        end if
        scale_grouping_seen = .true.
      case ('NLO_DECAY_COMBINATION')
        if (combination_mode_seen) then
          call fail_parameters('duplicate NLO_DECAY_COMBINATION record')
        end if
        read(line, *, iostat=ios) keyword, combination_mode
        if (ios == 0) then
          select case (trim(combination_mode))
          case ('ADDITIVE')
            nlo_combination_value = nlo_decay_additive
          case ('MULTIPLICATIVE')
            nlo_combination_value = nlo_decay_multiplicative
          case default
            call fail_parameters( &
                 'NLO decay combination must be ADDITIVE or MULTIPLICATIVE')
          end select
        end if
        combination_mode_seen = .true.
      case ('PRODUCTION_ORDER', 'DECAY_ORDER')
        read(line, *, iostat=ios) keyword, combination_mode
        if (ios == 0) then
          if (trim(combination_mode) /= 'LO' .and. &
              trim(combination_mode) /= 'NLO') then
            call fail_parameters('perturbative orders must be LO or NLO')
          end if
          if (trim(keyword) == 'PRODUCTION_ORDER') then
            if (production_order_seen) &
                 call fail_parameters('duplicate PRODUCTION_ORDER record')
            production_nlo = trim(combination_mode) == 'NLO'
            production_order_seen = .true.
          else
            if (decay_order_seen) &
                 call fail_parameters('duplicate DECAY_ORDER record')
            decays_nlo = trim(combination_mode) == 'NLO'
            decay_order_seen = .true.
          end if
        end if
      case ('DECAY_DYNAMICAL_SCALE_CHOICE', 'DECAY_WIDTH_SCALE_MODE', &
            'DECAY_PERTURBATIVE_ORDER')
        ! These indexed options are read in a separate pass so their order
        ! relative to DECAY_REN_SCALE is immaterial.
        continue
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
      case default
        call fail_parameters('unknown keyword '//trim(keyword))
      end select
      if (ios /= 0) call fail_parameters('malformed decay-card record')
    end do
    ! Keep the scale and width lookups on the same species domain. This is
    ! also required before indexed options and width variations are resolved.
    do scale_index = 1, size(scale_pdgs)
      if (find_pdg(scale_pdgs(scale_index), width_pdgs) == 0) &
           call fail_parameters('a decay scale has no physical width')
    end do
    do width_index = 1, number_of_width_species
      if (find_pdg(width_pdgs(width_index), scale_pdgs) == 0) &
           call fail_parameters('a physical width has no renormalisation scale')
    end do
    call read_indexed_options(unit_number)
    close(unit_number)
    if (nlo_combination_value == nlo_decay_multiplicative .and. &
        .not. has_nlo_contribution_bundle()) then
      call fail_parameters( &
           'multiplicative NLO decay combination requires a full NLO bundle')
    end if
    if (.not. momentum_mode_seen) then
      call fail_parameters('PRODUCTION_REN_SCALE_MOMENTA record is absent')
    end if
    if (production_w_system_value .and. use_decayed_production_momenta_value) &
         call fail_parameters('W_SYSTEM production scale grouping requires CORE momenta')
    if (production_current_sampling_value) then
      if (.not.sampling_mass_seen .or. .not.sampling_width_seen) &
           call fail_parameters('W_CURRENT requires explicit proposal mass and width')
      if (.not.ieee_is_finite(production_sampling_mass_value) .or. &
          .not.ieee_is_finite(production_sampling_width_value) .or. &
          production_sampling_mass_value <= 0d0 .or. production_sampling_width_value <= 0d0) &
           call fail_parameters('production proposal mass and width must be positive and finite')
    else if (sampling_mass_seen .or. sampling_width_seen) then
      call fail_parameters('proposal mass/width require W_CURRENT sampling')
    end if
    do width_index = 1, number_of_width_species
      if (.not. has_lo_width(width_index)) &
           call fail_parameters('an LO width is required for every species')
      if (lo_width_values(width_index) <= 0d0 .or. &
          .not. ieee_is_finite(lo_width_values(width_index))) then
        call fail_parameters('LO physical widths must be finite and positive')
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
          if (.not. has_nlo_width(width_index) .and. &
              species_order_is_nlo(node_pdg(node))) then
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
          if (.not. has_nlo_width(width_index) .and. &
              species_order_is_nlo(pdg)) then
            call fail_parameters(&
                 'the corrected species requires an NLO_DECAY_WIDTH record')
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
    call validate_scale_variations()
    initialized = .true.
  end subroutine initialize_decay_chain_parameters


  subroutine normalize_decay_card_record(line)
    character(len=*), intent(inout) :: line
    character(len=len(line)) :: key, value, indices
    integer :: pos, left, right, i, count, ios
    double precision :: factors(100)

    ! Use MG5 value = parameter syntax, with case-insensitive keywords and
    ! inline ! or # comments. Internally put the key and indices first.
    pos = scan(line, '!#')
    if (pos > 0) line(pos:) = ' '
    line = adjustl(line)
    if (len_trim(line) == 0) return
    pos = index(line, '=')
    if (pos <= 1 .or. pos == len_trim(line)) &
         call fail_parameters('expected value = parameter in decay_card.dat')
    value = adjustl(line(:pos-1))
    key = adjustl(line(pos+1:))
    if (index(key, '=') > 0) &
         call fail_parameters('expected one value = parameter assignment per line')
    call uppercase(key)
    indices = ' '
    left = index(key, '(')
    if (left > 0) then
      right = index(key, ')')
      if (right /= len_trim(key) .or. right <= left+1) &
           call fail_parameters('malformed indexed decay-card parameter')
      indices = key(left+1:right-1)
      key(left:) = ' '
    end if
    if (len_trim(key) == 0 .or. &
        verify(trim(key), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_') /= 0) &
         call fail_parameters('malformed decay-card parameter name')
    select case (trim(key))
    case ('LO_DECAY_WIDTH', 'NLO_DECAY_WIDTH', 'DECAY_REN_SCALE', &
          'DECAY_DYNAMICAL_SCALE_CHOICE', 'DECAY_WIDTH_SCALE_MODE', &
          'DECAY_PERTURBATIVE_ORDER', 'LO_DECAY_WIDTH_VARIATION', &
          'NLO_DECAY_WIDTH_VARIATION')
      if (len_trim(indices) == 0) &
           call fail_parameters('decay species parameters require a PDG index')
    case default
      if (len_trim(indices) > 0) &
           call fail_parameters('unexpected index on '//trim(key))
    end select
    do i = 1, len_trim(value)
      if (value(i:i) == "'" .or. value(i:i) == '"') value(i:i) = ' '
    end do
    line = trim(key)//' '//trim(indices)//' '//trim(adjustl(value))
    call uppercase(line)
    if (trim(key) == 'DECAY_SCALE_FACTORS') then
      ! Count list elements before allocating the scale-factor array.
      do i = 1, len_trim(value)
        if (scan(value(i:i), '[],') > 0) value(i:i) = ' '
      end do
      count = 0
      do i = 1, len_trim(value)
        if (value(i:i) == ' ' .or. value(i:i) == achar(9)) cycle
        if (i > 1) then
          if (value(i-1:i-1) /= ' ' .and. &
              value(i-1:i-1) /= achar(9)) cycle
        end if
        count = count + 1
      end do
      if (count < 1 .or. count > size(factors)) &
           call fail_parameters('invalid number of decay scale factors')
      read(value, *, iostat=ios) factors(1:count)
      if (ios /= 0) call fail_parameters('malformed decay scale factors')
      write(key, '(i0)') count
      line = 'DECAY_SCALE_FACTORS '//trim(key)//' '//trim(value)
    end if
  end subroutine normalize_decay_card_record


  subroutine uppercase(value)
    character(len=*), intent(inout) :: value
    integer :: i, code
    do i = 1, len_trim(value)
      code = iachar(value(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) &
           value(i:i) = achar(code + iachar('A') - iachar('a'))
    end do
  end subroutine uppercase


  subroutine read_indexed_options(unit_number)
    integer, intent(in) :: unit_number
    character(len=2048) :: line
    character(len=64) :: keyword, mode
    integer :: ios, pdg, idx, choice, node
    logical :: seen(3, size(scale_pdgs)), qcd_independent

    seen = .false.
    automatic_widths = .true.
    rewind(unit_number)
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_parameters('cannot read indexed decay options')
      call normalize_decay_card_record(line)
      if (skip_line(line)) cycle
      read(line, *) keyword
      select case (trim(keyword))
      case ('DECAY_DYNAMICAL_SCALE_CHOICE', 'DECAY_WIDTH_SCALE_MODE', &
            'DECAY_PERTURBATIVE_ORDER')
        read(line, *, iostat=ios) keyword, pdg, mode
        if (ios /= 0) call fail_parameters('malformed indexed decay option')
        idx = find_pdg(pdg, scale_pdgs)
        if (idx == 0) call fail_parameters('decay option has an unknown PDG')
        select case (trim(keyword))
        case ('DECAY_DYNAMICAL_SCALE_CHOICE')
          choice = 1
          read(mode, *, iostat=ios) dynamic_scale_choices(idx)
          if (ios /= 0) call fail_parameters('invalid decay scale choice')
          if (dynamic_scale_choices(idx) < -1 .or. &
              dynamic_scale_choices(idx) > 3) &
               call fail_parameters('decay scale choice must be -1, 0, 1, 2 or 3')
        case ('DECAY_WIDTH_SCALE_MODE')
          choice = 2
          if (trim(mode) /= 'AUTO' .and. trim(mode) /= 'EXPLICIT') &
               call fail_parameters('decay width scale mode must be AUTO or EXPLICIT')
          automatic_widths(idx) = trim(mode) == 'AUTO'
        case ('DECAY_PERTURBATIVE_ORDER')
          choice = 3
          if (trim(mode) /= 'LO' .and. trim(mode) /= 'NLO') &
               call fail_parameters('decay perturbative order must be LO or NLO')
          decay_orders(idx) = trim(mode) == 'NLO'
        end select
        if (seen(choice, idx)) call fail_parameters('duplicate indexed decay option')
        seen(choice, idx) = .true.
      end select
    end do
    do idx = 1, size(scale_pdgs)
      qcd_independent = .true.
      if (has_decay_chains()) then
        do node = 1, decay_node_count()
          if (abs(node_pdg(node)) == scale_pdgs(idx)) &
               qcd_independent = qcd_independent .and. node_qcd_order(node) == 0
        end do
      else
        do node = 1, nlo_decay_node_count()
          if (abs(nlo_decay_node_pdg(node)) == scale_pdgs(idx)) &
               qcd_independent = qcd_independent .and. nlo_decay_node_qcd_order(node) == 0
        end do
      end if
      if (automatic_widths(idx) .and. .not. qcd_independent) then
        if (seen(2, idx)) &
             call fail_parameters('AUTO widths require an alpha_s-independent Born decay')
        automatic_widths(idx) = .false.
      end if
      if (dynamic_scale_choices(idx) /= 0) then
        if (.not. qcd_independent) &
             call fail_parameters('event-by-event decay scales require an alpha_s-independent Born decay')
        if (.not. automatic_widths(idx)) &
             call fail_parameters('event-by-event decay scales require AUTO widths')
      end if
    end do
  end subroutine read_indexed_options


  subroutine set_decay_run_order(is_lo)
    logical, intent(in) :: is_lo
    if (lo_run .eqv. is_lo) return
    lo_run = is_lo
    if (initialized) then
      deallocate(scale_species_pdgs)
      call initialize_scale_species()
    end if
  end subroutine set_decay_run_order


  logical function species_order_is_nlo(pdg)
    integer, intent(in) :: pdg
    integer :: idx
    species_order_is_nlo = decays_nlo .and. .not. lo_run
    idx = find_pdg(pdg, scale_pdgs)
    if (idx > 0) species_order_is_nlo = species_order_is_nlo .and. decay_orders(idx)
  end function species_order_is_nlo


  logical function decay_species_nlo_enabled(pdg)
    integer, intent(in) :: pdg
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_species_nlo_enabled = species_order_is_nlo(pdg)
    if (has_nlo_contribution_bundle()) then
      decay_species_nlo_enabled = decay_species_nlo_enabled .and. bundle_species_is_nlo(pdg)
    else if (has_nlo_decay()) then
      decay_species_nlo_enabled = decay_species_nlo_enabled .and. &
           abs(pdg) == abs(corrected_parent_pdg())
    else
      decay_species_nlo_enabled = .false.
    end if
  end function decay_species_nlo_enabled


  logical function nlo_correction_enabled(contribution)
    integer, intent(in), optional :: contribution
    integer :: selected
    nlo_correction_enabled = .not. lo_run
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) return
    if (.not. initialized) call initialize_decay_chain_parameters()
    if (has_nlo_contribution_bundle()) then
      selected = active_nlo_contribution()
      if (present(contribution)) selected = contribution
      if (contribution_is_nlo_decay(selected)) then
        nlo_correction_enabled = decay_species_nlo_enabled(contribution_parent_pdg(selected))
        return
      end if
    else if (has_nlo_decay()) then
      nlo_correction_enabled = decay_species_nlo_enabled(corrected_parent_pdg())
      return
    end if
    nlo_correction_enabled = production_nlo .and. .not. lo_run
  end function nlo_correction_enabled


  integer function decay_dynamical_scale_choice(pdg)
    integer, intent(in) :: pdg
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_dynamical_scale_choice = dynamic_scale_choices(find_pdg(pdg, scale_pdgs))
  end function decay_dynamical_scale_choice


  logical function decay_automatic_width(pdg)
    integer, intent(in) :: pdg
    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_automatic_width = automatic_widths(find_pdg(pdg, scale_pdgs))
  end function decay_automatic_width


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
      decay_physical_width = decay_lo_width(pdg)
    else
      select_nlo_width = .false.
      if (present(use_nlo_width)) select_nlo_width = use_nlo_width
      if (select_nlo_width .and. has_nlo_width(width_index) .and. &
          decay_species_nlo_enabled(pdg)) then
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
    if (width_index == 0) then
      call fail_parameters('requested particle has no LO physical width')
    end if
    if (.not. has_lo_width(width_index)) then
      call fail_parameters('requested particle has no LO physical width')
    end if
    if (.not. present(factor_index)) then
      decay_lo_width = lo_width_values(width_index)
      return
    end if
    if (factor_index == 1) then
      decay_lo_width = lo_width_values(width_index)
      return
    end if
    call validate_factor_index(factor_index)
    variation_index = find_width_variation(&
         pdg, scale_factor_values(factor_index), &
         number_of_lo_width_variations, lo_variation_pdgs, &
         lo_variation_factors)
    if (variation_index == 0) then
      if (decay_automatic_width(pdg)) then
        decay_lo_width = lo_width_values(width_index)
        return
      end if
      call fail_parameters('requested LO width variation is absent')
    end if
    decay_lo_width = lo_variation_values(variation_index)
  end function decay_lo_width


  double precision function decay_nlo_width(pdg, factor_index, ren_scale)
    integer, intent(in) :: pdg
    integer, intent(in), optional :: factor_index
    double precision, intent(in), optional :: ren_scale
    integer :: width_index, variation_index, selected_index
    double precision :: selected_scale, reference_scale

    if (.not. initialized) call initialize_decay_chain_parameters()
    width_index = find_pdg(pdg, width_pdgs)
    if (width_index == 0) then
      call fail_parameters('requested particle has no NLO physical width')
    end if
    if (.not. has_nlo_width(width_index)) then
      call fail_parameters('requested particle has no NLO physical width')
    end if
    selected_index = 1
    if (present(factor_index)) selected_index = factor_index
    call validate_factor_index(selected_index)
    if (decay_automatic_width(pdg)) then
      reference_scale = decay_renormalization_scale(pdg)
      selected_scale = reference_scale
      if (present(ren_scale)) selected_scale = ren_scale
      selected_scale = selected_scale*decay_scale_factor(selected_index)
      ! For a QCD-independent Born width the complete NLO scale dependence
      ! is in its single power of alpha_s. The reference difference already
      ! includes that coupling: multiply by its ratio, never by alpha_s again.
      decay_nlo_width = lo_width_values(width_index) + &
           (nlo_width_values(width_index) - lo_width_values(width_index))* &
           alphas(selected_scale)/alphas(reference_scale)
      if (decay_nlo_width <= 0d0 .or. .not. ieee_is_finite(decay_nlo_width)) &
           call fail_parameters('the running NLO decay width is not positive and finite')
      return
    end if
    if (selected_index == 1) then
      decay_nlo_width = nlo_width_values(width_index)
      return
    end if
    variation_index = find_width_variation(&
         pdg, scale_factor_values(selected_index), &
         number_of_nlo_width_variations, nlo_variation_pdgs, &
         nlo_variation_factors)
    if (variation_index == 0) then
      call fail_parameters('requested NLO width variation is absent')
    end if
    decay_nlo_width = nlo_variation_values(variation_index)
  end function decay_nlo_width


  double precision function decay_width_expansion_coefficient(&
       factor_indices, node_scales)
    integer, intent(in), optional :: factor_indices(:)
    double precision, intent(in), optional :: node_scales(:)
    integer :: node, pdg, factor_index
    double precision :: ren_scale

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
      if (.not. decay_species_nlo_enabled(pdg)) cycle
      factor_index = selected_factor_index(pdg, factor_indices)
      ren_scale = decay_renormalization_scale(pdg)
      if (present(node_scales)) ren_scale = node_scales(node)
      decay_width_expansion_coefficient = &
           decay_width_expansion_coefficient - &
           (decay_nlo_width(pdg, factor_index, ren_scale) - &
            decay_lo_width(pdg, factor_index))/ &
           decay_lo_width(pdg, factor_index)
    end do
  end function decay_width_expansion_coefficient


  double precision function decay_width_denominator_rescaling(&
       factor_indices, node_scales)
    integer, intent(in), optional :: factor_indices(:)
    double precision, intent(in), optional :: node_scales(:)
    integer :: node, pdg, factor_index, node_count
    double precision :: ren_scale, selected_width
    logical :: corrected

    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_width_denominator_rescaling = 1d0
    if (has_nlo_decay()) then
      node_count = nlo_decay_node_count()
    else
      node_count = decay_node_count()
    end if
    do node = 1, node_count
      if (has_nlo_decay()) then
        pdg = nlo_decay_node_pdg(node)
        corrected = node == nlo_decay_corrected_node()
      else
        pdg = node_pdg(node)
        corrected = .false.
      end if
      factor_index = selected_factor_index(pdg, factor_indices)
      if (.not. has_nlo_contribution_bundle() .and. corrected .and. &
          decay_species_nlo_enabled(pdg)) then
        ren_scale = decay_renormalization_scale(pdg)
        if (present(node_scales)) ren_scale = node_scales(node)
        selected_width = decay_nlo_width(pdg, factor_index, ren_scale)
      else if (.not. has_lo_width(find_pdg(pdg, width_pdgs))) then
        selected_width = decay_nlo_width(pdg, factor_index)
      else
        selected_width = decay_lo_width(pdg, factor_index)
      end if
      decay_width_denominator_rescaling = &
           decay_width_denominator_rescaling*decay_physical_width(pdg, corrected)/selected_width
    end do
  end function decay_width_denominator_rescaling


  double precision function decay_multiplicative_width_rescaling( &
       factor_indices, node_scales)
    integer, intent(in), optional :: factor_indices(:)
    double precision, intent(in), optional :: node_scales(:)
    integer :: node, pdg, factor_index
    double precision :: selected_width, ren_scale

    if (.not. initialized) call initialize_decay_chain_parameters()
    decay_multiplicative_width_rescaling = 1d0
    if (.not. has_decay_chains()) return
    do node = 1, decay_node_count()
      pdg = node_pdg(node)
      factor_index = selected_factor_index(pdg, factor_indices)
      ren_scale = decay_renormalization_scale(pdg)
      if (present(node_scales)) ren_scale = node_scales(node)
      if (has_nlo_contribution_bundle() .and. &
          decay_species_nlo_enabled(pdg)) then
        selected_width = decay_nlo_width(pdg, factor_index, ren_scale)
      else
        selected_width = decay_lo_width(pdg, factor_index)
      end if
      ! The generated base measure is normalized with the central LO width.
      ! This single factor replaces it by the selected physical NLO width;
      ! no fixed-order width counterterm is part of a multiplicative leaf.
      decay_multiplicative_width_rescaling = &
           decay_multiplicative_width_rescaling* &
           decay_lo_width(pdg)/selected_width
    end do
  end function decay_multiplicative_width_rescaling


  double precision function decay_node_width_rescaling(pdg, factor_index, ren_scale)
    integer, intent(in) :: pdg, factor_index
    double precision, intent(in) :: ren_scale
    double precision :: selected_width
    selected_width = decay_lo_width(pdg, factor_index)
    if (decay_species_nlo_enabled(pdg)) &
         selected_width = decay_nlo_width(pdg, factor_index, ren_scale)
    decay_node_width_rescaling = decay_lo_width(pdg)/selected_width
  end function decay_node_width_rescaling


  integer function nlo_decay_combination_mode()
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) then
      nlo_decay_combination_mode = nlo_decay_additive
      return
    end if
    if (.not. initialized) call initialize_decay_chain_parameters()
    nlo_decay_combination_mode = nlo_combination_value
  end function nlo_decay_combination_mode


  logical function multiplicative_nlo_enabled()
    ! This query is used by shared fNLO code, including ordinary production
    ! processes with no decay metadata or decay_card.dat.  Such processes are
    ! necessarily on the ordinary additive path and must not initialize decay
    ! parameters merely to answer the mode question.
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) then
      multiplicative_nlo_enabled = .false.
      return
    end if
    if (.not. initialized) call initialize_decay_chain_parameters()
    multiplicative_nlo_enabled = &
         nlo_combination_value == nlo_decay_multiplicative .and. .not. lo_run
  end function multiplicative_nlo_enabled


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
    decay_scale_species_index = find_scale_pdg(pdg, scale_species_pdgs)
  end function decay_scale_species_index


  logical function use_decayed_production_ren_scale_momenta()
    if (.not. initialized) call initialize_decay_chain_parameters()
    use_decayed_production_ren_scale_momenta = &
         use_decayed_production_momenta_value
  end function use_decayed_production_ren_scale_momenta


  logical function production_w_system_scale()
    ! Undecayed generic processes need neither a decay card nor this option.
    production_w_system_scale = .false.
    if (.not. has_decay_chains() .and. .not. has_nlo_decay()) return
    if (.not. initialized) call initialize_decay_chain_parameters()
    production_w_system_scale = production_w_system_value
  end function production_w_system_scale


  subroutine production_current_proposal(enabled, mass, width)
    logical, intent(out) :: enabled
    double precision, intent(out) :: mass, width
    enabled = .false.
    mass = 0d0
    width = 0d0
    if (.not.has_decay_chains() .and. .not.has_nlo_decay()) return
    if (.not.initialized) call initialize_decay_chain_parameters()
    enabled = production_current_sampling_value
    mass = production_sampling_mass_value
    width = production_sampling_width_value
  end subroutine production_current_proposal


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
        pdg = node_pdg(node)
        scale_dependent = node_qcd_order(node) > 0
        if (has_nlo_contribution_bundle()) then
          scale_dependent = scale_dependent .or. &
               (bundle_species_is_nlo(pdg) .and. species_order_is_nlo(pdg))
        end if
      else
        pdg = nlo_decay_node_pdg(node)
        scale_dependent = nlo_decay_node_qcd_order(node) > 0 .or. &
             (abs(pdg) == abs(corrected_parent_pdg()) .and. species_order_is_nlo(pdg))
      end if
      if (.not. scale_dependent) cycle
      if (.not. signed_scale_axes) pdg = abs(pdg)
      if (find_scale_pdg(pdg, candidates) /= 0) cycle
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
        do while (previous >= 1)
          if (scale_species_pdgs(previous) <= pdg) exit
          scale_species_pdgs(previous + 1) = &
               scale_species_pdgs(previous)
          previous = previous - 1
        end do
        scale_species_pdgs(previous + 1) = pdg
      end do
    end if
    deallocate(candidates)
  end subroutine initialize_scale_species


  subroutine validate_scale_variations()
    integer :: index, factor_index, pdg

    if (scale_variation_mode_value == decay_scale_none) then
      if (number_of_lo_width_variations + number_of_nlo_width_variations > 0) &
           call fail_parameters('varied widths require decay scale reweighting')
    end if
    if (scale_variation_mode_value /= decay_scale_none .and. &
        number_of_scale_factors < 2) then
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
    if (scale_variation_mode_value == decay_scale_none) return

    call validate_variation_records(&
         number_of_lo_width_variations, lo_variation_pdgs, &
         lo_variation_factors, lo_variation_values, .false.)
    call validate_variation_records(&
         number_of_nlo_width_variations, nlo_variation_pdgs, &
         nlo_variation_factors, nlo_variation_values, .true.)
    do index = 1, number_of_scale_species
      pdg = scale_species_pdgs(index)
      if (automatic_widths(find_pdg(pdg, scale_pdgs))) cycle
      do factor_index = 2, number_of_scale_factors
        if (find_width_variation(&
             pdg, scale_factor_values(factor_index), &
             number_of_lo_width_variations, lo_variation_pdgs, &
             lo_variation_factors) == 0) then
          call fail_parameters(&
               'a varied decay species has no LO width at every factor')
        end if
        if (species_order_is_nlo(pdg) .and. &
            has_nlo_width(find_pdg(pdg, width_pdgs)) .and. &
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
      if (find_pdg(pdgs(index), scale_pdgs) == 0) then
        call fail_parameters(&
             'a width variation refers to a non-varied decay species')
      end if
      if (require_nlo_species .and. &
          .not. has_nlo_width(find_pdg(pdgs(index), width_pdgs))) then
        call fail_parameters(&
             'an NLO width variation refers to an uncorrected species')
      end if
      if (automatic_widths(find_pdg(pdgs(index), scale_pdgs))) &
           call fail_parameters('explicit width variations require EXPLICIT width scale mode')
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
    species_index = find_scale_pdg(pdg, scale_species_pdgs)
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


  integer function find_scale_pdg(pdg, pdgs)
    integer, intent(in) :: pdg, pdgs(:)
    integer :: index, key

    key = pdg
    if (.not. signed_scale_axes) key = abs(pdg)
    find_scale_pdg = 0
    do index = 1, size(pdgs)
      if (pdgs(index) == key) then
        find_scale_pdg = index
        return
      end if
    end do
  end function find_scale_pdg


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

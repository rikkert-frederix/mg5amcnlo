module multiplicative_fks_density
  use factorized_phase_space, only: factorized_radiation_state, &
       factorized_measure_state, fetch_factorized_radiation_state, &
       fetch_factorized_event_measure
  use fks_singular_module, only: evaluate_fks_sij, &
       record_soft_density_operator, sreal, sreal_deg, &
       fks_subtraction_shat
  use nlo_contribution_bundle, only: contribution_for_fks, &
       contribution_is_nlo_decay, contribution_corrected_node
  use fnlo_process_common, only: nfksprocess, i_fks, j_fks, &
       soft_counterevent, real_event, fkssymmetryfactor, &
       collinear_counterevent, soft_collinear_counterevent, &
       fkssymmetryfactordeg, xiscut_used, xicut_used, delta_used
  use spin_density_matrix_results, only: spin_density_real_insertion
  use multiplicative_density_terms, only: block_distribution_term, &
       density_scale_coefficient_count
  use density_operator_recorder, only: &
       begin_density_operator_recording, &
       recorded_density_operator_count, &
       recorded_density_operator_nonzero_count, &
       record_density_operator_primitive, &
       scale_recorded_density_operator, &
       cancel_density_operator_recording, &
       finish_density_operator_recording
  implicit none
  private

  public :: build_multiplicative_real_density_term
  public :: build_multiplicative_soft_density_term
  public :: build_multiplicative_collinear_density_term
  public :: build_multiplicative_soft_collinear_density_term

  interface
    integer function sdm_real_insertion_identifier(configuration)
      integer, intent(in) :: configuration
    end function sdm_real_insertion_identifier
  end interface

contains

  subroutine build_multiplicative_real_density_term( &
       contribution, term, available)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    type(factorized_radiation_state) :: radiation
    type(factorized_measure_state) :: measure
    complex(kind=8) :: coefficients(density_scale_coefficient_count)
    double precision :: damping, fks_partition, local_prefactor
    integer :: block, insertion_identifier
    logical :: radiation_available, measure_available

    available = .false.
    if (contribution_for_fks(nfksprocess) /= contribution) then
      call fail_multiplicative_fks_density( &
           'the requested contribution does not own the active FKS sector')
    end if
    block = contribution_physical_block(contribution)
    call fetch_factorized_radiation_state( &
         real_event, block, radiation, radiation_available)
    call fetch_factorized_event_measure( &
         real_event, block, measure, measure_available)
    if (.not. radiation_available .or. .not. measure_available) return
    if (radiation%jacobian <= 0d0 .or. measure%jacobian <= 0d0 .or. &
        measure%phase_space_weight <= 0d0) return
    if (radiation%xi <= 0d0 .or. 1d0 - radiation%y <= 0d0 .or. &
        radiation%xi_norm <= 0d0) return

    fks_partition = evaluate_fks_sij( &
         real_event, i_fks, j_fks, radiation%xi, radiation%y)
    if (fks_partition <= 0d0) return

    ! The phase-space Jacobian, phase-space weight, flux, incoming-channel
    ! weight, and eventually the single VEGAS weight are composed once for
    ! the complete Cartesian tuple.  This coefficient contains only the
    ! block-local FKS kernel multiplying the raw real density matrix.
    damping = radiation%xi**2*(1d0 - radiation%y)
    local_prefactor = radiation%xi_norm/radiation%xi/ &
         (1d0 - radiation%y)*fkssymmetryfactor
    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx( &
         damping*fks_partition*local_prefactor, 0d0, kind=8)
    insertion_identifier = sdm_real_insertion_identifier(nfksprocess)

    call begin_density_operator_recording(block, real_event, 1)
    call record_density_operator_primitive( &
         spin_density_real_insertion, insertion_identifier, 1, 0, &
         coefficients, .true.)
    call finish_density_operator_recording(term, 1)
    available = .true.
  end subroutine build_multiplicative_real_density_term


  subroutine build_multiplicative_soft_density_term( &
       contribution, term, available)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    type(factorized_radiation_state) :: real_radiation, soft_radiation
    type(factorized_measure_state) :: measure
    complex(kind=8) :: coefficients(density_scale_coefficient_count)
    double precision :: fks_partition, prefactor, endpoint_prefactor
    double precision :: cutoff
    integer :: block
    logical :: real_available, soft_available, measure_available

    available = .false.
    if (contribution_for_fks(nfksprocess) /= contribution) then
      call fail_multiplicative_fks_density( &
           'the requested contribution does not own the active FKS sector')
    end if
    block = contribution_physical_block(contribution)
    call fetch_factorized_radiation_state( &
         real_event, block, real_radiation, real_available)
    call fetch_factorized_radiation_state( &
         soft_counterevent, block, soft_radiation, soft_available)
    call fetch_factorized_event_measure( &
         soft_counterevent, block, measure, measure_available)
    if (.not. real_available .or. .not. soft_available .or. &
        .not. measure_available) return
    if (soft_radiation%jacobian <= 0d0 .or. measure%jacobian <= 0d0 .or. &
        measure%phase_space_weight <= 0d0) return
    if (real_radiation%xi <= 0d0 .or. &
        1d0 - real_radiation%y <= 0d0 .or. &
        real_radiation%xi_norm <= 0d0 .or. &
        real_radiation%xi_max <= 0d0 .or. &
        soft_radiation%xi_max <= 0d0) return
    if (real_radiation%xi_hat*soft_radiation%xi_max > xiscut_used .or. &
        real_radiation%xi > xiscut_used) return

    fks_partition = evaluate_fks_sij( &
         soft_counterevent, i_fks, j_fks, 0d0, real_radiation%y)
    if (fks_partition <= 0d0) return
    cutoff = min(real_radiation%xi_max, xiscut_used)
    if (cutoff <= 0d0) return
    prefactor = real_radiation%xi_norm/real_radiation%xi/ &
         (1d0 - real_radiation%y)
    endpoint_prefactor = real_radiation%xi_norm/cutoff* &
         log(xicut_used/cutoff)/(1d0 - real_radiation%y)
    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx( &
         fks_partition*(prefactor + endpoint_prefactor)* &
         fkssymmetryfactor, 0d0, kind=8)

    call begin_density_operator_recording( &
         block, soft_counterevent, 1)
    call record_soft_density_operator( &
         soft_counterevent, 0d0, real_radiation%y, &
         soft_radiation%sqrt_shat)
    if (recorded_density_operator_count() < 1) then
      call cancel_density_operator_recording()
      return
    end if
    call scale_recorded_density_operator(coefficients)
    call finish_density_operator_recording(term, -1)
    available = .true.
  end subroutine build_multiplicative_soft_density_term


  subroutine build_multiplicative_collinear_density_term( &
       contribution, term, available)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    type(factorized_radiation_state) :: real_radiation
    type(factorized_radiation_state) :: collinear_radiation
    type(factorized_measure_state) :: measure
    double precision :: fks_partition, collinear_prefactor
    double precision :: collinear_endpoint_prefactor, remnant_prefactor
    double precision :: dummy_weight, dummy_xi, dummy_lxi, subtraction_shat
    integer :: block, core_first, core_last
    integer :: xi_index, lxi_index, muf_index
    logical :: real_available, collinear_available, measure_available

    available = .false.
    if (contribution_for_fks(nfksprocess) /= contribution) then
      call fail_multiplicative_fks_density( &
           'the requested contribution does not own the active FKS sector')
    end if
    block = contribution_physical_block(contribution)
    call fetch_factorized_radiation_state( &
         real_event, block, real_radiation, real_available)
    call fetch_factorized_radiation_state( &
         collinear_counterevent, block, collinear_radiation, &
         collinear_available)
    call fetch_factorized_event_measure( &
         collinear_counterevent, block, measure, measure_available)
    if (.not. real_available .or. .not. collinear_available .or. &
        .not. measure_available) return
    if (collinear_radiation%jacobian <= 0d0 .or. &
        measure%jacobian <= 0d0 .or. measure%phase_space_weight <= 0d0) return
    if (real_radiation%y <= 0d0 .or. &
        1d0 - real_radiation%y <= 0d0 .or. &
        collinear_radiation%xi <= 0d0 .or. &
        collinear_radiation%xi_norm <= 0d0) return

    fks_partition = evaluate_fks_sij( &
         collinear_counterevent, i_fks, j_fks, &
         collinear_radiation%xi, 1d0)
    if (fks_partition <= 0d0) return
    collinear_prefactor = collinear_radiation%xi_norm/ &
         collinear_radiation%xi/(1d0 - real_radiation%y)
    collinear_endpoint_prefactor = collinear_radiation%xi_norm/ &
         collinear_radiation%xi*log(delta_used)/1d0
    subtraction_shat = fks_subtraction_shat(collinear_counterevent)
    if (subtraction_shat <= 0d0) return
    remnant_prefactor = collinear_radiation%xi_norm/ &
         collinear_radiation%xi/1d0/ &
         (subtraction_shat/(32d0*acos(-1d0)**2))* &
         fkssymmetryfactordeg

    call begin_density_operator_recording( &
         block, collinear_counterevent, 1)
    call sreal_deg( &
         collinear_counterevent, collinear_radiation%xi, dummy_xi, &
         dummy_lxi, xi_index, lxi_index, muf_index)
    call scale_one_recorded_primitive(xi_index, -remnant_prefactor)
    call scale_one_recorded_primitive( &
         lxi_index, -log(collinear_radiation%xi)*remnant_prefactor)
    call scale_one_recorded_primitive(muf_index, -remnant_prefactor)

    core_first = recorded_density_operator_count() + 1
    call sreal( &
         collinear_counterevent, collinear_radiation%xi, 1d0, dummy_weight)
    core_last = recorded_density_operator_count()
    if (core_last >= core_first) then
      call scale_recorded_range( &
           core_first, core_last, fks_partition* &
           (collinear_prefactor + collinear_endpoint_prefactor)* &
           fkssymmetryfactor)
    end if
    if (recorded_density_operator_nonzero_count() < 1) then
      call cancel_density_operator_recording()
      return
    end if
    call finish_density_operator_recording(term, -1)
    available = .true.
  end subroutine build_multiplicative_collinear_density_term


  subroutine build_multiplicative_soft_collinear_density_term( &
       contribution, term, available)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    type(factorized_radiation_state) :: real_radiation
    type(factorized_radiation_state) :: collinear_radiation
    type(factorized_radiation_state) :: soft_collinear_radiation
    type(factorized_measure_state) :: measure
    double precision :: fks_partition, prefactor_c, prefactor_coll
    double precision :: prefactor_cnt, prefactor_cnt_coll, prefactor_deg
    double precision :: prefactor_deg_sxi, prefactor_deg_slxi
    double precision :: local_fsc, local_fdsc1, local_fdsc2
    double precision :: local_fdsc3, local_fdsc4, cutoff
    double precision :: dummy_weight, dummy_xi, dummy_lxi, subtraction_shat
    integer :: block, core_first, core_last
    integer :: xi_index, lxi_index, muf_index
    logical :: real_available, collinear_available, sc_available
    logical :: measure_available

    available = .false.
    if (contribution_for_fks(nfksprocess) /= contribution) then
      call fail_multiplicative_fks_density( &
           'the requested contribution does not own the active FKS sector')
    end if
    block = contribution_physical_block(contribution)
    call fetch_factorized_radiation_state( &
         real_event, block, real_radiation, real_available)
    call fetch_factorized_radiation_state( &
         collinear_counterevent, block, collinear_radiation, &
         collinear_available)
    call fetch_factorized_radiation_state( &
         soft_collinear_counterevent, block, soft_collinear_radiation, &
         sc_available)
    call fetch_factorized_event_measure( &
         soft_collinear_counterevent, block, measure, measure_available)
    if (.not. real_available .or. .not. collinear_available .or. &
        .not. sc_available .or. .not. measure_available) return
    if (soft_collinear_radiation%jacobian <= 0d0 .or. &
        measure%jacobian <= 0d0 .or. measure%phase_space_weight <= 0d0) return
    if (real_radiation%y <= 0d0 .or. &
        1d0 - real_radiation%y <= 0d0 .or. &
        collinear_radiation%xi <= 0d0 .or. &
        collinear_radiation%xi_norm <= 0d0 .or. &
        collinear_radiation%xi_max <= 0d0) return
    if (real_radiation%xi_hat*collinear_radiation%xi_max >= &
        xiscut_used .or. collinear_radiation%xi >= xiscut_used) return

    fks_partition = evaluate_fks_sij( &
         soft_collinear_counterevent, i_fks, j_fks, 0d0, 1d0)
    if (fks_partition <= 0d0) return
    cutoff = min(collinear_radiation%xi_max, xiscut_used)
    if (cutoff <= 0d0) return
    prefactor_c = collinear_radiation%xi_norm/ &
         collinear_radiation%xi/(1d0 - real_radiation%y)
    prefactor_coll = collinear_radiation%xi_norm/ &
         collinear_radiation%xi*log(delta_used)
    prefactor_cnt = collinear_radiation%xi_norm/cutoff* &
         log(xicut_used/cutoff)/(1d0 - real_radiation%y)
    prefactor_cnt_coll = collinear_radiation%xi_norm/cutoff* &
         log(xicut_used/cutoff)*log(delta_used)
    prefactor_deg = collinear_radiation%xi_norm/ &
         collinear_radiation%xi
    prefactor_deg_sxi = collinear_radiation%xi_norm/cutoff* &
         log(xicut_used/cutoff)
    prefactor_deg_slxi = collinear_radiation%xi_norm/cutoff* &
         (log(xicut_used)**2 - log(cutoff)**2)/2d0
    subtraction_shat = fks_subtraction_shat( &
         soft_collinear_counterevent)
    if (subtraction_shat <= 0d0) return
    local_fsc = (prefactor_c + prefactor_coll + prefactor_cnt + &
         prefactor_cnt_coll)*fkssymmetryfactordeg
    local_fdsc1 = prefactor_deg/ &
         (subtraction_shat/(32d0*acos(-1d0)**2))* &
         fkssymmetryfactordeg
    local_fdsc2 = prefactor_deg_sxi*fkssymmetryfactordeg
    local_fdsc3 = prefactor_deg_slxi*fkssymmetryfactordeg
    local_fdsc4 = prefactor_deg_sxi*fkssymmetryfactordeg

    call begin_density_operator_recording( &
         block, soft_collinear_counterevent, 1)
    call sreal_deg( &
         soft_collinear_counterevent, 0d0, dummy_xi, dummy_lxi, &
         xi_index, lxi_index, muf_index)
    call scale_one_recorded_primitive( &
         xi_index, -(local_fdsc1 + local_fdsc2))
    call scale_one_recorded_primitive( &
         lxi_index, -(log(collinear_radiation%xi)*local_fdsc1 + &
                      local_fdsc3))
    call scale_one_recorded_primitive(muf_index, -local_fdsc4)

    core_first = recorded_density_operator_count() + 1
    call sreal( &
         soft_collinear_counterevent, 0d0, 1d0, dummy_weight)
    core_last = recorded_density_operator_count()
    if (core_last >= core_first) then
      call scale_recorded_range( &
           core_first, core_last, fks_partition*local_fsc)
    end if
    if (recorded_density_operator_nonzero_count() < 1) then
      call cancel_density_operator_recording()
      return
    end if
    call finish_density_operator_recording(term, 1)
    available = .true.
  end subroutine build_multiplicative_soft_collinear_density_term


  subroutine scale_one_recorded_primitive(primitive, factor)
    integer, intent(in) :: primitive
    double precision, intent(in) :: factor

    if (primitive == 0) return
    call scale_recorded_range(primitive, primitive, factor)
  end subroutine scale_one_recorded_primitive


  subroutine scale_recorded_range(first, last, factor)
    integer, intent(in) :: first, last
    double precision, intent(in) :: factor
    complex(kind=8) :: coefficients(density_scale_coefficient_count)

    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx(factor, 0d0, kind=8)
    call scale_recorded_density_operator(coefficients, first, last)
  end subroutine scale_recorded_range


  integer function contribution_physical_block(contribution)
    integer, intent(in) :: contribution

    if (contribution_is_nlo_decay(contribution)) then
      contribution_physical_block = &
           contribution_corrected_node(contribution)
    else
      contribution_physical_block = 0
    end if
  end function contribution_physical_block


  subroutine fail_multiplicative_fks_density(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_fks_density: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_fks_density
end module multiplicative_fks_density

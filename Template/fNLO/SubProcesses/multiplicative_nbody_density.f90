module multiplicative_nbody_density
  use factorized_phase_space, only: factorized_radiation_state, &
       factorized_measure_state, fetch_factorized_radiation_state, &
       fetch_factorized_event_measure
  use fks_singular_module, only: fks_subtraction_shat, &
       record_nbody_integrated_density_operator
  use nlo_contribution_bundle, only: contribution_for_fks, &
       contribution_has_virtual, contribution_has_fast_virtual, &
       contribution_is_nlo_decay, &
       contribution_corrected_node, contribution_representative_fks
  use fnlo_process_common, only: nfksprocess, soft_counterevent, real_event, &
       xibsvcut_used
  use spin_density_matrix_results, only: spin_density_no_insertion, &
       spin_density_virtual_insertion, &
       spin_density_fast_virtual_insertion
  use multiplicative_density_terms, only: block_distribution_term, &
       density_scale_coefficient_count
  use density_operator_recorder, only: &
       begin_density_operator_recording, &
       recorded_density_operator_nonzero_count, &
       record_density_operator_primitive, &
       scale_recorded_density_operator, &
       cancel_density_operator_recording, &
       finish_density_operator_recording
  implicit none
  private

  public :: build_multiplicative_lo_density_term
  public :: build_multiplicative_lo_only_density_term
  public :: build_multiplicative_nbody_density_term

  interface
    integer function sdm_multiplicative_born_qcd_power(block)
      integer, intent(in) :: block
    end function sdm_multiplicative_born_qcd_power
  end interface

contains

  subroutine build_multiplicative_lo_only_density_term(block, term)
    integer, intent(in) :: block
    type(block_distribution_term), intent(out) :: term
    complex(kind=8) :: coefficients(density_scale_coefficient_count)

    ! A QCD-inert or otherwise LO-only decay block has no auxiliary FKS
    ! variables to project out.  Its phase-space measure is part of the
    ! factorized base measure, so its operator coefficient is exactly one.
    coefficients = (0d0, 0d0)
    coefficients(1) = (1d0, 0d0)
    call begin_density_operator_recording(block, soft_counterevent, 0)
    call record_density_operator_primitive( &
         spin_density_no_insertion, 0, 0, 0, coefficients, .true.)
    call finish_density_operator_recording(term, 1, 0)
  end subroutine build_multiplicative_lo_only_density_term

  subroutine build_multiplicative_lo_density_term( &
       contribution, term, available, channel_partition)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    double precision, intent(in), optional :: channel_partition
    complex(kind=8) :: coefficients(density_scale_coefficient_count)
    double precision :: local_prefactor, partition
    integer :: block

    available = .false.
    call resolve_channel_partition(channel_partition, partition)
    call validate_active_contribution(contribution, block)
    call nbody_local_prefactor(block, local_prefactor, available)
    if (.not. available) return
    local_prefactor = local_prefactor*partition

    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx(local_prefactor, 0d0, kind=8)
    call begin_density_operator_recording( &
         block, soft_counterevent, 0)
    call record_density_operator_primitive( &
         spin_density_no_insertion, 0, 0, 0, coefficients, .true.)
    call finish_density_operator_recording( &
         term, 1, contribution_representative_fks(contribution))
    available = .true.
  end subroutine build_multiplicative_lo_density_term


  subroutine build_multiplicative_nbody_density_term( &
       contribution, term, available, channel_partition)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    double precision, intent(in), optional :: channel_partition
    complex(kind=8) :: coefficients(density_scale_coefficient_count)
    double precision :: local_prefactor, partition
    integer :: block, born_qcd_power

    available = .false.
    call resolve_channel_partition(channel_partition, partition)
    call validate_active_contribution(contribution, block)
    call nbody_local_prefactor(block, local_prefactor, available)
    if (.not. available) return
    local_prefactor = local_prefactor*partition
    born_qcd_power = sdm_multiplicative_born_qcd_power(block)
    if (born_qcd_power < 0) then
      call fail_multiplicative_nbody_density( &
           'a generated block has a negative Born QCD power')
    end if

    call begin_density_operator_recording( &
         block, soft_counterevent, 1)
    call record_nbody_integrated_density_operator( &
         soft_counterevent, born_qcd_power)
    if (contribution_has_virtual(contribution)) then
      coefficients = (0d0, 0d0)
      coefficients(1) = (1d0, 0d0)
      ! Rank one is the finite Laurent coefficient.  Ranks two and three of
      ! the generated loop provider are pole diagnostics and are never
      ! selected by a multiplicative primitive.
      call record_density_operator_primitive( &
           merge(spin_density_fast_virtual_insertion, &
                 spin_density_virtual_insertion, &
                 contribution_has_fast_virtual(contribution)), &
           contribution, 1, 0, &
           coefficients, .true.)
    end if
    if (recorded_density_operator_nonzero_count() < 1) then
      call cancel_density_operator_recording()
      available = .false.
      return
    end if
    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx(local_prefactor, 0d0, kind=8)
    call scale_recorded_density_operator(coefficients)
    call finish_density_operator_recording( &
         term, 1, contribution_representative_fks(contribution))
    available = .true.
  end subroutine build_multiplicative_nbody_density_term


  subroutine nbody_local_prefactor(block, prefactor, available)
    integer, intent(in) :: block
    double precision, intent(out) :: prefactor
    logical, intent(out) :: available
    double precision, parameter :: pi = 3.1415926535897932385d0
    type(factorized_radiation_state) :: real_radiation, soft_radiation
    type(factorized_measure_state) :: measure
    double precision :: cutoff, subtraction_shat
    logical :: real_available, soft_available, measure_available

    prefactor = 0d0
    available = .false.
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
    if (real_radiation%xi_norm <= 0d0 .or. &
        real_radiation%xi_max <= 0d0 .or. &
        soft_radiation%xi_max <= 0d0) return
    if (real_radiation%xi_hat*soft_radiation%xi_max > &
        xibsvcut_used) return
    cutoff = min(real_radiation%xi_max, xibsvcut_used)
    if (cutoff <= 0d0) return
    subtraction_shat = fks_subtraction_shat(soft_counterevent)
    if (subtraction_shat <= 0d0) return

    ! The selected tuple measure supplies the Born phase spaces, this
    ! block's dummy-radiation Jacobian/weight, the one physical flux, and the
    ! incoming integration weight.  This is only the local FKS projection
    ! factor that makes the three n-body radiation variables integrate to
    ! unity, matching COMPUTE_PREFACTORS_NBODY with its measure removed.
    ! The caller separately cancels the inverse probability of the sampled
    ! FKS sector for these sector-independent n-body atoms.  In particular,
    ! do not apply FKSSYMMETRYFACTORBORN here: that factor selects soft-gluon
    ! sectors in the additive sector sum, whereas every multiplicative draw
    ! must carry the same single Born/virtual atom.
    prefactor = real_radiation%xi_norm/ &
         (cutoff*subtraction_shat/(16d0*pi**2))
    available = prefactor > 0d0
  end subroutine nbody_local_prefactor


  subroutine resolve_channel_partition(requested, partition)
    double precision, intent(in), optional :: requested
    double precision, intent(out) :: partition

    partition = 1d0
    if (present(requested)) partition = requested
    if (partition <= 0d0 .or. partition > 1d0 .or. &
        partition /= partition) then
      call fail_multiplicative_nbody_density( &
           'the n-body integration-channel partition is invalid')
    end if
  end subroutine resolve_channel_partition


  subroutine validate_active_contribution(contribution, block)
    integer, intent(in) :: contribution
    integer, intent(out) :: block

    if (contribution_for_fks(nfksprocess) /= contribution) then
      call fail_multiplicative_nbody_density( &
           'the requested contribution does not own the active FKS sector')
    end if
    if (contribution_is_nlo_decay(contribution)) then
      block = contribution_corrected_node(contribution)
    else
      block = 0
    end if
  end subroutine validate_active_contribution


  subroutine fail_multiplicative_nbody_density(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_nbody_density: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_nbody_density
end module multiplicative_nbody_density

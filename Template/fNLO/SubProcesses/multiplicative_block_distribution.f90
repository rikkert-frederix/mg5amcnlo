module multiplicative_block_distribution
  use multiplicative_density_terms, only: block_distribution_term, &
       block_nlo_distribution, initialize_block_distribution, &
       finalize_block_distribution
  use multiplicative_nbody_density, only: &
       build_multiplicative_lo_only_density_term, &
       build_multiplicative_lo_density_term, &
       build_multiplicative_nbody_density_term
  use multiplicative_fks_density, only: &
       build_multiplicative_real_density_term, &
       build_multiplicative_soft_density_term, &
       build_multiplicative_collinear_density_term, &
       build_multiplicative_soft_collinear_density_term
  use multiplicative_scale_state, only: &
       activate_multiplicative_block_reference
  implicit none
  private

  integer, parameter :: maximum_block_terms = 6

  public :: build_multiplicative_block_nlo_distribution
  public :: build_multiplicative_lo_only_distribution

contains

  subroutine build_multiplicative_lo_only_distribution(block, distribution)
    integer, intent(in) :: block
    type(block_nlo_distribution), intent(out) :: distribution
    type(block_distribution_term) :: term

    call activate_multiplicative_block_reference(block)
    call build_multiplicative_lo_only_density_term(block, term)
    call initialize_block_distribution(distribution, block, 1)
    distribution%terms(1) = term
    call finalize_block_distribution(distribution)
  end subroutine build_multiplicative_lo_only_distribution

  subroutine build_multiplicative_block_nlo_distribution( &
       contribution, distribution, available, nbody_partition)
    integer, intent(in) :: contribution
    type(block_nlo_distribution), intent(out) :: distribution
    logical, intent(out) :: available
    double precision, intent(in), optional :: nbody_partition
    type(block_distribution_term) :: candidates(maximum_block_terms)
    type(block_distribution_term) :: candidate
    logical :: candidate_available
    integer :: candidate_count, term, block

    available = .false.
    candidate_count = 0

    ! Integrated kernels and any later loop calls must use exactly the same
    ! immutable reference scale and coupling for this block.
    call activate_multiplicative_block_reference_for_contribution( &
         contribution)

    call build_multiplicative_lo_density_term( &
         contribution, candidate, candidate_available, nbody_partition)
    if (candidate_available) call append_candidate(candidate)

    call build_multiplicative_nbody_density_term( &
         contribution, candidate, candidate_available, nbody_partition)
    if (candidate_available) call append_candidate(candidate)

    call build_multiplicative_real_density_term( &
         contribution, candidate, candidate_available)
    if (candidate_available) call append_candidate(candidate)

    call build_multiplicative_soft_density_term( &
         contribution, candidate, candidate_available)
    if (candidate_available) call append_candidate(candidate)

    call build_multiplicative_collinear_density_term( &
         contribution, candidate, candidate_available)
    if (candidate_available) call append_candidate(candidate)

    call build_multiplicative_soft_collinear_density_term( &
         contribution, candidate, candidate_available)
    if (candidate_available) call append_candidate(candidate)

    ! Some FKS directories (notably a q-qbar pair replacing a Born gluon)
    ! have FKSSYMMETRYFACTORBORN=0.  They are valid NLO-only Monte Carlo
    ! samples: omitting an LO atom in this one sampled sector is what makes
    ! the sector sum reproduce one Born contribution without duplication.
    if (candidate_count == 0) return

    block = candidates(1)%block
    call initialize_block_distribution( &
         distribution, block, candidate_count)
    do term = 1, candidate_count
      if (candidates(term)%block /= block) then
        call fail_block_distribution( &
             'one block distribution contains several physical blocks')
      end if
      if (.not. candidates(term)%finalized) then
        call fail_block_distribution( &
             'one block distribution contains an unfinished term')
      end if
      distribution%terms(term) = candidates(term)
    end do
    call finalize_block_distribution(distribution)
    available = .true.

  contains

    subroutine append_candidate(value)
      type(block_distribution_term), intent(in) :: value

      if (.not. value%finalized) then
        call fail_block_distribution( &
             'a block builder returned an unfinished term')
      end if
      if (candidate_count >= maximum_block_terms) then
        call fail_block_distribution( &
             'a block distribution exceeds its term capacity')
      end if
      candidate_count = candidate_count + 1
      candidates(candidate_count) = value
    end subroutine append_candidate
  end subroutine build_multiplicative_block_nlo_distribution


  subroutine activate_multiplicative_block_reference_for_contribution( &
       contribution)
    use nlo_contribution_bundle, only: contribution_is_nlo_decay, &
         contribution_corrected_node
    integer, intent(in) :: contribution
    integer :: block

    block = 0
    if (contribution_is_nlo_decay(contribution)) then
      block = contribution_corrected_node(contribution)
    end if
    call activate_multiplicative_block_reference(block)
  end subroutine activate_multiplicative_block_reference_for_contribution


  subroutine fail_block_distribution(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_block_distribution: '//trim(message)
    stop 1
  end subroutine fail_block_distribution
end module multiplicative_block_distribution

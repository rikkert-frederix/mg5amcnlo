module factorized_mint_policy
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       factorized_radiation_block, bundle_nlo_component
  use decay_chain_parameters, only: uses_multiplicative_nlo_combination
  implicit none
  private

  ! Keep process-specific training and channel-selection choices outside the
  ! generic MINT state machine.  MINT supplies only a dimension and its vector
  ! of diagnostic integrals; this adapter decides which positive proxy trains
  ! that dimension.
  public :: factorized_mint_grid_weight
  public :: factorized_mint_uses_uniform_channels
  public :: factorized_mint_shows_multiplicative_validation

contains

  subroutine factorized_mint_grid_weight( &
       total_dimension, dimension, integrals, lo_integral, &
       first_component_integral, grid_weight)
    integer, intent(in) :: total_dimension, dimension
    double precision, intent(in) :: integrals(:)
    integer, intent(in) :: lo_integral, first_component_integral
    double precision, intent(out) :: grid_weight
    integer :: radiation_block, component, component_integral

    if (dimension < 1 .or. dimension > total_dimension) then
      call fail_factorized_mint_policy('a MINT dimension is out of range')
    end if
    if (size(integrals) < 1) then
      call fail_factorized_mint_policy('the MINT integral vector is empty')
    end if
    grid_weight = integrals(1)
    if (uses_multiplicative_nlo_combination()) then
      if (lo_integral < 1 .or. lo_integral > size(integrals)) then
        call fail_factorized_mint_policy( &
             'the all-LO MINT proxy is out of range')
      end if
      grid_weight = abs(integrals(lo_integral))
    end if
    if (.not. has_nlo_contribution_bundle()) return

    radiation_block = factorized_radiation_block( &
         total_dimension, dimension)
    if (radiation_block <= 0) return
    component = bundle_nlo_component(radiation_block)
    component_integral = first_component_integral + component - 1
    if (component_integral < 1 .or. &
        component_integral > size(integrals)) then
      call fail_factorized_mint_policy( &
           'a block-local MINT proxy is out of range')
    end if
    grid_weight = integrals(component_integral)
  end subroutine factorized_mint_grid_weight


  logical function factorized_mint_uses_uniform_channels()
    factorized_mint_uses_uniform_channels = &
         uses_multiplicative_nlo_combination()
  end function factorized_mint_uses_uniform_channels


  logical function factorized_mint_shows_multiplicative_validation()
    factorized_mint_shows_multiplicative_validation = &
         has_nlo_contribution_bundle() .and. &
         uses_multiplicative_nlo_combination()
  end function factorized_mint_shows_multiplicative_validation


  subroutine fail_factorized_mint_policy(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in factorized_mint_policy: '//trim(message)
    stop 1
  end subroutine fail_factorized_mint_policy

end module factorized_mint_policy

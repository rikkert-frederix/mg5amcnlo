module fixed_order_user_hooks
  implicit none
  private

  public :: accept_dummy_cuts
  public :: fixed_user_scale
  public :: fixed_user_decay_scale

contains

  logical function accept_dummy_cuts()
    implicit none

    accept_dummy_cuts = .true.
  end function accept_dummy_cuts

  double precision function fixed_user_scale(mu_r_reference, scale_id)
    implicit none
    double precision, intent(in) :: mu_r_reference
    character(len=*), intent(out) :: scale_id

    fixed_user_scale = mu_r_reference
    scale_id = 'fixed scale'
  end function fixed_user_scale

  double precision function fixed_user_decay_scale(pdg, node, p, reference)
    integer, intent(in) :: pdg, node
    double precision, intent(in) :: p(0:, :), reference
    ! Used for decay_dynamical_scale_choice(PDG) = -1. P contains this
    ! decay block in the event frame: parent first, then immediate daughters
    ! (including the emitted parton for a resolved real event). NODE labels
    ! the physical resonance occurrence. Return a positive scale in GeV.
    ! Keep the prescription infrared safe so real and counterevents agree
    ! in unresolved limits. The default hook returns the card reference.
    fixed_user_decay_scale = reference
  end function fixed_user_decay_scale

end module fixed_order_user_hooks

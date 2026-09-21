! Dedicated direct-scale validation hook, NOT the main campaign default.
! Install only in a fresh export before compilation and archive its checksum.
module fixed_order_user_hooks
  implicit none
  private
  public :: accept_dummy_cuts, fixed_user_scale, fixed_user_decay_scale
contains
  logical function accept_dummy_cuts()
    accept_dummy_cuts = .true.
  end function accept_dummy_cuts

  double precision function fixed_user_scale(mu_r_reference, scale_id)
    double precision, intent(in) :: mu_r_reference
    character(len=*), intent(out) :: scale_id
    fixed_user_scale = mu_r_reference
    scale_id = 'fixed scale'
  end function fixed_user_scale

  double precision function fixed_user_decay_scale(pdg, node, p, reference)
    integer, intent(in) :: pdg, node
    double precision, intent(in) :: p(0:, :), reference
    ! Constant within each signed species: real/counterevent limits coincide.
    ! AUTO widths retain their coefficient reference at mt and evaluate their
    ! coupling at this local scale. No duplicate negative-PDG width is added.
    select case (pdg)
    case (6)
      fixed_user_decay_scale = 2d0*reference
    case (-6)
      fixed_user_decay_scale = .5d0*reference
    case default
      fixed_user_decay_scale = reference
    end select
  end function fixed_user_decay_scale
end module fixed_order_user_hooks

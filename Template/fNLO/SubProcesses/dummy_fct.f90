module fixed_order_user_hooks
  implicit none
  private

  public :: accept_dummy_cuts
  public :: fixed_user_scale

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

end module fixed_order_user_hooks

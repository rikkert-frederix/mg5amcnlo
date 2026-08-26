module binoth_lha_user_backend
  use process_dimensions, only: amp_split_size
  implicit none
  private

  public :: binoth_lha_eval

contains

  subroutine binoth_lha_eval(virt_wgt, amp_split, amp_split_finite)
    implicit none
    double precision, intent(out) :: virt_wgt
    double precision, intent(inout) :: amp_split(amp_split_size)
    double precision, intent(inout) :: amp_split_finite(amp_split_size)

    virt_wgt = 0d0
    amp_split = 0d0
    amp_split_finite = 0d0

    ! Replace this stub with the desired one-loop provider.  The provider
    ! must fill the finite virtual contribution and its split-order weights.
  end subroutine binoth_lha_eval

end module binoth_lha_user_backend

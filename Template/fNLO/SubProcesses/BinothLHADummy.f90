module binoth_lha_dummy
  implicit none
  private

  public :: BinothLHA

contains

  subroutine BinothLHA(virt_wgt)
    implicit none
    double precision, intent(out) :: virt_wgt

    virt_wgt = 0d0
  end subroutine BinothLHA

end module binoth_lha_dummy

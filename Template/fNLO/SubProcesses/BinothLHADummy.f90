module binoth_lha_dummy
  use process_dimensions, only: nexternal
  implicit none
  private

  public :: BinothLHA

contains

  subroutine BinothLHA(p_born, born_wgt, virt_wgt)
    implicit none
    double precision, intent(in) :: p_born(0:3, nexternal-1)
    double precision, intent(in) :: born_wgt
    double precision, intent(out) :: virt_wgt

    virt_wgt = 0d0
  end subroutine BinothLHA

end module binoth_lha_dummy

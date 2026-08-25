module binoth_lha_dummy
  use process_dimensions, only: nexternal
  implicit none
  private

  public :: BinothLHA, BinothLHAInit, ctsstatistics, sloopmatrix

contains

  subroutine BinothLHA(p_born, born_wgt, virt_wgt)
    implicit none
    double precision, intent(in) :: p_born(0:3, nexternal-1)
    double precision, intent(in) :: born_wgt
    double precision, intent(out) :: virt_wgt

    virt_wgt = 0d0
  end subroutine BinothLHA

  subroutine BinothLHAInit(filename)
    implicit none
    character(len=*), intent(in) :: filename
  end subroutine BinothLHAInit

  subroutine ctsstatistics(n_mp, n_disc)
    implicit none
    integer, intent(out) :: n_mp, n_disc

    n_mp = 0
    n_disc = 0
  end subroutine ctsstatistics

  subroutine sloopmatrix(p_born, virt_wgts)
    implicit none
    double precision, intent(in) :: p_born(0:3, nexternal-1)
    double precision, intent(out) :: virt_wgts(3)

    virt_wgts(1) = 0d0
    virt_wgts(2) = 0d0
    virt_wgts(3) = 0d0
  end subroutine sloopmatrix

end module binoth_lha_dummy

module pdfwrap_module
  implicit none
  private

  integer, parameter :: dp = kind(1.0d0)

  public :: configure_pdf
  public :: initpdfm_impl
  public :: initpdfsetbynamem_impl
  public :: numberpdfm_impl

  interface
    subroutine nnpdfdriver(gridfilename)
      character(len=*), intent(in) :: gridfilename
    end subroutine nnpdfdriver

    subroutine nninitpdf(irep)
      integer, intent(in) :: irep
    end subroutine nninitpdf
  end interface

contains

  subroutine configure_pdf(pdlabel, asmz, nloop)
    implicit none
    character(len=*), intent(in) :: pdlabel
    real(dp), intent(out) :: asmz
    integer, intent(out) :: nloop

    nloop = 2

    select case (trim(pdlabel))
    case ('nn23lo')
      call nnpdfdriver('NNPDF23_lo_as_0119_qed_mem0.grid')
      call nninitpdf(0)
      asmz = 0.119d0
    case ('nn23lo1')
      call nnpdfdriver('NNPDF23_lo_as_0130_qed_mem0.grid')
      call nninitpdf(0)
      asmz = 0.130d0
    case ('nn23nlo')
      call nnpdfdriver('NNPDF23nlo_as_0119_qed_mem0.grid')
      call nninitpdf(0)
      asmz = 0.119d0
    case default
      write (*, *) 'Unsupported bundled PDF label: ', trim(pdlabel)
      write (*, *) 'Use nn23lo, nn23lo1, nn23nlo, or LHAPDF 6'
      stop 1
    end select
  end subroutine configure_pdf


  subroutine numberpdfm_impl(idummy)
    implicit none
    integer, intent(in) :: idummy

    write (*, *) 'ERROR: numberPDFm requires the LHAPDF 6 backend', idummy
    stop 1
  end subroutine numberpdfm_impl


  subroutine initpdfm_impl(idummy1, idummy2)
    implicit none
    integer, intent(in) :: idummy1, idummy2

    write (*, *) 'ERROR: initPDFm requires the LHAPDF 6 backend', &
                 idummy1, idummy2
    stop 1
  end subroutine initpdfm_impl


  subroutine initpdfsetbynamem_impl(idummy, cdummy)
    implicit none
    integer, intent(in) :: idummy
    character(len=*), intent(in) :: cdummy

    write (*, *) 'ERROR: initPDFsetbynamem requires the LHAPDF 6 backend', &
                 idummy, trim(cdummy)
    stop 1
  end subroutine initpdfsetbynamem_impl

end module pdfwrap_module

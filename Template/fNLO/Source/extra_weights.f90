module extra_weights
  implicit none
  private

  public :: maxscales, maxPDFs, maxPDFsets, maxdynscales
  public :: doreweight, lscalevar, lpdfvar
  public :: lhaPDFid, nmemPDF, dyn_scale
  public :: QCD_power, orders_tag, amp_pos
  public :: scalevarR, scalevarF
  public :: LHAPDFsetname
  public :: initialize_extra_weights

  integer, parameter :: maxscales = 9
  integer, parameter :: maxPDFs = 200
  integer, parameter :: maxPDFsets = 25
  integer, parameter :: maxdynscales = 10

  logical :: doreweight
  logical :: lscalevar(maxdynscales), lpdfvar(maxPDFsets)
  integer :: lhaPDFid(0:maxPDFsets), nmemPDF(maxPDFsets)
  integer :: dyn_scale(0:maxdynscales)
  integer :: QCD_power, orders_tag, amp_pos
  double precision :: scalevarR(0:maxscales), scalevarF(0:maxscales)
  character(len=100) :: LHAPDFsetname(maxPDFsets)

  save

contains

  subroutine initialize_extra_weights()
    implicit none

    lscalevar = .false.
    lpdfvar = .false.
    lhaPDFid = 0
    nmemPDF = 0
    dyn_scale = 0
    scalevarR = 0d0
    scalevarF = 0d0
    LHAPDFsetname = ''
  end subroutine initialize_extra_weights

end module extra_weights

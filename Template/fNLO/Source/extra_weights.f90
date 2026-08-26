module extra_weights
  implicit none
  private

  public :: iwgtinfo, maxscales, maxPDFs, maxPDFsets, maxdynscales
  public :: maxorders, doreweight, lscalevar, lpdfvar
  public :: iwgtnumpartn, jwgtinfo, mexternal, lhaPDFid, nmemPDF
  public :: dyn_scale, n_ctr_found, n_mom_conf, QCD_power, orders_tag
  public :: amp_pos, wgtdegrem_xi, wgtdegrem_lxi, wgtdegrem_muF
  public :: wgtnstmp, wgtwnstmpmuf, wgtwnstmpmur, wgtnstmp_avgvirt
  public :: wgtref, scalevarR, scalevarF, wgtxsecmu, wgtxsecPDF
  public :: LHAPDFsetname
  public :: initialize_extra_weights

  integer, parameter :: iwgtinfo = -5
  integer, parameter :: maxscales = 9
  integer, parameter :: maxPDFs = 200
  integer, parameter :: maxPDFsets = 25
  integer, parameter :: maxdynscales = 10
  integer, parameter :: maxorders = 10

  logical :: doreweight
  logical :: lscalevar(maxdynscales), lpdfvar(maxPDFsets)
  integer :: iwgtnumpartn, jwgtinfo, mexternal
  integer :: lhaPDFid(0:maxPDFsets), nmemPDF(maxPDFsets)
  integer :: dyn_scale(0:maxdynscales)
  integer :: n_ctr_found, n_mom_conf
  integer :: QCD_power, orders_tag, amp_pos
  double precision :: wgtdegrem_xi, wgtdegrem_lxi, wgtdegrem_muF
  double precision :: wgtnstmp, wgtwnstmpmuf, wgtwnstmpmur
  double precision :: wgtnstmp_avgvirt, wgtref
  double precision :: scalevarR(0:maxscales), scalevarF(0:maxscales)
  double precision :: wgtxsecmu(0:maxorders, maxscales, maxscales, &
                                maxdynscales)
  double precision :: wgtxsecPDF(0:maxPDFs, maxPDFsets)
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
    wgtxsecmu = 0d0
    wgtxsecPDF = 0d0
    LHAPDFsetname = ''
  end subroutine initialize_extra_weights

end module extra_weights

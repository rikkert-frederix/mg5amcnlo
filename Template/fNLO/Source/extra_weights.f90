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
  public :: wgtbpower, wgtcpower, LHAPDFsetname
  public :: cpower_pos, runfac, ren_group_coeff_in
  public :: initialize_extra_weights, validate_extra_weights
  public :: finalize_extra_weights, extra_weights_are_initialized

  integer, parameter :: iwgtinfo = -5
  integer, parameter :: maxscales = 9
  integer, parameter :: maxPDFs = 200
  integer, parameter :: maxPDFsets = 25
  integer, parameter :: maxdynscales = 10
  integer, parameter :: maxorders = 10

  logical :: doreweight
  logical, allocatable :: lscalevar(:), lpdfvar(:)
  integer :: iwgtnumpartn, jwgtinfo, mexternal
  integer, allocatable :: lhaPDFid(:), nmemPDF(:)
  integer, allocatable :: dyn_scale(:)
  integer :: n_ctr_found, n_mom_conf
  integer :: QCD_power, orders_tag, amp_pos
  double precision :: wgtdegrem_xi, wgtdegrem_lxi, wgtdegrem_muF
  double precision :: wgtnstmp, wgtwnstmpmuf, wgtwnstmpmur
  double precision :: wgtnstmp_avgvirt, wgtref
  double precision, allocatable :: scalevarR(:), scalevarF(:)
  double precision, allocatable :: wgtxsecmu(:,:,:,:)
  double precision, allocatable :: wgtxsecPDF(:,:)
  double precision :: wgtbpower, wgtcpower
  character(len=100), allocatable :: LHAPDFsetname(:)

  logical :: extra_weights_initialized = .false.

  ! Position of cpower.
  integer, parameter :: cpower_pos = 0

  ! Switch for a running (1) or fixed (0) muR-dependent factor.
  integer, parameter :: runfac = 0

  ! If runfac is enabled, rwgt_muR_dep_fac in reweight_xsec.f and
  ! compute_cpower in fks_singular.f must include all factorizing
  ! muR-dependent overall factors other than alpha_s.
  !
  ! First-order coefficient of the renormalization-group equation of the
  ! muR-dependent factor.  For masses and Yukawa couplings this is
  ! gamma0 = 3/2*C_F.
  integer, parameter :: ren_group_coeff_in = 0

  save

contains

  subroutine initialize_extra_weights()
    implicit none

    if (extra_weights_initialized) then
      call validate_extra_weights()
      return
    end if

    allocate(lscalevar(maxdynscales), lpdfvar(maxPDFsets))
    allocate(lhaPDFid(0:maxPDFsets), nmemPDF(maxPDFsets))
    allocate(dyn_scale(0:maxdynscales))
    allocate(scalevarR(0:maxscales), scalevarF(0:maxscales))
    allocate(wgtxsecmu(0:maxorders, maxscales, maxscales, &
                       maxdynscales))
    allocate(wgtxsecPDF(0:maxPDFs, maxPDFsets))
    allocate(LHAPDFsetname(maxPDFsets))

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

    extra_weights_initialized = .true.
    call validate_extra_weights()
  end subroutine initialize_extra_weights


  subroutine validate_extra_weights()
    implicit none

    if (.not. extra_weights_are_initialized()) then
      write(*, *) 'extra_weights storage is not initialized correctly'
      stop 1
    end if
  end subroutine validate_extra_weights


  logical function extra_weights_are_initialized()
    implicit none

    extra_weights_are_initialized = extra_weights_initialized
    extra_weights_are_initialized = extra_weights_are_initialized .and. &
         allocated(lscalevar) .and. allocated(lpdfvar)
    extra_weights_are_initialized = extra_weights_are_initialized .and. &
         allocated(lhaPDFid) .and. allocated(nmemPDF)
    extra_weights_are_initialized = extra_weights_are_initialized .and. &
         allocated(dyn_scale)
    extra_weights_are_initialized = extra_weights_are_initialized .and. &
         allocated(scalevarR) .and. allocated(scalevarF)
    extra_weights_are_initialized = extra_weights_are_initialized .and. &
         allocated(wgtxsecmu) .and. allocated(wgtxsecPDF)
    extra_weights_are_initialized = extra_weights_are_initialized .and. &
         allocated(LHAPDFsetname)
  end function extra_weights_are_initialized


  subroutine finalize_extra_weights()
    implicit none

    if (allocated(lscalevar)) deallocate(lscalevar)
    if (allocated(lpdfvar)) deallocate(lpdfvar)
    if (allocated(lhaPDFid)) deallocate(lhaPDFid)
    if (allocated(nmemPDF)) deallocate(nmemPDF)
    if (allocated(dyn_scale)) deallocate(dyn_scale)
    if (allocated(scalevarR)) deallocate(scalevarR)
    if (allocated(scalevarF)) deallocate(scalevarF)
    if (allocated(wgtxsecmu)) deallocate(wgtxsecmu)
    if (allocated(wgtxsecPDF)) deallocate(wgtxsecPDF)
    if (allocated(LHAPDFsetname)) deallocate(LHAPDFsetname)

    extra_weights_initialized = .false.
  end subroutine finalize_extra_weights

end module extra_weights

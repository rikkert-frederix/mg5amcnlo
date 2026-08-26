module run_state
  implicit none
  private

  integer, parameter, public :: run_seed_kind = selected_int_kind(18)
  integer, parameter, public :: max_lhe_processes = 100

  public :: reset_run_state
  public :: scale, fixed_ren_scale, fixed_fac_scale
  public :: dynamical_scale_choice
  public :: fixed_qes_scale, mur_over_ref, muf1_over_ref, muf2_over_ref
  public :: qes_over_ref, mur_ref_fixed
  public :: muf1_ref_fixed, muf2_ref_fixed, qes_ref_fixed
  public :: mur2_current, muf12_current, muf22_current, qes2_current
  public :: lpp, ebeam, xbk, q2fact, bwcutoff
  public :: do_rwgt_scale, do_rwgt_pdf
  public :: pdg_cut, ptmin4pdg, ptmax4pdg, mxxmin4pdg
  public :: mxxpart_antipart
  public :: jetalgo, jetradius, ptj, etaj, ptl, etal, drll, drll_sf
  public :: mll, mll_sf, gamma_is_j, ptgmin, r0gamma, xn, epsgamma
  public :: etagamma, isoem, maxjetflavor
  public :: lhaid, pdfscheme, pdlabel
  public :: asmz, nloop, iseed, pineappl
  public :: idbmup, ebmup, pdfgup, pdfsup, idwtup, nprup
  public :: xsecup, xerrup, xmaxup, lprup

  double precision :: scale
  logical :: fixed_ren_scale
  logical :: fixed_fac_scale
  integer :: dynamical_scale_choice

  logical :: fixed_qes_scale

  double precision :: mur_over_ref
  double precision :: muf1_over_ref
  double precision :: muf2_over_ref
  double precision :: qes_over_ref

  double precision :: mur_ref_fixed
  double precision :: muf1_ref_fixed
  double precision :: muf2_ref_fixed
  double precision :: qes_ref_fixed

  double precision :: mur2_current
  double precision :: muf12_current
  double precision :: muf22_current
  double precision :: qes2_current

  integer :: lpp(2)
  double precision :: ebeam(2)
  double precision :: xbk(2)
  double precision :: q2fact(2)
  common /to_collider/ ebeam, xbk, q2fact, lpp

  double precision :: bwcutoff

  logical :: do_rwgt_scale
  logical :: do_rwgt_pdf

  integer :: pdg_cut(0:25)
  double precision :: ptmin4pdg(0:25)
  double precision :: ptmax4pdg(0:25)
  double precision :: mxxmin4pdg(0:25)
  logical :: mxxpart_antipart(1:25)

  double precision :: jetalgo
  double precision :: jetradius
  double precision :: ptj
  double precision :: etaj
  double precision :: ptl
  double precision :: etal
  double precision :: drll
  double precision :: drll_sf
  double precision :: mll
  double precision :: mll_sf

  logical :: gamma_is_j

  double precision :: ptgmin
  double precision :: r0gamma
  double precision :: xn
  double precision :: epsgamma
  double precision :: etagamma
  logical :: isoem

  integer :: maxjetflavor

  integer :: lhaid
  integer :: pdfscheme
  character(len=7) :: pdlabel
  common /to_pdf/ lhaid, pdfscheme, pdlabel

  double precision :: asmz
  integer :: nloop
  common /a_block/ asmz, nloop

  integer(kind=run_seed_kind) :: iseed
  common /to_seed/ iseed

  logical :: pineappl
  common /for_pineappl/ pineappl

  integer :: idbmup(2)
  double precision :: ebmup(2)
  integer :: pdfgup(2)
  integer :: pdfsup(2)
  integer :: idwtup
  integer :: nprup
  double precision :: xsecup(max_lhe_processes)
  double precision :: xerrup(max_lhe_processes)
  double precision :: xmaxup(max_lhe_processes)
  integer :: lprup(max_lhe_processes)
  common /heprup/ idbmup, ebmup, pdfgup, pdfsup, idwtup, nprup, &
       xsecup, xerrup, xmaxup, lprup

contains

  subroutine reset_run_state()
    implicit none

    scale = 0d0
    fixed_ren_scale = .false.
    fixed_fac_scale = .false.
    dynamical_scale_choice = 0
    fixed_qes_scale = .false.
    mur_over_ref = 0d0
    muf1_over_ref = 0d0
    muf2_over_ref = 0d0
    qes_over_ref = 0d0
    mur_ref_fixed = 0d0
    muf1_ref_fixed = 0d0
    muf2_ref_fixed = 0d0
    qes_ref_fixed = 0d0
    mur2_current = 0d0
    muf12_current = 0d0
    muf22_current = 0d0
    qes2_current = 0d0
    lpp = 0
    ebeam = 0d0
    xbk = 0d0
    q2fact = 0d0
    bwcutoff = 0d0
    do_rwgt_scale = .false.
    do_rwgt_pdf = .false.
    pdg_cut = 0
    ptmin4pdg = 0d0
    ptmax4pdg = 0d0
    mxxmin4pdg = 0d0
    mxxpart_antipart = .false.
    jetalgo = 0d0
    jetradius = 0d0
    ptj = 0d0
    etaj = 0d0
    ptl = 0d0
    etal = 0d0
    drll = 0d0
    drll_sf = 0d0
    mll = 0d0
    mll_sf = 0d0
    gamma_is_j = .false.
    ptgmin = 0d0
    r0gamma = 0d0
    xn = 0d0
    epsgamma = 0d0
    etagamma = 0d0
    isoem = .false.
    maxjetflavor = 0
    lhaid = 0
    pdfscheme = 0
    pdlabel = ''
    asmz = 0d0
    nloop = 0
    iseed = 0_run_seed_kind
    pineappl = .false.
    idbmup = 0
    ebmup = 0d0
    pdfgup = 0
    pdfsup = 0
    idwtup = 0
    nprup = 0
    xsecup = 0d0
    xerrup = 0d0
    xmaxup = 0d0
    lprup = 0
  end subroutine reset_run_state

end module run_state

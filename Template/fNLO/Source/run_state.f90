module run_state
  implicit none
  private

  integer, parameter, public :: run_seed_kind = selected_int_kind(18)
  integer, parameter, public :: max_lhe_processes = 100

  public :: reset_run_state
  public :: scale, scalefact, ellissextonfact, alpsfact
  public :: mue_ref_fixed, mue_over_ref, fixed_ren_scale
  public :: fixed_fac_scale, fixed_couplings, fixed_extra_scale
  public :: ickkw, nhmult, hmult, dynamical_scale_choice
  public :: fixed_qes_scale, mur_over_ref, muf1_over_ref, muf2_over_ref
  public :: qes_over_ref, muf_over_ref, mur_ref_fixed
  public :: muf1_ref_fixed, muf2_ref_fixed, qes_ref_fixed, muf_ref_fixed
  public :: mur2_current, muf12_current, muf22_current, qes2_current
  public :: lpp, ebeam, xbk, q2fact, bwcutoff
  public :: ktscheme, chcluster, pdfwgt
  public :: do_rwgt_scale, do_rwgt_pdf, store_rwgt_info
  public :: pdf_set_min, pdf_set_max, rw_fscale_down, rw_fscale_up
  public :: rw_rscale_down, rw_rscale_up, fo_lhe_weight_ratio
  public :: pdg_cut, ptmin4pdg, ptmax4pdg, mxxmin4pdg
  public :: mxxpart_antipart, alphascheme, nlep_run, nupq_run
  public :: ndnq_run, w_run, photons_from_lepton, has_bstrahl
  public :: nb_proton, nb_neutron, mass_ion, rai, xiai, flavour_bias
  public :: jetalgo, jetradius, ptj, etaj, ptl, etal, drll, drll_sf
  public :: mll, mll_sf, gamma_is_j, ptgmin, r0gamma, xn, epsgamma
  public :: etagamma, isoem, maxjetflavor, rphreco, etaphreco
  public :: lepphreco, quarkphreco, xmtc, xqcut
  public :: lhaid, pdfscheme, pdlabel, epa_label, pdsublabel
  public :: asmz, nloop, iseed, nevents, event_norm, pineappl
  public :: jet_distance_parameter
  public :: idbmup, ebmup, pdfgup, pdfsup, idwtup, nprup
  public :: xsecup, xerrup, xmaxup, lprup

  double precision :: scale
  double precision :: scalefact
  double precision :: ellissextonfact
  double precision :: alpsfact
  double precision :: mue_ref_fixed
  double precision :: mue_over_ref
  logical :: fixed_ren_scale
  logical :: fixed_fac_scale
  logical :: fixed_couplings
  logical :: fixed_extra_scale
  integer :: ickkw
  integer :: nhmult
  logical :: hmult
  integer :: dynamical_scale_choice
  common /to_scale/ scale, scalefact, ellissextonfact, alpsfact, &
       fixed_ren_scale, mue_ref_fixed, mue_over_ref, fixed_extra_scale, &
       fixed_fac_scale, fixed_couplings, ickkw, nhmult, hmult, &
       dynamical_scale_choice

  logical :: fixed_qes_scale
  common /cfxqes/ fixed_qes_scale

  double precision :: mur_over_ref
  double precision :: muf1_over_ref
  double precision :: muf2_over_ref
  double precision :: qes_over_ref
  double precision :: muf_over_ref
  common /cscales_fact/ mur_over_ref, muf1_over_ref, muf2_over_ref, &
       qes_over_ref, muf_over_ref

  double precision :: mur_ref_fixed
  double precision :: muf1_ref_fixed
  double precision :: muf2_ref_fixed
  double precision :: qes_ref_fixed
  double precision :: muf_ref_fixed
  common /cscales_fixed_values/ mur_ref_fixed, muf1_ref_fixed, &
       muf2_ref_fixed, qes_ref_fixed, muf_ref_fixed

  double precision :: mur2_current
  double precision :: muf12_current
  double precision :: muf22_current
  double precision :: qes2_current
  common /cscales_current_values/ mur2_current, muf12_current, &
       muf22_current, qes2_current

  integer :: lpp(2)
  double precision :: ebeam(2)
  double precision :: xbk(2)
  double precision :: q2fact(2)
  common /to_collider/ ebeam, xbk, q2fact, lpp

  double precision :: bwcutoff
  common /to_bwcutoff/ bwcutoff

  integer :: ktscheme
  logical :: chcluster
  logical :: pdfwgt
  common /to_cluster/ ktscheme, chcluster, pdfwgt

  logical :: do_rwgt_scale
  logical :: do_rwgt_pdf
  logical :: store_rwgt_info
  integer :: pdf_set_min
  integer :: pdf_set_max
  double precision :: rw_fscale_down
  double precision :: rw_fscale_up
  double precision :: rw_rscale_down
  double precision :: rw_rscale_up
  common /to_rwgt/ do_rwgt_scale, rw_fscale_down, rw_fscale_up, &
       rw_rscale_down, rw_rscale_up, do_rwgt_pdf, pdf_set_min, &
       pdf_set_max, store_rwgt_info

  double precision :: fo_lhe_weight_ratio

  integer :: pdg_cut(0:25)
  double precision :: ptmin4pdg(0:25)
  double precision :: ptmax4pdg(0:25)
  double precision :: mxxmin4pdg(0:25)
  logical :: mxxpart_antipart(1:25)
  common /to_pdg_specific_cut/ pdg_cut, ptmin4pdg, ptmax4pdg, &
       mxxmin4pdg, mxxpart_antipart

  integer :: alphascheme
  integer :: nlep_run
  integer :: nupq_run
  integer :: ndnq_run
  integer :: w_run
  common /to_alphascheme/ alphascheme, nlep_run, nupq_run, ndnq_run, w_run

  logical :: photons_from_lepton
  logical :: has_bstrahl
  common /to_afromee/ photons_from_lepton
  common /to_has_bs/ has_bstrahl

  integer :: nb_proton(2)
  integer :: nb_neutron(2)
  double precision :: mass_ion(2)
  double precision :: rai(2)
  double precision :: xiai(2)
  common /to_heavyion_pdg/ nb_proton, nb_neutron
  common /to_heavyion_mass/ mass_ion
  common /heavyion_upc_beams/ rai, xiai

  integer :: flavour_bias(0:2)
  common /c_flavour_bias/ flavour_bias

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
  common /to_new_auto_cuts/ jetalgo, jetradius, ptj, etaj, ptl, etal, &
       drll, drll_sf, mll, mll_sf

  logical :: gamma_is_j
  common /to_gamma_is_j/ gamma_is_j

  double precision :: ptgmin
  double precision :: r0gamma
  double precision :: xn
  double precision :: epsgamma
  double precision :: etagamma
  logical :: isoem
  common /to_isogamma_cuts/ ptgmin, r0gamma, xn, epsgamma, etagamma
  common /to_isogamma_em/ isoem

  integer :: maxjetflavor
  common /to_min_max_cuts/ maxjetflavor

  double precision :: rphreco
  double precision :: etaphreco
  logical :: lepphreco
  logical :: quarkphreco
  common /to_phreco/ rphreco, etaphreco, lepphreco, quarkphreco

  double precision :: xmtc
  double precision :: xqcut
  common /to_specxpt/ xmtc, xqcut

  integer :: lhaid
  integer :: pdfscheme
  character(len=7) :: pdlabel
  character(len=7) :: epa_label
  character(len=7) :: pdsublabel(2)
  common /to_pdf/ lhaid, pdfscheme, pdlabel, epa_label, pdsublabel

  double precision :: asmz
  integer :: nloop
  common /a_block/ asmz, nloop

  integer(kind=run_seed_kind) :: iseed
  common /to_seed/ iseed

  integer :: nevents
  character(len=7) :: event_norm
  common /event_normalisation/ event_norm

  logical :: pineappl
  common /for_pineappl/ pineappl

  double precision :: jet_distance_parameter
  common /to_dj/ jet_distance_parameter

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
    scalefact = 0d0
    ellissextonfact = 0d0
    alpsfact = 0d0
    mue_ref_fixed = 0d0
    mue_over_ref = 0d0
    fixed_ren_scale = .false.
    fixed_fac_scale = .false.
    fixed_couplings = .false.
    fixed_extra_scale = .false.
    ickkw = 0
    nhmult = 0
    hmult = .false.
    dynamical_scale_choice = 0
    fixed_qes_scale = .false.
    mur_over_ref = 0d0
    muf1_over_ref = 0d0
    muf2_over_ref = 0d0
    qes_over_ref = 0d0
    muf_over_ref = 0d0
    mur_ref_fixed = 0d0
    muf1_ref_fixed = 0d0
    muf2_ref_fixed = 0d0
    qes_ref_fixed = 0d0
    muf_ref_fixed = 0d0
    mur2_current = 0d0
    muf12_current = 0d0
    muf22_current = 0d0
    qes2_current = 0d0
    lpp = 0
    ebeam = 0d0
    xbk = 0d0
    q2fact = 0d0
    bwcutoff = 0d0
    ktscheme = 0
    chcluster = .false.
    pdfwgt = .false.
    do_rwgt_scale = .false.
    do_rwgt_pdf = .false.
    store_rwgt_info = .false.
    pdf_set_min = 0
    pdf_set_max = 0
    rw_fscale_down = 0d0
    rw_fscale_up = 0d0
    rw_rscale_down = 0d0
    rw_rscale_up = 0d0
    fo_lhe_weight_ratio = 0d0
    pdg_cut = 0
    ptmin4pdg = 0d0
    ptmax4pdg = 0d0
    mxxmin4pdg = 0d0
    mxxpart_antipart = .false.
    alphascheme = 0
    nlep_run = 0
    nupq_run = 0
    ndnq_run = 0
    w_run = 0
    photons_from_lepton = .false.
    has_bstrahl = .false.
    nb_proton = 0
    nb_neutron = 0
    mass_ion = 0d0
    rai = 0d0
    xiai = 0d0
    flavour_bias = 0
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
    rphreco = 0d0
    etaphreco = 0d0
    lepphreco = .false.
    quarkphreco = .false.
    xmtc = 0d0
    xqcut = 0d0
    lhaid = 0
    pdfscheme = 0
    pdlabel = ''
    epa_label = ''
    pdsublabel = ''
    asmz = 0d0
    nloop = 0
    iseed = 0_run_seed_kind
    nevents = 0
    event_norm = ''
    pineappl = .false.
    jet_distance_parameter = 0d0
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

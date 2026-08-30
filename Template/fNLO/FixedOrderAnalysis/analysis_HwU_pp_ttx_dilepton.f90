module analysis_hwu_pp_ttx_dilepton_module
  use process_dimensions, only: nincoming
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private

  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: tiny = 1d-12
  double precision, parameter :: jet_radius = 0.4d0
  double precision, parameter :: lepton_pt_lead_cut = 25d0
  double precision, parameter :: lepton_pt_trail_cut = 20d0
  double precision, parameter :: lepton_eta_cut = 2.5d0
  double precision, parameter :: bjet_pt_cut = 30d0
  double precision, parameter :: bjet_eta_cut = 2.5d0
  double precision, parameter :: analysis_jet_pt_cut = 25d0
  double precision, parameter :: analysis_jet_eta_cut = 4.5d0
  double precision, parameter :: lepton_jet_dr_cut = 0.4d0
  double precision, parameter :: dilepton_mass_cut = 20d0
  double precision, parameter :: missing_pt_cut = 20d0

  integer, parameter :: g_rate = 0
  integer, parameter :: g_top = 2
  integer, parameter :: g_lep = 28
  integer, parameter :: g_decay = 52
  integer, parameter :: g_neutrino = 76
  integer, parameter :: g_jet = 91
  integer, parameter :: g_fiducial = 121
  integer, parameter :: g_reco = 159
  integer, parameter :: g_spin = 190
  integer, parameter :: g_reco_spin = 216
  integer, parameter :: g_diagnostic = 242
  integer, parameter :: number_of_histograms = 247

  double precision :: top_mass_reference = 172.5d0
  double precision :: w_mass_reference = 80.379d0

  public :: analysis_begin, analysis_end, analysis_fill

  interface
    subroutine amcatnlo_fastjetppgenkt_etamax(pqcd, nn, rfj, sycut, &
                                               etamax, palg, pjet, njet, jet)
      implicit none
      integer, intent(in) :: nn
      double precision, intent(in) :: pqcd(0:3, nn)
      double precision, intent(in) :: rfj, sycut, etamax, palg
      double precision, intent(out) :: pjet(0:3, nn)
      integer, intent(out) :: njet, jet(nn)
    end subroutine amcatnlo_fastjetppgenkt_etamax

    double precision function get_mass_from_id(id)
      integer, intent(in) :: id
    end function get_mass_from_id
  end interface

contains

  subroutine analysis_begin(nwgt, weights_info)
    implicit none
    integer, intent(in) :: nwgt
    character(len=*), intent(in) :: weights_info(*)

    call HwU_inithist(nwgt, weights_info)
    top_mass_reference = abs(get_mass_from_id(6))
    w_mass_reference = abs(get_mass_from_id(24))
    if (top_mass_reference <= 0d0) top_mass_reference = 172.5d0
    if (w_mass_reference <= 0d0) w_mass_reference = 80.379d0
    call book_histograms()
  end subroutine analysis_begin


  subroutine analysis_end()
    use open_output_files_module, only: HwU_write_file
    implicit none
    call HwU_write_file
  end subroutine analysis_end


  subroutine analysis_fill(p, istatus, ipdg, wgts, ibody)
    implicit none
    double precision, intent(in) :: p(0:, :)
    integer, intent(in) :: istatus(:), ipdg(:)
    double precision, intent(in) :: wgts(*)
    integer, intent(in) :: ibody

    integer :: particle_count
    integer :: ilp, ilm, inu, ianu, ib, iab
    integer :: ibjet, iabjet, njet, njet25, njet40, nextra25
    integer :: lead_jet, trail_jet
    integer :: lead_extra_jet, assigned_radiation, pairing, i, stage
    integer :: jet_of_parton(size(p, 2))
    double precision :: pjet(0:3, size(p, 2))
    double precision :: lp(0:3), lm(0:3), nu(0:3), anu(0:3)
    double precision :: b(0:3), ab(0:3), wp(0:3), wm(0:3), ww(0:3)
    double precision :: top(0:3), atop(0:3), ttbar(0:3)
    double precision :: bjet(0:3), abjet(0:3), bbjet(0:3)
    double precision :: ll(0:3), visible(0:3), met(0:1)
    double precision :: reco_nu(0:3), reco_anu(0:3)
    double precision :: reco_wp(0:3), reco_wm(0:3)
    double precision :: reco_top(0:3), reco_atop(0:3), reco_ttbar(0:3)
    double precision :: truth_score, reco_score, mt2_ll, mt2_blbl
    double precision :: mlb_min, mlb_max, mlb_average
    double precision :: min_dr_lepton_jet, min_dr_lepton_b, max_dr_lepton_b
    double precision :: leading_extra_pt, leading_extra_y, extra_ht
    double precision :: jet_ht, visible_ht, st, threshold
    double precision :: lp_pt, lm_pt, lp_eta, lm_eta
    double precision :: bjet_pt, abjet_pt, bjet_eta, abjet_eta
    double precision :: spin_values(6), spin_products(9), spin_extra(11)
    double precision :: reco_spin_values(6), reco_spin_products(9)
    double precision :: reco_spin_extra(11)
    double precision :: wp_helicity, wm_helicity
    double precision :: reco_wp_helicity, reco_wm_helicity
    logical :: fiducial, reco_success, same_flavour_signal_region

    particle_count = size(p, 2)
    if (size(istatus) /= particle_count .or. size(ipdg) /= particle_count) then
      call fail_analysis('inconsistent momentum, status and PDG arrays')
    end if
    if (particle_count < nincoming + 6) then
      call fail_analysis('too few particles for a dileptonic ttbar event')
    end if

    call find_dilepton_objects(p, istatus, ipdg, ilp, ilm, inu, ianu, ib, iab)
    lp = p(0:3, ilp)
    lm = p(0:3, ilm)
    nu = p(0:3, inu)
    anu = p(0:3, ianu)
    b = p(0:3, ib)
    ab = p(0:3, iab)
    wp = lp + nu
    wm = lm + anu
    ww = wp + wm
    ll = lp + lm

    call reconstruct_truth_tops(p, istatus, ipdg, wp, wm, b, ab, &
         top, atop, truth_score, assigned_radiation)
    ttbar = top + atop

    call cluster_qcd_jets(p, istatus, ipdg, pjet, njet, jet_of_parton, &
         ibjet, iabjet)
    if (ibjet <= 0 .or. iabjet <= 0) then
      call fail_analysis('FastJet did not retain both bottom-flavoured jets')
    end if
    bjet = pjet(:, ibjet)
    abjet = pjet(:, iabjet)
    if (ibjet == iabjet) then
      bbjet = bjet
    else
      bbjet = bjet + abjet
    end if
    visible = ll + bbjet

    call missing_transverse_momentum(p, istatus, ipdg, met)
    call classify_jets(pjet, njet, ibjet, iabjet, njet25, njet40, &
         nextra25, lead_jet, trail_jet, lead_extra_jet, jet_ht, extra_ht)
    leading_extra_pt = 0d0
    leading_extra_y = 0d0
    if (lead_extra_jet > 0) then
      leading_extra_pt = transverse_momentum(pjet(:, lead_extra_jet))
      leading_extra_y = rapidity(pjet(:, lead_extra_jet))
    end if
    call lepton_b_pair_masses(lp, lm, bjet, abjet, &
         mlb_min, mlb_max, mlb_average)
    call lepton_jet_distances(lp, lm, pjet, njet, &
         min_dr_lepton_jet, min_dr_lepton_b, max_dr_lepton_b, &
         ibjet, iabjet)

    visible_ht = transverse_momentum(lp) + transverse_momentum(lm) + jet_ht
    st = visible_ht + vector_magnitude_2(met)

    call compute_spin_observables(top, atop, lp, lm, &
         spin_values, spin_products, spin_extra)
    call w_helicity_cosines(top, atop, wp, wm, lp, lm, &
         wp_helicity, wm_helicity)

    call fill_inclusive_truth(wgts, top, atop, ttbar, lp, lm, ll, &
         wp, wm, ww, nu, anu, bjet, abjet, bbjet, met, pjet, &
         njet25, njet40, nextra25, lead_jet, trail_jet, &
         leading_extra_pt, leading_extra_y, extra_ht, jet_ht, &
         visible_ht, st, visible, min_dr_lepton_b, max_dr_lepton_b, &
         spin_values, spin_products, spin_extra, wp_helicity, &
         wm_helicity)

    call HwU_fill(g_rate + 1, 1d0, wgts)
    call HwU_fill(g_rate + 2, 1d0, wgts)

    lp_pt = transverse_momentum(lp)
    lm_pt = transverse_momentum(lm)
    lp_eta = abs(pseudorapidity(lp))
    lm_eta = abs(pseudorapidity(lm))
    bjet_pt = transverse_momentum(bjet)
    abjet_pt = transverse_momentum(abjet)
    bjet_eta = abs(pseudorapidity(bjet))
    abjet_eta = abs(pseudorapidity(abjet))

    stage = 1
    if (max(lp_pt, lm_pt) > &
          lepton_pt_lead_cut .and. &
        min(lp_pt, lm_pt) > &
          lepton_pt_trail_cut .and. &
        lp_eta < lepton_eta_cut .and. lm_eta < lepton_eta_cut) then
      stage = 2
      call HwU_fill(g_rate + 2, 2d0, wgts)
    end if
    if (stage == 2 .and. invariant_mass(ll) > dilepton_mass_cut) then
      stage = 3
      call HwU_fill(g_rate + 2, 3d0, wgts)
    end if
    if (stage == 3 .and. ibjet /= iabjet .and. &
        bjet_pt > bjet_pt_cut .and. abjet_pt > bjet_pt_cut .and. &
        bjet_eta < bjet_eta_cut .and. abjet_eta < bjet_eta_cut) then
      stage = 4
      call HwU_fill(g_rate + 2, 4d0, wgts)
    end if
    if (stage == 4 .and. min_dr_lepton_jet > lepton_jet_dr_cut) then
      stage = 5
      call HwU_fill(g_rate + 2, 5d0, wgts)
    end if
    if (stage == 5 .and. vector_magnitude_2(met) > missing_pt_cut) then
      stage = 6
      call HwU_fill(g_rate + 2, 6d0, wgts)
    end if
    fiducial = stage == 6

    same_flavour_signal_region = .true.
    if (abs(ipdg(ilp)) == abs(ipdg(ilm))) then
      same_flavour_signal_region = &
           abs(invariant_mass(ll) - 91.1876d0) > 15d0 .and. &
           vector_magnitude_2(met) > 40d0
    end if
    if (fiducial .and. same_flavour_signal_region) then
      call HwU_fill(g_rate + 2, 7d0, wgts)
    end if

    reco_success = .false.
    reco_score = 0d0
    pairing = 0
    reco_nu = 0d0
    reco_anu = 0d0
    reco_wp = 0d0
    reco_wm = 0d0
    reco_top = 0d0
    reco_atop = 0d0
    reco_ttbar = 0d0
    mt2_ll = 0d0
    mt2_blbl = 0d0
    if (fiducial) then
      call reconstruct_neutrinos(lp, lm, bjet, abjet, met, &
           reco_nu, reco_anu, pairing, reco_score, reco_success)
      call compute_mt2(lp, lm, met, mt2_ll)
      call compute_mt2_blbl(lp, lm, bjet, abjet, met, mt2_blbl)
      if (reco_success) then
        reco_wp = lp + reco_nu
        reco_wm = lm + reco_anu
        if (pairing == 1) then
          reco_top = reco_wp + bjet
          reco_atop = reco_wm + abjet
        else
          reco_top = reco_wp + abjet
          reco_atop = reco_wm + bjet
        end if
        reco_ttbar = reco_top + reco_atop
        call compute_spin_observables(reco_top, reco_atop, lp, lm, &
             reco_spin_values, reco_spin_products, reco_spin_extra)
        call w_helicity_cosines(reco_top, reco_atop, reco_wp, reco_wm, &
             lp, lm, reco_wp_helicity, reco_wm_helicity)
        if (same_flavour_signal_region) &
             call HwU_fill(g_rate + 2, 8d0, wgts)
      end if
      call fill_fiducial(wgts, lp, lm, ll, bjet, abjet, bbjet, &
           met, visible, jet_ht, visible_ht, st, njet25, nextra25, &
           mlb_min, mlb_max, mlb_average, mt2_ll, mt2_blbl)
      if (reco_success) then
        call fill_reconstructed(wgts, reco_nu, reco_anu, reco_wp, reco_wm, &
             reco_top, reco_atop, reco_ttbar, top, atop, ttbar, &
             reco_score, pairing, reco_wp_helicity, reco_wm_helicity, &
             reco_spin_values, reco_spin_products, reco_spin_extra)
      end if
    end if

    do i = 1, 20
      threshold = 15d0 + 10d0*dble(i - 1)
      if (leading_extra_pt < threshold) &
           call HwU_fill(g_jet + 27, threshold, wgts)
      if (extra_ht < threshold) &
           call HwU_fill(g_jet + 28, threshold, wgts)
    end do

    call HwU_fill(g_diagnostic + 1, min(truth_score, 99.999d0), wgts)
    call HwU_fill(g_diagnostic + 2, dble(assigned_radiation), wgts)
    call HwU_fill(g_diagnostic + 3, &
         min(max(abs(invariant_mass(top) - top_mass_reference), &
                 abs(invariant_mass(atop) - top_mass_reference)), &
             49.999d0), wgts)
    if (reco_success) then
      call HwU_fill(g_diagnostic + 4, &
           min(max(abs(invariant_mass(reco_wp) - w_mass_reference), &
                   abs(invariant_mass(reco_wm) - w_mass_reference)), &
               49.999d0), wgts)
      call HwU_fill(g_diagnostic + 5, 1d0, wgts)
    end if

    ! ibody is deliberately not used to split these distributions.  Every
    ! real event and projected counterevent must be analysed with the same
    ! IR-safe measurement function.
    if (ibody < 1 .or. ibody > 3) then
      call fail_analysis('invalid fixed-order body type')
    end if
  end subroutine analysis_fill


  subroutine book_histograms()
    implicit none

    call book(g_rate + 1, 'inclusive rate', 1, 0.5d0, 1.5d0)
    call book(g_rate + 2, 'fiducial cutflow', 8, 0.5d0, 8.5d0)

    call book(g_top + 1, 'truth top pT', 50, 0d0, 1000d0)
    call book(g_top + 2, 'truth antitop pT', 50, 0d0, 1000d0)
    call book(g_top + 3, 'truth leading top pT', 50, 0d0, 1000d0)
    call book(g_top + 4, 'truth trailing top pT', 50, 0d0, 600d0)
    call book(g_top + 5, 'truth top rapidity', 50, -5d0, 5d0)
    call book(g_top + 6, 'truth antitop rapidity', 50, -5d0, 5d0)
    call book(g_top + 7, 'truth ttbar mass', 54, 300d0, 3000d0)
    call book(g_top + 8, 'truth ttbar pT', 50, 0d0, 1000d0)
    call book(g_top + 9, 'truth ttbar rapidity', 50, -5d0, 5d0)
    call book(g_top + 10, 'truth abs ttbar rapidity', 25, 0d0, 5d0)
    call book(g_top + 11, 'truth delta rapidity top antitop', 50, -5d0, 5d0)
    call book(g_top + 12, 'truth abs delta rapidity top antitop', 25, 0d0, 5d0)
    call book(g_top + 13, 'truth delta abs rapidity top antitop', 50, -5d0, 5d0)
    call book(g_top + 14, 'truth ttbar boost rapidity', 50, -5d0, 5d0)
    call book(g_top + 15, 'truth ttbar star rapidity', 50, -3d0, 3d0)
    call book(g_top + 16, 'truth top scattering cosine', 40, -1d0, 1d0)
    call book(g_top + 17, 'truth top speed in ttbar frame', 40, 0d0, 1d0)
    call book(g_top + 18, 'truth rho two mt over mtt', 40, 0d0, 1d0)
    call book(g_top + 19, 'truth delta phi top antitop', 32, 0d0, pi)
    call book(g_top + 20, 'truth delta R top antitop', 40, 0d0, 8d0)
    call book(g_top + 21, 'truth scalar top pT sum', 50, 0d0, 1500d0)
    call book(g_top + 22, 'truth top pT asymmetry', 40, -1d0, 1d0)
    call book(g_top + 23, 'truth ttbar transverse centrality', 40, 0d0, 2d0)
    call book(g_top + 24, 'truth top energy in ttbar frame', 50, 150d0, 1500d0)
    call book(g_top + 25, 'truth chi exp abs delta y', 49, 1d0, 50d0)
    call book(g_top + 26, 'truth ttbar longitudinal boost', 40, 0d0, 1d0)

    call book(g_lep + 1, 'truth positive lepton pT', 50, 0d0, 500d0)
    call book(g_lep + 2, 'truth negative lepton pT', 50, 0d0, 500d0)
    call book(g_lep + 3, 'truth leading lepton pT', 50, 0d0, 500d0)
    call book(g_lep + 4, 'truth trailing lepton pT', 50, 0d0, 300d0)
    call book(g_lep + 5, 'truth positive lepton eta', 50, -5d0, 5d0)
    call book(g_lep + 6, 'truth negative lepton eta', 50, -5d0, 5d0)
    call book(g_lep + 7, 'truth positive lepton rapidity', 50, -5d0, 5d0)
    call book(g_lep + 8, 'truth negative lepton rapidity', 50, -5d0, 5d0)
    call book(g_lep + 9, 'truth dilepton mass', 50, 0d0, 1000d0)
    call book(g_lep + 10, 'truth dilepton pT', 50, 0d0, 600d0)
    call book(g_lep + 11, 'truth dilepton rapidity', 50, -5d0, 5d0)
    call book(g_lep + 12, 'truth abs dilepton rapidity', 25, 0d0, 5d0)
    call book(g_lep + 13, 'truth dilepton delta phi', 32, 0d0, pi)
    call book(g_lep + 14, 'truth dilepton delta eta', 50, -5d0, 5d0)
    call book(g_lep + 15, 'truth abs dilepton delta eta', 25, 0d0, 5d0)
    call book(g_lep + 16, 'truth dilepton delta R', 40, 0d0, 8d0)
    call book(g_lep + 17, 'truth lepton opening cosine lab', 40, -1d0, 1d0)
    call book(g_lep + 18, 'truth dilepton cos theta star', 40, -1d0, 1d0)
    call book(g_lep + 19, 'truth lepton pT ratio', 40, 0d0, 1d0)
    call book(g_lep + 20, 'truth scalar lepton pT sum', 50, 0d0, 800d0)
    call book(g_lep + 21, 'truth delta abs lepton eta', 50, -5d0, 5d0)
    call book(g_lep + 22, 'truth dilepton acoplanarity', 32, 0d0, pi)
    call book(g_lep + 23, 'truth signed dilepton delta phi', 64, -pi, pi)
    call book(g_lep + 24, 'truth dilepton mass over lepton HT', 40, 0d0, 5d0)

    call book(g_decay + 1, 'truth Wplus pT', 50, 0d0, 600d0)
    call book(g_decay + 2, 'truth Wplus rapidity', 50, -5d0, 5d0)
    call book(g_decay + 3, 'truth Wplus mass', 50, 0d0, 150d0)
    call book(g_decay + 4, 'truth Wminus pT', 50, 0d0, 600d0)
    call book(g_decay + 5, 'truth Wminus rapidity', 50, -5d0, 5d0)
    call book(g_decay + 6, 'truth Wminus mass', 50, 0d0, 150d0)
    call book(g_decay + 7, 'truth WW mass', 50, 0d0, 1500d0)
    call book(g_decay + 8, 'truth WW pT', 50, 0d0, 800d0)
    call book(g_decay + 9, 'truth WW rapidity', 50, -5d0, 5d0)
    call book(g_decay + 10, 'truth WW delta phi', 32, 0d0, pi)
    call book(g_decay + 11, 'truth WW delta R', 40, 0d0, 8d0)
    call book(g_decay + 12, 'truth positive lepton b mass', 50, 0d0, 250d0)
    call book(g_decay + 13, 'truth negative lepton bbar mass', 50, 0d0, 250d0)
    call book(g_decay + 14, 'truth minimum correct lb mass', 50, 0d0, 250d0)
    call book(g_decay + 15, 'truth maximum correct lb mass', 50, 0d0, 250d0)
    call book(g_decay + 16, 'truth average correct lb mass', 50, 0d0, 250d0)
    call book(g_decay + 17, 'truth positive lepton energy top frame', 50, 0d0, 100d0)
    call book(g_decay + 18, 'truth negative lepton energy top frame', 50, 0d0, 100d0)
    call book(g_decay + 19, 'truth top recoil energy top frame', 50, 0d0, 120d0)
    call book(g_decay + 20, 'truth antitop recoil energy top frame', 50, 0d0, 120d0)
    call book(g_decay + 21, 'truth Wplus helicity cosine', 40, -1d0, 1d0)
    call book(g_decay + 22, 'truth Wminus helicity cosine', 40, -1d0, 1d0)
    call book(g_decay + 23, 'truth Wplus energy top frame', 50, 50d0, 150d0)
    call book(g_decay + 24, 'truth Wminus energy top frame', 50, 50d0, 150d0)

    call book(g_neutrino + 1, 'truth neutrino pT', 50, 0d0, 500d0)
    call book(g_neutrino + 2, 'truth neutrino eta', 50, -7d0, 7d0)
    call book(g_neutrino + 3, 'truth neutrino rapidity', 50, -7d0, 7d0)
    call book(g_neutrino + 4, 'truth antineutrino pT', 50, 0d0, 500d0)
    call book(g_neutrino + 5, 'truth antineutrino eta', 50, -7d0, 7d0)
    call book(g_neutrino + 6, 'truth antineutrino rapidity', 50, -7d0, 7d0)
    call book(g_neutrino + 7, 'truth dineutrino pT', 50, 0d0, 600d0)
    call book(g_neutrino + 8, 'truth dineutrino mass', 50, 0d0, 1000d0)
    call book(g_neutrino + 9, 'truth dineutrino rapidity', 50, -7d0, 7d0)
    call book(g_neutrino + 10, 'truth dineutrino delta phi', 32, 0d0, pi)
    call book(g_neutrino + 11, 'truth dineutrino delta R', 40, 0d0, 10d0)
    call book(g_neutrino + 12, 'visible missing pT', 50, 0d0, 600d0)
    call book(g_neutrino + 13, 'minimum delta phi missing pT lepton', 32, 0d0, pi)
    call book(g_neutrino + 14, 'maximum delta phi missing pT lepton', 32, 0d0, pi)
    call book(g_neutrino + 15, 'dilepton missing transverse mass', 50, 0d0, 1000d0)

    call book(g_jet + 1, 'truth b jet pT', 50, 0d0, 600d0)
    call book(g_jet + 2, 'truth bbar jet pT', 50, 0d0, 600d0)
    call book(g_jet + 3, 'truth leading b jet pT', 50, 0d0, 600d0)
    call book(g_jet + 4, 'truth trailing b jet pT', 50, 0d0, 400d0)
    call book(g_jet + 5, 'truth b jet eta', 50, -5d0, 5d0)
    call book(g_jet + 6, 'truth bbar jet eta', 50, -5d0, 5d0)
    call book(g_jet + 7, 'truth bb jet mass', 50, 0d0, 1200d0)
    call book(g_jet + 8, 'truth bb jet pT', 50, 0d0, 800d0)
    call book(g_jet + 9, 'truth bb jet delta phi', 32, 0d0, pi)
    call book(g_jet + 10, 'truth bb jet delta R', 40, 0d0, 8d0)
    call book(g_jet + 11, 'jet multiplicity pT25', 10, -0.5d0, 9.5d0)
    call book(g_jet + 12, 'jet multiplicity pT40', 10, -0.5d0, 9.5d0)
    call book(g_jet + 13, 'leading jet pT', 50, 0d0, 1000d0)
    call book(g_jet + 14, 'leading jet eta', 50, -5d0, 5d0)
    call book(g_jet + 15, 'subleading jet pT', 50, 0d0, 600d0)
    call book(g_jet + 16, 'subleading jet eta', 50, -5d0, 5d0)
    call book(g_jet + 17, 'additional jet multiplicity pT25', 8, -0.5d0, 7.5d0)
    call book(g_jet + 18, 'leading additional jet pT', 50, 0d0, 800d0)
    call book(g_jet + 19, 'leading additional jet rapidity', 50, -5d0, 5d0)
    call book(g_jet + 20, 'additional jet HT', 50, 0d0, 1000d0)
    call book(g_jet + 21, 'jet HT pT25', 50, 0d0, 1500d0)
    call book(g_jet + 22, 'visible HT leptons and jets', 50, 0d0, 2000d0)
    call book(g_jet + 23, 'visible ST including missing pT', 50, 0d0, 2500d0)
    call book(g_jet + 24, 'visible dilepton bb mass', 50, 0d0, 2000d0)
    call book(g_jet + 25, 'visible dilepton bb pT', 50, 0d0, 1000d0)
    call book(g_jet + 26, 'visible dilepton bb rapidity', 50, -5d0, 5d0)
    call book(g_jet + 27, 'additional jet gap Q0', 20, 10d0, 210d0)
    call book(g_jet + 28, 'additional jet gap Qsum', 20, 10d0, 210d0)
    call book(g_jet + 29, 'minimum lepton b jet delta R', 40, 0d0, 8d0)
    call book(g_jet + 30, 'maximum lepton b jet delta R', 40, 0d0, 8d0)

    call book_fiducial_histograms()
    call book_reconstructed_histograms()
    call book_spin_histograms(g_spin, 'truth spin')
    call book_spin_histograms(g_reco_spin, 'reco spin')

    call book(g_diagnostic + 1, 'diagnostic truth top assignment score', 50, 0d0, 100d0)
    call book(g_diagnostic + 2, 'diagnostic assigned radiation count', 6, -0.5d0, 5.5d0)
    call book(g_diagnostic + 3, 'diagnostic truth top mass residual', 50, 0d0, 50d0)
    call book(g_diagnostic + 4, 'diagnostic reco W mass residual', 50, 0d0, 50d0)
    call book(g_diagnostic + 5, 'diagnostic reconstruction rate', 1, 0.5d0, 1.5d0)

    if (g_diagnostic + 5 /= number_of_histograms) then
      call fail_analysis('internal histogram count is inconsistent')
    end if
  end subroutine book_histograms


  subroutine book_fiducial_histograms()
    implicit none
    call book(g_fiducial + 1, 'fiducial rate', 1, 0.5d0, 1.5d0)
    call book(g_fiducial + 2, 'fiducial leading lepton pT', 50, 0d0, 500d0)
    call book(g_fiducial + 3, 'fiducial trailing lepton pT', 50, 0d0, 300d0)
    call book(g_fiducial + 4, 'fiducial positive lepton eta', 40, -2.5d0, 2.5d0)
    call book(g_fiducial + 5, 'fiducial negative lepton eta', 40, -2.5d0, 2.5d0)
    call book(g_fiducial + 6, 'fiducial dilepton mass', 50, 0d0, 1000d0)
    call book(g_fiducial + 7, 'fiducial dilepton pT', 50, 0d0, 600d0)
    call book(g_fiducial + 8, 'fiducial dilepton rapidity', 40, -4d0, 4d0)
    call book(g_fiducial + 9, 'fiducial dilepton delta phi', 32, 0d0, pi)
    call book(g_fiducial + 10, 'fiducial abs dilepton delta eta', 25, 0d0, 5d0)
    call book(g_fiducial + 11, 'fiducial dilepton delta R', 40, 0d0, 8d0)
    call book(g_fiducial + 12, 'fiducial missing pT', 50, 0d0, 600d0)
    call book(g_fiducial + 13, 'fiducial dilepton missing transverse mass', 50, 0d0, 1000d0)
    call book(g_fiducial + 14, 'fiducial min delta phi missing lepton', 32, 0d0, pi)
    call book(g_fiducial + 15, 'fiducial leading b jet pT', 50, 0d0, 600d0)
    call book(g_fiducial + 16, 'fiducial trailing b jet pT', 50, 0d0, 400d0)
    call book(g_fiducial + 17, 'fiducial abs leading b jet eta', 25, 0d0, 2.5d0)
    call book(g_fiducial + 18, 'fiducial abs trailing b jet eta', 25, 0d0, 2.5d0)
    call book(g_fiducial + 19, 'fiducial bb jet mass', 50, 0d0, 1200d0)
    call book(g_fiducial + 20, 'fiducial bb jet pT', 50, 0d0, 800d0)
    call book(g_fiducial + 21, 'fiducial bb jet delta R', 40, 0d0, 8d0)
    call book(g_fiducial + 22, 'fiducial visible dilepton bb mass', 50, 0d0, 2000d0)
    call book(g_fiducial + 23, 'fiducial visible dilepton bb pT', 50, 0d0, 1000d0)
    call book(g_fiducial + 24, 'fiducial visible dilepton bb rapidity', 40, -4d0, 4d0)
    call book(g_fiducial + 25, 'fiducial jet HT', 50, 0d0, 1500d0)
    call book(g_fiducial + 26, 'fiducial visible HT', 50, 0d0, 2000d0)
    call book(g_fiducial + 27, 'fiducial ST', 50, 0d0, 2500d0)
    call book(g_fiducial + 28, 'fiducial jet multiplicity', 10, -0.5d0, 9.5d0)
    call book(g_fiducial + 29, 'fiducial additional jet multiplicity', 8, -0.5d0, 7.5d0)
    call book(g_fiducial + 30, 'fiducial minimum paired lb mass', 50, 0d0, 250d0)
    call book(g_fiducial + 31, 'fiducial maximum paired lb mass', 50, 0d0, 250d0)
    call book(g_fiducial + 32, 'fiducial average paired lb mass', 50, 0d0, 250d0)
    call book(g_fiducial + 33, 'fiducial mT2 dilepton', 50, 0d0, 250d0)
    call book(g_fiducial + 34, 'fiducial mT2 lb systems', 50, 0d0, 400d0)
    call book(g_fiducial + 35, 'fiducial delta abs lepton eta', 40, -4d0, 4d0)
    call book(g_fiducial + 36, 'fiducial lepton pT ratio', 40, 0d0, 1d0)
    call book(g_fiducial + 37, 'fiducial b jet pT ratio', 40, 0d0, 1d0)
    call book(g_fiducial + 38, 'fiducial lepton opening cosine lab', 40, -1d0, 1d0)
  end subroutine book_fiducial_histograms


  subroutine book_reconstructed_histograms()
    implicit none
    call book(g_reco + 1, 'reconstructed rate', 1, 0.5d0, 1.5d0)
    call book(g_reco + 2, 'reconstructed Wplus pT', 50, 0d0, 600d0)
    call book(g_reco + 3, 'reconstructed Wplus rapidity', 40, -4d0, 4d0)
    call book(g_reco + 4, 'reconstructed Wplus mass', 50, 0d0, 150d0)
    call book(g_reco + 5, 'reconstructed Wminus pT', 50, 0d0, 600d0)
    call book(g_reco + 6, 'reconstructed Wminus rapidity', 40, -4d0, 4d0)
    call book(g_reco + 7, 'reconstructed Wminus mass', 50, 0d0, 150d0)
    call book(g_reco + 8, 'reconstructed top pT', 50, 0d0, 1000d0)
    call book(g_reco + 9, 'reconstructed top rapidity', 40, -4d0, 4d0)
    call book(g_reco + 10, 'reconstructed top mass', 50, 100d0, 250d0)
    call book(g_reco + 11, 'reconstructed antitop pT', 50, 0d0, 1000d0)
    call book(g_reco + 12, 'reconstructed antitop rapidity', 40, -4d0, 4d0)
    call book(g_reco + 13, 'reconstructed antitop mass', 50, 100d0, 250d0)
    call book(g_reco + 14, 'reconstructed ttbar mass', 54, 300d0, 3000d0)
    call book(g_reco + 15, 'reconstructed ttbar pT', 50, 0d0, 1000d0)
    call book(g_reco + 16, 'reconstructed ttbar rapidity', 40, -4d0, 4d0)
    call book(g_reco + 17, 'reconstructed delta top rapidity', 40, -4d0, 4d0)
    call book(g_reco + 18, 'reconstructed delta abs top rapidity', 40, -4d0, 4d0)
    call book(g_reco + 19, 'reconstructed top scattering cosine', 40, -1d0, 1d0)
    call book(g_reco + 20, 'reconstructed top speed ttbar frame', 40, 0d0, 1d0)
    call book(g_reco + 21, 'reconstructed neutrino pz', 50, -1000d0, 1000d0)
    call book(g_reco + 22, 'reconstructed antineutrino pz', 50, -1000d0, 1000d0)
    call book(g_reco + 23, 'reconstructed dineutrino mass', 50, 0d0, 1500d0)
    call book(g_reco + 24, 'reconstructed dineutrino pz', 50, -1500d0, 1500d0)
    call book(g_reco + 25, 'reconstruction constraint score', 50, 0d0, 100d0)
    call book(g_reco + 26, 'reconstructed top mass difference', 50, -100d0, 100d0)
    call book(g_reco + 27, 'reconstructed b lepton pairing', 2, 0.5d0, 2.5d0)
    call book(g_reco + 28, 'reco truth average top pT ratio', 50, 0d0, 2.5d0)
    call book(g_reco + 29, 'reco truth ttbar mass residual', 50, -500d0, 500d0)
    call book(g_reco + 30, 'reconstructed Wplus helicity cosine', 40, -1d0, 1d0)
    call book(g_reco + 31, 'reconstructed Wminus helicity cosine', 40, -1d0, 1d0)
  end subroutine book_reconstructed_histograms


  subroutine book_spin_histograms(offset, prefix)
    implicit none
    integer, intent(in) :: offset
    character(len=*), intent(in) :: prefix

    call book(offset + 1, trim(prefix)//' cos plus k', 40, -1d0, 1d0)
    call book(offset + 2, trim(prefix)//' cos plus r', 40, -1d0, 1d0)
    call book(offset + 3, trim(prefix)//' cos plus n', 40, -1d0, 1d0)
    call book(offset + 4, trim(prefix)//' cos minus k', 40, -1d0, 1d0)
    call book(offset + 5, trim(prefix)//' cos minus r', 40, -1d0, 1d0)
    call book(offset + 6, trim(prefix)//' cos minus n', 40, -1d0, 1d0)
    call book(offset + 7, trim(prefix)//' product kk', 40, -1d0, 1d0)
    call book(offset + 8, trim(prefix)//' product kr', 40, -1d0, 1d0)
    call book(offset + 9, trim(prefix)//' product kn', 40, -1d0, 1d0)
    call book(offset + 10, trim(prefix)//' product rk', 40, -1d0, 1d0)
    call book(offset + 11, trim(prefix)//' product rr', 40, -1d0, 1d0)
    call book(offset + 12, trim(prefix)//' product rn', 40, -1d0, 1d0)
    call book(offset + 13, trim(prefix)//' product nk', 40, -1d0, 1d0)
    call book(offset + 14, trim(prefix)//' product nr', 40, -1d0, 1d0)
    call book(offset + 15, trim(prefix)//' product nn', 40, -1d0, 1d0)
    call book(offset + 16, trim(prefix)//' lepton opening cosine', 40, -1d0, 1d0)
    call book(offset + 17, trim(prefix)//' lab dilepton delta phi', 32, 0d0, pi)
    call book(offset + 18, trim(prefix)//' symmetric kr rk', 40, -2d0, 2d0)
    call book(offset + 19, trim(prefix)//' antisymmetric kr rk', 40, -2d0, 2d0)
    call book(offset + 20, trim(prefix)//' symmetric kn nk', 40, -2d0, 2d0)
    call book(offset + 21, trim(prefix)//' antisymmetric kn nk', 40, -2d0, 2d0)
    call book(offset + 22, trim(prefix)//' symmetric rn nr', 40, -2d0, 2d0)
    call book(offset + 23, trim(prefix)//' antisymmetric rn nr', 40, -2d0, 2d0)
    call book(offset + 24, trim(prefix)//' lepton triple product', 40, -1d0, 1d0)
    call book(offset + 25, trim(prefix)//' diagonal D estimator', 40, -1d0, 1d0)
    call book(offset + 26, trim(prefix)//' signed sine delta phi', 40, -1d0, 1d0)
  end subroutine book_spin_histograms


  subroutine book(label, title, bins, lower, upper)
    implicit none
    integer, intent(in) :: label, bins
    character(len=*), intent(in) :: title
    double precision, intent(in) :: lower, upper
    call HwU_book(label, trim(title)//' |T@NLO', bins, lower, upper)
  end subroutine book


  subroutine find_dilepton_objects(p, istatus, ipdg, &
                                    ilp, ilm, inu, ianu, ib, iab)
    implicit none
    double precision, intent(in) :: p(0:, :)
    integer, intent(in) :: istatus(:), ipdg(:)
    integer, intent(out) :: ilp, ilm, inu, ianu, ib, iab
    integer :: i, nlp, nlm, nnu, nanu, nb, nab

    ilp = 0
    ilm = 0
    inu = 0
    ianu = 0
    ib = 0
    iab = 0
    nlp = 0
    nlm = 0
    nb = 0
    nab = 0
    do i = nincoming + 1, size(p, 2)
      if (istatus(i) /= 1) cycle
      select case (abs(ipdg(i)))
      case (11, 13, 15)
        if (ipdg(i) < 0) then
          nlp = nlp + 1
          ilp = i
        else
          nlm = nlm + 1
          ilm = i
        end if
      case (5)
        if (ipdg(i) == 5) then
          nb = nb + 1
          ib = i
        else
          nab = nab + 1
          iab = i
        end if
      end select
    end do
    if (nlp /= 1 .or. nlm /= 1 .or. nb /= 1 .or. nab /= 1) then
      call fail_analysis('requires exactly lplus lminus b and bbar')
    end if

    nnu = 0
    nanu = 0
    do i = nincoming + 1, size(p, 2)
      if (istatus(i) /= 1) cycle
      if (ipdg(i) == abs(ipdg(ilp)) + 1) then
        nnu = nnu + 1
        inu = i
      end if
      if (ipdg(i) == -(abs(ipdg(ilm)) + 1)) then
        nanu = nanu + 1
        ianu = i
      end if
    end do
    if (nnu /= 1 .or. nanu /= 1) then
      call fail_analysis('charged leptons have no unique matching neutrinos')
    end if
  end subroutine find_dilepton_objects


  subroutine reconstruct_truth_tops(p, istatus, ipdg, wp, wm, b, ab, &
                                    top, atop, best_score, assigned_count)
    implicit none
    double precision, intent(in) :: p(0:, :), wp(0:3), wm(0:3)
    double precision, intent(in) :: b(0:3), ab(0:3)
    integer, intent(in) :: istatus(:), ipdg(:)
    double precision, intent(out) :: top(0:3), atop(0:3), best_score
    integer, intent(out) :: assigned_count
    integer :: extra_indices(size(p, 2)), nextra, i, code, work, digit
    integer :: combinations, count
    double precision :: candidate_top(0:3), candidate_atop(0:3), score

    nextra = 0
    do i = nincoming + 1, size(p, 2)
      if (istatus(i) /= 1 .or. p(0, i) <= tiny) cycle
      if (.not. (abs(ipdg(i)) <= 4 .or. ipdg(i) == 21)) cycle
      nextra = nextra + 1
      extra_indices(nextra) = i
    end do
    if (nextra > 10) then
      call fail_analysis('truth resonance assignment exceeds ten emissions')
    end if

    best_score = huge(1d0)
    top = wp + b
    atop = wm + ab
    assigned_count = 0
    combinations = 3**nextra
    do code = 0, combinations - 1
      candidate_top = wp + b
      candidate_atop = wm + ab
      work = code
      count = 0
      do i = 1, nextra
        digit = mod(work, 3)
        work = work/3
        if (digit == 1) then
          candidate_top = candidate_top + p(0:3, extra_indices(i))
          count = count + 1
        else if (digit == 2) then
          candidate_atop = candidate_atop + p(0:3, extra_indices(i))
          count = count + 1
        end if
      end do
      score = ((invariant_mass(candidate_top) - top_mass_reference)/5d0)**2 + &
              ((invariant_mass(candidate_atop) - top_mass_reference)/5d0)**2
      if (score < best_score) then
        best_score = score
        top = candidate_top
        atop = candidate_atop
        assigned_count = count
      end if
    end do
  end subroutine reconstruct_truth_tops


  subroutine cluster_qcd_jets(p, istatus, ipdg, pjet, njet, jet_of_parton, &
                              ibjet, iabjet)
    implicit none
    double precision, intent(in) :: p(0:, :)
    integer, intent(in) :: istatus(:), ipdg(:)
    double precision, intent(out) :: pjet(0:, :)
    integer, intent(out) :: njet, jet_of_parton(:)
    integer, intent(out) :: ibjet, iabjet
    double precision :: pqcd(0:3, size(p, 2))
    integer :: nqcd, i, qcd_b, qcd_ab

    pqcd = 0d0
    pjet = 0d0
    jet_of_parton = 0
    nqcd = 0
    qcd_b = 0
    qcd_ab = 0
    do i = nincoming + 1, size(p, 2)
      if (istatus(i) /= 1 .or. p(0, i) <= tiny) cycle
      if (.not. (abs(ipdg(i)) <= 5 .or. ipdg(i) == 21)) cycle
      nqcd = nqcd + 1
      pqcd(:, nqcd) = p(0:3, i)
      if (ipdg(i) == 5) qcd_b = nqcd
      if (ipdg(i) == -5) qcd_ab = nqcd
    end do
    if (nqcd < 2 .or. qcd_b == 0 .or. qcd_ab == 0) then
      call fail_analysis('cannot form the two bottom-flavoured jets')
    end if

    call amcatnlo_fastjetppgenkt_etamax(pqcd(:, 1:nqcd), nqcd, &
         jet_radius, 0d0, -1d0, -1d0, pjet(:, 1:nqcd), njet, &
         jet_of_parton(1:nqcd))
    ibjet = jet_of_parton(qcd_b)
    iabjet = jet_of_parton(qcd_ab)
  end subroutine cluster_qcd_jets


  subroutine classify_jets(pjet, njet, ibjet, iabjet, njet25, njet40, &
                           nextra25, lead, trail, lead_extra, jet_ht, extra_ht)
    implicit none
    double precision, intent(in) :: pjet(0:, :)
    integer, intent(in) :: njet, ibjet, iabjet
    integer, intent(out) :: njet25, njet40, nextra25
    integer, intent(out) :: lead, trail, lead_extra
    double precision, intent(out) :: jet_ht, extra_ht
    integer :: i
    double precision :: jet_pt, jet_eta

    njet25 = 0
    njet40 = 0
    nextra25 = 0
    lead = 0
    trail = 0
    lead_extra = 0
    jet_ht = 0d0
    extra_ht = 0d0
    do i = 1, njet
      jet_pt = transverse_momentum(pjet(:, i))
      jet_eta = abs(pseudorapidity(pjet(:, i)))
      if (jet_eta >= analysis_jet_eta_cut) cycle
      if (jet_pt > analysis_jet_pt_cut) then
        njet25 = njet25 + 1
        jet_ht = jet_ht + jet_pt
        if (lead == 0) then
          lead = i
        else if (trail == 0) then
          trail = i
        end if
        if (i /= ibjet .and. i /= iabjet) then
          nextra25 = nextra25 + 1
          extra_ht = extra_ht + jet_pt
          if (lead_extra == 0) lead_extra = i
        end if
      end if
      if (jet_pt > 40d0) njet40 = njet40 + 1
    end do
  end subroutine classify_jets


  subroutine missing_transverse_momentum(p, istatus, ipdg, met)
    implicit none
    double precision, intent(in) :: p(0:, :)
    integer, intent(in) :: istatus(:), ipdg(:)
    double precision, intent(out) :: met(0:1)
    integer :: i

    met = 0d0
    do i = nincoming + 1, size(p, 2)
      if (istatus(i) /= 1) cycle
      if (abs(ipdg(i)) == 12 .or. abs(ipdg(i)) == 14 .or. &
          abs(ipdg(i)) == 16) cycle
      met(0) = met(0) - p(1, i)
      met(1) = met(1) - p(2, i)
    end do
  end subroutine missing_transverse_momentum


  subroutine lepton_jet_distances(lp, lm, pjet, njet, min_all, &
                                  min_b, max_b, ibjet, iabjet)
    implicit none
    double precision, intent(in) :: lp(0:3), lm(0:3), pjet(0:, :)
    integer, intent(in) :: njet, ibjet, iabjet
    double precision, intent(out) :: min_all, min_b, max_b
    double precision :: distances(4), value
    integer :: i

    distances(1) = delta_r(lp, pjet(:, ibjet))
    distances(2) = delta_r(lp, pjet(:, iabjet))
    distances(3) = delta_r(lm, pjet(:, ibjet))
    distances(4) = delta_r(lm, pjet(:, iabjet))
    min_b = minval(distances)
    max_b = maxval(distances)
    min_all = huge(1d0)
    do i = 1, njet
      if (transverse_momentum(pjet(:, i)) <= analysis_jet_pt_cut) cycle
      if (abs(pseudorapidity(pjet(:, i))) >= analysis_jet_eta_cut) cycle
      value = min(delta_r(lp, pjet(:, i)), delta_r(lm, pjet(:, i)))
      min_all = min(min_all, value)
    end do
    if (min_all > 0.5d0*huge(1d0)) min_all = 99d0
  end subroutine lepton_jet_distances


  subroutine lepton_b_pair_masses(lp, lm, bjet, abjet, &
                                  minimum_mass, maximum_mass, average_mass)
    implicit none
    double precision, intent(in) :: lp(0:3), lm(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3)
    double precision, intent(out) :: minimum_mass, maximum_mass, average_mass
    double precision :: pair_a(2), pair_b(2), chosen(2)

    pair_a(1) = invariant_mass(lp + bjet)
    pair_a(2) = invariant_mass(lm + abjet)
    pair_b(1) = invariant_mass(lp + abjet)
    pair_b(2) = invariant_mass(lm + bjet)
    if (maxval(pair_a) <= maxval(pair_b)) then
      chosen = pair_a
    else
      chosen = pair_b
    end if
    minimum_mass = minval(chosen)
    maximum_mass = maxval(chosen)
    average_mass = 0.5d0*sum(chosen)
  end subroutine lepton_b_pair_masses


  subroutine reconstruct_neutrinos(lp, lm, bjet, abjet, met, &
                                   best_nu, best_anu, best_pairing, &
                                   best_score, success)
    implicit none
    double precision, intent(in) :: lp(0:3), lm(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), met(0:1)
    double precision, intent(out) :: best_nu(0:3), best_anu(0:3)
    integer, intent(out) :: best_pairing
    double precision, intent(out) :: best_score
    logical, intent(out) :: success
    double precision :: qx, qy, scale, step, trial_score
    double precision :: trial_nu(0:3), trial_anu(0:3)
    double precision :: candidate_qx, candidate_qy
    integer :: ix, iy, iteration, direction, trial_pairing
    double precision, parameter :: directions(2, 8) = reshape((/ &
         1d0, 0d0, -1d0, 0d0, 0d0, 1d0, 0d0, -1d0, &
         1d0, 1d0, 1d0, -1d0, -1d0, 1d0, -1d0, -1d0 /), (/2, 8/))

    scale = max(50d0, vector_magnitude_2(met), &
                transverse_momentum(lp), transverse_momentum(lm))
    best_score = huge(1d0)
    qx = 0.5d0*met(0)
    qy = 0.5d0*met(1)
    best_nu = 0d0
    best_anu = 0d0
    best_pairing = 0

    do ix = -2, 2
      do iy = -2, 2
        candidate_qx = 0.5d0*met(0) + 0.5d0*scale*dble(ix)
        candidate_qy = 0.5d0*met(1) + 0.5d0*scale*dble(iy)
        call neutrino_partition_score(lp, lm, bjet, abjet, met, &
             candidate_qx, candidate_qy, trial_score, trial_nu, trial_anu, &
             trial_pairing)
        if (trial_score < best_score) then
          best_score = trial_score
          qx = candidate_qx
          qy = candidate_qy
          best_nu = trial_nu
          best_anu = trial_anu
          best_pairing = trial_pairing
        end if
      end do
    end do

    step = 0.5d0*scale
    do iteration = 1, 12
      success = .false.
      do direction = 1, 8
        candidate_qx = qx + step*directions(1, direction)
        candidate_qy = qy + step*directions(2, direction)
        call neutrino_partition_score(lp, lm, bjet, abjet, met, &
             candidate_qx, candidate_qy, trial_score, trial_nu, trial_anu, &
             trial_pairing)
        if (trial_score < best_score) then
          best_score = trial_score
          qx = candidate_qx
          qy = candidate_qy
          best_nu = trial_nu
          best_anu = trial_anu
          best_pairing = trial_pairing
          success = .true.
        end if
      end do
      if (.not. success) step = 0.5d0*step
      if (step < 0.25d0) exit
    end do
    success = best_pairing > 0 .and. best_score < huge(1d0)/2d0
  end subroutine reconstruct_neutrinos


  subroutine neutrino_partition_score(lp, lm, bjet, abjet, met, qx, qy, &
                                      best_score, best_nu, best_anu, pairing)
    implicit none
    double precision, intent(in) :: lp(0:3), lm(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), met(0:1)
    double precision, intent(in) :: qx, qy
    double precision, intent(out) :: best_score, best_nu(0:3), best_anu(0:3)
    integer, intent(out) :: pairing
    double precision :: pz_plus(2), pz_minus(2), trial_nu(0:3), trial_anu(0:3)
    double precision :: wp(0:3), wm(0:3), top(0:3), atop(0:3), score
    integer :: roots_plus, roots_minus, i, j, pair

    call w_constraint_pz(lp, qx, qy, pz_plus, roots_plus)
    call w_constraint_pz(lm, met(0) - qx, met(1) - qy, &
         pz_minus, roots_minus)
    best_score = huge(1d0)
    best_nu = 0d0
    best_anu = 0d0
    pairing = 0
    do i = 1, roots_plus
      call massless_four_vector(qx, qy, pz_plus(i), trial_nu)
      wp = lp + trial_nu
      do j = 1, roots_minus
        call massless_four_vector(met(0) - qx, met(1) - qy, &
             pz_minus(j), trial_anu)
        wm = lm + trial_anu
        do pair = 1, 2
          if (pair == 1) then
            top = wp + bjet
            atop = wm + abjet
          else
            top = wp + abjet
            atop = wm + bjet
          end if
          score = ((invariant_mass(top) - top_mass_reference)/15d0)**2 + &
                  ((invariant_mass(atop) - top_mass_reference)/15d0)**2 + &
                  ((invariant_mass(wp) - w_mass_reference)/10d0)**2 + &
                  ((invariant_mass(wm) - w_mass_reference)/10d0)**2 + &
                  ((invariant_mass(top) - invariant_mass(atop))/30d0)**2
          if (score < best_score) then
            best_score = score
            best_nu = trial_nu
            best_anu = trial_anu
            pairing = pair
          end if
        end do
      end do
    end do
  end subroutine neutrino_partition_score


  subroutine w_constraint_pz(lepton, qx, qy, roots, number_of_roots)
    implicit none
    double precision, intent(in) :: lepton(0:3), qx, qy
    double precision, intent(out) :: roots(2)
    integer, intent(out) :: number_of_roots
    double precision :: denominator, mu, discriminant, root_disc

    denominator = lepton(0)**2 - lepton(3)**2
    if (denominator <= tiny) then
      roots = 0d0
      number_of_roots = 1
      return
    end if
    mu = 0.5d0*(w_mass_reference**2 - invariant_mass2(lepton)) + &
         lepton(1)*qx + lepton(2)*qy
    discriminant = mu**2 - denominator*(qx**2 + qy**2)
    if (discriminant >= 0d0) then
      root_disc = sqrt(discriminant)
      roots(1) = (mu*lepton(3) + lepton(0)*root_disc)/denominator
      roots(2) = (mu*lepton(3) - lepton(0)*root_disc)/denominator
      number_of_roots = 2
    else
      roots(1) = mu*lepton(3)/denominator
      roots(2) = roots(1)
      number_of_roots = 1
    end if
  end subroutine w_constraint_pz


  subroutine compute_mt2(first, second, met, result)
    implicit none
    double precision, intent(in) :: first(0:3), second(0:3), met(0:1)
    double precision, intent(out) :: result
    double precision :: qx, qy, scale, step, candidate, best
    double precision :: candidate_qx, candidate_qy
    integer :: ix, iy, iteration, direction
    logical :: improved
    double precision, parameter :: directions(2, 8) = reshape((/ &
         1d0, 0d0, -1d0, 0d0, 0d0, 1d0, 0d0, -1d0, &
         1d0, 1d0, 1d0, -1d0, -1d0, 1d0, -1d0, -1d0 /), (/2, 8/))

    scale = max(50d0, vector_magnitude_2(met), &
                transverse_momentum(first), transverse_momentum(second))
    best = huge(1d0)
    qx = 0.5d0*met(0)
    qy = 0.5d0*met(1)
    do ix = -2, 2
      do iy = -2, 2
        candidate_qx = 0.5d0*met(0) + 0.5d0*scale*dble(ix)
        candidate_qy = 0.5d0*met(1) + 0.5d0*scale*dble(iy)
        candidate = mt2_objective(first, second, met, &
                                  candidate_qx, candidate_qy)
        if (candidate < best) then
          best = candidate
          qx = candidate_qx
          qy = candidate_qy
        end if
      end do
    end do
    step = 0.5d0*scale
    do iteration = 1, 12
      improved = .false.
      do direction = 1, 8
        candidate_qx = qx + step*directions(1, direction)
        candidate_qy = qy + step*directions(2, direction)
        candidate = mt2_objective(first, second, met, &
                                  candidate_qx, candidate_qy)
        if (candidate < best) then
          best = candidate
          qx = candidate_qx
          qy = candidate_qy
          improved = .true.
        end if
      end do
      if (.not. improved) step = 0.5d0*step
      if (step < 0.25d0) exit
    end do
    result = best
  end subroutine compute_mt2


  subroutine compute_mt2_blbl(lp, lm, bjet, abjet, met, result)
    implicit none
    double precision, intent(in) :: lp(0:3), lm(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), met(0:1)
    double precision, intent(out) :: result
    double precision :: first_pairing, second_pairing
    call compute_mt2(lp + bjet, lm + abjet, met, first_pairing)
    call compute_mt2(lp + abjet, lm + bjet, met, second_pairing)
    result = min(first_pairing, second_pairing)
  end subroutine compute_mt2_blbl


  double precision function mt2_objective(first, second, met, qx, qy)
    implicit none
    double precision, intent(in) :: first(0:3), second(0:3), met(0:1)
    double precision, intent(in) :: qx, qy
    mt2_objective = max(transverse_mass(first, qx, qy), &
         transverse_mass(second, met(0) - qx, met(1) - qy))
  end function mt2_objective


  subroutine compute_spin_observables(top, atop, lp, lm, values, products, extra)
    implicit none
    double precision, intent(in) :: top(0:3), atop(0:3)
    double precision, intent(in) :: lp(0:3), lm(0:3)
    double precision, intent(out) :: values(6), products(9), extra(11)
    double precision :: ttbar(0:3), top_tt(0:3), atop_tt(0:3)
    double precision :: lp_tt(0:3), lm_tt(0:3), lp_rest(0:3), lm_rest(0:3)
    double precision :: beam(0:3), beam_tt(0:3)
    double precision :: k_axis(3), r_axis(3), n_axis(3)
    double precision :: lp_direction(3), lm_direction(3), cross_leptons(3)
    double precision :: normalized_axis(3)
    double precision :: cosine_scattering, transverse_norm, sign_factor

    values = 0d0
    products = 0d0
    extra = 0d0
    ttbar = top + atop
    call boost_to_rest(top, ttbar, top_tt)
    call boost_to_rest(atop, ttbar, atop_tt)
    call boost_to_rest(lp, ttbar, lp_tt)
    call boost_to_rest(lm, ttbar, lm_tt)
    call boost_to_rest(lp_tt, top_tt, lp_rest)
    call boost_to_rest(lm_tt, atop_tt, lm_rest)
    beam = (/1d0, 0d0, 0d0, 1d0/)
    call boost_to_rest(beam, ttbar, beam_tt)

    call unit_vector(top_tt(1:3), k_axis)
    call unit_vector(beam_tt(1:3), r_axis)
    cosine_scattering = dot_product(k_axis, r_axis)
    transverse_norm = sqrt(max(0d0, 1d0 - cosine_scattering**2))
    if (transverse_norm <= 1d-10) then
      r_axis = (/1d0, 0d0, 0d0/)
      if (abs(dot_product(r_axis, k_axis)) > 0.9d0) &
           r_axis = (/0d0, 1d0, 0d0/)
      r_axis = r_axis - dot_product(r_axis, k_axis)*k_axis
      call unit_vector(r_axis, normalized_axis)
      r_axis = normalized_axis
      call cross_product3(r_axis, k_axis, n_axis)
    else
      sign_factor = sign(1d0, cosine_scattering)
      r_axis = sign_factor*(r_axis - cosine_scattering*k_axis)/transverse_norm
      call cross_product3(beam_tt(1:3), top_tt(1:3), n_axis)
      call unit_vector(n_axis, normalized_axis)
      n_axis = normalized_axis
      n_axis = sign_factor*n_axis
    end if
    call unit_vector(lp_rest(1:3), lp_direction)
    call unit_vector(lm_rest(1:3), lm_direction)

    values(1) = clamp_cosine(dot_product(lp_direction, k_axis))
    values(2) = clamp_cosine(dot_product(lp_direction, r_axis))
    values(3) = clamp_cosine(dot_product(lp_direction, n_axis))
    values(4) = clamp_cosine(dot_product(lm_direction, k_axis))
    values(5) = clamp_cosine(dot_product(lm_direction, r_axis))
    values(6) = clamp_cosine(dot_product(lm_direction, n_axis))
    products = (/values(1)*values(4), values(1)*values(5), &
                 values(1)*values(6), values(2)*values(4), &
                 values(2)*values(5), values(2)*values(6), &
                 values(3)*values(4), values(3)*values(5), &
                 values(3)*values(6)/)
    extra(1) = clamp_cosine(dot_product(lp_direction, lm_direction))
    extra(2) = delta_phi(lp, lm)
    extra(3) = products(2) + products(4)
    extra(4) = products(2) - products(4)
    extra(5) = products(3) + products(7)
    extra(6) = products(3) - products(7)
    extra(7) = products(6) + products(8)
    extra(8) = products(6) - products(8)
    call cross_product3(lp_direction, lm_direction, cross_leptons)
    extra(9) = clamp_cosine(dot_product(cross_leptons, k_axis))
    extra(10) = (products(1) + products(5) + products(9))/3d0
    extra(11) = sin(signed_delta_phi(lp, lm))
  end subroutine compute_spin_observables


  subroutine w_helicity_cosines(top, atop, wp, wm, lp, lm, cos_plus, cos_minus)
    implicit none
    double precision, intent(in) :: top(0:3), atop(0:3), wp(0:3), wm(0:3)
    double precision, intent(in) :: lp(0:3), lm(0:3)
    double precision, intent(out) :: cos_plus, cos_minus
    double precision :: lp_w(0:3), lm_w(0:3), recoil_w(0:3)
    double precision :: direction_lepton(3), direction_recoil(3)

    call boost_to_rest(lp, wp, lp_w)
    call boost_to_rest(top - wp, wp, recoil_w)
    call unit_vector(lp_w(1:3), direction_lepton)
    call unit_vector(-recoil_w(1:3), direction_recoil)
    cos_plus = clamp_cosine(dot_product(direction_lepton, direction_recoil))

    call boost_to_rest(lm, wm, lm_w)
    call boost_to_rest(atop - wm, wm, recoil_w)
    call unit_vector(lm_w(1:3), direction_lepton)
    call unit_vector(-recoil_w(1:3), direction_recoil)
    cos_minus = clamp_cosine(dot_product(direction_lepton, direction_recoil))
  end subroutine w_helicity_cosines


  subroutine fill_inclusive_truth(wgts, top, atop, ttbar, lp, lm, ll, &
       wp, wm, ww, nu, anu, bjet, abjet, bbjet, met, pjet, &
       njet25, njet40, nextra25, lead_jet, trail_jet, leading_extra_pt, &
       leading_extra_y, extra_ht, jet_ht, visible_ht, st, visible, &
       min_dr_lepton_b, max_dr_lepton_b, spin_values, spin_products, &
       spin_extra, wp_helicity, wm_helicity)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: top(0:3), atop(0:3), ttbar(0:3)
    double precision, intent(in) :: lp(0:3), lm(0:3), ll(0:3)
    double precision, intent(in) :: wp(0:3), wm(0:3), ww(0:3)
    double precision, intent(in) :: nu(0:3), anu(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), bbjet(0:3)
    double precision, intent(in) :: met(0:1), pjet(0:, :)
    integer, intent(in) :: njet25, njet40, nextra25
    integer, intent(in) :: lead_jet, trail_jet
    double precision, intent(in) :: leading_extra_pt, leading_extra_y
    double precision, intent(in) :: extra_ht, jet_ht, visible_ht, st
    double precision, intent(in) :: visible(0:3), min_dr_lepton_b, max_dr_lepton_b
    double precision, intent(in) :: spin_values(6), spin_products(9)
    double precision, intent(in) :: spin_extra(11), wp_helicity, wm_helicity
    double precision :: ytop, yatop, yttbar, pt_top, pt_atop, mtt
    double precision :: top_tt(0:3), lp_top(0:3), lm_atop(0:3)
    double precision :: recoil_top(0:3), recoil_atop(0:3)
    double precision :: wp_top(0:3), wm_atop(0:3), nunu(0:3)
    double precision :: mlb_plus, mlb_minus, met_pt, dphi_met_lp, dphi_met_lm
    double precision :: denominator, scattering_cosine

    pt_top = transverse_momentum(top)
    pt_atop = transverse_momentum(atop)
    ytop = rapidity(top)
    yatop = rapidity(atop)
    yttbar = rapidity(ttbar)
    mtt = invariant_mass(ttbar)
    call top_scattering_variables(top, atop, scattering_cosine, denominator)
    call boost_to_rest(top, ttbar, top_tt)

    call HwU_fill(g_top + 1, pt_top, wgts)
    call HwU_fill(g_top + 2, pt_atop, wgts)
    call HwU_fill(g_top + 3, max(pt_top, pt_atop), wgts)
    call HwU_fill(g_top + 4, min(pt_top, pt_atop), wgts)
    call HwU_fill(g_top + 5, ytop, wgts)
    call HwU_fill(g_top + 6, yatop, wgts)
    call HwU_fill(g_top + 7, mtt, wgts)
    call HwU_fill(g_top + 8, transverse_momentum(ttbar), wgts)
    call HwU_fill(g_top + 9, yttbar, wgts)
    call HwU_fill(g_top + 10, abs(yttbar), wgts)
    call HwU_fill(g_top + 11, ytop - yatop, wgts)
    call HwU_fill(g_top + 12, abs(ytop - yatop), wgts)
    call HwU_fill(g_top + 13, abs(ytop) - abs(yatop), wgts)
    call HwU_fill(g_top + 14, 0.5d0*(ytop + yatop), wgts)
    call HwU_fill(g_top + 15, 0.5d0*(ytop - yatop), wgts)
    call HwU_fill(g_top + 16, scattering_cosine, wgts)
    call HwU_fill(g_top + 17, denominator, wgts)
    if (mtt > tiny) call HwU_fill(g_top + 18, 2d0*top_mass_reference/mtt, wgts)
    call HwU_fill(g_top + 19, delta_phi(top, atop), wgts)
    call HwU_fill(g_top + 20, delta_r(top, atop), wgts)
    call HwU_fill(g_top + 21, pt_top + pt_atop, wgts)
    call HwU_fill(g_top + 22, safe_ratio(pt_top - pt_atop, pt_top + pt_atop), wgts)
    call HwU_fill(g_top + 23, safe_ratio(2d0*(pt_top + pt_atop), mtt), wgts)
    call HwU_fill(g_top + 24, top_tt(0), wgts)
    call HwU_fill(g_top + 25, min(exp(abs(ytop - yatop)), 49.999d0), wgts)
    call HwU_fill(g_top + 26, &
         safe_ratio(abs(ttbar(3)), ttbar(0)), wgts)

    call HwU_fill(g_lep + 1, transverse_momentum(lp), wgts)
    call HwU_fill(g_lep + 2, transverse_momentum(lm), wgts)
    call HwU_fill(g_lep + 3, max(transverse_momentum(lp), transverse_momentum(lm)), wgts)
    call HwU_fill(g_lep + 4, min(transverse_momentum(lp), transverse_momentum(lm)), wgts)
    call HwU_fill(g_lep + 5, pseudorapidity(lp), wgts)
    call HwU_fill(g_lep + 6, pseudorapidity(lm), wgts)
    call HwU_fill(g_lep + 7, rapidity(lp), wgts)
    call HwU_fill(g_lep + 8, rapidity(lm), wgts)
    call HwU_fill(g_lep + 9, invariant_mass(ll), wgts)
    call HwU_fill(g_lep + 10, transverse_momentum(ll), wgts)
    call HwU_fill(g_lep + 11, rapidity(ll), wgts)
    call HwU_fill(g_lep + 12, abs(rapidity(ll)), wgts)
    call HwU_fill(g_lep + 13, delta_phi(lp, lm), wgts)
    call HwU_fill(g_lep + 14, pseudorapidity(lp) - pseudorapidity(lm), wgts)
    call HwU_fill(g_lep + 15, abs(pseudorapidity(lp) - pseudorapidity(lm)), wgts)
    call HwU_fill(g_lep + 16, delta_r(lp, lm), wgts)
    call HwU_fill(g_lep + 17, spatial_opening_cosine(lp, lm), wgts)
    call HwU_fill(g_lep + 18, &
         tanh(0.5d0*(pseudorapidity(lp) - pseudorapidity(lm))), wgts)
    call HwU_fill(g_lep + 19, safe_ratio(min(transverse_momentum(lp), &
         transverse_momentum(lm)), max(transverse_momentum(lp), &
         transverse_momentum(lm))), wgts)
    call HwU_fill(g_lep + 20, transverse_momentum(lp) + transverse_momentum(lm), wgts)
    call HwU_fill(g_lep + 21, abs(pseudorapidity(lp)) - abs(pseudorapidity(lm)), wgts)
    call HwU_fill(g_lep + 22, pi - delta_phi(lp, lm), wgts)
    call HwU_fill(g_lep + 23, signed_delta_phi(lp, lm), wgts)
    call HwU_fill(g_lep + 24, safe_ratio(invariant_mass(ll), &
         transverse_momentum(lp) + transverse_momentum(lm)), wgts)

    mlb_plus = invariant_mass(lp + bjet)
    mlb_minus = invariant_mass(lm + abjet)
    call boost_to_rest(lp, top, lp_top)
    call boost_to_rest(lm, atop, lm_atop)
    call boost_to_rest(top - wp, top, recoil_top)
    call boost_to_rest(atop - wm, atop, recoil_atop)
    call boost_to_rest(wp, top, wp_top)
    call boost_to_rest(wm, atop, wm_atop)
    call HwU_fill(g_decay + 1, transverse_momentum(wp), wgts)
    call HwU_fill(g_decay + 2, rapidity(wp), wgts)
    call HwU_fill(g_decay + 3, invariant_mass(wp), wgts)
    call HwU_fill(g_decay + 4, transverse_momentum(wm), wgts)
    call HwU_fill(g_decay + 5, rapidity(wm), wgts)
    call HwU_fill(g_decay + 6, invariant_mass(wm), wgts)
    call HwU_fill(g_decay + 7, invariant_mass(ww), wgts)
    call HwU_fill(g_decay + 8, transverse_momentum(ww), wgts)
    call HwU_fill(g_decay + 9, rapidity(ww), wgts)
    call HwU_fill(g_decay + 10, delta_phi(wp, wm), wgts)
    call HwU_fill(g_decay + 11, delta_r(wp, wm), wgts)
    call HwU_fill(g_decay + 12, mlb_plus, wgts)
    call HwU_fill(g_decay + 13, mlb_minus, wgts)
    call HwU_fill(g_decay + 14, min(mlb_plus, mlb_minus), wgts)
    call HwU_fill(g_decay + 15, max(mlb_plus, mlb_minus), wgts)
    call HwU_fill(g_decay + 16, 0.5d0*(mlb_plus + mlb_minus), wgts)
    call HwU_fill(g_decay + 17, lp_top(0), wgts)
    call HwU_fill(g_decay + 18, lm_atop(0), wgts)
    call HwU_fill(g_decay + 19, recoil_top(0), wgts)
    call HwU_fill(g_decay + 20, recoil_atop(0), wgts)
    call HwU_fill(g_decay + 21, wp_helicity, wgts)
    call HwU_fill(g_decay + 22, wm_helicity, wgts)
    call HwU_fill(g_decay + 23, wp_top(0), wgts)
    call HwU_fill(g_decay + 24, wm_atop(0), wgts)

    nunu = nu + anu
    met_pt = vector_magnitude_2(met)
    dphi_met_lp = delta_phi_vector(met, lp)
    dphi_met_lm = delta_phi_vector(met, lm)
    call HwU_fill(g_neutrino + 1, transverse_momentum(nu), wgts)
    call HwU_fill(g_neutrino + 2, pseudorapidity(nu), wgts)
    call HwU_fill(g_neutrino + 3, rapidity(nu), wgts)
    call HwU_fill(g_neutrino + 4, transverse_momentum(anu), wgts)
    call HwU_fill(g_neutrino + 5, pseudorapidity(anu), wgts)
    call HwU_fill(g_neutrino + 6, rapidity(anu), wgts)
    call HwU_fill(g_neutrino + 7, transverse_momentum(nunu), wgts)
    call HwU_fill(g_neutrino + 8, invariant_mass(nunu), wgts)
    call HwU_fill(g_neutrino + 9, rapidity(nunu), wgts)
    call HwU_fill(g_neutrino + 10, delta_phi(nu, anu), wgts)
    call HwU_fill(g_neutrino + 11, delta_r(nu, anu), wgts)
    call HwU_fill(g_neutrino + 12, met_pt, wgts)
    call HwU_fill(g_neutrino + 13, min(dphi_met_lp, dphi_met_lm), wgts)
    call HwU_fill(g_neutrino + 14, max(dphi_met_lp, dphi_met_lm), wgts)
    call HwU_fill(g_neutrino + 15, cluster_transverse_mass(ll, met), wgts)

    call fill_jet_truth(wgts, bjet, abjet, bbjet, pjet, njet25, &
         njet40, nextra25, lead_jet, trail_jet, leading_extra_pt, &
         leading_extra_y, extra_ht, jet_ht, visible_ht, st, visible, &
         min_dr_lepton_b, max_dr_lepton_b)
    call fill_spin(g_spin, spin_values, spin_products, spin_extra, wgts)
  end subroutine fill_inclusive_truth


  subroutine fill_jet_truth(wgts, bjet, abjet, bbjet, pjet, njet25, &
       njet40, nextra25, lead_jet, trail_jet, leading_extra_pt, &
       leading_extra_y, extra_ht, jet_ht, visible_ht, st, visible, &
       min_dr_lepton_b, max_dr_lepton_b)
    implicit none
    double precision, intent(in) :: wgts(*), bjet(0:3), abjet(0:3), bbjet(0:3)
    double precision, intent(in) :: pjet(0:, :), leading_extra_pt
    double precision, intent(in) :: leading_extra_y, extra_ht, jet_ht
    double precision, intent(in) :: visible_ht, st, visible(0:3)
    double precision, intent(in) :: min_dr_lepton_b, max_dr_lepton_b
    integer, intent(in) :: njet25, njet40, nextra25, lead_jet, trail_jet

    call HwU_fill(g_jet + 1, transverse_momentum(bjet), wgts)
    call HwU_fill(g_jet + 2, transverse_momentum(abjet), wgts)
    call HwU_fill(g_jet + 3, max(transverse_momentum(bjet), &
         transverse_momentum(abjet)), wgts)
    call HwU_fill(g_jet + 4, min(transverse_momentum(bjet), &
         transverse_momentum(abjet)), wgts)
    call HwU_fill(g_jet + 5, pseudorapidity(bjet), wgts)
    call HwU_fill(g_jet + 6, pseudorapidity(abjet), wgts)
    call HwU_fill(g_jet + 7, invariant_mass(bbjet), wgts)
    call HwU_fill(g_jet + 8, transverse_momentum(bbjet), wgts)
    call HwU_fill(g_jet + 9, delta_phi(bjet, abjet), wgts)
    call HwU_fill(g_jet + 10, delta_r(bjet, abjet), wgts)
    call HwU_fill(g_jet + 11, dble(njet25), wgts)
    call HwU_fill(g_jet + 12, dble(njet40), wgts)
    if (lead_jet > 0) then
      call HwU_fill(g_jet + 13, transverse_momentum(pjet(:, lead_jet)), wgts)
      call HwU_fill(g_jet + 14, pseudorapidity(pjet(:, lead_jet)), wgts)
    end if
    if (trail_jet > 0) then
      call HwU_fill(g_jet + 15, transverse_momentum(pjet(:, trail_jet)), wgts)
      call HwU_fill(g_jet + 16, pseudorapidity(pjet(:, trail_jet)), wgts)
    end if
    call HwU_fill(g_jet + 17, dble(nextra25), wgts)
    if (nextra25 > 0) then
      call HwU_fill(g_jet + 18, leading_extra_pt, wgts)
      call HwU_fill(g_jet + 19, leading_extra_y, wgts)
    end if
    call HwU_fill(g_jet + 20, extra_ht, wgts)
    call HwU_fill(g_jet + 21, jet_ht, wgts)
    call HwU_fill(g_jet + 22, visible_ht, wgts)
    call HwU_fill(g_jet + 23, st, wgts)
    call HwU_fill(g_jet + 24, invariant_mass(visible), wgts)
    call HwU_fill(g_jet + 25, transverse_momentum(visible), wgts)
    call HwU_fill(g_jet + 26, rapidity(visible), wgts)
    call HwU_fill(g_jet + 29, min_dr_lepton_b, wgts)
    call HwU_fill(g_jet + 30, max_dr_lepton_b, wgts)
  end subroutine fill_jet_truth


  subroutine fill_fiducial(wgts, lp, lm, ll, bjet, abjet, bbjet, met, &
       visible, jet_ht, visible_ht, st, njet25, nextra25, mlb_min, &
       mlb_max, mlb_average, mt2_ll, mt2_blbl)
    implicit none
    double precision, intent(in) :: wgts(*), lp(0:3), lm(0:3), ll(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), bbjet(0:3)
    double precision, intent(in) :: met(0:1), visible(0:3)
    double precision, intent(in) :: jet_ht, visible_ht, st
    double precision, intent(in) :: mlb_min, mlb_max, mlb_average
    double precision, intent(in) :: mt2_ll, mt2_blbl
    integer, intent(in) :: njet25, nextra25

    call HwU_fill(g_fiducial + 1, 1d0, wgts)
    call HwU_fill(g_fiducial + 2, max(transverse_momentum(lp), &
         transverse_momentum(lm)), wgts)
    call HwU_fill(g_fiducial + 3, min(transverse_momentum(lp), &
         transverse_momentum(lm)), wgts)
    call HwU_fill(g_fiducial + 4, pseudorapidity(lp), wgts)
    call HwU_fill(g_fiducial + 5, pseudorapidity(lm), wgts)
    call HwU_fill(g_fiducial + 6, invariant_mass(ll), wgts)
    call HwU_fill(g_fiducial + 7, transverse_momentum(ll), wgts)
    call HwU_fill(g_fiducial + 8, rapidity(ll), wgts)
    call HwU_fill(g_fiducial + 9, delta_phi(lp, lm), wgts)
    call HwU_fill(g_fiducial + 10, &
         abs(pseudorapidity(lp) - pseudorapidity(lm)), wgts)
    call HwU_fill(g_fiducial + 11, delta_r(lp, lm), wgts)
    call HwU_fill(g_fiducial + 12, vector_magnitude_2(met), wgts)
    call HwU_fill(g_fiducial + 13, cluster_transverse_mass(ll, met), wgts)
    call HwU_fill(g_fiducial + 14, min(delta_phi_vector(met, lp), &
         delta_phi_vector(met, lm)), wgts)
    call HwU_fill(g_fiducial + 15, max(transverse_momentum(bjet), &
         transverse_momentum(abjet)), wgts)
    call HwU_fill(g_fiducial + 16, min(transverse_momentum(bjet), &
         transverse_momentum(abjet)), wgts)
    if (transverse_momentum(bjet) >= transverse_momentum(abjet)) then
      call HwU_fill(g_fiducial + 17, abs(pseudorapidity(bjet)), wgts)
      call HwU_fill(g_fiducial + 18, abs(pseudorapidity(abjet)), wgts)
    else
      call HwU_fill(g_fiducial + 17, abs(pseudorapidity(abjet)), wgts)
      call HwU_fill(g_fiducial + 18, abs(pseudorapidity(bjet)), wgts)
    end if
    call HwU_fill(g_fiducial + 19, invariant_mass(bbjet), wgts)
    call HwU_fill(g_fiducial + 20, transverse_momentum(bbjet), wgts)
    call HwU_fill(g_fiducial + 21, delta_r(bjet, abjet), wgts)
    call HwU_fill(g_fiducial + 22, invariant_mass(visible), wgts)
    call HwU_fill(g_fiducial + 23, transverse_momentum(visible), wgts)
    call HwU_fill(g_fiducial + 24, rapidity(visible), wgts)
    call HwU_fill(g_fiducial + 25, jet_ht, wgts)
    call HwU_fill(g_fiducial + 26, visible_ht, wgts)
    call HwU_fill(g_fiducial + 27, st, wgts)
    call HwU_fill(g_fiducial + 28, dble(njet25), wgts)
    call HwU_fill(g_fiducial + 29, dble(nextra25), wgts)
    call HwU_fill(g_fiducial + 30, mlb_min, wgts)
    call HwU_fill(g_fiducial + 31, mlb_max, wgts)
    call HwU_fill(g_fiducial + 32, mlb_average, wgts)
    call HwU_fill(g_fiducial + 33, mt2_ll, wgts)
    call HwU_fill(g_fiducial + 34, mt2_blbl, wgts)
    call HwU_fill(g_fiducial + 35, &
         abs(pseudorapidity(lp)) - abs(pseudorapidity(lm)), wgts)
    call HwU_fill(g_fiducial + 36, safe_ratio(min(transverse_momentum(lp), &
         transverse_momentum(lm)), max(transverse_momentum(lp), &
         transverse_momentum(lm))), wgts)
    call HwU_fill(g_fiducial + 37, safe_ratio(min(transverse_momentum(bjet), &
         transverse_momentum(abjet)), max(transverse_momentum(bjet), &
         transverse_momentum(abjet))), wgts)
    call HwU_fill(g_fiducial + 38, spatial_opening_cosine(lp, lm), wgts)
  end subroutine fill_fiducial


  subroutine fill_reconstructed(wgts, nu, anu, wp, wm, top, atop, ttbar, &
       truth_top, truth_atop, truth_ttbar, score, pairing, wp_helicity, &
       wm_helicity, spin_values, spin_products, spin_extra)
    implicit none
    double precision, intent(in) :: wgts(*), nu(0:3), anu(0:3)
    double precision, intent(in) :: wp(0:3), wm(0:3), top(0:3), atop(0:3)
    double precision, intent(in) :: ttbar(0:3), truth_top(0:3)
    double precision, intent(in) :: truth_atop(0:3), truth_ttbar(0:3)
    double precision, intent(in) :: score, wp_helicity, wm_helicity
    integer, intent(in) :: pairing
    double precision, intent(in) :: spin_values(6), spin_products(9)
    double precision, intent(in) :: spin_extra(11)
    double precision :: scattering_cosine, speed, nunu(0:3)

    call top_scattering_variables(top, atop, scattering_cosine, speed)
    nunu = nu + anu
    call HwU_fill(g_reco + 1, 1d0, wgts)
    call HwU_fill(g_reco + 2, transverse_momentum(wp), wgts)
    call HwU_fill(g_reco + 3, rapidity(wp), wgts)
    call HwU_fill(g_reco + 4, invariant_mass(wp), wgts)
    call HwU_fill(g_reco + 5, transverse_momentum(wm), wgts)
    call HwU_fill(g_reco + 6, rapidity(wm), wgts)
    call HwU_fill(g_reco + 7, invariant_mass(wm), wgts)
    call HwU_fill(g_reco + 8, transverse_momentum(top), wgts)
    call HwU_fill(g_reco + 9, rapidity(top), wgts)
    call HwU_fill(g_reco + 10, invariant_mass(top), wgts)
    call HwU_fill(g_reco + 11, transverse_momentum(atop), wgts)
    call HwU_fill(g_reco + 12, rapidity(atop), wgts)
    call HwU_fill(g_reco + 13, invariant_mass(atop), wgts)
    call HwU_fill(g_reco + 14, invariant_mass(ttbar), wgts)
    call HwU_fill(g_reco + 15, transverse_momentum(ttbar), wgts)
    call HwU_fill(g_reco + 16, rapidity(ttbar), wgts)
    call HwU_fill(g_reco + 17, rapidity(top) - rapidity(atop), wgts)
    call HwU_fill(g_reco + 18, abs(rapidity(top)) - abs(rapidity(atop)), wgts)
    call HwU_fill(g_reco + 19, scattering_cosine, wgts)
    call HwU_fill(g_reco + 20, speed, wgts)
    call HwU_fill(g_reco + 21, nu(3), wgts)
    call HwU_fill(g_reco + 22, anu(3), wgts)
    call HwU_fill(g_reco + 23, invariant_mass(nunu), wgts)
    call HwU_fill(g_reco + 24, nunu(3), wgts)
    call HwU_fill(g_reco + 25, min(score, 99.999d0), wgts)
    call HwU_fill(g_reco + 26, invariant_mass(top) - invariant_mass(atop), wgts)
    call HwU_fill(g_reco + 27, dble(pairing), wgts)
    call HwU_fill(g_reco + 28, safe_ratio(transverse_momentum(top) + &
         transverse_momentum(atop), transverse_momentum(truth_top) + &
         transverse_momentum(truth_atop)), wgts)
    call HwU_fill(g_reco + 29, invariant_mass(ttbar) - &
         invariant_mass(truth_ttbar), wgts)
    call HwU_fill(g_reco + 30, wp_helicity, wgts)
    call HwU_fill(g_reco + 31, wm_helicity, wgts)
    call fill_spin(g_reco_spin, spin_values, spin_products, spin_extra, wgts)
  end subroutine fill_reconstructed


  subroutine fill_spin(offset, values, products, extra, wgts)
    implicit none
    integer, intent(in) :: offset
    double precision, intent(in) :: values(6), products(9), extra(11)
    double precision, intent(in) :: wgts(*)
    integer :: i
    do i = 1, 6
      call HwU_fill(offset + i, values(i), wgts)
    end do
    do i = 1, 9
      call HwU_fill(offset + 6 + i, products(i), wgts)
    end do
    do i = 1, 11
      call HwU_fill(offset + 15 + i, extra(i), wgts)
    end do
  end subroutine fill_spin


  subroutine top_scattering_variables(top, atop, cosine, speed)
    implicit none
    double precision, intent(in) :: top(0:3), atop(0:3)
    double precision, intent(out) :: cosine, speed
    double precision :: ttbar(0:3), top_tt(0:3), beam(0:3), beam_tt(0:3)
    double precision :: top_direction(3), beam_direction(3)

    ttbar = top + atop
    call boost_to_rest(top, ttbar, top_tt)
    beam = (/1d0, 0d0, 0d0, 1d0/)
    call boost_to_rest(beam, ttbar, beam_tt)
    call unit_vector(top_tt(1:3), top_direction)
    call unit_vector(beam_tt(1:3), beam_direction)
    cosine = clamp_cosine(dot_product(top_direction, beam_direction))
    speed = safe_ratio(spatial_magnitude(top_tt), top_tt(0))
  end subroutine top_scattering_variables


  subroutine boost_to_rest(vector, parent, result)
    implicit none
    double precision, intent(in) :: vector(0:3), parent(0:3)
    double precision, intent(out) :: result(0:3)
    double precision :: beta(3), beta2, gamma, beta_dot_p, factor

    if (parent(0) <= tiny) then
      result = vector
      return
    end if
    beta = parent(1:3)/parent(0)
    beta2 = dot_product(beta, beta)
    if (beta2 <= tiny) then
      result = vector
      return
    end if
    if (beta2 >= 1d0) beta = beta/sqrt(beta2)*(1d0 - 1d-14)
    beta2 = dot_product(beta, beta)
    gamma = 1d0/sqrt(1d0 - beta2)
    beta_dot_p = dot_product(beta, vector(1:3))
    result(0) = gamma*(vector(0) - beta_dot_p)
    factor = ((gamma - 1d0)*beta_dot_p/beta2 - gamma*vector(0))
    result(1:3) = vector(1:3) + factor*beta
  end subroutine boost_to_rest


  subroutine massless_four_vector(px, py, pz, vector)
    implicit none
    double precision, intent(in) :: px, py, pz
    double precision, intent(out) :: vector(0:3)
    vector(0) = sqrt(px**2 + py**2 + pz**2)
    vector(1) = px
    vector(2) = py
    vector(3) = pz
  end subroutine massless_four_vector


  double precision function transverse_momentum(vector)
    implicit none
    double precision, intent(in) :: vector(0:3)
    transverse_momentum = sqrt(max(0d0, vector(1)**2 + vector(2)**2))
  end function transverse_momentum


  double precision function invariant_mass2(vector)
    implicit none
    double precision, intent(in) :: vector(0:3)
    invariant_mass2 = vector(0)**2 - sum(vector(1:3)**2)
  end function invariant_mass2


  double precision function invariant_mass(vector)
    implicit none
    double precision, intent(in) :: vector(0:3)
    invariant_mass = sqrt(max(0d0, invariant_mass2(vector)))
  end function invariant_mass


  double precision function rapidity(vector)
    implicit none
    double precision, intent(in) :: vector(0:3)
    double precision :: plus, minus
    plus = vector(0) + vector(3)
    minus = vector(0) - vector(3)
    if (plus > tiny .and. minus > tiny) then
      rapidity = 0.5d0*log(plus/minus)
    else
      rapidity = sign(1d0, vector(3))*1d8
    end if
  end function rapidity


  double precision function pseudorapidity(vector)
    implicit none
    double precision, intent(in) :: vector(0:3)
    double precision :: momentum, cosine
    momentum = spatial_magnitude(vector)
    if (momentum <= tiny) then
      pseudorapidity = 0d0
      return
    end if
    cosine = max(-1d0 + 1d-15, min(1d0 - 1d-15, vector(3)/momentum))
    pseudorapidity = 0.5d0*log((1d0 + cosine)/(1d0 - cosine))
  end function pseudorapidity


  double precision function azimuth(vector)
    implicit none
    double precision, intent(in) :: vector(0:3)
    azimuth = atan2(vector(2), vector(1))
  end function azimuth


  double precision function signed_delta_phi(first, second)
    implicit none
    double precision, intent(in) :: first(0:3), second(0:3)
    signed_delta_phi = atan2(sin(azimuth(first) - azimuth(second)), &
                             cos(azimuth(first) - azimuth(second)))
  end function signed_delta_phi


  double precision function delta_phi(first, second)
    implicit none
    double precision, intent(in) :: first(0:3), second(0:3)
    delta_phi = abs(signed_delta_phi(first, second))
  end function delta_phi


  double precision function delta_phi_vector(transverse, vector)
    implicit none
    double precision, intent(in) :: transverse(0:1), vector(0:3)
    double precision :: difference
    difference = atan2(transverse(1), transverse(0)) - azimuth(vector)
    delta_phi_vector = abs(atan2(sin(difference), cos(difference)))
  end function delta_phi_vector


  double precision function delta_r(first, second)
    implicit none
    double precision, intent(in) :: first(0:3), second(0:3)
    delta_r = sqrt((pseudorapidity(first) - pseudorapidity(second))**2 + &
                   delta_phi(first, second)**2)
  end function delta_r


  double precision function spatial_magnitude(vector)
    implicit none
    double precision, intent(in) :: vector(0:3)
    spatial_magnitude = sqrt(max(0d0, sum(vector(1:3)**2)))
  end function spatial_magnitude


  double precision function spatial_opening_cosine(first, second)
    implicit none
    double precision, intent(in) :: first(0:3), second(0:3)
    spatial_opening_cosine = clamp_cosine(safe_ratio( &
         dot_product(first(1:3), second(1:3)), &
         spatial_magnitude(first)*spatial_magnitude(second)))
  end function spatial_opening_cosine


  double precision function vector_magnitude_2(vector)
    implicit none
    double precision, intent(in) :: vector(0:1)
    vector_magnitude_2 = sqrt(max(0d0, vector(0)**2 + vector(1)**2))
  end function vector_magnitude_2


  double precision function transverse_mass(visible, qx, qy)
    implicit none
    double precision, intent(in) :: visible(0:3), qx, qy
    double precision :: visible_et, invisible_pt, value
    visible_et = sqrt(max(0d0, invariant_mass2(visible) + &
                          transverse_momentum(visible)**2))
    invisible_pt = sqrt(qx**2 + qy**2)
    value = invariant_mass2(visible) + 2d0*(visible_et*invisible_pt - &
            visible(1)*qx - visible(2)*qy)
    transverse_mass = sqrt(max(0d0, value))
  end function transverse_mass


  double precision function cluster_transverse_mass(visible, met)
    implicit none
    double precision, intent(in) :: visible(0:3), met(0:1)
    double precision :: visible_et, value
    visible_et = sqrt(max(0d0, invariant_mass2(visible) + &
                          transverse_momentum(visible)**2))
    value = (visible_et + vector_magnitude_2(met))**2 - &
            (visible(1) + met(0))**2 - (visible(2) + met(1))**2
    cluster_transverse_mass = sqrt(max(0d0, value))
  end function cluster_transverse_mass


  double precision function safe_ratio(numerator, denominator)
    implicit none
    double precision, intent(in) :: numerator, denominator
    if (abs(denominator) <= tiny) then
      safe_ratio = 0d0
    else
      safe_ratio = numerator/denominator
    end if
  end function safe_ratio


  double precision function clamp_cosine(value)
    implicit none
    double precision, intent(in) :: value
    ! HwU bins are half-open.  Keep the physical endpoint inside the final bin.
    clamp_cosine = max(-1d0, min(1d0 - 1d-12, value))
  end function clamp_cosine


  subroutine unit_vector(vector, unit)
    implicit none
    double precision, intent(in) :: vector(3)
    double precision, intent(out) :: unit(3)
    double precision :: norm
    norm = sqrt(max(0d0, dot_product(vector, vector)))
    if (norm <= tiny) then
      unit = (/0d0, 0d0, 1d0/)
    else
      unit = vector/norm
    end if
  end subroutine unit_vector


  subroutine cross_product3(first, second, result)
    implicit none
    double precision, intent(in) :: first(3), second(3)
    double precision, intent(out) :: result(3)
    result(1) = first(2)*second(3) - first(3)*second(2)
    result(2) = first(3)*second(1) - first(1)*second(3)
    result(3) = first(1)*second(2) - first(2)*second(1)
  end subroutine cross_product3


  subroutine fail_analysis(message)
    implicit none
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in comprehensive ttbar dilepton analysis: '// &
                     trim(message)
    stop 1
  end subroutine fail_analysis

end module analysis_hwu_pp_ttx_dilepton_module

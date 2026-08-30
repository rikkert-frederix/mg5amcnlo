module analysis_hwu_pp_ttxw_trilepton_module
  use process_dimensions, only: nincoming
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private

  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: tiny = 1d-12
  double precision, parameter :: jet_radius = 0.4d0
  double precision, parameter :: opposite_sign_lepton_pt_cut = 10d0
  double precision, parameter :: same_sign_lepton_pt_cut = 20d0
  double precision, parameter :: electron_eta_cut = 2.47d0
  double precision, parameter :: muon_eta_cut = 2.5d0
  double precision, parameter :: bjet_pt_cut = 25d0
  double precision, parameter :: bjet_eta_cut = 2.5d0
  double precision, parameter :: analysis_jet_pt_cut = 25d0
  double precision, parameter :: analysis_jet_eta_cut = 2.5d0
  double precision, parameter :: ossf_mass_cut = 12d0
  double precision, parameter :: z_mass_reference = 91.1876d0
  double precision, parameter :: z_veto_half_width = 10d0

  integer, parameter :: g_rate = 0
  integer, parameter :: g_production = 4
  integer, parameter :: g_lepton = 34
  integer, parameter :: g_decay = 70
  integer, parameter :: g_neutrino = 94
  integer, parameter :: g_jet = 110
  integer, parameter :: g_fiducial = 138
  integer, parameter :: g_reco = 178
  integer, parameter :: g_spin = 214
  integer, parameter :: g_reco_spin = 240
  integer, parameter :: g_assoc_spin = 266
  integer, parameter :: g_diagnostic = 278
  integer, parameter :: number_of_histograms = 286

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
    integer :: ie1, ie2, imu, iep, iem, inue, ianue, inua, ib, iab
    integer :: qe1, qe2, qmu, qassoc
    integer :: ibjet, iabjet, njet, njet25, njet40, nextra25
    integer :: lead_jet, trail_jet, lead_extra_jet, lead_bjet
    integer :: assigned_radiation, pairing, stage, n_b_accepted
    integer :: jet_of_parton(size(p, 2))
    double precision :: pjet(0:3, size(p, 2))
    double precision :: e1(0:3), e2(0:3), ep(0:3), em(0:3), mu(0:3)
    double precision :: nue(0:3), anue(0:3), nua(0:3)
    double precision :: b(0:3), ab(0:3), wp(0:3), wm(0:3), wa(0:3)
    double precision :: www(0:3), top(0:3), atop(0:3), ttbar(0:3)
    double precision :: ttw(0:3), trilepton(0:3), trineutrino(0:3)
    double precision :: same_sign_first(0:3), same_sign_second(0:3)
    double precision :: opposite_sign_lepton(0:3)
    double precision :: same_sign_pair(0:3), opposite_sign_pair(0:3)
    double precision :: lead_lepton(0:3), sublead_lepton(0:3)
    double precision :: third_lepton(0:3)
    double precision :: bjet(0:3), abjet(0:3), bbjet(0:3)
    double precision :: visible(0:3), met(0:1), truth_met(0:1)
    double precision :: reco_nue(0:3), reco_anue(0:3), reco_nua(0:3)
    double precision :: reco_wp(0:3), reco_wm(0:3), reco_wa(0:3)
    double precision :: reco_top(0:3), reco_atop(0:3)
    double precision :: reco_ttbar(0:3), reco_ttw(0:3)
    double precision :: reco_lp(0:3), reco_lm(0:3), reco_la(0:3)
    double precision :: truth_score, reco_score, mt2_ee
    double precision :: min_dr_lepton_jet, min_dr_lepton_b, max_dr_lepton_b
    double precision :: lead_lepton_lead_b_dr
    double precision :: leading_extra_pt, leading_extra_y, extra_ht
    double precision :: jet_ht, lepton_ht, visible_ht, st
    double precision :: bjet_pt, abjet_pt, bjet_eta, abjet_eta
    double precision :: met_closure
    double precision :: spin_values(6), spin_products(9), spin_extra(11)
    double precision :: reco_spin_values(6), reco_spin_products(9)
    double precision :: reco_spin_extra(11)
    double precision :: assoc_spin_values(6), reco_assoc_spin_values(6)
    double precision :: wp_helicity, wm_helicity, wa_helicity
    double precision :: reco_wp_helicity, reco_wm_helicity
    double precision :: reco_wa_helicity
    logical :: fiducial, two_b_reconstruction, reco_success, pass_overlap
    logical :: theory_available, has_ossf_ee

    particle_count = size(p, 2)
    if (size(istatus) /= particle_count .or. size(ipdg) /= particle_count) then
      call fail_analysis('inconsistent momentum, status and PDG arrays')
    end if
    if (particle_count < nincoming + 8) then
      call fail_analysis('too few particles for a trileptonic ttW event')
    end if

    call find_visible_objects(p, istatus, ipdg, ie1, ie2, imu, ib, iab, &
         qe1, qe2, qmu, qassoc)
    e1 = p(0:3, ie1)
    e2 = p(0:3, ie2)
    if (transverse_momentum(e2) > transverse_momentum(e1)) then
      e1 = p(0:3, ie2)
      e2 = p(0:3, ie1)
      qe1 = lepton_charge(ipdg(ie2))
      qe2 = lepton_charge(ipdg(ie1))
    end if
    mu = p(0:3, imu)
    b = p(0:3, ib)
    ab = p(0:3, iab)
    trilepton = e1 + e2 + mu
    call classify_lepton_charges(e1, qe1, e2, qe2, mu, qmu, qassoc, &
         same_sign_first, same_sign_second, opposite_sign_lepton)
    same_sign_pair = same_sign_first + same_sign_second
    opposite_sign_pair = opposite_sign_lepton + same_sign_first
    call order_three_leptons(e1, e2, mu, lead_lepton, sublead_lepton, &
         third_lepton)
    has_ossf_ee = qe1*qe2 == -1

    call find_theory_objects(p, istatus, ipdg, qassoc, iep, iem, inue, &
         ianue, inua, theory_available)
    ep = 0d0
    em = 0d0
    nue = 0d0
    anue = 0d0
    nua = 0d0
    wp = 0d0
    wm = 0d0
    wa = 0d0
    www = 0d0
    trineutrino = 0d0
    top = 0d0
    atop = 0d0
    ttbar = 0d0
    ttw = 0d0
    truth_score = 0d0
    assigned_radiation = 0
    if (theory_available) then
      ep = p(0:3, iep)
      em = p(0:3, iem)
      nue = p(0:3, inue)
      anue = p(0:3, ianue)
      nua = p(0:3, inua)
      wp = ep + nue
      wm = em + anue
      wa = mu + nua
      www = wp + wm + wa
      trineutrino = nue + anue + nua
      call reconstruct_truth_tops(p, istatus, ipdg, wp, wm, b, ab, &
           top, atop, truth_score, assigned_radiation)
      ttbar = top + atop
      ttw = ttbar + wa
    end if

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
    visible = trilepton + bbjet

    call missing_transverse_momentum(p, istatus, ipdg, met)
    call classify_jets(pjet, njet, ibjet, iabjet, njet25, njet40, &
         nextra25, lead_jet, trail_jet, lead_extra_jet, jet_ht, extra_ht)
    leading_extra_pt = 0d0
    leading_extra_y = 0d0
    if (lead_extra_jet > 0) then
      leading_extra_pt = transverse_momentum(pjet(:, lead_extra_jet))
      leading_extra_y = rapidity(pjet(:, lead_extra_jet))
    end if
    call trilepton_jet_distances(e1, e2, mu, pjet, njet, &
         min_dr_lepton_jet, min_dr_lepton_b, max_dr_lepton_b, &
         ibjet, iabjet)
    call count_accepted_b_jets(bjet, abjet, ibjet, iabjet, &
         n_b_accepted, lead_bjet)
    if (lead_bjet == ibjet) then
      lead_lepton_lead_b_dr = delta_r(lead_lepton, bjet)
    else
      lead_lepton_lead_b_dr = delta_r(lead_lepton, abjet)
    end if

    lepton_ht = transverse_momentum(e1) + transverse_momentum(e2) + &
         transverse_momentum(mu)
    visible_ht = lepton_ht + jet_ht
    st = visible_ht + vector_magnitude_2(met)
    truth_met = 0d0
    met_closure = 0d0
    if (theory_available) then
      truth_met(0) = nue(1) + anue(1) + nua(1)
      truth_met(1) = nue(2) + anue(2) + nua(2)
      met_closure = vector_magnitude_2(met - truth_met)
      call compute_spin_observables(top, atop, ep, em, &
           spin_values, spin_products, spin_extra)
      call w_helicity_cosines(top, atop, wp, wm, ep, em, &
           wp_helicity, wm_helicity)
      call associated_w_helicity(ttw, wa, mu, wa_helicity)
      call compute_associated_spin_observables(ttw, wa, top, atop, ep, em, &
           mu, qassoc, assoc_spin_values)
      call fill_truth_production(wgts, top, atop, ttbar, wa, ttw, qassoc)
      call fill_truth_leptons(wgts, ep, em, mu, trilepton, &
           same_sign_pair, opposite_sign_pair, qassoc)
      call fill_truth_decays(wgts, top, atop, wp, wm, wa, www, ep, em, mu, &
           b, ab, qassoc, wp_helicity, wm_helicity, wa_helicity)
      call fill_truth_neutrinos(wgts, nue, anue, nua, trineutrino, met, &
           ep, em, mu, trilepton, met_closure)
      call fill_spin(g_spin, spin_values, spin_products, spin_extra, wgts)
      call fill_associated_spin(g_assoc_spin, assoc_spin_values, wgts)
    end if
    call fill_truth_jets(wgts, bjet, abjet, bbjet, pjet, njet25, njet40, &
         nextra25, lead_jet, trail_jet, leading_extra_pt, leading_extra_y, &
         extra_ht, jet_ht, lepton_ht, visible_ht, st, visible, &
         lead_lepton_lead_b_dr, min_dr_lepton_b)
    call HwU_fill(g_rate + 1, 1d0, wgts)
    call HwU_fill(g_rate + 2, dble(qassoc), wgts)
    call HwU_fill(g_rate + 3, 1d0, wgts)

    bjet_pt = transverse_momentum(bjet)
    abjet_pt = transverse_momentum(abjet)
    bjet_eta = abs(pseudorapidity(bjet))
    abjet_eta = abs(pseudorapidity(abjet))

    stage = 1
    if (transverse_momentum(same_sign_first) >= &
             same_sign_lepton_pt_cut .and. &
        transverse_momentum(same_sign_second) >= &
             same_sign_lepton_pt_cut .and. &
        transverse_momentum(opposite_sign_lepton) >= &
             opposite_sign_lepton_pt_cut .and. &
        abs(pseudorapidity(e1)) < electron_eta_cut .and. &
        abs(pseudorapidity(e2)) < electron_eta_cut .and. &
        abs(pseudorapidity(mu)) < muon_eta_cut) then
      stage = 2
      call HwU_fill(g_rate + 3, 2d0, wgts)
    end if
    if (stage == 2 .and. (.not. has_ossf_ee .or. &
        invariant_mass(e1 + e2) > ossf_mass_cut)) then
      stage = 3
      call HwU_fill(g_rate + 3, 3d0, wgts)
    end if
    if (stage == 3 .and. (.not. has_ossf_ee .or. &
        abs(invariant_mass(e1 + e2) - z_mass_reference) > &
             z_veto_half_width)) then
      stage = 4
      call HwU_fill(g_rate + 3, 4d0, wgts)
    end if
    if (stage == 4 .and. &
        abs(invariant_mass(trilepton) - z_mass_reference) > &
             z_veto_half_width) then
      stage = 5
      call HwU_fill(g_rate + 3, 5d0, wgts)
    end if
    if (stage == 5 .and. njet25 >= 1) then
      stage = 6
      call HwU_fill(g_rate + 3, 6d0, wgts)
    end if
    if (stage == 6 .and. n_b_accepted >= 1) then
      stage = 7
      call HwU_fill(g_rate + 3, 7d0, wgts)
    end if
    pass_overlap = passes_trilepton_jet_separation(e1, e2, mu, pjet, njet)
    if (stage == 7 .and. pass_overlap) then
      stage = 8
      call HwU_fill(g_rate + 3, 8d0, wgts)
    end if
    fiducial = stage >= 8
    two_b_reconstruction = fiducial .and. ibjet /= iabjet .and. &
         bjet_pt >= bjet_pt_cut .and. abjet_pt >= bjet_pt_cut .and. &
         bjet_eta < bjet_eta_cut .and. abjet_eta < bjet_eta_cut
    if (two_b_reconstruction) call HwU_fill(g_rate + 3, 9d0, wgts)
    if (fiducial) call HwU_fill(g_rate + 4, dble(qassoc), wgts)

    reco_success = .false.
    reco_score = 0d0
    pairing = 0
    reco_nue = 0d0
    reco_anue = 0d0
    reco_nua = 0d0
    reco_wp = 0d0
    reco_wm = 0d0
    reco_wa = 0d0
    reco_top = 0d0
    reco_atop = 0d0
    reco_ttbar = 0d0
    reco_ttw = 0d0
    mt2_ee = 0d0
    if (fiducial) then
      call compute_mt2(e1, e2, met, mt2_ee)
      call fill_ttw_fiducial(wgts, e1, e2, mu, trilepton, &
           same_sign_first, same_sign_second, same_sign_pair, &
           opposite_sign_lepton, opposite_sign_pair, lead_lepton, bjet, abjet, &
           bbjet, met, qassoc, n_b_accepted, njet25, nextra25, &
           jet_ht, lepton_ht, visible_ht, st, lead_lepton_lead_b_dr, mt2_ee)
    end if
    if (two_b_reconstruction) then
      call reconstruct_three_neutrinos(same_sign_first, same_sign_second, &
           opposite_sign_lepton, qassoc, bjet, abjet, met, reco_nue, &
           reco_anue, reco_nua, reco_lp, reco_lm, reco_la, pairing, &
           reco_score, reco_success)
      if (reco_success) then
        reco_wp = reco_lp + reco_nue
        reco_wm = reco_lm + reco_anue
        reco_wa = reco_la + reco_nua
        if (pairing == 1) then
          reco_top = reco_wp + bjet
          reco_atop = reco_wm + abjet
        else
          reco_top = reco_wp + abjet
          reco_atop = reco_wm + bjet
        end if
        reco_ttbar = reco_top + reco_atop
        reco_ttw = reco_ttbar + reco_wa
        call compute_spin_observables(reco_top, reco_atop, reco_lp, reco_lm, &
             reco_spin_values, reco_spin_products, reco_spin_extra)
        call w_helicity_cosines(reco_top, reco_atop, reco_wp, reco_wm, &
             reco_lp, reco_lm, reco_wp_helicity, reco_wm_helicity)
        call associated_w_helicity(reco_ttw, reco_wa, reco_la, &
             reco_wa_helicity)
        call compute_associated_spin_observables(reco_ttw, reco_wa, &
             reco_top, reco_atop, reco_lp, reco_lm, reco_la, qassoc, &
             reco_assoc_spin_values)
        call fill_ttw_reconstructed(wgts, reco_nue, reco_anue, reco_nua, &
             reco_wp, reco_wm, reco_wa, reco_top, reco_atop, reco_ttbar, &
             reco_ttw, qassoc, reco_score, pairing, reco_wp_helicity, &
             reco_wm_helicity, reco_wa_helicity)
        call fill_spin(g_reco_spin, reco_spin_values, reco_spin_products, &
             reco_spin_extra, wgts)
        call fill_associated_spin(g_assoc_spin + 6, &
             reco_assoc_spin_values, wgts)
        call HwU_fill(g_rate + 3, 10d0, wgts)
      end if
    end if

    if (theory_available) then
      call HwU_fill(g_diagnostic + 1, min(truth_score, 99.999d0), wgts)
      call HwU_fill(g_diagnostic + 2, dble(assigned_radiation), wgts)
      call HwU_fill(g_diagnostic + 3, &
           min(max(abs(invariant_mass(top) - top_mass_reference), &
                   abs(invariant_mass(atop) - top_mass_reference)), &
               49.999d0), wgts)
      call HwU_fill(g_diagnostic + 4, min(met_closure, 9.999d0), wgts)
    end if
    if (reco_success) then
      call HwU_fill(g_diagnostic + 5, min(reco_score, 99.999d0), wgts)
      call HwU_fill(g_diagnostic + 6, min(max( &
           abs(invariant_mass(reco_wp) - w_mass_reference), &
           abs(invariant_mass(reco_wm) - w_mass_reference), &
           abs(invariant_mass(reco_wa) - w_mass_reference)), 49.999d0), &
           wgts)
      call HwU_fill(g_diagnostic + 7, min(max( &
           abs(invariant_mass(reco_top) - top_mass_reference), &
           abs(invariant_mass(reco_atop) - top_mass_reference)), &
           49.999d0), wgts)
      call HwU_fill(g_diagnostic + 8, 1d0, wgts)
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
    call book(g_rate + 2, 'inclusive associated W charge', 3, -1.5d0, 1.5d0)
    call book(g_rate + 3, 'fiducial and reconstruction cutflow', 10, 0.5d0, 10.5d0)
    call book(g_rate + 4, 'fiducial charge yield summary', 3, -1.5d0, 1.5d0)

    call book(g_production + 1, 'truth top pT', 50, 0d0, 1000d0)
    call book(g_production + 2, 'truth antitop pT', 50, 0d0, 1000d0)
    call book(g_production + 3, 'truth top rapidity', 50, -5d0, 5d0)
    call book(g_production + 4, 'truth antitop rapidity', 50, -5d0, 5d0)
    call book(g_production + 5, 'truth ttbar mass', 54, 300d0, 3000d0)
    call book(g_production + 6, 'truth ttbar pT', 50, 0d0, 1000d0)
    call book(g_production + 7, 'truth ttbar rapidity', 50, -5d0, 5d0)
    call book(g_production + 8, 'truth delta rapidity top antitop', 50, -5d0, 5d0)
    call book(g_production + 9, 'truth delta abs rapidity top antitop', 50, -5d0, 5d0)
    call book(g_production + 10, 'truth delta phi top antitop', 32, 0d0, pi)
    call book(g_production + 11, 'truth delta R top antitop', 40, 0d0, 8d0)
    call book(g_production + 12, 'truth top scattering cosine', 40, -1d0, 1d0)
    call book(g_production + 13, 'truth top speed in ttbar frame', 40, 0d0, 1d0)
    call book(g_production + 14, 'truth rho two mt over mtt', 40, 0d0, 1d0)
    call book(g_production + 15, 'truth top pT asymmetry', 40, -1d0, 1d0)
    call book(g_production + 16, 'truth associated W pT', 50, 0d0, 800d0)
    call book(g_production + 17, 'truth associated W rapidity', 50, -5d0, 5d0)
    call book(g_production + 18, 'truth associated W mass', 50, 0d0, 150d0)
    call book(g_production + 19, 'truth ttW mass', 60, 400d0, 3400d0)
    call book(g_production + 20, 'truth ttW pT', 50, 0d0, 1000d0)
    call book(g_production + 21, 'truth ttW rapidity', 50, -5d0, 5d0)
    call book(g_production + 22, 'truth ttW threshold rho', 40, 0d0, 1d0)
    call book(g_production + 23, 'truth heavy system scalar pT sum', 50, 0d0, 1800d0)
    call book(g_production + 24, 'truth delta R associated W top', 40, 0d0, 8d0)
    call book(g_production + 25, 'truth delta R associated W antitop', 40, 0d0, 8d0)
    call book(g_production + 26, 'truth delta R associated W ttbar', 40, 0d0, 8d0)
    call book(g_production + 27, 'truth delta phi associated W ttbar', 32, 0d0, pi)
    call book(g_production + 28, 'truth charge signed associated W rapidity', 50, -5d0, 5d0)
    call book(g_production + 29, 'truth charge signed ttW rapidity', 50, -5d0, 5d0)
    call book(g_production + 30, 'truth charge signed delta abs top rapidity', 50, -5d0, 5d0)

    call book(g_lepton + 1, 'truth positron pT from top', 50, 0d0, 500d0)
    call book(g_lepton + 2, 'truth electron pT from antitop', 50, 0d0, 500d0)
    call book(g_lepton + 3, 'truth associated muon pT', 50, 0d0, 500d0)
    call book(g_lepton + 4, 'truth positron eta from top', 50, -5d0, 5d0)
    call book(g_lepton + 5, 'truth electron eta from antitop', 50, -5d0, 5d0)
    call book(g_lepton + 6, 'truth associated muon eta', 50, -5d0, 5d0)
    call book(g_lepton + 7, 'truth positron rapidity from top', 50, -5d0, 5d0)
    call book(g_lepton + 8, 'truth electron rapidity from antitop', 50, -5d0, 5d0)
    call book(g_lepton + 9, 'truth associated muon rapidity', 50, -5d0, 5d0)
    call book(g_lepton + 10, 'truth leading lepton pT', 50, 0d0, 500d0)
    call book(g_lepton + 11, 'truth subleading lepton pT', 50, 0d0, 400d0)
    call book(g_lepton + 12, 'truth third lepton pT', 50, 0d0, 300d0)
    call book(g_lepton + 13, 'truth trilepton mass', 50, 0d0, 1200d0)
    call book(g_lepton + 14, 'truth trilepton pT', 50, 0d0, 800d0)
    call book(g_lepton + 15, 'truth trilepton rapidity', 50, -5d0, 5d0)
    call book(g_lepton + 16, 'truth lepton HT', 50, 0d0, 1000d0)
    call book(g_lepton + 17, 'truth ee OSSF mass', 50, 0d0, 1000d0)
    call book(g_lepton + 18, 'truth ee OSSF pT', 50, 0d0, 600d0)
    call book(g_lepton + 19, 'truth ee OSSF delta phi', 32, 0d0, pi)
    call book(g_lepton + 20, 'truth ee OSSF abs delta eta', 25, 0d0, 5d0)
    call book(g_lepton + 21, 'truth ee OSSF delta R', 40, 0d0, 8d0)
    call book(g_lepton + 22, 'truth same sign pair mass', 50, 0d0, 1000d0)
    call book(g_lepton + 23, 'truth same sign pair pT', 50, 0d0, 600d0)
    call book(g_lepton + 24, 'truth same sign abs delta phi', 32, 0d0, pi)
    call book(g_lepton + 25, 'truth same sign abs delta eta', 25, 0d0, 5d0)
    call book(g_lepton + 26, 'truth same sign delta R', 40, 0d0, 8d0)
    call book(g_lepton + 27, 'truth opposite sign emu mass', 50, 0d0, 1000d0)
    call book(g_lepton + 28, 'truth opposite sign emu delta phi', 32, 0d0, pi)
    call book(g_lepton + 29, 'truth opposite sign emu abs delta eta', 25, 0d0, 5d0)
    call book(g_lepton + 30, 'truth opposite sign emu delta R', 40, 0d0, 8d0)
    call book(g_lepton + 31, 'truth minimum lepton pair mass', 50, 0d0, 600d0)
    call book(g_lepton + 32, 'truth maximum lepton pair mass', 50, 0d0, 1200d0)
    call book(g_lepton + 33, 'truth minimum lepton pair delta R', 40, 0d0, 8d0)
    call book(g_lepton + 34, 'truth maximum lepton pair delta R', 40, 0d0, 8d0)
    call book(g_lepton + 35, 'truth charge signed muon eta', 50, -5d0, 5d0)
    call book(g_lepton + 36, 'truth charge signed muon rapidity', 50, -5d0, 5d0)

    call book(g_decay + 1, 'truth top Wplus pT', 50, 0d0, 600d0)
    call book(g_decay + 2, 'truth top Wplus rapidity', 50, -5d0, 5d0)
    call book(g_decay + 3, 'truth top Wplus mass', 50, 0d0, 150d0)
    call book(g_decay + 4, 'truth antitop Wminus pT', 50, 0d0, 600d0)
    call book(g_decay + 5, 'truth antitop Wminus rapidity', 50, -5d0, 5d0)
    call book(g_decay + 6, 'truth antitop Wminus mass', 50, 0d0, 150d0)
    call book(g_decay + 7, 'truth associated W pT decay view', 50, 0d0, 800d0)
    call book(g_decay + 8, 'truth associated W rapidity decay view', 50, -5d0, 5d0)
    call book(g_decay + 9, 'truth associated W mass decay view', 50, 0d0, 150d0)
    call book(g_decay + 10, 'truth three W mass', 50, 0d0, 1800d0)
    call book(g_decay + 11, 'truth three W pT', 50, 0d0, 1000d0)
    call book(g_decay + 12, 'truth three W rapidity', 50, -5d0, 5d0)
    call book(g_decay + 13, 'truth correct positron b mass', 50, 0d0, 250d0)
    call book(g_decay + 14, 'truth correct electron bbar mass', 50, 0d0, 250d0)
    call book(g_decay + 15, 'truth maximum correct electron b mass', 50, 0d0, 250d0)
    call book(g_decay + 16, 'truth positron energy top frame', 50, 0d0, 120d0)
    call book(g_decay + 17, 'truth electron energy antitop frame', 50, 0d0, 120d0)
    call book(g_decay + 18, 'truth top Wplus helicity cosine', 40, -1d0, 1d0)
    call book(g_decay + 19, 'truth antitop Wminus helicity cosine', 40, -1d0, 1d0)
    call book(g_decay + 20, 'truth associated W helicity cosine', 40, -1d0, 1d0)
    call book(g_decay + 21, 'truth charge signed associated W helicity', 40, -1d0, 1d0)
    call book(g_decay + 22, 'truth top Wplus energy top frame', 50, 0d0, 180d0)
    call book(g_decay + 23, 'truth antitop Wminus energy antitop frame', 50, 0d0, 180d0)
    call book(g_decay + 24, 'truth muon energy associated W frame', 50, 0d0, 100d0)

    call book(g_neutrino + 1, 'truth electron neutrino pT', 50, 0d0, 500d0)
    call book(g_neutrino + 2, 'truth electron neutrino eta', 50, -7d0, 7d0)
    call book(g_neutrino + 3, 'truth electron antineutrino pT', 50, 0d0, 500d0)
    call book(g_neutrino + 4, 'truth electron antineutrino eta', 50, -7d0, 7d0)
    call book(g_neutrino + 5, 'truth associated neutrino pT', 50, 0d0, 500d0)
    call book(g_neutrino + 6, 'truth associated neutrino eta', 50, -7d0, 7d0)
    call book(g_neutrino + 7, 'truth trineutrino pT', 50, 0d0, 800d0)
    call book(g_neutrino + 8, 'truth trineutrino mass', 50, 0d0, 1200d0)
    call book(g_neutrino + 9, 'truth trineutrino rapidity', 50, -7d0, 7d0)
    call book(g_neutrino + 10, 'visible missing pT', 50, 0d0, 800d0)
    call book(g_neutrino + 11, 'minimum delta phi missing pT lepton', 32, 0d0, pi)
    call book(g_neutrino + 12, 'maximum delta phi missing pT lepton', 32, 0d0, pi)
    call book(g_neutrino + 13, 'trilepton missing cluster transverse mass', 50, 0d0, 1500d0)
    call book(g_neutrino + 14, 'leading lepton missing transverse mass', 50, 0d0, 800d0)
    call book(g_neutrino + 15, 'muon missing transverse mass', 50, 0d0, 800d0)
    call book(g_neutrino + 16, 'missing pT truth closure residual', 50, 0d0, 10d0)

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
    call book(g_jet + 11, 'jet multiplicity pT25 eta2p5', 10, -0.5d0, 9.5d0)
    call book(g_jet + 12, 'jet multiplicity pT40 eta2p5', 10, -0.5d0, 9.5d0)
    call book(g_jet + 13, 'leading jet pT', 50, 0d0, 1000d0)
    call book(g_jet + 14, 'leading jet eta', 50, -5d0, 5d0)
    call book(g_jet + 15, 'subleading jet pT', 50, 0d0, 600d0)
    call book(g_jet + 16, 'subleading jet eta', 50, -5d0, 5d0)
    call book(g_jet + 17, 'additional jet multiplicity pT25', 8, -0.5d0, 7.5d0)
    call book(g_jet + 18, 'leading additional jet pT', 50, 0d0, 800d0)
    call book(g_jet + 19, 'leading additional jet rapidity', 50, -5d0, 5d0)
    call book(g_jet + 20, 'additional jet HT', 50, 0d0, 1000d0)
    call book(g_jet + 21, 'jet HT pT25 eta2p5', 50, 0d0, 1800d0)
    call book(g_jet + 22, 'lepton HT', 50, 0d0, 1000d0)
    call book(g_jet + 23, 'visible HT leptons and jets', 50, 0d0, 2500d0)
    call book(g_jet + 24, 'visible ST including missing pT', 50, 0d0, 3000d0)
    call book(g_jet + 25, 'visible trilepton bb mass', 50, 0d0, 2500d0)
    call book(g_jet + 26, 'visible trilepton bb pT', 50, 0d0, 1200d0)
    call book(g_jet + 27, 'leading lepton leading b jet delta R', 40, 0d0, 8d0)
    call book(g_jet + 28, 'minimum lepton b jet delta R', 40, 0d0, 8d0)

    call book(g_fiducial + 1, 'fiducial rate', 1, 0.5d0, 1.5d0)
    call book(g_fiducial + 2, 'fiducial associated W charge distribution', 3, -1.5d0, 1.5d0)
    call book(g_fiducial + 3, 'fiducial leading lepton pT', 50, 0d0, 500d0)
    call book(g_fiducial + 4, 'fiducial subleading lepton pT', 50, 0d0, 400d0)
    call book(g_fiducial + 5, 'fiducial third lepton pT', 50, 0d0, 300d0)
    call book(g_fiducial + 6, 'fiducial leading electron flavor eta', 40, -2.5d0, 2.5d0)
    call book(g_fiducial + 7, 'fiducial subleading electron flavor eta', 40, -2.5d0, 2.5d0)
    call book(g_fiducial + 8, 'fiducial muon eta', 40, -2.5d0, 2.5d0)
    call book(g_fiducial + 9, 'fiducial trilepton mass', 50, 0d0, 1200d0)
    call book(g_fiducial + 10, 'fiducial trilepton pT', 50, 0d0, 800d0)
    call book(g_fiducial + 11, 'fiducial trilepton rapidity', 40, -4d0, 4d0)
    call book(g_fiducial + 12, 'fiducial lepton HT', 50, 0d0, 1000d0)
    call book(g_fiducial + 13, 'fiducial dielectron mass', 50, 0d0, 1000d0)
    call book(g_fiducial + 14, 'fiducial dielectron delta phi', 32, 0d0, pi)
    call book(g_fiducial + 15, 'fiducial dielectron delta R', 40, 0d0, 8d0)
    call book(g_fiducial + 16, 'fiducial same sign pair mass', 50, 0d0, 1000d0)
    call book(g_fiducial + 17, 'fiducial same sign abs delta phi', 32, 0d0, pi)
    call book(g_fiducial + 18, 'fiducial same sign abs delta eta', 25, 0d0, 5d0)
    call book(g_fiducial + 19, 'fiducial same sign delta R', 40, 0d0, 8d0)
    call book(g_fiducial + 20, 'fiducial leading opposite sign pair mass', 50, 0d0, 1000d0)
    call book(g_fiducial + 21, 'fiducial lead OS pair delta phi', 32, 0d0, pi)
    call book(g_fiducial + 22, 'fiducial leading opposite sign pair delta R', 40, 0d0, 8d0)
    call book(g_fiducial + 23, 'fiducial missing pT', 50, 0d0, 800d0)
    call book(g_fiducial + 24, 'fiducial min delta phi missing lepton', 32, 0d0, pi)
    call book(g_fiducial + 25, 'fiducial max delta phi missing lepton', 32, 0d0, pi)
    call book(g_fiducial + 26, 'fiducial leading lepton missing mT', 50, 0d0, 800d0)
    call book(g_fiducial + 27, 'fiducial muon missing mT', 50, 0d0, 800d0)
    call book(g_fiducial + 28, 'fiducial trilepton missing cluster mT', 50, 0d0, 1500d0)
    call book(g_fiducial + 29, 'fiducial leading b jet pT', 50, 0d0, 600d0)
    call book(g_fiducial + 30, 'fiducial trailing b jet pT', 50, 0d0, 400d0)
    call book(g_fiducial + 31, 'fiducial bb jet mass', 50, 0d0, 1200d0)
    call book(g_fiducial + 32, 'fiducial bb jet delta R', 40, 0d0, 8d0)
    call book(g_fiducial + 33, 'fiducial jet multiplicity', 10, -0.5d0, 9.5d0)
    call book(g_fiducial + 34, 'fiducial additional jet multiplicity', 8, -0.5d0, 7.5d0)
    call book(g_fiducial + 35, 'fiducial jet HT', 50, 0d0, 1800d0)
    call book(g_fiducial + 36, 'fiducial visible HT', 50, 0d0, 2500d0)
    call book(g_fiducial + 37, 'fiducial ST', 50, 0d0, 3000d0)
    call book(g_fiducial + 38, 'fiducial leading lepton leading b delta R', 40, 0d0, 8d0)
    call book(g_fiducial + 39, 'fiducial mT2 ee with total missing pT', 50, 0d0, 400d0)
    call book(g_fiducial + 40, 'fiducial accepted b jet multiplicity', 3, -0.5d0, 2.5d0)

    call book(g_reco + 1, 'reconstructed rate', 1, 0.5d0, 1.5d0)
    call book(g_reco + 2, 'reconstructed top Wplus pT', 50, 0d0, 600d0)
    call book(g_reco + 3, 'reconstructed top Wplus rapidity', 40, -4d0, 4d0)
    call book(g_reco + 4, 'reconstructed top Wplus mass', 50, 0d0, 150d0)
    call book(g_reco + 5, 'reconstructed antitop Wminus pT', 50, 0d0, 600d0)
    call book(g_reco + 6, 'reconstructed antitop Wminus rapidity', 40, -4d0, 4d0)
    call book(g_reco + 7, 'reconstructed antitop Wminus mass', 50, 0d0, 150d0)
    call book(g_reco + 8, 'reconstructed associated W pT', 50, 0d0, 800d0)
    call book(g_reco + 9, 'reconstructed associated W rapidity', 40, -4d0, 4d0)
    call book(g_reco + 10, 'reconstructed associated W mass', 50, 0d0, 150d0)
    call book(g_reco + 11, 'reconstructed top pT', 50, 0d0, 1000d0)
    call book(g_reco + 12, 'reconstructed top rapidity', 40, -4d0, 4d0)
    call book(g_reco + 13, 'reconstructed top mass', 50, 100d0, 250d0)
    call book(g_reco + 14, 'reconstructed antitop pT', 50, 0d0, 1000d0)
    call book(g_reco + 15, 'reconstructed antitop rapidity', 40, -4d0, 4d0)
    call book(g_reco + 16, 'reconstructed antitop mass', 50, 100d0, 250d0)
    call book(g_reco + 17, 'reconstructed ttbar mass', 54, 300d0, 3000d0)
    call book(g_reco + 18, 'reconstructed ttbar pT', 50, 0d0, 1000d0)
    call book(g_reco + 19, 'reconstructed ttbar rapidity', 40, -4d0, 4d0)
    call book(g_reco + 20, 'reconstructed ttW mass', 60, 400d0, 3400d0)
    call book(g_reco + 21, 'reconstructed ttW pT', 50, 0d0, 1000d0)
    call book(g_reco + 22, 'reconstructed ttW rapidity', 40, -4d0, 4d0)
    call book(g_reco + 23, 'reconstructed electron neutrino pT', 50, 0d0, 500d0)
    call book(g_reco + 24, 'reconstructed electron neutrino pz', 50, -1000d0, 1000d0)
    call book(g_reco + 25, 'reconstructed electron antineutrino pT', 50, 0d0, 500d0)
    call book(g_reco + 26, 'reconstructed electron antineutrino pz', 50, -1000d0, 1000d0)
    call book(g_reco + 27, 'reconstructed associated neutrino pT', 50, 0d0, 500d0)
    call book(g_reco + 28, 'reconstructed associated neutrino pz', 50, -1000d0, 1000d0)
    call book(g_reco + 29, 'reconstructed trineutrino mass', 50, 0d0, 1200d0)
    call book(g_reco + 30, 'reconstruction score', 50, 0d0, 100d0)
    call book(g_reco + 31, 'reconstruction b pairing', 2, 0.5d0, 2.5d0)
    call book(g_reco + 32, 'reconstructed top Wplus helicity', 40, -1d0, 1d0)
    call book(g_reco + 33, 'reconstructed antitop Wminus helicity', 40, -1d0, 1d0)
    call book(g_reco + 34, 'reconstructed associated W helicity', 40, -1d0, 1d0)
    call book(g_reco + 35, 'reco qW signed associated W helicity', 40, -1d0, 1d0)
    call book(g_reco + 36, 'reconstructed ttW threshold rho', 40, 0d0, 1d0)

    call book_spin_histograms(g_spin, 'truth ttbar spin')
    call book_spin_histograms(g_reco_spin, 'reco ttbar spin')

    call book(g_assoc_spin + 1, 'truth associated lepton beam cosine', 40, -1d0, 1d0)
    call book(g_assoc_spin + 2, 'truth qW associated lepton beam cosine', 40, -1d0, 1d0)
    call book(g_assoc_spin + 3, 'truth associated top SS lepton cosine', 40, -1d0, 1d0)
    call book(g_assoc_spin + 4, 'truth trilepton triple product', 40, -1d0, 1d0)
    call book(g_assoc_spin + 5, 'truth qW trilepton triple product', 40, -1d0, 1d0)
    call book(g_assoc_spin + 6, 'truth qW associated lepton top contrast', 40, -2d0, 2d0)
    call book(g_assoc_spin + 7, 'reco associated lepton beam cosine', 40, -1d0, 1d0)
    call book(g_assoc_spin + 8, 'reco qW associated lepton beam cosine', 40, -1d0, 1d0)
    call book(g_assoc_spin + 9, 'reco associated top SS lepton cosine', 40, -1d0, 1d0)
    call book(g_assoc_spin + 10, 'reco trilepton triple product', 40, -1d0, 1d0)
    call book(g_assoc_spin + 11, 'reco charge signed trilepton triple product', 40, -1d0, 1d0)
    call book(g_assoc_spin + 12, 'reco qW associated lepton top contrast', 40, -2d0, 2d0)

    call book(g_diagnostic + 1, 'diagnostic truth top assignment score', 50, 0d0, 100d0)
    call book(g_diagnostic + 2, 'diagnostic assigned radiation count', 6, -0.5d0, 5.5d0)
    call book(g_diagnostic + 3, 'diagnostic truth top mass residual', 50, 0d0, 50d0)
    call book(g_diagnostic + 4, 'diagnostic missing pT closure residual', 50, 0d0, 10d0)
    call book(g_diagnostic + 5, 'diagnostic reconstruction score', 50, 0d0, 100d0)
    call book(g_diagnostic + 6, 'diagnostic reco W mass residual', 50, 0d0, 50d0)
    call book(g_diagnostic + 7, 'diagnostic reco top mass residual', 50, 0d0, 50d0)
    call book(g_diagnostic + 8, 'diagnostic reconstruction rate', 1, 0.5d0, 1.5d0)

    if (g_diagnostic + 8 /= number_of_histograms) then
      call fail_analysis('internal histogram count is inconsistent')
    end if
  end subroutine book_histograms


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
    if (len_trim(title) > 43) then
      call fail_analysis('HwU base title exceeds 43 characters: '//trim(title))
    end if
    call HwU_book(label, trim(title)//' |T@NLO', bins, lower, upper)
  end subroutine book


  subroutine find_visible_objects(p, istatus, ipdg, ie1, ie2, imu, ib, iab, &
                                  qe1, qe2, qmu, qevent)
    implicit none
    double precision, intent(in) :: p(0:, :)
    integer, intent(in) :: istatus(:), ipdg(:)
    integer, intent(out) :: ie1, ie2, imu, ib, iab
    integer, intent(out) :: qe1, qe2, qmu, qevent
    integer :: i, nelectron, nmuon, nb, nab, nneutrino
    integer :: nsame, nopposite

    ie1 = 0
    ie2 = 0
    imu = 0
    ib = 0
    iab = 0
    qe1 = 0
    qe2 = 0
    qmu = 0
    nelectron = 0
    nmuon = 0
    nb = 0
    nab = 0
    nneutrino = 0
    do i = nincoming + 1, size(p, 2)
      if (istatus(i) /= 1) cycle
      if (abs(ipdg(i)) == 11) then
        nelectron = nelectron + 1
        if (nelectron == 1) then
          ie1 = i
          qe1 = lepton_charge(ipdg(i))
        else if (nelectron == 2) then
          ie2 = i
          qe2 = lepton_charge(ipdg(i))
        end if
      else if (abs(ipdg(i)) == 13) then
        nmuon = nmuon + 1
        imu = i
        qmu = lepton_charge(ipdg(i))
      else if (abs(ipdg(i)) == 12 .or. abs(ipdg(i)) == 14 .or. &
               abs(ipdg(i)) == 16) then
        nneutrino = nneutrino + 1
      else if (ipdg(i) == 5) then
        nb = nb + 1
        ib = i
      else if (ipdg(i) == -5) then
        nab = nab + 1
        iab = i
      end if
    end do
    if (nelectron /= 2 .or. nmuon /= 1 .or. nneutrino /= 3 .or. &
        nb /= 1 .or. nab /= 1) then
      call fail_analysis( &
           'requires two electrons or positrons, one muon, three neutrinos, b and bbar')
    end if
    qevent = qe1 + qe2 + qmu
    if (abs(qevent) /= 1) then
      call fail_analysis('the three charged leptons must have total charge plus or minus one')
    end if
    nsame = merge(1, 0, qe1 == qevent) + merge(1, 0, qe2 == qevent) + &
         merge(1, 0, qmu == qevent)
    nopposite = merge(1, 0, qe1 == -qevent) + &
         merge(1, 0, qe2 == -qevent) + merge(1, 0, qmu == -qevent)
    if (nsame /= 2 .or. nopposite /= 1) then
      call fail_analysis('invalid charge pattern for a trileptonic ttW event')
    end if
  end subroutine find_visible_objects


  subroutine find_theory_objects(p, istatus, ipdg, qassoc, iep, iem, inue, &
                                 ianue, inua, available)
    implicit none
    double precision, intent(in) :: p(0:, :)
    integer, intent(in) :: istatus(:), ipdg(:), qassoc
    integer, intent(out) :: iep, iem, inue, ianue, inua
    logical, intent(out) :: available
    integer :: i, nep, nem, nnue, nanue, nnua

    iep = 0
    iem = 0
    inue = 0
    ianue = 0
    inua = 0
    nep = 0
    nem = 0
    nnue = 0
    nanue = 0
    nnua = 0
    do i = nincoming + 1, size(p, 2)
      if (istatus(i) /= 1) cycle
      select case (ipdg(i))
      case (-11)
        nep = nep + 1
        iep = i
      case (11)
        nem = nem + 1
        iem = i
      case (12)
        nnue = nnue + 1
        inue = i
      case (-12)
        nanue = nanue + 1
        ianue = i
      end select
      if ((qassoc == 1 .and. ipdg(i) == 14) .or. &
          (qassoc == -1 .and. ipdg(i) == -14)) then
        nnua = nnua + 1
        inua = i
      end if
    end do
    available = nep == 1 .and. nem == 1 .and. nnue == 1 .and. &
         nanue == 1 .and. nnua == 1
  end subroutine find_theory_objects


  integer function lepton_charge(pdg)
    implicit none
    integer, intent(in) :: pdg
    if (abs(pdg) /= 11 .and. abs(pdg) /= 13) then
      call fail_analysis('requested the electric charge of a non-lepton')
    end if
    lepton_charge = -sign(1, pdg)
  end function lepton_charge


  subroutine classify_lepton_charges(e1, qe1, e2, qe2, mu, qmu, qevent, &
                                     same_first, same_second, opposite)
    implicit none
    double precision, intent(in) :: e1(0:3), e2(0:3), mu(0:3)
    integer, intent(in) :: qe1, qe2, qmu, qevent
    double precision, intent(out) :: same_first(0:3), same_second(0:3)
    double precision, intent(out) :: opposite(0:3)
    double precision :: leptons(0:3, 3), swap(0:3)
    integer :: charges(3), i, nsame

    leptons(:, 1) = e1
    leptons(:, 2) = e2
    leptons(:, 3) = mu
    charges = (/qe1, qe2, qmu/)
    nsame = 0
    opposite = 0d0
    do i = 1, 3
      if (charges(i) == qevent) then
        nsame = nsame + 1
        if (nsame == 1) same_first = leptons(:, i)
        if (nsame == 2) same_second = leptons(:, i)
      else
        opposite = leptons(:, i)
      end if
    end do
    if (transverse_momentum(same_second) > &
        transverse_momentum(same_first)) then
      swap = same_first
      same_first = same_second
      same_second = swap
    end if
  end subroutine classify_lepton_charges


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


  subroutine fill_truth_production(wgts, top, atop, ttbar, wa, ttw, qassoc)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: top(0:3), atop(0:3), ttbar(0:3)
    double precision, intent(in) :: wa(0:3), ttw(0:3)
    integer, intent(in) :: qassoc
    double precision :: pt_top, pt_atop, ytop, yatop, mtt, mttw
    double precision :: scattering_cosine, speed, threshold

    pt_top = transverse_momentum(top)
    pt_atop = transverse_momentum(atop)
    ytop = rapidity(top)
    yatop = rapidity(atop)
    mtt = invariant_mass(ttbar)
    mttw = invariant_mass(ttw)
    call top_scattering_variables(top, atop, scattering_cosine, speed)
    threshold = 2d0*top_mass_reference + w_mass_reference

    call HwU_fill(g_production + 1, pt_top, wgts)
    call HwU_fill(g_production + 2, pt_atop, wgts)
    call HwU_fill(g_production + 3, ytop, wgts)
    call HwU_fill(g_production + 4, yatop, wgts)
    call HwU_fill(g_production + 5, mtt, wgts)
    call HwU_fill(g_production + 6, transverse_momentum(ttbar), wgts)
    call HwU_fill(g_production + 7, rapidity(ttbar), wgts)
    call HwU_fill(g_production + 8, ytop - yatop, wgts)
    call HwU_fill(g_production + 9, abs(ytop) - abs(yatop), wgts)
    call HwU_fill(g_production + 10, delta_phi(top, atop), wgts)
    call HwU_fill(g_production + 11, delta_r(top, atop), wgts)
    call HwU_fill(g_production + 12, scattering_cosine, wgts)
    call HwU_fill(g_production + 13, speed, wgts)
    if (mtt > tiny) then
      call HwU_fill(g_production + 14, 2d0*top_mass_reference/mtt, wgts)
    end if
    call HwU_fill(g_production + 15, &
         safe_ratio(pt_top - pt_atop, pt_top + pt_atop), wgts)
    call HwU_fill(g_production + 16, transverse_momentum(wa), wgts)
    call HwU_fill(g_production + 17, rapidity(wa), wgts)
    call HwU_fill(g_production + 18, invariant_mass(wa), wgts)
    call HwU_fill(g_production + 19, mttw, wgts)
    call HwU_fill(g_production + 20, transverse_momentum(ttw), wgts)
    call HwU_fill(g_production + 21, rapidity(ttw), wgts)
    if (mttw > tiny) then
      call HwU_fill(g_production + 22, threshold/mttw, wgts)
    end if
    call HwU_fill(g_production + 23, pt_top + pt_atop + &
         transverse_momentum(wa), wgts)
    call HwU_fill(g_production + 24, delta_r(wa, top), wgts)
    call HwU_fill(g_production + 25, delta_r(wa, atop), wgts)
    call HwU_fill(g_production + 26, delta_r(wa, ttbar), wgts)
    call HwU_fill(g_production + 27, delta_phi(wa, ttbar), wgts)
    call HwU_fill(g_production + 28, dble(qassoc)*rapidity(wa), wgts)
    call HwU_fill(g_production + 29, dble(qassoc)*rapidity(ttw), wgts)
    call HwU_fill(g_production + 30, &
         dble(qassoc)*(abs(ytop) - abs(yatop)), wgts)
  end subroutine fill_truth_production


  subroutine fill_truth_leptons(wgts, ep, em, mu, trilepton, &
                                same_sign_pair, opposite_sign_pair, qassoc)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(in) :: trilepton(0:3)
    double precision, intent(in) :: same_sign_pair(0:3)
    double precision, intent(in) :: opposite_sign_pair(0:3)
    integer, intent(in) :: qassoc
    double precision :: lead(0:3), sublead(0:3), third(0:3)
    double precision :: same_sign_e(0:3), opposite_sign_e(0:3), ee(0:3)
    double precision :: masses(3), distances(3)

    call order_three_leptons(ep, em, mu, lead, sublead, third)
    if (qassoc > 0) then
      same_sign_e = ep
      opposite_sign_e = em
    else
      same_sign_e = em
      opposite_sign_e = ep
    end if
    ee = ep + em
    masses = (/invariant_mass(ee), invariant_mass(same_sign_pair), &
         invariant_mass(opposite_sign_pair)/)
    distances = (/delta_r(ep, em), delta_r(same_sign_e, mu), &
         delta_r(opposite_sign_e, mu)/)

    call HwU_fill(g_lepton + 1, transverse_momentum(ep), wgts)
    call HwU_fill(g_lepton + 2, transverse_momentum(em), wgts)
    call HwU_fill(g_lepton + 3, transverse_momentum(mu), wgts)
    call HwU_fill(g_lepton + 4, pseudorapidity(ep), wgts)
    call HwU_fill(g_lepton + 5, pseudorapidity(em), wgts)
    call HwU_fill(g_lepton + 6, pseudorapidity(mu), wgts)
    call HwU_fill(g_lepton + 7, rapidity(ep), wgts)
    call HwU_fill(g_lepton + 8, rapidity(em), wgts)
    call HwU_fill(g_lepton + 9, rapidity(mu), wgts)
    call HwU_fill(g_lepton + 10, transverse_momentum(lead), wgts)
    call HwU_fill(g_lepton + 11, transverse_momentum(sublead), wgts)
    call HwU_fill(g_lepton + 12, transverse_momentum(third), wgts)
    call HwU_fill(g_lepton + 13, invariant_mass(trilepton), wgts)
    call HwU_fill(g_lepton + 14, transverse_momentum(trilepton), wgts)
    call HwU_fill(g_lepton + 15, rapidity(trilepton), wgts)
    call HwU_fill(g_lepton + 16, transverse_momentum(ep) + &
         transverse_momentum(em) + transverse_momentum(mu), wgts)
    call HwU_fill(g_lepton + 17, invariant_mass(ee), wgts)
    call HwU_fill(g_lepton + 18, transverse_momentum(ee), wgts)
    call HwU_fill(g_lepton + 19, delta_phi(ep, em), wgts)
    call HwU_fill(g_lepton + 20, &
         abs(pseudorapidity(ep) - pseudorapidity(em)), wgts)
    call HwU_fill(g_lepton + 21, delta_r(ep, em), wgts)
    call HwU_fill(g_lepton + 22, invariant_mass(same_sign_pair), wgts)
    call HwU_fill(g_lepton + 23, transverse_momentum(same_sign_pair), wgts)
    call HwU_fill(g_lepton + 24, delta_phi(same_sign_e, mu), wgts)
    call HwU_fill(g_lepton + 25, &
         abs(pseudorapidity(same_sign_e) - pseudorapidity(mu)), wgts)
    call HwU_fill(g_lepton + 26, delta_r(same_sign_e, mu), wgts)
    call HwU_fill(g_lepton + 27, invariant_mass(opposite_sign_pair), wgts)
    call HwU_fill(g_lepton + 28, delta_phi(opposite_sign_e, mu), wgts)
    call HwU_fill(g_lepton + 29, &
         abs(pseudorapidity(opposite_sign_e) - pseudorapidity(mu)), wgts)
    call HwU_fill(g_lepton + 30, delta_r(opposite_sign_e, mu), wgts)
    call HwU_fill(g_lepton + 31, minval(masses), wgts)
    call HwU_fill(g_lepton + 32, maxval(masses), wgts)
    call HwU_fill(g_lepton + 33, minval(distances), wgts)
    call HwU_fill(g_lepton + 34, maxval(distances), wgts)
    call HwU_fill(g_lepton + 35, &
         dble(qassoc)*pseudorapidity(mu), wgts)
    call HwU_fill(g_lepton + 36, dble(qassoc)*rapidity(mu), wgts)
  end subroutine fill_truth_leptons


  subroutine fill_truth_decays(wgts, top, atop, wp, wm, wa, www, ep, em, &
                               mu, b, ab, qassoc, wp_helicity, wm_helicity, &
                               wa_helicity)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: top(0:3), atop(0:3)
    double precision, intent(in) :: wp(0:3), wm(0:3), wa(0:3), www(0:3)
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(in) :: b(0:3), ab(0:3)
    double precision, intent(in) :: wp_helicity, wm_helicity, wa_helicity
    integer, intent(in) :: qassoc
    double precision :: ep_top(0:3), em_atop(0:3), wp_top(0:3)
    double precision :: wm_atop(0:3), mu_wa(0:3), mlb_plus, mlb_minus

    call boost_to_rest(ep, top, ep_top)
    call boost_to_rest(em, atop, em_atop)
    call boost_to_rest(wp, top, wp_top)
    call boost_to_rest(wm, atop, wm_atop)
    call boost_to_rest(mu, wa, mu_wa)
    mlb_plus = invariant_mass(ep + b)
    mlb_minus = invariant_mass(em + ab)

    call HwU_fill(g_decay + 1, transverse_momentum(wp), wgts)
    call HwU_fill(g_decay + 2, rapidity(wp), wgts)
    call HwU_fill(g_decay + 3, invariant_mass(wp), wgts)
    call HwU_fill(g_decay + 4, transverse_momentum(wm), wgts)
    call HwU_fill(g_decay + 5, rapidity(wm), wgts)
    call HwU_fill(g_decay + 6, invariant_mass(wm), wgts)
    call HwU_fill(g_decay + 7, transverse_momentum(wa), wgts)
    call HwU_fill(g_decay + 8, rapidity(wa), wgts)
    call HwU_fill(g_decay + 9, invariant_mass(wa), wgts)
    call HwU_fill(g_decay + 10, invariant_mass(www), wgts)
    call HwU_fill(g_decay + 11, transverse_momentum(www), wgts)
    call HwU_fill(g_decay + 12, rapidity(www), wgts)
    call HwU_fill(g_decay + 13, mlb_plus, wgts)
    call HwU_fill(g_decay + 14, mlb_minus, wgts)
    call HwU_fill(g_decay + 15, max(mlb_plus, mlb_minus), wgts)
    call HwU_fill(g_decay + 16, ep_top(0), wgts)
    call HwU_fill(g_decay + 17, em_atop(0), wgts)
    call HwU_fill(g_decay + 18, wp_helicity, wgts)
    call HwU_fill(g_decay + 19, wm_helicity, wgts)
    call HwU_fill(g_decay + 20, wa_helicity, wgts)
    call HwU_fill(g_decay + 21, dble(qassoc)*wa_helicity, wgts)
    call HwU_fill(g_decay + 22, wp_top(0), wgts)
    call HwU_fill(g_decay + 23, wm_atop(0), wgts)
    call HwU_fill(g_decay + 24, mu_wa(0), wgts)
  end subroutine fill_truth_decays


  subroutine fill_truth_neutrinos(wgts, nue, anue, nua, trineutrino, met, &
                                  ep, em, mu, trilepton, met_closure)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: nue(0:3), anue(0:3), nua(0:3)
    double precision, intent(in) :: trineutrino(0:3), met(0:1)
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(in) :: trilepton(0:3), met_closure
    double precision :: lead(0:3), sublead(0:3), third(0:3), dphis(3)

    call order_three_leptons(ep, em, mu, lead, sublead, third)
    dphis = (/delta_phi_vector(met, ep), delta_phi_vector(met, em), &
         delta_phi_vector(met, mu)/)
    call HwU_fill(g_neutrino + 1, transverse_momentum(nue), wgts)
    call HwU_fill(g_neutrino + 2, pseudorapidity(nue), wgts)
    call HwU_fill(g_neutrino + 3, transverse_momentum(anue), wgts)
    call HwU_fill(g_neutrino + 4, pseudorapidity(anue), wgts)
    call HwU_fill(g_neutrino + 5, transverse_momentum(nua), wgts)
    call HwU_fill(g_neutrino + 6, pseudorapidity(nua), wgts)
    call HwU_fill(g_neutrino + 7, transverse_momentum(trineutrino), wgts)
    call HwU_fill(g_neutrino + 8, invariant_mass(trineutrino), wgts)
    call HwU_fill(g_neutrino + 9, rapidity(trineutrino), wgts)
    call HwU_fill(g_neutrino + 10, vector_magnitude_2(met), wgts)
    call HwU_fill(g_neutrino + 11, minval(dphis), wgts)
    call HwU_fill(g_neutrino + 12, maxval(dphis), wgts)
    call HwU_fill(g_neutrino + 13, &
         cluster_transverse_mass(trilepton, met), wgts)
    call HwU_fill(g_neutrino + 14, &
         transverse_mass(lead, met(0), met(1)), wgts)
    call HwU_fill(g_neutrino + 15, &
         transverse_mass(mu, met(0), met(1)), wgts)
    call HwU_fill(g_neutrino + 16, met_closure, wgts)
  end subroutine fill_truth_neutrinos


  subroutine fill_truth_jets(wgts, bjet, abjet, bbjet, pjet, njet25, &
                             njet40, nextra25, lead_jet, trail_jet, &
                             leading_extra_pt, leading_extra_y, extra_ht, &
                             jet_ht, lepton_ht, visible_ht, st, visible, &
                             lead_lepton_lead_b_dr, min_dr_lepton_b)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), bbjet(0:3)
    double precision, intent(in) :: pjet(0:, :)
    integer, intent(in) :: njet25, njet40, nextra25
    integer, intent(in) :: lead_jet, trail_jet
    double precision, intent(in) :: leading_extra_pt, leading_extra_y
    double precision, intent(in) :: extra_ht, jet_ht, lepton_ht
    double precision, intent(in) :: visible_ht, st, visible(0:3)
    double precision, intent(in) :: lead_lepton_lead_b_dr, min_dr_lepton_b
    double precision :: bpt, abpt

    bpt = transverse_momentum(bjet)
    abpt = transverse_momentum(abjet)
    call HwU_fill(g_jet + 1, bpt, wgts)
    call HwU_fill(g_jet + 2, abpt, wgts)
    call HwU_fill(g_jet + 3, max(bpt, abpt), wgts)
    call HwU_fill(g_jet + 4, min(bpt, abpt), wgts)
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
    call HwU_fill(g_jet + 22, lepton_ht, wgts)
    call HwU_fill(g_jet + 23, visible_ht, wgts)
    call HwU_fill(g_jet + 24, st, wgts)
    call HwU_fill(g_jet + 25, invariant_mass(visible), wgts)
    call HwU_fill(g_jet + 26, transverse_momentum(visible), wgts)
    call HwU_fill(g_jet + 27, lead_lepton_lead_b_dr, wgts)
    call HwU_fill(g_jet + 28, min_dr_lepton_b, wgts)
  end subroutine fill_truth_jets


  subroutine fill_ttw_fiducial(wgts, e1, e2, mu, trilepton, &
                               same_sign_first, same_sign_second, &
                               same_sign_pair, opposite_sign_lepton, &
                               opposite_sign_pair, &
                               lead_lepton, bjet, abjet, bbjet, met, &
                               qassoc, n_b_accepted, njet25, nextra25, &
                               jet_ht, lepton_ht, visible_ht, st, &
                               lead_lepton_lead_b_dr, mt2_ee)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: e1(0:3), e2(0:3), mu(0:3)
    double precision, intent(in) :: trilepton(0:3)
    double precision, intent(in) :: same_sign_first(0:3)
    double precision, intent(in) :: same_sign_second(0:3)
    double precision, intent(in) :: same_sign_pair(0:3)
    double precision, intent(in) :: opposite_sign_lepton(0:3)
    double precision, intent(in) :: opposite_sign_pair(0:3)
    double precision, intent(in) :: lead_lepton(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), bbjet(0:3)
    double precision, intent(in) :: met(0:1)
    integer, intent(in) :: qassoc, n_b_accepted, njet25, nextra25
    double precision, intent(in) :: jet_ht, lepton_ht, visible_ht, st
    double precision, intent(in) :: lead_lepton_lead_b_dr, mt2_ee
    double precision :: lead(0:3), sublead(0:3), third(0:3)
    double precision :: dphis(3)
    double precision :: bpt, abpt

    call order_three_leptons(e1, e2, mu, lead, sublead, third)
    dphis = (/delta_phi_vector(met, e1), delta_phi_vector(met, e2), &
         delta_phi_vector(met, mu)/)
    bpt = transverse_momentum(bjet)
    abpt = transverse_momentum(abjet)

    call HwU_fill(g_fiducial + 1, 1d0, wgts)
    call HwU_fill(g_fiducial + 2, dble(qassoc), wgts)
    call HwU_fill(g_fiducial + 3, transverse_momentum(lead), wgts)
    call HwU_fill(g_fiducial + 4, transverse_momentum(sublead), wgts)
    call HwU_fill(g_fiducial + 5, transverse_momentum(third), wgts)
    call HwU_fill(g_fiducial + 6, pseudorapidity(e1), wgts)
    call HwU_fill(g_fiducial + 7, pseudorapidity(e2), wgts)
    call HwU_fill(g_fiducial + 8, pseudorapidity(mu), wgts)
    call HwU_fill(g_fiducial + 9, invariant_mass(trilepton), wgts)
    call HwU_fill(g_fiducial + 10, transverse_momentum(trilepton), wgts)
    call HwU_fill(g_fiducial + 11, rapidity(trilepton), wgts)
    call HwU_fill(g_fiducial + 12, lepton_ht, wgts)
    call HwU_fill(g_fiducial + 13, invariant_mass(e1 + e2), wgts)
    call HwU_fill(g_fiducial + 14, delta_phi(e1, e2), wgts)
    call HwU_fill(g_fiducial + 15, delta_r(e1, e2), wgts)
    call HwU_fill(g_fiducial + 16, invariant_mass(same_sign_pair), wgts)
    call HwU_fill(g_fiducial + 17, &
         delta_phi(same_sign_first, same_sign_second), wgts)
    call HwU_fill(g_fiducial + 18, &
         abs(pseudorapidity(same_sign_first) - &
             pseudorapidity(same_sign_second)), wgts)
    call HwU_fill(g_fiducial + 19, &
         delta_r(same_sign_first, same_sign_second), wgts)
    call HwU_fill(g_fiducial + 20, invariant_mass(opposite_sign_pair), wgts)
    call HwU_fill(g_fiducial + 21, &
         delta_phi(opposite_sign_lepton, same_sign_first), wgts)
    call HwU_fill(g_fiducial + 22, &
         delta_r(opposite_sign_lepton, same_sign_first), wgts)
    call HwU_fill(g_fiducial + 23, vector_magnitude_2(met), wgts)
    call HwU_fill(g_fiducial + 24, minval(dphis), wgts)
    call HwU_fill(g_fiducial + 25, maxval(dphis), wgts)
    call HwU_fill(g_fiducial + 26, &
         transverse_mass(lead_lepton, met(0), met(1)), wgts)
    call HwU_fill(g_fiducial + 27, &
         transverse_mass(mu, met(0), met(1)), wgts)
    call HwU_fill(g_fiducial + 28, &
         cluster_transverse_mass(trilepton, met), wgts)
    call HwU_fill(g_fiducial + 29, max(bpt, abpt), wgts)
    call HwU_fill(g_fiducial + 30, min(bpt, abpt), wgts)
    call HwU_fill(g_fiducial + 31, invariant_mass(bbjet), wgts)
    call HwU_fill(g_fiducial + 32, delta_r(bjet, abjet), wgts)
    call HwU_fill(g_fiducial + 33, dble(njet25), wgts)
    call HwU_fill(g_fiducial + 34, dble(nextra25), wgts)
    call HwU_fill(g_fiducial + 35, jet_ht, wgts)
    call HwU_fill(g_fiducial + 36, visible_ht, wgts)
    call HwU_fill(g_fiducial + 37, st, wgts)
    call HwU_fill(g_fiducial + 38, lead_lepton_lead_b_dr, wgts)
    call HwU_fill(g_fiducial + 39, mt2_ee, wgts)
    call HwU_fill(g_fiducial + 40, dble(n_b_accepted), wgts)
  end subroutine fill_ttw_fiducial


  subroutine fill_ttw_reconstructed(wgts, nue, anue, nua, wp, wm, wa, &
                                    top, atop, ttbar, ttw, qassoc, score, &
                                    pairing, wp_helicity, wm_helicity, &
                                    wa_helicity)
    implicit none
    double precision, intent(in) :: wgts(*)
    double precision, intent(in) :: nue(0:3), anue(0:3), nua(0:3)
    double precision, intent(in) :: wp(0:3), wm(0:3), wa(0:3)
    double precision, intent(in) :: top(0:3), atop(0:3), ttbar(0:3)
    double precision, intent(in) :: ttw(0:3), score
    double precision, intent(in) :: wp_helicity, wm_helicity, wa_helicity
    integer, intent(in) :: qassoc, pairing
    double precision :: trineutrino(0:3), mttw

    trineutrino = nue + anue + nua
    mttw = invariant_mass(ttw)
    call HwU_fill(g_reco + 1, 1d0, wgts)
    call HwU_fill(g_reco + 2, transverse_momentum(wp), wgts)
    call HwU_fill(g_reco + 3, rapidity(wp), wgts)
    call HwU_fill(g_reco + 4, invariant_mass(wp), wgts)
    call HwU_fill(g_reco + 5, transverse_momentum(wm), wgts)
    call HwU_fill(g_reco + 6, rapidity(wm), wgts)
    call HwU_fill(g_reco + 7, invariant_mass(wm), wgts)
    call HwU_fill(g_reco + 8, transverse_momentum(wa), wgts)
    call HwU_fill(g_reco + 9, rapidity(wa), wgts)
    call HwU_fill(g_reco + 10, invariant_mass(wa), wgts)
    call HwU_fill(g_reco + 11, transverse_momentum(top), wgts)
    call HwU_fill(g_reco + 12, rapidity(top), wgts)
    call HwU_fill(g_reco + 13, invariant_mass(top), wgts)
    call HwU_fill(g_reco + 14, transverse_momentum(atop), wgts)
    call HwU_fill(g_reco + 15, rapidity(atop), wgts)
    call HwU_fill(g_reco + 16, invariant_mass(atop), wgts)
    call HwU_fill(g_reco + 17, invariant_mass(ttbar), wgts)
    call HwU_fill(g_reco + 18, transverse_momentum(ttbar), wgts)
    call HwU_fill(g_reco + 19, rapidity(ttbar), wgts)
    call HwU_fill(g_reco + 20, mttw, wgts)
    call HwU_fill(g_reco + 21, transverse_momentum(ttw), wgts)
    call HwU_fill(g_reco + 22, rapidity(ttw), wgts)
    call HwU_fill(g_reco + 23, transverse_momentum(nue), wgts)
    call HwU_fill(g_reco + 24, nue(3), wgts)
    call HwU_fill(g_reco + 25, transverse_momentum(anue), wgts)
    call HwU_fill(g_reco + 26, anue(3), wgts)
    call HwU_fill(g_reco + 27, transverse_momentum(nua), wgts)
    call HwU_fill(g_reco + 28, nua(3), wgts)
    call HwU_fill(g_reco + 29, invariant_mass(trineutrino), wgts)
    call HwU_fill(g_reco + 30, min(score, 99.999d0), wgts)
    call HwU_fill(g_reco + 31, dble(pairing), wgts)
    call HwU_fill(g_reco + 32, wp_helicity, wgts)
    call HwU_fill(g_reco + 33, wm_helicity, wgts)
    call HwU_fill(g_reco + 34, wa_helicity, wgts)
    call HwU_fill(g_reco + 35, dble(qassoc)*wa_helicity, wgts)
    if (mttw > tiny) then
      call HwU_fill(g_reco + 36, &
           (2d0*top_mass_reference + w_mass_reference)/mttw, wgts)
    end if
  end subroutine fill_ttw_reconstructed


  subroutine order_three_leptons(ep, em, mu, lead, sublead, third)
    implicit none
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(out) :: lead(0:3), sublead(0:3), third(0:3)
    double precision :: vectors(0:3, 3), pts(3), swap_vector(0:3)
    double precision :: swap_pt
    integer :: i, j

    vectors(:, 1) = ep
    vectors(:, 2) = em
    vectors(:, 3) = mu
    do i = 1, 3
      pts(i) = transverse_momentum(vectors(:, i))
    end do
    do i = 1, 2
      do j = i + 1, 3
        if (pts(j) > pts(i)) then
          swap_pt = pts(i)
          pts(i) = pts(j)
          pts(j) = swap_pt
          swap_vector = vectors(:, i)
          vectors(:, i) = vectors(:, j)
          vectors(:, j) = swap_vector
        end if
      end do
    end do
    lead = vectors(:, 1)
    sublead = vectors(:, 2)
    third = vectors(:, 3)
  end subroutine order_three_leptons


  subroutine trilepton_jet_distances(ep, em, mu, pjet, njet, min_all, &
                                     min_b, max_b, ibjet, iabjet)
    implicit none
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(in) :: pjet(0:, :)
    integer, intent(in) :: njet, ibjet, iabjet
    double precision, intent(out) :: min_all, min_b, max_b
    double precision :: leptons(0:3, 3), distance
    integer :: i, j

    leptons(:, 1) = ep
    leptons(:, 2) = em
    leptons(:, 3) = mu
    min_all = huge(1d0)
    min_b = huge(1d0)
    max_b = 0d0
    do i = 1, 3
      do j = 1, njet
        if (transverse_momentum(pjet(:, j)) < analysis_jet_pt_cut) cycle
        if (abs(pseudorapidity(pjet(:, j))) >= analysis_jet_eta_cut) cycle
        distance = delta_r(leptons(:, i), pjet(:, j))
        min_all = min(min_all, distance)
        if (j == ibjet .or. j == iabjet) then
          min_b = min(min_b, distance)
          max_b = max(max_b, distance)
        end if
      end do
    end do
    if (min_all > huge(1d0)/2d0) min_all = 99d0
    if (min_b > huge(1d0)/2d0) min_b = 99d0
  end subroutine trilepton_jet_distances


  subroutine count_accepted_b_jets(bjet, abjet, ibjet, iabjet, count, lead)
    implicit none
    double precision, intent(in) :: bjet(0:3), abjet(0:3)
    integer, intent(in) :: ibjet, iabjet
    integer, intent(out) :: count, lead
    logical :: accept_b, accept_ab

    accept_b = transverse_momentum(bjet) >= bjet_pt_cut .and. &
         abs(pseudorapidity(bjet)) < bjet_eta_cut
    accept_ab = transverse_momentum(abjet) >= bjet_pt_cut .and. &
         abs(pseudorapidity(abjet)) < bjet_eta_cut
    if (ibjet == iabjet) then
      count = merge(1, 0, accept_b .or. accept_ab)
      lead = ibjet
    else
      count = merge(1, 0, accept_b) + merge(1, 0, accept_ab)
      if (accept_b .and. .not. accept_ab) then
        lead = ibjet
      else if (accept_ab .and. .not. accept_b) then
        lead = iabjet
      else if (transverse_momentum(bjet) >= transverse_momentum(abjet)) then
        lead = ibjet
      else
        lead = iabjet
      end if
    end if
  end subroutine count_accepted_b_jets


  logical function passes_trilepton_jet_separation(ep, em, mu, pjet, njet)
    implicit none
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(in) :: pjet(0:, :)
    integer, intent(in) :: njet
    double precision :: leptons(0:3, 3), radius
    integer :: i, j

    leptons(:, 1) = ep
    leptons(:, 2) = em
    leptons(:, 3) = mu
    passes_trilepton_jet_separation = .true.
    do i = 1, 3
      radius = min(0.4d0, 0.04d0 + &
           10d0/max(transverse_momentum(leptons(:, i)), tiny))
      do j = 1, njet
        if (transverse_momentum(pjet(:, j)) < analysis_jet_pt_cut) cycle
        if (abs(pseudorapidity(pjet(:, j))) >= analysis_jet_eta_cut) cycle
        if (delta_r(leptons(:, i), pjet(:, j)) < radius) then
          passes_trilepton_jet_separation = .false.
          return
        end if
      end do
    end do
  end function passes_trilepton_jet_separation


  subroutine reconstruct_three_neutrinos(same_first, same_second, opposite, &
                                         qassoc, bjet, abjet, met, nue, &
                                         anue, nua, lp, lm, la, pairing, &
                                         best_score, success)
    implicit none
    double precision, intent(in) :: same_first(0:3), same_second(0:3)
    double precision, intent(in) :: opposite(0:3), bjet(0:3), abjet(0:3)
    double precision, intent(in) :: met(0:1)
    integer, intent(in) :: qassoc
    double precision, intent(out) :: nue(0:3), anue(0:3), nua(0:3)
    double precision, intent(out) :: lp(0:3), lm(0:3), la(0:3)
    integer, intent(out) :: pairing
    double precision, intent(out) :: best_score
    logical, intent(out) :: success
    double precision :: trial_lp(0:3), trial_lm(0:3), trial_la(0:3)
    double precision :: trial_nue(0:3), trial_anue(0:3), trial_nua(0:3)
    double precision :: trial_score
    integer :: assignment, trial_pairing
    logical :: trial_success

    best_score = huge(1d0)
    pairing = 0
    nue = 0d0
    anue = 0d0
    nua = 0d0
    lp = 0d0
    lm = 0d0
    la = 0d0
    do assignment = 1, 2
      if (assignment == 1) then
        trial_la = same_first
        if (qassoc > 0) then
          trial_lp = same_second
          trial_lm = opposite
        else
          trial_lp = opposite
          trial_lm = same_second
        end if
      else
        trial_la = same_second
        if (qassoc > 0) then
          trial_lp = same_first
          trial_lm = opposite
        else
          trial_lp = opposite
          trial_lm = same_first
        end if
      end if
      call reconstruct_neutrino_assignment(trial_lp, trial_lm, trial_la, &
           bjet, abjet, met, trial_nue, trial_anue, trial_nua, &
           trial_pairing, trial_score, trial_success)
      if (trial_success .and. trial_score < best_score) then
        best_score = trial_score
        pairing = trial_pairing
        nue = trial_nue
        anue = trial_anue
        nua = trial_nua
        lp = trial_lp
        lm = trial_lm
        la = trial_la
      end if
    end do
    success = pairing > 0 .and. best_score < huge(1d0)/2d0
  end subroutine reconstruct_three_neutrinos


  subroutine reconstruct_neutrino_assignment(ep, em, mu, bjet, abjet, met, &
                                              nue, anue, nua, pairing, &
                                              best_score, success)
    implicit none
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), met(0:1)
    double precision, intent(out) :: nue(0:3), anue(0:3), nua(0:3)
    integer, intent(out) :: pairing
    double precision, intent(out) :: best_score
    logical, intent(out) :: success
    double precision :: trial_nue(0:3), trial_anue(0:3), trial_nua(0:3)
    double precision :: q(4), best_q(4), trial_q(4), scale, step
    double precision :: score, directions(4, 8)
    integer :: i1, i2, i3, i4, direction, iteration, trial_pairing, dim
    logical :: improved

    directions = 0d0
    do dim = 1, 4
      directions(dim, 2*dim - 1) = 1d0
      directions(dim, 2*dim) = -1d0
    end do
    q = (/met(0)/3d0, met(1)/3d0, met(0)/3d0, met(1)/3d0/)
    scale = max(30d0, 0.5d0*vector_magnitude_2(met))
    best_score = huge(1d0)
    pairing = 0
    nue = 0d0
    anue = 0d0
    nua = 0d0

    do i1 = -1, 1
      do i2 = -1, 1
        do i3 = -1, 1
          do i4 = -1, 1
            trial_q = q + scale*dble((/i1, i2, i3, i4/))
            call three_neutrino_partition_score(ep, em, mu, bjet, abjet, &
                 met, trial_q, trial_nue, trial_anue, trial_nua, &
                 trial_pairing, score)
            if (score < best_score) then
              best_score = score
              best_q = trial_q
              nue = trial_nue
              anue = trial_anue
              nua = trial_nua
              pairing = trial_pairing
            end if
          end do
        end do
      end do
    end do

    step = scale
    do iteration = 1, 12
      improved = .false.
      do direction = 1, 8
        trial_q = best_q + step*directions(:, direction)
        call three_neutrino_partition_score(ep, em, mu, bjet, abjet, &
             met, trial_q, trial_nue, trial_anue, trial_nua, &
             trial_pairing, score)
        if (score < best_score) then
          best_score = score
          best_q = trial_q
          nue = trial_nue
          anue = trial_anue
          nua = trial_nua
          pairing = trial_pairing
          improved = .true.
        end if
      end do
      if (.not. improved) step = 0.5d0*step
    end do
    success = pairing > 0 .and. best_score < huge(1d0)/2d0
  end subroutine reconstruct_neutrino_assignment


  subroutine three_neutrino_partition_score(ep, em, mu, bjet, abjet, &
                                             met, q, best_nue, best_anue, &
                                             best_nua, best_pairing, &
                                             best_score)
    implicit none
    double precision, intent(in) :: ep(0:3), em(0:3), mu(0:3)
    double precision, intent(in) :: bjet(0:3), abjet(0:3), met(0:1)
    double precision, intent(in) :: q(4)
    double precision, intent(out) :: best_nue(0:3), best_anue(0:3)
    double precision, intent(out) :: best_nua(0:3), best_score
    integer, intent(out) :: best_pairing
    double precision :: roots_ep(2), roots_em(2), roots_mu(2)
    double precision :: cand_nue(0:3), cand_anue(0:3), cand_nua(0:3)
    double precision :: wp(0:3), wm(0:3), wa(0:3)
    double precision :: top(0:3), atop(0:3), score, pz_sum, pt_sum
    double precision :: qax, qay
    integer :: nep, nem, nmu, ie, im, ia, pair

    qax = met(0) - q(1) - q(3)
    qay = met(1) - q(2) - q(4)
    call w_constraint_pz(ep, q(1), q(2), roots_ep, nep)
    call w_constraint_pz(em, q(3), q(4), roots_em, nem)
    call w_constraint_pz(mu, qax, qay, roots_mu, nmu)
    best_score = huge(1d0)
    best_pairing = 0
    best_nue = 0d0
    best_anue = 0d0
    best_nua = 0d0
    if (nep == 0 .or. nem == 0 .or. nmu == 0) return

    do ie = 1, nep
      call massless_four_vector(q(1), q(2), roots_ep(ie), cand_nue)
      wp = ep + cand_nue
      do im = 1, nem
        call massless_four_vector(q(3), q(4), roots_em(im), cand_anue)
        wm = em + cand_anue
        do ia = 1, nmu
          call massless_four_vector(qax, qay, roots_mu(ia), cand_nua)
          wa = mu + cand_nua
          do pair = 1, 2
            if (pair == 1) then
              top = wp + bjet
              atop = wm + abjet
            else
              top = wp + abjet
              atop = wm + bjet
            end if
            pz_sum = abs(cand_nue(3)) + abs(cand_anue(3)) + &
                 abs(cand_nua(3))
            pt_sum = transverse_momentum(cand_nue) + &
                 transverse_momentum(cand_anue) + &
                 transverse_momentum(cand_nua)
            score = ((invariant_mass(top) - top_mass_reference)/15d0)**2 + &
                 ((invariant_mass(atop) - top_mass_reference)/15d0)**2 + &
                 ((invariant_mass(wp) - w_mass_reference)/10d0)**2 + &
                 ((invariant_mass(wm) - w_mass_reference)/10d0)**2 + &
                 ((invariant_mass(wa) - w_mass_reference)/10d0)**2 + &
                 ((invariant_mass(top) - invariant_mass(atop))/30d0)**2 + &
                 (pz_sum/1500d0)**2 + (pt_sum/2000d0)**2
            if (score < best_score) then
              best_score = score
              best_pairing = pair
              best_nue = cand_nue
              best_anue = cand_anue
              best_nua = cand_nua
            end if
          end do
        end do
      end do
    end do
  end subroutine three_neutrino_partition_score


  subroutine associated_w_helicity(ttw, wa, charged_lepton, cosine)
    implicit none
    double precision, intent(in) :: ttw(0:3), wa(0:3), charged_lepton(0:3)
    double precision, intent(out) :: cosine
    double precision :: wa_parent(0:3), mu_w(0:3)
    double precision :: axis(3), mu_axis(3)

    call boost_to_rest(wa, ttw, wa_parent)
    call boost_to_rest(charged_lepton, wa, mu_w)
    call unit_vector(wa_parent(1:3), axis)
    call unit_vector(mu_w(1:3), mu_axis)
    cosine = clamp_cosine(dot_product(axis, mu_axis))
  end subroutine associated_w_helicity


  subroutine compute_associated_spin_observables(ttw, wa, top, atop, &
                                                  lp, lm, la, qassoc, values)
    implicit none
    double precision, intent(in) :: ttw(0:3), wa(0:3), top(0:3), atop(0:3)
    double precision, intent(in) :: lp(0:3), lm(0:3), la(0:3)
    integer, intent(in) :: qassoc
    double precision, intent(out) :: values(6)
    double precision :: beam(0:3), beam_w(0:3), mu_w(0:3)
    double precision :: top_frame(0:3), atop_frame(0:3), mu_frame(0:3)
    double precision :: ep_axis(3), em_axis(3), mu_axis(3), beam_axis(3)
    double precision :: cross(3), triple
    double precision :: same_sign_top_lepton(0:3)

    beam = (/1d0, 0d0, 0d0, 1d0/)
    call boost_to_rest(beam, wa, beam_w)
    call boost_to_rest(la, wa, mu_w)
    call unit_vector(beam_w(1:3), beam_axis)
    call unit_vector(mu_w(1:3), mu_axis)
    values(1) = clamp_cosine(dot_product(beam_axis, mu_axis))
    values(2) = dble(qassoc)*values(1)

    if (qassoc > 0) then
      same_sign_top_lepton = lp
    else
      same_sign_top_lepton = lm
    end if
    values(3) = spatial_opening_cosine(la, same_sign_top_lepton)
    call unit_vector(lp(1:3), ep_axis)
    call unit_vector(lm(1:3), em_axis)
    call unit_vector(la(1:3), mu_axis)
    call cross_product3(ep_axis, em_axis, cross)
    triple = dot_product(cross, mu_axis)
    values(4) = max(-1d0, min(1d0, triple))
    values(5) = dble(qassoc)*values(4)

    call boost_to_rest(top, ttw, top_frame)
    call boost_to_rest(atop, ttw, atop_frame)
    call boost_to_rest(la, ttw, mu_frame)
    values(6) = dble(qassoc)*( &
         spatial_opening_cosine(mu_frame, top_frame) - &
         spatial_opening_cosine(mu_frame, atop_frame))
  end subroutine compute_associated_spin_observables


  subroutine fill_associated_spin(offset, values, wgts)
    implicit none
    integer, intent(in) :: offset
    double precision, intent(in) :: values(6), wgts(*)
    integer :: i

    do i = 1, 6
      call HwU_fill(offset + i, values(i), wgts)
    end do
  end subroutine fill_associated_spin


  subroutine fail_analysis(message)
    implicit none
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in comprehensive ttW trilepton analysis: '// &
                     trim(message)
    stop 1
  end subroutine fail_analysis

end module analysis_hwu_pp_ttxw_trilepton_module

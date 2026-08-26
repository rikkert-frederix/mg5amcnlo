module fks_weights_module
  use process_dimensions, only: nexternal, nincoming, maxproc, &
                                nsplitorders, qcd_pos, amp_split_size
  use run_state, only: q2fact, scale, dynamical_scale_choice
  use timing_state, only: t_as, tr_s, tr_pdf, t_plot
  use mint_module, only: nintegrals, virt_wgt_mint, born_wgt_mint
  use fks_metadata, only: fks_i_d, fks_j_d
  use setscales_module, only: set_ren_scale, set_fac_scale
  use split_orders, only: amp_split_pos_to_orders
  use chooser_functions_module, only: set_pdg_impl
  use madfks_plot_module, only: outfun_impl
  use fks_model_state_module, only: g => strong_coupling, external_masses
  use fnlo_process_common, only: nfksprocess, soft_counterevent, &
                                 collinear_counterevent, &
                                 soft_collinear_counterevent, real_event, &
                                 event_bjorken_x, ybst_til_tolab, qes2, &
                                 p_born, &
                                 stored_event_momenta => event_momenta, &
                                 idup, idup_d, &
                                 subproc_pd, subproc_iproc, &
                                 orders_tag_plot, virtual_over_born
  implicit none
  private

  integer, parameter, public :: real_contribution = 1
  integer, parameter, public :: born_contribution = 2
  integer, parameter, public :: integrated_contribution = 3
  integer, parameter, public :: soft_contribution = 4
  integer, parameter, public :: collinear_contribution = 5
  integer, parameter, public :: soft_collinear_contribution = 6
  integer, parameter, public :: virtual_contribution = 14
  integer, parameter, public :: averaged_born_contribution = 15

  public :: add_wgt, include_pdf_and_alphas
  public :: reweight_scale, reweight_pdf, get_wgt_no_nbody
  public :: fill_plots, fill_mint_function

  interface
    double precision function dlum(bjorken_x)
      implicit none
      double precision, intent(in) :: bjorken_x(2)
    end function dlum
  end interface

contains

  logical function pdg_equal(pdg1, pdg2)
! Returns .true. if the lists of PDG codes --'pdg1' and 'pdg2'-- are
! equal.
    implicit none
    integer i, pdg1(nexternal), pdg2(nexternal)
    pdg_equal = .true.
    do i = 1, nexternal
      if (pdg1(i) .ne. pdg2(i)) then
        pdg_equal = .false.
        return
      end if
    end do
  end function pdg_equal

  logical function momenta_equal(p1, p2)
! Returns .true. only if the momenta p1 and p2 are equal. To save time,
! it only checks the 0th and 3rd components (energy and z-direction).
    implicit none
    integer i, j
    double precision p1(0:3, nexternal), p2(0:3, nexternal), vtiny
    parameter(vtiny=1d-8)
    momenta_equal = .true.
    do i = 1, nexternal
      do j = 0, 3, 3
        if (p1(j, i) .eq. 0d0 .or. p2(j, i) .eq. 0d0) then
          if (abs(p1(j, i) - p2(j, i)) .gt. vtiny) then
            momenta_equal = .false.
            return
          end if
        else
          if (abs((p1(j, i) - p2(j, i))/max(abs(p1(j, i)), abs(p2(j, i)))) .gt. vtiny) then
            momenta_equal = .false.
            return
          end if
        end if
      end do
    end do
  end function momenta_equal

  logical function momenta_equal_uborn(p1, p2, jfks1, ifks1, jfks2, ifks2)
! Returns .true. only if the momenta p1 and p2 are equal, but with the
! momenta of i_fks and j_fks summed. To save time, it only checks the
! 0th and 3rd components (energy and z-direction).
    implicit none
    integer i, j, jfks1, ifks1, jfks2, ifks2
    double precision p1(0:3, nexternal), p2(0:3, nexternal), pb1(0:3, nexternal), pb2(0:3, nexternal)
! Fill the underlying Born momenta pb1 and pb2
    do i = 1, nexternal
      do j = 0, 3, 3 ! skip x and y components, since they are not used in
! the 'momenta_equal' function
        if (i .lt. ifks1) then
          pb1(j, i) = p1(j, i)
        elseif (i .eq. ifks1) then
! Sum the i_fks to the j_fks momenta (i_fks is always greater than
! j_fks, so this is fine: it will NOT be overwritten later in the
! do-loop)
          pb1(j, jfks1) = pb1(j, jfks1) + p1(j, i)
          pb1(j, nexternal) = 0d0 ! fill the final one with zero's
        else
          pb1(j, i - 1) = p1(j, i)   ! skip the i_fks momenta
        end if
        if (i .lt. ifks2) then
          pb2(j, i) = p2(j, i)
        elseif (i .eq. ifks2) then
          pb2(j, jfks2) = pb2(j, jfks2) + p2(j, i) ! sum i_fks to j_fks momenta
          pb2(j, nexternal) = 0d0 ! fill the final one with zero's
        else
          pb2(j, i - 1) = p2(j, i)   ! skip the i_fks momenta
        end if
      end do
    end do
! Check if they are equal
    momenta_equal_uborn = momenta_equal(pb1, pb2)
  end function momenta_equal_uborn

  subroutine add_wgt(event_slot, type, wgt1, wgt2, wgt3)
! Adds a contribution to the list in weight_lines. 'type' sets the type
! of the contribution and wgt1..wgt3 are the coefficients multiplying
! the logs. The arguments are:
!     type=1 : real-emission
!     type=2 : Born
!     type=3 : integrated counter terms
!     type=4 : soft counter-term
!     type=5 : collinear counter-term
!     type=6 : soft-collinear counter-term
!     type=14: virtual corrections
!     type=15: virt-trick: average born contribution
!     wgt1 : weight of the contribution not multiplying a scale log
!     wgt2 : coefficient of the weight multiplying the log(/mu_R^2/Q^2/)
!     wgt3 : coefficient of the weight multiplying the log(/mu_F^2/Q^2/)
! This subroutine increments the 'icontr' counter: each new call to this
! function makes sure that it's considered a new contribution. For each
! contribution, we save the
!     The type: itype(icontr)
!     The weights: wgt(1,icontr),wgt(2,icontr) and wgt(3,icontr) for
!         wgt1, wgt2 and wgt3, respectively.
!     The Bjorken x's: bjx(1,icontr), bjx(2,icontr)
!     The Ellis-Sexton scale squared used to compute the weight:
!        scales2(1,icontr)
!     The renormalisation scale squared used to compute the weight:
!        scales2(2,icontr)
!     The factorisation scale squared used to compute the weight:
!       scales2(3,icontr)
!     The value of the strong coupling: g_strong(icontr)
!     The FKS configuration: nFKS(icontr)
!     The boost to go from the momenta in the C.o.M. frame to the
!         laboratory frame: y_bst(icontr)
!     The power of the strong coupling (g_strong) for the current
!       weight: QCDpower(icontr)
!     The momenta: momenta(j,i,icontr). For the Born contribution, the
!        counter-term momenta are used. This is okay for any IR-safe
!        observables.
!     The PDG codes: pdg(i,icontr). Always the ones with length
!        'nexternal' are used, because the momenta are also the
!        'nexternal' ones. This is okay for IR-safe observables.
!     The PDG codes of the underlying Born process:
!        pdg_uborn(i,icontr). The PDGs of j_fks and i_fks are combined
!        to get the PDG code of the mother. The extra parton is given a
!        PDG=21 (gluon) code.
! Not set in this subroutine, but included in the weight_lines module
! are the
!     wgts(iwgt,icontr) : weights including scale/PDFs/logs. These are
!        normalised so that they can be used directly to compute cross
!        sections and fill plots. 'iwgt' goes from 1 to the maximum
!        number of weights obtained from scale and PDF reweighting, with
!        the iwgt=1 element being the central value.
!     plot_id(icontr) : =20 for Born, 11 for real-emission and 12 for
!        anything else.
!     plot_wgts(iwgt,icontr) : same as wgts(), but only non-zero for
!        unique contributions and non-unique are added to the unique
!        ones. 'Unique' here is defined that they would be identical in
!        an analysis routine (i.e. same momenta and PDG codes)
!     niproc(icontr) : number of combined subprocesses in parton_lum_*.f
!     parton_iproc(iproc,icontr) : value of the PDF for the iproc
!        contribution
!     parton_pdg(nexternal,iproc,icontr) : value of the PDG codes for
!     the iproc contribution
!     ipr(icontr): for separate_flavour_configs: the iproc of current
!        contribution
    use weight_lines
    use extra_weights
    use FKSParams
    implicit none
    integer event_slot, type, i, j
    logical foundIt
    double precision wgt1, wgt2, wgt3

    if (wgt1 .eq. 0d0 .and. wgt2 .eq. 0d0 .and. wgt3 .eq. 0d0) return
! Check for NaN's and INF's. Simply skip the contribution
    if (wgt1 .ne. wgt1) return
    if (wgt2 .ne. wgt2) return
    if (wgt3 .ne. wgt3) return

! Apply user-defined (in FKS_params.dat) contribution type filters if necessary

    if (SelectedContributionTypes(0) .gt. 0) then
      foundIt = .false.
      do i = 1, SelectedContributionTypes(0)
        if (type .eq. SelectedContributionTypes(i)) then
          foundIt = .true.
          exit
        end if
      end do
      if (.not. foundIt) then
! This contribution was not part of the user selection. Skip it.
        return
      end if
    end if

    icontr = icontr + 1
    call weight_lines_allocated(nexternal, icontr, max_wgt, max_iproc)
    itype(icontr) = type

    wgt(1, icontr) = wgt1
    wgt(2, icontr) = wgt2
    wgt(3, icontr) = wgt3

    bjx(1, icontr) = event_bjorken_x(1, event_slot)
    bjx(2, icontr) = event_bjorken_x(2, event_slot)
    scales2(1, icontr) = QES2
    scales2(2, icontr) = scale**2
    scales2(3, icontr) = q2fact(1)
    g_strong(icontr) = g
    nFKS(icontr) = nFKSprocess
    y_bst(icontr) = ybst_til_tolab(event_slot)
    qcdpower(icontr) = QCD_power
    orderstag(icontr) = orders_tag
    amppos(icontr) = amp_pos
    ipr(icontr) = 0
    call set_pdg_impl(icontr, nFKSprocess, idup)

    if (type .eq. real_contribution) then
! Real-emission contribution with n+1-body kinematics.
      do i = 1, nexternal
        do j = 0, 3
          momenta(j, i, icontr) = stored_event_momenta(j, i, real_event)
        end do
      end do
    elseif (type .ge. born_contribution .and. type .le. soft_collinear_contribution .or. &
            type .eq. virtual_contribution .or. &
            type .eq. averaged_born_contribution .or. &
            (type .ge. 20 .and. type .le. 22)) then
! Born, counter term, soft-virtual, or n-body real contributions.
      do i = 1, nexternal
        do j = 0, 3
          if (stored_event_momenta(0, 1, soft_counterevent) .gt. 0d0 &
              .and. type .ne. 5) then
            momenta(j, i, icontr) = &
              stored_event_momenta(j, i, soft_counterevent)
          elseif (stored_event_momenta(0, 1, collinear_counterevent) &
                  .gt. 0d0) then
            momenta(j, i, icontr) = &
              stored_event_momenta(j, i, collinear_counterevent)
          elseif (stored_event_momenta(0, 1, soft_collinear_counterevent) &
                  .gt. 0d0) then
            momenta(j, i, icontr) = &
              stored_event_momenta(j, i, soft_collinear_counterevent)
          elseif (i .lt. fks_i_d(nFKSprocess)) then
            momenta(j, i, icontr) = p_born(j, i)
          elseif (i .eq. fks_i_d(nFKSprocess)) then
            momenta(j, i, icontr) = 0d0
          else
            momenta(j, i, icontr) = p_born(j, i - 1)
          end if
        end do
      end do
    else
      write (*, *) 'ERROR: unknown type in add_wgt', type
      stop 1
    end if
    return
  end subroutine add_wgt
  subroutine include_PDF_and_alphas
! Multiply the saved wgt() info by the PDFs, alpha_S and the scale
! dependence and saves the weights in the wgts() array. The weights in
! this array are now correctly normalised to compute the cross section
! or to fill histograms.
    use weight_lines
    use extra_weights
    use mint_module
    use FKSParams
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer i, j, iamp, icontr_orig
    logical virt_found
    double precision xlum, mu2_r, mu2_f, mu2_q, wgt_wo_pdf, conv
    parameter(conv=389379660d0) ! conversion to picobarns
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
    virt_found = .false.
! number of contributions before they are (possibly) increased through a
! call to separate_flavour_config().
    icontr_orig = icontr
    i = 0
    do while (i .lt. icontr)
      i = i + 1
      nFKSprocess = nFKS(i)
      mu2_q = scales2(1, i)
      mu2_r = scales2(2, i)
      mu2_f = scales2(3, i)
      q2fact(1) = mu2_f
      q2fact(2) = mu2_f
! call the PDFs
      xlum = dlum(bjx(:, i))
! iwgt=1 is the central value (i.e. no scale/PDF reweighting).
      iwgt = 1
      call weight_lines_allocated(nexternal, max_contr, iwgt, subproc_iproc)
! set_pdg_codes fills the niproc, parton_iproc, parton_pdg and
! parton_pdg_uborn [Do only for the contributions that were already
! available as part of the input -- NOT the ones that are created
! through the call to separate_flavour_config(), since that will
! overwrite the relevant information.]
      if (i .le. icontr_orig) call set_pdg_codes(subproc_iproc, subproc_pd, nFKSprocess, i)
      if (separate_flavour_configs .and. ipr(i) .eq. 0) then
        call separate_flavour_config(i) ! this increases icontr
      end if
      if (separate_flavour_configs .and. ipr(i) .ne. 0) then
        if (nincoming .eq. 2) then
          xlum = subproc_pd(ipr(i))*conv
        else
          xlum = subproc_pd(ipr(i))
        end if
      end if
      wgt_wo_pdf = (wgt(1, i) + wgt(2, i)*log(mu2_r/mu2_q) &
                    + wgt(3, i)*log(mu2_f/mu2_q))*g_strong(i)**QCDpower(i)
      wgts(iwgt, i) = xlum*wgt_wo_pdf
      do j = 1, subproc_iproc
        parton_iproc(j, i) = parton_iproc(j, i)*wgt_wo_pdf
      end do
      if (itype(i) .eq. virtual_contribution .and. .not. virt_found) then
        virt_found = .true.
! Special for the soft-virtual needed for the virt-tricks. The
! *_wgt_mint variable should be directly passed to the mint-integrator
! and not be part of the plots nor computation of the cross section.
        virt_wgt_mint(0) = virt_wgt_mint(0)*xlum
        born_wgt_mint(0) = born_wgt_mint(0)*xlum
        do iamp = 1, amp_split_size
          call amp_split_pos_to_orders(iamp, orders)
          QCD_power = orders(qcd_pos)
          virt_wgt_mint(iamp) = virt_wgt_mint(iamp)*xlum
          born_wgt_mint(iamp) = born_wgt_mint(iamp)*xlum
        end do
      end if
    end do
    call cpu_time(tAfter)
    t_as = t_as + (tAfter - tBefore)
    return
  end subroutine include_PDF_and_alphas

  subroutine separate_flavour_config(ict)
    use weight_lines
    implicit none
    integer ict, i_add, i, j, k, ict_new
    if (niproc(ict) .eq. 1) return
    i_add = niproc(ict) - 1
    call weight_lines_allocated(nexternal, icontr + i_add, max_wgt, max_iproc)
    do i = 1, niproc(ict)
      if (i .eq. 1) then
        niproc(ict) = 1
        ipr(ict) = 1
        cycle
      end if
      ict_new = icontr + (i - 1)
      ipr(ict_new) = i
      itype(ict_new) = itype(ict)
      do j = 1, 3
        wgt(j, ict_new) = wgt(j, ict)
        scales2(j, ict_new) = scales2(j, ict)
      end do
      bjx(:, ict_new) = bjx(:, ict)
      g_strong(ict_new) = g_strong(ict)
      nFKS(ict_new) = nFKS(ict)
      y_bst(ict_new) = y_bst(ict)
      QCDpower(ict_new) = QCDpower(ict)
      orderstag(ict_new) = orderstag(ict)
      do k = 1, nexternal
        do j = 0, 3
          momenta(j, k, ict_new) = momenta(j, k, ict)
        end do
        pdg(k, ict_new) = parton_pdg(k, i, ict)
        pdg_uborn(k, ict_new) = parton_pdg_uborn(k, i, ict)
        parton_pdg(k, 1, ict_new) = parton_pdg(k, i, ict)
        parton_pdg_uborn(k, 1, ict_new) = parton_pdg_uborn(k, i, ict)
      end do
      niproc(ict_new) = 1
      parton_iproc(1, ict_new) = parton_iproc(i, ict)
    end do
    icontr = icontr + i_add
    return
  end subroutine separate_flavour_config

  subroutine set_pdg_codes(iproc, pd, iFKS, ict)
    use weight_lines
    implicit none
    integer j, k, iproc, ict, iFKS
    double precision pd(0:maxproc), conv
    parameter(conv=389379660d0) ! conversion to picobarns
! save also the separate contributions to the PDFs and the corresponding
! PDG codes
    niproc(ict) = iproc
    do j = 1, iproc
      if (nincoming .eq. 2) then
        parton_iproc(j, ict) = pd(j)*conv
      else
!           Keep GeV's for decay processes (no conv. factor needed)
        parton_iproc(j, ict) = pd(j)
      end if
      do k = 1, nexternal
        parton_pdg(k, j, ict) = idup_d(iFKS, k, j)
        if (k .lt. fks_j_d(iFKS)) then
          parton_pdg_uborn(k, j, ict) = idup_d(iFKS, k, j)
        elseif (k .eq. fks_j_d(iFKS)) then
          if (abs(idup_d(iFKS, fks_i_d(iFKS), j)) .eq. &
              abs(idup_d(iFKS, fks_j_d(iFKS), j)) .and. &
              abs(pdg(fks_i_d(iFKS), ict)) .ne. 21) then
            parton_pdg_uborn(k, j, ict) = 21
          elseif (abs(idup_d(iFKS, fks_i_d(iFKS), j)) .eq. 21) then
            parton_pdg_uborn(k, j, ict) = idup_d(iFKS, fks_j_d(iFKS), j)
          elseif (idup_d(iFKS, fks_j_d(iFKS), j) .eq. 21) then
            parton_pdg_uborn(k, j, ict) = -idup_d(iFKS, fks_i_d(iFKS), j)
          else
            write (*, *) 'set_pdg_codes ', 'ERROR#3 in PDG assigment for underlying Born'
            stop 1
          end if
        elseif (k .lt. fks_i_d(iFKS)) then
          parton_pdg_uborn(k, j, ict) = idup_d(iFKS, k, j)
        elseif (k .eq. nexternal) then
          parton_pdg_uborn(k, j, ict) = 21
        elseif (k .ge. fks_i_d(iFKS)) then
          parton_pdg_uborn(k, j, ict) = idup_d(iFKS, k + 1, j)
        end if
      end do
    end do
    return
  end subroutine set_pdg_codes

  subroutine reweight_scale
! Use the saved weight_lines info to perform scale reweighting. Extends the
! wgts() array to include the weights.
    use weight_lines
    use extra_weights
    use FKSParams
    use alfas_functions_module, only: alphas
    implicit none
    real :: tBefore, tAfter
    integer i, kr, kf, iwgt_save, dd
    double precision xlum(maxscales), pi, mu2_r(maxscales), c_mu2_r, c_mu2_f, mu2_f(maxscales), mu2_q, g(maxscales), conv
    parameter(pi=3.1415926535897932385d0)
    parameter(conv=389379660d0) ! conversion to picobarns
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
! currently we have 'iwgt' weights in the wgts() array.
    iwgt_save = iwgt
! loop over all the contributions in the weight lines module
    do i = 1, icontr
      iwgt = iwgt_save
      nFKSprocess = nFKS(i)
      mu2_q = scales2(1, i)
! Loop over the dynamical_scale_choices
      do dd = 1, dyn_scale(0)
! renormalisation scale variation (requires recomputation of the strong
! coupling)
        call set_mu_central(i, dd, c_mu2_r, c_mu2_f)
        do kr = 1, nint(scalevarR(0))
          if ((.not. lscalevar(dd)) .and. kr .ne. 1) exit
          mu2_r(kr) = c_mu2_r*scalevarR(kr)**2
          g(kr) = sqrt(4d0*pi*alphas(sqrt(mu2_r(kr))))
        end do
! factorisation scale variation (require recomputation of the PDFs)
        do kf = 1, nint(scalevarF(0))
          if ((.not. lscalevar(dd)) .and. kf .ne. 1) exit
          mu2_f(kf) = c_mu2_f*scalevarF(kf)**2
          q2fact(1) = mu2_f(kf)
          q2fact(2) = mu2_f(kf)
          xlum(kf) = dlum(bjx(:, i))
          if (separate_flavour_configs .and. ipr(i) .ne. 0) then
            if (nincoming .eq. 2) then
              xlum(kf) = subproc_pd(ipr(i))*conv
            else
              xlum(kf) = subproc_pd(ipr(i))
            end if
          end if
        end do
        do kf = 1, nint(scalevarF(0))
          if ((.not. lscalevar(dd)) .and. kf .ne. 1) exit
          do kr = 1, nint(scalevarR(0))
            if ((.not. lscalevar(dd)) .and. kr .ne. 1) exit
            iwgt = iwgt + 1   ! increment the iwgt for the wgts() array
            call weight_lines_allocated(nexternal, max_contr, iwgt, max_iproc)
! add the weights to the array
            wgts(iwgt, i) = xlum(kf) &
                            *(wgt(1, i) + wgt(2, i)*log(mu2_r(kr)/mu2_q) &
                              + wgt(3, i)*log(mu2_f(kf)/mu2_q)) &
                            *g(kr)**QCDpower(i)
          end do
        end do
      end do
    end do
    call cpu_time(tAfter)
    tr_s = tr_s + (tAfter - tBefore)
    return
  end subroutine reweight_scale

  subroutine reweight_pdf
! Use the saved weight_lines info to perform PDF reweighting. Extends the
! wgts() array to include the weights.
    use weight_lines
    use extra_weights
    use FKSParams
    use alfas_functions_module, only: alphas
    implicit none
    real :: tBefore, tAfter
    integer n, i, nn
    double precision xlum, pi, mu2_r, mu2_f, mu2_q, g, conv
    parameter(pi=3.1415926535897932385d0)
    parameter(conv=389379660d0) ! conversion to picobarns
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
    do nn = 1, lhaPDFid(0)
! Use as external loop the one over the PDF sets and as internal the one
! over the icontr. This reduces the number of calls to InitPDF and
! allows for better caching of the PDFs
      do n = 0, nmemPDF(nn)
        iwgt = iwgt + 1
        call weight_lines_allocated(nexternal, max_contr, iwgt, max_iproc)
        call InitPDFm(nn, n)
        do i = 1, icontr
          nFKSprocess = nFKS(i)
          mu2_q = scales2(1, i)
          mu2_r = scales2(2, i)
          mu2_f = scales2(3, i)
          q2fact(1) = mu2_f
          q2fact(2) = mu2_f
! Compute the luminosity
          xlum = dlum(bjx(:, i))
          if (separate_flavour_configs .and. ipr(i) .ne. 0) then
            if (nincoming .eq. 2) then
              xlum = subproc_pd(ipr(i))*conv
            else
              xlum = subproc_pd(ipr(i))
            end if
          end if
! Recompute the strong coupling: alpha_s in the PDF might change
          g = sqrt(4d0*pi*alphas(sqrt(mu2_r)))
! add the weights to the array
          wgts(iwgt, i) = xlum*(wgt(1, i) + wgt(2, i)*log(mu2_r/mu2_q) + wgt(3, i)*log(mu2_f/mu2_q))*g**QCDpower(i)
        end do
      end do
    end do
    call InitPDFm(1, 0)
    call cpu_time(tAfter)
    tr_pdf = tr_pdf + (tAfter - tBefore)
    return
  end subroutine reweight_pdf

  subroutine get_wgt_no_nbody(sig)
! Sums all the central weights that contribution to the cross section
! excluding the nbody contributions.
    use weight_lines
    implicit none
    double precision sig
    integer i
    sig = 0d0
    if (icontr .eq. 0) return
    do i = 1, icontr
      if (itype(i) .ne. born_contribution .and. &
          itype(i) .ne. integrated_contribution .and. &
          itype(i) .ne. virtual_contribution .and. &
          itype(i) .ne. averaged_born_contribution) then
        sig = sig + wgts(1, i)
      end if
    end do
    return
  end subroutine get_wgt_no_nbody

  subroutine fill_plots
! Calls the analysis routine (which fill plots) for all the
! contributions in the weight_lines module. Instead of really calling
! it for all, it first checks if weights can be summed (i.e. they have
! the same PDG codes and the same momenta) before calling the analysis
! to greatly reduce the calls to the analysis routines.
    use weight_lines
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer i, ii, j, max_weight
    double precision, allocatable :: www(:)
! stuff for plotting the different splitorders
    save max_weight
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
! fill the plots_wgts. Check if we can sum weights together before
! calling the analysis routines. This is the case if the PDG codes and
! the momenta are identical.
    do i = 1, icontr
      do j = 1, iwgt
        plot_wgts(j, i) = 0d0
      end do
! The following if lines have been changed with respect to the
! usual (with just 3 plot ids: 20, 11 and 12):
!  The kinematics of soft and collinear counterterms may
!  be different, for those processes without soft singularities
!  from initial(final)-state configurations when the
!  final(initial) confs are integrated (e.g. a a > e+ e-)
!  This gives no problem for normal histogramming (and in
!  fact plot_id 11 13 and 14 are merged into ibody=2 in
!  outfun, but separate counterterm types are useful for diagnostics.
!  Note that the separation between soft and soft-virtual
!  may not be needed in reality
      if (itype(i) .eq. born_contribution) then
        plot_id(i) = 20 ! Born
      elseif (itype(i) .eq. real_contribution) then
        plot_id(i) = 11 ! real-emission
      elseif (itype(i) .eq. collinear_contribution) then
        plot_id(i) = 13 ! collinear counter term
      elseif (itype(i) .eq. soft_collinear_contribution) then
        plot_id(i) = 14 ! soft collinear counter term
      else
        plot_id(i) = 12 ! soft-virtual and soft counter term
      end if
! Loop over all previous icontr. If the plot_id, PDGs and momenta are
! equal to a previous icountr, add the current weight to the plot_wgts
! of that contribution and exit the do-loop. This loop extends to 'i',
! so if the current weight cannot be summed to a previous one, the ii=i
! contribution makes sure that it is added as a new element.
      do ii = 1, i
        if (orderstag(ii) .ne. orderstag(i)) cycle
        if (plot_id(ii) .ne. plot_id(i)) cycle
        if (plot_id(i) .ne. 11) then
          if (.not. pdg_equal(pdg_uborn(1, ii), pdg_uborn(1, i))) cycle
        else
          if (.not. pdg_equal(pdg(1, ii), pdg(1, i))) cycle
        end if
        if (plot_id(i) .ne. 11) then
          if (.not. momenta_equal_uborn( &
              momenta(0, 1, ii), momenta(0, 1, i), &
              fks_j_d(nFKS(ii)), fks_i_d(nFKS(ii)), &
              fks_j_d(nFKS(i)), fks_i_d(nFKS(i)))) cycle
        else
          if (.not. momenta_equal(momenta(0, 1, ii), momenta(0, 1, i))) cycle
        end if
        do j = 1, iwgt
          plot_wgts(j, ii) = plot_wgts(j, ii) + wgts(j, i)
        end do
        exit
      end do
    end do
    do i = 1, icontr
      if (plot_wgts(1, i) .ne. 0d0) then
        if (.not. allocated(www)) then
          allocate (www(iwgt))
          max_weight = iwgt
        elseif (iwgt .ne. max_weight) then
          write (*, *) 'Error in fill_plots (fks_singular.f): '// &
            'number of weights should not vary between PS points', &
            iwgt, max_weight
          stop
        end if
        do j = 1, iwgt
          www(j) = plot_wgts(j, i)
        end do
! call the analysis/histogramming routines
        orders_tag_plot = orderstag(i)
        call outfun_impl(momenta(0, 1, i), y_bst(i), www, pdg(1, i), plot_id(i), &
                         external_masses)
      end if
    end do
    call cpu_time(tAfter)
    t_plot = t_plot + (tAfter - tBefore)
    return
  end subroutine fill_plots

  subroutine fill_mint_function(f)
! Fills the function that is returned to the MINT integrator
    use weight_lines
    use mint_module
    implicit none
    integer i, iamp, ithree, isix
    double precision f(nintegrals), sigint
    sigint = 0d0
    do i = 1, icontr
      sigint = sigint + wgts(1, i)
    end do
    f(1) = abs(sigint)
    f(2) = sigint
    f(4) = virtual_over_born    ! not used for anything
    do iamp = 0, amp_split_size
      if (iamp .eq. 0) then
        f(3) = 0d0
        f(6) = 0d0
        f(5) = 0d0
        do i = 1, amp_split_size
          f(3) = f(3) + virt_wgt_mint(i)
          f(6) = f(6) + born_wgt_mint(i)
        end do
        f(5) = abs(f(3))!v3.5.4, this fixes a wrong behaviour
      else
        ithree = 2*iamp + 5
        isix = 2*iamp + 6
        f(ithree) = virt_wgt_mint(iamp)
        f(isix) = born_wgt_mint(iamp)
      end if
    end do
    return
  end subroutine fill_mint_function



  subroutine set_mu_central(ic, dd, c_mu2_r, c_mu2_f)
    use weight_lines
    use extra_weights
    implicit none
    integer ic, dd, i, j
    double precision c_mu2_r, c_mu2_f, muR, muF(2), pp(0:3, nexternal)
    if (dd .eq. 1) then
      c_mu2_r = scales2(2, ic)
      c_mu2_f = scales2(3, ic)
    else
! need to recompute the scales using the momenta
      dynamical_scale_choice = dyn_scale(dd)
      do i = 1, nexternal
        do j = 0, 3
          pp(j, i) = momenta(j, i, ic)
        end do
      end do
      call set_ren_scale(pp, muR)
      c_mu2_r = muR**2
      call set_fac_scale(pp, muF)
      c_mu2_f = muF(1)**2
!     reset the default dynamical_scale_choice
      dynamical_scale_choice = dyn_scale(1)
    end if
    return
  end subroutine set_mu_central


end module fks_weights_module

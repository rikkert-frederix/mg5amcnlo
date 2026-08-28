module fks_contributions_module
  use process_dimensions, only: nexternal, nsplitorders, qcd_pos, &
                                amp_split_size
  use timing_state, only: tBorn, tIS, tReal, tCount, tf_nb, tf_all
  use mint_module, only: virt_wgt_mint, born_wgt_mint
  use split_orders, only: get_orders_tag, amp_split_pos_to_orders
  use madfks_plot_module, only: initplot_impl
  use fks_model_state_module, only: g => strong_coupling, external_masses
  use decay_chain_metadata, only: has_decay_chains
  use decay_chain_kinematics, only: fks_leg_mass
  use nlo_decay_metadata, only: has_nlo_decay
  use nlo_decay_kinematics, only: nlo_decay_fks_sister_mass
  use nlo_contribution_bundle, only: active_contribution_has_virtual, &
       active_virtual_grid_index
  use fks_singular_module, only: evaluate_fks_sij, sreal, sreal_deg, &
                                 bornsoftvirtual, fks_subtraction_shat
  use fks_weights_module, only: add_wgt, real_contribution, &
                                born_contribution, integrated_contribution, &
                                soft_contribution, collinear_contribution, &
                                soft_collinear_contribution, &
                                virtual_contribution, &
                                averaged_born_contribution
  use fnlo_process_common, only: soft_counterevent, &
                                 collinear_counterevent, &
                                 soft_collinear_counterevent, real_event, &
                                 event_xi, event_y, event_xi_hat, &
                                 event_xi_max, event_xi_norm, &
                                 event_shat, &
                                 stored_event_jacobian => event_jacobian, &
                                 stored_event_momenta => event_momenta, &
                                 p_born, nfksprocess, i_fks, j_fks, &
                                 f_b, f_nb, f_r, f_s, f_c, f_dc, f_sc, &
                                 f_dsc, &
                                 xiscut_used, xibsvcut_used, delta_used, &
                                 xicut_used, fkssymmetryfactor, &
                                 fkssymmetryfactorborn, &
                                 fkssymmetryfactordeg, nocntevents, &
                                 this_config, diagramsymmetryfactor, &
                                 config_map, amp2, amp_split, amp_split_virt, &
                                 amp_split_born_for_virt, amp_split_avv, &
                                 amp_split_wgtnstmp, &
                                 amp_split_wgtwnstmpmuf, &
                                 amp_split_wgtwnstmpmur, &
                                 amp_split_wgtdegrem_xi, &
                                 amp_split_wgtdegrem_lxi, &
                                 amp_split_wgtdegrem_muf
  implicit none
  private

  double precision, parameter :: deltas = 1d0

  public :: compute_born, compute_nbody_noborn, compute_real_emission
  public :: compute_soft_counter_term, compute_collinear_counter_term
  public :: compute_soft_collinear_ct_impl
  public :: compute_prefactors_nbody, compute_prefactors_n1body
  public :: include_multichannel_enhance

contains

  logical function fks_sister_is_massless()
    if (has_nlo_decay()) then
      fks_sister_is_massless = &
           nlo_decay_fks_sister_mass(nfksprocess) == 0d0
    else if (has_decay_chains()) then
      fks_sister_is_massless = &
           fks_leg_mass(nfksprocess, j_fks) == 0d0
    else
      fks_sister_is_massless = external_masses(j_fks) == 0d0
    end if
  end function fks_sister_is_massless

  subroutine compute_born
! This subroutine computes the Born matrix elements and adds its value
! to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp

    double precision wgt_c
    double precision wgt1
    call cpu_time(tBefore)
    if (f_b .eq. 0d0) return
    if (event_xi_hat(real_event)*event_xi_max(soft_counterevent) .gt. &
        xiBSVcut_used) return
    call sborn(p_born, wgt_c)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*f_b/g**(qcd_power)
      call add_wgt(soft_counterevent, born_contribution, wgt1, 0d0, 0d0)
    end do

    call cpu_time(tAfter)
    tBorn = tBorn + (tAfter - tBefore)
    return
  end subroutine compute_born

  subroutine compute_nbody_noborn
! This subroutine computes the soft-virtual matrix elements and adds its
! value to the list of weights using the add_wgt subroutine
    use extra_weights
    use mint_module
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp, virtual_grid

    double precision wgt1, wgt2, wgt3, bsv_wgt, virt_wgt, born_wgt, g22, wgt4
    call cpu_time(tBefore)
    if (f_nb .eq. 0d0) return
    if (event_xi_hat(real_event)*event_xi_max(soft_counterevent) .gt. &
        xiBSVcut_used) return
    call bornsoftvirtual(soft_counterevent, &
                         stored_event_momenta(:, :, soft_counterevent), &
                         bsv_wgt, virt_wgt, born_wgt)
    do iamp = 1, amp_split_size
      if (amp_split_wgtnstmp(iamp) .eq. 0d0 .and. &
          amp_split_wgtwnstmpmur(iamp) .eq. 0d0 .and. &
          amp_split_wgtwnstmpmuf(iamp) .eq. 0d0 .and. &
          amp_split_avv(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = amp_split_wgtnstmp(iamp)*f_nb/g22
      wgt2 = amp_split_wgtwnstmpmur(iamp)*f_nb/g22
      wgt3 = amp_split_wgtwnstmpmuf(iamp)*f_nb/g22
      wgt4 = amp_split_avv(iamp)*f_nb/g22
      call add_wgt(soft_counterevent, integrated_contribution, &
                   wgt1, wgt2, wgt3)
      call add_wgt(soft_counterevent, averaged_born_contribution, &
                   wgt4, 0d0, 0d0)
    end do
! Special for the soft-virtual needed for the virt-tricks. The
! *_wgt_mint variable should be directly passed to the mint-integrator
! and not be part of the plots nor computation of the cross section.
    if (active_contribution_has_virtual()) then
      virt_wgt_mint(0) = virt_wgt_mint(0) + virt_wgt*f_nb
      born_wgt_mint(0) = born_wgt_mint(0) + born_wgt*f_b
    end if
    do iamp = 1, amp_split_size
      if (.not. active_contribution_has_virtual()) exit
      if (amp_split_virt(iamp) .eq. 0d0) cycle
      virtual_grid = active_virtual_grid_index(iamp, amp_split_size)
      if (virtual_grid == 0) then
        write (*,*) 'ERROR: a virtual weight has no bundle grid'
        stop 1
      end if
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split_virt(iamp)*f_nb
      virt_wgt_mint(virtual_grid) = &
           virt_wgt_mint(virtual_grid) + wgt1
      born_wgt_mint(virtual_grid) = born_wgt_mint(virtual_grid) + &
           amp_split_born_for_virt(iamp)*f_nb
      wgt1 = wgt1/g**(QCD_power)
      call add_wgt(soft_counterevent, virtual_contribution, wgt1, 0d0, 0d0)
    end do

    call cpu_time(tAfter)
    tIS = tIS + (tAfter - tBefore)
    return
  end subroutine compute_nbody_noborn

  subroutine compute_real_emission(p)
! This subroutine computes the real-emission matrix elements and adds
! its value to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
    double precision s_ev, p(0:3, nexternal), wgt1, fx_ev
    call cpu_time(tBefore)
    if (f_r .eq. 0d0) return
    s_ev = evaluate_fks_sij(real_event, p, i_fks, j_fks, &
                            event_xi(real_event), event_y(real_event))
    if (s_ev .le. 0.d0) return
    call sreal(real_event, p, event_xi(real_event), &
               event_y(real_event), fx_ev)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*s_ev*f_r/g**(qcd_power)
      call add_wgt(real_event, real_contribution, wgt1, 0d0, 0d0)
    end do
    call cpu_time(tAfter)
    tReal = tReal + (tAfter - tBefore)
    return
  end subroutine compute_real_emission

  subroutine compute_soft_counter_term
! This subroutine computes the soft counter term and adds its value to
! the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
    double precision wgt1, s_s, fx_s, zero, g22
    parameter(zero=0d0)
    call cpu_time(tBefore)
    if (f_s .eq. 0d0) return
    if (event_xi_hat(real_event)*event_xi_max(soft_counterevent) .gt. &
        xiScut_used) return
    s_s = evaluate_fks_sij(soft_counterevent, &
            stored_event_momenta(:, :, soft_counterevent), &
            i_fks, j_fks, zero, event_y(real_event))
    if (s_s .le. 0d0) return
    call sreal(soft_counterevent, &
               stored_event_momenta(:, :, soft_counterevent), &
               0d0, event_y(real_event), fx_s)

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = 0d0
      if (event_xi(real_event) .le. xiScut_used) then
        wgt1 = -amp_split(iamp)*s_s*f_s/g22
      end if
      if (wgt1 .ne. 0d0) &
        call add_wgt(soft_counterevent, soft_contribution, wgt1, 0d0, 0d0)
    end do

    call cpu_time(tAfter)
    tCount = tCount + (tAfter - tBefore)
    return
  end subroutine compute_soft_counter_term

  subroutine compute_collinear_counter_term
! This subroutine computes the collinear counter term and adds its value
! to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
    double precision one, s_c, fx_c, deg_xi_c, deg_lxi_c, wgt1, wgt3, g22
    parameter(one=1d0)
    call cpu_time(tBefore)
    if (f_c .eq. 0d0 .and. f_dc .eq. 0d0) return
    if (event_y(real_event) .le. 1d0 - deltaS .or. &
        .not. fks_sister_is_massless()) return
    s_c = evaluate_fks_sij(collinear_counterevent, &
            stored_event_momenta(:, :, collinear_counterevent), &
            i_fks, j_fks, event_xi(collinear_counterevent), one)
    if (s_c .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg(collinear_counterevent, &
                   stored_event_momenta(:, :, collinear_counterevent), &
                   event_xi(collinear_counterevent), deg_xi_c, deg_lxi_c)
    call sreal(collinear_counterevent, &
               stored_event_momenta(:, :, collinear_counterevent), &
               event_xi(collinear_counterevent), one, fx_c)

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_xi(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_lxi(iamp) .eq. 0d0) cycle

      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = -amp_split(iamp)*s_c*f_c/g22
      wgt1 = wgt1 + (amp_split_wgtdegrem_xi(iamp) + &
                     amp_split_wgtdegrem_lxi(iamp)* &
                     log(event_xi(collinear_counterevent)))*f_dc/g22
      wgt3 = amp_split_wgtdegrem_muF(iamp)*f_dc/g22
      if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) &
        call add_wgt(collinear_counterevent, collinear_contribution, &
                     wgt1, 0d0, wgt3)
    end do

    call cpu_time(tAfter)
    tCount = tCount + (tAfter - tBefore)
    return
  end subroutine compute_collinear_counter_term

  subroutine compute_soft_collinear_ct_impl
! This subroutine computes the soft-collinear counter term and adds its
! value to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
    double precision zero, one, s_sc, fx_sc, wgt1, wgt3, deg_xi_sc, deg_lxi_sc, g22
    parameter(zero=0d0, one=1d0)
    call cpu_time(tBefore)
    if (f_sc .eq. 0d0 .and. f_dsc(1) .eq. 0d0 .and. f_dsc(2) .eq. 0d0 .and. f_dsc(3) .eq. 0d0 .and. f_dsc(4) .eq. 0d0) return
    if (event_xi_hat(real_event)*event_xi_max(collinear_counterevent) &
        .ge. xiScut_used .or. event_y(real_event) .le. 1d0 - deltaS &
        .or. .not. fks_sister_is_massless()) return
    s_sc = evaluate_fks_sij(soft_collinear_counterevent, &
             stored_event_momenta(:, :, soft_collinear_counterevent), &
             i_fks, j_fks, zero, one)
    if (s_sc .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg(soft_collinear_counterevent, &
      stored_event_momenta(:, :, soft_collinear_counterevent), &
      zero, deg_xi_sc, deg_lxi_sc)
    call sreal(soft_collinear_counterevent, &
               stored_event_momenta(:, :, soft_collinear_counterevent), &
               zero, one, fx_sc)

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_xi(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_lxi(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = 0d0
      wgt3 = 0d0
      if (event_xi(collinear_counterevent) .lt. xiScut_used) then
        wgt1 = amp_split(iamp)*s_sc*f_sc/g22
        wgt1 = wgt1 + ( &
               -(amp_split_wgtdegrem_xi(iamp) + &
                 amp_split_wgtdegrem_lxi(iamp)* &
                 log(event_xi(collinear_counterevent)))*f_dsc(1) &
               - (amp_split_wgtdegrem_xi(iamp)*f_dsc(2) + &
                  amp_split_wgtdegrem_lxi(iamp)*f_dsc(3)))/g22
        wgt3 = -amp_split_wgtdegrem_muF(iamp)*f_dsc(4)/g22
      end if
      if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) &
        call add_wgt(soft_collinear_counterevent, &
                     soft_collinear_contribution, wgt1, 0d0, wgt3)
    end do

    call cpu_time(tAfter)
    tCount = tCount + (tAfter - tBefore)
    return
  end subroutine compute_soft_collinear_ct_impl

  subroutine compute_prefactors_nbody(vegas_wgt)
! Compute all the relevant prefactors for the Born and the soft-virtual,
! i.e. all the nbody contributions. Also initialises the plots and
! bpower.
    use extra_weights
    use mint_module
    implicit none
    real :: tBefore, tAfter
    double precision pi, vegas_wgt, subtraction_shat
    logical firsttime
    data firsttime/.true./
    parameter(pi=3.1415926535897932385d0)
!
    call cpu_time(tBefore)
    if (firsttime) then
      call initplot_impl()
      firsttime = .false.
    end if
! f_* multiplication factors for Born and nbody
    subtraction_shat = fks_subtraction_shat(soft_counterevent)
    f_b = stored_event_jacobian(soft_counterevent)* &
          event_xi_norm(real_event)/ &
          (min(event_xi_max(real_event), xiBSVcut_used)* &
           subtraction_shat/(16*pi**2))* &
          fkssymmetryfactorBorn*vegas_wgt
    f_nb = f_b
    call cpu_time(tAfter)
    tf_nb = tf_nb + (tAfter - tBefore)
    return
  end subroutine compute_prefactors_nbody

  subroutine include_multichannel_enhance(imode)
    implicit none
    real :: tBefore, tAfter
    double precision xnoborn_cnt, xtot, wgt_c, enhance
    data xnoborn_cnt/0d0/
    integer inoborn_cnt, i, imode
    data inoborn_cnt/0/

    call cpu_time(tBefore)

! Compute the multi-channel enhancement factor 'enhance'.
    enhance = 1.d0
    if (p_born(0, 1) .gt. 0d0) then
      call sborn(p_born, wgt_c)
    elseif (p_born(0, 1) .lt. 0d0) then
      enhance = 0d0
    end if
    if (enhance .eq. 0d0) then
      xnoborn_cnt = xnoborn_cnt + 1.d0
      if (log10(xnoborn_cnt) .gt. inoborn_cnt) then
        write (*, *) 'WARNING: no Born momenta more than 10**', inoborn_cnt, 'times'
        inoborn_cnt = inoborn_cnt + 1
      end if
    else
      xtot = 0d0
      if (config_map(0, 0) .eq. 0) then
        write (*, *) 'Fatal error in compute_prefactor_nbody:'//' no Born diagrams ', config_map, '. Check bornfromreal.inc'
        write (*, *) 'Is fks_singular compiled correctly?'
        stop 1
      end if
      do i = 1, config_map(0, 0)
        xtot = xtot + amp2(config_map(i, 0))
      end do
      if (xtot .ne. 0d0) then
        enhance = amp2(config_map(this_config, 0))/xtot
        enhance = enhance*diagramsymmetryfactor
      else
        enhance = 0d0
      end if
    end if

    if (imode .eq. 1) then
      f_b = f_b*enhance
      f_nb = f_nb*enhance
    elseif (imode .eq. 2) then
      f_r = f_r*enhance
    elseif (imode .eq. 3) then
      f_s = f_s*enhance
      f_c = f_c*enhance
      f_dc = f_dc*enhance
      f_sc = f_sc*enhance
      f_dsc(1) = f_dsc(1)*enhance
      f_dsc(2) = f_dsc(2)*enhance
      f_dsc(3) = f_dsc(3)*enhance
      f_dsc(4) = f_dsc(4)*enhance
    end if
    call cpu_time(tAfter)
    tf_nb = tf_nb + (tAfter - tBefore)

    return
  end subroutine include_multichannel_enhance

  subroutine compute_prefactors_n1body(vegas_wgt, jac_ev)
! Compute all relevant prefactors for the real emission and counter
! terms.
    implicit none
    real :: tBefore, tAfter
    double precision vegas_wgt, prefact, prefact_cnt_ssc, prefact_deg
    double precision prefact_c, prefact_coll, jac_ev, pi
    double precision prefact_cnt_ssc_c, prefact_coll_c
    double precision prefact_deg_slxi, prefact_deg_sxi
    double precision collinear_shat, soft_collinear_shat
    integer i
    parameter(pi=3.1415926535897932385d0)
    call cpu_time(tBefore)
    collinear_shat = fks_subtraction_shat(collinear_counterevent)
    soft_collinear_shat = &
         fks_subtraction_shat(soft_collinear_counterevent)

! f_* multiplication factors for real-emission, soft counter, ... etc.
    prefact = event_xi_norm(real_event)/event_xi(real_event)/ &
              (1 - event_y(real_event))
    f_r = prefact*jac_ev*fkssymmetryfactor*vegas_wgt
    if (.not. nocntevents) then
      prefact_cnt_ssc = event_xi_norm(real_event)/ &
                        min(event_xi_max(real_event), xiScut_used)* &
                        log(xicut_used/min(event_xi_max(real_event), &
                                          xiScut_used))/ &
                        (1 - event_y(real_event))
      f_s = (prefact + prefact_cnt_ssc)* &
            stored_event_jacobian(soft_counterevent)* &
            fkssymmetryfactor*vegas_wgt
      if (fks_sister_is_massless()) then
! For the soft-collinear, these should be itwo. But they are always
! equal to ione, so no need to define separate factors.
        prefact_c = event_xi_norm(collinear_counterevent)/ &
                    event_xi(collinear_counterevent)/ &
                    (1 - event_y(real_event))
        prefact_coll = event_xi_norm(collinear_counterevent)/ &
                       event_xi(collinear_counterevent)* &
                       log(delta_used/deltaS)/deltaS
        f_c = (prefact_c + prefact_coll)* &
              stored_event_jacobian(collinear_counterevent)* &
              fkssymmetryfactor*vegas_wgt
        prefact_deg = event_xi_norm(collinear_counterevent)/ &
                      event_xi(collinear_counterevent)/deltaS
        prefact_cnt_ssc_c = event_xi_norm(collinear_counterevent)/ &
                            min(event_xi_max(collinear_counterevent), &
                                xiScut_used) &
                            *log(xicut_used/ &
                                 min(event_xi_max(collinear_counterevent), &
                                     xiScut_used)) &
                            /(1 - event_y(real_event))
        prefact_coll_c = event_xi_norm(collinear_counterevent)/ &
                         min(event_xi_max(collinear_counterevent), &
                             xiScut_used) &
                         *log(xicut_used/ &
                              min(event_xi_max(collinear_counterevent), &
                                  xiScut_used)) &
                         *log(delta_used/deltaS)/deltaS
        f_dc = stored_event_jacobian(collinear_counterevent)*prefact_deg/ &
               (collinear_shat/(32*pi**2))* &
               fkssymmetryfactorDeg*vegas_wgt
        f_sc = (prefact_c + prefact_coll + prefact_cnt_ssc_c + &
                prefact_coll_c)* &
               stored_event_jacobian(soft_collinear_counterevent)* &
               fkssymmetryfactorDeg*vegas_wgt
        prefact_deg_sxi = event_xi_norm(collinear_counterevent)/ &
                          min(event_xi_max(collinear_counterevent), &
                              xiScut_used)* &
                          log(xicut_used/ &
                              min(event_xi_max(collinear_counterevent), &
                                  xiScut_used))*1/deltaS
        prefact_deg_slxi = event_xi_norm(collinear_counterevent)/ &
                           min(event_xi_max(collinear_counterevent), &
                               xiScut_used) &
                           *(log(xicut_used)**2 &
                             - log(min(event_xi_max(collinear_counterevent), &
                                       xiScut_used))**2) &
                           /(2.d0*deltaS)
        f_dsc(1) = prefact_deg* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (soft_collinear_shat/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(2) = prefact_deg_sxi* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (soft_collinear_shat/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(3) = prefact_deg_slxi* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (soft_collinear_shat/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(4) = (prefact_deg + prefact_deg_sxi)* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (soft_collinear_shat/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
      else
        f_c = 0d0
        f_dc = 0d0
        f_sc = 0d0
        do i = 1, 4
          f_dsc(i) = 0d0
        end do
      end if
    else
      f_s = 0d0
      f_c = 0d0
      f_dc = 0d0
      f_sc = 0d0
      do i = 1, 4
        f_dsc(i) = 0d0
      end do
    end if
    call cpu_time(tAfter)
    tf_all = tf_all + (tAfter - tBefore)
    return
  end subroutine compute_prefactors_n1body


end module fks_contributions_module

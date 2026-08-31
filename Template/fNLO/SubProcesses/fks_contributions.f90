module fks_contributions_module
  use process_dimensions, only: nexternal, nsplitorders, qcd_pos, &
                                amp_split_size
  use timing_state, only: tBorn, tIS, tReal, tCount, tf_nb, tf_all
  use mint_module, only: virt_wgt_mint, born_wgt_mint
  use split_orders, only: get_orders_tag, orders_to_amp_split_pos, &
                          amp_split_pos_to_orders
  use madfks_plot_module, only: initplot_impl
  use fks_model_state_module, only: g => strong_coupling, external_masses
  use factorized_phase_space, only: factorized_radiation_state, &
       fetch_factorized_radiation_state, compose_factorized_block_measure
  use decay_chain_metadata, only: has_decay_chains
  use decay_chain_kinematics, only: fks_leg_mass
  use nlo_decay_metadata, only: has_nlo_decay, nlo_decay_corrected_node
  use nlo_decay_kinematics, only: nlo_decay_fks_sister_mass
  use nlo_contribution_bundle, only: active_contribution_has_virtual, &
       active_virtual_grid_index, active_contribution_is_production
  use fks_singular_module, only: evaluate_fks_sij, sreal, sreal_deg, &
                                 bornsoftvirtual, fks_subtraction_shat, &
                                 evaluate_born_matrix
  use spin_density_fks_matrices, only: &
       spin_density_fks_collection_enabled, &
       spin_density_active_open_size, &
       spin_density_active_component_position, &
       spin_density_born_amp_position, spin_density_nlo_amp_position, &
       get_spin_density_born_matrix, get_spin_density_reduced_matrix, &
       get_spin_density_degenerate_matrix, &
       get_spin_density_integrated_matrix, &
       get_spin_density_virtual_matrix, &
       spin_density_born_matrix_available, &
       spin_density_reduced_matrix_available, &
       spin_density_degenerate_matrix_available, &
       spin_density_integrated_matrix_available, &
       spin_density_virtual_matrix_available
  use spin_density_matrix_results, only: spin_density_bornlike_branch, &
       spin_density_real_branch
  use decay_chain_scales, only: active_block_qcd_squared_order
  use decay_chain_parameters, only: multiplicative_nlo_enabled
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

  ! The ordinary f_* factors contain the complete factorized event measure
  ! and remain untouched for the expanded calculation.  The density factors
  ! below contain only the active physical block.  Their product over all
  ! production/decay blocks therefore reconstructs the global measure once.
  double precision, save :: density_f_b = 0d0
  double precision, save :: density_f_nb = 0d0
  double precision, save :: density_f_r = 0d0
  double precision, save :: density_f_s = 0d0
  double precision, save :: density_f_c = 0d0
  double precision, save :: density_f_dc = 0d0
  double precision, save :: density_f_sc = 0d0
  double precision, save :: density_f_dsc(4) = 0d0

  public :: compute_born, compute_nbody_noborn, compute_real_emission
  public :: compute_decay_width_counterterm
  public :: compute_soft_counter_term, compute_collinear_counter_term
  public :: compute_soft_collinear_ct_impl
  public :: compute_prefactors_nbody, compute_prefactors_n1body
  public :: record_multiplicative_spectator_born_densities
  public :: include_multichannel_enhance

  interface
    integer function sdm_branch_component_count()
    end function sdm_branch_component_count

    integer function sdm_branch_max_open_size()
    end function sdm_branch_max_open_size

    integer function sdm_branch_component_id(position)
      integer, intent(in) :: position
    end function sdm_branch_component_id

  end interface

contains

  subroutine active_block_measure(event_slot, include_event, measure)
    integer, intent(in) :: event_slot
    logical, intent(in) :: include_event
    double precision, intent(out) :: measure
    double precision :: jacobian, phase_space_weight
    logical :: available

    call compose_factorized_block_measure( &
         active_radiation_block(), event_slot, include_event, &
         active_contribution_is_production(), jacobian, &
         phase_space_weight, available)
    if (.not. available) then
      write (*, '(a,i0,a,i0)') &
           'ERROR: block-local measure is unavailable for block ', &
           active_radiation_block(), ' and event slot ', event_slot
      stop 1
    end if
    measure = jacobian*phase_space_weight
  end subroutine active_block_measure


  double precision function active_block_vegas_weight(vegas_wgt)
    double precision, intent(in) :: vegas_wgt

    active_block_vegas_weight = 1d0
    if (active_contribution_is_production()) then
      active_block_vegas_weight = vegas_wgt
    end if
  end function active_block_vegas_weight

  logical function uses_factorized_radiation_state()
    if (has_nlo_decay()) then
      uses_factorized_radiation_state = .true.
    else
      uses_factorized_radiation_state = has_decay_chains()
    end if
  end function uses_factorized_radiation_state


  integer function active_radiation_block()
    if (has_nlo_decay()) then
      active_radiation_block = nlo_decay_corrected_node()
    else
      active_radiation_block = 0
    end if
  end function active_radiation_block


  subroutine load_radiation_state(event_slot, radiation)
    integer, intent(in) :: event_slot
    type(factorized_radiation_state), intent(out) :: radiation
    logical :: available

    if (uses_factorized_radiation_state()) then
      call fetch_factorized_radiation_state( &
           event_slot, active_radiation_block(), radiation, available)
      if (available) return
      ! Some massive mappings intentionally omit counterevents.  They retain
      ! negative legacy sentinels and are never evaluated as physical blocks.
      if (stored_event_jacobian(event_slot) > 0d0) then
        write (*, '(a,i0)') &
             'ERROR: block-local radiation state is unavailable for slot ', &
             event_slot
        stop 1
      end if
      ! A skipped massive counterevent has no independent mapping.  Retain
      ! its negative Jacobian sentinel, but inherit all unused kinematic
      ! fields from the physical real state of the same block.
      call fetch_factorized_radiation_state( &
           real_event, active_radiation_block(), radiation, available)
      if (.not. available) then
        write (*, '(a)') &
             'ERROR: real block-local radiation state is unavailable'
        stop 1
      end if
      radiation%jacobian = stored_event_jacobian(event_slot)
      return
    end if
    radiation = factorized_radiation_state()
    radiation%jacobian = stored_event_jacobian(event_slot)
    radiation%xi = event_xi(event_slot)
    radiation%y = event_y(event_slot)
    radiation%xi_hat = event_xi_hat(event_slot)
    radiation%xi_max = event_xi_max(event_slot)
    radiation%xi_norm = event_xi_norm(event_slot)
    radiation%shat = event_shat(event_slot)
    radiation%sqrt_shat = sqrt(max(0d0, event_shat(event_slot)))
  end subroutine load_radiation_state

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
    integer :: open_size, component, local_qcd_power

    double precision wgt_c
    double precision wgt1
    complex(kind=8), allocatable :: raw_density(:, :, :)
    complex(kind=8), allocatable :: density_coefficients(:, :, :)
    type(factorized_radiation_state) :: real_radiation, soft_radiation
    call cpu_time(tBefore)
    if (f_b .eq. 0d0) return
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(soft_counterevent, soft_radiation)
    if (real_radiation%xi_hat*soft_radiation%xi_max .gt. &
        xiBSVcut_used) return
    call evaluate_born_matrix(soft_counterevent, wgt_c)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*f_b/g**(qcd_power)
      if (spin_density_fks_collection_enabled() .and. &
          iamp == spin_density_born_amp_position() .and. &
          (.not. multiplicative_nlo_enabled() .or. &
           active_contribution_is_production())) then
        open_size = spin_density_active_open_size()
        component = spin_density_active_component_position()
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        allocate(raw_density(2, open_size, open_size))
        allocate(density_coefficients(3, open_size, open_size))
        call get_spin_density_born_matrix(raw_density)
        density_coefficients = (0d0, 0d0)
        density_coefficients(1, :, :) = raw_density(1, :, :)* &
             density_f_b/g**local_qcd_power
        call add_wgt(soft_counterevent, born_contribution, wgt1, 0d0, 0d0, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_bornlike_branch)
        deallocate(raw_density, density_coefficients)
      else
        call add_wgt(soft_counterevent, born_contribution, wgt1, 0d0, 0d0)
      end if
    end do

    call record_multiplicative_spectator_born_densities()
    call cpu_time(tAfter)
    tBorn = tBorn + (tAfter - tBefore)
    return
  end subroutine compute_born


  subroutine record_multiplicative_spectator_born_densities( &
       excluded_component)
    integer, intent(in), optional :: excluded_component
    integer :: component, component_count, open_size, maximum_open_size
    integer :: production_component
    integer :: qcd_power, scale_pdg, block
    double precision :: jacobian, phase_space_weight, measure
    logical :: available
    complex(kind=8), allocatable :: raw_density(:, :, :)
    complex(kind=8), allocatable :: density_coefficients(:, :, :)
    external :: sdm_component_born_density

    if (.not. spin_density_fks_collection_enabled() .or. &
        .not. multiplicative_nlo_enabled() .or. &
        .not. active_contribution_is_production()) return
    component_count = sdm_branch_component_count()
    maximum_open_size = sdm_branch_max_open_size()
    if (component_count < 1 .or. maximum_open_size < 1) then
      write (*, *) 'ERROR: invalid spectator density layout'
      stop 1
    end if
    production_component = spin_density_active_component_position()
    if (present(excluded_component)) production_component = &
         excluded_component
    if (production_component < 1 .or. &
        production_component > component_count) then
      write (*, *) 'ERROR: invalid spectator density exclusion', &
           production_component
      stop 1
    end if
    allocate(raw_density(2, maximum_open_size, maximum_open_size))
    do component = 1, component_count
      ! Production defines the common global Born embedding.  Record every
      ! other block here, including NLO-corrected decays; their own FKS
      ! contributions subsequently add only the local NLO increment.
      if (component == production_component) cycle
      raw_density = (0d0, 0d0)
      call sdm_component_born_density( &
           component, soft_counterevent, open_size, qcd_power, &
           scale_pdg, raw_density)
      if (open_size < 1 .or. open_size > maximum_open_size .or. &
          qcd_power < 0) then
        write (*, *) 'ERROR: invalid spectator Born density metadata'
        stop 1
      end if
      block = sdm_branch_component_id(component)
      call compose_factorized_block_measure( &
           block, soft_counterevent, .false., .false., jacobian, &
           phase_space_weight, available)
      if (.not. available) then
        write (*, *) 'ERROR: spectator Born measure is unavailable', block
        stop 1
      end if
      measure = jacobian*phase_space_weight
      allocate(density_coefficients(3, open_size, open_size))
      density_coefficients = (0d0, 0d0)
      density_coefficients(1, :, :) = &
           raw_density(1, 1:open_size, 1:open_size)*measure
      if (qcd_power > 0) then
        if (g <= 0d0 .or. scale_pdg == 0) then
          write (*, *) 'ERROR: invalid spectator QCD normalization'
          stop 1
        end if
        density_coefficients(1, :, :) = &
             density_coefficients(1, :, :)/g**qcd_power
      end if
      call add_wgt( &
           soft_counterevent, born_contribution, 0d0, 0d0, 0d0, &
           density_coefficients=density_coefficients, &
           density_component=component, &
           density_branch=spin_density_bornlike_branch, &
           density_qcd_power=qcd_power, density_scale_pdg=scale_pdg, &
           density_is_production=.false.)
      deallocate(density_coefficients)
    end do
    deallocate(raw_density)
  end subroutine record_multiplicative_spectator_born_densities


  subroutine compute_decay_width_counterterm
! Add the denominator term in the strict fixed-order expansion
!
!   product_i 1/Gamma_i = product_i 1/Gamma_i^(0)
!       * (1 - sum_i delta_Gamma_i/Gamma_i^(0)) + O(alpha_s^2).
!
! The supplied width difference already is an O(alpha_s) quantity.  It must
! therefore multiply the Born weight without changing its explicit coupling
! powers; doing otherwise would introduce an additional running coupling and
! would no longer represent the user-supplied NLO total width.
    use extra_weights
    implicit none
    integer :: born_orders(nsplitorders), correction_orders(nsplitorders)
    integer :: iamp
    double precision :: born_weight, weight
    type(factorized_radiation_state) :: real_radiation, soft_radiation

    if (multiplicative_nlo_enabled()) return
    if (.not. active_contribution_is_production()) return
    if (f_b .eq. 0d0) return
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(soft_counterevent, soft_radiation)
    if (real_radiation%xi_hat*soft_radiation%xi_max .gt. &
        xiBSVcut_used) return
    call evaluate_born_matrix(soft_counterevent, born_weight)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, born_orders)
      correction_orders = born_orders
      correction_orders(qcd_pos) = correction_orders(qcd_pos) + 2
      ! The numerical width difference already contains its alpha_s, so the
      ! coupling rescaling must remain that of the Born.  Nevertheless tag
      ! this line with the NLO squared order so plots and contribution
      ! bookkeeping cannot mistake the denominator counterterm for LO.
      QCD_power = born_orders(qcd_pos)
      orders_tag = get_orders_tag(correction_orders)
      amp_pos = orders_to_amp_split_pos(correction_orders)
      ! Store the raw Born line.  The width coefficient belongs to the
      ! selected scale point and is applied when the central/reweighted
      ! weight is evaluated.  This also permits a zero central coefficient
      ! with non-zero varied coefficients.
      weight = amp_split(iamp)*f_b/g**QCD_power
      call add_wgt(soft_counterevent, integrated_contribution, &
                   weight, 0d0, 0d0, is_width_counterterm=.true.)
    end do
  end subroutine compute_decay_width_counterterm

  subroutine compute_nbody_noborn
! This subroutine computes the soft-virtual matrix elements and adds its
! value to the list of weights using the add_wgt subroutine
    use extra_weights
    use mint_module
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp, virtual_grid
    integer :: open_size, component, local_qcd_power

    double precision wgt1, wgt2, wgt3, bsv_wgt, virt_wgt, born_wgt, g22, wgt4
    double precision :: averaged_virtual_factor
    logical :: born_density_available, integrated_density_available
    logical :: virtual_density_available
    complex(kind=8), allocatable :: raw_born_density(:, :, :)
    complex(kind=8), allocatable :: integrated_density(:, :, :)
    complex(kind=8), allocatable :: virtual_density(:, :, :)
    complex(kind=8), allocatable :: density_coefficients(:, :, :)
    type(factorized_radiation_state) :: real_radiation, soft_radiation
    call cpu_time(tBefore)
    if (f_nb .eq. 0d0) return
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(soft_counterevent, soft_radiation)
    if (real_radiation%xi_hat*soft_radiation%xi_max .gt. &
        xiBSVcut_used) return
    call bornsoftvirtual( &
         soft_counterevent, bsv_wgt, virt_wgt, born_wgt)
    born_density_available = spin_density_fks_collection_enabled() .and. &
         spin_density_born_matrix_available()
    integrated_density_available = &
         spin_density_fks_collection_enabled() .and. &
         spin_density_integrated_matrix_available()
    virtual_density_available = &
         spin_density_fks_collection_enabled() .and. &
         spin_density_virtual_matrix_available()
    if (born_density_available .or. integrated_density_available .or. &
        virtual_density_available) then
      open_size = spin_density_active_open_size()
      component = spin_density_active_component_position()
      allocate(density_coefficients(3, open_size, open_size))
    end if
    if (born_density_available) then
      allocate(raw_born_density(2, open_size, open_size))
      call get_spin_density_born_matrix(raw_born_density)
    end if
    if (integrated_density_available) then
      allocate(integrated_density(3, open_size, open_size))
      call get_spin_density_integrated_matrix(integrated_density)
    end if
    if (virtual_density_available) then
      allocate(virtual_density(3, open_size, open_size))
      call get_spin_density_virtual_matrix(virtual_density)
    end if
    do iamp = 1, amp_split_size
      if (amp_split_wgtnstmp(iamp) .eq. 0d0 .and. &
          amp_split_wgtwnstmpmur(iamp) .eq. 0d0 .and. &
          amp_split_wgtwnstmpmuf(iamp) .eq. 0d0 .and. &
          amp_split_avv(iamp) .eq. 0d0 .and. &
          .not. (integrated_density_available .and. &
                 iamp == spin_density_nlo_amp_position()) .and. &
          .not. (born_density_available .and. &
                 iamp == spin_density_nlo_amp_position())) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = amp_split_wgtnstmp(iamp)*f_nb/g22
      wgt2 = amp_split_wgtwnstmpmur(iamp)*f_nb/g22
      wgt3 = amp_split_wgtwnstmpmuf(iamp)*f_nb/g22
      wgt4 = amp_split_avv(iamp)*f_nb/g22
      if (integrated_density_available .and. &
          iamp == spin_density_nlo_amp_position()) then
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        density_coefficients = integrated_density*density_f_nb/ &
             g**local_qcd_power
        call add_wgt(soft_counterevent, integrated_contribution, &
             wgt1, wgt2, wgt3, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_bornlike_branch)
      else
        call add_wgt(soft_counterevent, integrated_contribution, &
                     wgt1, wgt2, wgt3)
      end if
      if (born_density_available .and. &
          iamp == spin_density_nlo_amp_position() .and. &
          amp_split_born_for_virt(iamp) /= 0d0) then
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        averaged_virtual_factor = &
             amp_split_avv(iamp)/amp_split_born_for_virt(iamp)
        density_coefficients = (0d0, 0d0)
        density_coefficients(1, :, :) = raw_born_density(1, :, :)* &
             averaged_virtual_factor*density_f_nb/g**local_qcd_power
        call add_wgt(soft_counterevent, averaged_born_contribution, &
             wgt4, 0d0, 0d0, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_bornlike_branch)
      else
        call add_wgt(soft_counterevent, averaged_born_contribution, &
                     wgt4, 0d0, 0d0)
      end if
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
      if (amp_split_virt(iamp) .eq. 0d0 .and. &
          .not. (virtual_density_available .and. &
                 iamp == spin_density_nlo_amp_position())) cycle
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
      if (virtual_density_available .and. &
          iamp == spin_density_nlo_amp_position()) then
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        density_coefficients = (0d0, 0d0)
        density_coefficients(1, :, :) = virtual_density(1, :, :)* &
             density_f_nb/g**local_qcd_power
        call add_wgt(soft_counterevent, virtual_contribution, &
             wgt1, 0d0, 0d0, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_bornlike_branch)
      else
        call add_wgt(soft_counterevent, virtual_contribution, &
                     wgt1, 0d0, 0d0)
      end if
    end do

    if (allocated(raw_born_density)) deallocate(raw_born_density)
    if (allocated(integrated_density)) deallocate(integrated_density)
    if (allocated(virtual_density)) deallocate(virtual_density)
    if (allocated(density_coefficients)) deallocate(density_coefficients)

    call cpu_time(tAfter)
    tIS = tIS + (tAfter - tBefore)
    return
  end subroutine compute_nbody_noborn

  subroutine compute_real_emission()
! This subroutine computes the real-emission matrix elements and adds
! its value to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
    integer :: open_size, component, local_qcd_power
    double precision s_ev, wgt1, fx_ev
    complex(kind=8), allocatable :: reduced_density(:, :)
    complex(kind=8), allocatable :: density_coefficients(:, :, :)
    type(factorized_radiation_state) :: real_radiation
    call cpu_time(tBefore)
    if (f_r .eq. 0d0) return
    call load_radiation_state(real_event, real_radiation)
    s_ev = evaluate_fks_sij(real_event, i_fks, j_fks, &
                            real_radiation%xi, real_radiation%y)
    if (s_ev .le. 0.d0) return
    call sreal(real_event, real_radiation%xi, real_radiation%y, fx_ev)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*s_ev*f_r/g**(qcd_power)
      if (spin_density_fks_collection_enabled() .and. &
          iamp == spin_density_nlo_amp_position() .and. &
          spin_density_reduced_matrix_available()) then
        open_size = spin_density_active_open_size()
        component = spin_density_active_component_position()
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        allocate(reduced_density(open_size, open_size))
        allocate(density_coefficients(3, open_size, open_size))
        call get_spin_density_reduced_matrix(reduced_density)
        density_coefficients = (0d0, 0d0)
        density_coefficients(1, :, :) = reduced_density*s_ev* &
             density_f_r/g**local_qcd_power
        call add_wgt(real_event, real_contribution, wgt1, 0d0, 0d0, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_real_branch)
        deallocate(reduced_density, density_coefficients)
      else
        call add_wgt(real_event, real_contribution, wgt1, 0d0, 0d0)
      end if
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
    integer :: open_size, component, local_qcd_power
    double precision wgt1, s_s, fx_s, zero, g22
    logical :: density_available
    complex(kind=8), allocatable :: reduced_density(:, :)
    complex(kind=8), allocatable :: density_coefficients(:, :, :)
    type(factorized_radiation_state) :: real_radiation, soft_radiation
    parameter(zero=0d0)
    call cpu_time(tBefore)
    if (f_s .eq. 0d0) return
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(soft_counterevent, soft_radiation)
    if (real_radiation%xi_hat*soft_radiation%xi_max .gt. &
        xiScut_used) return
    s_s = evaluate_fks_sij(soft_counterevent, i_fks, j_fks, &
                           zero, real_radiation%y)
    if (s_s .le. 0d0) return
    call sreal( &
         soft_counterevent, 0d0, real_radiation%y, fx_s)
    density_available = spin_density_fks_collection_enabled() .and. &
         spin_density_reduced_matrix_available() .and. &
         real_radiation%xi <= xiScut_used
    if (density_available) then
      open_size = spin_density_active_open_size()
      component = spin_density_active_component_position()
      allocate(reduced_density(open_size, open_size))
      allocate(density_coefficients(3, open_size, open_size))
      call get_spin_density_reduced_matrix(reduced_density)
    end if

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0 .and. &
          .not. (density_available .and. &
                 iamp == spin_density_nlo_amp_position())) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = 0d0
      if (real_radiation%xi .le. xiScut_used) then
        wgt1 = -amp_split(iamp)*s_s*f_s/g22
      end if
      if (density_available .and. &
          iamp == spin_density_nlo_amp_position()) then
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        density_coefficients = (0d0, 0d0)
        density_coefficients(1, :, :) = -reduced_density*s_s* &
             density_f_s/g**local_qcd_power
        call add_wgt(soft_counterevent, soft_contribution, wgt1, 0d0, 0d0, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_bornlike_branch)
      else if (wgt1 .ne. 0d0) then
        call add_wgt(soft_counterevent, soft_contribution, wgt1, 0d0, 0d0)
      end if
    end do
    if (density_available) then
      deallocate(reduced_density, density_coefficients)
    end if

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
    integer :: open_size, component, local_qcd_power
    double precision one, s_c, fx_c, deg_xi_c, deg_lxi_c, wgt1, wgt3, g22
    logical :: density_available
    complex(kind=8), allocatable :: reduced_density(:, :)
    complex(kind=8), allocatable :: degenerate_density(:, :, :)
    complex(kind=8), allocatable :: density_coefficients(:, :, :)
    type(factorized_radiation_state) :: real_radiation, collinear_radiation
    parameter(one=1d0)
    call cpu_time(tBefore)
    if (f_c .eq. 0d0 .and. f_dc .eq. 0d0) return
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(collinear_counterevent, collinear_radiation)
    if (real_radiation%y .le. 1d0 - deltaS .or. &
        .not. fks_sister_is_massless()) return
    s_c = evaluate_fks_sij(collinear_counterevent, i_fks, j_fks, &
                           collinear_radiation%xi, one)
    if (s_c .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg(collinear_counterevent, collinear_radiation%xi, &
                   deg_xi_c, deg_lxi_c)
    call sreal(collinear_counterevent, collinear_radiation%xi, one, fx_c)
    density_available = spin_density_fks_collection_enabled() .and. &
         spin_density_reduced_matrix_available()
    if (density_available) then
      open_size = spin_density_active_open_size()
      component = spin_density_active_component_position()
      allocate(reduced_density(open_size, open_size))
      allocate(degenerate_density(3, open_size, open_size))
      allocate(density_coefficients(3, open_size, open_size))
      call get_spin_density_reduced_matrix(reduced_density)
      degenerate_density = (0d0, 0d0)
      if (spin_density_degenerate_matrix_available()) then
        call get_spin_density_degenerate_matrix(degenerate_density)
      end if
    end if

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_xi(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_lxi(iamp) .eq. 0d0 .and. &
          .not. (density_available .and. &
                 iamp == spin_density_nlo_amp_position())) cycle

      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = -amp_split(iamp)*s_c*f_c/g22
      wgt1 = wgt1 + (amp_split_wgtdegrem_xi(iamp) + &
                     amp_split_wgtdegrem_lxi(iamp)* &
                     log(collinear_radiation%xi))*f_dc/g22
      wgt3 = amp_split_wgtdegrem_muF(iamp)*f_dc/g22
      if (density_available .and. &
          iamp == spin_density_nlo_amp_position()) then
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        density_coefficients = (0d0, 0d0)
        density_coefficients(1, :, :) = ( &
             -reduced_density*s_c*density_f_c + &
             (degenerate_density(1, :, :) + &
              degenerate_density(2, :, :)* &
              log(collinear_radiation%xi))*density_f_dc)/ &
             g**local_qcd_power
        density_coefficients(3, :, :) = &
             degenerate_density(3, :, :)*density_f_dc/ &
             g**local_qcd_power
        call add_wgt(collinear_counterevent, collinear_contribution, &
             wgt1, 0d0, wgt3, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_bornlike_branch)
      else if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) then
        call add_wgt(collinear_counterevent, collinear_contribution, &
                     wgt1, 0d0, wgt3)
      end if
    end do
    if (density_available) then
      deallocate(reduced_density, degenerate_density, density_coefficients)
    end if

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
    integer :: open_size, component, local_qcd_power
    double precision zero, one, s_sc, fx_sc, wgt1, wgt3, deg_xi_sc, deg_lxi_sc, g22
    logical :: density_available
    complex(kind=8), allocatable :: reduced_density(:, :)
    complex(kind=8), allocatable :: degenerate_density(:, :, :)
    complex(kind=8), allocatable :: density_coefficients(:, :, :)
    type(factorized_radiation_state) :: real_radiation
    type(factorized_radiation_state) :: collinear_radiation
    parameter(zero=0d0, one=1d0)
    call cpu_time(tBefore)
    if (f_sc .eq. 0d0 .and. f_dsc(1) .eq. 0d0 .and. f_dsc(2) .eq. 0d0 .and. f_dsc(3) .eq. 0d0 .and. f_dsc(4) .eq. 0d0) return
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(collinear_counterevent, collinear_radiation)
    if (real_radiation%xi_hat*collinear_radiation%xi_max &
        .ge. xiScut_used .or. real_radiation%y .le. 1d0 - deltaS &
        .or. .not. fks_sister_is_massless()) return
    s_sc = evaluate_fks_sij( &
         soft_collinear_counterevent, i_fks, j_fks, zero, one)
    if (s_sc .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg( &
         soft_collinear_counterevent, zero, deg_xi_sc, deg_lxi_sc)
    call sreal(soft_collinear_counterevent, zero, one, fx_sc)
    density_available = spin_density_fks_collection_enabled() .and. &
         spin_density_reduced_matrix_available() .and. &
         collinear_radiation%xi < xiScut_used
    if (density_available) then
      open_size = spin_density_active_open_size()
      component = spin_density_active_component_position()
      allocate(reduced_density(open_size, open_size))
      allocate(degenerate_density(3, open_size, open_size))
      allocate(density_coefficients(3, open_size, open_size))
      call get_spin_density_reduced_matrix(reduced_density)
      degenerate_density = (0d0, 0d0)
      if (spin_density_degenerate_matrix_available()) then
        call get_spin_density_degenerate_matrix(degenerate_density)
      end if
    end if

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_xi(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_lxi(iamp) .eq. 0d0 .and. &
          .not. (density_available .and. &
                 iamp == spin_density_nlo_amp_position())) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = 0d0
      wgt3 = 0d0
      if (collinear_radiation%xi .lt. xiScut_used) then
        wgt1 = amp_split(iamp)*s_sc*f_sc/g22
        wgt1 = wgt1 + ( &
               -(amp_split_wgtdegrem_xi(iamp) + &
                 amp_split_wgtdegrem_lxi(iamp)* &
                 log(collinear_radiation%xi))*f_dsc(1) &
               - (amp_split_wgtdegrem_xi(iamp)*f_dsc(2) + &
                  amp_split_wgtdegrem_lxi(iamp)*f_dsc(3)))/g22
        wgt3 = -amp_split_wgtdegrem_muF(iamp)*f_dsc(4)/g22
      end if
      if (density_available .and. &
          iamp == spin_density_nlo_amp_position()) then
        local_qcd_power = active_block_qcd_squared_order(QCD_power)
        density_coefficients = (0d0, 0d0)
        density_coefficients(1, :, :) = ( &
             reduced_density*s_sc*density_f_sc - &
             (degenerate_density(1, :, :) + &
              degenerate_density(2, :, :)* &
              log(collinear_radiation%xi))*density_f_dsc(1) - &
             degenerate_density(1, :, :)*density_f_dsc(2) - &
             degenerate_density(2, :, :)*density_f_dsc(3))/ &
             g**local_qcd_power
        density_coefficients(3, :, :) = &
             -degenerate_density(3, :, :)*density_f_dsc(4)/ &
             g**local_qcd_power
        call add_wgt(soft_collinear_counterevent, &
             soft_collinear_contribution, wgt1, 0d0, wgt3, &
             density_coefficients=density_coefficients, &
             density_component=component, &
             density_branch=spin_density_bornlike_branch)
      else if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) then
        call add_wgt(soft_collinear_counterevent, &
                     soft_collinear_contribution, wgt1, 0d0, wgt3)
      end if
    end do
    if (density_available) then
      deallocate(reduced_density, degenerate_density, density_coefficients)
    end if

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
    double precision :: block_measure, block_vegas
    type(factorized_radiation_state) :: real_radiation, soft_radiation
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
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(soft_counterevent, soft_radiation)
    subtraction_shat = fks_subtraction_shat(soft_counterevent)
    f_b = soft_radiation%jacobian* &
          real_radiation%xi_norm/ &
          (min(real_radiation%xi_max, xiBSVcut_used)* &
           subtraction_shat/(16*pi**2))* &
          fkssymmetryfactorBorn*vegas_wgt
    f_nb = f_b
    density_f_b = f_b
    density_f_nb = f_nb
    if (spin_density_fks_collection_enabled()) then
      call active_block_measure(soft_counterevent, .true., block_measure)
      block_vegas = active_block_vegas_weight(vegas_wgt)
      density_f_b = block_measure* &
           real_radiation%xi_norm/ &
           (min(real_radiation%xi_max, xiBSVcut_used)* &
            subtraction_shat/(16*pi**2))* &
           fkssymmetryfactorBorn*block_vegas
      density_f_nb = density_f_b
    end if
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
    logical single_channel_shortcut
    data inoborn_cnt/0/

    call cpu_time(tBefore)

! Compute the multi-channel enhancement factor 'enhance'.  With one Born
! channel the normalized partition is identically one, so evaluating the
! full helicity-summed Born solely to divide it by itself is wasted work.
    enhance = 1.d0
    single_channel_shortcut = &
         p_born(0, 1) > 0d0 .and. config_map(0, 0) == 1
    if (single_channel_shortcut) then
      enhance = diagramsymmetryfactor
    else if (p_born(0, 1) .gt. 0d0) then
      call evaluate_born_matrix(soft_counterevent, wgt_c)
    elseif (p_born(0, 1) .lt. 0d0) then
      enhance = 0d0
    end if
    if (enhance .eq. 0d0) then
      xnoborn_cnt = xnoborn_cnt + 1.d0
      if (log10(xnoborn_cnt) .gt. inoborn_cnt) then
        write (*, *) 'WARNING: no Born momenta more than 10**', inoborn_cnt, 'times'
        inoborn_cnt = inoborn_cnt + 1
      end if
    else if (config_map(0, 0) <= 0) then
      write (*, *) 'Fatal error in compute_prefactor_nbody:'// &
           ' no Born diagrams ', config_map, '. Check bornfromreal.inc'
      write (*, *) 'Is fks_singular compiled correctly?'
      stop 1
    else if (.not. single_channel_shortcut) then
      xtot = 0d0
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
      density_f_b = density_f_b*enhance
      density_f_nb = density_f_nb*enhance
    elseif (imode .eq. 2) then
      f_r = f_r*enhance
      density_f_r = density_f_r*enhance
    elseif (imode .eq. 3) then
      f_s = f_s*enhance
      f_c = f_c*enhance
      f_dc = f_dc*enhance
      f_sc = f_sc*enhance
      f_dsc(1) = f_dsc(1)*enhance
      f_dsc(2) = f_dsc(2)*enhance
      f_dsc(3) = f_dsc(3)*enhance
      f_dsc(4) = f_dsc(4)*enhance
      density_f_s = density_f_s*enhance
      density_f_c = density_f_c*enhance
      density_f_dc = density_f_dc*enhance
      density_f_sc = density_f_sc*enhance
      density_f_dsc = density_f_dsc*enhance
    end if
    call cpu_time(tAfter)
    tf_nb = tf_nb + (tAfter - tBefore)

    return
  end subroutine include_multichannel_enhance

  subroutine compute_prefactors_n1body(vegas_wgt)
! Compute all relevant prefactors for the real emission and counter
! terms.
    implicit none
    real :: tBefore, tAfter
    double precision vegas_wgt, prefact, prefact_cnt_ssc, prefact_deg
    double precision prefact_c, prefact_coll, pi
    double precision prefact_cnt_ssc_c, prefact_coll_c
    double precision prefact_deg_slxi, prefact_deg_sxi
    double precision collinear_shat, soft_collinear_shat
    double precision :: real_block_measure, soft_block_measure
    double precision :: collinear_block_measure
    double precision :: soft_collinear_block_measure, block_vegas
    type(factorized_radiation_state) :: real_radiation, soft_radiation
    type(factorized_radiation_state) :: collinear_radiation
    type(factorized_radiation_state) :: soft_collinear_radiation
    integer i
    parameter(pi=3.1415926535897932385d0)
    call cpu_time(tBefore)
    call load_radiation_state(real_event, real_radiation)
    call load_radiation_state(soft_counterevent, soft_radiation)
    call load_radiation_state(collinear_counterevent, collinear_radiation)
    call load_radiation_state(soft_collinear_counterevent, &
                              soft_collinear_radiation)
    collinear_shat = fks_subtraction_shat(collinear_counterevent)
    soft_collinear_shat = &
         fks_subtraction_shat(soft_collinear_counterevent)

! f_* multiplication factors for real-emission, soft counter, ... etc.
    prefact = real_radiation%xi_norm/real_radiation%xi/ &
              (1 - real_radiation%y)
    f_r = prefact*real_radiation%jacobian*fkssymmetryfactor*vegas_wgt
    density_f_r = f_r
    if (.not. nocntevents) then
      prefact_cnt_ssc = real_radiation%xi_norm/ &
                        min(real_radiation%xi_max, xiScut_used)* &
                        log(xicut_used/min(real_radiation%xi_max, &
                                          xiScut_used))/ &
                        (1 - real_radiation%y)
      f_s = (prefact + prefact_cnt_ssc)* &
            soft_radiation%jacobian* &
            fkssymmetryfactor*vegas_wgt
      if (fks_sister_is_massless()) then
! For the soft-collinear, these should be itwo. But they are always
! equal to ione, so no need to define separate factors.
        prefact_c = collinear_radiation%xi_norm/ &
                    collinear_radiation%xi/ &
                    (1 - real_radiation%y)
        prefact_coll = collinear_radiation%xi_norm/ &
                       collinear_radiation%xi* &
                       log(delta_used/deltaS)/deltaS
        f_c = (prefact_c + prefact_coll)* &
              collinear_radiation%jacobian* &
              fkssymmetryfactor*vegas_wgt
        prefact_deg = collinear_radiation%xi_norm/ &
                      collinear_radiation%xi/deltaS
        prefact_cnt_ssc_c = collinear_radiation%xi_norm/ &
                            min(collinear_radiation%xi_max, &
                                xiScut_used) &
                            *log(xicut_used/ &
                                 min(collinear_radiation%xi_max, &
                                     xiScut_used)) &
                            /(1 - real_radiation%y)
        prefact_coll_c = collinear_radiation%xi_norm/ &
                         min(collinear_radiation%xi_max, &
                             xiScut_used) &
                         *log(xicut_used/ &
                              min(collinear_radiation%xi_max, &
                                  xiScut_used)) &
                         *log(delta_used/deltaS)/deltaS
        f_dc = collinear_radiation%jacobian*prefact_deg/ &
               (collinear_shat/(32*pi**2))* &
               fkssymmetryfactorDeg*vegas_wgt
        f_sc = (prefact_c + prefact_coll + prefact_cnt_ssc_c + &
               prefact_coll_c)* &
               soft_collinear_radiation%jacobian* &
               fkssymmetryfactorDeg*vegas_wgt
        prefact_deg_sxi = collinear_radiation%xi_norm/ &
                          min(collinear_radiation%xi_max, &
                              xiScut_used)* &
                          log(xicut_used/ &
                              min(collinear_radiation%xi_max, &
                                  xiScut_used))*1/deltaS
        prefact_deg_slxi = collinear_radiation%xi_norm/ &
                           min(collinear_radiation%xi_max, &
                               xiScut_used) &
                           *(log(xicut_used)**2 &
                             - log(min(collinear_radiation%xi_max, &
                                       xiScut_used))**2) &
                           /(2.d0*deltaS)
        f_dsc(1) = prefact_deg* &
                   soft_collinear_radiation%jacobian/ &
                   (soft_collinear_shat/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(2) = prefact_deg_sxi* &
                   soft_collinear_radiation%jacobian/ &
                   (soft_collinear_shat/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(3) = prefact_deg_slxi* &
                   soft_collinear_radiation%jacobian/ &
                   (soft_collinear_shat/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(4) = (prefact_deg + prefact_deg_sxi)* &
                   soft_collinear_radiation%jacobian/ &
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
    density_f_s = f_s
    density_f_c = f_c
    density_f_dc = f_dc
    density_f_sc = f_sc
    density_f_dsc = f_dsc
    if (spin_density_fks_collection_enabled()) then
      block_vegas = active_block_vegas_weight(vegas_wgt)
      call active_block_measure(real_event, .true., real_block_measure)
      density_f_r = prefact*real_block_measure* &
           fkssymmetryfactor*block_vegas
      if (.not. nocntevents) then
        call active_block_measure( &
             soft_counterevent, .true., soft_block_measure)
        density_f_s = (prefact + prefact_cnt_ssc)*soft_block_measure* &
             fkssymmetryfactor*block_vegas
        if (fks_sister_is_massless()) then
          call active_block_measure( &
               collinear_counterevent, .true., collinear_block_measure)
          call active_block_measure(soft_collinear_counterevent, .true., &
                                    soft_collinear_block_measure)
          density_f_c = (prefact_c + prefact_coll)* &
               collinear_block_measure*fkssymmetryfactor*block_vegas
          density_f_dc = collinear_block_measure*prefact_deg/ &
               (collinear_shat/(32*pi**2))* &
               fkssymmetryfactorDeg*block_vegas
          density_f_sc = (prefact_c + prefact_coll + &
               prefact_cnt_ssc_c + prefact_coll_c)* &
               soft_collinear_block_measure*fkssymmetryfactorDeg* &
               block_vegas
          density_f_dsc(1) = prefact_deg* &
               soft_collinear_block_measure/ &
               (soft_collinear_shat/(32*pi**2))* &
               fkssymmetryfactorDeg*block_vegas
          density_f_dsc(2) = prefact_deg_sxi* &
               soft_collinear_block_measure/ &
               (soft_collinear_shat/(32*pi**2))* &
               fkssymmetryfactorDeg*block_vegas
          density_f_dsc(3) = prefact_deg_slxi* &
               soft_collinear_block_measure/ &
               (soft_collinear_shat/(32*pi**2))* &
               fkssymmetryfactorDeg*block_vegas
          density_f_dsc(4) = (prefact_deg + prefact_deg_sxi)* &
               soft_collinear_block_measure/ &
               (soft_collinear_shat/(32*pi**2))* &
               fkssymmetryfactorDeg*block_vegas
        end if
      end if
    end if
    call cpu_time(tAfter)
    tf_all = tf_all + (tAfter - tBefore)
    return
  end subroutine compute_prefactors_n1body


end module fks_contributions_module

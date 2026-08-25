module ewsudakov_dummy_module
  use kin_functions_module, only: sumdot => sumdot_impl
  implicit none
  private

  logical :: driver_initialized = .false.
  logical :: defaults_initialized = .false.

  public :: ewsudakov_f77_impl
  public :: get_lo2_orders_impl
  public :: initialize_ewsudakov_defaults
  public :: sdk_get_invariants_impl

  interface
    subroutine init_process_dimensions_bridge()
      implicit none
    end subroutine init_process_dimensions_bridge

    subroutine init_born_dimensions_bridge()
      implicit none
    end subroutine init_born_dimensions_bridge

    subroutine init_fks_metadata_bridge()
      implicit none
    end subroutine init_fks_metadata_bridge

    subroutine setpara(param_card_name)
      implicit none
      character(len=*), intent(in) :: param_card_name
    end subroutine setpara

    subroutine ewsud_dummy_set_coupling(gstr)
      implicit none
      double precision, intent(in) :: gstr
    end subroutine ewsud_dummy_set_coupling

    subroutine ewsud_dummy_evaluate(nexternal, &
         p_born, wgt_born, ewsud_lsc, ewsud_ssc, ewsud_xxc, ewsud_par)
      implicit none
      integer, intent(in) :: nexternal
      double precision, intent(in) :: p_born(0:3, nexternal-1)
      double precision, intent(out) :: wgt_born
      complex(kind=kind(0d0)), intent(out) :: ewsud_lsc, ewsud_ssc
      complex(kind=kind(0d0)), intent(out) :: ewsud_xxc, ewsud_par
    end subroutine ewsud_dummy_evaluate
  end interface

contains

  subroutine initialize_ewsudakov_defaults(sud_mod, sud_filter_hel, &
       sud_mc_hel, fav4, s_to_rij, cs_run, rij_ge_mw)
    implicit none
    integer, intent(inout) :: sud_mod
    logical, intent(inout) :: sud_filter_hel, sud_mc_hel, fav4
    logical, intent(inout) :: s_to_rij, cs_run, rij_ge_mw

    if (defaults_initialized) return
    sud_mod = 2
    sud_filter_hel = .true.
    sud_mc_hel = .true.
    fav4 = .false.
    s_to_rij = .true.
    cs_run = .false.
    rij_ge_mw = .true.
    defaults_initialized = .true.
  end subroutine initialize_ewsudakov_defaults

  subroutine get_lo2_orders_impl(nsplitorders, born_orders, qcd_pos, &
       qed_pos, lo2_orders)
    implicit none
    integer, intent(in) :: nsplitorders, qcd_pos, qed_pos
    integer, intent(in) :: born_orders(nsplitorders)
    integer, intent(out) :: lo2_orders(nsplitorders)

    ! This assumes that only one Born contribution is integrated, as the
    ! generated born.f checks.  LO2 trades two QCD powers for two QED powers.
    lo2_orders = born_orders
    lo2_orders(qcd_pos) = lo2_orders(qcd_pos) - 2
    lo2_orders(qed_pos) = lo2_orders(qed_pos) + 2
  end subroutine get_lo2_orders_impl


  subroutine sdk_get_invariants_impl(nlegs, p, iflist, mdl_mw, &
       rij_ge_mw, invariants)
    implicit none
    integer, intent(in) :: nlegs
    double precision, intent(in) :: p(0:3, nlegs), mdl_mw
    integer, intent(in) :: iflist(nlegs)
    logical, intent(in) :: rij_ge_mw
    double precision, intent(out) :: invariants(nlegs, nlegs)
    double precision :: mw_squared
    integer :: i, j

    mw_squared = mdl_mw**2
    do i = 1, nlegs
      do j = i, nlegs
        invariants(i, j) = sumdot(p(:, i), p(:, j), &
             dble(iflist(i) * iflist(j)))
        if (rij_ge_mw .and. abs(invariants(i, j)) < mw_squared) then
          invariants(i, j) = dsign(1d0, invariants(i, j)) * mw_squared
        end if
        invariants(j, i) = invariants(i, j)
      end do
    end do
  end subroutine sdk_get_invariants_impl


  subroutine ewsudakov_f77_impl(nexternal, p_born_in, gstr_in, results, &
       p_born, sud_mod, nfksprocess, sud_mc_hel, s_to_rij, rij_ge_mw)
    implicit none
    integer, intent(in) :: nexternal
    double precision, intent(in) :: p_born_in(0:3, nexternal-1)
    double precision, intent(in) :: gstr_in
    double precision, intent(out) :: results(6)
    double precision, intent(inout) :: p_born(0:3, nexternal-1)
    integer, intent(inout) :: sud_mod, nfksprocess
    logical, intent(inout) :: sud_mc_hel, s_to_rij, rij_ge_mw
    double precision :: wgt_born, wgt_sud

    nfksprocess = 1
    sud_mc_hel = .false.

    if (.not. driver_initialized) then
      call init_process_dimensions_bridge()
      call init_born_dimensions_bridge()
      call init_fks_metadata_bridge()
      call setpara('param_card.dat')
      driver_initialized = .true.
    end if

    call ewsud_dummy_set_coupling(gstr_in)
    p_born = p_born_in

    s_to_rij = .true.
    rij_ge_mw = .true.
    do sud_mod = 0, 1
      call evaluate_current_configuration(nexternal, p_born, wgt_born, &
           wgt_sud)
      results(1) = wgt_born
      results(2 + sud_mod) = wgt_sud
    end do

    ! Alternative prescriptions retained for comparison with the legacy
    ! standalone/Python interface.
    sud_mod = 1
    s_to_rij = .false.
    rij_ge_mw = .true.
    call evaluate_current_configuration(nexternal, p_born, wgt_born, &
         results(4))

    rij_ge_mw = .false.
    call evaluate_current_configuration(nexternal, p_born, wgt_born, &
         results(5))

    s_to_rij = .true.
    call evaluate_current_configuration(nexternal, p_born, wgt_born, &
         results(6))
  end subroutine ewsudakov_f77_impl


  subroutine evaluate_current_configuration(nexternal, p_born, &
       wgt_born, wgt_sud)
    implicit none
    integer, intent(in) :: nexternal
    double precision, intent(in) :: p_born(0:3, nexternal-1)
    double precision, intent(out) :: wgt_born, wgt_sud
    complex(kind=kind(0d0)) :: ewsud_lsc, ewsud_ssc
    complex(kind=kind(0d0)) :: ewsud_xxc, ewsud_par

    call ewsud_dummy_evaluate(nexternal, p_born, &
         wgt_born, ewsud_lsc, ewsud_ssc, ewsud_xxc, ewsud_par)
    wgt_sud = dble(2d0 * (ewsud_lsc + ewsud_ssc + ewsud_xxc + &
         ewsud_par))
  end subroutine evaluate_current_configuration

end module ewsudakov_dummy_module

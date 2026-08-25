module sa_ewsudakov_dummy_module
  use ranmar_module, only: ntuple
  implicit none
  private

  logical :: driver_initialized = .false.

  public :: ewsudakov_f77_impl
  public :: fill_needed_splittings_impl
  public :: ran2_impl

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

    subroutine sa_ewsud_dummy_set_coupling(gstr)
      implicit none
      double precision, intent(in) :: gstr
    end subroutine sa_ewsud_dummy_set_coupling

    subroutine sa_ewsud_dummy_evaluate(nexternal, p_born, &
         wgt_born, ewsud_lsc, ewsud_ssc, ewsud_xxc, ewsud_par)
      implicit none
      integer, intent(in) :: nexternal
      double precision, intent(in) :: p_born(0:3, nexternal-1)
      double precision, intent(out) :: wgt_born
      complex(kind=kind(0d0)), intent(out) :: ewsud_lsc, ewsud_ssc
      complex(kind=kind(0d0)), intent(out) :: ewsud_xxc, ewsud_par
    end subroutine sa_ewsud_dummy_evaluate
  end interface

contains

  subroutine fill_needed_splittings_impl()
    implicit none

    ! The standalone Sudakov build has no FKS splitting table to populate.
  end subroutine fill_needed_splittings_impl


  double precision function ran2_impl()
    implicit none
    double precision :: x, lower_bound, upper_bound
    integer :: integration_channel, ntuple_index

    lower_bound = 0d0
    upper_bound = 1d0
    ntuple_index = 0
    integration_channel = 1
    call ntuple(x, lower_bound, upper_bound, ntuple_index, &
         integration_channel)
    ran2_impl = x
  end function ran2_impl


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

    call sa_ewsud_dummy_set_coupling(gstr_in)
    p_born = p_born_in

    s_to_rij = .true.
    rij_ge_mw = .true.
    do sud_mod = 0, 1
      call evaluate_current_configuration(nexternal, p_born, wgt_born, &
           wgt_sud)
      results(1) = wgt_born
      results(2 + sud_mod) = wgt_sud
    end do

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

    call sa_ewsud_dummy_evaluate(nexternal, p_born, &
         wgt_born, ewsud_lsc, ewsud_ssc, ewsud_xxc, ewsud_par)
    wgt_sud = dble(2d0 * (ewsud_lsc + ewsud_ssc + ewsud_xxc + &
         ewsud_par))
  end subroutine evaluate_current_configuration

end module sa_ewsudakov_dummy_module

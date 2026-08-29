module top_decay_virtual_cdr
  !! Analytic renormalized QCD one-loop x Born interference for t -> b W+
  !! and t -> b l+ nu_l through an off-shell W.
  !!
  !! CDR/HV (loop_sm lhv=1).  In the two-body process the top and W are
  !! stable on-shell external particles.  In the three-body process the top
  !! is on shell and the internal W can use the same fixed-width propagator as
  !! the generated loop_sm amplitude.  No complex masses enter either
  !! calculation.
  !!
  !! Laurent order follows MadLoop: (finite, 1/eps, 1/eps^2).  Specialized
  !! loop-mode IXXXXX/OXXXXX wavefunctions and MadGraph's VXXXXX preserve the
  !! open-helicity phases of the fixed-order spin-density exporter.
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: tdv_dp = real64
  integer, parameter, public :: tdv_finite = 1
  integer, parameter, public :: tdv_single_pole = 2
  integer, parameter, public :: tdv_double_pole = 3
  integer, parameter, public :: tdv_top_minus = 1
  integer, parameter, public :: tdv_top_plus = 2
  integer, parameter, public :: tdv_w_minus = 1
  integer, parameter, public :: tdv_w_zero = 2
  integer, parameter, public :: tdv_w_plus = 3

  real(real64), parameter :: pi = acos(-1.0_real64)
  real(real64), parameter :: cf = 4.0_real64 / 3.0_real64
  complex(real64), parameter :: ci = (0.0_real64, 1.0_real64)

  type, public :: tdv_two_body_kernel
    private
    real(real64) :: mt, mb, mw, asfac
    real(real64) :: coeff(3,6)
    complex(real64) :: gc11
  end type tdv_two_body_kernel

  public :: tdv_prepare_two_body
  public :: tdv_virtual_rho_top
  public :: tdv_virtual_rho_top_prepared
  public :: tdv_virtual_rho_top_w
  public :: tdv_virtual_rho_top_w_prepared
  public :: tdv_virtual_scalar
  public :: tdv_virtual_scalar_prepared
  public :: tdv_virtual_rho_top_3body
  public :: tdv_virtual_scalar_3body

  interface
    subroutine vxxxxx(p, vmass, nhel, nsv, vc)
      import real64
      real(real64), intent(in) :: p(0:3), vmass
      integer, intent(in) :: nhel, nsv
      complex(real64), intent(out) :: vc(8)
    end subroutine vxxxxx
  end interface

contains

  subroutine tdv_prepare_two_body(kernel, mt, mb, mw, mur, alphas, gc11)
    !! Precompute all invariant-dependent data for repeated on-shell calls.
    type(tdv_two_body_kernel), intent(out) :: kernel
    real(real64), intent(in) :: mt, mb, mw, mur, alphas
    complex(real64), intent(in) :: gc11

    kernel%mt = mt
    kernel%mb = mb
    kernel%mw = mw
    kernel%asfac = alphas*cf/(4.0_real64*pi)
    kernel%gc11 = gc11
    call operator_coefficients(mt, mb, mw, mur, kernel%coeff)
  end subroutine tdv_prepare_two_body


  subroutine tdv_virtual_rho_top(pt, pb, pw, mt, mb, mw, mur, alphas, gc11, rho)
    !! Convenience path which prepares invariant data for this call.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    real(real64), intent(in) :: mt, mb, mw, mur, alphas
    complex(real64), intent(in) :: gc11
    complex(real64), intent(out) :: rho(3,2,2)
    type(tdv_two_body_kernel) :: kernel

    call tdv_prepare_two_body(kernel, mt, mb, mw, mur, alphas, gc11)
    call tdv_virtual_rho_top_prepared(pt, pb, pw, kernel, rho)
  end subroutine tdv_virtual_rho_top


  subroutine tdv_virtual_rho_top_prepared(pt, pb, pw, kernel, rho)
    !! Top helicity order is (-1,+1).  W and bottom helicities are traced
    !! immediately, without constructing the joint top-W density.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    type(tdv_two_body_kernel), intent(in) :: kernel
    complex(real64), intent(out) :: rho(3,2,2)
    real(real64) :: diagonal(3,2)
    complex(real64) :: ft(8,2), fb(8,2), fw(4,3)
    complex(real64) :: born(2), loop(3,2), off_diagonal(3)
    complex(real64) :: ptdot, pbdot
    integer :: it, ib, iw, k

    call incoming_spinor_pair(pt, kernel%mt, ft)
    call outgoing_spinors(pb, kernel%mb, fb)
    call vector_wavefunctions(pw, kernel%mw, fw)
    diagonal = 0.0_real64
    off_diagonal = (0.0_real64, 0.0_real64)

    if (kernel%mb == 0.0_real64) then
      do iw = 1, 3
        pbdot = pb(0)*fw(1,iw) - pb(1)*fw(2,iw) &
               - pb(2)*fw(3,iw) - pb(3)*fw(4,iw)
        do it = 1, 2
          call vertex_amplitudes_massless(ft(:,it), fb(:,1), fw(:,iw), &
                                          pbdot, kernel%coeff, kernel%gc11, &
                                          kernel%asfac, &
                                          born(it), loop(:,it))
        end do
        do k = 1, 3
          diagonal(k,1) = diagonal(k,1) &
                        + 2.0_real64*real(loop(k,1)*conjg(born(1)),real64)
          off_diagonal(k) = off_diagonal(k) &
                          + loop(k,1)*conjg(born(2)) &
                          + born(1)*conjg(loop(k,2))
          diagonal(k,2) = diagonal(k,2) &
                        + 2.0_real64*real(loop(k,2)*conjg(born(2)),real64)
        end do
      end do
    else
      do iw = 1, 3
        ptdot = pt(0)*fw(1,iw) - pt(1)*fw(2,iw) &
               - pt(2)*fw(3,iw) - pt(3)*fw(4,iw)
        pbdot = pb(0)*fw(1,iw) - pb(1)*fw(2,iw) &
               - pb(2)*fw(3,iw) - pb(3)*fw(4,iw)
        do ib = 1, 2
          do it = 1, 2
            call vertex_amplitudes_massive(ft(:,it), fb(:,ib), fw(:,iw), &
                                           ptdot, pbdot, kernel%coeff, &
                                           kernel%gc11, kernel%asfac, &
                                           born(it), loop(:,it))
          end do
          do k = 1, 2
            diagonal(k,1) = diagonal(k,1) &
                          + 2.0_real64*real(loop(k,1)*conjg(born(1)),real64)
            off_diagonal(k) = off_diagonal(k) &
                            + loop(k,1)*conjg(born(2)) &
                            + born(1)*conjg(loop(k,2))
            diagonal(k,2) = diagonal(k,2) &
                          + 2.0_real64*real(loop(k,2)*conjg(born(2)),real64)
          end do
        end do
      end do
    end if

    do k = 1, 3
      rho(k,1,1) = cmplx(diagonal(k,1), 0.0_real64, real64)
      rho(k,1,2) = off_diagonal(k)
      rho(k,2,1) = conjg(off_diagonal(k))
      rho(k,2,2) = cmplx(diagonal(k,2), 0.0_real64, real64)
    end do
  end subroutine tdv_virtual_rho_top_prepared


  subroutine tdv_virtual_rho_top_w(pt, pb, pw, mt, mb, mw, mur, alphas, gc11, rho)
    !! Convenience path which prepares invariant data for this call.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    real(real64), intent(in) :: mt, mb, mw, mur, alphas
    complex(real64), intent(in) :: gc11
    complex(real64), intent(out) :: rho(3,6,6)
    type(tdv_two_body_kernel) :: kernel

    call tdv_prepare_two_body(kernel, mt, mb, mw, mur, alphas, gc11)
    call tdv_virtual_rho_top_w_prepared(pt, pb, pw, kernel, rho)
  end subroutine tdv_virtual_rho_top_w


  subroutine tdv_virtual_rho_top_w_prepared(pt, pb, pw, kernel, rho)
    !! Open index = top_index + 2*(W_index-1), exactly as export_spin_density.
    !! Top order is (-1,+1), W order is (-1,0,+1); only bottom is summed.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    type(tdv_two_body_kernel), intent(in) :: kernel
    complex(real64), intent(out) :: rho(3,6,6)
    complex(real64) :: born(2,2,3), loop(3,2,2,3), value
    integer :: a, b, ia, ja, ib, iw, jw, k, nb

    call helicity_amplitudes_prepared(pt, pb, pw, kernel, born, loop)
    nb = merge(1, 2, kernel%mb == 0.0_real64)
    rho = (0.0_real64, 0.0_real64)
    do k = 1, 3
      do iw = 1, 3
        do ia = 1, 2
          a = ia + 2*(iw-1)
          do jw = iw, 3
            do ja = 1, 2
              b = ja + 2*(jw-1)
              if (b < a) cycle
              value = (0.0_real64, 0.0_real64)
              do ib = 1, nb
                value = value + loop(k,ia,ib,iw)*conjg(born(ja,ib,jw)) &
                              + born(ia,ib,iw)*conjg(loop(k,ja,ib,jw))
              end do
              if (a == b) value = cmplx(real(value,real64), 0.0_real64, real64)
              rho(k,a,b) = value
              rho(k,b,a) = conjg(value)
            end do
          end do
        end do
      end do
    end do
  end subroutine tdv_virtual_rho_top_w_prepared


  subroutine tdv_virtual_scalar(pt, pb, pw, mt, mb, mw, mur, alphas, gc11, born_me, virtual)
    !! Convenience path which prepares invariant data for this call.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    real(real64), intent(in) :: mt, mb, mw, mur, alphas
    complex(real64), intent(in) :: gc11
    real(real64), intent(out) :: born_me, virtual(3)
    type(tdv_two_body_kernel) :: kernel

    call tdv_prepare_two_body(kernel, mt, mb, mw, mur, alphas, gc11)
    call tdv_virtual_scalar_prepared(pt, pb, pw, kernel, born_me, virtual)
  end subroutine tdv_virtual_scalar


  subroutine tdv_virtual_scalar_prepared(pt, pb, pw, kernel, born_me, virtual)
    !! Spin/color averaged scalar convention returned by standalone MadLoop.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    type(tdv_two_body_kernel), intent(in) :: kernel
    real(real64), intent(out) :: born_me, virtual(3)
    complex(real64) :: ft(8,2), fb(8,2), fw(4,3)
    complex(real64) :: born(2), loop(3,2), ptdot, pbdot
    integer :: it, ib, iw, k

    call incoming_spinor_pair(pt, kernel%mt, ft)
    call outgoing_spinors(pb, kernel%mb, fb)
    call vector_wavefunctions(pw, kernel%mw, fw)
    born_me = 0.0_real64
    virtual = 0.0_real64

    if (kernel%mb == 0.0_real64) then
      do iw = 1, 3
        pbdot = pb(0)*fw(1,iw) - pb(1)*fw(2,iw) &
               - pb(2)*fw(3,iw) - pb(3)*fw(4,iw)
        do it = 1, 2
          call vertex_amplitudes_massless(ft(:,it), fb(:,1), fw(:,iw), &
                                          pbdot, kernel%coeff, kernel%gc11, &
                                          kernel%asfac, &
                                          born(it), loop(:,it))
          born_me = born_me + abs(born(it))**2
          do k = 1, 3
            virtual(k) = virtual(k) &
                       + 2.0_real64*real(loop(k,it)*conjg(born(it)),real64)
          end do
        end do
      end do
    else
      do iw = 1, 3
        ptdot = pt(0)*fw(1,iw) - pt(1)*fw(2,iw) &
               - pt(2)*fw(3,iw) - pt(3)*fw(4,iw)
        pbdot = pb(0)*fw(1,iw) - pb(1)*fw(2,iw) &
               - pb(2)*fw(3,iw) - pb(3)*fw(4,iw)
        do ib = 1, 2
          do it = 1, 2
            call vertex_amplitudes_massive(ft(:,it), fb(:,ib), fw(:,iw), &
                                           ptdot, pbdot, kernel%coeff, &
                                           kernel%gc11, kernel%asfac, &
                                           born(it), loop(:,it))
            born_me = born_me + abs(born(it))**2
            do k = 1, 2
              virtual(k) = virtual(k) &
                         + 2.0_real64*real(loop(k,it)*conjg(born(it)),real64)
            end do
          end do
        end do
      end do
    end if
    born_me = 0.5_real64*born_me
    virtual = 0.5_real64*virtual
  end subroutine tdv_virtual_scalar_prepared


  subroutine tdv_virtual_rho_top_3body(pt, pb, pl, pnu, mt, mb, mw, &
                                      mur, alphas, gc11, rho, ww)
    !! One massless flavor t -> b l+ nu_l.  Top order is (-1,+1);
    !! bottom is traced immediately and the unique nonzero lepton helicities
    !! are implicit.  The optional fixed W width defaults to zero for
    !! compatibility with the standalone topdecay interface.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pl(0:3), pnu(0:3)
    real(real64), intent(in) :: mt, mb, mw, mur, alphas
    complex(real64), intent(in) :: gc11
    complex(real64), intent(out) :: rho(3,2,2)
    real(real64), intent(in), optional :: ww
    real(real64) :: coeff(3,6), asfac, q(0:3), q2, qmass, wwidth
    real(real64) :: diagonal(3,2)
    complex(real64) :: ft(8,2), fb(8,2), fw(8)
    complex(real64) :: born(2), loop(3,2), off_diagonal(3)
    complex(real64) :: ptdot, pbdot
    integer :: it, ib, k

    wwidth = 0.0_real64
    if (present(ww)) wwidth = ww
    q = pl + pnu
    q2 = q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)
    qmass = sqrt(q2)
    call operator_coefficients(mt, mb, qmass, mur, coeff)
    asfac = alphas*cf/(4.0_real64*pi)
    call incoming_spinor_pair(pt, mt, ft)
    call outgoing_spinors(pb, mb, fb)
    call leptonic_w_current(pl, pnu, gc11, mw, wwidth, q2, fw)
    pbdot = pb(0)*fw(5) - pb(1)*fw(6) - pb(2)*fw(7) - pb(3)*fw(8)
    diagonal = 0.0_real64
    off_diagonal = (0.0_real64, 0.0_real64)

    if (mb == 0.0_real64) then
      do it = 1, 2
        call vertex_amplitudes_massless(ft(:,it), fb(:,1), fw(5:8), pbdot, &
                                        coeff, gc11, asfac, born(it), &
                                        loop(:,it))
      end do
      do k = 1, 3
        diagonal(k,1) = 2.0_real64*real(loop(k,1)*conjg(born(1)),real64)
        off_diagonal(k) = loop(k,1)*conjg(born(2)) &
                        + born(1)*conjg(loop(k,2))
        diagonal(k,2) = 2.0_real64*real(loop(k,2)*conjg(born(2)),real64)
      end do
    else
      ptdot = pt(0)*fw(5) - pt(1)*fw(6) - pt(2)*fw(7) - pt(3)*fw(8)
      do ib = 1, 2
        do it = 1, 2
          call vertex_amplitudes_massive(ft(:,it), fb(:,ib), fw(5:8), &
                                         ptdot, pbdot, coeff, gc11, &
                                         asfac, born(it), loop(:,it))
        end do
        do k = 1, 2
          diagonal(k,1) = diagonal(k,1) &
                        + 2.0_real64*real(loop(k,1)*conjg(born(1)),real64)
          off_diagonal(k) = off_diagonal(k) &
                          + loop(k,1)*conjg(born(2)) &
                          + born(1)*conjg(loop(k,2))
          diagonal(k,2) = diagonal(k,2) &
                        + 2.0_real64*real(loop(k,2)*conjg(born(2)),real64)
        end do
      end do
    end if

    do k = 1, 3
      rho(k,1,1) = cmplx(diagonal(k,1), 0.0_real64, real64)
      rho(k,1,2) = off_diagonal(k)
      rho(k,2,1) = conjg(off_diagonal(k))
      rho(k,2,2) = cmplx(diagonal(k,2), 0.0_real64, real64)
    end do
  end subroutine tdv_virtual_rho_top_3body


  subroutine tdv_virtual_scalar_3body(pt, pb, pl, pnu, mt, mb, mw, &
                                     mur, alphas, gc11, born_me, virtual, ww)
    !! Spin/color averaged MadLoop scalar for one massless lepton flavor.  The
    !! optional fixed W width defaults to zero.
    real(real64), intent(in) :: pt(0:3), pb(0:3), pl(0:3), pnu(0:3)
    real(real64), intent(in) :: mt, mb, mw, mur, alphas
    complex(real64), intent(in) :: gc11
    real(real64), intent(out) :: born_me, virtual(3)
    real(real64), intent(in), optional :: ww
    complex(real64) :: born(2,2), loop(3,2,2)
    real(real64) :: wwidth
    integer :: it, k, nk

    wwidth = 0.0_real64
    if (present(ww)) wwidth = ww
    call helicity_amplitudes_3body(pt, pb, pl, pnu, mt, mb, mw, &
                                  wwidth, mur, alphas, gc11, born, loop)
    nk = merge(3, 2, mb == 0.0_real64)
    born_me = 0.0_real64
    virtual = 0.0_real64
    do it = 1, 2
      born_me = born_me + abs(born(it,1))**2
      do k = 1, nk
        virtual(k) = virtual(k) + 2.0_real64*real(loop(k,it,1)*conjg(born(it,1)), real64)
      end do
    end do
    if (mb /= 0.0_real64) then
      do it = 1, 2
        born_me = born_me + abs(born(it,2))**2
        do k = 1, 2
          virtual(k) = virtual(k) + 2.0_real64*real(loop(k,it,2)*conjg(born(it,2)), real64)
        end do
      end do
    end if
    born_me = 0.5_real64*born_me
    virtual = 0.5_real64*virtual
  end subroutine tdv_virtual_scalar_3body


  subroutine helicity_amplitudes_prepared(pt, pb, pw, kernel, born, loop)
    real(real64), intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    type(tdv_two_body_kernel), intent(in) :: kernel
    complex(real64), intent(out) :: born(2,2,3), loop(3,2,2,3)
    complex(real64) :: ft(8,2), fb(8,2), fw(8,3)
    complex(real64) :: amp_loop(3), ptdot, pbdot
    logical :: massless_b
    integer, parameter :: hw(3) = [-1, 0, 1]
    integer :: it, ib, iw

    massless_b = (kernel%mb == 0.0_real64)
    call incoming_spinor_pair(pt, kernel%mt, ft)
    call outgoing_spinors(pb, kernel%mb, fb)
    do iw = 1, 3
      call vxxxxx(pw, kernel%mw, hw(iw), 1, fw(:,iw))
    end do

    born = (0.0_real64, 0.0_real64)
    loop = (0.0_real64, 0.0_real64)
    do iw = 1, 3
      pbdot = pb(0)*fw(5,iw) - pb(1)*fw(6,iw) - pb(2)*fw(7,iw) - pb(3)*fw(8,iw)
      if (massless_b) then
        do it = 1, 2
          call vertex_amplitudes_massless(ft(:,it), fb(:,1), fw(5:8,iw), &
                                          pbdot, kernel%coeff, kernel%gc11, &
                                          kernel%asfac, &
                                          born(it,1,iw), amp_loop)
          loop(:,it,1,iw) = amp_loop
        end do
      else
        ptdot = pt(0)*fw(5,iw) - pt(1)*fw(6,iw) - pt(2)*fw(7,iw) - pt(3)*fw(8,iw)
        do ib = 1, 2
          do it = 1, 2
            call vertex_amplitudes_massive(ft(:,it), fb(:,ib), fw(5:8,iw), &
                                           ptdot, pbdot, kernel%coeff, &
                                           kernel%gc11, kernel%asfac, &
                                           born(it,ib,iw), amp_loop)
            loop(:,it,ib,iw) = amp_loop
          end do
        end do
      end if
    end do
  end subroutine helicity_amplitudes_prepared


  subroutine helicity_amplitudes_3body(pt, pb, pl, pnu, mt, mb, mw, &
                                      ww, mur, alphas, gc11, born, loop)
    real(real64), intent(in) :: pt(0:3), pb(0:3), pl(0:3), pnu(0:3)
    real(real64), intent(in) :: mt, mb, mw, ww, mur, alphas
    complex(real64), intent(in) :: gc11
    complex(real64), intent(out) :: born(2,2), loop(3,2,2)
    real(real64) :: coeff(3,6), asfac, q(0:3), q2, qmass
    complex(real64) :: ft(8,2), fb(8,2), fw(8), amp_loop(3)
    complex(real64) :: ptdot, pbdot
    logical :: massless_b
    integer :: it, ib

    massless_b = (mb == 0.0_real64)
    q = pl + pnu
    q2 = q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)
    qmass = sqrt(q2)
    call operator_coefficients(mt, mb, qmass, mur, coeff)
    asfac = alphas*cf/(4.0_real64*pi)
    call incoming_spinor_pair(pt, mt, ft)
    call outgoing_spinors(pb, mb, fb)
    call leptonic_w_current(pl, pnu, gc11, mw, ww, q2, fw)

    born = (0.0_real64, 0.0_real64)
    loop = (0.0_real64, 0.0_real64)
    pbdot = pb(0)*fw(5) - pb(1)*fw(6) - pb(2)*fw(7) - pb(3)*fw(8)
    if (massless_b) then
      do it = 1, 2
        call vertex_amplitudes_massless(ft(:,it), fb(:,1), fw(5:8), pbdot, &
                                        coeff, gc11, asfac, born(it,1), amp_loop)
        loop(:,it,1) = amp_loop
      end do
    else
      ptdot = pt(0)*fw(5) - pt(1)*fw(6) - pt(2)*fw(7) - pt(3)*fw(8)
      do ib = 1, 2
        do it = 1, 2
          call vertex_amplitudes_massive(ft(:,it), fb(:,ib), fw(5:8), ptdot, &
                                         pbdot, coeff, gc11, asfac, &
                                         born(it,ib), amp_loop)
          loop(:,it,ib) = amp_loop
        end do
      end do
    end if
  end subroutine helicity_amplitudes_3body


  subroutine vertex_amplitudes_massless(ft, fb, fv, pbdot, coeff, gc11, asfac, born, loop)
    real(real64), intent(in) :: coeff(3,6), asfac
    complex(real64), intent(in) :: ft(8), fb(8), fv(4), gc11
    complex(real64), intent(in) :: pbdot
    complex(real64), intent(out) :: born, loop(3)
    complex(real64) :: gl, sr, prefactor

    gl = ft(5)*(fb(7)*(fv(1)+fv(4)) + fb(8)*(fv(2)+ci*fv(3))) &
       + ft(6)*(fb(7)*(fv(2)-ci*fv(3)) + fb(8)*(fv(1)-fv(4)))
    sr = ft(7)*fb(7) + ft(8)*fb(8)
    born = -ci*gc11*gl
    prefactor = -ci*gc11*asfac
    loop(1) = prefactor*(coeff(1,1)*gl + coeff(1,6)*pbdot*sr)
    loop(2) = prefactor*coeff(2,1)*gl
    loop(3) = prefactor*coeff(3,1)*gl
  end subroutine vertex_amplitudes_massless


  subroutine vertex_amplitudes_massive(ft, fb, fv, ptdot, pbdot, coeff, gc11, asfac, born, loop)
    real(real64), intent(in) :: coeff(3,6), asfac
    complex(real64), intent(in) :: ft(8), fb(8), fv(4), gc11, ptdot, pbdot
    complex(real64), intent(out) :: born, loop(3)
    complex(real64) :: gl, gr, sl, sr, current, prefactor

    gl = ft(5)*(fb(7)*(fv(1)+fv(4)) + fb(8)*(fv(2)+ci*fv(3))) &
       + ft(6)*(fb(7)*(fv(2)-ci*fv(3)) + fb(8)*(fv(1)-fv(4)))
    gr = ft(7)*(fb(5)*(fv(1)-fv(4)) - fb(6)*(fv(2)+ci*fv(3))) &
       + ft(8)*(fb(5)*(-fv(2)+ci*fv(3)) + fb(6)*(fv(1)+fv(4)))
    sl = ft(5)*fb(5) + ft(6)*fb(6)
    sr = ft(7)*fb(7) + ft(8)*fb(8)
    born = -ci*gc11*gl
    prefactor = -ci*gc11*asfac
    current = coeff(1,1)*gl + coeff(1,2)*gr &
            + (coeff(1,3)*ptdot + coeff(1,5)*pbdot)*sl &
            + (coeff(1,4)*ptdot + coeff(1,6)*pbdot)*sr
    loop(1) = prefactor*current
    loop(2) = prefactor*coeff(2,1)*gl
    loop(3) = (0.0_real64, 0.0_real64)
  end subroutine vertex_amplitudes_massive


  subroutine incoming_spinor_pair(p, mass, fi)
    !! Specialized loop-mode IXXXXX for a positive-mass incoming fermion,
    !! returning helicities (-1,+1) while sharing their kinematic square roots.
    real(real64), intent(in) :: p(0:3), mass
    complex(real64), intent(out) :: fi(8,2)
    real(real64) :: pp, pp3, omega1, omega2, root, denom
    complex(real64) :: chi1, chim, chip

    fi = (0.0_real64, 0.0_real64)
    pp = min(p(0), sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3)))
    if (pp == 0.0_real64) then
      root = sqrt(abs(mass))
      fi(6,1) = root
      fi(8,1) = root
      fi(5,2) = root
      fi(7,2) = root
      return
    end if

    omega1 = sqrt(p(0)+pp)
    omega2 = mass/omega1
    pp3 = max(pp+p(3), 0.0_real64)
    chi1 = cmplx(sqrt(0.5_real64*pp3/pp), 0.0_real64, real64)
    if (pp3 == 0.0_real64) then
      chim = (1.0_real64, 0.0_real64)
      chip = (-1.0_real64, 0.0_real64)
    else
      denom = sqrt(2.0_real64*pp*pp3)
      chim = cmplx(-p(1), p(2), real64)/denom
      chip = cmplx( p(1), p(2), real64)/denom
    end if
    fi(5,1) = omega1*chim
    fi(6,1) = omega1*chi1
    fi(7,1) = omega2*chim
    fi(8,1) = omega2*chi1
    fi(5,2) = omega2*chi1
    fi(6,2) = omega2*chip
    fi(7,2) = omega1*chi1
    fi(8,2) = omega1*chip
  end subroutine incoming_spinor_pair


  subroutine outgoing_spinors(p, mass, fo)
    !! Specialized loop-mode OXXXXX.  For mass=0 only the nonzero negative
    !! helicity is constructed; otherwise both (-1,+1) share their roots.
    real(real64), intent(in) :: p(0:3), mass
    complex(real64), intent(out) :: fo(8,2)
    real(real64) :: pp, pp3, sqp0p3, omega1, omega2, root, denom
    complex(real64) :: chi1, chim, chip

    fo = (0.0_real64, 0.0_real64)
    if (mass == 0.0_real64) then
      if (p(1) == 0.0_real64 .and. p(2) == 0.0_real64 .and. p(3) < 0.0_real64) then
        sqp0p3 = 0.0_real64
      else
        sqp0p3 = sqrt(max(p(0)+p(3), 0.0_real64))
      end if
      if (sqp0p3 == 0.0_real64) then
        fo(7,1) = sqrt(2.0_real64*p(0))
      else
        fo(7,1) = cmplx(-p(1), -p(2), real64)/sqp0p3
        fo(8,1) = sqp0p3
      end if
      return
    end if

    pp = min(p(0), sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3)))
    if (pp == 0.0_real64) then
      root = sqrt(abs(mass))
      fo(6,1) = root
      fo(8,1) = root
      fo(5,2) = root
      fo(7,2) = root
      return
    end if

    omega1 = sqrt(p(0)+pp)
    omega2 = mass/omega1
    pp3 = max(pp+p(3), 0.0_real64)
    chi1 = cmplx(sqrt(0.5_real64*pp3/pp), 0.0_real64, real64)
    if (pp3 == 0.0_real64) then
      chim = (1.0_real64, 0.0_real64)
      chip = (-1.0_real64, 0.0_real64)
    else
      denom = sqrt(2.0_real64*pp*pp3)
      chim = cmplx(-p(1), -p(2), real64)/denom
      chip = cmplx( p(1), -p(2), real64)/denom
    end if
    fo(5,1) = omega2*chim
    fo(6,1) = omega2*chi1
    fo(7,1) = omega1*chim
    fo(8,1) = omega1*chi1
    fo(5,2) = omega1*chi1
    fo(6,2) = omega1*chip
    fo(7,2) = omega2*chi1
    fo(8,2) = omega2*chip
  end subroutine outgoing_spinors


  subroutine vector_wavefunctions(p, mass, fv)
    !! Specialized loop-mode VXXXXX for a massive final-state vector.  All
    !! helicities (-1,0,+1) share the momentum norms and square roots.
    real(real64), intent(in) :: p(0:3), mass
    complex(real64), intent(out) :: fv(4,3)
    real(real64), parameter :: sqh = sqrt(0.5_real64)
    real(real64) :: pt2, pp, pt, emp, pzpt, transverse_sign

    fv = (0.0_real64, 0.0_real64)
    pt2 = p(1)*p(1) + p(2)*p(2)
    pp = min(p(0), sqrt(pt2+p(3)*p(3)))
    pt = min(pp, sqrt(pt2))

    if (pp == 0.0_real64) then
      fv(2,1) = sqh
      fv(3,1) = ci*sqh
      fv(4,2) = 1.0_real64
      fv(2,3) = -sqh
      fv(3,3) = ci*sqh
      return
    end if

    emp = p(0)/(mass*pp)
    fv(1,2) = pp/mass
    fv(2,2) = p(1)*emp
    fv(3,2) = p(2)*emp
    fv(4,2) = p(3)*emp
    fv(4,1) = -pt/pp*sqh
    fv(4,3) =  pt/pp*sqh
    if (pt /= 0.0_real64) then
      pzpt = p(3)/(pp*pt)*sqh
      fv(2,1) = cmplx( p(1)*pzpt, -p(2)/pt*sqh, real64)
      fv(3,1) = cmplx( p(2)*pzpt,  p(1)/pt*sqh, real64)
      fv(2,3) = cmplx(-p(1)*pzpt, -p(2)/pt*sqh, real64)
      fv(3,3) = cmplx(-p(2)*pzpt,  p(1)/pt*sqh, real64)
    else
      transverse_sign = sign(sqh, p(3))
      fv(2,1) = sqh
      fv(3,1) = ci*transverse_sign
      fv(2,3) = -sqh
      fv(3,3) = ci*transverse_sign
    end if
  end subroutine vector_wavefunctions


  subroutine leptonic_w_current(pl, pnu, gc11, mw, ww, q2, fw)
    !! ALOHA FFV2_3 with its fixed-width propagator.  The q_mu*q_nu term
    !! vanishes analytically against the massless lepton current.
    real(real64), intent(in) :: pl(0:3), pnu(0:3)
    complex(real64), intent(in) :: gc11
    real(real64), intent(in) :: mw, ww, q2
    complex(real64), intent(out) :: fw(8)
    real(real64) :: sl, sn
    complex(real64) :: l5, l6, n7, n8, denom

    if (pl(1) == 0.0_real64 .and. pl(2) == 0.0_real64 .and. pl(3) < 0.0_real64) then
      sl = 0.0_real64
    else
      sl = -sqrt(max(pl(0)+pl(3), 0.0_real64))
    end if
    if (sl == 0.0_real64) then
      l5 = cmplx(-sqrt(2.0_real64*pl(0)), 0.0_real64, real64)
      l6 = (0.0_real64, 0.0_real64)
    else
      l5 = cmplx(-pl(1), pl(2), real64)/sl
      l6 = cmplx(sl, 0.0_real64, real64)
    end if

    if (pnu(1) == 0.0_real64 .and. pnu(2) == 0.0_real64 .and. pnu(3) < 0.0_real64) then
      sn = 0.0_real64
    else
      sn = sqrt(max(pnu(0)+pnu(3), 0.0_real64))
    end if
    if (sn == 0.0_real64) then
      n7 = cmplx(sqrt(2.0_real64*pnu(0)), 0.0_real64, real64)
      n8 = (0.0_real64, 0.0_real64)
    else
      n7 = cmplx(-pnu(1), -pnu(2), real64)/sn
      n8 = cmplx(sn, 0.0_real64, real64)
    end if

    fw(1:4) = cmplx(pl+pnu, 0.0_real64, real64)
    denom = gc11/cmplx(q2-mw*mw, mw*ww, real64)
    fw(5) = -ci*denom*( n7*l5 + n8*l6)
    fw(6) = -ci*denom*(-n8*l5 - n7*l6)
    fw(7) = -ci*denom*(-ci*n8*l5 + ci*n7*l6)
    fw(8) = -ci*denom*(-n7*l5 + n8*l6)
  end subroutine leptonic_w_current


  subroutine operator_coefficients(mt, mb, mw, mur, coeff)
    !! Columns are gamma_L,gamma_R,pt_L,pt_R,pb_L,pb_R.
    !! The mb=0 current is Campbell--Ellis--Tramontano, eta=1 (CDR/HV).
    real(real64), intent(in) :: mt, mb, mw, mur
    real(real64), intent(out) :: coeff(3,6)
    real(real64) :: ff(6), pole_f1, z, lz, lmu, c0, c1

    coeff = 0.0_real64
    if (abs(mb) < tiny(1.0_real64)) then
      z = (mw/mt)**2
      lz = log_one_minus(z)
      lmu = 2.0_real64*log(mur/mt)
      c0 = -6.0_real64 - pi*pi/6.0_real64 - 2.0_real64*li2_01(z) &
           + 3.0_real64*lz - 2.0_real64*lz*lz - lz/z
      c1 = 2.0_real64*lz/z
      coeff(1,1) = c0 + lmu*(-2.5_real64 + 2.0_real64*lz) - 0.5_real64*lmu*lmu
      coeff(1,6) = c1/mt
      coeff(2,1) = -2.5_real64 + 2.0_real64*lz - lmu
      coeff(3,1) = -1.0_real64
      return
    end if

    call massive_form_factors(mt, mb, mw, mur, ff, pole_f1)
    coeff(1,1) = 0.5_real64*(ff(1) + ff(4))
    coeff(1,2) = 0.5_real64*(ff(1) - ff(4))
    coeff(1,3) = 0.5_real64*(ff(2) + ff(5))
    coeff(1,4) = 0.5_real64*(ff(2) - ff(5))
    coeff(1,5) = 0.5_real64*(ff(3) + ff(6))
    coeff(1,6) = 0.5_real64*(ff(3) - ff(6))
    coeff(2,1) = pole_f1
  end subroutine operator_coefficients


  subroutine massive_form_factors(mt, mb, mw, mur, ff, pole_f1)
    !! Fischer--Groote--Korner--Mauser Appendix-C form factors, with
    !! alpha_s*C_F/(4*pi) stripped.  ff=(dF1V,F2V,F3V,dF1A,F2A,F3A).
    !! log(mg^4/(mb^2*mt^2)) -> 2/eps+log(mur^4/(mb^2*mt^2)).
    real(real64), intent(in) :: mt, mb, mw, mur
    real(real64), intent(out) :: ff(6), pole_f1
    real(real64) :: root, w1, wmu, lw1, lwmu, lw, reglog, rcoef
    real(real64) :: mt2, mb2, q2, lmt, lmb, logmass
    real(real64) :: bracket, f1_common

    mt2 = mt*mt
    mb2 = mb*mb
    q2 = mw*mw
    root = kallen_root(mt, mb, mw)
    call w_values(mt, mb, mw, w1, wmu)
    lw1 = log(w1)
    lwmu = log(wmu)
    lw = lw1 + lwmu
    lmt = log(mt)
    lmb = log(mb)
    reglog = 2.0_real64*(2.0_real64*log(mur) - lmt - lmb)
    rcoef = -1.0_real64 - (mt2 + mb2 - q2)*lw/(2.0_real64*root)

    bracket = 2.0_real64*li2_01(1.0_real64-w1*w1) &
            - 2.0_real64*li2_01(1.0_real64-w1/wmu) &
            + 0.5_real64*reglog*lw &
            + (3.0_real64*lw1-lwmu) &
              *(lwmu + log_one_minus(w1*w1) - log(wmu-w1))
    logmass = 2.0_real64*(lmb-lmt)
    f1_common = -(mt2+mb2-q2)*bracket/root - reglog &
                - (mt2-mb2)*logmass/(2.0_real64*q2) - 4.0_real64 &
                + lw*root/(2.0_real64*q2)
    ff(1) = f1_common - lw*((mt+mb)**2-q2)/root
    ff(4) = f1_common - lw*((mt-mb)**2-q2)/root

    ff(2) = f2_vector_logs(mt, mb, mw, root, logmass, lw)
    ff(3) = f2_vector_logs(mb, mt, mw, root, -logmass, lw)
    ff(5) = f2_vector_logs(-mt, mb, mw, root, logmass, lw)
    ff(6) = f2_vector_logs(mb, -mt, mw, root, -logmass, lw)
    pole_f1 = 2.0_real64*rcoef
  end subroutine massive_form_factors


  pure real(real64) function f2_vector_logs(ma, mb, q, root, logm, logw) result(value)
    real(real64), intent(in) :: ma, mb, q, root, logm, logw
    real(real64) :: dm, termm, termw

    dm = ma - mb
    termm = (ma+2.0_real64*mb)/dm - (ma*ma-mb*mb)/(q*q)
    termw = root/(q*q) - mb/dm*(q*q + dm*(3.0_real64*ma+mb))/root
    value = dm/(q*q)*(2.0_real64 - termm*logm - termw*logw)
  end function f2_vector_logs


  pure real(real64) function kallen_root(a, b, q) result(root)
    real(real64), intent(in) :: a, b, q
    root = sqrt((a*a-(b+q)**2)*(a*a-(b-q)**2))
  end function kallen_root


  pure subroutine w_values(a, b, q, w1, wmu)
    real(real64), intent(in) :: a, b, q
    real(real64), intent(out) :: w1, wmu
    real(real64) :: aa, bb, root, anum, d0, dplus, dminus, product

    aa = abs(a)
    bb = abs(b)
    root = kallen_root(aa, bb, q)
    anum = 4.0_real64*aa*aa*bb*bb/(aa*aa+bb*bb-q*q+root)
    d0 = aa*aa + q*q - bb*bb
    dplus = d0 + root
    dminus = d0 - root
    product = 4.0_real64*aa*aa*q*q
    if (abs(dplus) < 0.5_real64*(abs(d0)+root)) dplus = product/dminus
    if (abs(dminus) < 0.5_real64*(abs(d0)+root)) dminus = product/dplus
    w1 = q/bb*anum/dplus
    wmu = q/bb*anum/dminus
  end subroutine w_values


  pure real(real64) function li2_01(xin) result(value)
    !! Real Li_2 on [0,1], the only range needed by physical t -> b W.
    !! The Bernoulli expansion in y=-log(1-x) needs only nine FMA-friendly
    !! terms after mapping x to [0,1/2].
    real(real64), intent(in) :: xin
    real(real64) :: x, t, y, y2, series

    x = min(1.0_real64, max(0.0_real64, xin))
    if (x <= 0.0_real64) then
      value = 0.0_real64
      return
    else if (x >= 1.0_real64) then
      value = pi*pi/6.0_real64
      return
    end if
    t = merge(1.0_real64-x, x, x > 0.5_real64)
    y = -log_one_minus(t)
    y2 = y*y
    series = 4.5189800296199182e-16_real64
    series = -1.9939295860721076e-14_real64 + y2*series
    series = 8.9216910204564526e-13_real64 + y2*series
    series = -4.0647616451442255e-11_real64 + y2*series
    series = 1.8978869988970999e-9_real64 + y2*series
    series = -9.1857730746619636e-8_real64 + y2*series
    series = 4.7241118669690098e-6_real64 + y2*series
    series = -2.7777777777777778e-4_real64 + y2*series
    series = 2.7777777777777778e-2_real64 + y2*series
    value = y + y2*(-0.25_real64 + y*series)
    if (x > 0.5_real64) then
      value = pi*pi/6.0_real64 - log(x)*log_one_minus(x) - value
    end if
  end function li2_01


  pure real(real64) function log_one_minus(x) result(value)
    real(real64), intent(in) :: x
    real(real64) :: term, sum, add
    integer :: n

    if (abs(x) >= 1.0e-4_real64) then
      value = log(1.0_real64-x)
      return
    end if
    term = x
    sum = x
    do n = 2, 20
      term = term*x
      add = term/real(n,real64)
      sum = sum + add
      if (abs(add) <= epsilon(1.0_real64)*abs(sum)) exit
    end do
    value = -sum
  end function log_one_minus

end module top_decay_virtual_cdr

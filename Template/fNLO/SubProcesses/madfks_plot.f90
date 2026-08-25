! Wrapper routines for the fixed-order analyses.
module madfks_plot_module
  use boostwdir2_module, only: boostwdir2
  use extra_weights
  use mint_module, only: itmax, ncalls0
  use process_dimensions, only: nexternal, nincoming, &
       validate_process_dimensions
  use run_state, only: do_rwgt_pdf, do_rwgt_scale, pineappl
  implicit none
  private

  public :: initplot_impl, topout_impl, outfun_impl

  integer, parameter :: pine_reset_action = 1
  integer, parameter :: pine_norm_action = 2
  integer, parameter :: pine_event_action = 3

  logical :: useitmax
  common /cuseitmax/ useitmax

  double precision :: xsecScale_acc(maxscales, maxscales, maxdynscales)
  double precision :: xsecPDFr_acc(0:maxPDFs, maxPDFsets)
  common /scale_pdf_print/ xsecScale_acc, xsecPDFr_acc

  interface
    subroutine analysis_begin(nwgt, weights_info)
      integer, intent(in) :: nwgt
      character(len=*), intent(in) :: weights_info(*)
    end subroutine analysis_begin

    subroutine analysis_end(xnorm)
      double precision, intent(in) :: xnorm
    end subroutine analysis_end

    subroutine analysis_fill(p, istatus, ipdg, wgts, ibody)
      use process_dimensions, only: nexternal
      double precision, intent(in) :: p(0:4, nexternal)
      integer, intent(in) :: istatus(nexternal), ipdg(nexternal)
      double precision, intent(in) :: wgts(*)
      integer, intent(in) :: ibody
    end subroutine analysis_fill

    subroutine APPL_delete_itype()
    end subroutine APPL_delete_itype

    subroutine initpdfsetbynamem(iset, setname)
      integer, intent(in) :: iset
      character(len=*), intent(in) :: setname
    end subroutine initpdfsetbynamem

    subroutine numberPDFm(iset, nmem)
      integer, intent(in) :: iset
      integer, intent(out) :: nmem
    end subroutine numberPDFm

    subroutine InitPDFm(iset, imem)
      integer, intent(in) :: iset, imem
    end subroutine InitPDFm

    subroutine plot_pine_bridge(action, norm, ibody, itype, www)
      integer, intent(in) :: action, ibody, itype
      double precision, intent(in) :: norm, www
    end subroutine plot_pine_bridge
  end interface

contains

  subroutine initplot_impl()
    implicit none

    character(len=50), allocatable :: weights_info(:)
    character(len=13) :: temp
    integer :: ii, jj, kk, n, nn, nwgt

    call validate_process_dimensions()
    call validate_extra_weights()

    ! Determine the PDF members before allocating the weight descriptions.
    nwgt = 1
    if (do_rwgt_scale) then
      do kk = 1, dyn_scale(0)
        if (lscalevar(kk)) then
          nwgt = nwgt + nint(scalevarF(0))*nint(scalevarR(0))
        else
          nwgt = nwgt + 1
        end if
      end do
    end if

    if (do_rwgt_pdf) then
      do nn = 1, lhaPDFid(0)
        if (lpdfvar(nn)) then
          write (*, *) 'Including central PDF with uncertainties for ', &
               trim(LHAPDFsetname(nn))
        else
          write (*, *) 'Including central PDF for ', &
               trim(LHAPDFsetname(nn))
        end if

        ! The first PDF set was already loaded by setrun.
        if (nn > 1) then
          call initpdfsetbynamem(nn, LHAPDFsetname(nn))
          if (lpdfvar(nn)) then
            call numberPDFm(nn, nmemPDF(nn))
            if (nmemPDF(nn) == 1) then
              nmemPDF(nn) = 0
              lpdfvar(nn) = .false.
            end if
          else
            nmemPDF(nn) = 0
          end if
        end if

        if (nmemPDF(nn) + 1 > maxPDFs) then
          write (*, *) 'Too many PDFs: increase maxPDFs in ', &
               'extra_weights.f90 to ', nmemPDF(nn) + 1
          stop 1
        end if

        if (lpdfvar(nn)) then
          nwgt = nwgt + nmemPDF(nn) + 1
        else
          nwgt = nwgt + 1
        end if
      end do
      call InitPDFm(1, 0)
    end if

    allocate(weights_info(nwgt))
    weights_info = ''
    weights_info(1) = 'central value'
    nwgt = 1

    if (do_rwgt_scale) then
      do kk = 1, dyn_scale(0)
        if (lscalevar(kk)) then
          do ii = 1, nint(scalevarF(0))
            do jj = 1, nint(scalevarR(0))
              nwgt = nwgt + 1
              write (weights_info(nwgt), &
                   '(a4,i4,1x,a4,f6.3,1x,a4,f6.3)') &
                   'dyn=', dyn_scale(kk), 'muR=', scalevarR(jj), &
                   'muF=', scalevarF(ii)
            end do
          end do
        else
          nwgt = nwgt + 1
          write (weights_info(nwgt), &
               '(a4,i4,1x,a4,f6.3,1x,a4,f6.3)') &
               'dyn=', dyn_scale(kk), 'muR=', scalevarR(1), &
               'muF=', scalevarF(1)
        end if
      end do
    end if

    if (do_rwgt_pdf) then
      do nn = 1, lhaPDFid(0)
        if (lpdfvar(nn)) then
          do n = 0, nmemPDF(nn)
            nwgt = nwgt + 1
            write (temp, '(a4,i8)') 'PDF=', lhaPDFid(nn) + n
            write (weights_info(nwgt), '(a)') &
                 trim(adjustl(temp)) // ' ' // &
                 trim(adjustl(LHAPDFsetname(nn)))
          end do
        else
          nwgt = nwgt + 1
          write (temp, '(a4,i8)') 'PDF=', lhaPDFid(nn)
          write (weights_info(nwgt), '(a)') &
               trim(adjustl(temp)) // ' ' // &
               trim(adjustl(LHAPDFsetname(nn)))
        end if
      end do
    end if

    call plot_pine_bridge(pine_reset_action, 0d0, 0, 0, 0d0)
    call analysis_begin(nwgt, weights_info)

    ! Keep track of accumulated scale and PDF results.
    do kk = 1, dyn_scale(0)
      if (lscalevar(kk)) then
        do ii = 1, nint(scalevarF(0))
          do jj = 1, nint(scalevarR(0))
            xsecScale_acc(jj, ii, kk) = 0d0
          end do
        end do
      else
        xsecScale_acc(1, 1, kk) = 0d0
      end if
    end do
    do nn = 1, lhaPDFid(0)
      if (lpdfvar(nn)) then
        do n = 0, nmemPDF(nn)
          xsecPDFr_acc(n, nn) = 0d0
        end do
      else
        xsecPDFr_acc(0, nn) = 0d0
      end if
    end do
  end subroutine initplot_impl


  subroutine topout_impl()
    implicit none

    double precision :: xnorm
    integer :: ii, jj, kk, n, nn

    call validate_extra_weights()

    xnorm = 1d0/float(ncalls0)
    if (useitmax) xnorm = xnorm/float(itmax)

    ! Normalization factor for the PineAPPL grids.
    call plot_pine_bridge(pine_norm_action, &
         1d0/(dble(ncalls0)*dble(itmax)), 0, 0, 0d0)
    call analysis_end(xnorm)

    open (unit=34, file='scale_pdf_dependence.dat', status='unknown')
    if (.not. useitmax) xnorm = xnorm/float(itmax)

    if (do_rwgt_scale) then
      write (34, *) 'scale variations:'
      do kk = 1, dyn_scale(0)
        if (lscalevar(kk)) then
          write (34, *) dyn_scale(kk), nint(scalevarR(0)), &
               nint(scalevarF(0))
          write (34, *) ((xsecScale_acc(jj, ii, kk)*xnorm, &
               jj=1, nint(scalevarR(0))), &
               ii=1, nint(scalevarF(0)))
        else
          write (34, *) dyn_scale(kk), 1, 1
          write (34, *) xsecScale_acc(1, 1, kk)*xnorm
        end if
      end do
    end if

    if (do_rwgt_pdf) then
      write (34, *) 'pdf variations:'
      do nn = 1, lhaPDFid(0)
        if (lpdfvar(nn)) then
          write (34, *) trim(adjustl(LHAPDFsetname(nn))), &
               nmemPDF(nn) + 1
          write (34, *) (xsecPDFr_acc(n, nn)*xnorm, &
               n=0, nmemPDF(nn))
        else
          write (34, *) LHAPDFsetname(nn), nmemPDF(nn) + 1
          write (34, *) xsecPDFr_acc(0, nn)*xnorm
        end if
      end do
    end if
    write (34, *) ' '
    close (34)
  end subroutine topout_impl


  subroutine outfun_impl(pp, ybst_til_tolab, www, ipdg, itype, pmass)
    ! MadFKS supplies pp in the reduced parton centre-of-mass frame.  A
    ! rapidity ybst_til_tolab transforms it according to
    ! y_lab = y_cm - ybst_til_tolab.
    implicit none

    double precision, intent(in) :: pp(0:3, nexternal)
    double precision, intent(in) :: ybst_til_tolab
    double precision, intent(in) :: www(*)
    integer, intent(in) :: ipdg(nexternal)
    integer, intent(in) :: itype
    double precision, intent(in) :: pmass(nexternal)

    double precision :: chybst, chybstmo, p(0:4, nexternal)
    double precision :: pplab(0:3, nexternal), shybst
    double precision, parameter :: xd(3) = (/ 0d0, 0d0, 1d0 /)
    integer :: i, ibody, ii, i_wgt, jj, kk, n, nn
    integer :: istatus(nexternal)

    call validate_process_dimensions()
    call validate_extra_weights()

    select case (itype)
    case (11)
      ibody = 1                 ! (n+1)-body
    case (12:14)
      ibody = 2                 ! n-body
    case (20)
      ibody = 3                 ! Born
    case default
      write (*, *) 'Error in outfun: unknown itype', itype
      stop 1
    end select

    chybst = cosh(ybst_til_tolab)
    shybst = sinh(ybst_til_tolab)
    chybstmo = chybst - 1d0
    do i = 1, nexternal
      call boostwdir2(chybst, shybst, chybstmo, xd, pp(0, i), &
           pplab(0, i))
    end do

    do i = 1, nexternal
      if (i <= nincoming) then
        istatus(i) = -1
      else
        istatus(i) = 1
      end if
      p(0:3, i) = pplab(0:3, i)
      p(4, i) = pmass(i)
    end do

    call plot_pine_bridge(pine_event_action, 0d0, ibody, itype, www(1))
    call analysis_fill(p, istatus, ipdg, www, ibody)
    if (pineappl) then
      ! PineAPPL combines contributions with identical kinematics while
      ! histograms are filled contribution by contribution.
      call APPL_delete_itype()
    end if

    i_wgt = 1
    if (do_rwgt_scale) then
      do kk = 1, dyn_scale(0)
        if (lscalevar(kk)) then
          do ii = 1, nint(scalevarF(0))
            do jj = 1, nint(scalevarR(0))
              i_wgt = i_wgt + 1
              xsecScale_acc(jj, ii, kk) = &
                   xsecScale_acc(jj, ii, kk) + www(i_wgt)
            end do
          end do
        else
          i_wgt = i_wgt + 1
          xsecScale_acc(1, 1, kk) = &
               xsecScale_acc(1, 1, kk) + www(i_wgt)
        end if
      end do
    end if

    if (do_rwgt_pdf) then
      do nn = 1, lhaPDFid(0)
        if (lpdfvar(nn)) then
          do n = 0, nmemPDF(nn)
            i_wgt = i_wgt + 1
            xsecPDFr_acc(n, nn) = xsecPDFr_acc(n, nn) + www(i_wgt)
          end do
        else
          i_wgt = i_wgt + 1
          xsecPDFr_acc(0, nn) = xsecPDFr_acc(0, nn) + www(i_wgt)
        end if
      end do
    end if
  end subroutine outfun_impl

end module madfks_plot_module

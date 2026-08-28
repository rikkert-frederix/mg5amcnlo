! Wrapper routines for the fixed-order analyses.
module madfks_plot_module
  use boostwdir2_module, only: boostwdir2
  use extra_weights
  use mint_module, only: itmax, ncalls0
  use process_dimensions, only: nexternal, event_capacity, nincoming, &
       validate_process_dimensions
  use run_state, only: do_rwgt_pdf, do_rwgt_scale, do_rwgt_decay_scale
  use fnlo_process_common, only: xsecPDFr_acc
  use fnlo_scale_variations, only: fnlo_scale_point_count, &
       fnlo_scale_max_point_count, fnlo_scale_point_label, &
       fnlo_scale_mode_name
  implicit none
  private

  double precision, allocatable, save :: xsec_scale_points(:, :)

  public :: initplot_impl, topout_impl, outfun_impl
  public :: outfun_multiplicative_impl

  interface
    subroutine analysis_begin(nwgt, weights_info)
      integer, intent(in) :: nwgt
      character(len=*), intent(in) :: weights_info(*)
    end subroutine analysis_begin

    subroutine analysis_end()
    end subroutine analysis_end

    subroutine analysis_fill(p, istatus, ipdg, wgts, ibody)
      use process_dimensions, only: event_capacity
      double precision, intent(in) :: p(0:4, event_capacity)
      integer, intent(in) :: istatus(event_capacity)
      integer, intent(in) :: ipdg(event_capacity)
      double precision, intent(in) :: wgts(*)
      integer, intent(in) :: ibody
    end subroutine analysis_fill

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

  end interface

contains

  subroutine initplot_impl()
    implicit none

    character(len=200), allocatable :: weights_info(:)
    character(len=13) :: temp
    integer :: kk, n, nn, nwgt, point

    call validate_process_dimensions()
    ! Determine the PDF members before allocating the weight descriptions.
    nwgt = 1
    if (do_rwgt_scale .or. do_rwgt_decay_scale) then
      do kk = 1, dyn_scale(0)
        nwgt = nwgt + fnlo_scale_point_count(kk)
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

    ! MINT can restart the integration and call initplot again after resetting
    ! the grids.  Some supported builds retain local variables between calls,
    ! so release the previous descriptions before rebuilding them.
    if (allocated(weights_info)) deallocate(weights_info)
    allocate(weights_info(nwgt))
    weights_info = ''
    weights_info(1) = 'central value'
    nwgt = 1

    if (do_rwgt_scale .or. do_rwgt_decay_scale) then
      do kk = 1, dyn_scale(0)
        do point = 1, fnlo_scale_point_count(kk)
          nwgt = nwgt + 1
          call fnlo_scale_point_label(&
               kk, point, weights_info(nwgt))
        end do
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

    call analysis_begin(nwgt, weights_info)

    ! Keep track of accumulated scale and PDF results.
    if (allocated(xsec_scale_points)) deallocate(xsec_scale_points)
    allocate(xsec_scale_points(&
         fnlo_scale_max_point_count(), dyn_scale(0)))
    xsec_scale_points = 0d0
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
    integer :: kk, n, nn, point
    character(len=16) :: scale_mode

    xnorm = 1d0/float(ncalls0)
    call analysis_end()

    open (unit=34, file='scale_pdf_dependence.dat', status='unknown')
    xnorm = xnorm/float(itmax)

    if (do_rwgt_scale .or. do_rwgt_decay_scale) then
      call fnlo_scale_mode_name(scale_mode)
      write (34, *) 'scale variations:'
      do kk = 1, dyn_scale(0)
        write (34, *) dyn_scale(kk), fnlo_scale_point_count(kk), 1, &
             trim(scale_mode)
        write (34, *) (xsec_scale_points(point, kk)*xnorm, &
             point=1, fnlo_scale_point_count(kk))
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

    double precision :: chybst, chybstmo, p(0:4, event_capacity)
    double precision :: pplab(0:3, event_capacity), shybst
    double precision, parameter :: xd(3) = (/ 0d0, 0d0, 1d0 /)
    integer :: i, ibody, i_wgt, kk, n, nn, point
    integer :: istatus(event_capacity), analysis_pdg(event_capacity)

    call validate_process_dimensions()
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
    p = 0d0
    pplab = 0d0
    istatus = 2
    analysis_pdg = 0
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
      analysis_pdg(i) = ipdg(i)
    end do

    call analysis_fill(p, istatus, analysis_pdg, www, ibody)

    call accumulate_analysis_weights(www)
  end subroutine outfun_impl


  subroutine outfun_multiplicative_impl( &
       pp, ybst_til_tolab, www, ipdg, origin_blocks)
    double precision, intent(in) :: pp(0:,:)
    double precision, intent(in) :: ybst_til_tolab, www(*)
    integer, intent(in) :: ipdg(:), origin_blocks(:)
    double precision :: chybst, chybstmo, shybst, mass_squared
    double precision :: p(0:4,event_capacity)
    double precision :: pplab(0:3,event_capacity)
    double precision, parameter :: direction(3) = (/0d0, 0d0, 1d0/)
    integer :: particle, istatus(event_capacity)
    integer :: analysis_pdg(event_capacity)

    if (size(pp,1) /= 4 .or. size(pp,2) /= event_capacity .or. &
        size(ipdg) /= event_capacity .or. &
        size(origin_blocks) /= event_capacity) then
      write (*, *) 'Invalid multiplicative analysis-event shape'
      stop 1
    end if
    p = 0d0
    pplab = 0d0
    istatus = 2
    analysis_pdg = ipdg
    chybst = cosh(ybst_til_tolab)
    shybst = sinh(ybst_til_tolab)
    chybstmo = chybst - 1d0
    do particle = 1, event_capacity
      call boostwdir2(chybst, shybst, chybstmo, direction, &
           pp(:,particle), pplab(:,particle))
      p(0:3,particle) = pplab(:,particle)
      mass_squared = (pp(0,particle) + pp(3,particle))* &
           (pp(0,particle) - pp(3,particle)) - &
           pp(1,particle)**2 - pp(2,particle)**2
      p(4,particle) = sqrt(max(0d0, mass_squared))
      if (particle <= nincoming) then
        istatus(particle) = -1
      else if (ipdg(particle) /= 0) then
        istatus(particle) = 1
      end if
      if (analysis_pdg(particle) == -21) analysis_pdg(particle) = 21
    end do
    ! Analyses always receive the complete event, including forced-decay
    ! products.  ORIGIN_BLOCKS is shape-checked here and intentionally not
    ! used as a cut mask.
    if (any(origin_blocks < -1)) then
      write (*, *) 'Invalid multiplicative analysis-event origin'
      stop 1
    end if
    call analysis_fill(p, istatus, analysis_pdg, www, 1)
    call accumulate_analysis_weights(www)
  end subroutine outfun_multiplicative_impl


  subroutine accumulate_analysis_weights(www)
    double precision, intent(in) :: www(*)
    integer :: i_wgt, kk, n, nn, point

    i_wgt = 1
    if (do_rwgt_scale .or. do_rwgt_decay_scale) then
      do kk = 1, dyn_scale(0)
        do point = 1, fnlo_scale_point_count(kk)
          i_wgt = i_wgt + 1
          xsec_scale_points(point, kk) = &
               xsec_scale_points(point, kk) + www(i_wgt)
        end do
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
  end subroutine accumulate_analysis_weights

end module madfks_plot_module

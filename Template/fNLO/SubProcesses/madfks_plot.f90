! Wrapper routines for the fixed-order analyses.
module madfks_plot_module
  use boostwdir2_module, only: boostwdir2
  use extra_weights
  use mint_module, only: itmax, ncalls0
  use process_dimensions, only: nexternal, nincoming, &
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
      use process_dimensions, only: nexternal
      double precision, intent(in) :: p(0:4, nexternal)
      integer, intent(in) :: istatus(nexternal), ipdg(nexternal)
      double precision, intent(in) :: wgts(*)
      integer, intent(in) :: ibody
    end subroutine analysis_fill

    subroutine analysis_fill_multiplicative( &
         p, particle_count, istatus, ipdg, wgts, ibody)
      integer, intent(in) :: particle_count
      double precision, intent(in) :: p(0:4, particle_count)
      integer, intent(in) :: istatus(particle_count)
      integer, intent(in) :: ipdg(particle_count)
      double precision, intent(in) :: wgts(*)
      integer, intent(in) :: ibody
    end subroutine analysis_fill_multiplicative

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

    double precision :: chybst, chybstmo, p(0:4, nexternal)
    double precision :: pplab(0:3, nexternal), shybst
    double precision, parameter :: xd(3) = (/ 0d0, 0d0, 1d0 /)
    integer :: i, ibody
    integer :: istatus(nexternal)

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
    do i = 1, nexternal
      call boostwdir2(chybst, shybst, chybstmo, xd, pp(0:3, i), &
           pplab(0:3, i))
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

    call analysis_fill(p, istatus, ipdg, www, ibody)

    call accumulate_analysis_weights(www)
  end subroutine outfun_impl


  subroutine outfun_multiplicative_impl( &
       pp, particle_count, ybst_til_tolab, www, ipdg, istatus, ibody)
    ! Multiplicative leaves can contain one resolved emission from every
    ! corrected block, so their particle count is not bounded by the legacy
    ! one-real ``nexternal`` layout.
    implicit none

    integer, intent(in) :: particle_count, ibody
    double precision, intent(in) :: pp(0:, :), ybst_til_tolab
    double precision, intent(in) :: www(*)
    integer, intent(in) :: ipdg(:), istatus(:)

    double precision :: chybst, chybstmo
    double precision :: p(0:4, particle_count)
    double precision :: pplab(0:3, particle_count), mass_squared, shybst
    double precision, parameter :: xd(3) = (/ 0d0, 0d0, 1d0 /)
    integer :: i

    call validate_process_dimensions()
    if (particle_count < nincoming .or. particle_count > size(pp, 2) .or. &
        particle_count > size(ipdg) .or. &
        particle_count > size(istatus)) then
      write (*, *) 'Error in multiplicative outfun: invalid particle count', &
           particle_count
      stop 1
    end if
    if (ibody /= 1 .and. ibody /= 2) then
      write (*, *) 'Error in multiplicative outfun: invalid body type', ibody
      stop 1
    end if

    chybst = cosh(ybst_til_tolab)
    shybst = sinh(ybst_til_tolab)
    chybstmo = chybst - 1d0
    do i = 1, particle_count
      call boostwdir2(chybst, shybst, chybstmo, xd, pp(0:3, i), &
           pplab(0:3, i))
      p(0:3, i) = pplab(0:3, i)
      mass_squared = pplab(0, i)**2 - sum(pplab(1:3, i)**2)
      p(4, i) = sqrt(max(0d0, mass_squared))
    end do

    call analysis_fill_multiplicative( &
         p, particle_count, istatus(1:particle_count), &
         ipdg(1:particle_count), www, ibody)

    call accumulate_analysis_weights(www)
  end subroutine outfun_multiplicative_impl


  subroutine accumulate_analysis_weights(www)
    implicit none
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

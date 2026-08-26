module binoth_lha_olp_backend
  use FKSParams, only: IRPoleCheckThreshold
  use process_dimensions, only: nexternal
  use fks_singular_module, only: getpoles
  use fnlo_process_common, only: isum_hel
  implicit none
  private

  logical, save :: firsttime_pole = .true.
  logical, save :: firsttime_init = .true.
  integer, save :: nbad = 0

  interface
    subroutine binoth_lha_update_couplings(mu_r_value, alpha_s)
      implicit none
      double precision, intent(out) :: mu_r_value, alpha_s
    end subroutine binoth_lha_update_couplings
  end interface

  public :: binoth_lha_eval

contains

  subroutine binoth_lha_eval(pin, born_wgt, virt_wgt, proc_label, pmass)
    implicit none
    double precision, intent(in) :: pin(0:3, nexternal-1)
    double precision, intent(in) :: born_wgt
    double precision, intent(out) :: virt_wgt
    integer, intent(in) :: proc_label
    double precision, intent(in) :: pmass(nexternal)

    integer, parameter :: nbadmax = 5
    double precision :: p(0:4, nexternal-1)
    double precision :: virt_wgts(4)
    double precision :: double_pole, single_pole, born
    double precision :: mu_r_value, alpha_s
    double precision :: tolerance, madfks_single, madfks_double
    integer :: i, j
    if (isum_hel /= 0) then
      write (*,*) 'Can only do explicit helicity sum' // &
           ' for Virtual corrections', isum_hel
    end if
    virt_wgt = 0d0

    call binoth_lha_update_couplings(mu_r_value, alpha_s)
    do i = 1, nexternal-1
      do j = 0, 3
        p(j,i) = pin(j,i)
      end do
      p(4,i) = pmass(i)
    end do

    if (firsttime_init) then
      call binoth_lha_init_impl()
      firsttime_init = .false.
    end if
    call OLP_EvalSubProcess(proc_label, p, mu_r_value, alpha_s, &
         virt_wgts)
    double_pole = virt_wgts(1)
    single_pole = virt_wgts(2)
    virt_wgt = virt_wgts(3)
    born = virt_wgts(4)

    if (firsttime_pole) then
      tolerance = IRPoleCheckThreshold
      call getpoles(pin, madfks_double, madfks_single)
      if (dabs(single_pole-madfks_single) < tolerance .and. &
          dabs(double_pole-madfks_double) < tolerance) then
        write (*,*) '---- POLES CANCELLED ----'
        firsttime_pole = .false.
      else
        write (*,*) 'POLES MISCANCELLATION, DIFFERENCE > ', tolerance
        write (*,*) ' BORN:'
        write (*,*) '       MadFKS: ', born_wgt, '          OLP: ', born
        write (*,*) ' COEFFICIENT DOUBLE POLE:'
        write (*,*) '       MadFKS: ', madfks_double, &
             '          OLP: ', double_pole
        write (*,*) ' COEFFICIENT SINGLE POLE:'
        write (*,*) '       MadFKS: ', madfks_single, &
             '          OLP: ', single_pole
        write (*,*) ' FINITE:'
        write (*,*) '          OLP: ', virt_wgt
        if (nbad < nbadmax) then
          nbad = nbad + 1
          write (*,*) ' Trying another PS point'
        else
          write (*,*) 'ERROR: TOO MANY FAILURES, QUITTING'
          stop
        end if
      end if
    end if
  end subroutine binoth_lha_eval


  subroutine binoth_lha_init_impl()
    implicit none
    character(len=13) :: filename
    integer :: ierr

    filename = 'OLE_order.olc'
    ierr = 0
    call OLP_Start(filename // char(0), ierr)
    if (ierr == 0) then
      write (*,*) 'ERROR in the BinothLHAInit process initialization'
      stop
    end if
  end subroutine binoth_lha_init_impl



end module binoth_lha_olp_backend

module setrun_module
  use extra_weights, only: dyn_scale, lhaPDFid, lpdfvar, lscalevar, &
       nmemPDF
  use run_state
  implicit none
  private

  public :: complete_setrun

contains

  subroutine complete_setrun()
    implicit none

    integer :: i
    double precision :: strong_coupling

    interface
      subroutine setrun_model_strong_coupling(value)
        implicit none
        double precision, intent(out) :: value
      end subroutine setrun_model_strong_coupling
    end interface

    ! Determine whether scale and/or PDF reweighting is needed.
    do_rwgt_scale = .false.
    do i = 1, dyn_scale(0)
      if (lscalevar(i) .or. dyn_scale(0) > 1) then
        do_rwgt_scale = .true.
        exit
      end if
    end do

    do_rwgt_pdf = .false.
    do i = 1, lhaPDFid(0)
      if (lpdfvar(i) .or. lhaPDFid(0) > 1) then
        do_rwgt_pdf = .true.
        exit
      end if
    end do

    ! Set the default scale and PDF choices used for the actual run.
    dynamical_scale_choice = dyn_scale(1)
    lhaid = lhaPDFid(1)

    scale = mur_ref_fixed
    q2fact(1) = muf1_ref_fixed**2
    q2fact(2) = muf2_ref_fixed**2

    ! Set alpha_s(mZ).  The model coupling is obtained through the fixed-form
    ! bridge, which is the only code in this unit that includes coupl.inc.
    if (lpp(1) /= 0 .or. lpp(2) /= 0) then
      write(*, *) 'A PDF is used, so alpha_s(MZ) is going to be modified'
      call setpara('param_card.dat')
      call setrun_model_strong_coupling(strong_coupling)
      asmz = strong_coupling**2 / (16d0 * atan(1d0))
      write(*, *) 'Old value of alpha_s from param_card: ', asmz
      call pdfwrap
      write(*, *) 'New value of alpha_s from PDF ', pdlabel, ':', asmz
    else
      call setpara('param_card.dat')
      call setrun_model_strong_coupling(strong_coupling)
      asmz = strong_coupling**2 / (16d0 * atan(1d0))
      nloop = 2
      pdlabel = 'none'
      write(*, *)
      write(*, *) 'No PDF is used, alpha_s(MZ) from param_card is used'
      write(*, *) 'Value of alpha_s from param_card: ', asmz
      write(*, *) 'The default order of alpha_s running is fixed to ', &
           nloop
    end if

    ! Fill the number of PDF error members using LHAPDF.
    if (lpdfvar(1) .and. (lpp(1) /= 0 .or. lpp(2) /= 0)) then
      call numberPDFm(1, nmemPDF(1))
      if (nmemPDF(1) == 1) then
        nmemPDF(1) = 0
        lpdfvar(1) = .false.
      end if
    else
      nmemPDF(1) = 0
    end if
  end subroutine complete_setrun

end module setrun_module

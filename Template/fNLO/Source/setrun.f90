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
    integer :: idum
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

    call fill_needed_splittings()

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

    if (pineappl) then
      call fill_lhe_initial_state(idum)
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


  subroutine fill_lhe_initial_state(idum)
    implicit none

    integer, intent(out) :: idum
    integer :: i
    integer :: read_status

    idum = 0
    do i = 1, 2
      if (lpp(i) == 1) then
        idbmup(i) = 2212
      else if (lpp(i) == -1) then
        idbmup(i) = -2212
      else if (lpp(i) == 0) then
        open(unit=71, status='old', file='initial_states_map.dat')
        read(71, *, iostat=read_status) idum, idum, idbmup(1), idbmup(2)
        if (read_status /= 0) then
          write(*, *) '"initial_states_map.dat" not found (or incorrect' // &
               ' format) by "Source/setrun"'
          stop 1
        end if
        close(71)
      else
        write(*, *) 'Unsupported fNLO beam type:', lpp(i)
        stop 1
      end if
      ebmup(i) = ebeam(i)
    end do

    if (abs(lpp(1)) == 1 .or. abs(lpp(2)) == 1) then
      call get_pdfup(pdlabel, pdfgup, pdfsup, lhaid)
    end if
  end subroutine fill_lhe_initial_state


  subroutine get_pdfup(pdfin, output_pdf_groups, output_pdf_sets, &
                       input_lhaid)
    implicit none

    character(len=*), intent(in) :: pdfin
    integer, intent(out) :: output_pdf_groups(2)
    integer, intent(out) :: output_pdf_sets(2)
    integer, intent(in) :: input_lhaid

    integer, parameter :: number_of_pdfs = 4
    character(len=7), parameter :: pdf_labels(number_of_pdfs) = (/ &
         'none   ', 'nn23lo ', 'nn23lo1', 'nn23nlo' /)
    integer, parameter :: pdf_numbers(number_of_pdfs) = (/ &
         0, 246800, 247000, 244800 /)

    integer :: i
    integer :: matched_pdf

    if (pdfin == 'lhapdf') then
      write(*, *) 'using LHAPDF'
      output_pdf_groups = -1
      output_pdf_sets = input_lhaid
      return
    end if

    matched_pdf = -1
    do i = 1, number_of_pdfs
      if (trim(pdfin) == trim(pdf_labels(i))) then
        matched_pdf = pdf_numbers(i)
      end if
    end do

    if (matched_pdf == -1) then
      write(*, *) 'ERROR: pdf ', pdfin, ' not implemented in get_pdfup.'
      write(*, *) 'known pdfs are'
      write(*, *) pdf_labels
      stop 1
    end if

    output_pdf_groups = -1
    output_pdf_sets = matched_pdf
  end subroutine get_pdfup

end module setrun_module

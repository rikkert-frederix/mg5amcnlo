module run_printout_module
  implicit none
  private

  public :: write_run_summary

contains

  subroutine write_run_summary(lpp, ebeam, pdlabel, asmz, nloop, &
                               fixed_ren_scale, scale, fixed_fac_scale, &
                               q2fact)
    implicit none
    integer, intent(in) :: lpp(2)
    double precision, intent(in) :: ebeam(2)
    character(len=*), intent(in) :: pdlabel
    double precision, intent(in) :: asmz
    integer, intent(in) :: nloop
    logical, intent(in) :: fixed_ren_scale
    double precision, intent(in) :: scale
    logical, intent(in) :: fixed_fac_scale
    double precision, intent(in) :: q2fact(2)

    character(len=2) :: beam_name(2)
    double precision :: energy
    integer :: i

    write(6, *)
    write(6, *) 'Collider parameters:'
    write(6, *) '--------------------'

    beam_name = '?'
    do i = 1, 2
      if (lpp(i) == 0) beam_name(i) = 'e'
      if (lpp(i) == 1) beam_name(i) = 'P'
      if (lpp(i) == -1) beam_name(i) = 'Pb'
      if (lpp(i) == 2) beam_name(i) = 'a2'
      if (lpp(i) == 3) beam_name(i) = 'e-'
      if (lpp(i) == -3) beam_name(i) = 'e+'
      if (lpp(i) == 4) beam_name(i) = 'm-'
      if (lpp(i) == -4) beam_name(i) = 'm+'
    end do

    energy = 2d0 * sqrt(ebeam(1) * ebeam(2))

    write(6, *)
    write(6, *) 'Running at ', beam_name(1), beam_name(2), &
                '  machine @ ', energy, ' GeV'
    write(6, *) 'PDF set = ', pdlabel
    write(6, '(1x,a12,1x,f6.4,a12,i1,a7)') &
         'alpha_s(Mz)=', asmz, ' running at ', nloop, ' loops.'
    if (lpp(1) /= 0 .or. lpp(2) /= 0) then
      write(6, '(1x,a12,1x,f6.4,a12,i1,a7)') &
           'alpha_s(Mz)=', asmz, ' running at ', nloop, &
           ' loops. Value tuned to the PDF set.'
    else
      write(6, '(1x,a12,1x,f6.4,a12,i1,a7)') &
           'alpha_s(Mz)=', asmz, ' running at ', nloop, &
           ' loops. Value set in param_card.dat'
    end if

    if (fixed_ren_scale) then
      write(6, *) 'Renormalization scale fixed @ ', scale
    else
      write(6, *) 'Renormalization scale set on event-by-event basis'
    end if
    if (fixed_fac_scale) then
      write(6, *) 'Factorization scales  fixed @ ', &
                  sqrt(q2fact(1)), sqrt(q2fact(2))
    else
      write(6, *) 'Factorization   scale set on event-by-event basis'
    end if

    write(6, *)
    write(6, *)
  end subroutine write_run_summary

end module run_printout_module

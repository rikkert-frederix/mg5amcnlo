program check_runtime_coupling
  use run_state, only: pdlabel, lhaid
  use alfas_functions_module, only: alphas
  implicit none
  double precision, parameter :: scales(5) = [80.385d0, 86.25d0, 91.1876d0, 172.5d0, 345d0]
  integer :: i

  ! Link the SAME compiled PDF setup and coupling objects as madevent_mintFO.
  ! The test changes no cards or source in the generated process.
  pdlabel = 'lhapdf'
  lhaid = 331700
  call pdfwrap()
  do i = 1, size(scales)
    write (*, '(a,2(1x,es24.16))') 'COUPLING', scales(i), alphas(scales(i))
  end do
end program check_runtime_coupling

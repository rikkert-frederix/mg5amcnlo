program ranmar_stage_driver
  use fnlo_runtime_common, only: iseed, random_offset_split
  use ranmar_module, only: ntuple
  implicit none
  integer :: configuration, i
  double precision :: value

  read(*, *) configuration, random_offset_split
  do i = 1, 8
    call ntuple(value, 0d0, 1d0, configuration)
    write(*, '(a,es24.16)') 'DRAW ', value
  end do
  write(*, '(a,i0)') 'BASE_SEED ', iseed
  open(unit=21, file='res.dat', status='replace')
  write(21, *) 'completed RNG fixture'
  close(21)
end program ranmar_stage_driver

program hwu_iteration_checks
  use HwU_module
  implicit none
  integer :: repeat, point
  character(len=80) :: labels(2)

  labels = [character(len=80) :: 'central value', 'variation']
  do repeat = 1, 2
    call HwU_inithist(2, labels)
    call HwU_book(1, 'iteration averaging', 4, 0d0, 4d0)
    ! Exercise reallocation of both the label and bin dimensions.
    call HwU_book(9, 'empty control', 6, 0d0, 6d0)

    ! A discarded iteration must not initialize the accumulated estimator.
    call HwU_fill(1, 0.5d0, [100d0, 200d0])
    call HwU_add_points()
    call HwU_accum_iter(.false., 4, [0d0, 0d0])

    ! First bin only in iteration one; third bin has exactly zero variance.
    do point = 1, 4
      if (point == 1) call HwU_fill(1, 0.5d0, [4d0, 8d0])
      call HwU_fill(1, 2.5d0, [1d0, 3d0])
      call HwU_fill(1, 3.5d0, [2d0, 4d0])
      call HwU_fill(1, 3.5d0, [-2d0, -4d0])
      call HwU_add_points()
    end do
    call HwU_accum_iter(.true., 4, [1d0, 0d0])

    ! Second bin only in iteration two. Global weights are 2/3 and 1/3.
    do point = 1, 4
      if (point == 1) call HwU_fill(1, 1.5d0, [4d0, 8d0])
      call HwU_fill(1, 2.5d0, [2d0, 6d0])
      call HwU_add_points()
    end do
    call HwU_accum_iter(.true., 4, [2d0, 1d0])
    call HwU_output(6, 1d0)

    ! An entirely empty included iteration also carries the global weight.
    call HwU_accum_iter(.true., 4, [3d0, 2d0/3d0])
    call HwU_output(6, 1d0)
  end do
end program hwu_iteration_checks

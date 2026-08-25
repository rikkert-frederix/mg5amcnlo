! Momentum reshuffling routines based on those used in MadSTR
! (arXiv:1907.04898).
module momentum_reshuffling
  use process_dimensions, only: nexternal, nincoming
  use kin_functions_module, only: dot => dot_impl, &
       sumdot => sumdot_impl, threedot => threedot_impl
  implicit none
  private

  public :: reshuffle_momenta

contains

  double precision function lambda_tr(x, y, z)
    double precision, intent(in) :: x, y, z

    lambda_tr = x**2 + y**2 + z**2 - 2d0*x*y - 2d0*x*z - 2d0*y*z
  end function lambda_tr


  double precision function lambda2(a, b, c)
    double precision, intent(in) :: a, b, c

    if (a <= 0d0 .or. abs(b + c) > abs(a) .or. &
        abs(b - c) > abs(a)) then
      write (*, *) 'Error #1 in lambda2: inputs not consistent', a, b, c
      stop 1
    end if
    lambda2 = sqrt(1d0 - (b + c)**2/a**2) * &
         sqrt(1d0 - (b - c)**2/a**2)
  end function lambda2


  subroutine reshuffle_momenta_old(p, q, iresh, pdg_old, pdg_new, pass)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(out) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh, pdg_old, pdg_new
    logical, intent(out) :: pass
    integer, parameter :: ihowresh = 1

    if (iresh > nincoming) then
      if (ihowresh == 1) then
        call reshuffle_initial(p, q, iresh, pdg_old, pdg_new, pass)
      else if (ihowresh == 2) then
        call reshuffle_final(p, q, iresh, pdg_old, pdg_new, pass)
      else
        write (*, *) 'ERROR: reshuffle momenta, wrong option', ihowresh
        stop 1
      end if
    else
      call reshuffle_initial_state(p, q, iresh, pdg_old, pdg_new, pass)
    end if
  end subroutine reshuffle_momenta_old


  ! Select the reshuffling strategy for as many as two particles.
  subroutine reshuffle_momenta(p, q, iresh, pdg_old, pdg_new, pass)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(out) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh(2), pdg_old(2), pdg_new(2)
    logical, intent(out) :: pass
    integer, parameter :: ihowresh = 1
    logical, parameter :: reshuffle_two = .true.
    integer :: i

    do i = 1, 2
      if (iresh(i) == 0 .or. pdg_old(i) == 0 .or. pdg_new(i) == 0) then
        if (iresh(i) /= 0 .or. pdg_old(i) /= 0 .or. pdg_new(i) /= 0) then
          write (*, *) 'ERROR in reshuffling, inconsistent values', &
               i, iresh, pdg_old, pdg_new
          stop 1
        end if
      end if
    end do

    if (reshuffle_two .and. iresh(1) > nincoming .and. &
        iresh(2) > nincoming) then
      call reshuffle_final_two(p, q, iresh, pdg_old, pdg_new, pass)
    else
      do i = 2, 1, -1
        ! Retain the original no-op CONTINUE for a zero padded entry.
        if (iresh(i) == 0) continue

        if (iresh(i) > nincoming) then
          if (ihowresh == 1) then
            call reshuffle_initial(p, q, iresh(i), pdg_old(i), &
                 pdg_new(i), pass)
          else if (ihowresh == 2) then
            call reshuffle_final(p, q, iresh(i), pdg_old(i), &
                 pdg_new(i), pass)
          else
            write (*, *) 'ERROR: reshuffle momenta, wrong option', ihowresh
            stop 1
          end if
        else
          call reshuffle_initial_state(p, q, iresh(i), pdg_old(i), &
               pdg_new(i), pass)
        end if
        if (.not. pass) return
      end do
    end if
  end subroutine reshuffle_momenta


  ! Change two final-state masses while conserving their combined momentum
  ! and invariant mass.
  subroutine reshuffle_final_two(p, q, iresh, pdg_old, pdg_new, pass)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(out) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh(2), pdg_old(2), pdg_new(2)
    logical, intent(out) :: pass
    integer :: i, j
    double precision :: invm2, ptot(0:3)
    double precision :: pcom(0:3, 2), qcom(0:3, 2)
    double precision :: pspac2, qspac2
    double precision :: mass_old(2), mass_new(2)
    double precision, external :: get_mass_from_id

    pass = .true.
    do i = 1, 2
      mass_old(i) = get_mass_from_id(pdg_old(i))
      mass_new(i) = get_mass_from_id(pdg_new(i))
    end do

    if (mass_old(1) == mass_new(1) .and. &
        mass_old(2) == mass_new(2)) then
      q = p
      return
    end if

    do j = 1, nexternal - 1
      if (j == iresh(1) .or. j == iresh(2)) cycle
      q(:, j) = p(:, j)
    end do

    ptot = p(:, iresh(1)) + p(:, iresh(2))
    invm2 = dot(ptot, ptot)
    if (sqrt(invm2) < mass_new(1) + mass_new(2)) then
      pass = .false.
      return
    end if

    do i = 1, 2
      call invboostx(p(0, iresh(i)), ptot, pcom(0, i))
    end do
    pspac2 = threedot(pcom(0, 1), pcom(0, 1))
    qspac2 = lambda_tr(invm2, mass_new(1)**2, mass_new(2)**2) / &
         4d0 / invm2
    do i = 1, 2
      do j = 1, 3
        qcom(j, i) = pcom(j, i) * dsqrt(qspac2 / pspac2)
      end do
      qcom(0, i) = dsqrt(mass_new(i)**2 + qspac2)
      call boostx(qcom(0, i), ptot, q(0, iresh(i)))
    end do

    call check_reshuffled_momenta_two(p, q, iresh, mass_new)
  end subroutine reshuffle_final_two


  ! Reshuffle an initial-state particle while keeping final-state momenta
  ! fixed.
  subroutine reshuffle_initial_state(p, q, iresh, pdg_old, pdg_new, pass)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(out) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh, pdg_old, pdg_new
    logical, intent(out) :: pass
    integer :: i
    double precision :: pcom(0:3)
    double precision :: qinboost(0:3, nincoming)
    double precision :: shat, m2_other
    double precision :: mass_old, mass_new
    double precision, external :: get_mass_from_id

    pass = .true.
    mass_old = get_mass_from_id(pdg_old)
    if (mass_old /= 0d0) then
      write (*, *) 'ERROR reshuffle initial state with mass_old!=0', &
           mass_old
    end if
    if (abs(p(3, 1) + p(3, 2)) / &
        (abs(p(3, 1)) + abs(p(3, 2))) > 1d-6) then
      write (*, *) 'ERROR reshuffle initial state, no com', &
           p(:, 1), p(:, 2)
    end if
    mass_new = get_mass_from_id(pdg_new)

    if (mass_old == mass_new) then
      q = p
      return
    end if

    shat = sumdot(p(0, 1), p(0, 2), 1d0)
    do i = 1, nincoming
      if (i == iresh) then
        qinboost(0, i) = mass_new
        qinboost(1:3, i) = 0d0
      else
        m2_other = dot(p(0, i), p(0, i))
        qinboost(0, i) = (shat - m2_other - mass_new**2) / &
             2d0 / mass_new
        qinboost(1:2, i) = 0d0
        if (qinboost(0, i)**2 - m2_other < 0d0 .or. &
            qinboost(0, i) < 0d0) then
          pass = .false.
          return
        end if
        qinboost(3, i) = sign(dsqrt(qinboost(0, i)**2 - m2_other), &
             p(3, i))
      end if
    end do

    pcom = qinboost(:, 1) + qinboost(:, 2)
    do i = 1, nincoming
      call invboostx(qinboost(0, i), pcom, q(0, i))
    end do
    do i = nincoming + 1, nexternal - 1
      q(:, i) = p(:, i)
    end do

    call check_reshuffled_momenta(p, q, iresh, mass_new)
  end subroutine reshuffle_initial_state


  ! Recoil a changed final-state mass against the other final-state
  ! particles.
  subroutine reshuffle_final(p, q, iresh, pdg_old, pdg_new, pass)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(out) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh, pdg_old, pdg_new
    logical, intent(out) :: pass
    integer :: i, nu
    double precision :: preco(0:3), qreco(0:3), ptmp(0:3)
    double precision :: shat, msq_reco, totmass
    double precision :: mass_old, mass_new
    double precision, external :: get_mass_from_id

    pass = .true.
    mass_old = get_mass_from_id(pdg_old)
    mass_new = get_mass_from_id(pdg_new)
    if (mass_old == mass_new) then
      q = p
      return
    end if

    totmass = 0d0
    do i = nincoming + 1, nexternal - 1
      totmass = totmass + dsqrt(max(0d0, dot(p(0, i), p(0, i))))
    end do
    shat = 2d0 * dot(p(0, 1), p(0, 2))
    if (sqrt(shat) < totmass + mass_new - mass_old) then
      pass = .false.
      return
    end if

    preco = 0d0
    do i = nincoming + 1, nexternal - 1
      if (i == iresh) cycle
      do nu = 0, 3
        preco(nu) = preco(nu) + p(nu, i)
      end do
    end do
    msq_reco = dot(preco, preco)

    q(0, iresh) = dsqrt(shat) / 2d0 * &
         (1d0 + (mass_new**2 - msq_reco) / shat)
    qreco(0) = dsqrt(shat) / 2d0 * &
         (1d0 - (mass_new**2 - msq_reco) / shat)
    do nu = 1, 3
      q(nu, iresh) = p(nu, iresh) / &
           dsqrt(threedot(p(0, iresh), p(0, iresh))) * &
           dsqrt(q(0, iresh)**2 - mass_new**2)
      qreco(nu) = preco(nu) / dsqrt(threedot(preco, preco)) * &
           dsqrt(qreco(0)**2 - msq_reco)
    end do

    do i = 1, nexternal - 1
      if (i == iresh) cycle
      call invboostx(p(0, i), preco, ptmp)
      call boostx(ptmp, qreco, q(0, i))
    end do
    do i = 1, nincoming
      q(:, i) = p(:, i)
    end do

    call check_reshuffled_momenta(p, q, iresh, mass_new)
  end subroutine reshuffle_final


  ! Recoil a changed final-state mass against the incoming momenta.
  subroutine reshuffle_initial(p, q, iresh, pdg_old, pdg_new, pass)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(out) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh, pdg_old, pdg_new
    logical, intent(out) :: pass
    integer :: j
    double precision :: etot, ztot
    double precision :: mass_old, mass_new
    double precision, external :: get_mass_from_id

    pass = .true.
    mass_old = get_mass_from_id(pdg_old)
    mass_new = get_mass_from_id(pdg_new)
    if (mass_old == mass_new) then
      q = p
      return
    end if

    etot = 0d0
    ztot = 0d0
    q(1:3, iresh) = p(1:3, iresh)
    q(0, iresh) = dsqrt(mass_new**2 + &
         threedot(q(0, iresh), q(0, iresh)))

    do j = nincoming + 1, nexternal - 1
      if (j /= iresh) q(:, j) = p(:, j)
      etot = etot + q(0, j)
      ztot = ztot + q(3, j)
    end do

    q(0, 1) = (etot + ztot) / 2d0
    q(1:2, 1) = 0d0
    q(3, 1) = dsign(q(0, 1), p(3, 1))
    q(0, 2) = (etot - ztot) / 2d0
    q(1:2, 2) = 0d0
    q(3, 2) = dsign(q(0, 2), p(3, 2))

    call check_reshuffled_momenta(p, q, iresh, mass_new)
  end subroutine reshuffle_initial


  subroutine check_reshuffled_momenta(p, q, iresh, mass_new)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(in) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh
    double precision, intent(in) :: mass_new
    double precision :: a, b
    integer :: i, j

    if (nincoming /= 2) then
      write (*, *) 'ERROR IN OS_CHECK_MOMENTA: nincoming != 2 not ' // &
           'implemented', nincoming
      stop
    end if

    do i = 1, nexternal - 1
      if (i /= iresh) then
        if (dabs(dot(q(0, i), q(0, i)) - dot(p(0, i), p(0, i))) > &
            1d-3 * max(dot(p(0, i), p(0, i)), p(0, i)**2)) then
          write (*, *) 'ERROR IN CHECK_RESHUFFLED_MOMENTA: NOT ON SHELL', i
          write (*, *) 'MSQ before', dot(p(0, i), p(0, i))
          write (*, *) 'MSQ after ', dot(q(0, i), q(0, i))
          stop
        end if
      else
        if (dabs(dot(q(0, i), q(0, i)) - mass_new**2) > &
            1d-3 * max(dot(q(0, i), q(0, i)), q(0, i)**2)) then
          write (*, *) 'ERROR IN CHECK_RESHUFFLED_MOMENTA: NOT ON SHELL', i
          write (*, *) 'MSQ (iresh)', mass_new**2
          write (*, *) 'MSQ after ', dot(q(0, i), q(0, i))
          stop
        end if
      end if
    end do

    do i = 0, 3
      a = 0d0
      b = 0d0
      do j = 1, nexternal - 1
        b = max(b, dabs(q(i, j)))
        if (j <= nincoming) then
          a = a - q(i, j)
        else
          a = a + q(i, j)
        end if
      end do
      if (dabs(a) / b > 1d-6) then
        write (*, *) 'ERROR IN CHECK_RESHUFFLED_MOMENTA: MOM. CONS', &
             i, dabs(a), b
        do j = 1, nexternal - 1
          write (*, *) q(0, j), q(1, j), q(2, j), q(3, j), &
               dsqrt(dot(q(0, j), q(0, j)))
        end do
        stop
      end if
    end do
  end subroutine check_reshuffled_momenta


  subroutine check_reshuffled_momenta_two(p, q, iresh, mass_new)
    double precision, intent(in) :: p(0:3, nexternal - 1)
    double precision, intent(in) :: q(0:3, nexternal - 1)
    integer, intent(in) :: iresh(2)
    double precision, intent(in) :: mass_new(2)
    double precision :: a, b
    integer :: i, j, location

    if (nincoming /= 2) then
      write (*, *) 'ERROR IN OS_CHECK_MOMENTA: nincoming != 2 not ' // &
           'implemented', nincoming
      stop
    end if

    do j = 1, 2
      do i = 1, nexternal - 1
        location = filoc_int(iresh, size(iresh), i)
        if (location == 0) then
          if (dabs(dot(q(0, i), q(0, i)) - dot(p(0, i), p(0, i))) > &
              1d-3 * max(dot(p(0, i), p(0, i)), p(0, i)**2)) then
            write (*, *) &
                 'ERROR IN CHECK_RESHUFFLED_MOMENTA: NOT ON SHELL', i
            write (*, *) 'MSQ before', dot(p(0, i), p(0, i))
            write (*, *) 'MSQ after ', dot(q(0, i), q(0, i))
            stop
          end if
        else
          if (dabs(dot(q(0, i), q(0, i)) - mass_new(location)**2) > &
              1d-3 * max(dot(q(0, i), q(0, i)), q(0, i)**2)) then
            write (*, *) &
                 'ERROR IN CHECK_RESHUFFLED_MOMENTA: NOT ON SHELL', i
            write (*, *) 'MSQ (iresh)', mass_new(j)**2
            write (*, *) 'MSQ after ', dot(q(0, i), q(0, i))
            stop
          end if
        end if
      end do
    end do

    do i = 0, 3
      a = 0d0
      b = 0d0
      do j = 1, nexternal - 1
        b = max(b, dabs(q(i, j)))
        if (j <= nincoming) then
          a = a - q(i, j)
        else
          a = a + q(i, j)
        end if
      end do
      if (dabs(a) / b > 1d-6) then
        write (*, *) 'ERROR IN CHECK_RESHUFFLED_MOMENTA: MOM. CONS', &
             i, dabs(a), b
        do j = 1, nexternal - 1
          write (*, *) q(0, j), q(1, j), q(2, j), q(3, j), &
               dsqrt(dot(q(0, j), q(0, j)))
        end do
        stop
      end if
    end do
  end subroutine check_reshuffled_momenta_two


  integer function filoc_int(array, n, val)
    integer, intent(in) :: n, val
    integer, intent(in) :: array(n)
    integer :: i

    filoc_int = 0
    do i = 1, n
      if (array(i) == val) filoc_int = i
    end do
  end function filoc_int


  subroutine invboostx(p, q, pboost)
    double precision, intent(in) :: p(0:3), q(0:3)
    double precision, intent(out) :: pboost(0:3)
    double precision, parameter :: zero = 0d0
    double precision :: pq, qq, mass, lf

    qq = q(1)**2 + q(2)**2 + q(3)**2
    if (qq /= zero) then
      pq = p(1)*q(1) + p(2)*q(2) + p(3)*q(3)
      mass = dsqrt(max(q(0)**2 - qq, 1d-99))
      lf = (-(q(0) - mass)*pq/qq + p(0)) / mass
      pboost(0) = (p(0)*q(0) - pq) / mass
      pboost(1) = p(1) - q(1)*lf
      pboost(2) = p(2) - q(2)*lf
      pboost(3) = p(3) - q(3)*lf
    else
      pboost = p
    end if
  end subroutine invboostx


  subroutine write_momenta(p)
    double precision, intent(in) :: p(0:3, nexternal)
    integer :: i

    do i = 1, nexternal
      write (*, *) i, p(0, i), p(1, i), p(2, i), p(3, i)
    end do
  end subroutine write_momenta


  subroutine write_momenta4(p)
    double precision, intent(in) :: p(0:4, nexternal)
    integer :: i

    do i = 1, nexternal
      write (*, *) i, p(0, i), p(1, i), p(2, i), p(3, i), p(4, i)
    end do
  end subroutine write_momenta4

end module momentum_reshuffling

module ranmar_module
  implicit none
  private

  integer, parameter :: ranmar_real_kind = selected_real_kind(15, 307)
  integer, parameter :: ranmar_seed_kind = selected_int_kind(18)

  real(kind=ranmar_real_kind) :: random_values(97) = 0.0_ranmar_real_kind
  real(kind=ranmar_real_kind) :: carry = 0.0_ranmar_real_kind
  real(kind=ranmar_real_kind) :: carry_decrement = 0.0_ranmar_real_kind
  real(kind=ranmar_real_kind) :: carry_modulus = 0.0_ranmar_real_kind
  integer :: first_index = 0
  integer :: second_index = 0

  logical :: ntuple_needs_initialization = .true.
  integer :: ntuple_ij = 0
  integer :: ntuple_kl = 0

  integer(kind=ranmar_seed_kind) :: base_seed_common
  common /to_seed/ base_seed_common
  integer :: random_offset_split_common
  common /c_random_offset_split/ random_offset_split_common

  public :: ntuple

contains

  subroutine ntuple(x, a, b, jconfig)
    ! Front end to RANMAR which lets the caller choose a seed.
    real(kind=ranmar_real_kind), intent(out) :: x
    real(kind=ranmar_real_kind), intent(in) :: a
    real(kind=ranmar_real_kind), intent(in) :: b
    integer, intent(in) :: jconfig

    integer :: ioffset
    integer :: joffset
    integer(kind=ranmar_seed_kind) :: iseed

    if (ntuple_needs_initialization) then
      ntuple_needs_initialization = .false.
      call get_offset(ioffset)

      ! Always read the seed from randinit.  The Python run interface updates
      ! that file at the start of a run, including resumed integrations.
      call get_base(iseed)
      base_seed_common = iseed
      call get_moffset(joffset)

      joffset = joffset * 3157
      iseed = iseed * 31300_ranmar_seed_kind
      ntuple_ij = 1802 + jconfig + int(mod(iseed, 30081_ranmar_seed_kind))
      ntuple_kl = 9373 + int(iseed / 30081_ranmar_seed_kind) + ioffset + joffset

      write (*, '(a,i6,a3,i6,a3,i6)') &
           'Using random seed offsets:', jconfig, ' , ', ioffset, ' , ', joffset
      write (*, *) ' with seed', iseed / 31300_ranmar_seed_kind

      do while (ntuple_ij > 31328)
        ntuple_ij = ntuple_ij - 31328
      end do
      do while (ntuple_kl > 30081)
        ntuple_kl = ntuple_kl - 30081
      end do
      call rmarin(ntuple_ij, ntuple_kl)
    end if

    call ranmar(x)
    do while (x < 1.0e-16_ranmar_real_kind)
      call ranmar(x)
    end do
    x = a + x * (b - a)
  end subroutine ntuple


  subroutine get_base(iseed)
    ! Read the base seed from randinit in this or a parent directory.
    integer(kind=ranmar_seed_kind), intent(out) :: iseed

    integer, parameter :: lun = 22
    character(len=60) :: fname
    logical :: done
    integer :: i
    integer :: io_status
    integer :: level

    fname = 'randinit'
    done = .false.
    level = 1
    do while (.not. done .and. level < 5)
      open(unit=lun, file=fname, status='old', iostat=io_status)
      done = io_status == 0
      level = level + 1
      call prepend_parent_directory(fname)
    end do

    if (done) then
      read(lun, '(a)', iostat=io_status) fname
      if (io_status == 0) then
        i = index(fname, '=')
        if (i > 0) fname = fname(i + 1:)
        read(fname, *, iostat=io_status) iseed
      end if
      close(lun)
      if (io_status == 0) return
    end if

    iseed = 0_ranmar_seed_kind
  end subroutine get_base


  subroutine get_offset(iseed)
    ! Read the subprocess offset from iproc.dat.
    integer, intent(out) :: iseed

    integer, parameter :: lun = 22
    integer :: io_status

    open(unit=lun, file='./iproc.dat', status='old', iostat=io_status)
    if (io_status == 0) then
      read(lun, *, iostat=io_status) iseed
      close(lun)
      if (io_status == 0) return
    end if

    open(unit=lun, file='../iproc.dat', status='old', iostat=io_status)
    if (io_status == 0) then
      read(lun, *, iostat=io_status) iseed
      close(lun)
      if (io_status == 0) return
    end if

    iseed = 0
  end subroutine get_offset


  subroutine get_moffset(iseed)
    ! Read the multi-run offset, falling back to the driver-provided split.
    integer, intent(out) :: iseed

    integer, parameter :: lun = 22
    integer :: io_status
    open(unit=lun, file='./moffset.dat', status='old', iostat=io_status)
    if (io_status == 0) then
      read(lun, *, iostat=io_status) iseed
      if (io_status == 0) write (*, *) 'Got moffset', iseed
      close(lun)
      if (io_status == 0) return
    end if

    iseed = random_offset_split_common
  end subroutine get_moffset


  subroutine ranmar(rvec)
    ! Universal random-number generator by Marsaglia and Zaman.
    real(kind=ranmar_real_kind), intent(out) :: rvec

    real(kind=ranmar_real_kind) :: uniform_value

    uniform_value = random_values(first_index) - random_values(second_index)
    if (uniform_value < 0.0_ranmar_real_kind) uniform_value = uniform_value + 1.0_ranmar_real_kind
    random_values(first_index) = uniform_value
    first_index = first_index - 1
    second_index = second_index - 1
    if (first_index == 0) first_index = 97
    if (second_index == 0) second_index = 97

    carry = carry - carry_decrement
    if (carry < 0.0_ranmar_real_kind) carry = carry + carry_modulus
    uniform_value = uniform_value - carry
    if (uniform_value < 0.0_ranmar_real_kind) uniform_value = uniform_value + 1.0_ranmar_real_kind
    rvec = uniform_value
  end subroutine ranmar


  subroutine rmarin(ij, kl)
    ! Initialize RANMAR.  Valid inputs are 0 <= ij <= 31328 and
    ! 0 <= kl <= 30081.
    integer, intent(in) :: ij
    integer, intent(in) :: kl

    character(len=30) :: filename
    logical :: file_exists
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: ii
    integer :: jj
    real(kind=ranmar_real_kind) :: s
    real(kind=ranmar_real_kind) :: t

    write (*, *) 'Ranmar initialization seeds', ij, kl

    if (ij < 0 .or. ij > 31328 .or. kl < 0 .or. kl > 30081) then
      filename = '../../error'
      inquire(file='../../RunWeb', exist=file_exists)
      if (.not. file_exists) call prepend_parent_directory(filename)
      open(unit=26, file=filename, status='unknown')
      if (ij < 0 .or. ij > 31328) then
        write (26, *) 'Bad initialization value of ij in rmarin ', ij
        write (*, *) 'Bad initialization value of ij in rmarin ', ij
      else if (kl < 0 .or. kl > 30081) then
        write (26, *) 'Bad initialization value of kl in rmarin ', kl
        write (*, *) 'Bad initialization value of kl in rmarin ', kl
      end if
      stop
    end if

    i = mod(ij / 177, 177) + 2
    j = mod(ij, 177) + 2
    k = mod(kl / 169, 178) + 1
    l = mod(kl, 169)
    do ii = 1, 97
      s = 0.0_ranmar_real_kind
      t = 0.5_ranmar_real_kind
      do jj = 1, 24
        m = mod(mod(i * j, 179) * k, 179)
        i = j
        j = k
        k = m
        l = mod(53 * l + 1, 169)
        if (mod(l * m, 64) >= 32) s = s + t
        t = 0.5_ranmar_real_kind * t
      end do
      random_values(ii) = s
    end do

    carry = 362436.0_ranmar_real_kind / 16777216.0_ranmar_real_kind
    carry_decrement = 7654321.0_ranmar_real_kind / 16777216.0_ranmar_real_kind
    carry_modulus = 16777213.0_ranmar_real_kind / 16777216.0_ranmar_real_kind
    first_index = 97
    second_index = 33
  end subroutine rmarin


  subroutine prepend_parent_directory(path)
    character(len=*), intent(inout) :: path

    integer :: used_length

    used_length = min(len_trim(path), len(path) - 3)
    if (used_length > 0) path(4:used_length + 3) = path(1:used_length)
    path(1:3) = '../'
    if (used_length + 3 < len(path)) path(used_length + 4:) = ' '
  end subroutine prepend_parent_directory
end module ranmar_module

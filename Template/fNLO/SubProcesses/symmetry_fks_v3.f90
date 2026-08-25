module symmetry_fks_module
  use process_dimensions, only: nexternal, nincoming, ngraphs, &
       fks_configs, validate_process_dimensions
  use fks_metadata, only: validate_fks_metadata, fks_i_d, fks_j_d
  use genps_fks, only: generate_momenta
  use kin_functions_module, only: switchmom => switchmom_impl
  use mint_module, only: ndim, iconfig, ichan, iconfigs, new_point
  use run_state, only: lpp
  use setscales_module, only: set_alphas
  use cuts_module, only: passcuts
  implicit none
  private

  integer, allocatable :: born_mapconfig(:)
  logical :: symmetry_data_initialized = .false.

  public :: initialize_symmetry_data
  public :: finalize_symmetry_data
  public :: run_symmetry

  interface
    subroutine fks_inc_chooser()
    end subroutine fks_inc_chooser

    subroutine leshouche_inc_chooser()
    end subroutine leshouche_inc_chooser

    subroutine setrun()
    end subroutine setrun

    subroutine setpara(card_name)
      character(len=*), intent(in) :: card_name
    end subroutine setpara

    subroutine setcuts()
    end subroutine setcuts

    subroutine sync_cuts_bridge_state()
    end subroutine sync_cuts_bridge_state

    subroutine printout()
    end subroutine printout

    subroutine run_printout()
    end subroutine run_printout

    subroutine fill_configurations_common()
    end subroutine fill_configurations_common

    subroutine setfksfactor()
    end subroutine setfksfactor

    double precision function ran2()
    end function ran2

    subroutine set_cms_stuff(iconfiguration)
      integer, intent(in) :: iconfiguration
    end subroutine set_cms_stuff

    subroutine sborn(p, wgt)
      double precision, intent(in) :: p(0:3, *)
      complex(kind=kind(0d0)), intent(out) :: wgt(2)
    end subroutine sborn

    subroutine open_bash_file(unit_number, file_name, name_length)
      integer, intent(in) :: unit_number
      character(len=*), intent(in) :: file_name
      integer, intent(in) :: name_length
    end subroutine open_bash_file

    subroutine close_bash_file(unit_number)
      integer, intent(in) :: unit_number
    end subroutine close_bash_file
  end interface

contains

  subroutine initialize_symmetry_data(mapconfig_input)
    implicit none
    integer, intent(in) :: mapconfig_input(0:)

    call validate_process_dimensions(require_born=.true.)
    call validate_fks_metadata()
    if (size(mapconfig_input) < 2) then
      call fail_symmetry('the Born configuration table is empty')
    end if
    if (mapconfig_input(0) < 1 .or. &
        mapconfig_input(0) > ubound(mapconfig_input, 1)) then
      call fail_symmetry('the Born configuration count is invalid')
    end if
    if (any(mapconfig_input(1:mapconfig_input(0)) < 1) .or. &
        any(mapconfig_input(1:mapconfig_input(0)) > ngraphs)) then
      call fail_symmetry('a Born configuration index is invalid')
    end if

    if (symmetry_data_initialized) then
      if (size(born_mapconfig) /= size(mapconfig_input) .or. &
          any(born_mapconfig /= mapconfig_input)) then
        call fail_symmetry('Born configurations changed after initialization')
      end if
      return
    end if

    allocate(born_mapconfig(0:ubound(mapconfig_input, 1)))
    born_mapconfig = mapconfig_input
    symmetry_data_initialized = .true.
  end subroutine initialize_symmetry_data


  subroutine finalize_symmetry_data()
    implicit none

    if (allocated(born_mapconfig)) deallocate(born_mapconfig)
    symmetry_data_initialized = .false.
  end subroutine finalize_symmetry_data


  subroutine run_symmetry(run_mode, nndim, amp2, p_born, &
       nfksprocess, calculated_born, multi_channel, nbody, is_aorg, &
       idup)
    implicit none
    character(len=*), intent(in) :: run_mode
    integer, intent(inout) :: nndim
    double precision, intent(inout) :: amp2(:)
    double precision, intent(inout) :: p_born(0:, :)
    integer, intent(inout) :: nfksprocess
    logical, intent(inout) :: calculated_born
    logical, intent(inout) :: multi_channel
    logical, intent(inout) :: nbody
    logical, intent(inout) :: is_aorg(:)
    integer, intent(inout) :: idup(:, :)
    logical :: mtc, even, force_one_job
    integer :: j, k, nmatch, ibase, ntry
    integer, allocatable :: icb(:), inverse_permutation(:)
    integer, allocatable :: use_config(:)
    double precision :: diff, rwgt, wgt
    double precision, allocatable :: p(:, :), x(:)
    double precision, allocatable :: p_born1(:, :), p_born_save(:, :)
    double precision, allocatable :: saved_amplitudes(:)
    complex(kind=kind(0d0)) :: wgt1(2)

    call validate_symmetry_state(amp2, p_born, is_aorg, idup)
    if (len_trim(run_mode) < 2) then
      call fail_symmetry('unknown run_mode in gensym')
    end if
    if (run_mode(1:3) == 'NLO' .or. run_mode(1:2) == 'LO') then
      force_one_job = .false.
    else
      call fail_symmetry('unknown run_mode in gensym')
    end if

    write (*,*) 'run_mode given is: ', run_mode
    allocate(icb(nexternal - 1), inverse_permutation(nexternal))
    allocate(p(0:3, nexternal))
    allocate(p_born1(0:3, nexternal - 1))
    allocate(p_born_save(0:3, nexternal - 1))
    allocate(saved_amplitudes(ngraphs))
    allocate(use_config(0:born_mapconfig(0)))

    multi_channel = .true.
    nbody = .true.

    do nfksprocess = 1, fks_configs
      call fks_inc_chooser()
      if (is_aorg(fks_i_d(nfksprocess))) exit
    end do
    if (nfksprocess > fks_configs) nfksprocess = 1
    call leshouche_inc_chooser()
    call setrun()
    call setpara('param_card.dat')
    call setcuts()
    call sync_cuts_bridge_state()
    call printout()
    call run_printout()
    call fill_configurations_common()
    iconfig = 1
    ichan = 1
    iconfigs(1) = iconfig
    call setfksfactor()

    ndim = 3 * (nexternal - nincoming) - 4
    if (abs(lpp(1)) >= 1) ndim = ndim + 1
    if (abs(lpp(2)) >= 1) ndim = ndim + 1
    if (ndim < 1) call fail_symmetry('invalid integration dimension')
    allocate(x(ndim))
    nndim = ndim
    use_config = 1
    use_config(0) = 0

    ntry = 1
    do j = 1, ndim
      x(j) = ran2()
    end do
    new_point = .true.
    wgt = 1d0
    call generate_momenta(ndim, iconfig, wgt, x, p)
    call set_cms_stuff(-100)
    do while ((.not. passcuts(p, rwgt) .or. wgt < 0d0 .or. &
         p(0, 1) <= 0d0 .or. p_born(0, 1) <= 0d0) .and. &
         ntry < 10000)
      do j = 1, ndim
        x(j) = ran2()
      end do
      new_point = .true.
      wgt = 1d0
      call generate_momenta(ndim, iconfig, wgt, x, p)
      call set_cms_stuff(-100)
      ntry = ntry + 1
    end do
    write (*,*) 'ntry', ntry
    call set_alphaS(p)

    calculated_born = .false.
    call sborn(p_born, wgt1)
    call sborn(p_born, wgt1)
    do j = 1, born_mapconfig(0)
      saved_amplitudes(born_mapconfig(j)) = amp2(born_mapconfig(j))
    end do
    write (*,*) 'born momenta'
    p_born_save = p_born
    do j = 1, nexternal - 1
      write (*,'(i4,4e15.5)') j, p_born(:, j)
    end do

    do k = 1, nexternal - 1
      icb(k) = k
    end do
    nmatch = 0
    mtc = .false.
    call nexper(nexternal - 3, icb(3:), mtc, even)
    do while (mtc)
      call nexper(nexternal - 3, icb(3:), mtc, even)
      icb(3:nexternal - 1) = icb(3:nexternal - 1) + 2
      if (check_swap(icb, idup)) then
        write (*,*) 'Good swap', icb
        call switchmom(p_born_save, p_born1, icb, &
             inverse_permutation, nexternal - 1)
        p_born = p_born1
        calculated_born = .false.
        call sborn(p_born, wgt1)
        do j = 2, born_mapconfig(0)
          do k = 1, j - 1
            diff = abs((amp2(born_mapconfig(j)) - &
                 saved_amplitudes(born_mapconfig(k))) / &
                 (amp2(born_mapconfig(j)) + 1d-99))
            if (diff > 1d-8) cycle
            if (use_config(j) < 0) exit
            nmatch = nmatch + 1
            if (use_config(k) > 0) then
              use_config(k) = use_config(k) + use_config(j)
              use_config(j) = -k
            else
              ibase = -use_config(k)
              use_config(ibase) = use_config(ibase) + use_config(j)
              use_config(j) = -ibase
            end if
          end do
        end do
      else
        write (*,*) 'Bad swap', icb
      end if
      icb(3:nexternal - 1) = icb(3:nexternal - 1) - 2
    end do
    write (*,*) 'Found ', nmatch, ' matches. ', &
         born_mapconfig(0) - nmatch, ' channels remain for integration.'
    call write_bash(use_config, force_one_job)
  end subroutine run_symmetry


  logical function check_swap(permutation, idup)
    implicit none
    integer, intent(in) :: permutation(:)
    integer, intent(in) :: idup(:, :)
    integer :: particle

    check_swap = .true.
    do particle = 1, nexternal - 1
      if (particle == permutation(particle)) cycle
      if (idup(particle, 1) /= idup(permutation(particle), 1)) then
        check_swap = .false.
        return
      end if
      if (is_fks_particle(particle) .or. &
          is_fks_particle(permutation(particle))) then
        check_swap = .false.
        return
      end if
    end do
  end function check_swap


  logical function is_fks_particle(particle)
    implicit none
    integer, intent(in) :: particle
    integer :: configuration

    is_fks_particle = .false.
    do configuration = 1, fks_configs
      if (particle == fks_i_d(configuration) .or. &
          particle == fks_j_d(configuration)) then
        is_fks_particle = .true.
        return
      end if
    end do
  end function is_fks_particle


  subroutine nexper(n, values, mtc, even)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: values(*)
    logical, intent(inout) :: mtc
    logical, intent(inout) :: even
    integer :: difference, i, i1, ia, j, l, m, nm3, sum_value

    if (mtc) goto 10
    nm3 = n - 3
    do i = 1, n
      values(i) = i
    end do
    mtc = .true.
    even = .true.
    if (n == 1) goto 8
6   if (values(n) /= 1 .or. &
        values(1) /= 2 + mod(n, 2)) return
    if (n <= 3) goto 8
    do i = 1, nm3
      if (values(i + 1) /= values(i) + 1) return
    end do
8   mtc = .false.
    return
10  if (n == 1) goto 27
    if (.not. even) goto 20
    ia = values(1)
    values(1) = values(2)
    values(2) = ia
    even = .false.
    goto 6
20  sum_value = 0
    do i1 = 2, n
      ia = values(i1)
      i = i1 - 1
      difference = 0
      do j = 1, i
        if (values(j) > ia) difference = difference + 1
      end do
      sum_value = difference + sum_value
      if (difference /= i * mod(sum_value, 2)) goto 35
    end do
27  values(1) = 0
    goto 8
35  m = mod(sum_value + 1, 2) * (n + 1)
    l = 1
    do j = 1, i
      if (isign(1, values(j) - ia) == &
          isign(1, values(j) - m)) cycle
      m = values(j)
      l = j
    end do
    values(l) = ia
    values(i1) = m
    even = .true.
  end subroutine nexper


  subroutine write_bash(use_config, force_one_job)
    implicit none
    integer, intent(in) :: use_config(0:)
    logical, intent(in) :: force_one_job
    integer :: configuration, name_length, record_length
    character(len=30) :: file_name
    character(len=2) :: postfix
    character(len=14 * max(1, size(use_config) - 1)) :: channel_record
    logical :: initial_fks_partner, final_fks_partner, two_jobs

    initial_fks_partner = .false.
    final_fks_partner = .false.
    do configuration = 1, fks_configs
      if (fks_j_d(configuration) <= nincoming) then
        initial_fks_partner = .true.
      end if
      if (fks_j_d(configuration) > nincoming) then
        final_fks_partner = .true.
      end if
    end do
    two_jobs = .not. force_one_job .and. initial_fks_partner .and. &
         final_fks_partner

    file_name = 'ajob'
    name_length = 4
    call open_bash_file(26, file_name, name_length)
    call close_bash_file(26)
    channel_record = ' '
    record_length = 0
    do configuration = 1, born_mapconfig(0)
      if (use_config(configuration) <= 0) cycle
      if (two_jobs) then
        postfix = '.1'
      else
        postfix = '.0'
      end if
      do
        call append_channel(channel_record, record_length, &
             born_mapconfig(configuration), postfix)
        if (postfix /= '.1') exit
        postfix = '.2'
      end do
    end do
    if (record_length < 1) then
      call fail_symmetry('no independent integration channel was found')
    end if
    open(unit=26, file='channels.txt', status='replace', &
         access='direct', form='formatted', recl=record_length)
    write (26,'(a)', rec=1) channel_record(1:record_length)
    close(26)

    if (born_mapconfig(0) > 9999) then
      write (*,*) 'ERROR: only writing first 9999 jobs', &
           born_mapconfig(0)
      stop 1
    end if

    open(unit=26, file='symfact.dat', status='unknown')
    do configuration = 1, born_mapconfig(0)
      if (use_config(configuration) > 0) then
        if (two_jobs) then
          write (26,'(i6,a2,i6)') born_mapconfig(configuration), &
               '.1', use_config(configuration)
          write (26,'(i6,a2,i6)') born_mapconfig(configuration), &
               '.2', use_config(configuration)
        else
          write (26,'(i6,a2,i6)') born_mapconfig(configuration), &
               '.0', use_config(configuration)
        end if
      else
        if (two_jobs) then
          write (26,'(i6,a2,i6)') born_mapconfig(configuration), &
               '.1', -born_mapconfig(-use_config(configuration))
          write (26,'(i6,a2,i6)') born_mapconfig(configuration), &
               '.2', -born_mapconfig(-use_config(configuration))
        else
          write (26,'(i6,a2,i6)') born_mapconfig(configuration), &
               '.0', -born_mapconfig(-use_config(configuration))
        end if
      end if
    end do
    close(26)
  end subroutine write_bash


  subroutine append_channel(record, record_length, configuration, postfix)
    implicit none
    character(len=*), intent(inout) :: record
    integer, intent(inout) :: record_length
    integer, intent(in) :: configuration
    character(len=2), intent(in) :: postfix
    integer :: first, last

    first = record_length + 1
    if (configuration < 10) then
      last = record_length + 4
      write (record(first:last),'(1x,i1,a2)') &
           configuration, postfix
    else if (configuration < 100) then
      last = record_length + 5
      write (record(first:last),'(1x,i2,a2)') &
           configuration, postfix
    else if (configuration < 1000) then
      last = record_length + 6
      write (record(first:last),'(1x,i3,a2)') &
           configuration, postfix
    else if (configuration < 10000) then
      last = record_length + 7
      write (record(first:last),'(1x,i4,a2)') &
           configuration, postfix
    else
      call fail_symmetry('a channel identifier exceeds four digits')
    end if
    record_length = last
  end subroutine append_channel


  subroutine validate_symmetry_state(amp2, p_born, is_aorg, idup)
    implicit none
    double precision, intent(in) :: amp2(:)
    double precision, intent(in) :: p_born(0:, :)
    logical, intent(in) :: is_aorg(:)
    integer, intent(in) :: idup(:, :)

    call validate_process_dimensions(require_born=.true.)
    call validate_fks_metadata()
    if (.not. symmetry_data_initialized) then
      call fail_symmetry('generated symmetry data are not initialized')
    end if
    if (size(amp2) < ngraphs) then
      call fail_symmetry('the Born amplitude array is too small')
    end if
    if (ubound(p_born, 1) < 3 .or. &
        size(p_born, 2) < nexternal - 1) then
      call fail_symmetry('the Born momentum array has the wrong shape')
    end if
    if (size(is_aorg) < nexternal) then
      call fail_symmetry('the FKS-origin array is too small')
    end if
    if (size(idup, 1) < nexternal .or. size(idup, 2) < 1) then
      call fail_symmetry('the Les Houches identity array is too small')
    end if
  end subroutine validate_symmetry_state


  subroutine fail_symmetry(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*,'(a)') 'symmetry_fks: ' // trim(message)
    stop 1
  end subroutine fail_symmetry

end module symmetry_fks_module

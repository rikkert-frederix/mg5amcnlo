module iproc_map_module
  implicit none
  private

  integer, allocatable, save :: current_ids(:, :)
  integer, allocatable, save :: first_ids(:, :)
  integer, allocatable, save :: equal_to(:)
  integer, allocatable, save :: equal_to_inverse(:)
  integer, save :: first_maxproc = 0

  integer, allocatable, save :: lumi_nproc(:)
  integer, allocatable, save :: lumi_pdgs(:, :, :)
  logical, allocatable, save :: found_appl(:)
  logical, allocatable, save :: found_mg(:)

  logical, parameter :: flavour_map_debug = .false.

  public :: initialize_iproc_map_workspace
  public :: map_iproc_configuration
  public :: finalize_iproc_map_workspace
  public :: initialize_flavour_workspace
  public :: read_initial_states_map
  public :: match_flavour_configuration
  public :: validate_flavour_map
  public :: finalize_flavour_workspace

contains

  subroutine initialize_iproc_map_workspace(nexternal, maxproc)
    implicit none
    integer, intent(in) :: nexternal, maxproc
    integer :: nborn

    call finalize_iproc_map_workspace()

    if (nexternal < 2 .or. maxproc < 1) then
      write(*,*) 'Invalid dimensions in initialize_iproc_map_workspace', &
           nexternal, maxproc
      stop 1
    end if

    nborn = nexternal - 1
    allocate(current_ids(nborn, maxproc))
    allocate(first_ids(nborn, maxproc))
    allocate(equal_to(maxproc))
    allocate(equal_to_inverse(maxproc))

    current_ids = 0
    first_ids = 0
    equal_to = 0
    equal_to_inverse = 0
    first_maxproc = 0
  end subroutine initialize_iproc_map_workspace


  subroutine finalize_iproc_map_workspace()
    implicit none

    if (allocated(current_ids)) deallocate(current_ids)
    if (allocated(first_ids)) deallocate(first_ids)
    if (allocated(equal_to)) deallocate(equal_to)
    if (allocated(equal_to_inverse)) deallocate(equal_to_inverse)
    first_maxproc = 0
  end subroutine finalize_iproc_map_workspace


  subroutine map_iproc_configuration(nfksprocess, qcd_pos, qed_pos, &
       split_type, idup, i_fks, j_fks, iproc, iproc_save, eto, etoi, &
       maxproc_found)
    implicit none
    integer, intent(in) :: nfksprocess, qcd_pos, qed_pos
    logical, intent(in) :: split_type(:)
    integer, intent(in) :: idup(:, :)
    integer, intent(in) :: i_fks, j_fks, iproc
    integer, intent(inout) :: iproc_save(:)
    integer, intent(inout) :: eto(:, :), etoi(:, :)
    integer, intent(out) :: maxproc_found

    integer :: i, ii, j, jj, nequal
    integer :: mother_position, removed_position
    logical :: qcd_split, qed_split

    call check_iproc_workspace(size(idup, 1), size(idup, 2))

    if (nfksprocess < 1 .or. nfksprocess > size(iproc_save)) then
      write(*,*) 'Invalid FKS process in iproc_map:', nfksprocess
      stop 1
    end if
    if (nfksprocess > size(eto, 2) .or. &
         nfksprocess > size(etoi, 2)) then
      write(*,*) 'FKS process exceeds process-map dimensions:', nfksprocess
      stop 1
    end if
    if (iproc < 1 .or. iproc > size(idup, 2) .or. &
         iproc > size(eto, 1) .or. iproc > size(etoi, 1)) then
      write(*,*) 'Invalid IPROC in iproc_map:', iproc
      stop 1
    end if
    if (i_fks < 1 .or. i_fks > size(idup, 1) .or. &
         j_fks < 1 .or. j_fks > size(idup, 1) .or. i_fks == j_fks) then
      write(*,*) 'Invalid FKS pair in iproc_map:', i_fks, j_fks
      stop 1
    end if

    qcd_split = split_order_is_active(qcd_pos, split_type)
    qed_split = split_order_is_active(qed_pos, split_type)
    if (qcd_split .and. qed_split) then
      write(*,*) 'IPROCMAP: NOT IMPLEMENTED'
      stop
    end if

    current_ids = 0
    equal_to = 0
    equal_to_inverse = 0
    maxproc_found = 0
    iproc_save(nfksprocess) = iproc
    mother_position = min(j_fks, i_fks)
    removed_position = max(j_fks, i_fks)

    do j = 1, iproc
      do i = 1, size(idup, 1) - 1
        if (i == mother_position) then
          if (abs(idup(i_fks, j)) == abs(idup(j_fks, j))) then
            if (qcd_split) then
              current_ids(i, j) = 21
            else if (qed_split) then
              current_ids(i, j) = 22
            else
              write(*,*) 'No active splitting order in iproc_map', &
                   nfksprocess
              stop
            end if
          else if (abs(idup(i_fks, j)) == 21 .or. &
               idup(i_fks, j) == 22) then
            current_ids(i, j) = idup(j_fks, j)
          else if (idup(j_fks, j) == 21 .or. &
               idup(j_fks, j) == 22) then
            current_ids(i, j) = -idup(i_fks, j)
          else
            write(*,*) 'Error #1 in iproc_map', nfksprocess, &
                 idup(i_fks, j), idup(j_fks, j)
            stop
          end if
        else if (i < removed_position) then
          current_ids(i, j) = idup(i, j)
        else
          current_ids(i, j) = idup(i + 1, j)
        end if
      end do

      if (j == 1) then
        maxproc_found = 1
        equal_to(j) = 1
        equal_to_inverse(1) = j
      else
        nequal = 0
        do jj = 1, maxproc_found
          nequal = count(current_ids(:, j) == &
               current_ids(:, equal_to_inverse(jj)))
          if (nequal == size(idup, 1) - 1) then
            equal_to(j) = jj
            exit
          end if
        end do
        if (nequal /= size(idup, 1) - 1) then
          maxproc_found = maxproc_found + 1
          equal_to(j) = maxproc_found
          equal_to_inverse(maxproc_found) = j
        end if
      end if
    end do

    if (nfksprocess == 1) then
      first_maxproc = maxproc_found
      do j = 1, iproc
        if (j <= maxproc_found) then
          first_ids(:, j) = current_ids(:, equal_to_inverse(j))
          eto(j, nfksprocess) = equal_to(j)
          etoi(j, nfksprocess) = equal_to_inverse(j)
        else
          eto(j, nfksprocess) = equal_to(j)
        end if
      end do
    else
      if (maxproc_found /= first_maxproc) then
        write(*,*) 'Number of unique IPROCs not identical among ', &
             'nFKSprocesses', nfksprocess, maxproc_found, first_maxproc
        stop
      end if

      do j = 1, iproc
        do jj = 1, maxproc_found
          nequal = count(current_ids(:, j) == first_ids(:, jj))
          if (nequal == size(idup, 1) - 1) then
            eto(j, nfksprocess) = jj
            etoi(jj, nfksprocess) = j
          end if
        end do
      end do

      do j = 1, maxproc_found
        do i = 1, size(idup, 1) - 1
          if (current_ids(i, etoi(j, nfksprocess)) /= first_ids(i, j)) then
            write(*,*) 'Particle IDs not equal (inverse)', j, &
                 nfksprocess, maxproc_found, iproc
            do jj = 1, maxproc_found
              write(*,*) jj, etoi(jj, nfksprocess), ' current:', &
                   (current_ids(ii, etoi(jj, nfksprocess)), &
                    ii = 1, size(idup, 1) - 1)
              write(*,*) jj, jj, ' saved  :', &
                   (first_ids(ii, jj), ii = 1, size(idup, 1) - 1)
            end do
            stop
          end if
        end do
      end do

      do j = 1, iproc
        do i = 1, size(idup, 1) - 1
          if (current_ids(i, j) /= first_ids(i, eto(j, nfksprocess))) then
            write(*,*) 'Particle IDs not equal', j, nfksprocess, &
                 maxproc_found, iproc
            do jj = 1, iproc
              write(*,*) jj, jj, ' current:', &
                   (current_ids(ii, jj), ii = 1, size(idup, 1) - 1)
              write(*,*) jj, jj, ' saved  :', &
                   (first_ids(ii, eto(jj, nfksprocess)), &
                    ii = 1, size(idup, 1) - 1)
            end do
            stop
          end if
        end do
      end do
    end if

    if (nfksprocess == 1) write(*,*) '================================'
    if (nfksprocess == 1) then
      write(*,*) 'process combination map (specified per FKS dir):'
    end if
    write(*, '(i3)', advance='no') nfksprocess
    write(*, '(a)', advance='no') ' map     '
    do j = 1, iproc
      write(*, '(i4)', advance='no') eto(j, nfksprocess)
    end do
    write(*, '(a)') ''
    write(*, '(i3)', advance='no') nfksprocess
    write(*, '(a)', advance='no') ' inv. map'
    do j = 1, maxproc_found
      write(*, '(i4)', advance='no') etoi(j, nfksprocess)
    end do
    write(*, '(a)') ''
    if (nfksprocess == size(iproc_save)) then
      write(*,*) '================================'
    end if
  end subroutine map_iproc_configuration


  subroutine check_iproc_workspace(nexternal, maxproc)
    implicit none
    integer, intent(in) :: nexternal, maxproc

    if (.not. allocated(current_ids)) then
      write(*,*) 'IPROC map workspace has not been initialized'
      stop 1
    end if
    if (size(current_ids, 1) /= nexternal - 1 .or. &
         size(current_ids, 2) /= maxproc) then
      write(*,*) 'IPROC map workspace has inconsistent dimensions'
      stop 1
    end if
  end subroutine check_iproc_workspace


  logical function split_order_is_active(position, split_type)
    implicit none
    integer, intent(in) :: position
    logical, intent(in) :: split_type(:)

    split_order_is_active = .false.
    if (position >= 1 .and. position <= size(split_type)) then
      split_order_is_active = split_type(position)
    end if
  end function split_order_is_active


  subroutine initialize_flavour_workspace(mxpdflumi, max_nproc, &
       maxproc)
    implicit none
    integer, intent(in) :: mxpdflumi, max_nproc, maxproc

    call finalize_flavour_workspace()

    if (mxpdflumi < 1 .or. max_nproc < 1 .or. maxproc < 1) then
      write(*,*) 'Invalid dimensions in initialize_flavour_map_workspace', &
           mxpdflumi, max_nproc, maxproc
      stop 1
    end if

    allocate(lumi_nproc(mxpdflumi))
    allocate(lumi_pdgs(2, max_nproc, mxpdflumi))
    allocate(found_appl(max_nproc))
    allocate(found_mg(maxproc))

    lumi_nproc = 0
    lumi_pdgs = 0
    found_appl = .false.
    found_mg = .false.
  end subroutine initialize_flavour_workspace


  subroutine finalize_flavour_workspace()
    implicit none

    if (allocated(lumi_nproc)) deallocate(lumi_nproc)
    if (allocated(lumi_pdgs)) deallocate(lumi_pdgs)
    if (allocated(found_appl)) deallocate(found_appl)
    if (allocated(found_mg)) deallocate(found_mg)
  end subroutine finalize_flavour_workspace


  subroutine read_initial_states_map(appl_lumimap, appl_nproc, &
       appl_nlumi, npdflumi)
    implicit none
    integer, intent(out) :: appl_lumimap(:, :, :)
    integer, intent(out) :: appl_nproc(:)
    integer, intent(out) :: appl_nlumi, npdflumi

    character(len=200) :: buffer
    integer :: i, ilumi, ios, j, kpdflumi, ncomponents

    call check_flavour_workspace(size(appl_nproc), &
         size(appl_lumimap, 2), size(found_mg))
    if (size(appl_lumimap, 1) /= 2 .or. &
         size(appl_lumimap, 3) /= size(appl_nproc)) then
      write(*,*) 'Inconsistent PineAPPL luminosity-map dimensions'
      stop 1
    end if

    lumi_nproc = 0
    lumi_pdgs = 0
    appl_lumimap = 0
    appl_nproc = 0
    appl_nlumi = 0
    npdflumi = 0

    open(unit=71, status='old', file='initial_states_map.dat')
    do
      read(71, '(a)', iostat=ios) buffer
      if (ios /= 0) exit
      write(*,*) buffer

      read(buffer, *, iostat=ios) kpdflumi, ncomponents
      if (ios /= 0) exit
      if (kpdflumi < 1 .or. kpdflumi > size(lumi_nproc)) then
        write(*,*) 'ERROR in iproc_map.f90, too many PDF luminosities'
        write(*,*) 'increase mxpdflumi in pineappl_maxproc.inc'
        write(*,*) 'and __max_nproc__ in pineappl_interface.cc'
        write(*,*) 'Make sure to assign all variables the same value!'
        stop 1
      end if
      if (ncomponents < 0 .or. ncomponents > size(lumi_pdgs, 2)) then
        write(*,*) 'ERROR in iproc_map.f90, too many processes:', &
             ncomponents
        write(*,*) 'increase max_nproc in pineappl_maxproc.inc'
        write(*,*) 'and __max_nproc__ in pineappl_interface.cc'
        write(*,*) 'Make sure to assign all variables the same value!'
        stop 1
      end if

      lumi_nproc(kpdflumi) = ncomponents
      if (ncomponents > 0) then
        read(buffer, *, iostat=ios) ilumi, lumi_nproc(kpdflumi), &
             ((lumi_pdgs(i, j, kpdflumi), i = 1, 2), &
              j = 1, ncomponents)
        if (ios /= 0 .or. ilumi /= kpdflumi) then
          write(*,*) 'Malformed entry in initial_states_map.dat:', buffer
          stop 1
        end if
      end if
      appl_nproc(kpdflumi) = ncomponents
      appl_nlumi = kpdflumi
    end do
    close(71)

    if (appl_nlumi < 1) then
      write(*,*) 'No luminosities found in initial_states_map.dat'
      stop 1
    end if

    do ilumi = 1, appl_nlumi
      do j = 1, appl_nproc(ilumi)
        do i = 1, 2
          if (lumi_pdgs(i, j, ilumi) == 21) then
            appl_lumimap(i, j, ilumi) = 0
          else
            appl_lumimap(i, j, ilumi) = lumi_pdgs(i, j, ilumi)
          end if
        end do
      end do
    end do

    if (lumi_nproc(appl_nlumi) == 0) then
      npdflumi = appl_nlumi - 1
    else
      npdflumi = appl_nlumi
    end if

    if (flavour_map_debug) then
      write(6,*) 'kpdflumi = ', appl_nlumi
      write(6,*) 'npdflumi = ', npdflumi
      do kpdflumi = 1, npdflumi
        write(6,*) kpdflumi, lumi_nproc(kpdflumi), &
             ((lumi_pdgs(i, j, kpdflumi), i = 1, 2), &
              j = 1, lumi_nproc(kpdflumi))
      end do
    end if
  end subroutine read_initial_states_map


  subroutine match_flavour_configuration(nfksprocess, npdflumi, niprocs, &
       idup, flavour_map, nmatch_total)
    implicit none
    integer, intent(in) :: nfksprocess, npdflumi, niprocs
    integer, intent(in) :: idup(:, :)
    integer, intent(inout) :: flavour_map(:)
    integer, intent(inout) :: nmatch_total

    integer :: found_a, found_m, kpdflumi, l, ll

    call check_flavour_workspace(size(lumi_nproc), size(lumi_pdgs, 2), &
         size(idup, 2))
    if (nfksprocess < 1 .or. nfksprocess > size(flavour_map)) then
      write(*,*) 'Invalid FKS process in setup_flavourmap:', nfksprocess
      stop 1
    end if
    if (npdflumi < 1 .or. npdflumi > size(lumi_nproc)) then
      write(*,*) 'Invalid luminosity count in setup_flavourmap:', npdflumi
      stop 1
    end if
    if (niprocs < 1 .or. niprocs > size(idup, 2) .or. &
         niprocs > size(found_mg)) then
      write(*,*) 'Invalid subprocess count in setup_flavourmap:', niprocs
      stop 1
    end if
    if (size(idup, 1) < 2) then
      write(*,*) 'Missing incoming flavours in setup_flavourmap'
      stop 1
    end if

    flavour_map(nfksprocess) = 0

    if (flavour_map_debug) then
      write(6,*) 'nFKSprocess = ', nfksprocess
      write(6,*) 'niprocs = ', niprocs
      do l = 1, niprocs
        write(6,*) l, idup(1, l), idup(2, l)
      end do
    end if

    do kpdflumi = 1, npdflumi
      found_appl = .false.
      found_mg = .false.

      do l = 1, lumi_nproc(kpdflumi)
        do ll = 1, niprocs
          if (lumi_pdgs(1, l, kpdflumi) == idup(1, ll) .and. &
               lumi_pdgs(2, l, kpdflumi) == idup(2, ll)) then
            found_appl(l) = .true.
            found_mg(ll) = .true.
          end if
        end do
      end do

      found_a = count(found_appl(1:lumi_nproc(kpdflumi)))
      found_m = count(found_mg(1:niprocs))
      if (found_a == lumi_nproc(kpdflumi) .and. &
           found_m == niprocs) then
        flavour_map(nfksprocess) = kpdflumi
        nmatch_total = nmatch_total + 1
      end if
    end do

    if (nmatch_total /= nfksprocess) then
      write(6,*) 'Problem with setup_flavourmap in iproc_map.f90'
      write(6,*) 'nFKSprocess = ', nfksprocess
      write(6,*) 'flavour_map(nFKSprocess) = ', &
           flavour_map(nfksprocess)
      stop
    end if
  end subroutine match_flavour_configuration


  subroutine validate_flavour_map(flavour_map, npdflumi)
    implicit none
    integer, intent(in) :: flavour_map(:), npdflumi
    integer :: nfksprocess

    do nfksprocess = 1, size(flavour_map)
      if (flavour_map(nfksprocess) < 1 .or. &
           flavour_map(nfksprocess) > npdflumi) then
        write(6,*) 'Problem with flavor map, stopping'
        write(6,*) 'flavour_map(nFKSprocess) = ', &
             flavour_map(nfksprocess)
        stop
      end if
    end do

    if (flavour_map_debug) then
      write(*,*) 'flavour map found:'
      do nfksprocess = 1, size(flavour_map)
        write(*,*) nfksprocess, flavour_map(nfksprocess)
      end do
    end if
  end subroutine validate_flavour_map


  subroutine check_flavour_workspace(mxpdflumi, max_nproc, maxproc)
    implicit none
    integer, intent(in) :: mxpdflumi, max_nproc, maxproc

    if (.not. allocated(lumi_nproc) .or. &
         .not. allocated(lumi_pdgs) .or. &
         .not. allocated(found_appl) .or. &
         .not. allocated(found_mg)) then
      write(*,*) 'Flavour-map workspace has not been initialized'
      stop 1
    end if
    if (size(lumi_nproc) /= mxpdflumi .or. &
         size(lumi_pdgs, 2) /= max_nproc .or. &
         size(found_mg) < maxproc) then
      write(*,*) 'Flavour-map workspace has inconsistent dimensions'
      stop 1
    end if
  end subroutine check_flavour_workspace

end module iproc_map_module

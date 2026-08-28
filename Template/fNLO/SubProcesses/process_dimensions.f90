module process_dimensions
  implicit none
  private

  integer, parameter, public :: order_name_length = 64

  integer, public :: nexternal = 0
  integer, public :: event_capacity = 0
  integer, public :: nincoming = 0
  integer, public :: max_particles = 0
  integer, public :: max_branch = 0
  integer, public :: lmaxconfigs = 0
  integer, public :: maxproc = 0
  integer, public :: ngraphs = 0
  integer, public :: ncolor = 0
  integer, public :: maxflow = 0
  integer, public :: fks_configs = 0

  integer, public :: nsplitorders = 0
  integer, public :: qcd_pos = 0
  integer, public :: amp_split_size = 0
  character(len=order_name_length), allocatable, public :: order_names(:)
  integer, allocatable, public :: amp_split_orders(:,:)

  integer, public :: max_bhel = 0
  integer, public :: max_bcol = 0
  integer, public :: maxamps = 0
  integer, public :: born_maxflow = 0
  integer, public :: born_maxproc = 0
  integer, public :: maxsproc = 0
  logical, public :: process_dimensions_initialized = .false.
  logical :: born_dimensions_initialized = .false.

  public :: initialize_process_dimensions
  public :: configure_event_capacity
  public :: initialize_born_dimensions
  public :: validate_process_dimensions
  public :: validate_process_and_born_dimensions

contains

  subroutine initialize_process_dimensions(nexternal_in, nincoming_in, &
       max_particles_in, max_branch_in, lmaxconfigs_in, maxproc_in, &
       ngraphs_in, ncolor_in, maxflow_in, fks_configs_in, &
       nsplitorders_in, qcd_pos_in, amp_split_size_in, &
       order_names_in, amp_split_orders_in)
    implicit none
    integer, intent(in) :: nexternal_in
    integer, intent(in) :: nincoming_in
    integer, intent(in) :: max_particles_in
    integer, intent(in) :: max_branch_in
    integer, intent(in) :: lmaxconfigs_in
    integer, intent(in) :: maxproc_in
    integer, intent(in) :: ngraphs_in
    integer, intent(in) :: ncolor_in
    integer, intent(in) :: maxflow_in
    integer, intent(in) :: fks_configs_in
    integer, intent(in) :: nsplitorders_in
    integer, intent(in) :: qcd_pos_in
    integer, intent(in) :: amp_split_size_in
    character(len=*), intent(in) :: order_names_in(:)
    integer, intent(in) :: amp_split_orders_in(:,:)

    call validate_initial_values(nexternal_in, nincoming_in, &
         max_particles_in, max_branch_in, lmaxconfigs_in, maxproc_in, &
         ngraphs_in, ncolor_in, maxflow_in, fks_configs_in, &
         nsplitorders_in, qcd_pos_in, amp_split_size_in, &
         order_names_in, amp_split_orders_in)

    if (process_dimensions_initialized) then
      if (.not. same_process_dimensions(nexternal_in, nincoming_in, &
           max_particles_in, max_branch_in, lmaxconfigs_in, maxproc_in, &
           ngraphs_in, ncolor_in, maxflow_in, fks_configs_in, &
           nsplitorders_in, qcd_pos_in, amp_split_size_in, &
           order_names_in, amp_split_orders_in)) then
        call fail_validation('process dimensions were reinitialized ' // &
             'with different values')
      end if
      return
    end if

    nexternal = nexternal_in
    event_capacity = nexternal_in
    nincoming = nincoming_in
    max_particles = max_particles_in
    max_branch = max_branch_in
    lmaxconfigs = lmaxconfigs_in
    maxproc = maxproc_in
    ngraphs = ngraphs_in
    ncolor = ncolor_in
    maxflow = maxflow_in
    fks_configs = fks_configs_in
    nsplitorders = nsplitorders_in
    qcd_pos = qcd_pos_in
    amp_split_size = amp_split_size_in

    allocate(order_names(nsplitorders))
    allocate(amp_split_orders(amp_split_size, nsplitorders))
    order_names = order_names_in
    amp_split_orders = amp_split_orders_in

    process_dimensions_initialized = .true.
    call validate_process_dimensions()
  end subroutine initialize_process_dimensions


  subroutine configure_event_capacity(capacity)
    integer, intent(in) :: capacity

    call validate_process_dimensions()
    if (capacity < nexternal) then
      call fail_validation('the event capacity is smaller than NEXTERNAL')
    end if
    event_capacity = capacity
  end subroutine configure_event_capacity


  subroutine initialize_born_dimensions(max_bhel_in, max_bcol_in, &
       maxamps_in, maxflow_in, maxproc_in, maxsproc_in)
    implicit none
    integer, intent(in) :: max_bhel_in
    integer, intent(in) :: max_bcol_in
    integer, intent(in) :: maxamps_in
    integer, intent(in) :: maxflow_in
    integer, intent(in) :: maxproc_in
    integer, intent(in) :: maxsproc_in

    if (max_bhel_in < 1) then
      call fail_validation('MAX_BHEL must be positive')
    end if
    if (max_bcol_in < 1) then
      call fail_validation('MAX_BCOL must be positive')
    end if
    if (maxamps_in < 1) then
      call fail_validation('MAXAMPS must be positive')
    end if
    if (maxflow_in < 1) then
      call fail_validation('the Born MAXFLOW must be positive')
    end if
    if (maxproc_in < 1) then
      call fail_validation('the Born MAXPROC must be positive')
    end if
    if (maxsproc_in < 1 .or. maxsproc_in > maxproc_in) then
      call fail_validation('MAXSPROC is outside the Born process range')
    end if
    if (born_dimensions_initialized) then
      if (max_bhel /= max_bhel_in .or. max_bcol /= max_bcol_in .or. &
           maxamps /= maxamps_in .or. born_maxflow /= maxflow_in .or. &
           born_maxproc /= maxproc_in .or. maxsproc /= maxsproc_in) then
        call fail_validation('Born dimensions were reinitialized ' // &
             'with different values')
      end if
      return
    end if

    max_bhel = max_bhel_in
    max_bcol = max_bcol_in
    maxamps = maxamps_in
    born_maxflow = maxflow_in
    born_maxproc = maxproc_in
    maxsproc = maxsproc_in
    born_dimensions_initialized = .true.

    call validate_born_dimensions()
  end subroutine initialize_born_dimensions


  subroutine validate_process_dimensions()
    implicit none

    if (.not. process_dimensions_are_valid()) then
      call fail_validation('the process dimensions are not valid')
    end if
    if (born_dimensions_initialized) call validate_born_dimensions()
  end subroutine validate_process_dimensions


  subroutine validate_process_and_born_dimensions()
    implicit none

    call validate_process_dimensions()
    call validate_born_dimensions()
  end subroutine validate_process_and_born_dimensions


  logical function process_dimensions_are_valid()
    implicit none

    process_dimensions_are_valid = process_dimensions_initialized
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         nexternal > 0
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         event_capacity >= nexternal
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         nincoming >= 1 .and. nincoming <= 2
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         nincoming <= nexternal
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         max_particles >= nexternal
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         max_branch == max_particles - 1
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         lmaxconfigs > 0 .and. maxproc > 0 .and. ngraphs > 0
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         ncolor > 0 .and. maxflow > 0 .and. fks_configs > 0
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         nsplitorders > 0 .and. amp_split_size > 0
    process_dimensions_are_valid = process_dimensions_are_valid .and. &
         allocated(order_names) .and. allocated(amp_split_orders)
    if (allocated(order_names)) then
      process_dimensions_are_valid = process_dimensions_are_valid .and. &
           size(order_names) == nsplitorders
    end if
    if (allocated(amp_split_orders)) then
      process_dimensions_are_valid = process_dimensions_are_valid .and. &
           size(amp_split_orders, 1) == amp_split_size .and. &
           size(amp_split_orders, 2) == nsplitorders
    end if
  end function process_dimensions_are_valid


  subroutine validate_initial_values(nexternal_in, nincoming_in, &
       max_particles_in, max_branch_in, lmaxconfigs_in, maxproc_in, &
       ngraphs_in, ncolor_in, maxflow_in, fks_configs_in, &
       nsplitorders_in, qcd_pos_in, amp_split_size_in, &
       order_names_in, amp_split_orders_in)
    implicit none
    integer, intent(in) :: nexternal_in
    integer, intent(in) :: nincoming_in
    integer, intent(in) :: max_particles_in
    integer, intent(in) :: max_branch_in
    integer, intent(in) :: lmaxconfigs_in
    integer, intent(in) :: maxproc_in
    integer, intent(in) :: ngraphs_in
    integer, intent(in) :: ncolor_in
    integer, intent(in) :: maxflow_in
    integer, intent(in) :: fks_configs_in
    integer, intent(in) :: nsplitorders_in
    integer, intent(in) :: qcd_pos_in
    integer, intent(in) :: amp_split_size_in
    character(len=*), intent(in) :: order_names_in(:)
    integer, intent(in) :: amp_split_orders_in(:,:)
    integer :: i

    if (nexternal_in < 1) then
      call fail_validation('NEXTERNAL must be positive')
    end if
    if (nincoming_in < 1 .or. nincoming_in > 2 .or. &
         nincoming_in > nexternal_in) then
      call fail_validation('NINCOMING is outside the supported range')
    end if
    if (max_particles_in < nexternal_in) then
      call fail_validation('MAX_PARTICLES is smaller than NEXTERNAL')
    end if
    if (max_branch_in /= max_particles_in - 1) then
      call fail_validation('MAX_BRANCH is inconsistent with MAX_PARTICLES')
    end if
    if (lmaxconfigs_in < 1 .or. maxproc_in < 1 .or. &
         ngraphs_in < 1 .or. ncolor_in < 1 .or. maxflow_in < 1 .or. &
         fks_configs_in < 1) then
      call fail_validation('one or more process dimensions are not positive')
    end if
    if (nsplitorders_in < 1) then
      call fail_validation('NSPLITORDERS must be positive')
    end if
    if (size(order_names_in) /= nsplitorders_in) then
      call fail_validation('the order metadata has inconsistent dimensions')
    end if
    if (len(order_names_in) > order_name_length) then
      call fail_validation('an order name exceeds ORDER_NAME_LENGTH')
    end if
    do i = 1, nsplitorders_in
      if (len_trim(order_names_in(i)) == 0) then
        call fail_validation('an order name is empty')
      end if
    end do
    if (qcd_pos_in < 1 .or. qcd_pos_in > nsplitorders_in) then
      call fail_validation('QCD_POS is outside the order metadata')
    end if
    if (amp_split_size_in < 1) then
      call fail_validation('AMP_SPLIT_SIZE must be positive')
    end if
    if (size(amp_split_orders_in, 1) /= amp_split_size_in .or. &
         size(amp_split_orders_in, 2) /= nsplitorders_in) then
      call fail_validation('AMP_SPLIT_ORDERS has inconsistent dimensions')
    end if
    if (any(amp_split_orders_in < 0)) then
      call fail_validation('AMP_SPLIT_ORDERS contains a negative order')
    end if
  end subroutine validate_initial_values


  logical function same_process_dimensions(nexternal_in, nincoming_in, &
       max_particles_in, max_branch_in, lmaxconfigs_in, maxproc_in, &
       ngraphs_in, ncolor_in, maxflow_in, fks_configs_in, &
       nsplitorders_in, qcd_pos_in, amp_split_size_in, &
       order_names_in, amp_split_orders_in)
    implicit none
    integer, intent(in) :: nexternal_in
    integer, intent(in) :: nincoming_in
    integer, intent(in) :: max_particles_in
    integer, intent(in) :: max_branch_in
    integer, intent(in) :: lmaxconfigs_in
    integer, intent(in) :: maxproc_in
    integer, intent(in) :: ngraphs_in
    integer, intent(in) :: ncolor_in
    integer, intent(in) :: maxflow_in
    integer, intent(in) :: fks_configs_in
    integer, intent(in) :: nsplitorders_in
    integer, intent(in) :: qcd_pos_in
    integer, intent(in) :: amp_split_size_in
    character(len=*), intent(in) :: order_names_in(:)
    integer, intent(in) :: amp_split_orders_in(:,:)

    same_process_dimensions = nexternal == nexternal_in .and. &
         nincoming == nincoming_in .and. &
         max_particles == max_particles_in .and. &
         max_branch == max_branch_in .and. &
         lmaxconfigs == lmaxconfigs_in .and. maxproc == maxproc_in .and. &
         ngraphs == ngraphs_in .and. ncolor == ncolor_in .and. &
         maxflow == maxflow_in .and. fks_configs == fks_configs_in .and. &
         nsplitorders == nsplitorders_in .and. qcd_pos == qcd_pos_in .and. &
         amp_split_size == amp_split_size_in
    same_process_dimensions = same_process_dimensions .and. &
         allocated(order_names) .and. allocated(amp_split_orders)
    if (.not. same_process_dimensions) return
    same_process_dimensions = size(order_names) == nsplitorders_in .and. &
         size(amp_split_orders, 1) == amp_split_size_in .and. &
         size(amp_split_orders, 2) == nsplitorders_in
    if (.not. same_process_dimensions) return
    same_process_dimensions = all(order_names == order_names_in) .and. &
         all(amp_split_orders == amp_split_orders_in)
  end function same_process_dimensions

  subroutine validate_born_dimensions()
    implicit none

    if (.not. born_dimensions_initialized) then
      call fail_validation('the Born dimensions are not initialized')
    end if
    if (max_bhel < 1 .or. max_bcol < 1 .or. maxamps < 1 .or. &
         born_maxflow < 1 .or. born_maxproc < 1 .or. maxsproc < 1 .or. &
         maxsproc > born_maxproc) then
      call fail_validation('the Born dimensions are not valid')
    end if
  end subroutine validate_born_dimensions


  subroutine fail_validation(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*,*) 'Invalid fNLO process dimensions: ', trim(message)
    stop 1
  end subroutine fail_validation

end module process_dimensions

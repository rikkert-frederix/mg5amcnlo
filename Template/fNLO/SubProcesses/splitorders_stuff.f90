module split_orders
  use process_dimensions, only: nsplitorders, qcd_pos, &
       amp_split_size, amp_split_size_born, order_names, &
       amp_split_orders, validate_process_dimensions
  implicit none
  private

  logical :: first_tag_message = .true.
  logical, allocatable :: first_contribution_message(:)
  logical :: split_orders_initialized = .false.

  public :: initialize_split_orders
  public :: finalize_split_orders
  public :: get_orders_tag
  public :: get_orders_tag_from_amp_pos
  public :: orders_to_amp_split_pos
  public :: amp_split_pos_to_orders
  public :: lo_qcd_to_amp_pos
  public :: nlo_qcd_to_amp_pos
  public :: check_amp_split
  public :: orders_equal

contains

  subroutine initialize_split_orders()
    implicit none

    call validate_process_dimensions()
    if (split_orders_initialized) then
      if (.not. allocated(first_contribution_message)) then
        call fail_split_orders('the contribution state is not allocated')
      end if
      if (size(first_contribution_message) /= amp_split_size) then
        call fail_split_orders('the contribution state has the wrong size')
      end if
      return
    end if

    allocate(first_contribution_message(amp_split_size))
    first_contribution_message = .true.
    first_tag_message = .true.
    split_orders_initialized = .true.
  end subroutine initialize_split_orders


  subroutine finalize_split_orders()
    implicit none

    if (allocated(first_contribution_message)) then
      deallocate(first_contribution_message)
    end if
    first_tag_message = .true.
    split_orders_initialized = .false.
  end subroutine finalize_split_orders


  integer function get_orders_tag(ord)
    implicit none
    integer, intent(in) :: ord(:)
    integer :: i
    integer :: position
    integer :: step

    call ensure_split_orders_initialized()
    call validate_order_vector(ord)

    if (first_tag_message) then
      write (*, fmt='(a)', advance='no') &
           'INFO: orders_tag_plot is computed as:'
    end if

    get_orders_tag = 0
    step = 1
    do i = 1, nsplitorders
      if (first_tag_message) then
        write (*, fmt='(3a,i8)', advance='no') &
             '         + ', trim(order_names(i)), ' * ', step
      end if
      get_orders_tag = get_orders_tag + step * ord(i)
      step = step * 100
    end do
    if (first_tag_message) then
      write (*,*)
      first_tag_message = .false.
    end if

    position = orders_to_amp_split_pos(ord)
    if (first_contribution_message(position)) then
      write (*,*) 'orders_tag_plot= ', get_orders_tag, ' for ', &
           (trim(order_names(i)), ',', i=1, nsplitorders), ' = ', &
           (ord(i), ',', i=1, nsplitorders)
      first_contribution_message(position) = .false.
    end if
  end function get_orders_tag


  integer function get_orders_tag_from_amp_pos(iamp)
    implicit none
    integer, intent(in) :: iamp
    integer, allocatable :: orders(:)

    call ensure_split_orders_initialized()
    allocate(orders(nsplitorders))
    call amp_split_pos_to_orders(iamp, orders)
    get_orders_tag_from_amp_pos = get_orders_tag(orders)
    deallocate(orders)
  end function get_orders_tag_from_amp_pos


  integer function orders_to_amp_split_pos(ord)
    implicit none
    integer, intent(in) :: ord(:)
    integer :: i

    call ensure_split_orders_initialized()
    call validate_order_vector(ord)

    do i = 1, amp_split_size
      if (all(amp_split_orders(i,:) == ord)) then
        orders_to_amp_split_pos = i
        return
      end if
    end do

    write (*,*) 'ERROR:: Stopping function orders_to_amp_split_pos'
    write (*,*) 'Could not find orders ', ord
    stop 1
  end function orders_to_amp_split_pos


  subroutine amp_split_pos_to_orders(position, orders)
    implicit none
    integer, intent(in) :: position
    integer, intent(out) :: orders(:)

    call ensure_split_orders_initialized()
    if (size(orders) /= nsplitorders) then
      call fail_split_orders('an order vector has the wrong size')
    end if
    if (position < 1 .or. position > amp_split_size) then
      write (*,*) 'ERROR in amp_split_pos_to_orders'
      write (*,*) 'Invalid position', position, amp_split_size
      stop 1
    end if

    orders = amp_split_orders(position,:)
  end subroutine amp_split_pos_to_orders


  integer function lo_qcd_to_amp_pos(qcdpower)
    implicit none
    integer, intent(in) :: qcdpower
    integer :: position
    integer, allocatable :: orders(:)

    call ensure_qcd_order_available()
    allocate(orders(nsplitorders))
    do position = 1, amp_split_size_born
      call amp_split_pos_to_orders(position, orders)
      if (orders(qcd_pos) == qcdpower) exit
    end do
    lo_qcd_to_amp_pos = position
    deallocate(orders)
  end function lo_qcd_to_amp_pos


  integer function nlo_qcd_to_amp_pos(qcdpower)
    implicit none
    integer, intent(in) :: qcdpower
    integer :: position
    integer, allocatable :: orders(:)

    call ensure_qcd_order_available()
    allocate(orders(nsplitorders))
    do position = amp_split_size_born + 1, amp_split_size
      call amp_split_pos_to_orders(position, orders)
      if (orders(qcd_pos) == qcdpower) exit
    end do
    nlo_qcd_to_amp_pos = position
    deallocate(orders)
  end function nlo_qcd_to_amp_pos


  subroutine check_amp_split()
    implicit none
    integer :: i
    integer :: position
    integer, allocatable :: orders(:)

    call ensure_split_orders_initialized()
    allocate(orders(nsplitorders))
    do i = 1, amp_split_size
      call amp_split_pos_to_orders(i, orders)
      position = orders_to_amp_split_pos(orders)

      if (position /= i) then
        write (*,*) 'ERROR#1 in check amp_split', position, i
        write (*,*) 'ORD is ', orders
        stop 1
      end if

      if (get_orders_tag(orders) /= get_orders_tag_from_amp_pos(i)) then
        write (*,*) 'ERROR#2 in check amp_split', &
             get_orders_tag(orders), get_orders_tag_from_amp_pos(i)
        write (*,*) 'I, ORD ', i, orders
        stop 1
      end if

      write (*,*) 'AMP_SPLIT: ', i, 'correspond to S.O.', orders
    end do
    deallocate(orders)
  end subroutine check_amp_split


  logical function orders_equal(orders1, orders2)
    implicit none
    integer, intent(in) :: orders1(:)
    integer, intent(in) :: orders2(:)

    call ensure_split_orders_initialized()
    call validate_order_vector(orders1)
    call validate_order_vector(orders2)
    orders_equal = all(orders1 == orders2)
  end function orders_equal


  subroutine ensure_split_orders_initialized()
    implicit none

    if (.not. split_orders_initialized) call initialize_split_orders()
  end subroutine ensure_split_orders_initialized


  subroutine ensure_qcd_order_available()
    implicit none

    call ensure_split_orders_initialized()
    if (qcd_pos < 1 .or. qcd_pos > nsplitorders) then
      call fail_split_orders('the process has no QCD split order')
    end if
  end subroutine ensure_qcd_order_available


  subroutine validate_order_vector(orders)
    implicit none
    integer, intent(in) :: orders(:)

    if (size(orders) /= nsplitorders) then
      call fail_split_orders('an order vector has the wrong size')
    end if
  end subroutine validate_order_vector


  subroutine fail_split_orders(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*,*) 'Invalid fNLO split-order state: ', trim(message)
    stop 1
  end subroutine fail_split_orders

end module split_orders

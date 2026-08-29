module multiplicative_lambda_validation
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  ! Bookkeeping for the formal block parameters
  !
  !   B_b(lambda_b) = B_b^(0) + lambda_b Delta B_b^(1).
  !
  ! Every entry accumulated here is still one separately measured
  ! Cartesian momentum atom.  Only its scalar contribution to a given
  ! formal-lambda coefficient is added after cuts and observables have seen
  ! that atom.  In particular, this type must never be used to combine FKS
  ! real and counterevent momenta before measurement.
  type, public :: multiplicative_lambda_accumulator
    integer :: block_count = 0
    double precision :: exact_weight = 0d0
    double precision :: absolute_weight = 0d0
    double precision, allocatable :: order_coefficients(:)
    double precision, allocatable :: single_coefficients(:)
    double precision, allocatable :: pair_coefficients(:, :)
  end type multiplicative_lambda_accumulator

  public :: initialize_multiplicative_lambda_accumulator
  public :: accumulate_multiplicative_lambda_atom
  public :: formal_lambda_global_weight
  public :: formal_lambda_lo_weight
  public :: formal_lambda_linear_correction
  public :: formal_lambda_additive_weight
  public :: formal_lambda_block_linear_correction
  public :: formal_lambda_block_mixed_coefficient
  public :: require_formal_lambda_closure
  public :: require_formal_lambda_linear_closure

contains

  subroutine initialize_multiplicative_lambda_accumulator( &
       accumulator, block_count)
    type(multiplicative_lambda_accumulator), intent(inout) :: accumulator
    integer, intent(in) :: block_count

    if (block_count < 1) then
      call fail_lambda_validation('the block count is not positive')
    end if
    accumulator%block_count = block_count
    accumulator%exact_weight = 0d0
    accumulator%absolute_weight = 0d0
    if (allocated(accumulator%order_coefficients)) then
      if (lbound(accumulator%order_coefficients, 1) /= 0 .or. &
          ubound(accumulator%order_coefficients, 1) /= block_count) &
           deallocate(accumulator%order_coefficients)
    end if
    if (.not. allocated(accumulator%order_coefficients)) &
         allocate(accumulator%order_coefficients(0:block_count))
    if (allocated(accumulator%single_coefficients)) then
      if (size(accumulator%single_coefficients) /= block_count) &
           deallocate(accumulator%single_coefficients)
    end if
    if (.not. allocated(accumulator%single_coefficients)) &
         allocate(accumulator%single_coefficients(block_count))
    if (allocated(accumulator%pair_coefficients)) then
      if (size(accumulator%pair_coefficients, 1) /= block_count .or. &
          size(accumulator%pair_coefficients, 2) /= block_count) &
           deallocate(accumulator%pair_coefficients)
    end if
    if (.not. allocated(accumulator%pair_coefficients)) &
         allocate(accumulator%pair_coefficients(block_count, block_count))
    accumulator%order_coefficients = 0d0
    accumulator%single_coefficients = 0d0
    accumulator%pair_coefficients = 0d0
  end subroutine initialize_multiplicative_lambda_accumulator


  subroutine accumulate_multiplicative_lambda_atom( &
       accumulator, block_orders, weight)
    type(multiplicative_lambda_accumulator), intent(inout) :: accumulator
    integer, intent(in) :: block_orders(:)
    double precision, intent(in) :: weight
    integer :: block, first_block, second_block, order

    call require_accumulator_shape(accumulator)
    if (size(block_orders) /= accumulator%block_count) then
      call fail_lambda_validation( &
           'a block-order vector has the wrong size')
    end if
    if (any(block_orders < 0) .or. any(block_orders > 1)) then
      call fail_lambda_validation( &
           'a block-order vector is not multi-affine')
    end if
    if (.not. ieee_is_finite(weight)) then
      call fail_lambda_validation('an accumulated atom is not finite')
    end if

    order = sum(block_orders)
    accumulator%exact_weight = accumulator%exact_weight + weight
    accumulator%absolute_weight = accumulator%absolute_weight + abs(weight)
    accumulator%order_coefficients(order) = &
         accumulator%order_coefficients(order) + weight

    if (order == 1) then
      do block = 1, accumulator%block_count
        if (block_orders(block) == 0) cycle
        accumulator%single_coefficients(block) = &
             accumulator%single_coefficients(block) + weight
        exit
      end do
    else if (order == 2) then
      first_block = 0
      second_block = 0
      do block = 1, accumulator%block_count
        if (block_orders(block) == 0) cycle
        if (first_block == 0) then
          first_block = block
        else
          second_block = block
          exit
        end if
      end do
      accumulator%pair_coefficients(first_block, second_block) = &
           accumulator%pair_coefficients(first_block, second_block) + weight
      accumulator%pair_coefficients(second_block, first_block) = &
           accumulator%pair_coefficients(first_block, second_block)
    end if
  end subroutine accumulate_multiplicative_lambda_atom


  double precision function formal_lambda_global_weight( &
       accumulator, lambda, lo_widths, nlo_widths)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    double precision, intent(in) :: lambda
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    double precision :: numerator
    integer :: order

    call require_accumulator_shape(accumulator)
    if (.not. ieee_is_finite(lambda)) then
      call fail_lambda_validation('the formal lambda is not finite')
    end if
    numerator = 0d0
    do order = 0, accumulator%block_count
      numerator = numerator + &
           accumulator%order_coefficients(order)*lambda**order
    end do
    formal_lambda_global_weight = numerator* &
         formal_lambda_width_reweight(lambda, lo_widths, nlo_widths)
  end function formal_lambda_global_weight


  double precision function formal_lambda_lo_weight( &
       accumulator, lo_widths, nlo_widths)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)

    call require_accumulator_shape(accumulator)
    formal_lambda_lo_weight = accumulator%order_coefficients(0)* &
         formal_lambda_width_reweight(0d0, lo_widths, nlo_widths)
  end function formal_lambda_lo_weight


  double precision function formal_lambda_linear_correction( &
       accumulator, lo_widths, nlo_widths)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    double precision :: width_slope

    call require_accumulator_shape(accumulator)
    call require_width_shape(lo_widths, nlo_widths)
    width_slope = sum((nlo_widths - lo_widths)/lo_widths)
    formal_lambda_linear_correction = &
         formal_lambda_width_reweight(0d0, lo_widths, nlo_widths)* &
         (accumulator%order_coefficients(1) - &
          width_slope*accumulator%order_coefficients(0))
  end function formal_lambda_linear_correction


  double precision function formal_lambda_additive_weight( &
       accumulator, lo_widths, nlo_widths)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)

    formal_lambda_additive_weight = &
         formal_lambda_lo_weight(accumulator, lo_widths, nlo_widths) + &
         formal_lambda_linear_correction( &
         accumulator, lo_widths, nlo_widths)
  end function formal_lambda_additive_weight


  double precision function formal_lambda_block_linear_correction( &
       accumulator, block, width_blocks, lo_widths, nlo_widths)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    integer, intent(in) :: block
    integer, intent(in) :: width_blocks(:)
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    double precision :: slope

    call require_accumulator_shape(accumulator)
    call require_block(block, accumulator%block_count)
    call require_width_blocks( &
         width_blocks, lo_widths, nlo_widths, accumulator%block_count)
    slope = width_slope_for_block( &
         block, width_blocks, lo_widths, nlo_widths)
    formal_lambda_block_linear_correction = &
         formal_lambda_width_reweight(0d0, lo_widths, nlo_widths)* &
         (accumulator%single_coefficients(block) - &
          slope*accumulator%order_coefficients(0))
  end function formal_lambda_block_linear_correction


  double precision function formal_lambda_block_mixed_coefficient( &
       accumulator, first_block, second_block, width_blocks, &
       lo_widths, nlo_widths)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    integer, intent(in) :: first_block, second_block
    integer, intent(in) :: width_blocks(:)
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    double precision :: first_slope, second_slope

    call require_accumulator_shape(accumulator)
    call require_block(first_block, accumulator%block_count)
    call require_block(second_block, accumulator%block_count)
    if (first_block == second_block) then
      call fail_lambda_validation( &
           'a mixed coefficient needs two distinct blocks')
    end if
    call require_width_blocks( &
         width_blocks, lo_widths, nlo_widths, accumulator%block_count)
    first_slope = width_slope_for_block( &
         first_block, width_blocks, lo_widths, nlo_widths)
    second_slope = width_slope_for_block( &
         second_block, width_blocks, lo_widths, nlo_widths)
    formal_lambda_block_mixed_coefficient = &
         formal_lambda_width_reweight(0d0, lo_widths, nlo_widths)* &
         (accumulator%pair_coefficients(first_block, second_block) - &
          first_slope*accumulator%single_coefficients(second_block) - &
          second_slope*accumulator%single_coefficients(first_block) + &
          first_slope*second_slope* &
          accumulator%order_coefficients(0))
  end function formal_lambda_block_mixed_coefficient


  subroutine require_formal_lambda_closure( &
       accumulator, exact_weight, lo_widths, nlo_widths, tolerance)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    double precision, intent(in) :: exact_weight
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    double precision, intent(in), optional :: tolerance
    double precision :: allowed, reconstructed, relative_tolerance, scale

    call require_accumulator_shape(accumulator)
    relative_tolerance = 1d-11
    if (present(tolerance)) relative_tolerance = tolerance
    if (.not. ieee_is_finite(relative_tolerance) .or. &
        relative_tolerance < 0d0) then
      call fail_lambda_validation('the closure tolerance is invalid')
    end if
    reconstructed = formal_lambda_global_weight( &
         accumulator, 1d0, lo_widths, nlo_widths)
    scale = max(abs(exact_weight), accumulator%absolute_weight)
    allowed = relative_tolerance*scale
    if (abs(accumulator%exact_weight - exact_weight) > allowed) then
      call fail_lambda_validation( &
           'the accumulated exact weight does not close')
    end if
    if (abs(reconstructed - exact_weight) > allowed) then
      call fail_lambda_validation( &
           'the lambda=1 weight does not reproduce the exact product')
    end if
  end subroutine require_formal_lambda_closure


  subroutine require_formal_lambda_linear_closure( &
       accumulator, width_blocks, lo_widths, nlo_widths, tolerance)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator
    integer, intent(in) :: width_blocks(:)
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    double precision, intent(in), optional :: tolerance
    double precision :: allowed, block_sum, global_coefficient
    double precision :: relative_tolerance, scale, value
    integer :: block

    call require_accumulator_shape(accumulator)
    call require_width_blocks( &
         width_blocks, lo_widths, nlo_widths, accumulator%block_count)
    relative_tolerance = 1d-11
    if (present(tolerance)) relative_tolerance = tolerance
    if (.not. ieee_is_finite(relative_tolerance) .or. &
        relative_tolerance < 0d0) then
      call fail_lambda_validation('the linear-closure tolerance is invalid')
    end if
    block_sum = 0d0
    scale = 0d0
    do block = 1, accumulator%block_count
      value = formal_lambda_block_linear_correction( &
           accumulator, block, width_blocks, lo_widths, nlo_widths)
      block_sum = block_sum + value
      scale = scale + abs(value)
    end do
    global_coefficient = formal_lambda_linear_correction( &
         accumulator, lo_widths, nlo_widths)
    ! The global and block sums traverse the same signed atoms in different
    ! orders.  Their physical result can be much smaller than the individual
    ! single-block coefficients and their underlying signed atoms, so use
    ! both as the round-off scale as well as the final result.  Otherwise a
    ! legitimate FKS cancellation can make the nominal relative tolerance
    ! smaller than double-precision summation noise.
    scale = max(scale, abs(global_coefficient), &
                abs(accumulator%order_coefficients(1)), &
                sum(abs(accumulator%single_coefficients)), &
                accumulator%absolute_weight)
    allowed = relative_tolerance*scale
    if (abs(block_sum - global_coefficient) > allowed) then
      call fail_lambda_validation( &
           'the block-linear coefficients do not reproduce the derivative')
    end if
  end subroutine require_formal_lambda_linear_closure


  double precision function formal_lambda_width_reweight( &
       lambda, lo_widths, nlo_widths)
    double precision, intent(in) :: lambda
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    double precision :: interpolated_width
    integer :: occurrence

    call require_width_shape(lo_widths, nlo_widths)
    formal_lambda_width_reweight = 1d0
    do occurrence = 1, size(lo_widths)
      interpolated_width = lo_widths(occurrence) + lambda* &
           (nlo_widths(occurrence) - lo_widths(occurrence))
      if (.not. ieee_is_finite(interpolated_width) .or. &
          interpolated_width <= 0d0) then
        call fail_lambda_validation( &
             'an interpolated physical width is not positive')
      end if
      formal_lambda_width_reweight = formal_lambda_width_reweight* &
           nlo_widths(occurrence)/interpolated_width
    end do
  end function formal_lambda_width_reweight


  double precision function width_slope_for_block( &
       block, width_blocks, lo_widths, nlo_widths)
    integer, intent(in) :: block, width_blocks(:)
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)
    integer :: occurrence

    width_slope_for_block = 0d0
    do occurrence = 1, size(width_blocks)
      if (width_blocks(occurrence) /= block) cycle
      width_slope_for_block = width_slope_for_block + &
           (nlo_widths(occurrence) - lo_widths(occurrence))/ &
           lo_widths(occurrence)
    end do
  end function width_slope_for_block


  subroutine require_accumulator_shape(accumulator)
    type(multiplicative_lambda_accumulator), intent(in) :: accumulator

    if (accumulator%block_count < 1 .or. &
        .not. allocated(accumulator%order_coefficients) .or. &
        .not. allocated(accumulator%single_coefficients) .or. &
        .not. allocated(accumulator%pair_coefficients)) then
      call fail_lambda_validation('the lambda accumulator is uninitialized')
    end if
    if (lbound(accumulator%order_coefficients, 1) /= 0 .or. &
        ubound(accumulator%order_coefficients, 1) /= &
        accumulator%block_count .or. &
        size(accumulator%single_coefficients) /= accumulator%block_count .or. &
        size(accumulator%pair_coefficients, 1) /= accumulator%block_count .or. &
        size(accumulator%pair_coefficients, 2) /= accumulator%block_count) then
      call fail_lambda_validation('the lambda accumulator has invalid storage')
    end if
  end subroutine require_accumulator_shape


  subroutine require_width_shape(lo_widths, nlo_widths)
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)

    if (size(lo_widths) /= size(nlo_widths)) then
      call fail_lambda_validation('the LO/NLO width vectors disagree')
    end if
    if (any(.not. ieee_is_finite(lo_widths)) .or. &
        any(.not. ieee_is_finite(nlo_widths)) .or. &
        any(lo_widths <= 0d0) .or. any(nlo_widths <= 0d0)) then
      call fail_lambda_validation('a physical width is not positive and finite')
    end if
  end subroutine require_width_shape


  subroutine require_width_blocks( &
       width_blocks, lo_widths, nlo_widths, block_count)
    integer, intent(in) :: width_blocks(:), block_count
    double precision, intent(in) :: lo_widths(:), nlo_widths(:)

    call require_width_shape(lo_widths, nlo_widths)
    if (size(width_blocks) /= size(lo_widths)) then
      call fail_lambda_validation('the width-owner vector has the wrong size')
    end if
    if (any(width_blocks < 1) .or. any(width_blocks > block_count)) then
      call fail_lambda_validation('a width owner is not a physical block')
    end if
  end subroutine require_width_blocks


  subroutine require_block(block, block_count)
    integer, intent(in) :: block, block_count

    if (block < 1 .or. block > block_count) then
      call fail_lambda_validation('a validation block is out of range')
    end if
  end subroutine require_block


  subroutine fail_lambda_validation(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_lambda_validation: '//trim(message)
    stop 1
  end subroutine fail_lambda_validation

end module multiplicative_lambda_validation

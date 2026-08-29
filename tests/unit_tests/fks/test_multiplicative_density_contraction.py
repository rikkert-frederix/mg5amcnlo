from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeDensityContractionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.compiler = shutil.which('gfortran')
        cls.repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        cls.sources = [os.path.join(
            cls.repository, 'Template', 'fNLO', 'SubProcesses', name)
            for name in (
                'multiplicative_density_terms.f90',
                'multiplicative_density_contraction.f90')]

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_signed_cartesian_terms_are_contracted_before_summing(self):
        stubs = r'''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 4
contains
  subroutine validate_process_dimensions()
  end subroutine validate_process_dimensions
end module process_dimensions

module fnlo_process_common
  implicit none
  integer, parameter :: soft_counterevent = 0
  integer, parameter :: real_event = 3
end module fnlo_process_common

module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_no_insertion = 0
  integer, parameter :: spin_density_virtual_insertion = 3
  integer, parameter :: spin_density_color_insertion = 4
  integer, parameter :: spin_density_fast_virtual_insertion = 5
end module spin_density_matrix_results

module multiplicative_kinematics
  implicit none
contains
  subroutine realize_factorized_event_tuple(event_slots, pass)
    integer, intent(in) :: event_slots(0:)
    logical, intent(out) :: pass
    pass = .true.
  end subroutine realize_factorized_event_tuple
end module multiplicative_kinematics

module multiplicative_generated_metadata
  implicit none
  integer, parameter :: multiplicative_block_count = 2
  integer, parameter :: multiplicative_physical_blocks(2) = (/0, 1/)
contains
  integer function multiplicative_component_position(block)
    integer, intent(in) :: block
    multiplicative_component_position = block + 1
  end function multiplicative_component_position
end module multiplicative_generated_metadata
'''
        generated = r'''
integer function sdm_multiplicative_block_count()
  sdm_multiplicative_block_count = 2
end function sdm_multiplicative_block_count

integer function sdm_multiplicative_physical_block(position)
  integer, intent(in) :: position
  sdm_multiplicative_physical_block = position - 1
end function sdm_multiplicative_physical_block

integer function sdm_multiplicative_component_position(block)
  integer, intent(in) :: block
  sdm_multiplicative_component_position = block + 1
end function sdm_multiplicative_component_position

module basis_stub
  implicit none
  integer :: saved_slots(2), saved_counts(2), saved_kinds(8,2)
  integer :: evaluate_calls = 0
end module basis_stub

subroutine sdm_multiplicative_prepare_basis( &
     event_slots, maximum_primitives, primitive_counts, insertion_kinds, &
     insertion_ids, insertion_ranks, correlation_legs, include_virtual, &
     precision_asked, precision_found, return_code)
  use basis_stub
  integer, intent(in) :: event_slots(*), maximum_primitives
  integer, intent(in) :: primitive_counts(*)
  integer, intent(in) :: insertion_kinds(maximum_primitives, *)
  integer, intent(in) :: insertion_ids(maximum_primitives, *)
  integer, intent(in) :: insertion_ranks(maximum_primitives, *)
  integer, intent(in) :: correlation_legs(maximum_primitives, *)
  logical, intent(in) :: include_virtual
  double precision, intent(in) :: precision_asked
  double precision, intent(out) :: precision_found
  integer, intent(out) :: return_code
  integer :: position

  saved_slots = event_slots(1:2)
  saved_counts = primitive_counts(1:2)
  saved_kinds = 0
  do position = 1, 2
    saved_kinds(1:saved_counts(position), position) = &
         insertion_kinds(1:saved_counts(position), position)
  end do
  precision_found = precision_asked
  return_code = 0
end subroutine sdm_multiplicative_prepare_basis

subroutine sdm_multiplicative_evaluate_basis( &
     maximum_primitives, primitive_counts, coefficients, result)
  use basis_stub
  integer, intent(in) :: maximum_primitives, primitive_counts(*)
  complex(kind=8), intent(in) :: coefficients(maximum_primitives, *)
  complex(kind=8), intent(out) :: result
  integer :: position, primitive
  complex(kind=8) :: effective

  evaluate_calls = evaluate_calls + 1
  result = (1d0, 0d0)
  do position = 1, 2
    effective = (0d0, 0d0)
    do primitive = 1, primitive_counts(position)
      effective = effective + coefficients(primitive, position)* &
           dble(1 + saved_slots(position) + &
                saved_kinds(primitive, position))
    end do
    result = result*effective
  end do
end subroutine sdm_multiplicative_evaluate_basis
'''
        program = r'''
program test_density_contraction
  use basis_stub
  use multiplicative_density_terms
  use multiplicative_density_contraction
  implicit none
  type(block_nlo_distribution) :: distributions(2)
  type(multiplicative_density_tuple) :: tuple
  type(density_tuple_schedule) :: schedule
  type(multiplicative_density_basis) :: basis
  complex(kind=8) :: coefficients(3), result, direct_result
  double precision :: logs(0:4), varied_logs_r(0:4), varied_logs_f(0:4)
  double precision :: coupling_rescaling(0:4,0:1)
  integer :: distribution
  integer :: term_slots(2), term_signs(2), term_kinds(2)
  double precision :: term_coefficients(2)
  double precision :: expected(4)

  term_slots = (/3, 0/)
  term_signs = (/1, -1/)
  term_kinds = (/2, 1/)
  expected = (/84d0, -180d0, -168d0, 360d0/)
  logs = 0d0
  coupling_rescaling = 1d0
  coupling_rescaling(0, 1) = 2d0
  coupling_rescaling(1, 1) = 3d0
  do distribution = 1, 2
    call initialize_block_distribution( &
         distributions(distribution), distribution - 1, 2)
    term_coefficients = merge((/2d0, 3d0/), (/5d0, 7d0/), &
                              distribution == 1)
    call fill_term(distributions(distribution), 1, term_slots(1), &
                   term_signs(1), term_kinds(1), term_coefficients(1))
    call fill_term(distributions(distribution), 2, term_slots(2), &
                   term_signs(2), term_kinds(2), term_coefficients(2))
    call finalize_block_distribution(distributions(distribution))
  end do

  call initialize_density_tuple_schedule(distributions, schedule)
  do distribution = 1, 4
    call decode_scheduled_density_tuple( &
         distributions, schedule, distribution, tuple)
    call prepare_multiplicative_density_basis( &
         distributions, tuple, 1d-6, basis)
    call evaluate_multiplicative_density_basis( &
         basis, logs, logs, result)
    if (abs(dble(result) - expected(distribution)) > 1d-12) stop 11
    if (abs(aimag(result)) > 1d-12) stop 12
    if (basis%nlo_order /= 2 .or. basis%return_code /= 0) stop 13
    if (abs(basis%precision_found - 1d-6) > 1d-15) stop 14
    call evaluate_multiplicative_density_basis( &
         basis, logs, logs, result, coupling_rescaling)
    if (abs(dble(result) - 6d0*expected(distribution)) > 1d-12) stop 15
  end do

  call decode_scheduled_density_tuple(distributions, schedule, 1, tuple)
  call prepare_multiplicative_density_basis( &
       distributions, tuple, 1d-6, basis, .true.)
  varied_logs_r = 0d0
  varied_logs_f = 0d0
  varied_logs_r(0) = 0.2d0
  varied_logs_r(1) = -0.3d0
  varied_logs_f(0) = 0.4d0
  varied_logs_f(1) = 0.1d0
  evaluate_calls = 0
  call evaluate_multiplicative_density_basis( &
       basis, varied_logs_r, varied_logs_f, direct_result, &
       coupling_rescaling)
  call configure_multiplicative_scale_evaluation(basis, 7)
  if (basis%scale_polynomial_enabled) stop 16
  evaluate_calls = 0
  call evaluate_multiplicative_scale_polynomial( &
       basis, varied_logs_r, varied_logs_f, result, coupling_rescaling)
  if (abs(result - direct_result) > 1d-12 .or. evaluate_calls /= 1) stop 17
  call configure_multiplicative_scale_evaluation(basis, 10)
  if (.not. basis%scale_polynomial_enabled) stop 18
  call evaluate_multiplicative_scale_polynomial( &
       basis, varied_logs_r, varied_logs_f, result, coupling_rescaling)
  if (basis%scale_monomial_count /= 9 .or. evaluate_calls /= 10) stop 19
  call evaluate_multiplicative_scale_polynomial( &
       basis, varied_logs_r, varied_logs_f, result, coupling_rescaling)
  if (evaluate_calls /= 10) stop 20

contains

  subroutine fill_term(distribution, term, event_slot, sign, kind, value)
    type(block_nlo_distribution), intent(inout) :: distribution
    integer, intent(in) :: term, event_slot, sign, kind
    double precision, intent(in) :: value

    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx(value, 0d0, kind=8)
    coefficients(2) = cmplx(0.5d0*value, 0d0, kind=8)
    coefficients(3) = cmplx(-0.25d0*value, 0d0, kind=8)
    call initialize_block_distribution_term( &
         distribution, term, event_slot, sign, 1, 1)
    call set_density_primitive( &
         distribution%terms(term), 1, kind, term, 1, 0, 1, &
         coefficients, .true.)
    call finalize_block_distribution_term(distribution%terms(term))
  end subroutine fill_term
end program test_density_contraction
'''
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for name, source in (
                    ('stubs.f90', stubs), ('generated.f90', generated),
                    ('test.f90', program)):
                path = os.path.join(directory, name)
                with open(path, 'w') as stream:
                    stream.write(source)
                paths.append(path)
            executable = os.path.join(directory, 'test')
            subprocess.check_call([
                self.compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory, paths[0]] +
                self.sources + paths[1:] + ['-o', executable])
            subprocess.check_call([executable])

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_same_momentum_primitives_keep_loop_times_loop_term(self):
        stubs = r'''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 2
contains
  subroutine validate_process_dimensions()
  end subroutine validate_process_dimensions
end module process_dimensions

module fnlo_process_common
  implicit none
  integer, parameter :: soft_counterevent = 0
  integer, parameter :: real_event = 3
end module fnlo_process_common

module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_no_insertion = 0
  integer, parameter :: spin_density_virtual_insertion = 3
  integer, parameter :: spin_density_color_insertion = 4
  integer, parameter :: spin_density_fast_virtual_insertion = 5
end module spin_density_matrix_results

module multiplicative_kinematics
  implicit none
contains
  subroutine realize_factorized_event_tuple(event_slots, pass)
    integer, intent(in) :: event_slots(0:)
    logical, intent(out) :: pass
    pass = .true.
  end subroutine realize_factorized_event_tuple
end module multiplicative_kinematics

module multiplicative_generated_metadata
  implicit none
  integer, parameter :: multiplicative_block_count = 2
  integer, parameter :: multiplicative_physical_blocks(2) = (/0, 1/)
contains
  integer function multiplicative_component_position(block)
    integer, intent(in) :: block
    multiplicative_component_position = block + 1
  end function multiplicative_component_position
end module multiplicative_generated_metadata

module contraction_trace
  implicit none
  integer :: calls(0:4, 2) = 0
  integer :: effective_calls = 0
end module contraction_trace
'''
        generated = r'''
integer function sdm_multiplicative_block_count()
  sdm_multiplicative_block_count = 2
end function sdm_multiplicative_block_count

integer function sdm_multiplicative_physical_block(position)
  integer, intent(in) :: position
  sdm_multiplicative_physical_block = position - 1
end function sdm_multiplicative_physical_block

integer function sdm_multiplicative_component_position(block)
  integer, intent(in) :: block
  sdm_multiplicative_component_position = block + 1
end function sdm_multiplicative_component_position

module traced_basis_stub
  implicit none
  integer :: saved_counts(2), saved_kinds(8,2)
end module traced_basis_stub

subroutine sdm_multiplicative_prepare_basis( &
     event_slots, maximum_primitives, primitive_counts, insertion_kinds, &
     insertion_ids, insertion_ranks, correlation_legs, include_virtual, &
     precision_asked, precision_found, return_code)
  use contraction_trace
  use traced_basis_stub
  integer, intent(in) :: event_slots(*), maximum_primitives
  integer, intent(in) :: primitive_counts(*)
  integer, intent(in) :: insertion_kinds(maximum_primitives, *)
  integer, intent(in) :: insertion_ids(maximum_primitives, *)
  integer, intent(in) :: insertion_ranks(maximum_primitives, *)
  integer, intent(in) :: correlation_legs(maximum_primitives, *)
  logical, intent(in) :: include_virtual
  double precision, intent(in) :: precision_asked
  double precision, intent(out) :: precision_found
  integer, intent(out) :: return_code

  integer :: position, primitive

  saved_counts = primitive_counts(1:2)
  saved_kinds = 0
  do position = 1, 2
    do primitive = 1, primitive_counts(position)
      if (.not. include_virtual .and. &
          insertion_kinds(primitive, position) == 3) cycle
      calls(insertion_kinds(primitive, position), position) = &
           calls(insertion_kinds(primitive, position), position) + 1
      saved_kinds(primitive, position) = &
           insertion_kinds(primitive, position)
    end do
  end do
  precision_found = precision_asked
  return_code = 0
end subroutine sdm_multiplicative_prepare_basis

subroutine sdm_multiplicative_evaluate_basis( &
     maximum_primitives, primitive_counts, coefficients, result)
  use contraction_trace
  use traced_basis_stub
  integer, intent(in) :: maximum_primitives, primitive_counts(*)
  complex(kind=8), intent(in) :: coefficients(maximum_primitives, *)
  complex(kind=8), intent(out) :: result
  integer :: position, primitive
  complex(kind=8) :: effective(2)

  effective_calls = effective_calls + 1
  effective = (0d0, 0d0)
  do position = 1, 2
    do primitive = 1, primitive_counts(position)
      effective(position) = effective(position) + &
           coefficients(primitive, position)* &
           dble(saved_kinds(primitive, position))* &
           merge(10d0, 1d0, position == 1)
    end do
  end do
  result = effective(1)*effective(2)
end subroutine sdm_multiplicative_evaluate_basis
'''
        program = r'''
program test_loop_times_loop
  use contraction_trace
  use multiplicative_density_terms
  use multiplicative_density_contraction
  implicit none
  type(block_nlo_distribution) :: distributions(2)
  type(multiplicative_density_tuple) :: tuple
  type(density_tuple_schedule) :: schedule
  type(multiplicative_density_basis) :: basis
  complex(kind=8) :: coefficients(3), result
  double precision :: logs(0:2)
  integer :: distribution, primitive

  coefficients = (0d0, 0d0)
  coefficients(1) = (1d0, 0d0)
  logs = 0d0
  do distribution = 1, 2
    call initialize_block_distribution( &
         distributions(distribution), distribution - 1, 1)
    call initialize_block_distribution_term( &
         distributions(distribution), 1, 0, 1, 1, 2)
    do primitive = 1, 2
      ! Kinds two and three stand for real and virtual insertions.  Both
      ! blocks therefore expose a virtual primitive on the same Born point.
      call set_density_primitive( &
           distributions(distribution)%terms(1), primitive, &
           primitive + 1, primitive, 1, 0, 1, coefficients, .true.)
    end do
    call finalize_block_distribution_term(distributions(distribution)%terms(1))
    call finalize_block_distribution(distributions(distribution))
  end do

  call initialize_density_tuple_schedule(distributions, schedule)
  call decode_scheduled_density_tuple(distributions, schedule, 1, tuple)
  call prepare_multiplicative_density_basis( &
       distributions, tuple, 1d-6, basis, .true.)
  call evaluate_multiplicative_density_basis( &
       basis, logs, logs, result)
  if (basis%nlo_order /= 2 .or. basis%return_code /= 0) stop 11
  if (abs(result - (250d0, 0d0)) > 1d-12) stop 12
  if (effective_calls /= 1 .or. sum(calls) /= 4) stop 13
  if (calls(2, 1) /= 1 .or. calls(3, 1) /= 1 .or. &
      calls(2, 2) /= 1 .or. calls(3, 2) /= 1) stop 14

  calls = 0
  effective_calls = 0
  call prepare_multiplicative_density_basis( &
       distributions, tuple, 1d-6, basis, .true., .false.)
  call evaluate_multiplicative_density_basis( &
       basis, logs, logs, result, include_virtual=.false.)
  if (abs(result - (40d0, 0d0)) > 1d-12) stop 15
  if (effective_calls /= 1 .or. sum(calls) /= 2) stop 16
  if (any(calls(3, :) /= 0)) stop 17
end program test_loop_times_loop
'''
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for name, source in (
                    ('stubs.f90', stubs), ('generated.f90', generated),
                    ('test.f90', program)):
                path = os.path.join(directory, name)
                with open(path, 'w') as stream:
                    stream.write(source)
                paths.append(path)
            executable = os.path.join(directory, 'test')
            subprocess.check_call([
                self.compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory, paths[0]] +
                self.sources + paths[1:] + ['-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()

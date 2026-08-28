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
  integer, parameter :: spin_density_color_insertion = 4
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

subroutine sdm_multiplicative_contraction( &
     event_slots, insertion_kinds, insertion_ids, insertion_ranks, &
     correlation_legs, result, precision_asked, precision_found, &
     return_code)
  integer, intent(in) :: event_slots(*), insertion_kinds(*)
  integer, intent(in) :: insertion_ids(*), insertion_ranks(*)
  integer, intent(in) :: correlation_legs(*)
  complex(kind=8), intent(out) :: result
  double precision, intent(in) :: precision_asked
  double precision, intent(out) :: precision_found
  integer, intent(out) :: return_code
  integer :: position

  result = (1d0, 0d0)
  do position = 1, 2
    result = result*dble(1 + event_slots(position) + &
                         insertion_kinds(position))
  end do
  precision_found = precision_asked
  return_code = 0
end subroutine sdm_multiplicative_contraction
'''
        program = r'''
program test_density_contraction
  use multiplicative_density_terms
  use multiplicative_density_contraction
  implicit none
  type(block_nlo_distribution) :: distributions(2)
  complex(kind=8) :: coefficients(3), result
  double precision :: logs(0:4), precision, coupling_rescaling(0:4,0:1)
  integer :: slots(0:4), order, return_code, distribution
  integer :: term_slots(2), term_signs(2), term_kinds(2)
  double precision :: term_coefficients(2)
  double precision :: expected(4)

  term_slots = (/3, 0/)
  term_signs = (/1, -1/)
  term_kinds = (/2, 1/)
  expected = (/360d0, -180d0, -168d0, 84d0/)
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

  do distribution = 1, 4
    call contract_multiplicative_density_tuple( &
         distributions, distribution, logs, logs, 1d-6, result, order, &
         slots, precision, return_code)
    if (abs(dble(result) - expected(distribution)) > 1d-12) stop 11
    if (abs(aimag(result)) > 1d-12) stop 12
    if (order /= 2 .or. return_code /= 0) stop 13
    if (abs(precision - 1d-6) > 1d-15) stop 14
    call contract_multiplicative_density_tuple( &
         distributions, distribution, logs, logs, 1d-6, result, order, &
         slots, precision, return_code, coupling_rescaling)
    if (abs(dble(result) - 6d0*expected(distribution)) > 1d-12) stop 15
  end do

contains

  subroutine fill_term(distribution, term, event_slot, sign, kind, value)
    type(block_nlo_distribution), intent(inout) :: distribution
    integer, intent(in) :: term, event_slot, sign, kind
    double precision, intent(in) :: value

    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx(value, 0d0, kind=8)
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
  integer, parameter :: spin_density_color_insertion = 4
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

module contraction_trace
  implicit none
  integer :: calls(0:4, 0:4) = 0
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

subroutine sdm_multiplicative_contraction( &
     event_slots, insertion_kinds, insertion_ids, insertion_ranks, &
     correlation_legs, result, precision_asked, precision_found, &
     return_code)
  use contraction_trace
  integer, intent(in) :: event_slots(*), insertion_kinds(*)
  integer, intent(in) :: insertion_ids(*), insertion_ranks(*)
  integer, intent(in) :: correlation_legs(*)
  complex(kind=8), intent(out) :: result
  double precision, intent(in) :: precision_asked
  double precision, intent(out) :: precision_found
  integer, intent(out) :: return_code

  calls(insertion_kinds(1), insertion_kinds(2)) = &
       calls(insertion_kinds(1), insertion_kinds(2)) + 1
  result = cmplx(10*insertion_kinds(1) + insertion_kinds(2), &
                 0, kind=8)
  precision_found = precision_asked
  return_code = 0
end subroutine sdm_multiplicative_contraction
'''
        program = r'''
program test_loop_times_loop
  use contraction_trace
  use multiplicative_density_terms
  use multiplicative_density_contraction
  implicit none
  type(block_nlo_distribution) :: distributions(2)
  complex(kind=8) :: coefficients(3), result
  double precision :: logs(0:2), precision
  integer :: slots(0:2), order, return_code, distribution, primitive

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

  call contract_multiplicative_density_tuple( &
       distributions, 1, logs, logs, 1d-6, result, order, slots, &
       precision, return_code)
  if (order /= 2 .or. return_code /= 0) stop 11
  if (abs(result - (110d0, 0d0)) > 1d-12) stop 12
  if (sum(calls) /= 4) stop 13
  if (calls(2, 2) /= 1 .or. calls(2, 3) /= 1 .or. &
      calls(3, 2) /= 1 .or. calls(3, 3) /= 1) stop 14
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

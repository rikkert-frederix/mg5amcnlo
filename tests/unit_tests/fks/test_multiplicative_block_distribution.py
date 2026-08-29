from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeBlockDistributionTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_distribution_retains_all_signed_momentum_terms(self):
        compiler = shutil.which('gfortran')
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        density_terms = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_density_terms.f90')
        block_builder = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_block_distribution.f90')
        base_stubs = r'''
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

module multiplicative_scale_state
  implicit none
contains
  subroutine activate_multiplicative_block_reference(block)
    integer, intent(in) :: block
  end subroutine activate_multiplicative_block_reference
end module multiplicative_scale_state

module nlo_contribution_bundle
  implicit none
contains
  logical function contribution_is_nlo_decay(contribution)
    integer, intent(in) :: contribution
    contribution_is_nlo_decay = .true.
  end function contribution_is_nlo_decay
  integer function contribution_corrected_node(contribution)
    integer, intent(in) :: contribution
    contribution_corrected_node = 2
  end function contribution_corrected_node
end module nlo_contribution_bundle
'''
        builders = r'''
module multiplicative_nbody_density
  use multiplicative_density_terms
  implicit none
contains
  subroutine build_multiplicative_lo_only_density_term(block, term)
    integer, intent(in) :: block
    type(block_distribution_term), intent(out) :: term
    complex(kind=8) :: coefficients(3)
    type(block_nlo_distribution) :: holder
    coefficients = (0d0, 0d0)
    coefficients(1) = (1d0, 0d0)
    call initialize_block_distribution(holder, block, 1)
    call initialize_block_distribution_term(holder, 1, 0, 1, 0, 1)
    call set_density_primitive( &
         holder%terms(1), 1, 0, 0, 0, 0, 0, coefficients, .true.)
    call finalize_block_distribution_term(holder%terms(1))
    term = holder%terms(1)
  end subroutine build_multiplicative_lo_only_density_term

  subroutine build_multiplicative_lo_density_term( &
       contribution, term, available, channel_partition)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    double precision, intent(in), optional :: channel_partition
    call make_term(term, 0, 1, 0, 1d0)
    available = .true.
  end subroutine build_multiplicative_lo_density_term

  subroutine build_multiplicative_nbody_density_term( &
       contribution, term, available, channel_partition)
    integer, intent(in) :: contribution
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    double precision, intent(in), optional :: channel_partition
    call make_term(term, 0, 1, 1, 2d0)
    available = .true.
  end subroutine build_multiplicative_nbody_density_term

  subroutine make_term(term, slot, sign, order, value)
    type(block_distribution_term), intent(out) :: term
    integer, intent(in) :: slot, sign, order
    double precision, intent(in) :: value
    complex(kind=8) :: coefficients(3)
    type(block_nlo_distribution) :: holder
    coefficients = (0d0, 0d0)
    coefficients(1) = cmplx(value, 0d0, kind=8)
    call initialize_block_distribution(holder, 2, 1)
    call initialize_block_distribution_term( &
         holder, 1, slot, sign, order, 1)
    call set_density_primitive( &
         holder%terms(1), 1, merge(0, 1, order == 0), &
         merge(0, 2, order == 0), merge(0, 1, order == 0), 0, order, &
         coefficients, .true.)
    call finalize_block_distribution_term(holder%terms(1))
    term = holder%terms(1)
  end subroutine make_term
end module multiplicative_nbody_density

module multiplicative_fks_density
  use multiplicative_density_terms
  implicit none
contains
  subroutine build_multiplicative_real_density_term(c, term, available)
    integer, intent(in) :: c
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    call make_term(term, 3, 1)
    available = .true.
  end subroutine build_multiplicative_real_density_term
  subroutine build_multiplicative_soft_density_term(c, term, available)
    integer, intent(in) :: c
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    call make_term(term, 0, -1)
    available = .true.
  end subroutine build_multiplicative_soft_density_term
  subroutine build_multiplicative_collinear_density_term(c, term, available)
    integer, intent(in) :: c
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    call make_term(term, 1, -1)
    available = .true.
  end subroutine build_multiplicative_collinear_density_term
  subroutine build_multiplicative_soft_collinear_density_term( &
       c, term, available)
    integer, intent(in) :: c
    type(block_distribution_term), intent(out) :: term
    logical, intent(out) :: available
    call make_term(term, 2, 1)
    available = .true.
  end subroutine build_multiplicative_soft_collinear_density_term

  subroutine make_term(term, slot, sign)
    type(block_distribution_term), intent(out) :: term
    integer, intent(in) :: slot, sign
    complex(kind=8) :: coefficients(3)
    type(block_nlo_distribution) :: holder
    coefficients = (0d0, 0d0)
    coefficients(1) = (1d0, 0d0)
    call initialize_block_distribution(holder, 2, 1)
    call initialize_block_distribution_term(holder, 1, slot, sign, 1, 1)
    call set_density_primitive( &
         holder%terms(1), 1, 1, 2, 1, 0, 1, coefficients, .true.)
    call finalize_block_distribution_term(holder%terms(1))
    term = holder%terms(1)
  end subroutine make_term
end module multiplicative_fks_density
'''
        program = r'''
program test_block_distribution
  use multiplicative_density_terms
  use multiplicative_block_distribution
  implicit none
  type(block_nlo_distribution) :: distribution
  type(block_nlo_distribution) :: lo_only
  integer :: expected_slots(6), expected_signs(6), expected_orders(6), term
  logical :: available

  expected_slots = (/0, 0, 3, 0, 1, 2/)
  expected_signs = (/1, 1, 1, -1, -1, 1/)
  expected_orders = (/0, 1, 1, 1, 1, 1/)
  call build_multiplicative_block_nlo_distribution( &
       2, distribution, available)
  if (.not. available .or. .not. distribution%finalized) stop 11
  if (distribution%block /= 2 .or. distribution%term_count /= 6) stop 12
  do term = 1, 6
    if (distribution%terms(term)%event_slot /= expected_slots(term)) stop 13
    if (distribution%terms(term)%sign /= expected_signs(term)) stop 14
    if (distribution%terms(term)%nlo_order /= expected_orders(term)) stop 15
  end do
  if (abs(dble(distribution%terms(1)%primitives(1)% &
          scale_coefficients(1)) - 1d0) > 1d-12) stop 16
  if (abs(dble(distribution%terms(2)%primitives(1)% &
          scale_coefficients(1)) - 2d0) > 1d-12) stop 17
  call build_multiplicative_lo_only_distribution(4, lo_only)
  if (.not. lo_only%finalized .or. lo_only%block /= 4) stop 18
  if (lo_only%term_count /= 1 .or. &
      lo_only%terms(1)%event_slot /= 0 .or. &
      lo_only%terms(1)%nlo_order /= 0) stop 19
end program test_block_distribution
'''
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for name, source in (
                    ('base.f90', base_stubs), ('builders.f90', builders),
                    ('test.f90', program)):
                path = os.path.join(directory, name)
                with open(path, 'w') as stream:
                    stream.write(source)
                paths.append(path)
            executable = os.path.join(directory, 'test')
            subprocess.check_call([
                compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory, paths[0], density_terms,
                paths[1], block_builder, paths[2], '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()

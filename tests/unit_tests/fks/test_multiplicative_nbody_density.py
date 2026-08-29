from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeNBodyDensityTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.compiler = shutil.which('gfortran')
        cls.repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        cls.sources = [os.path.join(
            cls.repository, 'Template', 'fNLO', 'SubProcesses', name)
            for name in (
                'multiplicative_density_terms.f90',
                'density_operator_recorder.f90',
                'multiplicative_nbody_density.f90')]

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_lo_and_finite_nlo_share_born_slot_but_remain_distinct_terms(self):
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
  integer :: nfksprocess = 5
  double precision :: fkssymmetryfactorborn = 11d0
  double precision :: xibsvcut_used = 0.5d0
end module fnlo_process_common

module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_no_insertion = 0
  integer, parameter :: spin_density_born_insertion = 1
  integer, parameter :: spin_density_virtual_insertion = 3
  integer, parameter :: spin_density_color_insertion = 4
  integer, parameter :: spin_density_fast_virtual_insertion = 5
end module spin_density_matrix_results

module factorized_phase_space
  implicit none
  type :: factorized_radiation_state
    double precision :: jacobian = 13d0
    double precision :: xi_hat = 0.1d0
    double precision :: xi_max = 0.6d0
    double precision :: xi_norm = 0.4d0
  end type factorized_radiation_state
  type :: factorized_measure_state
    double precision :: jacobian = 17d0
    double precision :: phase_space_weight = 19d0
  end type factorized_measure_state
contains
  subroutine fetch_factorized_radiation_state(slot, block, state, available)
    integer, intent(in) :: slot, block
    type(factorized_radiation_state), intent(out) :: state
    logical, intent(out) :: available
    state = factorized_radiation_state()
    if (slot == 0) state%xi_max = 0.4d0
    available = (slot == 0 .or. slot == 3) .and. block == 2
  end subroutine fetch_factorized_radiation_state

  subroutine fetch_factorized_event_measure(slot, block, state, available)
    integer, intent(in) :: slot, block
    type(factorized_measure_state), intent(out) :: state
    logical, intent(out) :: available
    state = factorized_measure_state()
    available = slot == 0 .and. block == 2
  end subroutine fetch_factorized_event_measure
end module factorized_phase_space

module nlo_contribution_bundle
  implicit none
contains
  integer function contribution_for_fks(configuration)
    integer, intent(in) :: configuration
    contribution_for_fks = 2
  end function contribution_for_fks
  logical function contribution_has_virtual(contribution)
    integer, intent(in) :: contribution
    contribution_has_virtual = contribution == 2
  end function contribution_has_virtual
  logical function contribution_has_fast_virtual(contribution)
    integer, intent(in) :: contribution
    contribution_has_fast_virtual = contribution == 2
  end function contribution_has_fast_virtual
  logical function contribution_is_nlo_decay(contribution)
    integer, intent(in) :: contribution
    contribution_is_nlo_decay = contribution == 2
  end function contribution_is_nlo_decay
  integer function contribution_corrected_node(contribution)
    integer, intent(in) :: contribution
    contribution_corrected_node = 2
  end function contribution_corrected_node
  integer function contribution_representative_fks(contribution)
    integer, intent(in) :: contribution
    contribution_representative_fks = 7
  end function contribution_representative_fks
end module nlo_contribution_bundle
'''
        fks_stub = r'''
module fks_singular_module
  use density_operator_recorder, only: record_density_operator_primitive
  implicit none
contains
  double precision function fks_subtraction_shat(slot)
    integer, intent(in) :: slot
    fks_subtraction_shat = 100d0
  end function fks_subtraction_shat

  subroutine record_nbody_integrated_density_operator(slot, born_qcd_power)
    integer, intent(in) :: slot, born_qcd_power
    complex(kind=8) :: coefficients(3)
    if (slot /= 0 .or. born_qcd_power /= 4) stop 90
    coefficients = (/ (2d0, 0d0), (3d0, 0d0), (5d0, 0d0) /)
    call record_density_operator_primitive(1, 2, 1, 0, &
                                           coefficients, .true.)
    coefficients = (0d0, 0d0)
    coefficients(1) = (7d0, 0d0)
    call record_density_operator_primitive(4, 9, 1, 0, &
                                           coefficients, .true.)
  end subroutine record_nbody_integrated_density_operator
end module fks_singular_module
'''
        generated = r'''
integer function sdm_multiplicative_born_qcd_power(block)
  integer, intent(in) :: block
  if (block /= 2) stop 91
  sdm_multiplicative_born_qcd_power = 4
end function sdm_multiplicative_born_qcd_power
'''
        program = r'''
program test_nbody_density
  use multiplicative_density_terms
  use multiplicative_nbody_density
  implicit none
  type(block_distribution_term) :: lo_term, nlo_term
  logical :: available
  double precision :: pi, prefactor

  pi = acos(-1d0)
  ! The shared n-body atom is independent of the sampled FKS sector.  The
  ! caller cancels that sector's sampling probability, so the additive
  ! fkssymmetryfactorborn=11 stub must not enter this prefactor.
  prefactor = 0.4d0/(0.5d0*100d0/(16d0*pi**2))

  call build_multiplicative_lo_density_term( &
       2, lo_term, available, 0.5d0)
  if (.not. available .or. .not. lo_term%finalized) stop 11
  if (lo_term%block /= 2 .or. lo_term%event_slot /= 0) stop 12
  if (lo_term%nlo_order /= 0 .or. lo_term%sign /= 1) stop 13
  if (lo_term%luminosity_configuration /= 7) stop 35
  if (lo_term%primitive_count /= 1) stop 14
  if (lo_term%primitives(1)%insertion_kind /= 0) stop 15
  if (abs(dble(lo_term%primitives(1)%scale_coefficients(1)) - &
          0.5d0*prefactor) > 1d-12) stop 16

  call build_multiplicative_nbody_density_term( &
       2, nlo_term, available, 0.5d0)
  if (.not. available .or. .not. nlo_term%finalized) stop 21
  if (nlo_term%block /= 2 .or. nlo_term%event_slot /= 0) stop 22
  if (nlo_term%nlo_order /= 1 .or. nlo_term%sign /= 1) stop 23
  if (nlo_term%luminosity_configuration /= 7) stop 36
  if (nlo_term%primitive_count /= 3) stop 24
  if (nlo_term%primitives(1)%insertion_kind /= 1) stop 25
  if (nlo_term%primitives(2)%insertion_kind /= 4) stop 26
  if (nlo_term%primitives(3)%insertion_kind /= 5) stop 27
  if (nlo_term%primitives(3)%insertion_rank /= 1) stop 28
  if (abs(dble(nlo_term%primitives(1)%scale_coefficients(1)) - &
          1d0*prefactor) > 1d-12) stop 29
  if (abs(dble(nlo_term%primitives(1)%scale_coefficients(2)) - &
          1.5d0*prefactor) > 1d-12) stop 30
  if (abs(dble(nlo_term%primitives(1)%scale_coefficients(3)) - &
          2.5d0*prefactor) > 1d-12) stop 31
  if (abs(dble(nlo_term%primitives(2)%scale_coefficients(1)) - &
          3.5d0*prefactor) > 1d-12) stop 32
  if (abs(dble(nlo_term%primitives(3)%scale_coefficients(1)) - &
          0.5d0*prefactor) > 1d-12) stop 33
  ! Neither the event Jacobian nor its phase-space weight belongs here.
  if (abs(dble(lo_term%primitives(1)%scale_coefficients(1)) - &
          prefactor*17d0*19d0) < 1d-6) stop 34
end program test_nbody_density
'''
        with tempfile.TemporaryDirectory() as directory:
            paths = {}
            for name, source in (
                    ('stubs.f90', stubs), ('fks_stub.f90', fks_stub),
                    ('generated.f90', generated), ('test.f90', program)):
                path = os.path.join(directory, name)
                with open(path, 'w') as stream:
                    stream.write(source)
                paths[name] = path
            executable = os.path.join(directory, 'test')
            subprocess.check_call([
                self.compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory, paths['stubs.f90'],
                self.sources[0], self.sources[1], paths['fks_stub.f90'],
                self.sources[2], paths['generated.f90'],
                paths['test.f90'], '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()

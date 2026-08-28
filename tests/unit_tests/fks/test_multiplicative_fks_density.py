from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeFKSDensityTest(unittest.TestCase):

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
                'multiplicative_fks_density.f90')]

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_real_operator_retains_real_event_slot_and_local_kernel(self):
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
  integer, parameter :: collinear_counterevent = 1
  integer, parameter :: soft_collinear_counterevent = 2
  integer, parameter :: real_event = 3
  integer :: nfksprocess = 5
  integer :: i_fks = 3
  integer :: j_fks = 4
  double precision :: fkssymmetryfactor = 7d0
  double precision :: fkssymmetryfactordeg = 11d0
  double precision :: xiscut_used = 0.5d0
  double precision :: xicut_used = 0.8d0
  double precision :: delta_used = 0.25d0
end module fnlo_process_common

module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_no_insertion = 0
  integer, parameter :: spin_density_real_insertion = 2
  integer, parameter :: spin_density_color_insertion = 4
end module spin_density_matrix_results

module factorized_phase_space
  implicit none
  type :: factorized_radiation_state
    double precision :: jacobian = 13d0
    double precision :: xi = 0.2d0
    double precision :: y = 0.25d0
    double precision :: xi_norm = 0.4d0
    double precision :: xi_hat = 0.2d0
    double precision :: xi_max = 0.6d0
    double precision :: sqrt_shat = 10d0
  end type factorized_radiation_state
  type :: factorized_measure_state
    double precision :: jacobian = 11d0
    double precision :: phase_space_weight = 17d0
  end type factorized_measure_state
contains
  subroutine fetch_factorized_radiation_state(slot, block, state, available)
    integer, intent(in) :: slot, block
    type(factorized_radiation_state), intent(out) :: state
    logical, intent(out) :: available
    state = factorized_radiation_state()
    available = slot == 3 .and. block == 2
  end subroutine fetch_factorized_radiation_state

  subroutine fetch_factorized_event_measure(slot, block, state, available)
    integer, intent(in) :: slot, block
    type(factorized_measure_state), intent(out) :: state
    logical, intent(out) :: available
    state = factorized_measure_state()
    available = slot == 3 .and. block == 2
  end subroutine fetch_factorized_event_measure
end module factorized_phase_space

module fks_singular_module
  implicit none
contains
  double precision function evaluate_fks_sij(slot, i, j, xi, y)
    integer, intent(in) :: slot, i, j
    double precision, intent(in) :: xi, y
    evaluate_fks_sij = 5d0
  end function evaluate_fks_sij
  subroutine record_soft_density_operator(slot, xi, y, sqrt_shat)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi, y, sqrt_shat
  end subroutine record_soft_density_operator
  subroutine sreal(slot, xi, y, weight)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi, y
    double precision, intent(out) :: weight
    weight = 0d0
  end subroutine sreal
  subroutine sreal_deg(slot, xi, value_xi, value_lxi, &
       xi_index, lxi_index, muf_index)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi
    double precision, intent(out) :: value_xi, value_lxi
    integer, intent(out), optional :: xi_index, lxi_index, muf_index
    value_xi = 0d0
    value_lxi = 0d0
    if (present(xi_index)) xi_index = 0
    if (present(lxi_index)) lxi_index = 0
    if (present(muf_index)) muf_index = 0
  end subroutine sreal_deg
  double precision function fks_subtraction_shat(slot)
    integer, intent(in) :: slot
    fks_subtraction_shat = 1d0
  end function fks_subtraction_shat
end module fks_singular_module

module nlo_contribution_bundle
  implicit none
contains
  integer function contribution_for_fks(configuration)
    integer, intent(in) :: configuration
    contribution_for_fks = 2
  end function contribution_for_fks
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
        generated = r'''
integer function sdm_real_insertion_identifier(configuration)
  integer, intent(in) :: configuration
  sdm_real_insertion_identifier = 23
end function sdm_real_insertion_identifier
'''
        program = r'''
program test_real_density
  use multiplicative_density_terms
  use multiplicative_fks_density
  implicit none
  type(block_distribution_term) :: term
  logical :: available
  double precision :: expected

  call build_multiplicative_real_density_term(2, term, available)
  if (.not. available) stop 11
  if (.not. term%finalized) stop 12
  if (term%block /= 2 .or. term%event_slot /= 3) stop 13
  if (term%sign /= 1 .or. term%nlo_order /= 1) stop 14
  if (term%primitive_count /= 1) stop 15
  if (term%primitives(1)%insertion_kind /= 2) stop 16
  if (term%primitives(1)%insertion_identifier /= 23) stop 17
  if (term%primitives(1)%insertion_rank /= 1) stop 18
  expected = (0.2d0**2*(1d0 - 0.25d0))*5d0* &
       (0.4d0/0.2d0/(1d0 - 0.25d0)*7d0)
  if (abs(dble(term%primitives(1)%scale_coefficients(1)) - &
          expected) > 1d-12) stop 19
  ! The 11*17 phase-space measure must not be hidden in this coefficient.
  if (abs(expected - 2.8d0) > 1d-12) stop 20
end program test_real_density
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
    def test_soft_operator_is_one_distinct_signed_counterevent_term(self):
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
  integer, parameter :: collinear_counterevent = 1
  integer, parameter :: soft_collinear_counterevent = 2
  integer, parameter :: real_event = 3
  integer :: nfksprocess = 5
  integer :: i_fks = 3
  integer :: j_fks = 4
  double precision :: fkssymmetryfactor = 7d0
  double precision :: fkssymmetryfactordeg = 11d0
  double precision :: xiscut_used = 0.5d0
  double precision :: xicut_used = 0.8d0
  double precision :: delta_used = 0.25d0
end module fnlo_process_common

module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_no_insertion = 0
  integer, parameter :: spin_density_real_insertion = 2
  integer, parameter :: spin_density_color_insertion = 4
end module spin_density_matrix_results

module factorized_phase_space
  implicit none
  type :: factorized_radiation_state
    double precision :: jacobian = 13d0
    double precision :: xi = 0.2d0
    double precision :: y = 0.25d0
    double precision :: xi_norm = 0.4d0
    double precision :: xi_hat = 0.2d0
    double precision :: xi_max = 0.6d0
    double precision :: sqrt_shat = 10d0
  end type factorized_radiation_state
  type :: factorized_measure_state
    double precision :: jacobian = 11d0
    double precision :: phase_space_weight = 17d0
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
        fks_stub = r'''
module fks_singular_module
  use density_operator_recorder, only: record_density_operator_primitive
  implicit none
contains
  double precision function evaluate_fks_sij(slot, i, j, xi, y)
    integer, intent(in) :: slot, i, j
    double precision, intent(in) :: xi, y
    evaluate_fks_sij = 5d0
  end function evaluate_fks_sij

  subroutine record_soft_density_operator(slot, xi, y, sqrt_shat)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi, y, sqrt_shat
    complex(kind=8) :: coefficients(3)
    coefficients = (0d0, 0d0)
    coefficients(1) = (2d0, 0d0)
    call record_density_operator_primitive(4, 9, 1, 0, &
                                           coefficients, .true.)
  end subroutine record_soft_density_operator
  subroutine sreal(slot, xi, y, weight)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi, y
    double precision, intent(out) :: weight
    weight = 0d0
  end subroutine sreal
  subroutine sreal_deg(slot, xi, value_xi, value_lxi, &
       xi_index, lxi_index, muf_index)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi
    double precision, intent(out) :: value_xi, value_lxi
    integer, intent(out), optional :: xi_index, lxi_index, muf_index
    value_xi = 0d0
    value_lxi = 0d0
    if (present(xi_index)) xi_index = 0
    if (present(lxi_index)) lxi_index = 0
    if (present(muf_index)) muf_index = 0
  end subroutine sreal_deg
  double precision function fks_subtraction_shat(slot)
    integer, intent(in) :: slot
    fks_subtraction_shat = 1d0
  end function fks_subtraction_shat
end module fks_singular_module
'''
        generated = r'''
integer function sdm_real_insertion_identifier(configuration)
  integer, intent(in) :: configuration
  sdm_real_insertion_identifier = 23
end function sdm_real_insertion_identifier
'''
        program = r'''
program test_soft_density
  use multiplicative_density_terms
  use multiplicative_fks_density
  implicit none
  type(block_distribution_term) :: term
  logical :: available
  double precision :: expected, prefactor, endpoint

  call build_multiplicative_soft_density_term(2, term, available)
  if (.not. available) stop 11
  if (.not. term%finalized) stop 12
  if (term%block /= 2 .or. term%event_slot /= 0) stop 13
  if (term%sign /= -1 .or. term%nlo_order /= 1) stop 14
  if (term%primitive_count /= 1) stop 15
  if (term%primitives(1)%insertion_kind /= 4) stop 16
  if (term%primitives(1)%insertion_identifier /= 9) stop 17
  prefactor = 0.4d0/0.2d0/(1d0 - 0.25d0)
  endpoint = 0.4d0/0.5d0*log(0.8d0/0.5d0)/(1d0 - 0.25d0)
  expected = 2d0*5d0*(prefactor + endpoint)*7d0
  if (abs(dble(term%primitives(1)%scale_coefficients(1)) - &
          expected) > 1d-12) stop 18
  if (abs(dble(term%primitives(1)%scale_coefficients(1)) - &
          expected*11d0*17d0) < 1d-6) stop 19
end program test_soft_density
'''
        with tempfile.TemporaryDirectory() as directory:
            paths = {}
            for name, source in (
                    ('base_stubs.f90', base_stubs),
                    ('fks_stub.f90', fks_stub),
                    ('generated.f90', generated), ('test.f90', program)):
                path = os.path.join(directory, name)
                with open(path, 'w') as stream:
                    stream.write(source)
                paths[name] = path
            executable = os.path.join(directory, 'test')
            subprocess.check_call([
                self.compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory,
                paths['base_stubs.f90'], self.sources[0], self.sources[1],
                paths['fks_stub.f90'], self.sources[2],
                paths['generated.f90'], paths['test.f90'],
                '-o', executable])
            subprocess.check_call([executable])

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_collinear_family_keeps_c_and_sc_as_distinct_terms(self):
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
  integer, parameter :: collinear_counterevent = 1
  integer, parameter :: soft_collinear_counterevent = 2
  integer, parameter :: real_event = 3
  integer :: nfksprocess = 5
  integer :: i_fks = 3
  integer :: j_fks = 4
  double precision :: fkssymmetryfactor = 7d0
  double precision :: fkssymmetryfactordeg = 11d0
  double precision :: xiscut_used = 0.5d0
  double precision :: xicut_used = 0.8d0
  double precision :: delta_used = 0.25d0
end module fnlo_process_common

module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_no_insertion = 0
  integer, parameter :: spin_density_born_insertion = 1
  integer, parameter :: spin_density_real_insertion = 2
  integer, parameter :: spin_density_color_insertion = 4
end module spin_density_matrix_results

module factorized_phase_space
  implicit none
  type :: factorized_radiation_state
    double precision :: jacobian = 13d0
    double precision :: xi = 0.2d0
    double precision :: y = 0.25d0
    double precision :: xi_norm = 0.4d0
    double precision :: xi_hat = 0.1d0
    double precision :: xi_max = 0.6d0
    double precision :: sqrt_shat = 10d0
  end type factorized_radiation_state
  type :: factorized_measure_state
    double precision :: jacobian = 11d0
    double precision :: phase_space_weight = 17d0
  end type factorized_measure_state
contains
  subroutine fetch_factorized_radiation_state(slot, block, state, available)
    integer, intent(in) :: slot, block
    type(factorized_radiation_state), intent(out) :: state
    logical, intent(out) :: available
    state = factorized_radiation_state()
    if (slot == 1) then
      state%xi = 0.25d0
      state%xi_norm = 0.5d0
      state%xi_max = 0.4d0
      state%y = 1d0
    else if (slot == 2) then
      state%xi = 0d0
      state%xi_max = 0.4d0
      state%y = 1d0
    end if
    available = slot >= 1 .and. slot <= 3 .and. block == 2
  end subroutine fetch_factorized_radiation_state

  subroutine fetch_factorized_event_measure(slot, block, state, available)
    integer, intent(in) :: slot, block
    type(factorized_measure_state), intent(out) :: state
    logical, intent(out) :: available
    state = factorized_measure_state()
    available = (slot == 1 .or. slot == 2) .and. block == 2
  end subroutine fetch_factorized_event_measure
end module factorized_phase_space

module nlo_contribution_bundle
  implicit none
contains
  integer function contribution_for_fks(configuration)
    integer, intent(in) :: configuration
    contribution_for_fks = 2
  end function contribution_for_fks
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
        fks_stub = r'''
module fks_singular_module
  use density_operator_recorder, only: &
       record_density_operator_primitive, recorded_density_operator_count
  implicit none
contains
  double precision function evaluate_fks_sij(slot, i, j, xi, y)
    integer, intent(in) :: slot, i, j
    double precision, intent(in) :: xi, y
    evaluate_fks_sij = merge(5d0, 6d0, slot == 1)
  end function evaluate_fks_sij

  double precision function fks_subtraction_shat(slot)
    integer, intent(in) :: slot
    fks_subtraction_shat = merge(100d0, 200d0, slot == 1)
  end function fks_subtraction_shat

  subroutine record_soft_density_operator(slot, xi, y, sqrt_shat)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi, y, sqrt_shat
  end subroutine record_soft_density_operator

  subroutine sreal_deg(slot, xi, collrem_xi, collrem_lxi, &
       xi_index, lxi_index, muf_index)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi
    double precision, intent(out) :: collrem_xi, collrem_lxi
    integer, intent(out), optional :: xi_index, lxi_index, muf_index
    complex(kind=8) :: coefficients(3)
    collrem_xi = 0d0
    collrem_lxi = 0d0
    coefficients = (0d0, 0d0)
    coefficients(1) = (2d0, 0d0)
    call record_density_operator_primitive(1, 2, 1, 4, coefficients, .true.)
    if (present(xi_index)) xi_index = recorded_density_operator_count()
    coefficients = (0d0, 0d0)
    coefficients(1) = (3d0, 0d0)
    call record_density_operator_primitive(1, 2, 1, 4, coefficients, .true.)
    if (present(lxi_index)) lxi_index = recorded_density_operator_count()
    coefficients = (0d0, 0d0)
    coefficients(3) = (4d0, 0d0)
    call record_density_operator_primitive(1, 2, 1, 4, coefficients, .true.)
    if (present(muf_index)) muf_index = recorded_density_operator_count()
  end subroutine sreal_deg

  subroutine sreal(slot, xi, y, weight)
    integer, intent(in) :: slot
    double precision, intent(in) :: xi, y
    double precision, intent(out) :: weight
    complex(kind=8) :: coefficients(3)
    weight = 0d0
    coefficients = (0d0, 0d0)
    coefficients(1) = (5d0, 0d0)
    call record_density_operator_primitive(1, 2, 1, 4, coefficients, .true.)
    coefficients(1) = (7d0, 0d0)
    call record_density_operator_primitive(1, 2, 2, 4, coefficients, .true.)
  end subroutine sreal
end module fks_singular_module
'''
        generated = r'''
integer function sdm_real_insertion_identifier(configuration)
  integer, intent(in) :: configuration
  sdm_real_insertion_identifier = 23
end function sdm_real_insertion_identifier
'''
        program = r'''
program test_collinear_density
  use multiplicative_density_terms
  use multiplicative_fks_density
  implicit none
  type(block_distribution_term) :: cterm, scterm
  logical :: available
  double precision :: pi, fc, fdc, core_scale
  double precision :: cutoff, fsc, fdsc1, fdsc2, fdsc3, fdsc4
  double precision :: pref_c, pref_coll, pref_cnt, pref_cnt_coll

  pi = acos(-1d0)
  call build_multiplicative_collinear_density_term(2, cterm, available)
  if (.not. available .or. .not. cterm%finalized) stop 11
  if (cterm%event_slot /= 1 .or. cterm%sign /= -1) stop 12
  if (cterm%primitive_count /= 5) stop 13
  fdc = (0.5d0/0.25d0)/(100d0/(32d0*pi**2))*11d0
  fc = (0.5d0/0.25d0/0.75d0 + &
        0.5d0/0.25d0*log(0.25d0))*7d0
  core_scale = 5d0*fc
  call assert_close(cterm%primitives(1)%scale_coefficients(1), -2d0*fdc, 21)
  call assert_close(cterm%primitives(2)%scale_coefficients(1), &
                    -3d0*log(0.25d0)*fdc, 22)
  call assert_close(cterm%primitives(3)%scale_coefficients(3), -4d0*fdc, 23)
  call assert_close(cterm%primitives(4)%scale_coefficients(1), &
                    5d0*core_scale, 24)
  call assert_close(cterm%primitives(5)%scale_coefficients(1), &
                    7d0*core_scale, 25)

  call build_multiplicative_soft_collinear_density_term( &
       2, scterm, available)
  if (.not. available .or. .not. scterm%finalized) stop 31
  if (scterm%event_slot /= 2 .or. scterm%sign /= 1) stop 32
  if (scterm%primitive_count /= 5) stop 33
  cutoff = 0.4d0
  pref_c = 0.5d0/0.25d0/0.75d0
  pref_coll = 0.5d0/0.25d0*log(0.25d0)
  pref_cnt = 0.5d0/cutoff*log(0.8d0/cutoff)/0.75d0
  pref_cnt_coll = 0.5d0/cutoff*log(0.8d0/cutoff)*log(0.25d0)
  fsc = (pref_c + pref_coll + pref_cnt + pref_cnt_coll)*11d0
  fdsc1 = (0.5d0/0.25d0)/(200d0/(32d0*pi**2))*11d0
  fdsc2 = 0.5d0/cutoff*log(0.8d0/cutoff)*11d0
  fdsc3 = 0.5d0/cutoff*(log(0.8d0)**2 - &
                        log(cutoff)**2)/2d0*11d0
  fdsc4 = fdsc2
  call assert_close(scterm%primitives(1)%scale_coefficients(1), &
                    -2d0*(fdsc1 + fdsc2), 41)
  call assert_close(scterm%primitives(2)%scale_coefficients(1), &
                    -3d0*(log(0.25d0)*fdsc1 + fdsc3), 42)
  call assert_close(scterm%primitives(3)%scale_coefficients(3), &
                    -4d0*fdsc4, 43)
  call assert_close(scterm%primitives(4)%scale_coefficients(1), &
                    5d0*6d0*fsc, 44)
  call assert_close(scterm%primitives(5)%scale_coefficients(1), &
                    7d0*6d0*fsc, 45)

contains
  subroutine assert_close(actual, expected, code)
    complex(kind=8), intent(in) :: actual
    double precision, intent(in) :: expected
    integer, intent(in) :: code
    if (abs(actual - cmplx(expected, 0d0, kind=8)) > &
        1d-10*max(1d0, abs(expected))) then
      write (*, *) 'assertion failed', code, actual, expected
      stop 99
    end if
  end subroutine assert_close
end program test_collinear_density
'''
        with tempfile.TemporaryDirectory() as directory:
            paths = {}
            for name, source in (
                    ('base_stubs.f90', base_stubs),
                    ('fks_stub.f90', fks_stub),
                    ('generated.f90', generated), ('test.f90', program)):
                path = os.path.join(directory, name)
                with open(path, 'w') as stream:
                    stream.write(source)
                paths[name] = path
            executable = os.path.join(directory, 'test')
            subprocess.check_call([
                self.compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory,
                paths['base_stubs.f90'], self.sources[0], self.sources[1],
                paths['fks_stub.f90'], self.sources[2],
                paths['generated.f90'], paths['test.f90'],
                '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()

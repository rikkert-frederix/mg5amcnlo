from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class SpinDensityMatrixResultsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.compiler = shutil.which('gfortran')
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        cls.source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'spin_density_matrix_results.f90')

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_insertions_are_cached_by_momenta_and_operator(self):
        stubs = r'''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 2
end module process_dimensions

module fnlo_process_common
  implicit none
  integer, parameter :: soft_counterevent = 0
  integer, parameter :: real_event = 3
end module fnlo_process_common

module factorized_phase_space
  implicit none
  integer(kind=8) :: revisions(0:2,0:3) = 0_8
contains
  integer(kind=8) function factorized_block_momentum_revision(event_slot, block)
    integer, intent(in) :: event_slot, block
    factorized_block_momentum_revision = revisions(block,event_slot)
  end function factorized_block_momentum_revision
end module factorized_phase_space

module multiplicative_scale_state
  implicit none
  integer(kind=8) :: active_context = 1_8
contains
  integer(kind=8) function multiplicative_active_coupling_context()
    multiplicative_active_coupling_context = active_context
  end function multiplicative_active_coupling_context
end module multiplicative_scale_state
'''
        program = r'''
program test_insertion_cache
  use factorized_phase_space, only: revisions
  use multiplicative_scale_state, only: active_context
  use spin_density_matrix_results
  implicit none
  type(spin_density_block_result) :: stored, loaded
  type(spin_density_cache_statistics) :: statistics
  complex(kind=8) :: density(2,2,2)
  logical :: available
  double precision :: precision_found
  integer :: return_code, rank, left, right

  do rank = 1, 2
    do left = 1, 2
      do right = 1, 2
        density(rank,left,right) = cmplx( &
             dble(100*rank + 10*left + right), dble(rank-left), kind=8)
      end do
    end do
  end do
  revisions(1,0) = 17_8
  call initialize_spin_density_block(stored, 0, 1, 2)
  call load_cached_spin_density_insertion( &
       stored, spin_density_real_insertion, 7, 0, 1d-5, available, &
       precision_found, return_code)
  if (available) stop 11
  call record_spin_density_insertion( &
       stored, spin_density_real_insertion, 1, 7, 0, 1d-5, 2d-7, 3, &
       density)
  if (.not. stored%has_insertion .or. &
      any(stored%insertion /= density)) stop 12

  call initialize_spin_density_block(loaded, 0, 1, 2)
  call load_cached_spin_density_insertion( &
       loaded, spin_density_real_insertion, 7, 0, 1d-5, available, &
       precision_found, return_code)
  if (.not. available .or. any(loaded%insertion /= density)) stop 13
  if (loaded%insertion_order /= 1 .or. &
      loaded%insertion_kind /= spin_density_real_insertion) stop 14
  if (precision_found /= 2d-7 .or. return_code /= 3) stop 15

  active_context = 2_8
  call initialize_spin_density_block(loaded, 0, 1, 2)
  call load_cached_spin_density_insertion( &
       loaded, spin_density_real_insertion, 7, 0, 1d-5, available, &
       precision_found, return_code)
  if (available) stop 21
  active_context = 1_8

  call initialize_spin_density_block(loaded, 0, 1, 2)
  call load_cached_spin_density_insertion( &
       loaded, spin_density_real_insertion, 8, 0, 1d-5, available, &
       precision_found, return_code)
  if (available) stop 16
  revisions(1,0) = 18_8
  call initialize_spin_density_block(loaded, 0, 1, 2)
  call load_cached_spin_density_insertion( &
       loaded, spin_density_real_insertion, 7, 0, 1d-5, available, &
       precision_found, return_code)
  if (available) stop 17

  revisions(1,0) = 19_8
  call initialize_spin_density_block(stored, 0, 1, 2)
  call record_spin_density_insertion( &
       stored, spin_density_real_insertion, 1, 7, 0, 1d-5, 2d-7, 3, &
       density)
  call reset_spin_density_caches()
  call initialize_spin_density_block(loaded, 0, 1, 2)
  call load_cached_spin_density_insertion( &
       loaded, spin_density_real_insertion, 7, 0, 1d-5, available, &
       precision_found, return_code)
  if (available) stop 18
  call fetch_spin_density_cache_statistics(statistics)
  if (statistics%insertion_hits /= 1_8 .or. &
      statistics%insertion_misses /= 5_8) stop 19
  if (statistics%insertion_provider_evaluations /= 2_8) stop 20

  revisions(1,0) = 20_8
  call initialize_spin_density_block(stored, 0, 1, 2)
  call record_spin_density_insertion( &
       stored, spin_density_fast_virtual_insertion, 1, 9, 0, 1d-5, &
       0d0, 0, density)
  if (stored%insertion_kind /= spin_density_fast_virtual_insertion) stop 22
  call record_raw_amplitude_cache_miss()
  call record_raw_amplitude_cache_hit()
  call fetch_spin_density_cache_statistics(statistics)
  if (statistics%raw_amplitude_hits /= 1_8 .or. &
      statistics%raw_amplitude_misses /= 1_8 .or. &
      statistics%raw_amplitude_provider_evaluations /= 1_8) stop 23
  call record_direct_virtual_reconstruction(.true., 2)
  call record_direct_virtual_reconstruction(.false., 8)
  call record_virtual_tomography_reconstruction(7)
  call fetch_spin_density_cache_statistics(statistics)
  if (statistics%virtual_madloop_evaluations /= 17_8 .or. &
      statistics%virtual_direct_reconstructions /= 1_8 .or. &
      statistics%virtual_direct_fallbacks /= 1_8 .or. &
      statistics%virtual_tomography_reconstructions /= 1_8) stop 24
  call record_exact_family_candidate()
  call record_exact_family_cut_rejection()
  call record_exact_family_candidate()
  call record_exact_family_acceptance()
  call record_density_contraction()
  call record_scale_reweight_evaluation()
  call record_pdf_member_initialization()
  call record_pdf_luminosity_evaluation()
  call record_histogram_family_fill()
  call fetch_spin_density_cache_statistics(statistics)
  if (statistics%exact_family_candidates /= 2_8 .or. &
      statistics%exact_family_cut_rejections /= 1_8 .or. &
      statistics%exact_family_acceptances /= 1_8 .or. &
      statistics%density_contractions /= 1_8 .or. &
      statistics%scale_reweight_evaluations /= 1_8 .or. &
      statistics%pdf_member_initializations /= 1_8 .or. &
      statistics%pdf_luminosity_evaluations /= 1_8 .or. &
      statistics%histogram_family_fills /= 1_8) stop 25
end program test_insertion_cache
'''
        with tempfile.TemporaryDirectory() as directory:
            stubs_path = os.path.join(directory, 'stubs.f90')
            program_path = os.path.join(directory, 'test.f90')
            with open(stubs_path, 'w') as stream:
                stream.write(stubs)
            with open(program_path, 'w') as stream:
                stream.write(program)
            executable = os.path.join(directory, 'test')
            subprocess.check_call([
                self.compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory, stubs_path, self.source,
                program_path, '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()

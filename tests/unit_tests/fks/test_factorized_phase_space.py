from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class FactorizedPhaseSpaceTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.compiler = shutil.which('gfortran')
        cls.repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        cls.source = os.path.join(
            cls.repository, 'Template', 'fNLO', 'SubProcesses',
            'factorized_phase_space.f90')

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_matrix_momenta_follow_each_provider_leg_order(self):
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
'''
        program = r'''
program test_ordered_matrix_momenta
  use factorized_phase_space
  use fnlo_process_common, only: soft_counterevent, real_event
  implicit none
  double precision :: stored(0:3,4), born(0:3,3), real(0:3,4)
  integer :: pdgs(4), target_kinds(4), target_ids(4)
  integer :: born_pdgs(3), born_final(3), real_final(4)
  integer :: production_pdgs(4)
  logical :: is_final(4), available
  integer :: particle, component
  integer(kind=8) :: revision, repeated_revision

  do particle = 1, 4
    do component = 0, 3
      stored(component, particle) = 100d0*particle + component
    end do
  end do
  ! The stored real-decay layout is parent,b,W,g, whereas the reusable
  ! underlying-Born provider generated from ``t > W b`` expects parent,W,b.
  pdgs = (/ 6, 5, 24, 21 /)
  is_final = (/ .false., .true., .true., .true. /)
  target_kinds = (/ factorized_no_target, factorized_visible_target, &
                    factorized_visible_target, factorized_visible_target /)
  target_ids = (/ 0, 1, 2, 3 /)
  call store_factorized_local_momenta(soft_counterevent, 1, 4, stored)
  call store_factorized_local_layout(soft_counterevent, 1, 4, pdgs, is_final, &
                                     target_kinds, target_ids)
  call store_factorized_block_momenta(soft_counterevent, 1, 3, stored(:,1:3))

  born_pdgs = (/ 6, 24, 5 /)
  born_final = (/ 0, 1, 1 /)
  call fetch_factorized_ordered_matrix_momenta( &
       soft_counterevent, 1, 3, born_pdgs, born_final, born, available)
  if (.not. available) stop 11
  if (any(born(:,1) /= stored(:,1))) stop 12
  if (any(born(:,2) /= stored(:,3))) stop 13
  if (any(born(:,3) /= stored(:,2))) stop 14

  ! A non-radiating spectator is only re-boosted into another event slot.
  ! Its momenta must come from that slot while its immutable layout may be
  ! inherited from the underlying-Born slot.
  call store_factorized_block_momenta(1, 1, 3, 2d0*stored(:,1:3))
  call fetch_factorized_ordered_matrix_momenta( &
       1, 1, 3, born_pdgs, born_final, born, available)
  if (.not. available) stop 15
  if (any(born(:,1) /= 2d0*stored(:,1))) stop 16
  if (any(born(:,2) /= 2d0*stored(:,3))) stop 17
  if (any(born(:,3) /= 2d0*stored(:,2))) stop 18

  real_final = (/ 0, 1, 1, 1 /)
  call store_factorized_local_momenta(real_event, 1, 4, stored)
  call store_factorized_local_layout(real_event, 1, 4, pdgs, is_final, &
                                     target_kinds, target_ids)
  call store_factorized_block_momenta(real_event, 1, 4, stored)
  call fetch_factorized_ordered_matrix_momenta( &
       real_event, 1, 4, pdgs, real_final, real, available)
  if (.not. available .or. any(real /= stored)) stop 21

  ! The live real context of an ISR soft projection may start q,g,... even
  ! though the cached reduced matrix block is g,g,... .  Updating that live
  ! layout must not mutate the matrix block's captured Born identity.
  production_pdgs = (/ 21, 21, 6, -6 /)
  is_final = (/ .false., .false., .true., .true. /)
  target_kinds = (/ factorized_no_target, factorized_no_target, &
                    factorized_visible_target, factorized_visible_target /)
  target_ids = (/ 0, 0, 1, 2 /)
  call store_factorized_local_momenta(soft_counterevent, 0, 4, stored)
  call store_factorized_local_layout(soft_counterevent, 0, 4, &
       production_pdgs, is_final, target_kinds, target_ids)
  call store_factorized_block_momenta(soft_counterevent, 0, 4, stored)
  pdgs = (/ -1, 21, 6, -6 /)
  call store_factorized_local_layout(soft_counterevent, 0, 4, pdgs, &
       is_final, target_kinds, target_ids)
  real_final = (/ 0, 0, 1, 1 /)
  call fetch_factorized_ordered_matrix_momenta( &
       soft_counterevent, 0, 4, production_pdgs, real_final, real, available)
  if (.not. available .or. any(real /= stored)) stop 22

  ! A reduced ISR collinear atom keeps its real-context qg identity in the
  ! local event layout, but its density provider needs the underlying gg
  ! identity.  Only the identity falls back; momenta stay in event slot 1.
  call store_factorized_local_momenta(1, 0, 4, 2d0*stored)
  call store_factorized_local_layout(1, 0, 4, pdgs, is_final, &
                                     target_kinds, target_ids)
  call store_factorized_block_momenta(1, 0, 4, 2d0*stored)
  call fetch_factorized_ordered_matrix_momenta( &
       1, 0, 4, production_pdgs, real_final, real, available)
  if (.not. available .or. any(real /= 2d0*stored)) stop 25

  ! Re-realizing an identical tuple must preserve the matrix-input identity,
  ! while an actual four-vector change must invalidate density caches.
  revision = factorized_block_momentum_revision( &
       soft_counterevent, 0)
  call store_factorized_block_momenta( &
       soft_counterevent, 0, 4, stored)
  repeated_revision = factorized_block_momentum_revision( &
       soft_counterevent, 0)
  if (repeated_revision /= revision) stop 23
  stored(0,4) = stored(0,4) + 1d0
  call store_factorized_block_momenta( &
       soft_counterevent, 0, 4, stored)
  if (factorized_block_momentum_revision(soft_counterevent, 0) == &
      repeated_revision) stop 24

  born_pdgs(2) = -24
  call fetch_factorized_ordered_matrix_momenta( &
       soft_counterevent, 1, 3, born_pdgs, born_final, born, available)
  if (available .or. any(born /= 0d0)) stop 31
end program test_ordered_matrix_momenta
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

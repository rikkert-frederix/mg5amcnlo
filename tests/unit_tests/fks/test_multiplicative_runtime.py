from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeRuntimeTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_event_is_materialized_only_after_density_contraction(self):
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_runtime.f90')
        stubs = r'''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 4
end module process_dimensions

module runtime_order_state
  implicit none
  integer :: stage = 0
end module runtime_order_state

module multiplicative_density_terms
  implicit none
  type :: block_nlo_distribution
    integer :: block = 0
  end type block_nlo_distribution
end module multiplicative_density_terms

module fnlo_process_common
  implicit none
  integer, parameter :: soft_counterevent = 0
end module fnlo_process_common

module nlo_contribution_bundle
  implicit none
contains
  integer function multiplicative_event_capacity()
    multiplicative_event_capacity = 6
  end function multiplicative_event_capacity
end module nlo_contribution_bundle

module factorized_phase_space
  use runtime_order_state
  implicit none
  type :: factorized_radiation_state
    double precision :: bjorken_x(2) = -1d0
    double precision :: y_to_lab = 0d0
  end type factorized_radiation_state
contains
  subroutine compose_factorized_tuple_measure(slots, jac, weight, available)
    integer, intent(in) :: slots(0:)
    double precision, intent(out) :: jac, weight
    logical, intent(out) :: available
    if (stage /= 2) stop 81
    stage = 3
    jac = 2d0 + 0d0*slots(0)
    weight = 3d0
    available = .true.
  end subroutine compose_factorized_tuple_measure
  subroutine fetch_factorized_radiation_state(slot, block, value, available)
    integer, intent(in) :: slot, block
    type(factorized_radiation_state), intent(out) :: value
    logical, intent(out) :: available
    if (stage /= 3) stop 82
    value%bjorken_x = (/0.2d0, 0.3d0/)
    value%y_to_lab = 0.4d0
    available = slot == 3 .and. block == 0
  end subroutine fetch_factorized_radiation_state
end module factorized_phase_space

module multiplicative_density_contraction
  use runtime_order_state
  use multiplicative_density_terms
  implicit none
contains
  subroutine contract_multiplicative_density_tuple( &
       distributions, tuple_index, logs_r, logs_f, precision_asked, &
       result, nlo_order, slots, precision_found, return_code, &
       coupling_rescaling, already_realized)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    integer, intent(in) :: tuple_index
    double precision, intent(in) :: logs_r(0:), logs_f(0:)
    double precision, intent(in) :: precision_asked
    complex(kind=8), intent(out) :: result
    integer, intent(out) :: nlo_order, slots(0:), return_code
    double precision, intent(out) :: precision_found
    double precision, intent(in) :: coupling_rescaling(0:,0:)
    logical, intent(in) :: already_realized
    if (stage /= 0) stop 83
    stage = 1
    slots = 0
    slots(0) = 3
    result = (5d0, 0d0) + 0d0*tuple_index + &
         0d0*sum(logs_r) + 0d0*sum(logs_f) + &
         0d0*sum(coupling_rescaling) + 0d0*size(distributions)
    if (already_realized) stop 85
    nlo_order = 2
    precision_found = precision_asked
    return_code = 0
  end subroutine contract_multiplicative_density_tuple
end module multiplicative_density_contraction

module multiplicative_kinematics
  use runtime_order_state
  implicit none
contains
  subroutine materialize_factorized_event_tuple( &
       slots, capacity, count, momenta, pdgs, pass, origins)
    integer, intent(in) :: slots(0:), capacity
    integer, intent(out) :: count, pdgs(capacity)
    double precision, intent(out) :: momenta(0:3, capacity)
    logical, intent(out) :: pass
    integer, intent(out) :: origins(capacity)
    if (stage /= 1) stop 84
    stage = 2
    count = capacity
    momenta = 0d0
    momenta(0, :) = 1d0
    pdgs = 21
    origins = 0
    pass = slots(0) == 3
  end subroutine materialize_factorized_event_tuple
end module multiplicative_kinematics
'''
        program = r'''
program test_runtime
  use runtime_order_state
  use multiplicative_density_terms
  use multiplicative_runtime
  implicit none
  type(block_nlo_distribution) :: distributions(1)
  type(multiplicative_event_evaluation) :: evaluation
  double precision :: logs(0:4), rescaling(0:4,0:1)

  logs = 0d0
  rescaling = 1d0
  call evaluate_multiplicative_event_tuple( &
       distributions, 1, logs, logs, rescaling, 7d0, 1d-6, evaluation)
  if (.not. evaluation%available .or. stage /= 3) stop 1
  if (abs(dble(evaluation%partonic_weight) - 210d0) > 1d-12) stop 2
  if (evaluation%visible_count /= 6 .or. evaluation%nlo_order /= 2) stop 3
  if (maxval(abs(evaluation%bjorken_x - (/0.2d0, 0.3d0/))) > &
      1d-12) stop 4
end program test_runtime
'''
        with tempfile.TemporaryDirectory() as directory:
            stub_path = os.path.join(directory, 'stubs.f90')
            program_path = os.path.join(directory, 'test.f90')
            executable = os.path.join(directory, 'test')
            with open(stub_path, 'w') as stream:
                stream.write(stubs)
            with open(program_path, 'w') as stream:
                stream.write(program)
            subprocess.check_call([
                shutil.which('gfortran'), '-std=f2008',
                '-ffree-line-length-none', '-J', directory, '-I', directory,
                stub_path, source, program_path, '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()

from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeScaleStateTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_each_block_owns_an_activated_reference(self):
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_scale_state.f90')
        stubs = r'''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 5
contains
  subroutine validate_process_dimensions()
  end subroutine validate_process_dimensions
end module process_dimensions

module alfas_functions_module
  implicit none
contains
  double precision function alphas(value)
    double precision, intent(in) :: value
    alphas = 0.1d0 + value/10000d0
  end function alphas
end module alfas_functions_module

module factorized_phase_space
  implicit none
  type :: factorized_radiation_state
    double precision :: sqrt_shat = 0d0
  end type factorized_radiation_state
contains
  subroutine fetch_factorized_radiation_state(slot, block, value, available)
    integer, intent(in) :: slot, block
    type(factorized_radiation_state), intent(out) :: value
    logical, intent(out) :: available
    value%sqrt_shat = 100d0
    ! Exercise the massive-sector fallback where no counterevent exists.
    available = slot == 3 .and. block == 0
  end subroutine fetch_factorized_radiation_state
end module factorized_phase_space

module decay_chain_parameters
  implicit none
contains
  double precision function decay_renormalization_scale(pdg, factor_index)
    integer, intent(in) :: pdg
    integer, intent(in), optional :: factor_index
    integer :: selected
    selected = 1
    if (present(factor_index)) selected = factor_index
    decay_renormalization_scale = 20d0*dble(selected) + 0d0*pdg
  end function decay_renormalization_scale
  logical function use_decayed_production_ren_scale_momenta()
    use_decayed_production_ren_scale_momenta = .false.
  end function use_decayed_production_ren_scale_momenta
end module decay_chain_parameters

module run_state
  implicit none
  double precision :: scale = 0d0
  double precision :: q2fact(2) = 0d0
  double precision :: mur2_current = 0d0
  double precision :: muf12_current = 0d0
  double precision :: muf22_current = 0d0
  double precision :: qes2_current = 0d0
  logical :: fixed_ren_scale = .false.
  logical :: fixed_fac_scale = .false.
  double precision :: mur_over_ref = 1d0
  double precision :: muf1_over_ref = 0.5d0
  double precision :: mur_ref_fixed = 91d0
  double precision :: muf1_ref_fixed = 91d0
  integer :: dynamical_scale_choice = 2
end module run_state

module fixed_order_user_hooks
  implicit none
contains
  double precision function fixed_user_scale(reference, scale_id)
    double precision, intent(in) :: reference
    character(len=*), intent(out) :: scale_id
    fixed_user_scale = reference
    scale_id = 'fixed'
  end function fixed_user_scale
end module fixed_order_user_hooks

module multiplicative_kinematics
  implicit none
contains
  subroutine factorized_production_scale_sums(slots, sum_et, sum_mt)
    integer, intent(in) :: slots(0:)
    double precision, intent(out) :: sum_et, sum_mt
    sum_et = 10d0
    sum_mt = 20d0 + 0d0*slots(0)
  end subroutine factorized_production_scale_sums
  subroutine factorized_visible_scale_sums(slots, sum_et, sum_mt)
    integer, intent(in) :: slots(0:)
    double precision, intent(out) :: sum_et, sum_mt
    sum_et = 30d0
    sum_mt = 40d0 + 0d0*slots(0)
  end subroutine factorized_visible_scale_sums
end module multiplicative_kinematics

module fnlo_process_common
  implicit none
  integer, parameter :: soft_counterevent = 0
  integer, parameter :: real_event = 3
end module fnlo_process_common

module bridge_state
  implicit none
  double precision :: bridge_mur = 0d0
  double precision :: bridge_g = 0d0
  double precision :: bridge_qes2 = 0d0
end module bridge_state

integer function sdm_multiplicative_block_count()
  sdm_multiplicative_block_count = 3
end function sdm_multiplicative_block_count

integer function sdm_multiplicative_physical_block(position)
  integer, intent(in) :: position
  integer, parameter :: values(3) = (/0, 1, 2/)
  sdm_multiplicative_physical_block = values(position)
end function sdm_multiplicative_physical_block

integer function sdm_multiplicative_block_pdg(block)
  integer, intent(in) :: block
  integer, parameter :: values(0:2) = (/0, 6, 24/)
  sdm_multiplicative_block_pdg = values(block)
end function sdm_multiplicative_block_pdg

integer function sdm_multiplicative_born_qcd_power(block)
  integer, intent(in) :: block
  integer, parameter :: values(0:2) = (/4, 2, 0/)
  sdm_multiplicative_born_qcd_power = values(block)
end function sdm_multiplicative_born_qcd_power

subroutine set_model_ren_scale_bridge(mur, g_value)
  use bridge_state
  double precision, intent(in) :: mur, g_value
  bridge_mur = mur
  bridge_g = g_value
end subroutine set_model_ren_scale_bridge

subroutine set_model_qes_scale_bridge(qes_squared)
  use bridge_state
  double precision, intent(in) :: qes_squared
  bridge_qes2 = qes_squared
end subroutine set_model_qes_scale_bridge
'''
        program = r'''
program test_scale_state
  use run_state
  use bridge_state
  use multiplicative_scale_state
  implicit none
  double precision :: value, logs_r(0:5), logs_f(0:5)
  double precision :: rescaling(0:5,0:1), production_mu2_r
  double precision :: production_mu2_f
  integer :: slots(0:5), factor_indices(0:5)

  call initialize_multiplicative_scale_references()
  if (abs(multiplicative_reference_scale_squared(0) - 10000d0) > &
      1d-12) stop 1
  if (abs(multiplicative_reference_scale_squared(1) - 400d0) > &
      1d-12) stop 2
  if (multiplicative_reference_coupling(0) <= 0d0 .or. &
      multiplicative_reference_coupling(1) <= 0d0) stop 3
  call activate_multiplicative_block_reference(1)
  if (abs(scale - 20d0) > 1d-12) stop 4
  if (maxval(abs(q2fact - 400d0)) > 1d-12) stop 5
  if (abs(bridge_mur - 20d0) > 1d-12 .or. &
      abs(bridge_qes2 - 400d0) > 1d-12) stop 6
  if (abs(multiplicative_block_scale_logarithm(1, 1600d0) - &
      log(4d0)) > 1d-12) stop 7
  value = multiplicative_block_coupling_rescaling(2, 1, 40d0)
  if (value <= 0d0 .or. abs(value - 1d0) < 1d-8) stop 8
  slots = 0
  factor_indices = 1
  factor_indices(1) = 2
  call build_multiplicative_scale_tables( &
       slots, logs_r, logs_f, rescaling, production_mu2_r, &
       production_mu2_f, decay_factor_indices=factor_indices)
  if (abs(production_mu2_r - 400d0) > 1d-12 .or. &
      abs(production_mu2_f - 100d0) > 1d-12) stop 9
  if (abs(logs_r(0) - log(400d0/10000d0)) > 1d-12) stop 10
  if (abs(logs_f(0) - log(100d0/10000d0)) > 1d-12) stop 11
  if (abs(logs_r(1) - log(1600d0/400d0)) > 1d-12) stop 12
  if (abs(rescaling(1,0) - &
      multiplicative_block_coupling_rescaling(1,0,40d0)) > 1d-12) stop 13
  call build_multiplicative_scale_tables( &
       slots, logs_r, logs_f, rescaling, production_mu2_r, &
       production_mu2_f, production_dynamic_choice=1)
  if (abs(production_mu2_r - 100d0) > 1d-12 .or. &
      abs(production_mu2_f - 25d0) > 1d-12) stop 14
end program test_scale_state
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

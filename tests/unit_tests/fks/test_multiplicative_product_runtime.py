################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of MadGraph5_aMC@NLO.
#
################################################################################

"""Compiled checks for the lazy multiplicative-product Fortran runtime."""

from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile

import tests.unit_tests as unittest

from madgraph import MG5DIR


class TestMultiplicativeProductRuntime(unittest.TestCase):

    def test_three_real_tensor_event_is_mapped_and_contracted(self):
        if not shutil.which('gfortran'):
            self.skipTest('gfortran is unavailable')
        metadata = """\
FORMAT 1
PRESCRIPTION STAGEWISE_NLO_PRODUCT
ENUMERATION CARTESIAN_LAZY
COUNTEREVENTS TENSOR_PRODUCT
STAGES 3
SECTORS 18
FIRST_ORDER_SECTORS 5
MAX_RADIATIONS 3
STAGE 1 PRODUCTION PRODUCTION 0 0 0 1 1 3 1
VIRTUAL_ORDER 1 1 3 0
CHOICE 1 1 BORN 0 0 0 0 0 0 0
CHOICE 1 2 FINITE 0 0 0 0 0 0 0
CHOICE 1 3 REAL 1 1 6 5 5 1 1
STAGE 2 DECAY_1 NLO_DECAY 6 1 1 0 1 2 0
CHOICE 2 1 BORN 0 0 0 0 0 0 0
CHOICE 2 2 REAL 1 1 4 3 3 1 0
STAGE 3 DECAY_2 NLO_DECAY -6 1 2 1 1 3 1
VIRTUAL_ORDER 3 1 1 2
CHOICE 3 1 BORN 0 0 0 0 0 0 0
CHOICE 3 2 FINITE 0 0 0 0 0 0 0
CHOICE 3 3 REAL 1 1 5 4 4 0 1
END
"""
        dimensions = """\
module process_dimensions
  implicit none
  integer :: nsplitorders = 2
contains
  subroutine validate_process_dimensions()
  end subroutine validate_process_dimensions
end module process_dimensions
"""
        driver = r"""
program test_product_runtime
  use iso_fortran_env, only: int64, real64
  use multiplicative_product
  implicit none

  type(product_event_descriptor) :: event
  integer, allocatable :: choices(:), slots(:)
  integer :: sign, orders(2), map_calls, map_order(3)
  real(real64) :: coordinates(11)
  real(real64) :: seed_momenta(0:3, 6), mapped_momenta(0:3, 6)
  real(real64) :: seed_masses(6), mapped_masses(6)
  real(real64) :: jacobian, carrier_value, kernels(3), weight
  logical :: pass

  map_calls = 0
  map_order = 0
  if (command_argument_count() > 0) then
    call initialize_multiplicative_product()
    if (has_multiplicative_product()) error stop 8
    if (multiplicative_product_stage_count() /= 0 .or. &
        multiplicative_product_sector_count() /= 0_int64) error stop 9
    stop
  end if
  call initialize_multiplicative_product()
  if (.not. has_multiplicative_product()) error stop 10
  if (multiplicative_product_stage_count() /= 3) error stop 11
  if (multiplicative_product_sector_count() /= 18_int64) error stop 12
  if (multiplicative_product_first_order_sector_count() /= 5_int64) &
       error stop 13
  if (multiplicative_product_max_radiations() /= 3) error stop 14

  call multiplicative_product_virtual_order(3, 1, orders)
  if (any(orders /= (/1, 2/))) error stop 15
  call decode_multiplicative_sector(18_int64, choices)
  if (any(choices /= (/3, 2, 3/))) error stop 16
  if (multiplicative_counterevent_count(18_int64) /= 16_int64) &
       error stop 17
  call decode_multiplicative_counterevent(18_int64, 16_int64, slots, sign)
  if (any(slots /= (/product_slot_soft_collinear, product_slot_soft, &
                     product_slot_collinear/))) error stop 18
  if (sign /= 1) error stop 19

  coordinates = (/10._real64, 11._real64, &
       .2_real64, .3_real64, .4_real64, &
       .5_real64, .6_real64, .7_real64, &
       .8_real64, .9_real64, 1._real64/)
  if (multiplicative_phase_space_dimension(18_int64, 2) /= 11) &
       error stop 20
  call build_multiplicative_event(&
       18_int64, 16_int64, 2, coordinates, event)
  if (event%perturbative_order /= 3 .or. event%real_order /= 3 .or. &
      event%finite_order /= 0 .or. event%inclusion_sign /= 1) error stop 21
  if (any(event%born_coordinates /= (/10._real64, 11._real64/))) &
       error stop 22
  if (event%stages(1)%xi /= 0._real64 .or. &
      event%stages(1)%y /= 1._real64 .or. &
      event%stages(2)%xi /= 0._real64 .or. &
      event%stages(2)%y /= .6_real64 .or. &
      event%stages(3)%xi /= .8_real64 .or. &
      event%stages(3)%y /= 1._real64) error stop 23

  seed_momenta = 0._real64
  seed_momenta(0, 1) = 100._real64
  seed_masses = 0._real64
  call evaluate_product_counterevent(event, seed_momenta, seed_masses, &
       map_stage, evaluate_carrier, evaluate_kernel, mapped_momenta, &
       mapped_masses, jacobian, carrier_value, kernels, weight, pass)
  if (.not. pass) error stop 24
  if (map_calls /= 3 .or. any(map_order /= (/1, 2, 3/))) error stop 25
  if (jacobian /= 6._real64 .or. carrier_value /= 2._real64) error stop 26
  if (any(kernels /= (/14._real64, 22._real64, 33._real64/))) &
       error stop 27
  if (weight /= 121968._real64) error stop 28
  if (mapped_momenta(0, 1) /= 101._real64 .or. &
      mapped_momenta(0, 2) /= 2._real64 .or. &
      mapped_momenta(0, 3) /= 3._real64) error stop 29

  map_calls = 0
  call evaluate_product_counterevent(event, seed_momenta, seed_masses, &
       fail_second_map, evaluate_carrier, evaluate_kernel, mapped_momenta, &
       mapped_masses, jacobian, carrier_value, kernels, weight, pass)
  if (pass) error stop 30
  if (any(mapped_momenta /= seed_momenta) .or. &
      any(mapped_masses /= seed_masses)) error stop 31
  if (jacobian /= 0._real64 .or. carrier_value /= 0._real64 .or. &
      weight /= 0._real64) error stop 32

contains

  subroutine map_stage(stage, momenta, masses, local_jacobian, local_pass)
    type(product_stage_event), intent(in) :: stage
    real(real64), intent(inout) :: momenta(0:, :), masses(:)
    real(real64), intent(out) :: local_jacobian
    logical, intent(out) :: local_pass
    map_calls = map_calls + 1
    map_order(map_calls) = stage%stage_id
    momenta(0, stage%stage_id) = &
         momenta(0, stage%stage_id) + real(stage%stage_id, real64)
    masses(stage%stage_id) = real(stage%stage_id, real64)/10._real64
    local_jacobian = real(stage%stage_id, real64)
    local_pass = .true.
  end subroutine map_stage

  subroutine fail_second_map(stage, momenta, masses, local_jacobian, &
                             local_pass)
    type(product_stage_event), intent(in) :: stage
    real(real64), intent(inout) :: momenta(0:, :), masses(:)
    real(real64), intent(out) :: local_jacobian
    logical, intent(out) :: local_pass
    call map_stage(stage, momenta, masses, local_jacobian, local_pass)
    if (stage%stage_id == 2) local_pass = .false.
  end subroutine fail_second_map

  subroutine evaluate_carrier(active_event, momenta, masses, value, &
                              local_pass)
    type(product_event_descriptor), intent(in) :: active_event
    real(real64), intent(in) :: momenta(0:, :), masses(:)
    real(real64), intent(out) :: value
    logical, intent(out) :: local_pass
    value = 2._real64
    local_pass = active_event%real_order == 3 .and. &
         all(momenta(0, 1:3) > 0._real64) .and. &
         all(masses(1:3) > 0._real64)
  end subroutine evaluate_carrier

  subroutine evaluate_kernel(stage, momenta, masses, value, local_pass)
    type(product_stage_event), intent(in) :: stage
    real(real64), intent(in) :: momenta(0:, :), masses(:)
    real(real64), intent(out) :: value
    logical, intent(out) :: local_pass
    value = real(10*stage%stage_id + stage%slot, real64)
    local_pass = momenta(0, stage%stage_id) > 0._real64 .and. &
         masses(stage%stage_id) > 0._real64
  end subroutine evaluate_kernel

end program test_product_runtime
"""

        runtime = os.path.join(
            MG5DIR, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_product.f90')
        compiler = shutil.which('gfortran')
        with tempfile.TemporaryDirectory() as workdir:
            dimensions_path = os.path.join(workdir, 'process_dimensions.f90')
            driver_path = os.path.join(workdir, 'test_product_runtime.f90')
            metadata_path = os.path.join(
                workdir, 'multiplicative_product_info.dat')
            executable = os.path.join(workdir, 'test_product_runtime')
            with open(dimensions_path, 'w') as stream:
                stream.write(dimensions)
            with open(driver_path, 'w') as stream:
                stream.write(driver)
            with open(metadata_path, 'w') as stream:
                stream.write(metadata)
            subprocess.check_call([
                compiler, '-std=f2008', '-Wall', '-Wextra',
                '-Wno-compare-reals', '-fcheck=all',
                '-J', workdir, dimensions_path, runtime, driver_path,
                '-o', executable], cwd=workdir)
            subprocess.check_call([executable], cwd=workdir)
            empty_dir = os.path.join(workdir, 'without_metadata')
            os.mkdir(empty_dir)
            subprocess.check_call([executable, 'absent'], cwd=empty_dir)
            invalid_dir = os.path.join(workdir, 'invalid_metadata')
            os.mkdir(invalid_dir)
            with open(os.path.join(
                    invalid_dir, 'multiplicative_product_info.dat'),
                    'w') as stream:
                stream.write(metadata.replace('SECTORS 18', 'SECTORS 17'))
            with open(os.devnull, 'w') as devnull:
                status = subprocess.call(
                    [executable], cwd=invalid_dir,
                    stdout=devnull, stderr=devnull)
            self.assertNotEqual(status, 0)


if __name__ == '__main__':
    unittest.main()

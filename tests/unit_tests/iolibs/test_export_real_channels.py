"""Production real-diagram proposals retain undecayed core topologies."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph.fks.fks_helas_objects import FKSHelasMultiProcess
from madgraph.interface.master_interface import MasterCmd
from madgraph.iolibs.export_fks import ProcessExporterFortranFKS


def generate(process):
    command = MasterCmd()
    command.exec_cmd('import model loop_sm', printcmd=False, precmd=True)
    command.exec_cmd('generate ' + process, printcmd=False, precmd=True)
    return command


class TestProductionRealChannels(unittest.TestCase):

    def test_undecayed_core_and_lo_capability(self):
        command = generate('u d~ > w+ [real=QCD]')
        matrix_element = FKSHelasMultiProcess(
            command._fks_multi_proc)['matrix_elements'][0]
        self.assertIsNone(matrix_element.spin_density_plan)
        groups = ProcessExporterFortranFKS.production_real_channel_data(
            matrix_element)
        self.assertEqual(len(groups), 3)
        self.assertTrue(all(group['nexternal'] == 4 for group in groups))

        command = generate('u d~ > w+ [LOonly=QCD]')
        matrix_element = FKSHelasMultiProcess(
            command._fks_multi_proc)['matrix_elements'][0]
        self.assertEqual(ProcessExporterFortranFKS.production_real_channel_data(
            matrix_element), [])
        exporter = ProcessExporterFortranFKS.__new__(ProcessExporterFortranFKS)
        previous = os.getcwd()
        with tempfile.TemporaryDirectory() as temporary:
            try:
                os.chdir(temporary)
                exporter.write_production_real_channels(matrix_element, None)
                self.assertEqual(Path('hard_production_channels.dat').read_text(), '')
                self.assertEqual(Path('hard_symfact.dat').read_text(), '')
                source = Path('spin_density_real_channels.f').read_text()
                self.assertIn('SDM_HARD_CHANNEL_COUNT=0', source)
                self.assertIn('SDM_REAL_CORE_PROCESS=0', source)
                compiler = shutil.which('gfortran')
                if compiler:
                    subprocess.run([compiler, '-c', '-ffixed-line-length-none',
                        'spin_density_real_channels.f'], check=True,
                        capture_output=True, text=True)
            finally:
                os.chdir(previous)

    def test_core_diagrams_exclude_decay_insertions_and_duplicate_sectors(self):
        command = generate(
            'u u~ > t t~ [real=QCD], (t > w+ b [real=QCD])')
        matrix_element = FKSHelasMultiProcess(
            command._fks_multi_proc)['matrix_elements'][0]
        groups = ProcessExporterFortranFKS.production_real_channel_data(
            matrix_element)
        self.assertEqual([group['matrix_id'] for group in groups], [1, 2, 3])
        self.assertEqual([group['configurations'] for group in groups],
                         [[1, 2, 3, 4], [5], [6]])
        self.assertEqual([group['nexternal'] for group in groups], [5, 5, 5])
        for group in groups:
            self.assertEqual(group['ninitial'], 2)
            self.assertEqual(len(group['channels']), 5)
            self.assertNotIn(24, [leg['id'] for leg in
                group['matrix_element']['processes'][0]['legs']])
            for channel in group['channels']:
                self.assertEqual(channel['nbranch'], 3)
                self.assertEqual(sorted(branch['branch'] for branch in
                                        channel['branches']), [1, 2, 3])
                for branch in channel['branches']:
                    self.assertEqual(len(branch['daughters']), 2)
                    self.assertTrue(all(-3 <= leg <= 5 and leg != 0
                                        for leg in branch['daughters']))
                terminal = channel['branches'][-1]
                self.assertEqual(terminal['branch'], 3)
                self.assertEqual(terminal['propagator'], 0)

    def test_exported_channel_tables_and_topology_abi(self):
        command = generate('u d~ > w+ [real=QCD], (w+ > e+ ve)')
        previous = os.getcwd()
        with tempfile.TemporaryDirectory() as temporary:
            try:
                process = Path(temporary)/'PROC'
                command.exec_cmd('output fNLO %s' % process,
                                 printcmd=False, precmd=True)
                directory = next((process/'SubProcesses').glob('P*'))
                for name in ['hard_real_kinematics.f90',
                             'production_hard_split.f90']:
                    self.assertTrue((directory/name).is_symlink())
                    self.assertTrue((directory/name).is_file())
                rows = [list(map(int, row.split())) for row in
                    (directory/'hard_production_channels.dat').read_text().splitlines()]
                self.assertEqual(rows, [[1, 1, 1], [2, 1, 2],
                                        [3, 3, 1], [4, 3, 2],
                                        [5, 4, 1], [6, 4, 2]])
                self.assertEqual((directory/'hard_symfact.dat').read_text(),
                                 ''.join('%d.3 1\n' % i for i in range(1, 7)))
                # The existing runtime creates symfact.dat from its Born
                # configurations; exporting hard maps must not create it.
                self.assertFalse((directory/'symfact.dat').exists())
                source = (directory/'spin_density_real_channels.f').read_text()
                self.assertIn('SDM_HARD_CHANNEL_COUNT=6', source)
                self.assertIn('PARAMETER (NEXTERNAL=4,NCOMB=', source)
                self.assertIn('AMP2(1)=AMP2(1)+AMPLITUDES(1)*', source)
                self.assertNotIn('RHO(', source)
                compiler = shutil.which('gfortran')
                if compiler is None:
                    self.skipTest('gfortran is required for the generated ABI check')
                harness = Path(temporary)/'check.f90'
                harness.write_text('''program check
  implicit none
  integer :: sdm_hard_channel_count,sdm_hard_channel_fks
  integer :: sdm_hard_channel_config,sdm_real_core_process
  integer :: sdm_real_core_channel_count,sdm_real_core_leg_count
  integer :: tree(2,9),sp(9),tp(9),pow(9),nleg,ninitial,nb,id
  double precision :: mass(9),width(9)
  if(sdm_hard_channel_count()/=6) stop 1
  if(sdm_real_core_process(1)/=sdm_real_core_process(2)) stop 2
  if(sdm_real_core_channel_count(1)/=2) stop 3
  if(sdm_real_core_leg_count(1)/=4) stop 4
  do id=1,6
    call sdm_real_core_channel_topology(sdm_hard_channel_fks(id), &
      sdm_hard_channel_config(id),nleg,ninitial,nb,tree,sp,tp,mass,width,pow)
    if(nleg/=4.or.ninitial/=2.or.nb/=2) stop 5
    if(any(tree(:,1:nb)==0)) stop 6
    if(mass(nb)/=0d0.or.width(nb)/=0d0.or.tp(nb)/=0) stop 7
  enddo
  if(sdm_real_core_process(99)/=0) stop 8
  print *, 'PASS'
end program
double precision function get_mass_from_id(pdg)
  integer pdg
  get_mass_from_id=abs(dble(pdg))
end function
double precision function get_width_from_id(pdg)
  integer pdg
  get_width_from_id=0d0
end function
''')
                obj = Path(temporary)/'channels.o'
                subprocess.run([compiler, '-c', '-ffixed-line-length-none',
                    '-ffunction-sections', '-I.', '-I../../Source/MODEL',
                    'spin_density_real_channels.f', '-o', str(obj)],
                    cwd=directory, check=True, capture_output=True, text=True)
                executable = Path(temporary)/'check'
                subprocess.run([compiler, '-Wl,--gc-sections', str(harness),
                    str(obj), '-o', str(executable)], check=True,
                    capture_output=True, text=True)
                result = subprocess.run([str(executable)], check=True,
                                        capture_output=True, text=True)
                self.assertIn('PASS', result.stdout)
            finally:
                os.chdir(previous)


if __name__ == '__main__':
    unittest.main()

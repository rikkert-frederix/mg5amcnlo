"""The optional current proposal changes only integration coordinates."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.fks.fks_decay import decay_card_text
import tests.unit_tests.fks.test_fnlo_decay_card as card_fixture


class TestCurrentProposal(unittest.TestCase):
    def test_current_map_and_rejected_cores(self):
        compiler=shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root=Path(MG5DIR)
        fixture=root/'tests/input_files/fks_decay'
        with tempfile.TemporaryDirectory(prefix='fnlo_current_proposal_') as directory:
            executable=Path(directory)/'check'
            sources=[fixture/'phase_space_test_dimensions.f90',root/'Template/fNLO/Source/kin_functions.f90',
                     root/'Template/fNLO/SubProcesses/phase_space_kinematics.f90',
                     root/'Template/fNLO/SubProcesses/factorized_block_kinematics.f90',
                     fixture/'current_proposal_checks.f90']
            result=subprocess.run([compiler,'-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
                                   '-o',str(executable),*map(str,sources)],cwd=directory,
                                  capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            result=subprocess.run([str(executable)],cwd=directory,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('TRANSFORMED_POINT_CHECKS 80',result.stdout)
            self.assertIn('ORDERING_CHECKS 24',result.stdout)
            for mode,message in [('bad_flavour','mismatched charged-current'),
                                 ('bad_charge','mismatched charged-current'),
                                 ('bad_duplicate','ambiguous W_CURRENT'),
                                 ('bad_mass','massless current daughters'),
                                 ('bad_count','four-particle production core')]:
                result=subprocess.run([str(executable),mode],cwd=directory,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn(message,result.stdout)

    def test_card_defaults_validation_and_runtime(self):
        card_fixture.TestFNLODecayCard.setUpClass()
        fixture=card_fixture.TestFNLODecayCard()
        try:
            base=fixture.card()
            explicit_flat=fixture.card(production_phase_space_sampling='FLAT')
            self.assertEqual(base,explicit_flat)
            default=fixture.run_card(base,'proposal')['PROPOSAL']
            self.assertEqual(default,['F','0.0000000000000000E+00','0.0000000000000000E+00'])
            options=dict(production_phase_space_sampling='W_CURRENT',production_sampling_mass=80.385,
                         production_sampling_width=2.097673562052797)
            card=fixture.card(**options)
            row=fixture.run_card(card,'proposal')['PROPOSAL']
            self.assertEqual(row[0],'T')
            self.assertEqual(list(map(float,row[1:])),[80.385,2.097673562052797])
            self.assertEqual(fixture.run_card(base),fixture.run_card(card))
            for changes in ({'production_phase_space_sampling':'BAD'},
                            {'production_sampling_mass':None},{'production_sampling_width':0.},
                            {'production_sampling_width':float('nan')},
                            {'production_phase_space_sampling':'FLAT'}):
                with self.assertRaises(ValueError):
                    fixture.card(**dict(options,**changes))
            for text,message in [(card+'\nW_CURRENT = production_phase_space_sampling','duplicate'),
                                 (card+'\n3. = production_sampling_width','duplicate'),
                                 (base+'\nW_CURRENT = production_phase_space_sampling','requires explicit'),
                                 (base+'\n2. = production_sampling_width','require W_CURRENT'),
                                 (card.replace('W_CURRENT =','BAD ='),'FLAT or W_CURRENT'),
                                 (card.replace('2.0976735620527971e+00 = production_sampling_width',
                                               '-1. = production_sampling_width'),'positive and finite')]:
                fixture.run_card(text,'proposal',error=message)
        finally:
            card_fixture.TestFNLODecayCard.doClassCleanups()


if __name__=='__main__':
    unittest.main()

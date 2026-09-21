import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from campaign import digest, import_grid


class TestGridRecovery(unittest.TestCase):
    def fixture(self,root):
        study = root/'study'
        old,new = [study/'processes'/name for name in ('old','new')]
        manifest = dict(variant='S',w_treatment='onshell',production_scale='core-w-ht-half',
                        decay_bottom_mass=0.,decay_scale_grouping='separate',
                        topology_hashes={'P':'same'},internal_width_hashes={'P':'same'},
                        settings={'iseed':1,'lhaid':331700,'ebeam1':6500.},hashes={})
        for process in (old,new):
            archive = process/'study_cards'/'S'
            archive.mkdir(parents=True)
            for name in ('param_card.dat','decay_card.dat','FO_analyse_card.dat'):
                (archive/name).write_text(name)
                manifest['hashes'][name] = digest(archive/name)
            for relative in ('FixedOrderAnalysis/HwU.f90',
                             'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90',
                             'SubProcesses/mint_module.f90',
                             'SubProcesses/phase_space_kinematics.f90',
                             'SubProcesses/factorized_block_kinematics.f90',
                             'SubProcesses/decay_chain_parameters.f90',
                             'SubProcesses/decay_chain_kinematics.f90',
                             'SubProcesses/nlo_decay_kinematics.f90'):
                path = process/relative
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(relative)
            (archive/'manifest.json').write_text(json.dumps(manifest))
        (old/'study_cards/S/execution.json').write_text(json.dumps({'source_hashes':{}}))
        (old/'SubProcesses/job_status.pkl').write_bytes(b'trusted local job fixture')
        grid = old/'SubProcesses/P0_test/all_G1'
        grid.mkdir(parents=True)
        for name in ('mint_grids','grid.MC_integer','res_0.dat','res.dat','results.dat',
                     'contribution_results_0.dat','contribution_results.dat','log.txt','MADatNLO.HwU'):
            (grid/name).write_text(name)
        return study,old,new,grid

    def test_complete_grid_and_component_ledgers_are_transferred(self):
        with tempfile.TemporaryDirectory() as temporary:
            study,old,new,grid = self.fixture(Path(temporary))
            with patch('campaign.STUDY',study):
                result = import_grid(new,new/'study_cards/S',old/'study_cards/S/manifest.json')
            self.assertEqual(len(result['files']),9)
            self.assertIn('SubProcesses/P0_test/all_G1/contribution_results.dat',result['files'])
            self.assertIn('SubProcesses/P0_test/all_G1/results.dat',result['files'])
            self.assertFalse((new/'SubProcesses/P0_test/all_G1/MADatNLO.HwU').exists())
            for relative,checksum in result['files'].items():
                self.assertEqual(digest(old/relative),checksum)
                self.assertEqual(digest(new/relative),checksum)

    def test_missing_summary_is_rejected_before_any_grid_copy(self):
        with tempfile.TemporaryDirectory() as temporary:
            study,old,new,grid = self.fixture(Path(temporary))
            (grid/'contribution_results.dat').unlink()
            with patch('campaign.STUDY',study),self.assertRaises(ValueError):
                import_grid(new,new/'study_cards/S',old/'study_cards/S/manifest.json')
            self.assertFalse((new/'SubProcesses/job_status.pkl').exists())

    def test_changed_inputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            study,old,new,grid = self.fixture(Path(temporary))
            (new/'study_cards/S/param_card.dat').write_text('changed weak inputs')
            with patch('campaign.STUDY',study),self.assertRaises(ValueError):
                import_grid(new,new/'study_cards/S',old/'study_cards/S/manifest.json')
            self.assertFalse((new/'SubProcesses/job_status.pkl').exists())

    def test_changed_map_and_private_runtime_are_rejected(self):
        for mode in ('map', 'parameters', 'private_runtime'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as temporary:
                study,old,new,grid = self.fixture(Path(temporary))
                path = new/'study_cards/S/manifest.json'
                manifest = json.loads(path.read_text())
                if mode == 'map':
                    manifest['production_sampling'] = 'w-current'
                elif mode == 'parameters':
                    manifest['production_sampling_parameters'] = dict(mass_GeV=80.,width_GeV=2.)
                else:
                    (new/'SubProcesses/factorized_block_kinematics.f90').write_text('new map')
                path.write_text(json.dumps(manifest))
                with patch('campaign.STUDY',study), self.assertRaisesRegex(ValueError, 'sampling-|runtime changed'):
                    import_grid(new,new/'study_cards/S',old/'study_cards/S/manifest.json')
                self.assertFalse((new/'SubProcesses/job_status.pkl').exists())


if __name__ == '__main__':
    unittest.main()

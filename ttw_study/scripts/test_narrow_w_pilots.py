import contextlib
import copy
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from campaign import ROOT, STUDY, digest
from narrow_w_inputs import COMMON_SCALE, rescaling
from run_narrow_w_pilots import cases, configure, decay_card, width_for


class TestNarrowWPilots(unittest.TestCase):
    def inputs(self):
        path=STUDY/'inputs/narrow_w_common_inputs_v1.json'
        return path,json.loads(path.read_text())

    def test_three_current_normalizations_are_removed_once(self):
        for factor in (1.,.2,.05):
            self.assertAlmostEqual(rescaling(factor)/factor**3,1.)
            self.assertAlmostEqual((.108345431330546/factor)**3*rescaling(factor),.108345431330546**3)
        for bad in (0.,-1.,2.,float('nan'),float('inf')):
            with self.assertRaises(ValueError):
                rescaling(bad)

    def test_full_mode_width_order_grid_and_distinct_seeds(self):
        rows=cases()
        self.assertEqual(len(rows),21)
        self.assertEqual({r['seed'] for r in rows},set(range(65001,65022)))
        for variant in ('LO','S','Pi'):
            selected=[r for r in rows if r['variant']==variant]
            self.assertEqual(len(selected),7)
            self.assertEqual(selected[0]['w_treatment'],'onshell')
            for mode in ('top-bw','all-bw'):
                self.assertEqual({r['width_factor'] for r in selected if r['w_treatment']==mode},{1.,.2,.05})

    def test_common_decay_and_width_reference_are_not_left_at_mt(self):
        _,inputs=self.inputs()
        for row in cases():
            width=width_for(row,inputs)
            card=decay_card(row,width)
            self.assertIn('%.16e = decay_ren_scale(6)' % COMMON_SCALE,card)
            self.assertIn('%.16e = nlo_decay_width(6)' % width['gamma_nlo'],card)
            self.assertIn('AUTO = decay_width_scale_mode(6)',card)
            self.assertEqual('W_CURRENT = production_phase_space_sampling' in card,row['w_treatment']=='all-bw')
            self.assertEqual(' = lo_decay_width(24)' in card,row['w_treatment']!='all-bw')
        bad=copy.deepcopy(inputs)
        bad['rows'].append(copy.deepcopy(bad['rows'][0]))
        with self.assertRaises(ValueError):
            width_for(cases()[0],bad)

    def test_actual_setup_writes_distinct_complete_common_scale_archives(self):
        import tests.unit_tests.fks.test_ttw_product_tools as fixtures
        from madgraph.various.banner import RunCardNLO
        from models.check_param_card import ParamCard
        inputs_path,inputs=self.inputs()
        source=STUDY/'processes/TTWplus_onshell_eemu_mb0p0_both_rng_v1/study_cards/S_onshell_core-w-ht-half_separate_57003'
        with tempfile.TemporaryDirectory(prefix='narrow_cards_') as temporary:
            for row in cases():
                process=Path(temporary)/str(row['seed'])
                cards=process/'Cards'
                cards.mkdir(parents=True)
                (process/'Events').mkdir()
                sub=process/'SubProcesses'
                sub.mkdir()
                (process/'FixedOrderAnalysis').mkdir()
                for name in ('run_card.dat','decay_card.dat','param_card.dat','FO_analyse_card.dat'):
                    shutil.copy2(source/name,cards/name)
                shutil.copy2(source/'run_card.dat',cards/'run_card_default.dat')
                for name in ('decay_chain_parameters.f90','setscales.f90','decay_chain_scales.f90',
                             'factorized_block_kinematics.f90','phase_space_kinematics.f90',
                             'decay_chain_kinematics.f90','nlo_decay_kinematics.f90'):
                    shutil.copy2(ROOT/'Template/fNLO/SubProcesses'/name,sub/name)
                shutil.copy2(ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90',
                             process/'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90')
                fixtures.TestTTWProductTools.topology_fixture(process,row['w_treatment'])
                grid_points=1000 if row['seed']==65001 else 200
                with contextlib.redirect_stdout(io.StringIO()):
                    name,archive=configure(process,row,inputs,inputs_path,.03,grid_points=grid_points)
                manifest=json.loads((archive/'manifest.json').read_text())
                base=Path(manifest['base_configuration'])
                self.assertEqual(digest(base),manifest['base_configuration_sha256'])
                self.assertFalse(list((process/'Events').iterdir()))
                self.assertEqual(manifest['top_width_reference_scale'],COMMON_SCALE)
                self.assertEqual(manifest['generated_event_weight_rescaling'],1.)
                self.assertEqual(manifest['comparison_cross_section_rescaling'],row['width_factor']**3)
                card=RunCardNLO(str(archive/'run_card.dat'))
                self.assertEqual(card['npoints_fo_grid'],grid_points)
                for flag in ('fixed_ren_scale','fixed_fac_scale','fixed_qes_scale'):
                    self.assertTrue(card[flag])
                for value in ('mur_ref_fixed','muf_ref_fixed','qes_ref_fixed'):
                    self.assertEqual(card[value],COMMON_SCALE)
                self.assertEqual(list(card['lhaid']),[331700])
                param=ParamCard(str(archive/'param_card.dat'))
                width=width_for(row,inputs)
                self.assertEqual(param['decay'].get((24,)).value,width['w_width_GeV'])
                self.assertEqual(param['decay'].get((6,)).value,width['gamma_nlo'])
                self.assertEqual((archive/'decay_card.dat').read_text(),decay_card(row,width))
                for filename in ('param_card.dat','run_card.dat','decay_card.dat','FO_analyse_card.dat'):
                    self.assertEqual(digest(archive/filename),manifest['hashes'][filename])


if __name__=='__main__':
    unittest.main()

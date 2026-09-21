import copy
import unittest

from campaign import STUDY
from production_invariance import compare_subprocess, core_key, fixed_statements, production_fks, production_owner


def fks_fixture():
    lines=['      DATA FKS_I_D / 7,7,6 /','      DATA FKS_J_D / 1,2,3 /',
           '      DATA NEED_COLOR_LINKS_D / .TRUE.,.TRUE.,',
           '     $ .TRUE. /']
    for i in range(1,4):
        lines.extend(['      DATA FKS_J_FROM_I_D(%d,1,0) / 0 /' % i,
                      '      DATA (PARTICLE_TYPE_D(%d,IPOS),IPOS=1,NEXTERNAL) / 3,-3 /' % i,
                      '      DATA (PDG_TYPE_D(%d,IPOS),IPOS=1,NEXTERNAL) / 5,-5 /' % i])
    return '\n'.join(lines)


class TestProductionInvariance(unittest.TestCase):
    def test_core_identity_uses_ordered_pdgs_not_particle_aliases(self):
        source=('FORMAT 4\nCONTEXT 1 BORN 1 5 10\n'
                'CORE_LEG 1 1 -1 I\nCORE_LEG 1 2 2 I\nCORE_LEG 1 3 6 F\n'
                'CORE_LEG 1 4 -6 F\nCORE_LEG 1 5 24 F\n')
        self.assertEqual(core_key(source),((-1,'I'),(2,'I'),(6,'F'),(-6,'F'),(24,'F')))
        self.assertNotEqual(core_key(source),core_key(source.replace('1 1 -1','1 1 -3')))
        for broken in (source.replace('CORE_LEG 1 5 24 F\n',''),
                       source.replace('1 5 24 F','1 5 24 I'),source.replace('1 2 2 I','1 1 2 I')):
            with self.assertRaises(ValueError):
                core_key(broken)

    def test_production_ownership_is_not_inferred_from_total_fks_count(self):
        text=('FORMAT 3\nCOUNT 3\nVIRTUAL_GRIDS 3\n'
              'CONTRIBUTION 1 PRODUCTION 1 6 1 1 0 0 0\nVIRTUAL_GRID 1 1 6 12\n'
              'CONTRIBUTION 2 NLO_DECAY 7 7 7 1 6 1 2\nEND\n')
        owner=production_owner(text)
        self.assertEqual((owner['first'],owner['last']),(1,6))
        for changed in (text.replace('FORMAT 3','FORMAT 2'),
                        text.replace('1 PRODUCTION 1','1 PRODUCTION 2'),
                        text.replace(' 1 1 0 0 0',' 1 1 6 0 0')):
            with self.assertRaises(ValueError):
                production_owner(changed)

    def test_decay_region_changes_are_allowed_but_production_changes_are_visible(self):
        original=fks_fixture()
        reference=production_fks(original,1,2)
        decay_change=original.replace('7,7,6','7,7,9').replace('FKS_J_FROM_I_D(3,1,0) / 0',
                                                               'FKS_J_FROM_I_D(3,1,0) / 4')
        self.assertEqual(reference,production_fks(decay_change,1,2))
        for changed in (original.replace('7,7,6','8,7,6'),
                        original.replace('FKS_J_FROM_I_D(1,1,0) / 0','FKS_J_FROM_I_D(1,1,0) / 4'),
                        original.replace('PDG_TYPE_D(1,IPOS),IPOS=1,NEXTERNAL) / 5,-5',
                                         'PDG_TYPE_D(1,IPOS),IPOS=1,NEXTERNAL) / 4,-4')):
            self.assertNotEqual(reference,production_fks(changed,1,2))

    def test_parser_fails_closed_on_incomplete_or_unknown_records(self):
        original=fks_fixture()
        for changed in (original.replace('FKS_I_D','UNKNOWN'),
                        original.replace('      DATA FKS_I_D / 7,7,6 /',''),
                        original.replace('PDG_TYPE_D(2,','MISSING(2,'),
                        original+'\n      DATA FKS_I_D / 7,7,6 /'):
            with self.assertRaises(ValueError):
                production_fks(changed,1,2)
        with self.assertRaises(ValueError):
            fixed_statements('     $ orphan')

    def test_actual_massless_subprocess_self_comparison(self):
        # This fixture is the preserved physical export, not a generated toy
        # amplitude. It tests inventory/ownership parsing but not mass effects.
        process=STUDY/'processes/TTWplus_all-bw_eemu_mb0p0_both_rng_v1'
        if not process.exists():
            self.skipTest('Physical pilot export unavailable')
        directories=sorted(p for p in (process/'SubProcesses').glob('P*') if p.is_dir())
        self.assertEqual(len(directories),2)
        for directory in directories:
            result=compare_subprocess(directory,directory)
            self.assertGreater(len(result['production_sources_sha256']),30)
            self.assertEqual(result['production_ownership']['last'],6)

    def test_private_model_preserves_all_production_couplings_at_three_scales(self):
        from madgraph.fks.fks_decay_masses import massive_bottom_decay_model
        from models import import_ufo, model_reader
        from models.check_param_card import ParamCard
        production=import_ufo.import_model('loop_sm-no_b_mass')
        unchanged=copy.deepcopy(production)
        massive_bottom_decay_model(production,4.8)
        readers=[model_reader.ModelReader(model) for model in (unchanged,production)]
        source=STUDY/'processes/TTWplus_all-bw_eemu_mb0p0_both_rng_v1/Cards/param_card.dat'
        cards=[ParamCard(str(source)),ParamCard(str(source))]
        cards[1].add_param('decaymass',[5],4.8)
        for scale in (86.25,172.5,345.):
            for reader,card in zip(readers,cards):
                reader.set_parameters_and_couplings(copy.deepcopy(card),scale=scale)
            for name,value in readers[0]['coupling_dict'].items():
                self.assertEqual(value,readers[1]['coupling_dict'][name],(scale,name))
            self.assertEqual(readers[1].get_mass(5),0.)
            self.assertEqual(readers[1]['parameter_dict']['dc_mdl_MB'],4.8)
            self.assertEqual(production.get_nflav(),5)


if __name__=='__main__':
    unittest.main()

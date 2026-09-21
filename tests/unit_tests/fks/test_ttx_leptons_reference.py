from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.various.FO_analyse_card import FOAnalyseCard


class TestTTXLeptonsReference(unittest.TestCase):
    def test_extra_bottom_and_leptonic_measurement(self):
        if not shutil.which('gfortran'):
            self.skipTest('gfortran required')
        name = 'analysis_HwU_pp_ttx_leptons_reference'
        card = FOAnalyseCard('FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE='+name+'.o\n'
                             'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n',testing=True)
        self.assertIn(name+'.o',card.write_card('',fixed_order_only=True))
        root = Path(MG5DIR)
        fixtures = root/'tests/input_files/fks_decay'
        with tempfile.TemporaryDirectory(prefix='ttx_reference_') as build:
            flags = ['gfortran','-O0','-g','-fcheck=all','-ffree-line-length-none',
                     '-ffpe-trap=invalid,zero,overflow']
            for command in (
                flags+['-c',fixtures/'ttw_product_analysis_driver.f90'],
                flags+['-c',root/'Template/fNLO/FixedOrderAnalysis'/(name+'.f90')],
                flags+['-o','check',fixtures/'ttx_leptons_reference_checks.f90',
                       'ttw_product_analysis_driver.o',name+'.o'],
                [Path(build)/'check']):
                result = subprocess.run(list(map(str,command)),cwd=build,text=True,
                                        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                self.assertEqual(result.returncode,0,result.stdout)
            self.assertIn('PASS: lepton-only reference',result.stdout)

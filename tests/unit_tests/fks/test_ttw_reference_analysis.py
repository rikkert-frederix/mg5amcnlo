"""Compile the independent literature selection with the actual FastJet wrapper."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.various.FO_analyse_card import FOAnalyseCard


class TestTTWReferenceAnalysis(unittest.TestCase):
    def test_standalone_abi_and_measurement(self):
        for program in ('gfortran','c++','fastjet-config'):
            if not shutil.which(program):
                self.skipTest(program+' required')
        card = FOAnalyseCard('FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE=analysis_HwU_pp_ttxw_reference.o\n'
                             'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n',testing=True)
        content = card.write_card('',fixed_order_only=True)
        self.assertIn('analysis_HwU_pp_ttxw_reference.o',content)
        self.assertNotIn('analysis_HwU_pp_ttxw_reference_bridge.o',content)
        root = Path(MG5DIR)
        fixture = root/'tests/input_files/fks_decay'
        source = root/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_reference.f90'
        wrapper = root/'Template/fNLO/SubProcesses/fastjetfortran_madfks_full.cc'
        cflags = shlex.split(subprocess.check_output(['fastjet-config','--cxxflags'],text=True))
        libs = shlex.split(subprocess.check_output(['fastjet-config','--libs','--plugins'],text=True))
        libs += ['-lc++' if os.uname().sysname == 'Darwin' else '-lstdc++']
        with tempfile.TemporaryDirectory(prefix='ttw_reference_analysis_') as build:
            def run(command):
                result = subprocess.run(list(map(str,command)),cwd=build,text=True,
                                        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                self.assertEqual(result.returncode,0,result.stdout)
                return result.stdout
            run(['c++','-std=c++11','-c',wrapper,'-o','fastjet.o']+cflags)
            flags = ['gfortran','-O0','-g','-fcheck=all','-ffree-line-length-none',
                     '-ffpe-trap=invalid,zero,overflow']
            run(flags+['-c',fixture/'ttw_product_analysis_driver.f90'])
            run(flags+['-c',source])
            run(flags+['-o','check',fixture/'ttw_reference_analysis_checks.f90',
                       'ttw_product_analysis_driver.o','analysis_HwU_pp_ttxw_reference.o','fastjet.o']+libs)
            self.assertIn('PASS: reference cuts',run([Path(build)/'check']))


if __name__ == '__main__':
    unittest.main()

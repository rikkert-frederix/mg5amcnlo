"""Compile the paper measurement function against the real FastJet wrapper."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.various.FO_analyse_card import FOAnalyseCard


class TestTTWProductAnalysis(unittest.TestCase):
    def test_bridge_selection(self):
        card = FOAnalyseCard('''FO_ANALYSIS_FORMAT=HwU
FO_ANALYSE=analysis_HwU_pp_ttxw_product.o
FO_EXTRALIBS=
FO_EXTRAPATHS=
FO_INCLUDEPATHS=
''', testing=True)
        self.assertIn('analysis_HwU_pp_ttxw_product_bridge.o',
                      card.write_card('', fixed_order_only=True))

    def test_measurement_function(self):
        for program in ('gfortran', 'c++', 'fastjet-config'):
            if not shutil.which(program):
                self.skipTest('%s is required' % program)
        root = Path(MG5DIR)
        fixture = root / 'tests/input_files/fks_decay'
        analysis = root / 'Template/fNLO/FixedOrderAnalysis'
        wrapper = root / 'Template/fNLO/SubProcesses/fastjetfortran_madfks_full.cc'
        cflags = shlex.split(subprocess.check_output(
            ['fastjet-config', '--cxxflags'], text=True))
        libs = shlex.split(subprocess.check_output(
            ['fastjet-config', '--libs', '--plugins'], text=True))
        # Match the platform C++ runtime used by the shipped wrapper.
        libs += ['-lc++' if os.uname().sysname == 'Darwin' else '-lstdc++']
        with tempfile.TemporaryDirectory(prefix='ttw_product_analysis_') as build:
            def run(args, success=True):
                result = subprocess.run(list(map(str, args)), cwd=build,
                                        text=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
                if success:
                    self.assertEqual(result.returncode, 0, result.stdout)
                return result
            run(['c++', '-std=c++11', '-c', wrapper, '-o', 'fastjet.o'] + cflags)
            flags = ['gfortran', '-O0', '-g', '-fcheck=all',
                     '-ffree-line-length-none', '-ffpe-trap=invalid,zero,overflow']
            run(flags + ['-c', fixture / 'ttw_product_analysis_driver.f90'])
            run(flags + ['-c', analysis / 'analysis_HwU_pp_ttxw_product.f90'])
            run(flags + ['-c', analysis / 'analysis_HwU_pp_ttxw_product_bridge.f'])
            run(flags + ['-o', 'check', fixture / 'ttw_product_analysis_checks.f90',
                         'ttw_product_analysis_driver.o',
                         'analysis_HwU_pp_ttxw_product.o',
                         'analysis_HwU_pp_ttxw_product_bridge.o', 'fastjet.o'] + libs)
            self.assertIn('PASS: 210', run([Path(build) / 'check']).stdout)
            bad = run([Path(build) / 'check', 'bad'], success=False)
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn('requires 3 e/mu, 3 neutrinos', bad.stdout)


if __name__ == '__main__':
    unittest.main()

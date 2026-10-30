"""Compile the actual trilepton analysis and bridges with actual FastJet."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.various.FO_analyse_card import FOAnalyseCard


class TestTTWTrileptonAnalysis(unittest.TestCase):
    def test_analysis_selection_includes_both_bridges(self):
        card = FOAnalyseCard('''FO_ANALYSIS_FORMAT=HwU
FO_ANALYSE=analysis_HwU_pp_ttxw_trilepton.o
FO_EXTRALIBS=
FO_EXTRAPATHS=
FO_INCLUDEPATHS=
''', testing=True)
        self.assertIn('analysis_HwU_pp_ttxw_trilepton_bridge.o',
                      card.write_card('', fixed_order_only=True))

    def test_visible_measurements_and_unresolved_limits(self):
        for program in ('gfortran', 'c++', 'fastjet-config'):
            if not shutil.which(program):
                self.skipTest('%s is required' % program)
        root = Path(MG5DIR)
        fixture = root/'tests/input_files/fks_decay'
        analysis = root/'Template/fNLO/FixedOrderAnalysis'
        source = analysis/'analysis_HwU_pp_ttxw_trilepton.f90'
        self.assertNotIn('reconstruct_three_neutrinos', source.read_text())
        cflags = shlex.split(subprocess.check_output(
            ['fastjet-config', '--cxxflags'], text=True))
        libs = shlex.split(subprocess.check_output(
            ['fastjet-config', '--libs', '--plugins'], text=True))
        libs += ['-lc++' if os.uname().sysname == 'Darwin' else '-lstdc++']
        with tempfile.TemporaryDirectory(prefix='ttw_trilepton_analysis_') as build:
            def run(args, success=True):
                result = subprocess.run(list(map(str, args)), cwd=build,
                    text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    timeout=120)
                if success:
                    self.assertEqual(result.returncode, 0, result.stdout)
                return result
            wrapper = root/'Template/fNLO/SubProcesses/fastjetfortran_madfks_full.cc'
            run(['c++', '-std=c++11', '-c', wrapper, '-o', 'fastjet.o']+cflags)
            flags = ['gfortran', '-O0', '-g', '-fcheck=all',
                     '-ffree-line-length-none', '-ffpe-trap=invalid,zero,overflow']
            run(flags+['-c', fixture/'ttw_trilepton_analysis_driver.f90'])
            run(flags+['-c', source])
            run(flags+['-c', analysis/'analysis_HwU_pp_ttxw_trilepton_bridge.f'])
            run(flags+['-o', 'check', fixture/'ttw_trilepton_analysis_checks.f90',
                'ttw_trilepton_analysis_driver.o', 'analysis_HwU_pp_ttxw_trilepton.o',
                'analysis_HwU_pp_ttxw_trilepton_bridge.o', 'fastjet.o']+libs)
            self.assertIn('PASS: trilepton categories', run([Path(build)/'check']).stdout)
            invalid = run([Path(build)/'check', 'bad'], success=False)
            self.assertNotEqual(invalid.returncode, 0)
            self.assertNotIn('invalid visible content was accepted', invalid.stdout)


if __name__ == '__main__':
    unittest.main()

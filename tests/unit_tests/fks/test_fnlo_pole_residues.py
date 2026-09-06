"""The fNLO SDE partition needs resonance numerators, not regulated poles."""

import math
import sys
from types import SimpleNamespace
from unittest import mock

import tests.unit_tests as unittest
import aloha
from aloha import create_aloha
from aloha.template_files import wavefunctions


class TestFNLOPoleResidues(unittest.TestCase):
    def test_channel_ratios_and_on_shell_finiteness(self):
        # The special tag is used only by tagged, flattened fNLO connectors.
        # Include both fermion flows and polarized massive-vector projectors.
        scalar = SimpleNamespace(name='SSS1', spins=[1, 1, 1], structure='1')
        vector = SimpleNamespace(name='FFV1', spins=[2, 2, 3],
                                 structure='Gamma(3,2,1)')
        cases = [(scalar, 3, []), (vector, 1, []), (vector, 2, []),
                 (vector, 3, []), (vector, 3, ['P1L']),
                 (vector, 3, ['P1T'])]
        with mock.patch.object(aloha, 'loop_mode', False), \
                mock.patch.object(aloha, 'unitary_gauge', True), \
                mock.patch.object(aloha, 'complex_mass', False), \
                mock.patch.dict(sys.modules, wavefunctions=wavefunctions):
            for lorentz, outgoing, tags in cases:
                with self.subTest(spins=lorentz.spins, outgoing=outgoing, tags=tags):
                    builder = create_aloha.AbstractRoutineBuilder(lorentz)
                    ordinary = builder.compute_routine(outgoing, tags)
                    residue = builder.compute_routine(outgoing, ['NWA'] + tags)
                    self.assertEqual(residue.denominator, 1)
                    # Numerators, including mass/projector terms, are unchanged.
                    self.assertEqual(str(ordinary.expr), str(residue.expr))
                    namespace = {}
                    exec(ordinary.write(None, 'Python'), namespace)
                    exec(residue.write(None, 'Python'), namespace)
                    full = namespace[lorentz.name + ''.join(tags) + '_%d' % outgoing]
                    stripped = namespace[lorentz.name + 'NWA' + ''.join(tags) + '_%d' % outgoing]
                    ordinary_weights, residue_weights = [], []
                    for sample in range(1, 4):
                        first = wavefunctions.WaveFunction(size=6)
                        second = wavefunctions.WaveFunction(size=6)
                        first[:] = [6.+1.j, 1.+.5j] + [complex(sample+i, i) for i in range(4)]
                        second[:] = [5.+2.j, 3.-1.5j] + [complex(i-sample, sample*i) for i in range(4)]
                        full_value = full(first, second, 1., 3., 0.)
                        residue_value = stripped(first, second, 1., 3., 0.)
                        ordinary_weights.append(sum(abs(v)**2 for v in full_value[2:]))
                        residue_weights.append(sum(abs(v)**2 for v in residue_value[2:]))
                        # E=11, px=4, py=-1, pz=3: the parent is exactly on shell.
                        on_shell = stripped(first, second, 1., math.sqrt(95.), 0.)
                        self.assertTrue(all(math.isfinite(abs(v)) for v in on_shell))
                        self.assertGreater(sum(abs(v)**2 for v in on_shell[2:]), 0.)
                    for full_weight, stripped_weight in zip(ordinary_weights, residue_weights):
                        self.assertAlmostEqual(full_weight/sum(ordinary_weights),
                                               stripped_weight/sum(residue_weights), places=13)

import math
import unittest

from bw_branching import GF, MW, MT, convolution_density, current_density, partial_width


class TestBWBranching(unittest.TestCase):
    def test_fixed_coupling_current_identity_massless_and_massive(self):
        for mb in (0.,4.8):
            for gw in (2.,.1):
                b = GF*MW**3/(6*math.pi*math.sqrt(2)*gw)
                for s in (1.e-5,MW**2/3,MW**2,.99*(MT-mb)**2):
                    self.assertAlmostEqual(current_density(s,MT,MW,mb,gw)/
                                           convolution_density(s,MT,MW,mb,gw),b,places=13)

    def test_semileptonic_narrow_width_limit_has_inverse_width_scaling(self):
        gw = 1.e-3
        value,error = partial_width(MT,MW,0.,gw)
        b = GF*MW**3/(6*math.pi*math.sqrt(2)*gw)
        x = (MW/MT)**2
        g0 = GF*MT**3/(8*math.pi*math.sqrt(2))*(1-x)**2*(1+2*x)
        self.assertLess(abs(value/b/g0-1.),2.e-5)
        self.assertLess(error/value,1.e-8)

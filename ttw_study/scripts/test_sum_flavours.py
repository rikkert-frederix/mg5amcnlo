import copy
import unittest

import numpy as np

from sum_flavours import FLAVOURS, combine
import test_load_results


class FullFlavourChecks(unittest.TestCase):
    def samples(self):
        prototype=test_load_results.TestMatchedIdentity().samples()['LO']
        results=[]
        for seed,flavours in enumerate(sorted(FLAVOURS),1):
            row=copy.deepcopy(prototype)
            row['report']['manifest']['settings']['iseed']=seed
            row['report']['execution']['export_generation']['args']['flavours']=list(flavours)
            results.append(row)
        return results

    def test_exact_disjoint_sum_and_no_input_mutation(self):
        samples=self.samples()
        before=copy.deepcopy(samples[0]['report'])
        result=combine(samples)
        np.testing.assert_allclose(result['values'],8.)
        np.testing.assert_allclose(result['errors'],np.sqrt(8.)*.1)
        self.assertEqual(samples[0]['report'],before)
        self.assertEqual(result['flavours'],sorted(FLAVOURS))

    def test_missing_repeated_and_mismatched_samples_rejected(self):
        with self.assertRaises(ValueError):
            combine(self.samples()[:-1])
        for change in ('flavour','seed','mass','charge','variant'):
            samples=self.samples()
            report=samples[-1]['report']
            generation=report['execution']['export_generation']['args']
            if change=='flavour':
                generation['flavours']=samples[0]['report']['execution']['export_generation']['args']['flavours']
            elif change=='seed':
                report['manifest']['settings']['iseed']=1
            elif change=='mass':
                report['manifest']['hashes']['param_card.dat']='changed mass'
            elif change=='charge':
                generation['charge']='minus'
            else:
                report['variant']='S'
            with self.assertRaises(ValueError):
                combine(samples)


if __name__ == '__main__':
    unittest.main()

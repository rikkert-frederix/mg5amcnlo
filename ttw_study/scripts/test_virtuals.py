import unittest

from audit_virtuals import check_log


class TestVirtualAudit(unittest.TestCase):
    def log(self):
        return '\n'.join([
            *('INFO: analytic top-decay virtual validation %d/3 for contribution %d; '
              'relative difference 1.2345E-13' % (i, c)
              for c in (2, 3) for i in (1, 2, 3)),
            *('INFO: analytic top-decay virtual provider for contribution %d validated;' % c
              for c in (2, 3)), 'Time spent in Total :'])

    def test_complete_both_branches(self):
        self.assertEqual(check_log(self.log()), {'2': 1.2345e-13, '3': 1.2345e-13})

    def test_reject_incomplete_duplicate_failed_or_nonfinite(self):
        for log in (self.log().replace('3/3', '2/3', 1),
                    self.log().replace('1.2345E-13', '1.1E-8', 1),
                    self.log().replace('1.2345E-13', 'NaN', 1),
                    self.log().replace('contribution 3 validated;', 'contribution 4 validated;'),
                    self.log().replace('Time spent in Total :', '')):
            with self.assertRaises(ValueError):
                check_log(log)

import hashlib
import unittest

from card_writer_transition import AFTER, BEFORE, old_source


class TestCardWriterTransition(unittest.TestCase):
    def test_exact_single_line_reversal_requires_frozen_hash(self):
        before='prefix\n'+BEFORE+'suffix\n'
        after='prefix\n'+AFTER+'suffix\n'
        checksum=hashlib.sha256(before.encode()).hexdigest()
        self.assertEqual(old_source(after,checksum),before)
        for changed in (after+'extra\n',after+AFTER,after.replace(AFTER,BEFORE)):
            with self.assertRaises(ValueError):
                old_source(changed,checksum)
        with self.assertRaises(ValueError):
            old_source(after,'0'*64)


if __name__=='__main__':
    unittest.main()

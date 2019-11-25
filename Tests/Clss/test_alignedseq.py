from Clss.Model.AlignedSeq import AlignedSeq
import unittest


class TestAlignedSeq(unittest.TestCase):
    def test_content_calc(self):
        seq = "--GA--ACTG--"
        start, end = AlignedSeq.content_calc(seq)
        self.assertEqual(2, start)
        self.assertEqual(9, end)

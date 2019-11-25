import unittest
from Clss.Model.Alignment import Alignment
from Clss.Model.AlignedSeq import AlignedSeq


class TestAligns(unittest.TestCase):
    def test_count(self):
        counts = {"A": [0, 1, 0], "T": [0, 2, 1], "G": [1, 0, 0], "C": [0, 0, 0]}
        Alignment.count("T", counts, 0)
        expected_counts = {"A": [0, 1, 0], "T": [1, 2, 1], "G": [1, 0, 0], "C": [0, 0, 0]}
        self.assertEqual(expected_counts, counts)

    def test_get_counts(self):
        # setup data
        seq1 = AlignedSeq("ACTCG")
        seq2 = AlignedSeq("AT-GG")
        seq3 = AlignedSeq("-CTGG")
        alignment = Alignment([seq1, seq2, seq3])

        # testing function
        counts = alignment.get_counts(full_length=True)

        expected_counts = {
            "A": [2, 0, 0, 0, 0], "C": [0, 2, 0, 1, 0], "G": [0, 0, 0, 2, 3], "T": [0, 1, 2, 0, 0], "-": [1, 0, 1, 0, 0]
        }
        # assertion
        self.assertEqual(expected_counts, counts)

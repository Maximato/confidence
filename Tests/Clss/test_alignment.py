import unittest
from Clss.Model.Alignment import Alignment
from Clss.Model.AlignedSeq import AlignedSeq


class TestAligns(unittest.TestCase):
    def test_count(self):
        counts = {"A": [0, 1, 0], "T": [0, 2, 1], "G": [1, 0, 0], "C": [0, 0, 0]}
        Alignment.count("T", counts, 0)
        expected_counts = {"A": [0, 1, 0], "T": [1, 2, 1], "G": [1, 0, 0], "C": [0, 0, 0]}
        self.assertEqual(expected_counts, counts)

    def test_get_counts_full_length(self):
        # setup data
        seq1 = AlignedSeq("ACTC")
        seq2 = AlignedSeq("AT-G")
        seq3 = AlignedSeq("-CTG")
        alignment = Alignment([seq1, seq2, seq3])

        # testing function
        counts = alignment.get_counts(full_length=True)

        # assertion
        expected_counts = {"A": [2, 0, 0, 0], "C": [0, 2, 0, 1], "G": [0, 0, 0, 2], "T": [0, 1, 2, 0], "-": [1, 0, 1, 0]}
        self.assertEqual(expected_counts, counts)

    def test_get_counts_part_length(self):
        # setup data
        seq1 = AlignedSeq("AC--")
        seq2 = AlignedSeq("AT-G")
        seq3 = AlignedSeq("-CTG")
        alignment = Alignment([seq1, seq2, seq3])

        # testing function
        counts = alignment.get_counts(full_length=False)

        # assertion
        expected_counts = {"A": [2, 0, 0, 0], "C": [0, 2, 0, 0], "G": [0, 0, 0, 2], "T": [0, 1, 1, 0], "-": [0, 0, 1, 0]}
        self.assertEqual(expected_counts, counts)

import unittest
from Clss.Model.Alignment import Alignment
from Clss.Model.AlignedSeq import AlignedSeq
from Clss.Model.Counts import Counts


class TestAlignmemt(unittest.TestCase):
    def setUp(self):
        # setup data
        seq1 = AlignedSeq("ACTC")
        seq2 = AlignedSeq("AT-G")
        seq3 = AlignedSeq("-CTG")
        self.alignment = Alignment([seq1, seq2, seq3])

    def test_get_counts_full_length(self):
        # testing function
        counts = self.alignment._Alignment__get_counts(full_length=True)

        # assertion
        expected_counts = {"A": [2, 0, 0, 0], "C": [0, 2, 0, 1], "G": [0, 0, 0, 2], "T": [0, 1, 2, 0], "-": [1, 0, 1, 0]}
        self.assertEqual(expected_counts, counts.counts)

    def test_get_counts_part_length(self):
        # testing function
        counts = self.alignment._Alignment__get_counts(full_length=False)

        # assertion
        expected_counts = {"A": [2, 0, 0, 0], "C": [0, 2, 0, 1], "G": [0, 0, 0, 2], "T": [0, 1, 2, 0], "-": [0, 0, 1, 0]}
        self.assertEqual(expected_counts, counts.counts)

    def test_get_consensus(self):
        # test function
        consensus = self.alignment.get_consensus()

        # assertion
        self.assertIsNotNone(consensus)

    def test_consensus_with(self):
        # setup data
        seqs = ["ATCGAC", "ACTAGC"]

        # test function
        consensus = self.alignment.consensus_with(seqs)

        # assertion
        self.assertIsNotNone(consensus)

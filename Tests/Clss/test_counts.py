from Clss.Model.Counts import Counts
from unittest import TestCase


class TestCounts(TestCase):
    def test_get_consensus(self):
        # setup data
        counts = Counts({"A": [2, 0, 0], "C": [0, 2, 0], "G": [0, 0, 0], "T": [0, 1, 2], "-": [1, 0, 1]})

        # test function
        consensus = counts.get_consensus()
        print(consensus)

        # assertion
        self.assertIsNotNone(consensus)

    def test_recount_with(self):
        # setup data
        counts = Counts({"A": [2, 0, 0], "C": [0, 2, 0], "G": [0, 0, 0], "T": [0, 1, 2], "-": [1, 0, 1]})

        # test function
        recount = counts.recount_with(["AT"])
        print(recount)
        # assertion
        self.assertIsNotNone(recount)

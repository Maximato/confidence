from Clss.Model.Counts import Counts
from unittest import TestCase


class TestCounts(TestCase):
    def test_count(self):
        # setup data
        counts = Counts(3)

        # testing function
        counts.count("A", 1)
        counts.count("-", 0)

        # assertion
        expected_counts = {"A": [0, 1, 0], "T": [0, 0, 0], "G": [0, 0, 0], "C": [0, 0, 0], "-": [1, 0, 0]}
        self.assertEqual(expected_counts, counts.counts)

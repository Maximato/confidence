from Clss.Model.Clusterization import Clusterization
from Bio.SeqRecord import SeqRecord
from unittest import TestCase


class TestClusterisation(TestCase):

    def setUp(self):
        self.seq1 = SeqRecord("CGGG", id="o1_1", description="organism1")
        self.seq2 = SeqRecord("ACGA", id="o1_2", description="organism1")
        self.seq3 = SeqRecord("ATGATCGA", id="o3_1", description="organism2")
        self.sequences = Clusterization([self.seq1, self.seq2, self.seq3])

    def test_create_dist_matrix(self):
        # test function
        dist_matrix = self.sequences.create_dist_matrix()
        print(dist_matrix)

        # assertion
        self.assertEqual(3, len(dist_matrix))

    def test_clusterize(self):
        # test function
        db = self.sequences.clusterize()

        # assertion
        self.assertIsNotNone(db)

    def test_get_clusters(self):
        # test function
        self.sequences.clusterize()
        clusters = self.sequences.get_clusters()

        # assertion
        self.assertIsNotNone(clusters)

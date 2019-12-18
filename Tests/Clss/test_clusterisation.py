from Clss.Model.Clusterization import Clusterization
from unittest import TestCase


class TestClusterisation(TestCase):

    def setUp(self):
        dm = [[0, 2], [2, 0]]
        self.cl = Clusterization(dm)

    def test_get_clusters_indexes(self):
        # test function
        indexes = self.cl.get_clusters_indexes()

        # assertion
        self.assertIsNotNone(indexes)

    def test_get_coordinates(self):
        # test function
        coordinates = self.cl.get_coordinates()

        # assertion
        self.assertIsNotNone(coordinates)

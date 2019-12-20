from Clss.Controller.ClusterizationController import ClusterisationController
from definitions import *
from os.path import join
import pytest

cc = ClusterisationController()
test_data = join(TEST_IN_DATA, "seqs.fasta")


@pytest.mark.skip
def test_clusterisation_without_dm():
    cc.clusterization(test_data, join(TEST_OUT_DATA, "clusters"), 0.4, 2)


def test_clusterisation_with_gadm():
    dm_path = join(TEST_IN_DATA, "seqs_ga.dm")
    cc.clusterization(test_data, join(TEST_OUT_DATA, "clusters_ga"), 0.4, 2, dm_path)


def test_clusterisation_with_clustdm():
    dm_path = join(TEST_IN_DATA, "seqs_clustalo.dm")
    cc.clusterization(test_data, join(TEST_OUT_DATA, "clusters_clustalo"), 0.2, 2, dm_path)

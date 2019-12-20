from Clss.Controller.AlignController import AlignController
from definitions import *
from os.path import join


ac = AlignController()
test_data = join(TEST_IN_DATA, "align.fasta")


def test_convert_to():
    ac.convert_to(test_data, join(TEST_OUT_DATA, "consensus.fasta"), False, "c", "fasta")


def test_convert_in_all_combinations():
    ac.convert_in_all_combinations(test_data, join(TEST_OUT_DATA, "cns"), "consensus")


def test_consensus_with():
    seqs_file = join(TEST_DATA_PATH, "in", "short_seqs.fasta")
    ac.consensus_with(test_data, seqs_file, join(TEST_OUT_DATA, "new_consensus.html"))


def test_unite_aligns():
    data_dir = join(TEST_IN_DATA, "aligns")
    ac.unite_aligns(data_dir, join(TEST_OUT_DATA, "aligns_consensus"), False)

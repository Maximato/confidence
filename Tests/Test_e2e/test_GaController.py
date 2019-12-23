from Clss.Controller.RunAlignController import GaController
from definitions import *
from os.path import join

ga = GaController()


def test_align():
    test_data = join(TEST_IN_DATA, "seqs.fasta")
    ga.align(test_data, join(TEST_OUT_DATA, "aligned_seq.fasta"))


def test_align_groups():
    test_dir = join(TEST_IN_DATA, "groups")
    ga.align_groups(test_dir, join(TEST_OUT_DATA, "aligned_groups"))

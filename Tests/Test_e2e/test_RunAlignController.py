from Clss.Controller.RunAlignController import RunAlignController
from definitions import *
from os.path import join

rac = RunAlignController()


def test_align():
    test_data = join(TEST_IN_DATA, "seqs.fasta")
    rac.align(test_data, "muscle", join(TEST_OUT_DATA, "aligned_seq.fasta"))


def test_align_groups():
    test_dir = join(TEST_IN_DATA, "groups")
    rac.align_groups(test_dir, "muscle", join(TEST_OUT_DATA, "aligned_groups"))

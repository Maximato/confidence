from Clss.Controller.GaController import GaController
from definitions import *
from os.path import join


ga = GaController()
test_data = join(TEST_DATA_PATH, "seqs.fasta")
ga.align(test_data, join(TEST_DATA_PATH, "aligned_seq.fasta"))

test_dir = join(TEST_DATA_PATH, "groups")
ga.align_groups(test_dir, join(TEST_DATA_PATH, "aligned_groups"))

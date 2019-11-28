from Clss.Controller.AlignController import AlignController
from definitions import *
from os.path import join


ac = AlignController()
test_data = join(TEST_DATA_PATH, "align.fasta")
ac.convert_to(test_data, join(TEST_DATA_PATH, "consensus.fasta"), False, "c", "html")
ac.convert_in_all_combinations(test_data, join(TEST_DATA_PATH, "cns"), "consensus")

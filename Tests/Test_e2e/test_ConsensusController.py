from Clss.Controller.HtmlConsensusController import HtmlConsensusController
from definitions import *
from os.path import join


cc = HtmlConsensusController()
test_data = join(TEST_IN_DATA, "consensus.html")


def test_convert_to_mutations():
    cc.convert_to_mutations(test_data, join(TEST_OUT_DATA, "cwith_mut.fasta"), ["c90", "c80"])

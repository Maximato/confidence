from Clss.Controller.RecordsController import RecordsController
from definitions import *
from os.path import join


rc = RecordsController()
test_data = join(TEST_IN_DATA, "seqs.fasta")


def test_grouping():
    out_dir = join(TEST_OUT_DATA, "groups")
    rc.grouping(test_data, out_dir, 300, 1000)


def test_filtrating():
    organism = "tick-borne encephalitis virus"
    rc.filtrating(test_data, join(TEST_OUT_DATA, "filtr.fasta"), organism, 250, 300)

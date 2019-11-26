from Clss.Controller.RecordsController import RecordsController
from definitions import *
from os.path import join


rc = RecordsController()
test_data = join(TEST_DATA_PATH, "seqs.fasta")
out_dir = join(TEST_DATA_PATH, "outdir")
rc.grouping(test_data, 500, out_dir)

organism = "tick-borne encephalitis virus"
rc.filtrating(test_data, join(TEST_DATA_PATH, "filtr.fasta"), organism, 250, 300)

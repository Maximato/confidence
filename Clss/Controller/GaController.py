from Clss.Model.Alignment import Alignment
from Bio import SeqIO
from Clss.Model.Consensus import Consensus
from Clss.Model.Records import Records
from Clss.GramAlign import GramAlign
import os
from os.path import join, basename, dirname


class GaController:
    @staticmethod
    def align_groups(groups_dir, align_dir):
        # running GramAlign for all groped sequences
        ga = GramAlign()
        # ga.run_gram_align("file.fasta")
        ga.run_for_all_in(groups_dir, align_dir)

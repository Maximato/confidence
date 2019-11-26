from Clss.Model.Alignment import Alignment
from Clss.Model.Consensus import Consensus
from Clss.Model.Records import Records
from Clss.Controller.ConsensusController import ConsensusController
from Clss.GramAlign import GramAlign
import os
from os.path import join, basename, dirname


class AlignController:

    def __init__(self, alignment):
        self.alignment = alignment

    def write_html_to(self, filename, full_length, coloring):
        counts = self.alignment.get_counts(full_length)
        consensus = Consensus.get_consensus(counts)
        cc = ConsensusController(consensus)
        cc.write_html_to(filename, coloring)

    def write_seqrec_to(self, filename, full_length):
        counts = self.alignment.get_counts(full_length)
        consensus = Consensus.get_consensus(counts)
        cc = ConsensusController(consensus)
        cc.write_seqrec_to(filename)

    def write_all_to(self, prefix):
        self.write_html_to(prefix + "TC.html", True, "c")
        self.write_html_to(prefix + "TD.html", True, "d")
        self.write_html_to(prefix + "FC.html", False, "c")
        self.write_html_to(prefix + "FD.html", False, "d")

        self.write_seqrec_to(prefix + "T.txt", True)
        self.write_seqrec_to(prefix + "F.txt", False)

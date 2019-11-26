from Bio import SeqIO
from Clss.Model.Alignment import Alignment
from Clss.Model.AlignedSeq import AlignedSeq
from os.path import join, abspath, basename, dirname
from definitions import *


CONF_HTML_HEADER = join(PROJECT_PATH, "Clss", "Extractor", "header_for_confidence.txt")
DEEP_HTML_HEADER = join(PROJECT_PATH, "Clss", "Extractor", "header_for_deeps.txt")


class Extractor:
    @staticmethod
    def get_html_header(coloring):
        if coloring == "c":
            filename = CONF_HTML_HEADER
        elif coloring == "d":
            filename = DEEP_HTML_HEADER
        else:
            raise ValueError("Argument should be 'c' or 'd'")

        with open(filename, "r") as f:
            header = f.read()
        return header

    @staticmethod
    def recs_extractor(filename):
        return list(SeqIO.parse(filename, "fasta"))

    @staticmethod
    def align_extractor(filename):
        recs = Extractor.recs_extractor(filename)
        aligned_seqs = [AlignedSeq(rec.seq) for rec in recs]
        align = Alignment(aligned_seqs)
        return align

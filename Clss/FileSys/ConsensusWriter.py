from Clss.Extractor.Extractor import Extractor
from Clss.Model.Consensus import Consensus
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from os.path import join, basename, abspath, dirname, isdir
import os


class ConsensusWriter:
    def __init__(self, consensus):
        self.consensus = consensus

    def write_html_to(self, filename, coloring="c"):
        html_header = Extractor.get_html_header(coloring)
        html_body = Consensus.get_html_body(self.consensus, coloring)

        dname = dirname(filename)
        if not isdir(dname):
            os.mkdir(dname)

        with open(filename, "w") as f:
            f.write(html_header + html_body)

    def write_seqrec_to(self, filename):
        seq = Consensus.get_seq_consensus(self.consensus)
        seq_rec = SeqRecord(seq)
        print(seq_rec)
        SeqIO.write(seq_rec, filename, "fasta")

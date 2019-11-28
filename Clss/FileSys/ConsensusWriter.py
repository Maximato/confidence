from Clss.Extractor.Extractor import Extractor
from Clss.Model.Consensus import Consensus
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from os.path import join, basename, abspath, dirname, isdir
import os


class ConsensusWriter:
    def __init__(self, consensus):
        self.consensus = consensus

    def __write_html_to(self, filename, coloring="c"):
        html_header = Extractor.get_html_header(coloring)
        html_body = Consensus.get_html_body(self.consensus, coloring)

        dname = dirname(filename)
        if not isdir(dname):
            os.mkdir(dname)

        with open(filename, "w") as f:
            f.write(html_header + html_body)

    def __write_seqrec_to(self, filename):
        seq = Consensus.get_seq_consensus(self.consensus)
        seq_rec = SeqRecord(seq)
        print(seq_rec)
        SeqIO.write(seq_rec, filename, "fasta")

    def write(self, filename, coloring="c", fmt="html"):
        if fmt not in ["html", "fasta"]:
            fmt = "html"
            print("Wrong format! Should be fasta or html.")

        if fmt == "html":
            self.__write_html_to(filename, coloring)
        elif fmt == "fasta":
            self.__write_seqrec_to(filename)

    def write_all(self, outdir, prefix):
        self.write(join(outdir, f"{prefix}_TC.html"), "c", "html")
        self.write(join(outdir, f"{prefix}_TC.html"), "d", "html")
        self.write(join(outdir, f"{prefix}_TC.html"), fmt="fasta")

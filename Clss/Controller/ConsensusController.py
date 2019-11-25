from Clss.htmlWriter import htmlWriter
from Clss.Extractor.Extractor import Extractor
from Clss.Model.Consensus import Consensus
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write


class ConsensusController:
    def __init__(self, consensus):
        self.consensus = consensus

    def write_html_to(self, filename, coloring="c"):
        html_header = Extractor.get_html_header(coloring)
        html_body = Consensus.get_html_body(self.consensus, coloring)

        html_writer = htmlWriter(html_header + html_body)
        html_writer.write_to(filename)

    def write_seqrec_to(self, filename):
        seq = Consensus.get_seq_consensus(self.consensus)
        seq_rec = SeqRecord(seq)
        write(seq_rec, filename, "fasta")

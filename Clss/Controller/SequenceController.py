from Clss.Model.Alignment import Alignment
from Bio import SeqIO
from Clss.Model.Consensus import Consensus
from Clss.Model.Sequences import Sequences
from Clss.GramAlign import GramAlign
import os
from os.path import join, basename, dirname


class SequencesController:

    def __init__(self, sequences):
        self.sequences = sequences

    def write_groups(self, esize, directory):
        groups = self.sequences.group(esize)
        if not os.path.isdir(directory):
            os.mkdir(directory)
        for key in groups:
            SeqIO.write(groups[key], join(directory, key + ".fasta"), "fasta")

    def write_filtrated(self, out_file, organism, minsize, maxsize):
        fseqs = self.sequences.filtr_organizm_by_size(organism, minsize, maxsize)
        SeqIO.write(fseqs, out_file, "fasta")

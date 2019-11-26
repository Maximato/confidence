from os.path import join, basename, abspath, dirname, isdir
from Bio import SeqIO
import os


class RecordsWriter:

    def __init__(self, records):
        self.records = records

    def write_to(self, filename):
        dname = dirname(filename)
        if not isdir(dname):
            os.mkdir(dname)
        SeqIO.write(self.records, filename, "fasta")

    def write_to_dir(self, basename, dirname):
        if not isdir(dirname):
            os.mkdir(dirname)
        SeqIO.write(self.records, join(dirname, basename), "fasta")

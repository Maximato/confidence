from Clss.FileSys.check_directory import check_dir
from os.path import join, isdir
from Bio import SeqIO
import os


class RecordsWriter:

    def __init__(self, records):
        self.records = records

    def write_to(self, filename):
        check_dir(filename)
        SeqIO.write(self.records, filename, "fasta")

    def write_to_dir(self, basename, dirname):
        if not isdir(dirname):
            os.mkdir(dirname)
        SeqIO.write(self.records, join(dirname, basename), "fasta")

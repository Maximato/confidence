from Bio import SeqIO
from definitions import *
from os.path import abspath, join
from os import listdir


CONF_HTML_HEADER = join(PROJECT_PATH, "Clss", "Extractor", "header_for_confidence.txt")
DEEP_HTML_HEADER = join(PROJECT_PATH, "Clss", "Extractor", "header_for_deeps.txt")


class Extractor:
    @staticmethod
    def extract_html_header(coloring):
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
    def extract_records(filename):
        return list(SeqIO.parse(filename, "fasta"))

    @staticmethod
    def extract_filenames(dirname):
        files = listdir(dirname)
        full_paths = [join(dirname, fn) for fn in files]
        return full_paths

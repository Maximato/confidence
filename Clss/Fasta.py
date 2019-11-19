from Bio import SeqIO


class Fasta:
    @staticmethod
    def extract_from(filename):
        return list(SeqIO.parse(filename, "fasta"))

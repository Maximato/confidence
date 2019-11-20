from Bio import SeqIO
from Clss.Fasta import Fasta
import os
from os.path import join


class Sequences(Fasta):

    @staticmethod
    def group_seqs(lseqs, estimated_size=10000):
        groups = {"cds": []}
        for record in lseqs:
            seq_id = record.id
            group_id = seq_id[0:5]
            length = len(record.seq)
            if length > estimated_size:
                groups["cds"].append(record)
            else:
                if group_id in groups.keys():
                    groups[group_id].append(record)
                else:
                    groups[group_id] = [record]
        return groups

    @staticmethod
    def write_groups(groups, directory):
        if not os.path.isdir(directory):
            os.mkdir(directory)
        for key in groups:
            SeqIO.write(groups[key], join(directory, key + ".fasta"), "fasta")

    @staticmethod
    def filtr_by_organizm(lseqs, organizm):
        fseqs = []
        for record in lseqs:
            if organizm in record.description:
                fseqs.append(record)
        return fseqs

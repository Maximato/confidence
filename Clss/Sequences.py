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
            group_id = seq_id[0:2]
            length = len(record.seq)
            if length > estimated_size:
                groups["cds"].append(record)
            else:
                if group_id in groups.keys():
                    groups[group_id].append(record)
                else:
                    groups[group_id] = [record]
        print("Number of groups: ", len(groups))
        return groups

    @staticmethod
    def write_groups(groups, directory):
        if not os.path.isdir(directory):
            os.mkdir(directory)
        for key in groups:
            SeqIO.write(groups[key], join(directory, key + ".fasta"), "fasta")

    @staticmethod
    def filtr_by(lseqs, organizm, minsize=100):
        print("Initially number of sequences: ", len(lseqs))
        fseqs = []
        for record in lseqs:
            description = record.description.lower()
            if organizm in description and len(record.seq) >= minsize:
                if ("chimeric" not in description) and ("chimera" not in description):
                    fseqs.append(record)
        print("Number of sequences after filtrating: ", len(fseqs))
        return fseqs

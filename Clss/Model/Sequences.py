from Bio import SeqIO
import os
from os.path import join


class Sequences:
    def __init__(self, seqs):
        self.seqs = seqs

    def group(self, estimated_size=10000):
        groups = {"cds": []}
        for record in self.seqs:
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

    def filtr_organizm_by_size(self, organizm, minsize, maxsize):
        print("Initially number of sequences: ", len(self.seqs))
        fseqs = []
        for record in self.seqs:
            description = record.description.lower()
            if organizm in description and (minsize <= len(record.seq) <= maxsize):
                if ("chimeric" not in description) and ("chimera" not in description):
                    fseqs.append(record)
        print("Number of sequences after filtrating: ", len(fseqs))
        return fseqs

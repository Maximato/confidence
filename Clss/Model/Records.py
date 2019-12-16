import numpy as np
from Bio import pairwise2
import time


class Records:
    def __init__(self, records):
        self.records = records
        self.size = len(records)

    def get_seqs(self):
        seqs = []
        for record in self.records:
            seqs.append(record.seq)
        return seqs

    def group(self, genome_min=9000, genome_max=12000):
        groups = {"cds": []}
        for record in self.records:
            rec_id = record.id
            group_id = rec_id[0:2]
            length = len(record.seq)
            if genome_min < length < genome_max:
                groups["cds"].append(record)
            else:
                if group_id in groups.keys():
                    groups[group_id].append(record)
                else:
                    groups[group_id] = [record]
        print("Number of groups: ", len(groups))
        return groups

    def filtr_organism_by_size(self, organizm, minsize, maxsize):
        print("Initially number of sequences: ", len(self.records))
        fseqs = []
        for record in self.records:
            description = record.description.lower()
            if organizm in description and (minsize <= len(record.seq) <= maxsize):
                if ("chimeric" not in description) and ("chimera" not in description):
                    fseqs.append(record)
        print("Number of sequences after filtrating: ", len(fseqs))
        return fseqs

    def create_dist_matrix(self):
        dm = np.zeros((self.size, self.size))
        for i in range(0, self.size):
            for j in range(i+1, self.size):
                align = pairwise2.align.globalxx(self.records[i].seq, self.records[j].seq, one_alignment_only=1)
                dm[i, j] = 1-align[0][2]/len(align[0][0])
                print(time.process_time())
                # raw.append(align[0][2])
        dm = dm + dm.T
        return dm

    def get_clusters(self, indexes):
        clusters = {}
        for ids in indexes:
            clusters[ids] = [self.records[i] for i in indexes[ids]]
        return clusters

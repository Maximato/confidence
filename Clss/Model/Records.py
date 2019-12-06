from sklearn.cluster import DBSCAN
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


class Records:
    def __init__(self, records):
        self.records = records

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
        dm = []
        for record1 in self.records:
            raw = []
            for record2 in self.records:
                align = pairwise2.align.globalxx(record1.seq, record2.seq, one_alignment_only=1)
                # print((align[0]))
                raw.append(align[0][2])
            dm.append(raw)
        return dm

    def pca(self, eps=0.825, ms=5, dist_matrix=None):
        if dist_matrix is None:
            dm = self.create_dist_matrix()
        else:
            dm = dist_matrix

        X = np.array(dm)
        db = DBSCAN(eps=eps, min_samples=ms).fit(X)
        print(db)
        print(db.labels_)
        return db

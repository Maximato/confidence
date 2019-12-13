import numpy as np
from Bio import pairwise2
from sklearn.cluster import DBSCAN
import time
from sklearn.preprocessing import StandardScaler


class Clusterization:
    def __init__(self, records, dm=None):
        super()
        self.records = records
        self.dm = dm
        # self.X = None
        self.db = None

    def create_dist_matrix(self):
        dm = []
        for record1 in self.records:
            raw = []
            for record2 in self.records:
                align = pairwise2.align.globalxx(record1.seq, record2.seq, one_alignment_only=1)
                raw.append(1-align[0][2]/len(align[0][0]))
                print(time.process_time())
                #raw.append(align[0][2])
            dm.append(raw)
        dm = np.array(dm)
        self.dm = dm
        return dm

    def clusterize(self, eps=1, ms=2):
        if self.dm is None:
            self.create_dist_matrix()
        db = DBSCAN(eps=eps, min_samples=ms, metric="precomputed").fit(self.dm)
        self.db = db
        return db

    def get_clusters(self):
        if self.db is None:
            raise ValueError("Db should not be None")

        labels = self.db.labels_
        clusters = {}
        for i, lb in enumerate(labels):
            if lb not in clusters.keys():
                clusters[lb] = [self.records[i]]
            else:
                clusters[lb].append(self.records[i])
        return clusters

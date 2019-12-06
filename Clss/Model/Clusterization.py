from sklearn.cluster import DBSCAN
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from sklearn.datasets import load_iris
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
import matplotlib.pyplot as pl


class Clusterization:
    def __init__(self, reords):
        super()
        self.records = reords
        self.dm = None
        self.db = None

    def create_dist_matrix(self):
        dm = []
        for record1 in self.records:
            raw = []
            for record2 in self.records:
                align = pairwise2.align.globalxx(record1.seq, record2.seq, one_alignment_only=1)
                # print((align[0]))
                raw.append(align[0][2])
            dm.append(raw)
        self.dm = dm
        return dm

    def clusterize(self, eps=1, ms=2):
        if self.dm is None:
            X = np.array(self.create_dist_matrix())
        else:
            X = np.array(self.dm)
        db = DBSCAN(eps=eps, min_samples=ms).fit(X)
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

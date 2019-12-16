from sklearn import manifold
from sklearn.cluster import DBSCAN


class Clusterization:
    def __init__(self, dm):
        super()
        self.dm = dm

    def get_clusters_indexes(self, eps=1, ms=2):
        db = DBSCAN(eps=eps, min_samples=ms, metric="precomputed").fit(self.dm)

        labels = db.labels_
        clusters_indexes = {}
        for i, lb in enumerate(labels):
            if lb not in clusters_indexes.keys():
                clusters_indexes[lb] = [i]
            else:
                clusters_indexes[lb].append(i)
        return clusters_indexes

    def get_coordinates(self):
        mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
        results = mds.fit(self.dm)
        return results.embedding_

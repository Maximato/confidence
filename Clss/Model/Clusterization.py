from sklearn import manifold
from sklearn.cluster import DBSCAN


class Clusterization:
    def __init__(self, dm):
        super()
        self.dm = dm

    def get_clusters_indexes(self, eps, ms):
        """
        Clusterisation based on distance matrix

        :param eps: float, the maximum distance between two samples
        :param ms: int, the number of samples (or total weight) in a neighborhood for a point to be considered
        as a core point. This includes the point itself
        :return: dict, keys - clusters, elements - indexes of cluster
        """
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
        """
        Calculation of coordinates of sequences (mark as blobs) on 2D expanse
        :return: list, coordinates
        """
        mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
        results = mds.fit(self.dm)
        return results.embedding_

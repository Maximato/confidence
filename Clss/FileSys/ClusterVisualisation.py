import numpy as np
import matplotlib.pyplot as plt
from os.path import join, isdir
import os


class ClusterVisualisation:
    def __init__(self, coordinates, dm, indexes):
        self.crds = coordinates
        self.dm = dm
        self.indexes = indexes

    def write_dm(self, outdir):
        if not isdir(outdir):
            os.mkdir(outdir)
        with open(join(outdir, "dist_matrix"), "w") as f:
            for x in self.dm:
                line = " ".join(map(str, x))
                f.write(line + "\n")

    def visualize(self, outdir):

        unique_labels = self.indexes.keys()
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]

        for cluster_index, color in zip(self.indexes, colors):
            cluster = self.crds[self.indexes[cluster_index]]
            ms = 14
            if cluster_index == -1:
                # Black used for noise.
                color = [0, 0, 0, 1]
                ms = 6
            plt.plot(cluster[:, 0], cluster[:, 1], 'o', markerfacecolor=tuple(color),
                     markeredgecolor='k', markersize=ms)

        n_clusters_ = len(unique_labels) - (1 if -1 in unique_labels else 0)
        plt.title('Estimated number of clusters: %d' % n_clusters_)
        plt.savefig(join(outdir, "clusters"))
        plt.show()

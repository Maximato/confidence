from sklearn.datasets import make_blobs
import numpy as np
from sklearn.preprocessing import StandardScaler
# X, y = make_blobs(n_samples=10, centers=3, n_features=3,random_state=0)


centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(n_samples=750, centers=centers, cluster_std=0.4,
                            random_state=0)
print(X)
X = StandardScaler().fit_transform(X)
print(X)

dist_m = [[1, 2, 2, 1], [2, 1, 2, 1], [3, 2, 4, 1], [1, 2, 3, 1]]
dist_m = np.array(dist_m)
print(dist_m)
dist_m = StandardScaler().fit_transform(dist_m)
print(dist_m)

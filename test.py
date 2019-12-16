from sklearn.datasets import make_blobs
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
import matplotlib.pyplot as plt
# X, y = make_blobs(n_samples=10, centers=3, n_features=3,random_state=0)


centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(n_samples=750, centers=centers, cluster_std=0.4,
                            random_state=0)
print(X)
X = StandardScaler().fit_transform(X)
print(X)

dist_m = [[1, 2, 2, 1], [2, 1, 2, 1], [3, 2, 4, 1], [43, 32, 1, 8]]
dist_m = np.array(dist_m)
print(dist_m)
#dist_m = StandardScaler().fit_transform(dist_m)
#print(dist_m)

labels = np.array([1, -1, 1, 2])
core_samples_mask = np.array([True, False, True, True])
class_member_mask = (labels == 2)
print("class member mask")
print(class_member_mask)

print("mask")
mask = (class_member_mask & core_samples_mask)
print(mask)

print("new")
xy = dist_m[mask]
print(xy)
print(xy[:, 0])
print(xy[:, 1])

adist = np.array([[0, 4, 23], [4, 0, 32], [23, 32, 0]])
amax = np.amax(adist)
print(amax)
adist = adist / amax
print(adist)

mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
results = mds.fit(adist)

coords = results.embedding_
print(coords)

plt.subplots_adjust(bottom=0.1)
plt.scatter(
    coords[:, 0], coords[:, 1], marker='o'
    )
plt.show()

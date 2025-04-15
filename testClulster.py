#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 09:30:01 2025

load in arrays and 
@author: allen
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn import metrics
from sklearn.cluster import AffinityPropagation
from sklearn.datasets import make_blobs

import pylab as plt

centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(
    n_samples=300, centers=centers, cluster_std=0.5, random_state=0
)

affmat =   X[ :, np.newaxis, : ] - X[ np.newaxis, :, : ]
affmatR = -np.square(np.abs(affmat).sum(axis=2))
affmatR2 = -(affmat*affmat).sum(axis=2)

plt.plot( X[:,0],X[:,1] , '.')

af = AffinityPropagation(preference=-50, random_state=0, affinity='precomputed').fit(affmatR)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_

for i in cluster_centers_indices:
    plt.plot(X[i,0],X[i,1],'+k')

plt.scatter( X[:,0],X[:,1] ,  c=labels)

n_clusters_ = len(cluster_centers_indices)

print("Estimated number of clusters: %d" % n_clusters_)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
print("Adjusted Rand Index: %0.3f" % metrics.adjusted_rand_score(labels_true, labels))
print(
    "Adjusted Mutual Information: %0.3f"
    % metrics.adjusted_mutual_info_score(labels_true, labels)
)
print(
    "Silhouette Coefficient: %0.3f"
    % metrics.silhouette_score(X, labels, metric="sqeuclidean")
)
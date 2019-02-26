#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
from sklearn.mixture import GaussianMixture

def gmm(encodings, **kw):
    nclusters = int(kw['nclusters']) if kw['nclusters'] != None else 3
    if kw['sof'] == 'sample':
        encodings = np.array(encodings)[1:]
    else:
        encodings = np.array(encodings).T[1:]

    if len(encodings) < nclusters:
        return 0, 'sample number should be greater than n_clusters.'

    data = encodings[:, 1:]
    shape = data.shape
    data = np.reshape(data, shape[0] * shape[1])
    data = np.reshape([float(i) for i in data], shape)
    e = ''
    try:
        gmm = GaussianMixture(n_components = nclusters).fit(data)
        cluster_pred = gmm.predict(data)
    except ValueError as e:
        print(e)
        return 0, e

    res = []
    for i in range(len(data)):
        res.append([encodings[i][0], cluster_pred[i]])
    return res, e

#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
from sklearn.decomposition import LatentDirichletAllocation


def lda(encodings, labels, n_components = 2):
    encodings = np.array(encodings)[1:]
    data = encodings[:, 1:]
    shape = data.shape
    data = np.reshape(data, shape[0] * shape[1])
    data = np.reshape([float(i) for i in data], shape)
    newData = LatentDirichletAllocation(n_components=n_components, max_iter=5,
                                learning_method='batch',
                                learning_offset=50.,
                                random_state=0).fit_transform(data, labels)
    lda = []
    for i in range(len(data)):
        lda.append([encodings[i][0]] + list(newData[i]))
    return lda
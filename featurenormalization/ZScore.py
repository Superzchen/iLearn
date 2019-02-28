#!/usr/bin/env python
# _*_coding:utf-8_*_

import numpy as np

def ZScore(encodings, labels):
    normalized_vector = np.zeros((len(encodings), len(encodings[0]) + 1)).astype(str)
    normalized_vector[0, 2:] = encodings[0][1:]
    normalized_vector[:, 0] = encodings[:, 0]
    normalized_vector[:, 1] = ['label'] + labels

    data = np.array(encodings[1:, 1:]).astype(float)
    std_array = np.std(data, axis=0)
    mean_array = np.mean(data, axis=0)

    e = ''
    try:
        for i in range(len(mean_array)):
            data[:, i] = (data[:, i] - mean_array[i]) / std_array[i]
    except ZeroDivisionError as e:
        return 0, e

    normalized_vector[1:, 2:] = data.astype(str)

    return normalized_vector.tolist(), e

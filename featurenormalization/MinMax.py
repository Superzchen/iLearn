#!/usr/bin/env python
# _*_coding:utf-8_*_

import numpy as np


def MinMax(encodings, labels):
    normalized_vector = np.zeros((len(encodings), len(encodings[0]) + 1)).astype(str)
    normalized_vector[0, 2:] = encodings[0][1:]
    normalized_vector[:, 0] = encodings[:, 0]
    normalized_vector[:, 1] = ['label'] + labels

    data = np.array(encodings[1:, 1:]).astype(float)

    e = ''
    for i in range(len(data[0])):
        maxValue, minValue = max(data[:, i]), min(data[:, i])
        try:
            data[:, i] = (data[:, i] - minValue) / maxValue
        except ZeroDivisionError as e:
            return 0, e

    normalized_vector[1:, 2:] = data.astype(str)

    return normalized_vector.tolist(), e

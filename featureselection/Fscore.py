#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
import math

def Calculate_Fscore(array, labels):
    if len(array) != len(labels):
        print('Error. inconsistent data shape with sample number')
        return 0

    array_po = []
    array_ne = []
    for i in range(len(labels)):
        if labels[i] == 1:
            array_po.append(array[i])
        else:
            array_ne.append(array[i])

    mean_po = sum(array_po) / len(array_po)
    mean_ne = sum(array_ne) / len(array_ne)
    mean = sum(array) / len(array)


    score_1 = ((mean_po - mean) ** 2 + (mean_ne - mean) ** 2)
    score_2 = sum([(i-mean_po) ** 2 for i in array_po]) / (len(array_po) - 1)
    score_3 = sum([(i-mean_ne) ** 2 for i in array_ne]) / (len(array_ne) - 1)
    f_score = score_1 / (score_2 + score_3)
    return f_score

def Fscore(encodings, labels):
    features = encodings[0][1:]
    encodings = np.array(encodings)[1:]
    data = encodings[:, 1:].astype(float)
    shape = data.shape

    e = ''
    if shape[0] < 5 or shape[1] < 2:
        return 0, e

    dataShape = data.shape

    if dataShape[1] != len(features):
        print('Error: inconsistent data shape with feature number.')
        return 0, 'Error: inconsistent data shape with feature number.'
    if dataShape[0] != len(labels):
        print('Error: inconsistent data shape with sample number.')
        return 0, 'Error: inconsistent data shape with sample number.'

    myFea = {}
    for i in range(len(features)):
        array = list(data[:, i])
        try:
            myFea[features[i]] = Calculate_Fscore(array, labels)
        except (ValueError, RuntimeWarning) as e:
            return 0, e

    res = []
    res.append(['feature', 'F-score'])
    for key in sorted(myFea.items(), key=lambda item: item[1], reverse=True):
        res.append([key[0], '{0:.3f}'.format(myFea[key[0]])])
    return res, e

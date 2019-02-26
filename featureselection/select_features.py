#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import re
import numpy as np

def select_features(encodings, labels, selected_features, selected_feature_number):
    rows = []
    for i in range(1, len(selected_features)):
        searchObj = re.search('f\.(\d+)', selected_features[i][0])
        if searchObj:
            rows.append(int(searchObj.group(1)))

    feature_number = selected_feature_number
    if len(rows) < selected_feature_number:
        feature_number = len(rows)
    rows = [0] + rows
    selected_feature_vectors = np.zeros((len(encodings), feature_number + 2)).astype(str)
    encodings = np.array(encodings)
    selected_feature_vectors[:, 0] = encodings[:, 0]
    labels = np.array(['label'] + labels)
    selected_feature_vectors[:, 1] =labels

    for i in range(1, feature_number + 1):
        selected_feature_vectors[:, i+1] = encodings[:, rows[i]]

    return selected_feature_vectors


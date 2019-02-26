#!/usr/bin/env python
#_*_coding:utf-8_*_

import argparse
from featureselection import *
from pubscripts import read_code, save_file
import sys, os, re
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="feature selection")
    parser.add_argument("--file", required=True, help="input encoding file format")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm',
                        help="the encoding type")
    parser.add_argument("--method", required=True,
                        choices=['CHI2', 'IG', 'MIC', 'pearsonr', 'Fscore'], help="select feature selection method")
    parser.add_argument("--num", type=int, default=100, help="the number of selected features, default is 100")
    parser.add_argument("--out", default="selected_features.txt", help="file with selected features")
    parser.add_argument("--rank", help="feature rank file")
    args = parser.parse_args()
    encodings, labels = read_code.read_code(args.file, format=args.format)
    cmd = args.method + '.' + args.method + '(encodings, labels)'
    print('Feature selection method: ' + args.method)
    selected_features, e = eval(cmd)
    featureRank = args.rank if args.rank != None else 'featureRank.txt'
    save_file.save_FS_result(selected_features, e, args.method, featureRank)

    # rows = []
    # for i in range(1, len(selected_features)):
    #     searchObj = re.search('f\.(\d+)', selected_features[i][0])
    #     if searchObj:
    #         rows.append(int(searchObj.group(1)))
    #
    # if len(rows) < args.num:
    #     args.num = len(rows)
    # rows = [0] + rows
    # selected_feature_vectors = np.zeros((len(encodings), args.num + 2)).astype(str)
    # encodings = np.array(encodings)
    # selected_feature_vectors[:, 0] = encodings[:, 0]
    # labels = np.array(['label'] + labels)
    # selected_feature_vectors[:, 1] =labels
    #
    # for i in range(1, args.num + 1):
    #     selected_feature_vectors[:, i+1] = encodings[:, rows[i]]

    selected_feature_vectors = select_features.select_features(encodings, labels, selected_features, args.num)
    save_file.save_file(selected_feature_vectors.tolist(), args.format, args.out)
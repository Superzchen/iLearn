#!/usr/bin/env python
#_*_coding:utf-8_*_

import argparse
from featurenormalization import *
from pubscripts import read_code, save_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="feature vector normalization")
    parser.add_argument("--file", required=True, help="input encoding file format")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm',
                        help="the encoding type")
    parser.add_argument("--method", required=True,
                        choices=['ZScore', 'MinMax'], help="select feature normalization method")
    parser.add_argument("--out", default="normalized_features.txt", help="file with normalized features vectors")

    args = parser.parse_args()
    encodings, labels = read_code.read_code(args.file, format=args.format)
    cmd = args.method + '.' + args.method + '(encodings, labels)'
    print('Feature normalization method: ' + args.method)
    normalized_feature_vectors, e = eval(cmd)
    save_file.save_file(normalized_feature_vectors, args.format, args.out)


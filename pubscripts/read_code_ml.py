#!/usr/bin/env python
# _*_coding:utf-8_*_

import os, sys, re
import numpy as np

def read_svm(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.readlines()

    for line in records:
        line = re.sub('\d+:', '', line)
        array = line.strip().split() if line.strip() != '' else None
        encodings.append(array[1:])
        labels.append(int(array[0]))

    return np.array(encodings).astype(float), np.array(labels).astype(int)

def read_tsv(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.readlines()

    for line in records:
        array = line.strip().split('\t') if line.strip() != '' else None
        encodings.append(array[1:])
        labels.append(int(array[0]))

    return np.array(encodings).astype(float), np.array(labels).astype(int)

def read_csv(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.readlines()

    sample = 1
    for line in records:
        array = line.strip().split(',') if line.strip() != '' else None
        encodings.append(array[1:])
        labels.append(int(array[0]))

    return np.array(encodings).astype(float), np.array(labels).astype(int)

def read_weka(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.read().strip().split('data\n')

    labels = []
    tmp = records[1].strip().split('\n')
    for i in tmp:
        tmp_arr = i.strip().split(',') if i.strip() != '' else None
        encodings.append(tmp_arr[0:-1])
        if tmp_arr[-1] == 'yes':
            labels.append(1)
        else:
            labels.append(0)

    return np.array(encodings).astype(float), np.array(labels).astype(int)

def read_code(file, format='svm'):
    encodings = []
    labels = []
    if not os.path.exists(file):
        print('Error: file does not exist.')
        sys.exit(1)
    if format == 'svm':
        encodings, labels = read_svm(file)
    if format == 'tsv':
        encodings, labels = read_tsv(file)
    if format == 'csv':
        encodings, labels = read_csv(file)
    if format == 'weka':
        encodings, labels = read_weka(file)
    return encodings, labels

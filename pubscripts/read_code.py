#!/usr/bin/env python
# _*_coding:utf-8_*_

import os, sys, re
import numpy as np


def read_tsv_1(file):
    encodings = []
    if not os.path.exists(file):
        print('Error: file does not exist.')
        sys.exit(1)
    with open(file) as f:
        records = f.readlines()
    for i in records:
        array = i.rstrip().split('\t') if i.strip() != '' else None
        encodings.append(array)
    return np.array(encodings)


def read_svm(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.readlines()

    ##
    feature = 1
    header = ['#']
    for i in range(1, len(records[0].split())):
        header.append('f.%d' %feature)
        feature = feature + 1
    encodings.append(header)

    ##
    sample = 1
    for line in records:
        line = re.sub('\d+:', '', line)
        array = line.strip().split() if line.strip() != '' else None
        encodings.append(['s.%d' % sample] + array[1:])
        labels.append(int(array[0]))
        sample = sample + 1
    return np.array(encodings), labels

def read_tsv(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.readlines()

    ##
    feature = 1
    header = ['#']
    for i in range(1, len(records[0].split())):
        header.append('f.%d' % feature)
        feature = feature + 1
    encodings.append(header)

    ##
    sample = 1
    for line in records:
        array = line.strip().split('\t') if line.strip() != '' else None
        encodings.append(['s.%d' % sample] + array[1:])
        labels.append(int(array[0]))
        sample = sample + 1
    return np.array(encodings), labels

def read_csv(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.readlines()

    ##
    feature = 1
    header = ['#']
    for i in range(1, len(records[0].split(','))):
        header.append('f.%d' % feature)
        feature = feature + 1
    encodings.append(header)

    ##
    sample = 1
    for line in records:
        array = line.strip().split(',') if line.strip() != '' else None
        encodings.append(['s.%d' % sample] + array[1:])
        labels.append(int(array[0]))
        sample = sample + 1
    return np.array(encodings), labels

def read_weka(file):
    encodings = []
    labels = []
    with open(file) as f:
        records = f.read()
    array = records.split('@data\n')
    tmp = array[0].split('\n')
    header = ['#']
    for i in tmp:
        tmp_arr = i.split()
        if len(tmp_arr) == 3 and tmp_arr[1] != 'play':
            header.append(tmp_arr[1])
    encodings.append(header)

    labels = []
    tmp = array[1].strip().split('\n')
    sample = 1
    for i in tmp:
        tmp_arr = i.strip().split(',') if i.strip() != '' else None
        encodings.append(['s.%d' %sample] + tmp_arr[0:-1])
        if tmp_arr[-1] == 'yes':
            labels.append(1)
        else:
            labels.append(0)
        sample = sample + 1

    return np.array(encodings), labels

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

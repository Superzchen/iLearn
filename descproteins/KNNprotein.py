#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import math
import numpy as np
import sys, os, platform
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import check_sequences
import needleman_wunsch

def CalculateSimilarity(sequence1, sequence2):
    blosumFile = re.sub('descproteins$', '', os.path.split(os.path.realpath(__file__))[
        0]) + r'\data\blosum62.txt' if platform.system() == 'Windows' else sys.path[0] + '/data/blosum62.txt'
    gap = [-10, -1]
    f = open(blosumFile)
    raw_matrix = [line.split() for line in f]
    f.close()
    raw_dicts = [dict() for x in range(len(raw_matrix[0]))]
    for i in range(len(raw_matrix[0])):
        raw_dicts[i] = dict(zip(raw_matrix[0], map(int, raw_matrix[i + 1])))

    s_matrix = dict()
    for i in range(len(raw_matrix[0])):
        s_matrix[raw_matrix[0][i]] = raw_dicts[i]
    #alignment = pairwise2.align.globalds(sequence1, sequence2, blosum62, -10, -0.5)[0]
    alignment = needleman_wunsch.matrix_filling_NW([sequence1, sequence2], s_matrix, gap)
    sum = 0
    for i in range(len(alignment[0])):
        if alignment[0][i] == alignment[1][i]:
            sum = sum + 1
    return 2 * sum / (len(sequence1) + len(sequence2))

def CalculateContent(mySimilarity, j, myLabelSets):
    content = []
    myDict = {}
    for i in myLabelSets:
        myDict[i] = 0
    for i in range(j):
        myDict[mySimilarity[i][0]] = myDict[mySimilarity[i][0]] + 1
    for i in myLabelSets:
        content.append(myDict[myLabelSets[i]] / j)
    return content

def KNNprotein(fastas, **kw):
    trainData = []
    myLabel = {}

    for i in fastas:
        if i[3] == 'training':
            trainData.append(i)
            myLabel[i[0]] = int(i[2])
    myLabelSets = list(set(myLabel.values()))

    kValues = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15,
               0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30]
    kNum = []
    for i in kValues:
        kNum.append(math.ceil(len(trainData) * i))

    encodings = []
    header = ['#', 'label']
    for k in kValues:
        for l in myLabelSets:
            header.append('Top' + str(k) + '.label' + str(l))
    encodings.append(header)

    for i in fastas:
        if i[3] == 'testing':
            name, sequence, label = i[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', i[1]), i[2]
            code = [name, label]
            mySimilarity = []
            for j in range(len(trainData)):
                if name != trainData[j][0]:
                    mySimilarity.append([myLabel[trainData[j][0]], CalculateSimilarity(re.sub('[^ARNDCQEGHILKMFPSTWYV]', '', trainData[j][1]), sequence)])
            mySimilarity = np.array(mySimilarity)
            mySimilarity = mySimilarity[np.lexsort(-mySimilarity.T)]
            for j in kNum:
                code = code + CalculateContent(mySimilarity, j, myLabelSets)
            encodings.append(code)
    return encodings


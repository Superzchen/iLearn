#!/usr/bin/env python
# _*_coding:utf-8_*_

import re
import itertools

myDiIndex = {
    'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3,
    'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
    'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11,
    'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15
}
myTriIndex = {
    'AAA': 0, 'AAC': 1, 'AAG': 2, 'AAT': 3,
    'ACA': 4, 'ACC': 5, 'ACG': 6, 'ACT': 7,
    'AGA': 8, 'AGC': 9, 'AGG': 10, 'AGT': 11,
    'ATA': 12, 'ATC': 13, 'ATG': 14, 'ATT': 15,
    'CAA': 16, 'CAC': 17, 'CAG': 18, 'CAT': 19,
    'CCA': 20, 'CCC': 21, 'CCG': 22, 'CCT': 23,
    'CGA': 24, 'CGC': 25, 'CGG': 26, 'CGT': 27,
    'CTA': 28, 'CTC': 29, 'CTG': 30, 'CTT': 31,
    'GAA': 32, 'GAC': 33, 'GAG': 34, 'GAT': 35,
    'GCA': 36, 'GCC': 37, 'GCG': 38, 'GCT': 39,
    'GGA': 40, 'GGC': 41, 'GGG': 42, 'GGT': 43,
    'GTA': 44, 'GTC': 45, 'GTG': 46, 'GTT': 47,
    'TAA': 48, 'TAC': 49, 'TAG': 50, 'TAT': 51,
    'TCA': 52, 'TCC': 53, 'TCG': 54, 'TCT': 55,
    'TGA': 56, 'TGC': 57, 'TGG': 58, 'TGT': 59,
    'TTA': 60, 'TTC': 61, 'TTG': 62, 'TTT': 63
}

baseSymbol = 'ACGT'


def get_kmer_frequency(sequence, kmer):
    myFrequency = {}
    for pep in [''.join(i) for i in list(itertools.product(baseSymbol, repeat=kmer))]:
        myFrequency[pep] = 0
    for i in range(len(sequence) - kmer + 1):
        myFrequency[sequence[i: i + kmer]] = myFrequency[sequence[i: i + kmer]] + 1
    for key in myFrequency:
        myFrequency[key] = myFrequency[key] / (len(sequence) - kmer + 1)
    return myFrequency


def correlationFunction(pepA, pepB, myIndex, myPropertyName, myPropertyValue):
    CC = 0
    for p in myPropertyName:
        CC = CC + (float(myPropertyValue[p][myIndex[pepA]]) - float(myPropertyValue[p][myIndex[pepB]])) ** 2
    return CC / len(myPropertyName)


def correlationFunction_type2(pepA, pepB, myIndex, myPropertyName, myPropertyValue):
    CC = 0
    for p in myPropertyName:
        CC = CC + float(myPropertyValue[p][myIndex[pepA]]) * float(myPropertyValue[p][myIndex[pepB]])
    return CC


def get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, kmer):
    thetaArray = []
    for tmpLamada in range(lamadaValue):
        theta = 0
        for i in range(len(sequence) - tmpLamada - kmer):
            theta = theta + correlationFunction(sequence[i:i + kmer],
                                                sequence[i + tmpLamada + 1: i + tmpLamada + 1 + kmer], myIndex,
                                                myPropertyName, myPropertyValue)
        thetaArray.append(theta / (len(sequence) - tmpLamada - kmer))
    return thetaArray


def get_theta_array_type2(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, kmer):
    thetaArray = []
    for tmpLamada in range(lamadaValue):
        for p in myPropertyName:
            theta = 0
            for i in range(len(sequence) - tmpLamada - kmer):
                theta = theta + correlationFunction_type2(sequence[i:i + kmer],
                                                          sequence[i + tmpLamada + 1: i + tmpLamada + 1 + kmer],
                                                          myIndex,
                                                          [p], myPropertyValue)
            thetaArray.append(theta / (len(sequence) - tmpLamada - kmer))
    return thetaArray


def make_PseDNC_vector(fastas, myPropertyName, myPropertyValue, lamadaValue, weight):
    encodings = []
    myIndex = myDiIndex
    header = ['#', 'label']
    for pair in sorted(myIndex):
        header.append(pair)
    for k in range(1, lamadaValue + 1):
        header.append('lamada_' + str(k))
    encodings.append(header)
    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        dipeptideFrequency = get_kmer_frequency(sequence, 2)
        thetaArray = get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, 2)
        for pair in sorted(myIndex.keys()):
            code.append(dipeptideFrequency[pair] / (1 + weight * sum(thetaArray)))
        for k in range(17, 16 + lamadaValue + 1):
            code.append((weight * thetaArray[k - 17]) / (1 + weight * sum(thetaArray)))
        encodings.append(code)
    return encodings


def make_PseKNC_vector(fastas, myPropertyName, myPropertyValue, lamadaValue, weight, kmer):
    encodings = []
    myIndex = myDiIndex
    header = ['#', 'label']
    header = header + sorted([''.join(i) for i in list(itertools.product(baseSymbol, repeat=kmer))])
    for k in range(1, lamadaValue + 1):
        header.append('lamada_' + str(k))
    encodings.append(header)
    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        kmerFreauency = get_kmer_frequency(sequence, kmer)
        thetaArray = get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, 2)
        for pep in sorted([''.join(j) for j in list(itertools.product(baseSymbol, repeat=kmer))]):
            code.append(kmerFreauency[pep] / (1 + weight * sum(thetaArray)))
        for k in range(len(baseSymbol) ** kmer + 1, len(baseSymbol) ** kmer + lamadaValue + 1):
            code.append((weight * thetaArray[k - (len(baseSymbol) ** kmer + 1)]) / (1 + weight * sum(thetaArray)))
        encodings.append(code)
    return encodings


def make_PCPseTNC_vector(fastas, myPropertyName, myPropertyValue, lamadaValue, weight):
    encodings = []
    myIndex = myTriIndex

    header = ['#', 'label']
    for tripeptide in sorted(myIndex):
        header.append(tripeptide)
    for k in range(1, lamadaValue + 1):
        header.append('lamada_' + str(k))
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        tripeptideFrequency = get_kmer_frequency(sequence, 3)
        thetaArray = get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, 3)
        for pep in sorted(myIndex.keys()):
            code.append(tripeptideFrequency[pep] / (1 + weight * sum(thetaArray)))
        for k in range(65, 64 + lamadaValue + 1):
            code.append((weight * thetaArray[k - 65]) / (1 + weight * sum(thetaArray)))
        encodings.append(code)
    return encodings


def make_SCPseDNC_vector(fastas, myPropertyName, myPropertyValue, lamadaValue, weight):
    encodings = []
    myIndex = myDiIndex
    header = ['#', 'label']
    for pair in sorted(myIndex):
        header.append(pair)
    for k in range(1, lamadaValue * len(myPropertyName) + 1):
        header.append('lamada_' + str(k))

    encodings.append(header)
    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        dipeptideFrequency = get_kmer_frequency(sequence, 2)
        thetaArray = get_theta_array_type2(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, 2)
        for pair in sorted(myIndex.keys()):
            code.append(dipeptideFrequency[pair] / (1 + weight * sum(thetaArray)))
        for k in range(17, 16 + lamadaValue * len(myPropertyName) + 1):
            code.append((weight * thetaArray[k - 17]) / (1 + weight * sum(thetaArray)))
        encodings.append(code)
    return encodings


def make_SCPseTNC_vector(fastas, myPropertyName, myPropertyValue, lamadaValue, weight):
    encodings = []
    myIndex = myTriIndex
    header = ['#', 'label']
    for pep in sorted(myIndex):
        header.append(pep)
    for k in range(1, lamadaValue * len(myPropertyName) + 1):
        header.append('lamada_' + str(k))
    encodings.append(header)
    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        tripeptideFrequency = get_kmer_frequency(sequence, 3)
        thetaArray = get_theta_array_type2(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, 3)
        for pep in sorted(myIndex.keys()):
            code.append(tripeptideFrequency[pep] / (1 + weight * sum(thetaArray)))
        for k in range(65, 64 + lamadaValue * len(myPropertyName) + 1):
            code.append((weight * thetaArray[k - 65]) / (1 + weight * sum(thetaArray)))
        encodings.append(code)
    return encodings

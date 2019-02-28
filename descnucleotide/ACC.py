#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import sys

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

def generatePropertyPairs(myPropertyName):
    pairs = []
    for i in range(len(myPropertyName)):
        for j in range(i+1, len(myPropertyName)):
            pairs.append([myPropertyName[i], myPropertyName[j]])
            pairs.append([myPropertyName[j], myPropertyName[i]])

    return pairs

def make_ac_vector(fastas, myPropertyName, myPropertyValue, lag, kmer):
    encodings = []
    myIndex = myDiIndex if kmer == 2 else myTriIndex
    header = ['#', 'label']
    for p in myPropertyName:
        for l in range(1, lag + 1):
            header.append('%s.lag%d' %(p, l))
    encodings.append(header)
    
    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        
        for p in myPropertyName:
            meanValue = 0
            #for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                meanValue = meanValue + float(myPropertyValue[p][myIndex[sequence[j: j+kmer]]])
            #meanValue = meanValue / (len(sequence) - kmer)
            meanValue = meanValue / (len(sequence) - kmer + 1)
            
            for l in range(1, lag + 1):
                acValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    #acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j+kmer]]]) - meanValue) * (float(myPropertyValue[p][myIndex[sequence[j+l:j+l+kmer]]]))
                    acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j + kmer]]]) - meanValue) * (
                    float(myPropertyValue[p][myIndex[sequence[j + l:j + l + kmer]]]) - meanValue)
                acValue = acValue / (len(sequence) - kmer - l + 1)
                #print(acValue)
                code.append(acValue)
        encodings.append(code)
    return encodings

def make_cc_vector(fastas, myPropertyName, myPropertyValue, lag, kmer):
    encodings = []
    myIndex = myDiIndex if kmer == 2 else myTriIndex
    if len(myPropertyName) < 2:
        print('Error: two or more property are needed for cross covariance (i.e. DCC and TCC) descriptors')
        sys.exit(1)
    propertyPairs = generatePropertyPairs(myPropertyName)
    header = ['#', 'label'] + [n[0] + '-' + n[1] + '-lag.' + str(l) for n in propertyPairs for l in range(1, lag + 1)]
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]

        for pair in propertyPairs:
            meanP1 = 0
            meanP2 = 0
            #for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                meanP1 = meanP1 + float(myPropertyValue[pair[0]][myIndex[sequence[j: j+kmer]]])
                meanP2 = meanP2 + float(myPropertyValue[pair[1]][myIndex[sequence[j: j+kmer]]])
            #meanP1 = meanP1 / (len(sequence) - kmer)
            #meanP2 = meanP2 / (len(sequence) - kmer)
            meanP1 = meanP1 / (len(sequence) - kmer + 1)
            meanP2 = meanP2 / (len(sequence) - kmer + 1)

            for l in range(1, lag + 1):
                ccValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    ccValue = ccValue + (float(myPropertyValue[pair[0]][myIndex[sequence[j: j + kmer]]]) - meanP1) * (
                        float(myPropertyValue[pair[1]][myIndex[sequence[j + l:j + l + kmer]]]) - meanP2)
                ccValue = ccValue / (len(sequence) - kmer - l + 1)
                code.append(ccValue)
        encodings.append(code)
    return encodings

def make_acc_vector(fastas, myPropertyName, myPropertyValue, lag, kmer):
    encodings = []
    myIndex = myDiIndex if kmer == 2 else myTriIndex
    if len(myPropertyName) < 2:
        print('Error: two or more property are needed for cross covariance (i.e. DCC and TCC) descriptors')
        sys.exit(1)

    header = ['#', 'label']
    for p in myPropertyName:
        for l in range(1, lag + 1):
            header.append('%s.lag%d' %(p, l))
    propertyPairs = generatePropertyPairs(myPropertyName)
    header = header + [n[0] + '-' + n[1] + '-lag.' + str(l) for n in propertyPairs for l in range(1, lag + 1)]
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        ## Auto covariance
        for p in myPropertyName:
            meanValue = 0
            # for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                meanValue = meanValue + float(myPropertyValue[p][myIndex[sequence[j: j + kmer]]])
            # meanValue = meanValue / (len(sequence) - kmer)
            meanValue = meanValue / (len(sequence) - kmer + 1)

            for l in range(1, lag + 1):
                acValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    # acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j+kmer]]]) - meanValue) * (float(myPropertyValue[p][myIndex[sequence[j+l:j+l+kmer]]]))
                    acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j + kmer]]]) - meanValue) * (
                        float(myPropertyValue[p][myIndex[sequence[j + l:j + l + kmer]]]) - meanValue)
                acValue = acValue / (len(sequence) - kmer - l + 1)
                # print(acValue)
                code.append(acValue)

        ## Cross covariance
        for pair in propertyPairs:
            meanP1 = 0
            meanP2 = 0
            #for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                meanP1 = meanP1 + float(myPropertyValue[pair[0]][myIndex[sequence[j: j+kmer]]])
                meanP2 = meanP2 + float(myPropertyValue[pair[1]][myIndex[sequence[j: j+kmer]]])
            #meanP1 = meanP1 / (len(sequence) - kmer)
            #meanP2 = meanP2 / (len(sequence) - kmer)
            meanP1 = meanP1 / (len(sequence) - kmer + 1)
            meanP2 = meanP2 / (len(sequence) - kmer + 1)

            for l in range(1, lag + 1):
                ccValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    ccValue = ccValue + (float(myPropertyValue[pair[0]][myIndex[sequence[j: j + kmer]]]) - meanP1) * (
                        float(myPropertyValue[pair[1]][myIndex[sequence[j + l:j + l + kmer]]]) - meanP2)
                ccValue = ccValue / (len(sequence) - kmer - l + 1)
                code.append(ccValue)
        encodings.append(code)
    return encodings
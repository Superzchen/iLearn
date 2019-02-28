#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def DNC(fastas, **kw):
    base = 'ACGT'

    encodings = []
    dinucleotides = [n1 + n2 for n1 in base for n2 in base]
    header = ['#', 'label'] + dinucleotides
    encodings.append(header)

    AADict = {}
    for i in range(len(base)):
        AADict[base[i]] = i

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        tmpCode = [0] * 16
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings
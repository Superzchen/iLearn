#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys, os, re

def TriNcleotideComposition(sequence, base):
    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    tnc_dict = {}
    for triN in trincleotides:
        tnc_dict[triN] = 0
    for i in range(len(sequence) - 2):
        tnc_dict[sequence[i:i + 3]] += 1
    for key in tnc_dict:
       tnc_dict[key] /= (len(sequence) - 2)
    return tnc_dict

def PseEIIP(fastas, **kw):
    for i in fastas:
        if re.search('[^ACGT-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT-" are allowed by this PseEIIP scheme.')
            return 0

    base = 'ACGT'

    EIIP_dict = {
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'T': 0.1335,
    }

    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    EIIPxyz = {}
    for triN in trincleotides:
        EIIPxyz[triN] = EIIP_dict[triN[0]] + EIIP_dict[triN[1]] + EIIP_dict[triN[2]]

    encodings = []
    header = ['#', 'label'] + trincleotides
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        trincleotide_frequency = TriNcleotideComposition(sequence, base)
        code = code + [EIIPxyz[triN] * trincleotide_frequency[triN] for triN in trincleotides]
        encodings.append(code)
    return encodings

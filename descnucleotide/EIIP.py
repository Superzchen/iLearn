#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, platform
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import check_sequences

def EIIP(fastas, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "EIIP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGT'

    EIIP_dict ={
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'T': 0.1335,
        '-': 0,
    }

    encodings = []
    header = ['#', 'label']
    for i in range(1, len(fastas[0][1]) + 1):
        header.append('F'+str(i))
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], i[1], i[2]
        code = [name, label]
        for aa in sequence:
            code.append(EIIP_dict.get(aa, 0))
        encodings.append(code)
    return encodings
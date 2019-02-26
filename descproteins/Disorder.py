#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re, platform
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import check_sequences

def Disorder(fastas, **kw):
    if check_sequences.check_fasta_with_equal_length(fastas) == False:
        print('Error: for "Disorder" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    encodings = []
    header = ['#', 'label']
    for p in range(1, len(fastas[0][1])+1):
        header.append('disorder.F' + str(p))
    encodings.append(header)

    disDir = kw['path']
    if disDir == None:
        print('Error: please specify the directory of predicted protein disorder files by "--path"')
        return 0
    for i in fastas:
        name, sequence, label = i[0], i[1], i[2]
        code = [name, label]
        if os.path.exists(disDir + '/' + name + '.dis') == False:
            print('Error: the predicted disorder information file (.dis) for protein ' + name + ' does not exist.')
            return 0

        with open(disDir + '/' + name + '.dis') as f:
            records = f.readlines()
        tag = 0
        for i in range(len(records)):
            if re.search('^-------', records[i]):
                tag = i
                break
        records = records[tag+1:-1]

        proteinSeq = ''
        disValue = []
        for line in records:
            array = line.rstrip().split() if line.rstrip() != '' else None
            proteinSeq = proteinSeq + array[1]
            disValue.append(array[2])
        pos = proteinSeq.find(sequence)
        if pos == -1:
            print('Warning: could not find the peptide in proteins.\n\n')
        else:
            for p in range(pos, pos+len(sequence)):
                code.append(disValue[p])
        encodings.append(code)

    return encodings
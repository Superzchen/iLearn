#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys, platform, os, re
import numpy as np
import argparse

pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import read_fasta_sequences
import save_file
import check_sequences


def SOCNumber(fastas, nlag=30, **kw):
    if check_sequences.get_min_sequence_length_1(fastas) < nlag + 1:
        print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
        return 0

    dataFile = re.sub('descproteins$', '', os.path.split(os.path.realpath(__file__))[
        0]) + r'\data\Schneider-Wrede.txt' if platform.system() == 'Windows' else re.sub('descproteins$', '', os.path.split(
        os.path.realpath(__file__))[0]) + r'/data/Schneider-Wrede.txt'
    dataFile1 = re.sub('descproteins$', '', os.path.split(os.path.realpath(__file__))[
        0]) + r'\data\Grantham.txt' if platform.system() == 'Windows' else re.sub('descproteins$', '', os.path.split(
        os.path.realpath(__file__))[0]) + r'/data/Grantham.txt'
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    AA1 = 'ARNDCQEGHILKMFPSTWYV'

    DictAA = {}
    for i in range(len(AA)):
        DictAA[AA[i]] = i

    DictAA1 = {}
    for i in range(len(AA1)):
        DictAA1[AA1[i]] = i

    with open(dataFile) as f:
        records = f.readlines()[1:]
    AADistance = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance.append(array)
    AADistance = np.array(
        [float(AADistance[i][j]) for i in range(len(AADistance)) for j in range(len(AADistance[i]))]).reshape((20, 20))

    with open(dataFile1) as f:
        records = f.readlines()[1:]
    AADistance1 = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance1.append(array)
    AADistance1 = np.array(
        [float(AADistance1[i][j]) for i in range(len(AADistance1)) for j in range(len(AADistance1[i]))]).reshape(
        (20, 20))

    encodings = []
    header = ['#', 'label']
    for n in range(1, nlag + 1):
        header.append('Schneider.lag' + str(n))
    for n in range(1, nlag + 1):
        header.append('gGrantham.lag' + str(n))
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        for n in range(1, nlag + 1):
            code.append(sum(
                [AADistance[DictAA[sequence[j]]][DictAA[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]) / (
                            len(sequence) - n))

        for n in range(1, nlag + 1):
            code.append(sum([AADistance1[DictAA1[sequence[j]]][DictAA1[sequence[j + n]]] ** 2 for j in
                             range(len(sequence) - n)]) / (len(sequence) - n))
        encodings.append(code)
    return encodings


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating SOCNumber feature vector for nucleotide sequences")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--lag", type=int, default=30, help="the lag value for SOCNumber descriptor")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()
    output = args.out if args.out != None else 'encoding.txt'
    kw = {'order': 'ACDEFGHIKLMNPQRSTVWY'}
    fastas = read_fasta_sequences.read_protein_sequences(args.file)
    encodings = SOCNumber(fastas, nlag=args.lag, **kw)
    save_file.save_file(encodings, args.format, output)

#!/usr/bin/env python
# _*_coding:utf-8_*_

import re, sys, os, platform
import math
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


def Rvalue(aa1, aa2, AADict, Matrix):
    return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)


def PAAC(fastas, lambdaValue=30, w=0.05, **kw):
    if check_sequences.get_min_sequence_length_1(fastas) < lambdaValue + 1:
        print(
            'Error: all the sequence length should be larger than the lambdaValue+1: ' + str(lambdaValue + 1) + '\n\n')
        return 0

    dataFile = re.sub('descproteins$', '', os.path.split(os.path.realpath(__file__))[
        0]) + r'\data\PAAC.txt' if platform.system() == 'Windows' else re.sub('descproteins$', '',
                                                                              os.path.split(os.path.realpath(__file__))[
                                                                                  0]) + r'/data/PAAC.txt'

    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])

    encodings = []
    header = ['#', 'label']
    for aa in AA:
        header.append('Xc1.' + aa)
    for n in range(1, lambdaValue + 1):
        header.append('Xc2.lambda' + str(n))
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        theta = []
        for n in range(1, lambdaValue + 1):
            theta.append(
                sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (
                    len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
    return encodings

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating PAAC feature vector for nucleotide sequences")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--lamada", type=int, default=30, help="the lamada value for SOCNumber descriptor")
    parser.add_argument("--weight", type=float, default=0.05, help="the weight value for SOCNumber descriptor")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()
    output = args.out if args.out != None else 'encoding.txt'
    kw = {'order': 'ACDEFGHIKLMNPQRSTVWY'}
    fastas = read_fasta_sequences.read_protein_sequences(args.file)
    encodings = PAAC(fastas, lambdaValue=args.lamada, w=args.weight, **kw)
    save_file.save_file(encodings, args.format, output)
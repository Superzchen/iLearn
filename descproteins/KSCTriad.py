#!/usr/bin/env python
# _*_coding:utf-8_*_

import re, sys, os, platform
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


def CalculateKSCTriad(sequence, gap, features, AADict):
    res = []
    for g in range(gap + 1):
        myDict = {}
        for f in features:
            myDict[f] = 0

        for i in range(len(sequence)):
            if i + g + 1 < len(sequence) and i + 2 * g + 2 < len(sequence):
                fea = AADict[sequence[i]] + '.' + AADict[sequence[i + g + 1]] + '.' + AADict[
                    sequence[i + 2 * g + 2]]
                myDict[fea] = myDict[fea] + 1

        maxValue, minValue = max(myDict.values()), min(myDict.values())
        for f in features:
            res.append((myDict[f] - minValue) / maxValue)

    return res


def KSCTriad(fastas, gap=0, **kw):
    AAGroup = {
        'g1': 'AGV',
        'g2': 'ILFP',
        'g3': 'YMTS',
        'g4': 'HNQW',
        'g5': 'RK',
        'g6': 'DE',
        'g7': 'C'
    }

    myGroups = sorted(AAGroup.keys())

    AADict = {}
    for g in myGroups:
        for aa in AAGroup[g]:
            AADict[aa] = g

    features = [f1 + '.' + f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]

    encodings = []
    header = ['#', 'label']
    for g in range(gap + 1):
        for f in features:
            header.append(f + '.gap' + str(g))
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        if len(sequence) < 2 * gap + 3:
            print('Error: for "KSCTriad" encoding, the input fasta sequences should be greater than (2*gap+3). \n\n')
            return 0
        code = code + CalculateKSCTriad(sequence, gap, features, AADict)
        encodings.append(code)

    return encodings


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating KSCTriad feature vector for nucleotide sequences")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--gap", type=int, default=5, help="the k-space value for KSCTriad descriptor")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()
    output = args.out if args.out != None else 'encoding.txt'
    kw = {}
    fastas = read_fasta_sequences.read_protein_sequences(args.file)
    encodings = KSCTriad(fastas, gap=args.gap, **kw)
    save_file.save_file(encodings, args.format, output)

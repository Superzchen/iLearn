#!/usr/bin/env python
# _*_coding:utf-8_*_

import re, sys, os, platform
from collections import Counter
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


def ENAC(fastas, window=5, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "ENAC" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    if window < 1:
        print('Error: the sliding window should be greater than zero' + '\n\n')
        return 0

    if check_sequences.get_min_sequence_length(fastas) < window:
        print('Error: all the sequence length should be larger than the sliding window :' + str(window) + '\n\n')
        return 0

    AA = kw['order'] if kw['order'] != None else 'ACGT'
    encodings = []
    header = ['#', 'label']
    for w in range(1, len(fastas[0][1]) - window + 2):
        for aa in AA:
            header.append('SW.' + str(w) + '.' + aa)
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], i[1], i[2]
        code = [name, label]
        for j in range(len(sequence)):
            if j < len(sequence) and j + window <= len(sequence):
                count = Counter(sequence[j:j + window])
                for key in count:
                    count[key] = count[key] / len(sequence[j:j + window])
                for aa in AA:
                    code.append(count[aa])
        encodings.append(code)
    return encodings

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating ENAC feature vector for nucleotide sequences")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--slwindow", type=int, default=5, help="the sliding window of ENAC descriptor")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()
    output = args.out if args.out != None else 'encoding.txt'
    kw = {'order': 'ACGT'}
    fastas = read_fasta_sequences.read_nucleotide_sequences(args.file)
    encodings = ENAC(fastas, window=args.slwindow, **kw)
    save_file.save_file(encodings, args.format, output)

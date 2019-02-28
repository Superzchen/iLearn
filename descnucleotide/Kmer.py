#!/usr/bin/env python
#_*_coding:utf-8_*_

import argparse
import re, sys, os, platform
import itertools
from collections import Counter

pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import read_fasta_sequences
import save_file


def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer


def Kmer(fastas, k=2, type="DNA", upto=False, normalize=True, **kw):
    encoding = []
    header = ['#', 'label']
    NA = 'ACGT'
    if type in ("DNA", 'RNA'):
        NA = 'ACGT'
    else:
        NA = 'ACDEFGHIKLMNPQRSTVWY'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            for kmer in itertools.product(NA, repeat=tmpK):
                header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            code = [name, label]
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        for kmer in itertools.product(NA, repeat=k):
            header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            kmers = kmerArray(sequence, k)
            count = Counter()
            count.update(kmers)
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            code = [name, label]
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    return encoding


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating Kmer feature vector for nucleotide sequences")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--kmer", type=int, default=2, help="the k value of kmer")
    parser.add_argument("--type", choices=['DNA', 'RNA', 'Protein'], default='DNA', help="input sequence type")
    parser.add_argument("--upto", action='store_true',
                        help="will generate all the kmers: 1mer, 2mer, ..., kmer")
    parser.add_argument("--normalize", action='store_true',
                        help="with this option, the final feature vector will be normalized based on the total occurrences of all kmers")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()
    kmer = args.kmer
    output = args.out if args.out != None else 'encoding.txt'
    kw = {}
    fastas = read_fasta_sequences.read_nucleotide_sequences(args.file)
    encodings = Kmer(fastas, kmer, args.type, args.upto, args.normalize, **kw)
    save_file.save_file(encodings, args.format, output)


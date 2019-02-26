#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
from descnucleotide import *
from pubscripts import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating various numerical representation schemes for nucleotide sequences.")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--method", required=True,
                        choices=['Kmer', 'RCKmer', 'NAC', 'DNC', 'TNC', 'ANF', 'ENAC', 'binary', 'CKSNAP', 'NCP',
                                 'PSTNPss', 'PSTNPds', 'EIIP', 'PseEIIP'],
                        help="the encoding type")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka', 'tsv_1'], default='svm',
                        help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()
    fastas = read_fasta_sequences.read_nucleotide_sequences(args.file)
    kw = {'order': 'ACGT', }
    cmd = args.method + '.' + args.method + '(fastas, **kw)'
    print('Descriptor method: ' + args.method)
    encodings = eval(cmd)
    out_file = args.out if args.out != None else 'encoding.txt'
    save_file.save_file(encodings, args.format, out_file)

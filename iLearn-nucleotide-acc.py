#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
from descnucleotide import *
from pubscripts import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating auto-correlation encoding for nucleotide sequences.")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--method", required=True,
                        choices=['DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC'],
                        help="the encoding method")
    parser.add_argument("--type", choices=['DNA', 'RNA'], default='DNA', help="the nucleotide, default: DNA.")
    parser.add_argument('--lag', type=int, default=2, help="The value of lag.")
    parser.add_argument('--index',
                        help="The indices file user choose.\n"
                             "Default indices:\n"
                             "DNA dinucleotide: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                             "DNA trinucleotide: Dnase I, Bendability (DNAse).\n"
                             "RNA: Rise, Roll, Shift, Slide, Tilt, Twist.\n")
    parser.add_argument("--udi", help="The user-defined indices file.")
    parser.add_argument('--all_index', action='store_true', help="Choose all physico-chemical indices, default: False.")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka', 'tsv_1'], default='svm',
                        help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file.")
    args = parser.parse_args()
    my_property_name, my_property_value, kmer = check_parameters.check_acc_arguments(args)
    fastas = read_fasta_sequences.read_nucleotide_sequences(args.file)
    encodings = []
    if args.method == 'DAC' or args.method == 'TAC':
        encodings = ACC.make_ac_vector(fastas, my_property_name, my_property_value, args.lag, kmer)
    if args.method == 'DCC' or args.method == 'TCC':
        encodings = ACC.make_cc_vector(fastas, my_property_name, my_property_value, args.lag, kmer)
    if args.method == 'DACC' or args.method == 'TACC':
        encodings = ACC.make_acc_vector(fastas, my_property_name, my_property_value, args.lag, kmer)
    out_file = args.out if args.out != None else 'encoding.txt'
    save_file.save_file(encodings, args.format, out_file)

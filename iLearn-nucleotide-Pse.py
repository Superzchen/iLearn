#!/usr/bin/env python
# _*_coding:utf-8_*_

import argparse
from descnucleotide import *
from pubscripts import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating pseudo nucleic acid composition encoding for nucleotide sequences.")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--method", required=True,
                        choices=['PseDNC', 'PseKNC', 'PCPseDNC', 'PCPseTNC', 'SCPseDNC', 'SCPseTNC'],
                        help="the encoding type")
    parser.add_argument("--type", choices=['DNA', 'RNA'], help="the nucleotide, default: DNA.")
    parser.add_argument('--lamada', dest='lamadaValue', type=int, default=2, help="The value of lamada; default: 2")
    parser.add_argument('--weight', type=float, default=0.1, help="The value of weight; default: 0.1")
    parser.add_argument('--kmer', type=int, default=3, help="The value of kmer; it works only with PseKNC method.")
    parser.add_argument('--index',
                        help="The indices file user choose.\n"
                             "Default indices:\n"
                             "DNA dinucleotide: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                             "DNA trinucleotide: Dnase I, Bendability (DNAse).\n"
                             "RNA: Rise, Roll, Shift, Slide, Tilt, Twist.\n")
    parser.add_argument("--udi", help="The user-defined indices file.")
    parser.add_argument('--all_index', action='store_true', help="Choose all physico-chemical indices, default: False.")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka', 'tsv_1'], default='svm',
                        help="the output format")
    parser.add_argument("--out", help="the generated descriptor file.")
    parser.set_defaults(type='DNA')
    args = parser.parse_args()

    fastas = read_fasta_sequences.read_nucleotide_sequences(args.file)
    my_property_name, my_property_value, lamada_value, weight, kmer = check_parameters.check_Pse_arguments(args, fastas)
    encodings = []
    if args.method == 'PseDNC':
        encodings = Pse.make_PseDNC_vector(fastas, my_property_name, my_property_value, lamada_value, weight)
    if args.method == 'PseKNC':
        encodings = Pse.make_PseKNC_vector(fastas, my_property_name, my_property_value, lamada_value, weight, kmer)
    if args.method == 'PCPseDNC':  # PseDNC is identical to PC-PseDNC
        encodings = Pse.make_PseDNC_vector(fastas, my_property_name, my_property_value, lamada_value, weight)
    if args.method == 'PCPseTNC':
        encodings = Pse.make_PCPseTNC_vector(fastas, my_property_name, my_property_value, lamada_value, weight)
    if args.method == 'SCPseDNC':
        encodings = Pse.make_SCPseDNC_vector(fastas, my_property_name, my_property_value, lamada_value, weight)
    if args.method == 'SCPseTNC':
        encodings = Pse.make_SCPseTNC_vector(fastas, my_property_name, my_property_value, lamada_value, weight)
    out_file = args.out if args.out != None else 'encoding.txt'
    save_file.save_file(encodings, args.format, out_file)

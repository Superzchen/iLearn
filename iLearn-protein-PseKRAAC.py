#!/usr/bin/env python
# _*_coding:utf-8_*_

import argparse, sys
from descproteins import *
from pubscripts import *

USAGE = """The 'raactype' value for each subtype descriptor could be chosen from:
    type1    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type2    [2, 3, 4, 5, 6,    8,                        15,                 20]
    type3A   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type3B   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type4    [         5,       8, 9,     11,     13,                         20]
    type5    [   3, 4,          8,    10,                 15,                 20]
    type6A   [      4, 5,                                                     20]
    type6B   [         5,                                                       ]
    type6C   [         5,                                                       ]
    type7    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type8    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type9    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type10   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type11   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type12   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,     20]
    type13   [      4,                        12,                 17,         20]
    type14   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    type15   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20]
    type16   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20]    
"""

USAGEHASH = {
    'type1': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type2': [2, 3, 4, 5, 6, 8, 15, 20],
    'type3A': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type3B': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type4': [5, 8, 9, 11, 13, 20],
    'type5': [3, 4, 8, 10, 15, 20],
    'type6A': [4, 5, 20],
    'type6B': [5, ],
    'type6C': [5, ],
    'type7': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type8': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type9': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type10': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type11': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type12': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20],
    'type13': [4, 12, 17, 20],
    'type14': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'type15': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20],
    'type16': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20],
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating PseKRAAC descriptors for protein sequences/peptides:")
    parser.add_argument("--file", help="input fasta file")
    parser.add_argument("--method",
                        choices=['type1', 'type2', 'type3A', 'type3B', 'type4', 'type5', 'type6A', 'type6B', 'type6C',
                                 'type7', 'type8', 'type9', 'type10', 'type11', 'type12', 'type13', 'type14', 'type15',
                                 'type16'], help="the descriptor type")
    parser.add_argument("--model", default='g-gap', choices=['g-gap', 'lambda-correlation'],
                        help="the model of the descriptor method, default is 'g-gap'")
    parser.add_argument("--ktuple", default=2, choices=[1, 2, 3], help="k-tuple peptide, default is 2", type=int)
    parser.add_argument("--gap_lambda", choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                        help="the gap value or lambda value for the 'g-gap' model or 'lambda-correlation' model",
                        type=int)
    parser.add_argument("--type", help="the reduced amino acids cluster type", type=int)
    parser.add_argument("--show", help="show detatiled available '--type' value for each type",
                        action="store_true")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka', 'tsv_1'], default='svm',
                        help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()

    if args.show:
        print(USAGE)
        sys.exit(1)

    if args.file == None:
        print('The following arguments are required: --file\n')
        sys.exit(1)
    if args.method == None:
        print('The following arguments are required: --medhod\n')
        sys.exit(1)
    if args.type == None:
        print('The following arguments are required: --type\n')
        sys.exit(1)
    if args.gap_lambda == None:
        print('The following arguments are required: --gap_lambda\n')
        sys.exit(1)

    if args.type not in USAGEHASH[args.method]:
        print('The "type" value error. For detailed parameter, please see:\n\n')
        print(USAGE)
        sys.exit(1)

    fastas = read_fasta_sequences.read_protein_sequences(args.file)
    cmd = args.method + '.type1' + '(fastas, args.model, args.type, args.ktuple, args.gap_lambda)'
    print(cmd)
    print('Descriptor method: ' + args.method)
    print('Descriptor model: ' + args.model)
    print('Reduced amino acids cluster type: ' + str(args.type))
    print('k_tuple peptide number: ' + str(args.ktuple))
    print('gap or lambda value: ' + str(args.gap_lambda) + '\n\n')
    encodings = eval(cmd)
    out_file = args.out if args.out != None else 'encoding.txt'
    save_file.save_file(encodings, args.format, out_file)

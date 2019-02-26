#!/usr/bin/env python
# _*_coding:utf-8_*_

import argparse
import sys, os, re, platform

pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import check_sequences
import read_fasta_sequences
import save_file


def AAINDEX(fastas, props=None, **kw):
    if check_sequences.check_fasta_with_equal_length(fastas) == False:
        print('Error: for "AAINDEX" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ARNDCQEGHILKMFPSTWYV'

    fileAAindex = re.sub('descproteins$', '', os.path.split(os.path.realpath(__file__))[
        0]) + r'\data\AAindex.txt' if platform.system() == 'Windows' else re.sub('descproteins$', '', os.path.split(
        os.path.realpath(__file__))[0]) + r'/data/AAindex.txt'
    with open(fileAAindex) as f:
        records = f.readlines()[1:]

    AAindex = []
    AAindexName = []
    for i in records:
        AAindex.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
        AAindexName.append(i.rstrip().split()[0] if i.rstrip() != '' else None)

    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i

    #  use the user inputed properties
    if props:
        tmpIndexNames = []
        tmpIndex = []
        for p in props:
            if AAindexName.index(p) != -1:
                tmpIndexNames.append(p)
                tmpIndex.append(AAindex[AAindexName.index(p)])
        if len(tmpIndexNames) != 0:
            AAindexName = tmpIndexNames
            AAindex = tmpIndex
    
    encodings = []
    header = ['#', 'label']
    for pos in range(1, len(fastas[0][1]) + 1):
        for idName in AAindexName:
            header.append('SeqPos.' + str(pos) + '.' + idName)
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], i[1], i[2]
        code = [name, label]
        for aa in sequence:
            if aa == '-':
                for j in AAindex:
                    code.append(0)
                continue
            for j in AAindex:
                code.append(j[index[aa]])
        encodings.append(code)
    return encodings

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="AAINDEX descriptor")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--props", help="phy-chemical property list")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()

    fastas = read_fasta_sequences.read_protein_sequences(args.file)
    props = args.props.split(':') if args.props != None else None    
    output = args.out if args.out != None else 'encoding.txt'
    encodings = AAINDEX(fastas, props)
    save_file.save_file(encodings, args.format, output)

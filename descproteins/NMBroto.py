#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys, platform, os, re
import argparse
import numpy as np

pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import read_fasta_sequences
import save_file
import check_sequences


def NMBroto(fastas,
            props=['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102', 'CHOC760101', 'BIGC670101', 'CHAM810101',
                   'DAYM780201'], nlag=30, **kw):
    if check_sequences.get_min_sequence_length_1(fastas) < nlag + 1:
        print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
        return 0

    AA = 'ARNDCQEGHILKMFPSTWYV'
    fileAAidx = re.sub('descproteins$', '', os.path.split(os.path.realpath(__file__))[
        0]) + r'\data\AAidx.txt' if platform.system() == 'Windows' else sys.path[0] + r'/data/AAidx.txt'
    with open(fileAAidx) as f:
        records = f.readlines()[1:]
    myDict = {}
    for i in records:
        array = i.rstrip().split('\t')
        myDict[array[0]] = array[1:]

    AAidx = []
    AAidxName = []
    for i in props:
        if i in myDict:
            AAidx.append(myDict[i])
            AAidxName.append(i)
        else:
            print('"' + i + '" properties not exist.')
            return None

    AAidx1 = np.array([float(j) for i in AAidx for j in i])
    AAidx = AAidx1.reshape((len(AAidx), 20))
    pstd = np.std(AAidx, axis=1)
    pmean = np.average(AAidx, axis=1)

    for i in range(len(AAidx)):
        for j in range(len(AAidx[i])):
            AAidx[i][j] = (AAidx[i][j] - pmean[i]) / pstd[i]

    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i

    encodings = []
    header = ['#', 'label']
    for p in props:
        for n in range(1, nlag + 1):
            header.append(p + '.lag' + str(n))
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        N = len(sequence)
        for prop in range(len(props)):
            for n in range(1, nlag + 1):
                if len(sequence) > nlag:
                    # if key is '-', then the value is 0
                    rn = sum(
                        [AAidx[prop][index.get(sequence[j], 0)] * AAidx[prop][index.get(sequence[j + n], 0)] for j in
                         range(len(sequence) - n)]) / (N - n)
                else:
                    rn = 'NA'
                code.append(rn)
        encodings.append(code)
    return encodings


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Moran descriptor")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--props", help="input fasta file")
    parser.add_argument("--nlag", help="input fasta file")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the encoding type")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()

    fastas = read_fasta_sequences.read_protein_sequences(args.file)
    props = args.props.split(':') if args.props != None else ['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
                                                              'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201']
    nlag = int(args.nlag) if args.nlag != None else 30
    output = args.out if args.out != None else 'encoding.txt'
    encodings = NMBroto(fastas, props, nlag)
    save_file.save_file(encodings, args.format, output)


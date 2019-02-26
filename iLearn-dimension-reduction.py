#!/usr/bin/env python
# _*_coding:utf-8_*_

import argparse
from pubscripts import read_code, save_file, draw_plot
from dimreduction import *
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="dimension reduction")
    parser.add_argument("--file", required=True, help="input encoding file format")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm',
                        help="the encoding type")
    parser.add_argument("--method", required=True,
                        choices=['pca', 'lda', 'tsne'], help="select dimension reduction method")
    parser.add_argument("--ncomponents", default=2, type=int, help="number of n components, default 2")
    parser.add_argument("--out", default="dimension_reduction_res.txt", help="output file")
    args = parser.parse_args()

    encodings, labels = read_code.read_code(args.file, format=args.format)
    reduced_data = []
    if args.method == 'pca':
        reduced_data = pca.pca(encodings, n_components=args.ncomponents)
    if args.method == 'lda':
        reduced_data = lda.lda(encodings, labels, n_components=args.ncomponents)
    if args.method == 'tsne':
        reduced_data = tsne.tsne(encodings[1:, 1:].astype(float), no_dims=args.ncomponents)
        new_data = np.zeros((reduced_data.shape[0], reduced_data.shape[1] + 1))
        new_data[:, 1:] = reduced_data
        new_data = new_data.astype(str)
        new_data[:, 0] = np.array(encodings[1:])[:, 0]
        reduced_data = new_data.tolist()
    save_file.save_reduction_result(reduced_data, file=args.out)
    draw_plot.plot_2d(reduced_data, labels, file=args.out+'_2d.png')
    if args.ncomponents >= 3:
        draw_plot.plot_3d(reduced_data, labels, file=args.out+'_3d.png')


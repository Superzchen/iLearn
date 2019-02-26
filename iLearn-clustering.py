#!/usr/bin/env python
# _*_coding:utf-8_*_

import argparse
from pubscripts import read_code, save_file, draw_plot
from clusters import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="cluster for the generated numerical representation")
    parser.add_argument("--file", required=True, help="input encoding file")
    parser.add_argument("--method", required=True,
                        choices=['kmeans', 'gmm', 'hcluster', 'apc', 'meanshift', 'dbscan'], help="select cluster method")
    parser.add_argument("--sof", default='sample', choices=['sample', 'feature'],
                        help="cluster for sample or feature, default: sample")
    parser.add_argument("--nclusters", help="specify the cluster number for kmeans cluster method. default: 3")
    parser.add_argument("--out", help="output file")
    args = parser.parse_args()

    kw = {'nclusters': args.nclusters if args.nclusters != None else 3,
          'sof': args.sof if args.sof != None else 'sample'}
    output = args.out if args.out != None else '%s_cluster.txt' % args.method
    encodings = read_code.read_tsv_1(args.file)
    cmd = args.method + '.' + args.method + '(encodings, **kw)'
    print('Cluster method: ' + args.method)
    myCluster, e = eval(cmd)
    save_file.save_cluster_result(myCluster, e, output)
    ## t-sne plot
    draw_plot.plot_clustering_2d(encodings, myCluster, output, **kw)
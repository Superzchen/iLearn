#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os, platform, sys
import numpy as np
import pandas as pd
from itertools import cycle
from mpl_toolkits.mplot3d import Axes3D
from scipy import interp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve
pPath = os.path.split(os.path.realpath(__file__))[0]
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\clusters' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/clusters'
sys.path.append(father_path)
import tsne, pca


def plot_roc_cv(data, out, label_column=0, score_column=2):
    tprs = []
    aucs = []
    fprArray = []
    tprArray = []
    thresholdsArray = []
    mean_fpr = np.linspace(0, 1, 100)

    for i in range(len(data)):
        fpr, tpr, thresholds = roc_curve(data[i][:, label_column], data[i][:, score_column])
        fprArray.append(fpr)
        tprArray.append(tpr)
        thresholdsArray.append(thresholds)
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)

    colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'blueviolet', 'deeppink', 'cyan'])
    ## ROC plot for CV
    fig = plt.figure(0)
    for i, color in zip(range(len(fprArray)), colors):
        plt.plot(fprArray[i], tprArray[i], lw=1, alpha=0.7, color=color,
                 label='ROC fold %d (AUC = %0.2f)' % (i + 1, aucs[i]))
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Random', alpha=.8)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='blue',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.9)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')
    plt.xlim([0, 1.0])
    plt.ylim([0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig(out)
    plt.close(0)
    return mean_auc

def plot_mean_roc_cv(data_dict, out, label_column=0, score_column=2):
    fig = plt.figure(0)
    colors = cycle(['red', 'darkorange', 'aqua', 'cornflowerblue', 'blueviolet', 'deeppink', 'cyan'])
    line_styles = cycle(['-', '--', '-.'])
    for ls, color, key in zip(line_styles, colors, data_dict):
        data = data_dict[key]
        tprs = []
        aucs = []
        fprArray = []
        tprArray = []
        thresholdsArray = []
        mean_fpr = np.linspace(0, 1, 100)

        for i in range(len(data)):
            fpr, tpr, thresholds = roc_curve(data[i][:, label_column], data[i][:, score_column])
            fprArray.append(fpr)
            tprArray.append(tpr)
            thresholdsArray.append(thresholds)
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        plt.plot(mean_fpr, mean_tpr, color=color, linestyle=ls, label=r'%s: %.2f ' %(key, mean_auc), lw=2, alpha=.9)
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='b', label='Random', alpha=.8)
    plt.xlim([0, 1.0])
    plt.ylim([0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig(out)
    plt.close(0)


def plot_roc_ind(data, out, label_column=0, score_column=2):
    fprIndep, tprIndep, thresholdsIndep = roc_curve(data[:, label_column], data[:, score_column])
    ind_auc = auc(fprIndep, tprIndep)
    fig = plt.figure(0)
    plt.plot(fprIndep, tprIndep, lw=2, alpha=0.7, color='red',
             label='ROC curve (area = %0.2f)' % ind_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig(out)
    plt.close(0)
    return ind_auc

def plot_roc_muti_ind(data_dict, out, label_column=0, score_column=2):
    fig = plt.figure(0)
    colors = cycle(['red', 'darkorange', 'aqua', 'cornflowerblue', 'blueviolet', 'deeppink', 'cyan'])
    line_styles = cycle(['-', '--', '-.'])
    for ls, color, key in zip(line_styles, colors, data_dict):
        data = data_dict[key]
        fprIndep, tprIndep, thresholdsIndep = roc_curve(data[:, label_column], data[:, score_column])
        ind_auc = auc(fprIndep, tprIndep)

        plt.plot(fprIndep, tprIndep, lw=2, alpha=0.7, color=color, linestyle=ls,
                 label=r'%s: %.2f ' % (key, ind_auc))
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='b', label='Random', alpha=.8)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig(out)
    plt.close(0)


def plot_prc_CV(data, out, label_column=0, score_column=2):
    precisions = []
    aucs = []
    recall_array = []
    precision_array = []
    mean_recall = np.linspace(0, 1, 100)

    for i in range(len(data)):
        precision, recall, _ = precision_recall_curve(data[i][:, label_column], data[i][:, score_column])
        recall_array.append(recall)
        precision_array.append(precision)
        precisions.append(interp(mean_recall, recall[::-1], precision[::-1])[::-1])
        roc_auc = auc(recall, precision)
        aucs.append(roc_auc)

    colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'blueviolet', 'deeppink', 'cyan'])
    ## ROC plot for CV
    fig = plt.figure(0)
    for i, color in zip(range(len(recall_array)), colors):
        plt.plot(recall_array[i], precision_array[i], lw=1, alpha=0.7, color=color,
                 label='PRC fold %d (AUPRC = %0.2f)' % (i + 1, aucs[i]))
    mean_precision = np.mean(precisions, axis=0)
    mean_recall = mean_recall[::-1]
    mean_auc = auc(mean_recall, mean_precision)
    std_auc = np.std(aucs)

    plt.plot(mean_recall, mean_precision, color='blue',
             label=r'Mean PRC (AUPRC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.9)
    std_precision = np.std(precisions, axis=0)
    precision_upper = np.minimum(mean_precision + std_precision, 1)
    precision_lower = np.maximum(mean_precision - std_precision, 0)
    plt.fill_between(mean_recall, precision_lower, precision_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')
    plt.xlim([0, 1.0])
    plt.ylim([0, 1.0])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc="lower left")
    plt.savefig(out)
    plt.close(0)
    return mean_auc


def plot_mean_prc_CV(data_dict, out, label_column=0, score_column=2):
    fig = plt.figure(0)
    colors = cycle(['red', 'darkorange', 'aqua', 'cornflowerblue', 'blueviolet', 'deeppink', 'cyan'])
    line_styles = cycle(['-', '--', '-.'])
    for ls, color, key in zip(line_styles, colors, data_dict):
        data = data_dict[key]
        precisions = []
        aucs = []
        recall_array = []
        precision_array = []
        mean_recall = np.linspace(0, 1, 100)

        for i in range(len(data)):
            precision, recall, _ = precision_recall_curve(data[i][:, label_column], data[i][:, score_column])
            recall_array.append(recall)
            precision_array.append(precision)
            precisions.append(interp(mean_recall, recall[::-1], precision[::-1])[::-1])
            roc_auc = auc(recall, precision)
            aucs.append(roc_auc)

        mean_precision = np.mean(precisions, axis=0)
        mean_recall = mean_recall[::-1]
        mean_auc = auc(mean_recall, mean_precision)

        plt.plot(mean_recall, mean_precision, color=color, linestyle=ls, label=r'%s: %.2f ' % (key, mean_auc),
                 lw=2, alpha=.9)
    plt.xlim([0, 1.0])
    plt.ylim([0, 1.0])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc="lower left")
    plt.savefig(out)
    plt.close(0)


def plot_prc_ind(data, out, label_column=0, score_column=2):
    precision, recall, _ = precision_recall_curve(data[:, label_column], data[:, score_column])
    ind_auc = auc(recall, precision)
    fig = plt.figure(0)
    plt.plot(recall, precision, lw=2, alpha=0.7, color='red',
             label='PRC curve (area = %0.2f)' % ind_auc)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc="lower left")
    plt.savefig(out)
    plt.close(0)
    return ind_auc

def plot_prc_muti_ind(data_dict, out, label_column=0, score_column=2):
    fig = plt.figure(0)
    colors = cycle(['red', 'darkorange', 'aqua', 'cornflowerblue', 'blueviolet', 'deeppink', 'cyan'])
    line_styles = cycle(['-', '--', '-.'])
    for ls, color, key in zip(line_styles, colors, data_dict):
        data = data_dict[key]
        precision, recall, _ = precision_recall_curve(data[:, label_column], data[:, score_column])
        ind_auc = auc(recall, precision)
        plt.plot(recall, precision, lw=2, alpha=0.7, color=color, linestyle=ls, label=r'%s: %.2f ' % (key, ind_auc))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc="lower left")
    plt.savefig(out)
    plt.close(0)

# def plot_2d(data, labels, file='scatter.png'):
#     data = np.array(data)[:]
#     data = data[:, 1:].astype(float)
#
#     fig = plt.figure(0)
#     if len(labels) == 0:
#         plt.scatter(data[:, 0], data[:, 1], 20, c='r')
#     else:
#         df = pd.DataFrame({'X': data[:, 0], 'Y': data[:, 1], 'L': labels})
#         mySet = set(labels)
#         for l in mySet:
#             newData = df.loc[df.loc[:, "L"] == l, :]
#             plt.scatter(np.array(newData.X), np.array(newData.Y), 20, label="%s" % l)
#             plt.legend(loc='best')
#     plt.savefig(file)
#     plt.close(0)
#     return None


def plot_2d(data, labels, file='scatter.png'):
    data = np.array(data)[:]
    data = data[:, 1:].astype(float)

    color_sets = cycle(['dodgerblue', 'coral', 'limegreen', 'violet', 'mediumslateblue'])
    color_set = []
    label_set = list(set(labels))
    for i, j in zip(label_set, color_sets):
        color_set.append(j)

    my_dict = {}
    for i in range(len(label_set)):
        my_dict[label_set[i]] = color_set[i]

    fig = plt.figure(0)
    if len(labels) == 0:
        plt.scatter(data[:, 0], data[:, 1], 20, c='r')
    else:
        df = pd.DataFrame({'X': data[:, 0], 'Y': data[:, 1], 'L': labels})
        mySet = set(labels)
        for l in mySet:
            newData = df.loc[df.loc[:, "L"] == l, :]
            plt.scatter(np.array(newData.X), np.array(newData.Y), 20, c=my_dict[l], label="%s" % l)
            plt.legend(loc='best')
    plt.xlabel('pc.1')
    plt.ylabel('pc.2')
    plt.savefig(file)
    plt.close(0)
    return None


def plot_3d(data, labels, file='scatter_3d.png'):
    data = np.array(data)[:]
    data = data[:, 1:].astype(float)

    # mark_sets = cycle(['o', '^', '+', ','])
    mark_sets = cycle(['o', 'o'])
    color_sets = cycle(['dodgerblue', 'coral', 'limegreen', 'violet', 'mediumslateblue'])
    label_set = list(set(labels))
    my_dict = {}
    m = 0
    for i in label_set:
        my_dict[i] = m
        m = m + 1

    mark_set = []
    color_set = []
    for i, j, k in zip(label_set, mark_sets, color_sets):
        mark_set.append(j)
        color_set.append(k)

    mc = np.zeros((len(labels), 2)).astype(str)
    for i in range(len(labels)):
        mc[i][0], mc[i][1] = mark_set[my_dict[labels[i]]], color_set[my_dict[labels[i]]]

    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(data)):
        ax.scatter(data[i][0], data[i][1], data[i][2], c=mc[i][1], marker=mc[i][0])
    ax.set_xlabel('pc.1')
    ax.set_ylabel('pc.2')
    ax.set_zlabel('pc.3')
    plt.savefig(file)
    plt.close(0)
    return None

def plot_clustering_2d(encodings, myCluster, output, **kw):
    if myCluster != 0:
        if kw['sof'] == 'sample':
            data = np.array(encodings)[1:, 1:].astype(float)
        else:
            data = np.array(encodings).T[1:, 1:].astype(float)
        labels = np.array(myCluster)[0:, 1:].reshape(-1, )
        e = ''
        try:
            Y = tsne.tsne(data, 2, 50, 20.0)
        except RuntimeWarning as e:
            Y = pca.pca(data, n_components=2)

        df = pd.DataFrame({'X': Y[:, 0], 'Y': Y[:, 1], 'L': labels})

        fig = plt.figure(0)
        mySet = set(labels)
        if len(mySet) > 5:
            plt.scatter(Y[:, 0], Y[:, 1], 20, labels)
        else:
            for l in mySet:
                newData = df.loc[df.loc[:, "L"] == l, :]
                plt.scatter(np.array(newData.X), np.array(newData.Y), 20, label="Cluster_%s" % l)
        plt.legend(loc='best')
        plt.savefig('%s.png' % output)
        plt.close(0)

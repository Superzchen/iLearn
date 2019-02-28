#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import math
import sys


def calculate_metrics(labels, scores, cutoff=0.5, po_label=1):
    my_metrics = {
        'Sensitivity': 'NA',
        'Specificity': 'NA',
        'Accuracy': 'NA',
        'MCC': 'NA',
        'Recall': 'NA',
        'Precision': 'NA',
        'F1-score': 'NA',
        'Cutoff': cutoff,
    }

    tp, tn, fp, fn = 0, 0, 0, 0
    for i in range(len(scores)):
        if labels[i] == po_label:
            if scores[i] >= cutoff:
                tp = tp + 1
            else:
                fn = fn + 1
        else:
            if scores[i] < cutoff:
                tn = tn + 1
            else:
                fp = fp + 1

    my_metrics['Sensitivity'] = tp / (tp + fn) if (tp + fn) != 0 else 'NA'
    my_metrics['Specificity'] = tn / (fp + tn) if (fp + tn) != 0 else 'NA'
    my_metrics['Accuracy'] = (tp + tn) / (tp + fn + tn + fp)
    my_metrics['MCC'] = (tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) if (tp + fp) * (
        tp + fn) * (tn + fp) * (tn + fn) != 0 else 'NA'
    my_metrics['Precision'] = tp / (tp + fp) if (tp + fp) != 0 else 'NA'
    my_metrics['Recall'] = my_metrics['Sensitivity']
    my_metrics['F1-score'] = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) != 0 else 'NA'
    return my_metrics


def calculate_metrics_cv(cv_res, label_column=0, score_column=2, cutoff=0.5, po_label=1):
    metrics_list = []
    for i in cv_res:
        metrics_list.append(calculate_metrics(i[:, label_column], i[:, score_column], cutoff=cutoff, po_label=po_label))
    return metrics_list

# calculate accuracy for each class
def calculate_metrics_cv_muti(cv_res, classes, label_column=0, cutoff=0.5):
    acc_list = []
    for i in cv_res:
        tmp_dict = {}
        for c in range(len(classes)):
            tmp_dict[classes[c]] = \
            calculate_metrics(i[:, label_column], i[:, c + 1], cutoff=cutoff, po_label=classes[c])['Accuracy']
        acc_list.append(tmp_dict)
    return acc_list


def calculate_metrics_ind_muti(ind_res, classes, label_column=0, cutoff=0.5):
    acc_dict = {}
    for c in range(len(classes)):
        acc_dict[classes[c]] = \
        calculate_metrics(ind_res[:, label_column], ind_res[:, c + 1], cutoff=cutoff, po_label=classes[c])['Accuracy']
    return acc_dict

def calculate_acc(res, classes, label_column):
    # column-1: label
    my_dict = {}
    num = 0
    for i in sorted(classes):
        my_dict[num] = int(i)
        num += 1
    max_column = np.argmax(res[:, 1:], axis=1)
    max_label = np.zeros((len(max_column,)))
    for i in range(len(max_column)):
        max_label[i] = my_dict[max_column[i]]

    corrected_pred = 0
    for i in range(len(max_label)):
        if max_label[i] == res[:, label_column][i]:
            corrected_pred += 1
    return corrected_pred / len(max_label)


# calculate accuracy for all class: Accuracy = Num(pred=label) / Num(pred)
def calculate_metrics_cv_muti_V2(cv_res, classes, label_column=0):
    accuracy = np.zeros((len(cv_res)))
    for i in range(len(cv_res)):
        accuracy[i] = calculate_acc(cv_res[i], classes, label_column)
        # print(accuracy[i])
    return np.mean(accuracy)

def calculate_metrics_ind_muti_V2(ind_res, classes, label_column=0):
    return calculate_acc(ind_res, classes, label_column)

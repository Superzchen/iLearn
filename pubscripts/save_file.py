#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import pandas as pd


def write_to_svm(encodings, file):
    with open(file, 'w') as f:
        for line in encodings[1:]:
            line = line[1:]
            f.write('%s' % line[0])
            for i in range(1, len(line)):
                f.write('  %d:%s' % (i, line[i]))
            f.write('\n')


def write_to_tsv(encodings, file):
    with open(file, 'w') as f:
        for line in encodings[1:]:
            line = line[1:]
            f.write('%s' % line[0])
            for i in range(1, len(line)):
                f.write('\t%s' % line[i])
            f.write('\n')


def write_to_tsv_1(encodings, file):
    with open(file, 'w') as f:
        for line in encodings:
            f.write('%s' % line[0])
            for i in range(2, len(line)):
                f.write('\t%s' % line[i])
            f.write('\n')


def write_to_csv(encodings, file):
    with open(file, 'w') as f:
        for line in encodings[1:]:
            line = line[1:]
            f.write('%s' % line[0])
            for i in range(1, len(line)):
                f.write(',%s' % line[i])
            f.write('\n')


def write_to_weka(encodings, file):
    with open('%s.weka' % file, 'w') as f:
        f.write('@relation descriptor\n\n')
        for i in range(1, len(encodings[0][2:]) + 1):
            f.write('@attribute f.%d numeric\n' % i)
        f.write('@attribute play {yes, no}\n\n')
        f.write('@data\n')
        for line in encodings[1:]:
            line = line[1:]
            for fea in line[1:]:
                f.write('%s,' % fea)
            if line[0] == '1':
                f.write('yes\n')
            else:
                f.write('no\n')


def save_file(encodings, format='svm', file='encodings.txt'):
    if encodings == 0:
        with open(file, 'w') as f:
            f.write('An error encountered.')
    else:
        if format == 'svm':
            write_to_svm(encodings, file)
        if format == 'tsv':
            write_to_tsv(encodings, file)
        if format == 'csv':
            write_to_csv(encodings, file)
        if format == 'weka':
            write_to_weka(encodings, file)
        if format == 'tsv_1':
            write_to_tsv_1(encodings, file)


def save_cluster_result(cluster, e, file='Clusters.txt'):
    with open(file, 'w') as f:
        if cluster == 0:
            f.write(str(e))
        else:
            myCluster = np.array(cluster)
            df = pd.DataFrame({'name': myCluster[:, 0], 'cluster': myCluster[:, 1]})
            mySet = set(np.array(df.cluster).tolist())
            f.write('# The sample/feature can be clustered into %d clusters:\n' % len(mySet))
            for l in sorted(mySet):
                newData = np.array(df.loc[df.loc[:, "cluster"] == l, :].name)
                f.write('Cluster_%s:\t' % l)
                for i in newData:
                    f.write(i + '\t')
                f.write('\n')
            f.write('\n==============================================================\n')

            f.write('Protein/Feature\tcluster\n')
            for i in cluster:
                f.write(i[0] + '\t' + str(i[1]) + '\n')
    return None


def save_FS_result(feature, e, method, file='featureRank.txt'):
    with open(file, 'w') as f:
        if feature == 0:
            f.write(str(e))
        else:
            f.write('# Feature selection method: %s\n' % method)
            f.write(
                '# The features were ranked according to their importance, the topper the more important the feature is\n');
            f.write('=======================\n')
            for i in feature:
                f.write(i[0] + '\t' + str(i[1]) + '\n')
    return None


def save_reduction_result(reduced_data, file='dimension_reduction.txt'):
    with open(file, 'w') as f:
        f.write('sample')
        for i in range(1, len(reduced_data[0])):
            f.write('\tpc.' + str(i))
        f.write('\n')
        for i in reduced_data:
            f.write(i[0])
            for j in range(1, len(i)):
                f.write('\t' + str(i[j]))
            f.write('\n')
    return None


def save_CV_result_binary(data, out, info=None):
    with open(out, 'w') as f:
        if info:
            f.write('%s\n' % info)
        for i in range(len(data)):
            f.write('# result for fold %d\n' %(i + 1))
            for j in range(len(data[i])):
                f.write('%d\t%s\n' % (data[i][j][0], data[i][j][2]))
    return None

def save_CV_result(data, classes, out, info=None):
    with open(out, 'w') as f:
        if info:
            f.write('%s\n' % info)
        for i in range(len(data)):
            f.write('result for fold %d\n' %(i + 1))
            f.write('label')
            for k in classes:
                f.write('\t%s' %k)
            f.write('\n')
            for j in range(len(data[i])):
                f.write('%d' % data[i][j][0])
                for k in range(1, len(data[i][j])):
                    f.write('\t%s' %data[i][j][k])
                f.write('\n')
    return None


def save_IND_result_binary(data, out, info=None):
    with open(out, 'w') as f:
        if info:
            f.write('%s\n' % info)
        for i in data:
            f.write('%d\t%s\n' % (i[0], i[2]))
    return None

def save_IND_result(data, classes, out, info=None):
    with open(out, 'w') as f:
        if info:
            f.write('%s\n' % info)
        f.write('label')
        for k in classes:
            f.write('\t%s' % k)
        f.write('\n')
        for i in data:
            f.write('%d' %i[0])
            for j in range(1, len(i)):
                f.write('\t%s' %i[j])
            f.write('\n')
    return None


def save_prediction_metrics_ind(m_dict, out):
    with open(out, 'w') as f:
        f.write('#')
        for key in m_dict:
            f.write('\t%s' %key)
        f.write('\n')
        f.write('Indep')
        for key in m_dict:
            f.write('\t%s' %m_dict[key])
        f.write('\n')
    return None


def save_prediction_metrics_ind_muti(m_dict, classes, out):
    with open(out, 'w') as f:
        f.write('#')
        for c in classes:
            f.write('\tclass_%s_acc' %c)
        f.write('\n')
        f.write('Indep')
        for c in classes:
            f.write('\t%s' % m_dict[c])
        f.write('\n')
    return None


def save_prediction_metrics_cv(m_list, out):
    with open(out, 'w') as f:
        f.write('Fold')
        for key in m_list[0]:
            f.write('\t%s' %key)
        f.write('\n')
        for i in range(len(m_list)):
            f.write('%d' %(i + 1))
            for key in m_list[i]:
                f.write('\t%s' %m_list[i][key])
            f.write('\n')
    return None

def save_prediction_metrics_2_classes(m_dict, out):
    with open(out, 'w') as f:
        f.write('#')
        for m in ['Sensitivity', 'Specificity', 'Accuracy', 'MCC', 'Recall', 'Precision', 'F1-score']:
            f.write('\t%s' %m)
        f.write('\n')
        for key in m_dict:
            f.write('%s' %key)
            for m in ['Sensitivity', 'Specificity', 'Accuracy', 'MCC', 'Recall', 'Precision', 'F1-score']:
                f.write('\t%.3f' %m_dict[key][m])
            f.write('\n')
    return None

def save_prediction_metrics_cv_muti(m_list, classes, out):
    with open(out, 'w') as f:
        f.write('Fold')
        for c in classes:
            f.write('\tclass_%s_acc' %c)
        f.write('\n')
        for fold in range(len(m_list)):
            f.write('fold_%s' %(fold + 1))
            for c in classes:
                f.write('\t%s' %m_list[fold][c])
            f.write('\n')
    return None

def save_prediction_metrics_muti_V2(muti_record_metrics, out):
    with open(out, 'w') as f:
        f.write('#\tCV_accuracy\tIND_accuracy\n')
        for key in muti_record_metrics:
            f.write('%s\t%.4f\t%.4f\n' %(key, muti_record_metrics[key][0], muti_record_metrics[key][1]))

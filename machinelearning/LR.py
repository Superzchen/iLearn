#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import numpy as np
# from pubscripts import read_code_ml, save_file, draw_plot, calculate_prediction_metrics
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression


def LR_Classifier_binary(X, y, indep=None, fold=5, out='LR_output'):
    classes = sorted(list(set(y)))
    prediction_result_cv = []
    indep_out = np.zeros([])
    if indep.shape[0] != 0:
        indep_out = np.zeros((indep.shape[0], len(classes) + 1))
        indep_out[:, 0] = indep[:, 0]

    folds = StratifiedKFold(fold).split(X, y)
    for trained, valided in folds:
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = LogisticRegression(C=1.0, random_state=0).fit(train_X, train_y)
        scores = model.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)
        # independent
        if indep.shape[0] != 0:
            indep_out[:, 1:] += model.predict_proba(indep[:, 1:])
    if indep.shape[0] != 0:
        indep_out[:, 1:] /= fold
    header = ''
    return header, prediction_result_cv, indep_out

'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="training LR model.")
    parser.add_argument("--train", required=True, help="input training coding file")
    parser.add_argument("--indep", help="independent coding file")
    parser.add_argument("--format", choices=['tsv', 'svm', 'csv', 'weka'], default='tsv',
                        help="input file format (default tsv format)")
    parser.add_argument("--fold", type=int, default=5,
                        help="n-fold cross validation mode (default 5-fold cross-validation, 1 means jack-knife cross-validation)")
    parser.add_argument("--out", default="LR_output", help="set prefix for output score file")
    args = parser.parse_args()

    X, y, independent = 0, 0, np.array([])
    X, y = read_code_ml.read_code(args.train, format='%s' % args.format)
    if args.indep:
        ind_X, ind_y = read_code_ml.read_code(args.indep, format='%s' % args.format)
        independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
        independent[:, 0], independent[:, 1:] = ind_y, ind_X

    para_info, cv_res, ind_res = LR_Classifier_binary(X, y, indep=independent, fold=args.fold, out=args.out)
    # save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % args.out, para_info)
    # save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % args.out, para_info)
    # mean_auc = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % args.out, label_column=0, score_column=2)
    # mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % args.out, label_column=0, score_column=2)
    # cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2, )
    # save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % args.out)
    #
    # if args.indep:
    #     ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % args.out, label_column=0, score_column=2)
    #     ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % args.out, label_column=0, score_column=2)
    #     ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
    #     save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' % args.out)
    classes = sorted(list(set(y)))
    if len(classes) == 2:
        save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % args.out, para_info)
        mean_auc = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % args.out, label_column=0, score_column=2)
        mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % args.out, label_column=0, score_column=2)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2, )
        save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % args.out)

        if args.indep:
            save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % args.out, para_info)
            ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % args.out, label_column=0, score_column=2)
            ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % args.out, label_column=0, score_column=2)
            ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
            save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' % args.out)
    if len(classes) > 2:
        save_file.save_CV_result(cv_res, classes, '%s_CV.txt' % args.out, para_info)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv_muti(cv_res, classes, label_column=0)
        save_file.save_prediction_metrics_cv_muti(cv_metrics, classes, '%s_metrics_CV.txt' % args.out)

        if args.indep:
            save_file.save_IND_result(ind_res, classes, '%s_IND.txt' % args.out, para_info)
            ind_metrics = calculate_prediction_metrics.calculate_metrics_ind_muti(ind_res, classes, label_column=0)
            save_file.save_prediction_metrics_ind_muti(ind_metrics, classes, '%s_metrics_IND.txt' % args.out)
'''
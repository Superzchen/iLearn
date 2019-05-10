#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import numpy as np
from pubscripts import read_code_ml, save_file, draw_plot, calculate_prediction_metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold


def RF_Classifier(X, y, indep=None, fold=5, n_trees=100, out='RF_output'):
    """
    Parameters:
    ----------
    :param X: 2-D ndarray
    :param y: 1-D ndarray
    :param indep: 2-D ndarray, the first column is labels and the rest are feature values
    :param fold: int, default 5
    :param n_trees: int, number of trees, default: 5
    :param out:
    :return:
        info: str, the model parameters
        cross-validation result: list with element is ndarray
        independent result: ndarray, the first column is labels and the rest are prediction scores.
    """
    classes = sorted(list(set(y)))
    if indep.shape[0] != 0:
        indep_out = np.zeros((indep.shape[0], len(classes) + 1))
        indep_out[:, 0] = indep[:, 0]

    prediction_result_cv = []
    prediction_result_ind = np.array([])
    if indep.shape[0] != 0:
        prediction_result_ind = np.zeros((len(indep), len(classes) + 1))
        prediction_result_ind[:, 0] = indep[:, 0]

    folds = StratifiedKFold(fold).split(X, y)
    for i, (trained, valided) in enumerate(folds):
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = RandomForestClassifier(n_estimators=n_trees, bootstrap=False)
        rfc = model.fit(train_X, train_y)
        scores = rfc.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)
        # independent
        if indep.shape[0] != 0:
            prediction_result_ind[:, 1:] += rfc.predict_proba(indep[:, 1:])
    if indep.shape[0] != 0:
        prediction_result_ind[:, 1:] /= fold
    header = 'n_trees: %d' % n_trees
    return header, prediction_result_cv, prediction_result_ind


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.", description="training RF model.")
    parser.add_argument("--train", required=True, help="input training coding file")
    parser.add_argument("--indep", help="independent coding file")
    parser.add_argument("--format", choices=['tsv', 'svm', 'csv', 'weka'], default='tsv',
                        help="input file format (default tab format)")
    parser.add_argument("--n_trees", type=int, default=100, help="the number of trees in the forest (default 100)")
    parser.add_argument("--fold", type=int, default=5,
                        help="n-fold cross validation mode (default 5-fold cross-validation, 1 means jack-knife cross-validation)")
    parser.add_argument("--out", default="RF_output", help="set prefix for output score file")
    args = parser.parse_args()

    X, y, independent = 0, 0, np.array([])
    X, y = read_code_ml.read_code(args.train, format='%s' % args.format)
    if args.indep:
        ind_X, ind_y = read_code_ml.read_code(args.indep, format='%s' % args.format)
        independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
        independent[:, 0], independent[:, 1:] = ind_y, ind_X

    para_info, cv_res, ind_res = RF_Classifier(X, y, indep=independent, fold=args.fold, n_trees=args.n_trees,
                                               out=args.out)
    classes = sorted(list(set(y)))
    if len(classes) == 2:
        save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % args.out, para_info)
        mean_auc = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % args.out, label_column=0, score_column=2)
        mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % args.out, label_column=0, score_column=2)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2,)
        save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % args.out)

        if args.indep:
            save_file.save_IND_result_binary(ind_res,'%s_IND.txt' % args.out, para_info)
            ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % args.out, label_column=0, score_column=2)
            ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % args.out, label_column=0, score_column=2)
            ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
            save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' %args.out)
    if len(classes) > 2:
        save_file.save_CV_result(cv_res, classes, '%s_CV.txt' % args.out, para_info)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv_muti(cv_res, classes, label_column=0)
        save_file.save_prediction_metrics_cv_muti(cv_metrics, classes, '%s_metrics_CV.txt' % args.out)

        if args.indep:
            save_file.save_IND_result(ind_res, classes, '%s_IND.txt' % args.out, para_info)
            ind_metrics = calculate_prediction_metrics.calculate_metrics_ind_muti(ind_res, classes, label_column=0)
            save_file.save_prediction_metrics_ind_muti(ind_metrics, classes, '%s_metrics_IND.txt' % args.out)


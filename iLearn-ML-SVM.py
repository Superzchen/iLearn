#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import math
import numpy as np
from pubscripts import read_code_ml, save_file, draw_plot, calculate_prediction_metrics
from sklearn import svm
from sklearn.model_selection import StratifiedKFold, GridSearchCV


def SVM_Classifier(X, y, indep=None, fold=5, batch=None, auto=False, kernel='rbf', degree=3, gamma='auto',
                          coef0=0, C=1.0):
    default_params = {'degree': degree, 'gamma': gamma, 'coef0': coef0, 'C': C}
    if auto:
        data = np.zeros((X.shape[0], X.shape[1] + 1))
        data[:, 0] = y
        data[:, 1:] = X
        np.random.shuffle(data)
        X1 = data[:, 1:]
        y1 = data[:, 0]
        parameters = {'kernel': ['linear'], 'C': [1, 15]} if kernel == 'linear' else {'kernel': [kernel],
                                                                                      'C': [1, 15],
                                                                                      'gamma': 2.0 ** np.arange(-10, 4)}
        optimizer = GridSearchCV(svm.SVC(probability=True), parameters)
        optimizer = optimizer.fit(X1[0:math.ceil(batch * X1.shape[0]), ],
                                  y1[0:math.ceil(batch * y1.shape[0]), ]) if batch else optimizer.fit(X, y)
        params = optimizer.best_params_
        default_params['C'] = params['C']
        if kernel != 'linear':
            default_params['gamma'] = params['gamma']

    classes = sorted(list(set(y)))
    svms = []
    cvs = np.zeros((X.shape[0], len(classes) + 1))
    folds = StratifiedKFold(fold).split(X, y)

    inds = np.array([])
    if indep.shape[0] != 0:
        inds = np.zeros((len(indep), len(classes) + 1))
        inds[:, 0] = indep[:, 0]

    prediction_result_cv = []
    for trained, valided in folds:
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = svm.SVC(C=default_params['C'], kernel=kernel, degree=default_params['degree'],
                        gamma=default_params['gamma'], coef0=default_params['coef0'], probability=True,
                        random_state=1)
        svc = model.fit(train_X, train_y)
        svms.append(svc)
        proba_ = svc.predict_proba(valid_X)
        cvs[valided, 0], cvs[valided, 1:] = valid_y, proba_
        # save the sample label and prediction result to ndarray
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, proba_
        prediction_result_cv.append(tmp_result)

        # independent
        if indep.shape[0] != 0:
            inds[:, 1:] += svc.predict_proba(indep[:, 1:])

    header = 'C=%f\tgamma=%s' % (default_params['C'], default_params['gamma'])
    if indep.shape[0] != 0:
        inds[:, 1:] /= fold
    return header, prediction_result_cv, inds

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Training SVM model.")
    parser.add_argument("--train", required=True, help="input training coding file")
    parser.add_argument("--indep", help="independent coding file")
    parser.add_argument("--format", choices=['tsv', 'svm', 'csv', 'weka'], default='tsv',
                        help="input file format (default tab format)")
    parser.add_argument("--kernel", choices=['linear', 'poly', 'rbf', 'sigmoid'], default='rbf',
                        help="SVM kernel type (default rbf kernel)")
    parser.add_argument("--auto", action='store_true', help="auto optimize parameters (defalult: False)")
    parser.add_argument("--batch", type=float, default=1,
                        help="random select part (batch * samples) samples for parameters optimization")
    parser.add_argument("--degree", type=int, default=3, help="set degree in polynomial kernel function (default 3)")
    parser.add_argument("--gamma", type=float, help="set gamma in polynomial/rbf/sigmoid kernel function (default 1/k)")
    parser.add_argument("--coef0", type=float, default=0,
                        help="set coef0 in polynomial/rbf/sigmoid kernel function (default 0)")
    parser.add_argument("--cost", type=float, default=1, help="set the parameter cost value (default 1)")
    parser.add_argument("--fold", type=int, default=5,
                        help="n-fold cross validation mode (default 5-fold cross-validation, 1 means jack-knife cross-validation)")
    parser.add_argument("--out", default="SVM_output", help="set prefix for output score file")
    args = parser.parse_args()
    if args.gamma == None:
        args.gamma = 'auto'

    X, y, independent = 0, 0, np.array([])
    X, y = read_code_ml.read_code(args.train, format='%s' % args.format)
    if args.indep:
        ind_X, ind_y = read_code_ml.read_code(args.indep, format='%s' % args.format)
        independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
        independent[:, 0], independent[:, 1:] = ind_y, ind_X

    para_info, cv_res, ind_res = SVM_Classifier(X, y, indep=independent, fold=args.fold, batch=args.batch,
                                                       auto=args.auto, kernel=args.kernel, degree=args.degree,
                                                       gamma=args.gamma, coef0=args.coef0, C=args.cost)


    classes = sorted(list(set(y)))
    if len(classes) == 2:
        save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % args.out, para_info)
        mean_auc = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % args.out, label_column=0, score_column=2)
        mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % args.out, label_column=0, score_column=2)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2,)
        save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % args.out)

        if args.indep:
            save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % args.out, para_info)
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




#!/usr/bin/env python

import argparse
import pickle
import numpy as np
import pandas as pd
from sklearn import metrics, tree, ensemble, model_selection
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='ensemble classifier fitting (training data) or application (test data)')
parser.add_argument('-m', '--mode', required=True, help='train | predict')
parser.add_argument('-ft1', '--filetrainpos', required=False, default='featspace_train_pos.txt',
                                                help='feature space representation of positive training examples')
parser.add_argument('-ft2', '--filetrainneg', required=False, default='featspace_train_neg.txt',
                                                help='feature space representation of negative training examples')
parser.add_argument('-fp', '--filepredict', required=False, default='featspace_pred.txt',
                                                help='feature space representation of test set sequences')
parser.add_argument('-clf', '--classifier', required=False, default='../Out/classifier.pkl',
                                                help='trained classifier serialized as a pickle file')
parser.add_argument('-n1', '--Nseqpos', required=False, default=1000, help='use N of the 10k positive examples for training')
parser.add_argument('-n2', '--Nseqneg', required=False, default=10000, help='use N of the 10k negative examples for training')
args = parser.parse_args()

if args.mode == 'predict':
    with open('../Output/classifier.pkl', 'rb') as fh:
        clf = pickle.load(fh)
    featspace_pred = pd.read_csv(args.filepredict, header=0)
    X_pred = np.array(featspace_pred)
    y_pred = clf.predict_proba(X_pred)[:,1]
    print(y_pred[:100])

elif args.mode == 'train':
    topN_pos = int(args.Nseqpos)
    topN_neg = int(args.Nseqneg)
    featspace_train_pos = pd.read_csv(args.filetrainpos, header=0)
    featspace_train_neg = pd.read_csv(args.filetrainneg, header=0)

    X_train = np.concatenate((np.array(featspace_train_pos[:topN_pos]),
                              np.array(featspace_train_neg[-topN_neg:])), axis=0)
    y_train = np.array([1]*topN_pos + [0]*topN_neg)
    print(np.shape(X_train), np.shape(y_train))

    numfolds = 10
    skf = model_selection.StratifiedKFold(n_splits=numfolds)
    skf.get_n_splits(X_train, y_train)
    GB_list = []
    scores = np.zeros((numfolds,))
    boosting_rounds = 100
    depth = 3
    LR = 0.1
    train_losses = np.zeros((numfolds, boosting_rounds))
    test_losses = np.zeros((numfolds, boosting_rounds))
    feature_scores = np.zeros((numfolds, X_train.shape[1]))
    counter = 0
    for train_index, test_index in skf.split(X_train, y_train):
        trainX, trainY = X_train[train_index], y_train[train_index]
        testX, testY = X_train[test_index], y_train[test_index]
        GB_list.append(ensemble.GradientBoostingClassifier(n_estimators=boosting_rounds, learning_rate=LR, max_depth=depth))
        GB_list[counter].fit(trainX, trainY)
        train_losses[counter,:] = GB_list[counter].train_score_
        feature_scores[counter,:] = GB_list[counter].feature_importances_
        for i, predY in enumerate(GB_list[counter].staged_decision_function(testX)):
            test_losses[counter,i] = GB_list[counter].loss_(testY, predY)
        predY = GB_list[counter].predict_proba(testX)[:,1]
        fpr, tpr, thresholds = metrics.roc_curve(testY, predY, pos_label=1)
        roc_auc = metrics.auc(fpr, tpr)
        scores[counter] = roc_auc
        counter += 1
        print(counter)

    clf = ensemble.GradientBoostingClassifier(max_depth=depth, n_estimators=boosting_rounds, learning_rate=LR)
    clf.fit(X_train, y_train)
    with open('../Output/classifier.pkl', 'wb') as fh:
        pickle.dump(clf, fh)

    print(scores)
    print(scores.mean(), scores.std())

    x = np.arange(boosting_rounds) + 1
    train_mean = np.mean(train_losses, axis=0)
    train_std = np.std(train_losses, axis=0)
    train_upper = train_mean + train_std
    train_lower = train_mean - train_std
    test_mean = np.mean(test_losses, axis=0)
    test_std = np.std(test_losses, axis=0)
    test_upper = test_mean + test_std
    test_lower = test_mean - test_std

    fig1 = plt.figure()
    plt.plot(x, train_mean, color='blue', linewidth=2.0, label='train split')
    plt.plot(x, test_mean, color='red', linewidth=2.0, label='test split')
    plt.fill_between(x, train_upper, train_lower, color='blue', alpha='0.1')
    plt.fill_between(x, test_upper, test_lower, color='red', alpha='0.1')
    fig1.gca().set_xlim([0, boosting_rounds + 1])
    plt.legend(loc='upper right')
    plt.xlabel('Boosting Iterations')
    plt.ylabel('Binomial Deviance')
    plt.grid('off')
    plt.tight_layout()
    plt.savefig('../Output/boosting_loss_func.pdf')

    feat_mean = np.mean(feature_scores, axis=0)
    feat_std = np.std(feature_scores, axis=0)

    fig2 = plt.figure()
    sorted_idx = np.argsort(feat_mean)[-10:]
    pos = np.arange(sorted_idx.shape[0]) + .5
    plt.barh(pos, feat_mean[sorted_idx], xerr=feat_std[sorted_idx], color='#a7d3ff', ecolor='#23337c', align='center')
    plt.yticks(pos, featspace_train_pos.columns[:][sorted_idx])
    plt.xlabel('Relative Feature Importance')
    plt.grid('off')
    plt.tight_layout()
    plt.savefig('../Output/boosting_importances.pdf')

else:
    print('-m flag must be either train | predict')

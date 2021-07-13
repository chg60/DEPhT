import numpy as np
import pandas as pd

from sklearn.model_selection import KFold
from sklearn.naive_bayes import GaussianNB

from prophicient.classes.prophage_classifier import ProphageClassifier
from prophicient.functions.statistics import average, mcc


def _score(real_classes, predict_classes):
    """
    Scores predicted classes against known classes to count how many
    tps, tns, fps, fns there were.

    :param real_classes: known classes
    :type real_classes: list
    :param predict_classes: predicted classes
    :type predict_classes: list
    :return: tps, tns, fps, fns
    """
    tps, tns, fps, fns = 0, 0, 0, 0

    for real_value, predict_value in zip(real_classes, predict_classes):
        if real_value and predict_value:
            tps += 1
        elif real_value and not predict_value:
            fns += 1
        elif predict_value and not real_value:
            fps += 1
        else:
            tns += 1

    return tps, tns, fps, fns


def _mcc_score(real_classes, predict_classes):
    """
    Helper function to score class predictions and then return
    the mcc.
    """
    tps, tns, fps, fns = _score(real_classes, predict_classes)
    print(tps, tns, fps, fns)
    return mcc(tps, fns, tns, fps)


def train_prophage_classifier(prophage_data, bacteria_data):
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    prophage_splits = [x for x in kf.split(prophage_data)]
    bacteria_splits = [x for x in kf.split(bacteria_data)]

    p_trains = [x[0] for x in prophage_splits]
    p_tests = [x[1] for x in prophage_splits]
    b_trains = [x[0] for x in bacteria_splits]
    b_tests = [x[1] for x in bacteria_splits]

    mcc_scores = list()

    clf = ProphageClassifier()

    zipper = zip(p_trains, p_tests, b_trains, b_tests)
    for p_train, p_test, b_train, b_test in zipper:
        train_feats = pd.concat((prophage_data.iloc[list(p_train), :-1],
                                 bacteria_data.iloc[list(b_train), :-1]),
                                axis=0)
        train_labels = pd.concat((prophage_data.iloc[list(p_train), -1],
                                  bacteria_data.iloc[list(b_train), -1]),
                                 axis=0)
        test_feats = pd.concat((prophage_data.iloc[list(p_test), :-1],
                                bacteria_data.iloc[list(b_test), :-1]), axis=0)
        test_labels = pd.concat((prophage_data.iloc[list(p_test), -1],
                                 bacteria_data.iloc[list(b_test), -1]), axis=0)

        clf.fit(train_feats, train_labels)

        predictions = clf.predict(test_feats)
        mcc_scores.append(_mcc_score(test_labels, predictions))

    mean_mcc = average(mcc_scores)

    all_feats = pd.concat((prophage_data.iloc[:, :-1],
                           bacteria_data.iloc[:, :-1]), axis=0)
    all_labels = pd.concat((prophage_data.iloc[:, -1],
                            bacteria_data.iloc[:, -1]), axis=0)

    clf.fit(all_feats, all_labels)

    return clf, mean_mcc


def train_bayes_classifier(prophage_data, bacteria_data):
    """
    Train a Naive Bayes classifier to descriminate between prophage
    and bacterial genes.

    :param prophage_data: prophage gene features
    :type prophage_data: pandas.DataFrame
    :param bacteria_data: bacterial gene features
    :type bacteria_data: pandas.DataFrame
    :return: best_model
    """
    # Cast pandas DataFrames to numpy Arrays
    prophage_data = np.array(prophage_data)
    bacteria_data = np.array(bacteria_data)

    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    prophage_splits = [x for x in kf.split(prophage_data)]
    bacteria_splits = [x for x in kf.split(bacteria_data)]

    p_trains = [x[0] for x in prophage_splits]
    p_tests = [x[1] for x in prophage_splits]
    b_trains = [x[0] for x in bacteria_splits]
    b_tests = [x[1] for x in bacteria_splits]

    mcc_scores = list()

    clf = GaussianNB()

    zipper = zip(p_trains, p_tests, b_trains, b_tests)
    for p_train, p_test, b_train, b_test in zipper:
        train_feats = np.concatenate((prophage_data[p_train, :-1],
                                      bacteria_data[b_train, :-1]), axis=0)
        train_labels = np.concatenate((prophage_data[p_train, -1],
                                       bacteria_data[b_train, -1]), axis=0)
        test_feats = np.concatenate((prophage_data[p_test, :-1],
                                     bacteria_data[b_test, :-1]), axis=0)
        test_labels = np.concatenate((prophage_data[p_test, -1],
                                      bacteria_data[b_test, -1]), axis=0)

        ratio = float(len(b_train))/len(p_train)
        sample_weights = [ratio] * len(p_train)
        sample_weights.extend([1] * len(b_train))

        clf.fit(train_feats, train_labels, sample_weight=sample_weights)

        predictions = clf.predict(test_feats)
        mcc_scores.append(_mcc_score(test_labels, predictions))

    mean_mcc = average(mcc_scores)

    all_feats = np.concatenate((prophage_data[:, :-1],
                                bacteria_data[:, :-1]), axis=0)
    all_labels = np.concatenate((prophage_data[:, -1],
                                 bacteria_data[:, -1]), axis=0)

    clf.fit(all_feats, all_labels)

    return clf, mean_mcc

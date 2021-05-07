import numpy as np

from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier

from prophicient.functions.statistics import average, \
    matthews_correlation_coefficient

# Seed RNG: the answer to the question of life, the universe, and everything
RANDOM_STATE = 42


def _score(real_classes, predict_classes):
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
    tps, tns, fps, fns = _score(real_classes, predict_classes)
    return matthews_correlation_coefficient(tps, fns, tns, fps)


def train_random_forest_classifier(prophage_data, bacteria_data):
    """
    Train a random forest classifier, to achieve maximal training
    accuracy in discriminating between prophage genes (`prophage_data`)
    and bacterial genes (`bacteria_data`).

    The model with the highest average accuracy (across training
    datasets) is re-trained against all the data, then returned
    as `best_model`.

    :param prophage_data: prophage gene features
    :type prophage_data: pandas.DataFrame
    :param bacteria_data: bacterial gene features
    :type bacteria_data: pandas.DataFrame
    :return: best_model
    """
    # Cast pandas DataFrames to numpy Arrays
    prophage_data = np.array(prophage_data)
    bacteria_data = np.array(bacteria_data)

    kf = KFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)

    prophage_splits = [x for x in kf.split(prophage_data)]
    bacteria_splits = [x for x in kf.split(bacteria_data)]

    p_trains = [x[0] for x in prophage_splits]
    p_tests = [x[1] for x in prophage_splits]
    b_trains = [x[0] for x in bacteria_splits]
    b_tests = [x[1] for x in bacteria_splits]

    mcc_scores = list()
    clf = RandomForestClassifier(n_estimators=250, criterion="gini",
                                 bootstrap=False, n_jobs=-1, class_weight="balanced")

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

        clf.fit(train_feats, train_labels)
        predictions = clf.predict(test_feats)
        mcc_scores.append(_mcc_score(test_labels, predictions))
    mean_mcc = average(mcc_scores)

    # Now we've got the best model parameters, let's re-train on all data
    print(f"RandomForest classifier got average MCC = {mean_mcc:.3f}")

    all_feats = np.concatenate((prophage_data[:, :-1],
                                bacteria_data[:, :-1]), axis=0)
    all_labels = np.concatenate((prophage_data[:, -1],
                                 bacteria_data[:, -1]), axis=0)
    clf.fit(all_feats, all_labels)
    return clf

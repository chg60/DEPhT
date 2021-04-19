import itertools

import numpy as np

from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

from prophicient.functions.statistics import average

# Seed RNG: the answer to the question of life, the universe, and everything
RANDOM_STATE = 42

# Default value ranges for model training
N_ESTIMATORS = (100, 200, 250)
CRITERION = ("entropy", "gini")
MAX_FEATURES = ("log2", "sqrt", None)
CLASS_WEIGHT = ("balanced", None)
PARAMS = itertools.product(N_ESTIMATORS, CRITERION, MAX_FEATURES, CLASS_WEIGHT)


def train_random_forest_classifier(prophage_data, bacteria_data):
    """
    Trains a random forest classifier

    :param prophage_data:
    :param bacteria_data:
    :return:
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

    best_accuracy, best_param, best_model = 0, None, None
    for i, params in enumerate(PARAMS):
        n_estimators, criterion, max_features, class_weight = params
        accuracies = list()
        clf = RandomForestClassifier(n_estimators=n_estimators,
                                     criterion=criterion,
                                     max_features=max_features,
                                     n_jobs=-1,
                                     random_state=RANDOM_STATE,
                                     oob_score=True,
                                     class_weight=class_weight)

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
            accuracies.append(accuracy_score(test_labels, predictions))
        mean_accuracy = average(accuracies)
        if mean_accuracy > best_accuracy:
            # print(f"Local winner ({mean_accuracy*100:.3f}%)")
            best_accuracy = mean_accuracy
            best_param = params
            best_model = clf

    # Now we've got the best model parameters, let's re-train on all data
    estimators, criterion, max_features, class_weight = best_param
    print(f"Best model achieves {best_accuracy * 100:.3f}% training accuracy.")
    all_feats = np.concatenate((prophage_data[:, :-1],
                                bacteria_data[:, :-1]), axis=0)
    all_labels = np.concatenate((prophage_data[:, -1],
                                 bacteria_data[:, -1]), axis=0)
    best_model.fit(all_feats, all_labels)
    return best_model

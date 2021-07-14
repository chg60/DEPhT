import random
import pandas as pd

from prophicient.classes.prophage_classifier import ProphageClassifier
from prophicient.functions.statistics import average, mcc


def mcc_score(real_classes, predict_classes):
    """
    Helper function to score class predictions and then return
    the mcc.
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

    return mcc(tps, fns, tns, fps)


class KFold:
    """
    Slim implementation of sklearn.model_selection.KFold functionality.
    """
    def __init__(self, n_splits=5, shuffle=False, random_state=None):
        if random_state and not shuffle:
            raise ValueError("random_state cannot be used without shuffle")

        if random_state:
            random.seed(random_state)

        self.n_splits = n_splits
        self.shuffle = shuffle

    def split(self, x):
        """
        Generator that yields the train/test indices for the indicated
        n_splits.

        :param x: the length of the data to be split
        :type x: int
        """
        indices = list(range(x))

        # Shuffle indices if told to do so
        if self.shuffle:
            random.shuffle(indices)

        # Now start creating the splits
        for i in range(self.n_splits):
            train, test = list(), list()
            for x, index in enumerate(indices):
                if x % self.n_splits == i:
                    test.append(index)
                else:
                    train.append(index)
            yield sorted(train), sorted(test)


def train_prophage_classifier(prophage_data, bacteria_data):
    """
    Trains a ProphageClassifier on the given prophage_data and
    bacteria_data. This is effectively a simple implementation of
    a Naive Bayes
    """
    # Seed RNG for reproducibility - use the "Answer to the Ultimate
    # Question of Life, the Universe, and Everything." (Douglas Adams)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    # Get train/test split indices for both groups of genes
    prophage_splits = [x for x in kf.split(len(prophage_data))]
    bacteria_splits = [x for x in kf.split(len(bacteria_data))]

    mcc_scores = list()
    clf = ProphageClassifier()

    for prophage_split, bacteria_split in zip(prophage_splits, bacteria_splits):
        p_train, p_test = prophage_split
        b_train, b_test = bacteria_split

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
        mcc_scores.append(mcc_score(test_labels, predictions))

    mean_mcc = average(mcc_scores)

    all_feats = pd.concat((prophage_data.iloc[:, :-1],
                           bacteria_data.iloc[:, :-1]), axis=0)
    all_labels = pd.concat((prophage_data.iloc[:, -1],
                            bacteria_data.iloc[:, -1]), axis=0)

    clf.fit(all_feats, all_labels)

    return clf, mean_mcc

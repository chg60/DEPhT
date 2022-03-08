"""Utility script for training a PhageClassifier."""

import argparse
import pathlib

import pandas as pd

from depht.classes.prophage_classifier import ProphageClassifier
from depht_train.classes.kfold import KFold


def parse_args():
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("phage_csv", type=pathlib.Path,
                   help="path to CSV file containing phage data")
    p.add_argument("bacteria_csv", type=pathlib.Path,
                   help="path to CSV file containing bacterial data")
    p.add_argument("classifier_path", type=pathlib.Path,
                   help="path where classifier should be stored")

    p.add_argument("-k", "--k-fold", type=int, default=5,
                   help=f"perform [k]-fold cross validation [default: {5}]")
    p.add_argument("-w", "--optimize-weights", action="store_true",
                   help="optimize feature weights for the trained model")
    p.add_argument("-t", "--optimize-threshold", action="store_true",
                   help="optimize threshold for the model to predict 0/1")

    return p.parse_args()


def even_weights(n, k=100):
    """Returns the feature weights that weigh `n` features evenly.

    :param n: the number of features to weigh evenly
    :type n: int
    :param k: value the weights should sum to
    :type k: int
    """
    return [round(float(k) / n, 2)] * n


def multiset(n, k=100, k_step=5):
    """Returns the multiset of `n` integers that sum to `k` given a
    step size of `k_step`.

    Useful for generating feature weight parameters to test, given a
    modest `n`, `k`=100 and `k_step` an even divisor of 100 (e.g. 1,
    2, 5, 10, 20, 25, 50).

    :param n: the number of features being trained against
    :type n: int
    :param k: sum of feature weights must equal this number
    :type k: int
    :param k_step: granularity for changing feature weights
    :type k_step: int
    """
    # Neither n, nor k can be negative
    if n < 0 or k < 0:
        raise ValueError(f"map_feature_weights requires (k, n) >= (0, 0): "
                         f"({n}, {k}) invalid")

    if not n:
        return []

    if not k:
        return [[0] * n]

    if n == 1:
        return [[k]]

    return [[0] + v for v in multiset(n - 1, k)] + \
           [[v[0] + k_step] + v[1:] for v in multiset(n, k - k_step)]


def mcc(true_pos, false_neg, true_neg, false_pos):
    """Calculate the Matthews correlation coefficient.

    NOTE: this can be thought of as a balanced accuracy metric
    that isn't skewed by differences in class size.

    :param true_pos: number of correctly identified positives
    :type true_pos: int
    :param false_neg: number of incorrectly identified negatives
    :type false_neg: int
    :param true_neg: number of correctly identified negatives
    :type true_neg: int
    :param false_pos: number of incorrectly identified positives
    :type false_pos: int
    :return: mcc
    """
    numer = (true_pos * true_neg) + (false_pos * false_neg)
    denom = ((true_pos + false_pos) * (true_pos + false_neg) *
             (true_neg + false_pos) * (true_neg + false_neg)) ** 0.5

    return float(numer)/denom


def score(real_labels, predict_labels):
    """
    Helper function to score class predictions and then return
    the mcc and f1_scores.
    """
    tps, tns, fps, fns = 0, 0, 0, 0

    for real_value, predict_value in zip(real_labels, predict_labels):
        if real_value and predict_value:
            tps += 1
        elif real_value and not predict_value:
            fns += 1
        elif predict_value and not real_value:
            fps += 1
        else:
            tns += 1

    print(f"tps: {tps}\tfns: {fns}\ttns: {tns}\tfps: {fps}")

    return mcc(tps, fns, tns, fps)


def train_classifier(df, k=5, verbose=False):
    """Train a ProphageClassifier.

    :param df: concatenated phage genomes dataframe
    :type df: pandas.DataFrame
    :param k: fold validation to perform
    :type k: int
    :param verbose:
    :type verbose: bool
    :return: classifier, feature_weights, threshold
    """
    phg_df, bct_df = df[df["class"] == 1], df[df["class"] == 0]

    clf = ProphageClassifier()

    # Seed RNG for reproducibility - use the "Answer to the Ultimate
    # Question of Life, the Universe, and Everything." (Douglas Adams)
    kf = KFold(n_splits=k, shuffle=True, random_state=42)

    phg_splits, bct_splits = kf.split(len(phg_df)), kf.split(len(bct_df))

    evaluations = dict()
    for phg_split, bct_split in zip(phg_splits, bct_splits):
        p_train, p_test = phg_split
        b_train, b_test = bct_split

        train_feats = pd.concat((phg_df.iloc[list(p_train), :-1],
                                 bct_df.iloc[list(b_train), :-1]), axis=0)
        train_labels = pd.concat((phg_df.iloc[list(p_train), -1],
                                  bct_df.iloc[list(b_train), -1]), axis=0)
        test_feats = pd.concat((phg_df.iloc[list(p_test), :-1],
                                bct_df.iloc[list(b_test), :-1]), axis=0)
        test_labels = pd.concat((phg_df.iloc[list(p_test), -1],
                                 bct_df.iloc[list(b_test), -1]), axis=0)

        clf.fit(train_feats, train_labels)

        # Optimize feature weights
        for fw in multiset(n=len(clf), k=100, k_step=5):
            fw = [x / 100 for x in fw]
            predictions = clf.predict(test_feats, feature_weights=fw, alpha=0.5)
            mcc_score = score(test_labels, predictions)

            key = tuple(fw)
            if key in evaluations:
                evaluations[key].append(mcc_score)
            else:
                evaluations[key] = [mcc_score]

            if verbose:
                print(f"{key}: {mcc_score}")

    ranked_params = list()
    for key, values in evaluations.items():
        value = sum(values) / k
        ranked_params.append((value, key))
    ranked_params = sorted(ranked_params)

    best_params = ranked_params[-1]
    opt_mcc = best_params[0]
    opt_weights = best_params[1]

    if verbose:
        print(f"The best weights ({opt_weights}) get average MCC = "
              f"{round(opt_mcc, 5)}")

    all_feats = pd.concat((phg_df.iloc[:, :-1],
                           bct_df.iloc[:, :-1]), axis=0)
    all_labels = pd.concat((phg_df.iloc[:, -1],
                            bct_df.iloc[:, -1]), axis=0)

    figs = clf.fit(all_feats, all_labels, plot=True)

    return clf, figs

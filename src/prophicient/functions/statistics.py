"""
Some functions useful for statistics - evaluating models, program
performance, transforming data, etc.
"""
import math


def average(values, kind="arithmetic"):
    if kind == "arithmetic":
        numerator, denominator = sum(values), len(values)
    elif kind == "geometric":
        numerator, denominator = math.prod(values), len(values)
    elif kind == "harmonic":
        numerator, denominator = len(values), sum([1.0/x for x in values])
    else:
        raise ValueError(f"'{kind}' is not a supported kind of average")
    return float(numerator)/denominator


def stdev(values, mean=None):
    if not mean:
        mean = average(values)
    var = [(x - mean) ** 2 for x in values]
    return average(var) ** 0.5


def true_positive_rate(true_pos, false_neg):
    """
    Calculate the true positive rate.

    NOTE: same as sensitivity, recall, and hit rate

    :param true_pos: number of correctly identified positives
    :type true_pos: int
    :param false_neg: number of missed positives
    :type false_neg: int
    :return: tpr
    """
    return float(true_pos)/(true_pos + false_neg)


def false_negative_rate(true_pos, false_neg):
    """
    Calculate the false negative rate.

    NOTE: same as miss rate

    :param true_pos: number of correctly identified positives
    :type true_pos: int
    :param false_neg: number of missed positives
    :type false_neg: int
    :return: fnr
    """
    return 1 - true_positive_rate(true_pos, false_neg)


def true_negative_rate(true_neg, false_pos):
    """
    Calculate the true negative rate.

    NOTE: same as specificity, selectivity

    :param true_neg: number of correctly identified negatives
    :type true_neg: int
    :param false_pos: number of incorrectly identified positives
    :type false_pos: int
    :return: tnr
    """
    return float(true_neg)/(true_neg + false_pos)


def false_positive_rate(true_neg, false_pos):
    """
    Calculate the false positive rate.

    NOTE: same as fall-out

    :param true_neg: number of correctly identified negatives
    :type true_neg: int
    :param false_pos: number of incorrectly identified positives
    :type false_pos: int
    :return: fpr
    """
    return 1 - true_negative_rate(true_neg, false_pos)


def positive_predictive_value(true_pos, false_pos):
    """
    Calculate the positive predictive value.

    NOTE: same as precision

    :param true_pos: number of correctly identified positives
    :type true_pos: int
    :param false_pos: number of incorrectly identified positives
    :type false_pos: int
    :return: ppv
    """
    return float(true_pos)/(true_pos + false_pos)


def false_discovery_rate(true_pos, false_pos):
    """
    Calculate the false discovery rate.

    :param true_pos: number of correctly identified positives
    :type true_pos: int
    :param false_pos: number of incorrectly identified positives
    :type false_pos: int
    :return: fdr
    """
    return 1 - positive_predictive_value(true_pos, false_pos)


def negative_predictive_value(true_neg, false_neg):
    """
    Calculate the negative predictive value.

    :param true_neg: number of correctly identified negatives
    :type true_neg: int
    :param false_neg: number of incorrectly identified negatives
    :type false_neg: int
    :return: npv
    """
    return float(true_neg)/(true_neg + false_neg)


def false_omission_rate(true_neg, false_neg):
    """
    Calculate the false omission rate.

    :param true_neg: number of correctly identified negatives
    :type true_neg: int
    :param false_neg: number of incorrectly identified negatives
    :type false_neg: int
    :return: for
    """
    return 1 - negative_predictive_value(true_neg, false_neg)


def f1_score(true_pos, false_pos, false_neg):
    """
    Calculates the F1 score.

    NOTE: harmonic mean of TPR and PPV

    :param true_pos: number of correctly identified positives
    :type true_pos: int
    :param false_pos: number of incorrectly identified positives
    :type false_pos: int
    :param false_neg: number of incorrectly identified negatives
    :type false_neg: int
    :return: f1
    """
    tpr = true_positive_rate(true_pos, false_neg)
    ppv = positive_predictive_value(true_pos, false_pos)
    return 2 * float(tpr * ppv)/(tpr + ppv)


def accuracy(true_pos, false_neg, true_neg, false_pos):
    """
    Calculates the accuracy.

    :param true_pos: number of correctly identified positives
    :type true_pos: int
    :param false_neg: number of incorrectly identified negatives
    :type false_neg: int
    :param true_neg: number of correctly identified negatives
    :type true_neg: int
    :param false_pos: number of incorrectly identified positives
    :type false_pos: int
    :return: acc
    """
    numer = true_pos + true_neg
    denom = true_pos + false_neg + true_neg + false_pos
    return float(numer)/denom


def matthews_correlation_coefficient(true_pos, false_neg, true_neg, false_pos):
    """
    Calculates the Matthews correlation coefficient.

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


def minmax(values):
    """
    Identifies the minimum and maximum values in a list of floats in a
    single pass (as opposed to calling min() and max() separately).

    :param values: a list of floats to find min() and max() values in
    :type values: list of float
    :return: minimum, maximum
    """
    minimum = maximum = None
    for value in values:
        if minimum is None or value < minimum:
            minimum = value
        if maximum is None or value > maximum:
            maximum = value
    return minimum, maximum


def transform(values, min_t=0, max_t=0):
    min_o, max_o = minmax(values)

    try:
        for i, value in enumerate(values):
            values[i] = (max_t - min_t) * float(value - min_o)/(max_o - min_o) + min_t
    except ZeroDivisionError:
        raise ValueError(f"cannot transform to range [{min_t},{max_t}] because "
                         f"all values are identical")

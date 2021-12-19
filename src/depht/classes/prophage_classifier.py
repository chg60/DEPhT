import pandas as pd
import plotly.express as px


blue, red = "#0008ff", "#ff0800"
colormap = {"prophage": blue, "bacteria": red}


class ProphageClassifier:
    """
    Classifier that uses a Naive Bayes-type strategy to construct
    probability distributions for an n-feature space
    """
    def __init__(self):
        """
        Constructor for instances of ProphageClassifier.
        """
        self.distributions_ = dict()

    @property
    def n_features(self):
        return len(self.distributions_)

    def fit(self, x, y, plot=False):
        """
        Fits a model that attempts to differentiate between

        :param x: training data to fit a model against
        :type x: pandas.DataFrame
        :param y: training data class labels
        :type y: list or numpy.ndarray or pd.DataFrame
        :param plot: show distribution of feature space?
        :type plot: bool
        """
        # Build a probability distribution for each feature
        prophage_df = x.iloc[[i for i, value in enumerate(y) if value == 1], :]
        bacteria_df = x.iloc[[i for i, value in enumerate(y) if value == 0], :]

        ratio = float(len(bacteria_df)) / len(prophage_df)

        figs = list()
        for i, feature in enumerate(x.columns):
            prophage_series = prophage_df.loc[:, feature]
            bacteria_series = bacteria_df.loc[:, feature]

            prophage_hist = Histogram(prophage_series)
            bacteria_hist = Histogram(bacteria_series, prophage_hist.bin_width)

            if plot:
                pdf = prophage_hist.as_dataframe()
                pdf["class"] = ["prophage"] * len(pdf)
                pdf["values"] = [float(x)/prophage_hist.n_samples * 100 for x
                                 in pdf["values"]]
                pdf["color"] = [blue] * len(pdf)

                bdf = bacteria_hist.as_dataframe()
                bdf["class"] = ["bacteria"] * len(bdf)
                bdf["values"] = [float(x)/bacteria_hist.n_samples * 100 for x
                                 in bdf["values"]]
                bdf["color"] = [red] * len(bdf)

                df = pd.concat([pdf, bdf])
                fig = px.bar(data_frame=df, x="bins", y="values",
                             color="class", barmode="overlay",
                             color_discrete_map=colormap,
                             template="simple_white")
                figs.append((feature, fig))

            dist = ProbabilityDistribution(prophage_hist, bacteria_hist,
                                           weights=[1, ratio])

            self.distributions_[feature] = dist

        return figs

    def predict_proba(self, x, feature_weights=None):
        """
        Uses the known probability distributions for input features
        that each input row belongs to a prophage.

        NOTE: sum(feature_weights) must equal 1.

        :param x: the data rows to make predictions for
        :type x: pandas.DataFrame
        :param feature_weights: how much weigh to place on each feature
        :type feature_weights: list of float
        :return: temp_proba
        """
        if len(x.columns) != self.n_features:
            raise ValueError(f"input dataframe must have {self.n_features} "
                             f"columns")

        if feature_weights and len(feature_weights) != self.n_features:
            raise ValueError(f"feature_weights should have {self.n_features} "
                             f"elements")

        if feature_weights and sum(feature_weights) != 1:
            raise ValueError("feature_weights must sum to 1")

        if not feature_weights:
            feature_weights = [1.0 / self.n_features] * self.n_features

        temp_proba = list()
        for i, row in x.iterrows():
            probabilities = list()
            for feature in row.keys().values:
                value = row[feature]
                proba = self.distributions_[feature].get_probability(value)
                probabilities.append(proba)
            row_proba = 0
            for weight, proba in zip(feature_weights, probabilities):
                row_proba += weight * proba
            temp_proba.append(row_proba)

        return temp_proba

    def predict(self, x, feature_weights=None, alpha=0.5):
        """
        Get binary class predictions by calling predict_proba and
        setting values to 0 if < alpha, else 1.

        :param x: the data rows to make predictions for
        :type x: pandas.DataFrame
        :param feature_weights: how much weigh to place on each feature
        :type feature_weights: list of float
        :param alpha: threshold for class 1
        :type alpha: float
        :return: temp_predict
        """
        proba = self.predict_proba(x, feature_weights)
        temp_predict = list()
        for p in proba:
            if p >= alpha:
                temp_predict.append(1)
            else:
                temp_predict.append(0)

        return temp_predict

    def __len__(self):
        return len(self.distributions_)


class Histogram:
    def __init__(self, data, bin_width=None):
        """
        Constructor for instances of the Histogram class.

        :param data: data to build the histogram from
        :type data: pandas.Series
        :param bin_width: the desired histogram bin width
        :type bin_width: int or float
        """
        self.hist = dict()
        self.n_samples = len(data)

        if bin_width:
            self.bin_width = bin_width
        else:
            data_range = data.max() - data.min()

            if 0 < data_range <= 1:         # For very small numbers
                self.bin_width = 0.01
            elif 0 < data_range <= 100:     # For smallish numbers
                self.bin_width = 1
            else:                           # Everything else
                self.bin_width = 10

        for value in data:
            index = round(value // self.bin_width * self.bin_width, 2)
            if index in self.hist:
                self.hist[index] += 1
            else:
                self.hist[index] = 1

    def as_dataframe(self):
        return pd.DataFrame({"bins": self.hist.keys(),
                             "values": self.hist.values()})


class ProbabilityDistribution:
    def __init__(self, prophage_hist, bacteria_hist, weights=None):
        """

        :param prophage_hist: prophage histogram
        :type prophage_hist: Histogram
        :param bacteria_hist: bacteria histogram
        :type bacteria_hist: Histogram
        :param weights: adjust weights to account for uneven sampling
        :type weights: list of int or list of float or None
        """
        if not prophage_hist.bin_width == bacteria_hist.bin_width:
            raise ValueError("input histograms must have the same bin_width")

        self.distribution = dict()

        self.minimum = None
        self.maximum = None

        self.bin_width = prophage_hist.bin_width

        if weights:
            prophage_weight, bacteria_weight = weights
        else:
            prophage_weight, bacteria_weight = 1, 1

        prophage_count = prophage_hist.n_samples * prophage_weight
        bacteria_count = bacteria_hist.n_samples * bacteria_weight

        temp_distribution = dict()

        for prophage_key, prophage_value in prophage_hist.hist.items():
            numer = float(prophage_value) / prophage_count * prophage_weight
            if prophage_key in bacteria_hist.hist:
                bacteria_value = bacteria_hist.hist[prophage_key]
                bacteria_proportion = float(bacteria_value) / bacteria_count
                bacteria_proportion *= bacteria_weight
            else:
                bacteria_proportion = 0.0
            denom = numer + bacteria_proportion
            prophage_probability = round(float(numer) / denom, 2)
            temp_distribution[prophage_key] = prophage_probability

        for bacteria_key, bacteria_value in bacteria_hist.hist.items():
            if bacteria_key not in temp_distribution:
                temp_distribution[bacteria_key] = 0.0

        sorted_keys = sorted(temp_distribution.keys())
        self.minimum, self.maximum = sorted_keys[0], sorted_keys[-1]
        for key in sorted_keys:
            self.distribution[key] = temp_distribution[key]

        self._smooth()

    def _smooth(self):
        expect_steps = (self.maximum - self.minimum + self.bin_width) // \
                       self.bin_width

        # Fill missing values with their left neighbor
        for i in range(round(expect_steps)):
            key = round(self.minimum + (i * self.bin_width), 2)
            if key not in self.distribution:
                prior_key = round(key - self.bin_width, 2)
                self.distribution[key] = self.distribution[prior_key]

    def get_probability(self, value):
        if value < self.minimum:
            probability = 1.0
        elif value > self.maximum:
            probability = 0.0
        else:
            key = round(value // self.bin_width * self.bin_width, 2)
            probability = self.distribution[key]

        return probability

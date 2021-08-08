import random


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

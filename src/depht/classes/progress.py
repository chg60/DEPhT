"""Class for easy visualization of progress level for multiprocessing
or multithreading."""


class ProgressBar:
    """Inline-updatable progress bar."""
    def __init__(self, current, end, width=50):
        self._percent = int(float(current) / end * 100)
        self._width = int(width)
        self._ratio = int(100 / width)
        self._multi = int(self._percent / self._ratio)
        self._pad = self._width - self._multi

    def show(self):
        """Print ProgressBar to the console."""
        if self._percent < 100:
            print("\r" + str(self), end="")
        else:
            print("\r" + str(self))

    def __str__(self):
        s = f"[{'#' * self._multi}{' ' * self._pad}] {self._percent}%"
        return s


def show_progress(current, end, width=50):
    """Create an instance of ProgressBar and print it in-line.

    :param current: current step (1 through n)
    :type current: int
    :param end: number of steps (n)
    :type end: int
    :param width: character width for the progressbar
    :type width: int
    :return: progressbar
    """
    ProgressBar(current, end, width).show()

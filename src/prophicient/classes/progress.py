class Progress:
    def __init__(self, current, end, width=50):
        self._percent = int(float(current) / end * 100)
        self._width = int(width)
        self._width_ratio = int(100 / width)
        self._multiplier = int(self._percent / self._width_ratio)
        self._padding = self._width - self._multiplier

    def __str__(self):
        s = f"[{'#' * self._multiplier}{' ' * self._padding}] {self._percent}%"
        return s

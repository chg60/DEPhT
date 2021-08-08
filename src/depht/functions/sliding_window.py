def leading_window(window_size, range_size):
    """
    Generator function that yields left/cursor/right index values for
    the prescribed `window_size` ahead of the cursor (inclusive) over
    the range `range_size`.

    :param window_size: size of window (including cursor)
    :type window_size: int
    :param range_size: size of the index range
    :type range_size: int
    """
    for cursor in range(range_size):
        left, right = cursor, cursor + window_size - 1
        yield left, cursor, right


def lagging_window(window_size, range_size):
    """
    Generator function that yields left/cursor/right index values for
    the prescribed `window_size` behind the cursor (inclusive) over
    the range `range_size`.

    :param window_size: size of window (including cursor)
    :type window_size: int
    :param range_size: size of the index range
    :type range_size: int
    """
    for cursor in range(range_size):
        left, right = cursor - window_size + 1, cursor
        yield left, cursor, right


def center_window(window_size, range_size):
    """
    Generator function that yields left/cursor/right index values for
    the prescribed `window_size` around the cursor (inclusive) over the
    range `range_size`.

    NOTE: for even values of `window_size`, windows will be asymmetric,
    with the right side of the window longer than the left side by 1
    index.

    :param window_size: size of window (including cursor)
    :type window_size: int
    :param range_size: size of the index range
    :type range_size: int
    """
    for cursor in range(range_size):
        half = (window_size - 1) // 2
        left, right = cursor - half, cursor + half
        if window_size % 2 == 0:
            right += 1
        yield left, cursor, right

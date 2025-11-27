from math import isclose
from typing import List
import numpy as np
from thermal_time.calcs import calc_thermal_time_range
from .units import PiecewiseFunction, TimeUnit


def PLF_value(points: List[List[float]], x: float) -> float:
    """Calculate the value of a piecewise linear function at a particular x value.

    Given an array of points (2xN) which describe a piecewise linear function,
    this function gets the y value for the given x value.  Values before/after
    the range of the function are set to the first/last y value respectively.
    Within the range of the function, zero-width segments are ignored so that
    discontinuous functions can be defined - be aware that the the value at the
    shared x value will come from the first of the two segments that share it.

    real, dimension(: , : ), intent( in ) : : points
    real, intent(in ) : : x
    real: : y
    """
    # TODO: sanity-check points: should be size 2 in first dimension, and
    #       x values should never decrease.

    n = len(points[1]) - 1
    if x < points[0][0]:
        y = points[1][0]
    elif x > points[0][n]:
        y = points[1][n]
    else:
        bx = points[0][0]
        by = points[1][0]

        for i in range(1, n + 1):
            ax = bx
            ay = by
            bx = points[0][i]
            by = points[1][i]

            # Skip zero-width pieces(this should be equivalent to an # equality check,
            # but checking floating point equality is evil # and the compiler warns about it)
            if isclose(abs(ax - bx), 0.0):
                continue

            if x <= bx:
                y = ay + (by - ay) * ((x - ax) / (bx - ax))
                break
    return y


def get_PLF_value(
    x_values: List[TimeUnit],
    y_values: List[float],
    x: float,
) -> float:
    """Calculate the value of a piecewise linear function at a particular x value.
    Uses numpy interp for simplicity.

    """

    return float(np.interp(x, x_values, y_values))


def offset(arr, zero, wrap) -> list:
    """offset an array to be relative to a particular *zero* value.  If *wrap*
    is supplied, new values less than 0 are offset by the *wrap* value.
    Note: formally known as reindex
    """
    out_a = [i - zero for i in arr]
    out_b = [(i + wrap) if i < 0.0 else i for i in out_a] if wrap else out_a
    return out_b


def generate_hours_data(day_count: int) -> List[int]:
    return np.array([[hr for hr in range(24)] for _ in range(day_count)]).reshape(day_count * 24)


def generate_days_data(day_count: int) -> List[int]:
    return np.array([[dd for _ in range(24)] for dd in range(day_count)]).reshape(day_count * 24)


def generate_example_td_data(day_count: int, T_b) -> List[float]:
    hr_data = generate_hours_data(day_count)
    demo_temp_data = [24 - abs(hr - 12) for hr in hr_data]
    tsc = demo_temp_data
    td_data = calc_thermal_time_range(tsc, t_b=T_b)
    return td_data


def get_day_from_td(target_td, dd_data, td_range, INVALID_VAL=999) -> int:
    """Helper to get the julian day from a thermal time value.

    Parameters
    ----------
    target_td : [type]
        The thermal time value at target day
    dd_data : [type]
        Array of dd data
    td_range : [type]
        Array of thermal time data
    INVALID_VAL : int, optional
        Value to return if day not found, by default 999

    Returns
    -------
    int
        the Day that matches thermal time value.

    """
    try:
        i, d = next((i, dd) for i, (td, dd) in enumerate(zip(td_range, dd_data)) if td >= target_td)
    except StopIteration:
        d = INVALID_VAL
        i = INVALID_VAL
    return i, d


def get_td_from_day(target_dd, dd_data, td_range, INVALID_VAL=999) -> int:
    """Helper to get the julian day from a thermal time value.

    Parameters
    ----------
    target_dd : [type]
        The day to get thermal time at
    dd_data : [type]
        Array of dd data
    td_range : [type]
        Array of thermal time data
    INVALID_VAL : int, optional
        Value to return if day not found, by default 999

    Returns
    -------
    int
        the Day that matches thermal time value.

    """
    try:
        i, d = next((i, td) for i, (td, dd) in enumerate(zip(td_range, dd_data)) if dd >= target_dd)
    except StopIteration:
        d = INVALID_VAL
        i = INVALID_VAL
    return i, d

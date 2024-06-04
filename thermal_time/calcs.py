"""Functions related to thermal time."""


import warnings
from math import ceil
from typing import List
import numpy as np


def calc_daily_thermal_time(
    effective_temperature: float,
    previous_thermal_time: float,
) -> float:
    return previous_thermal_time + effective_temperature


def calc_effective_thermal_time_range(
    Ts_C: List[float],
    dd_list: List[int] = None,
    t_b: float = 0,
    t_o: float = 20,
    t_m: float = 30,
    thermal_time_start: int = 0,
    thermal_time_offset: float = 0,
) -> List[float]:
    dd_count = ceil(len(Ts_C) / 24)
    padding = 24 - len(Ts_C) % 24
    Ts_C_padded = np.concatenate([Ts_C, np.zeros(padding)]) if padding < 24 else Ts_C
    t_at_zero_day = 0

    dd_list = dd_list if dd_list is not None else [
        j for j in range(dd_count * 24) for _ in range(24)]

    acc_temp_array = np.array(Ts_C_padded).reshape(dd_count, 24).sum(axis=1)
    effective_temp_array = np.array([calc_effective_temperature(
        t, t_b=t_b, t_o=t_o, t_m=t_m)for t in acc_temp_array])

    td_data = np.add.accumulate(effective_temp_array)
    td_data_flat = (td_data[:, np.newaxis] + np.zeros(24)).reshape(24 * dd_count)
    t_at_zero_day = next(t for t, dd in zip(td_data_flat, dd_list) if dd >= thermal_time_start)
    td_data_offset = td_data_flat - t_at_zero_day + thermal_time_offset
    return td_data_offset


def calc_thermal_time_range(
    Ts_C: List[float],
    dd_list: List[int] = None,
    t_b: float = 0,
    thermal_time_start: int = 0,
    thermal_time_offset: float = 0,
) -> List[float]:
    """Calculate the thermal temperature at each hour from start day for day count.

    Leaves are assumed to grow and senesce at a rate determined by the time integral
    of the mean daily temperature (thermal time) above a base temperature of 1 °C ((°C day) tb=1)
     (Porter 1984)

    We use a base temperature of 0.

    This is currently just the accumulated average daily temperature.

    Parameters
    ----------
    Ts_C : List[float]
        Ambient temperature row data where each row represents an hour [deg C]
    dd_list: List[int]
        list of days
    t_b : float, optional
        base temperature, by default 0
    thermal_time_start: int = 0
        day at which thermal time is 0
    thermal_time_offset: float = 0
        offset thermal time at thermal_time_start

    Returns
    -------
    Thermal_time: List[float]
        The thermal time as a list between start day to day count [deg days]

    """
    mean_temperature_accumulated = 0
    dd_count = ceil(len(Ts_C) / 24)
    mean_temp_array = np.zeros(dd_count * 24)
    t_at_zero_day = 0
    day_zero = dd_list is not None and dd_list[0] or 0

    dd_list = dd_list if dd_list is not None else [
        j for j in range(dd_count * 24) for _ in range(24)]

    for i, dd in enumerate(range(day_zero, day_zero + dd_count)):
        # Get the row id
        start_row = i * 24
        end_row = (i + 1) * 24
        temperature_array = Ts_C[start_row: end_row]
        mean_temp_day = sum(temperature_array) / max(1, len(temperature_array))
        # TODO: Should base temp be taken from each day
        # TODO: What happens if temperature is below base temp?
        mean_temp_above_base_temp = mean_temp_day - t_b
        mean_temperature_accumulated += max(0, mean_temp_above_base_temp)

        start_row = i * 24
        end_row = (i + 1) * 24
        mean_temp_array[start_row:end_row] = [mean_temperature_accumulated for _ in range(24)]
        if dd == thermal_time_start:
            t_at_zero_day = mean_temp_array[start_row]
    if t_at_zero_day != 0 or thermal_time_offset != 0:
        mean_temp_array = mean_temp_array - t_at_zero_day + thermal_time_offset

    # TODO: implement season vars
    # if (season_vars%accumulation_period_method == 'thermal time') then
    #     call calc_seasonal_config(mean_temps_array, start_day, species, gsto_config, season_vars)
    # end if
    return mean_temp_array[0:len(Ts_C)]


def calc_thermal_time(
    Ts_C: List[float],
    td_prev: float,
    t_b: float = 0,
) -> float:
    """Calculate the thermal temperature at each day.

    Leaves are assumed to grow and senesce at a rate determined by the time integral
    of the mean daily temperature (thermal time) above a base temperature of 1 °C ((°C day) tb=1)
     (Porter 1984)

    We use a base temperature of 0.

    This is just the accumulated average daily temperature.

    Parameters
    ----------
    Ts_C : List[float]
        Ambient temperature data over 24 hours [deg C]
    td_prev: float
        Previous day thermal time [deg C days]
    t_b : float, optional
        base temperature, by default 0 [deg C]

    Returns
    -------
    Thermal_time: float
        The thermal time [deg C days]
    """
    mean_temperature = np.average(Ts_C)
    mean_above_base = mean_temperature - t_b if mean_temperature > t_b else 0
    td = td_prev + mean_above_base
    return td


def get_thermal_time_at_day(dd: int, td: List[float], start_day: int = 0) -> float:
    """Get the td at a specific day.

    Parameters
    ----------
    dd: int
        Day to check thermal time [jd]
    td: list(int)
        Should be a list of td for each hour
    start_day : int, optional
        first julian day, by default 0
    Returns
    -------
    td_dd: float
        thermal time at dd
    """
    try:
        if dd >= start_day:
            return td[int((dd - start_day) * 24)]
        else:
            # raise ValueError(f" {dd} is before start day {start_day}")
            # We should extrapolate td backwards here
            return td[0]
    except IndexError:
        raise ValueError(
            f"""{dd} is before or after day range {start_day} -> {start_day + len(td)/24}.
            {int((dd - start_day) * 24)} not in td data.""")


def get_day_at_thermal_time(td: float, td_list: List[float], dd_list: List[int]) -> int:
    """Get the td at a specific day.

    Parameters
    ----------
    td: float
        Thermal time at which to find date
    td_list: list(int)
        Should be a list of td for each hour
    dd_list
        list of day ids

    Returns
    -------
    dd: int
        Julian day when thermal time is greater than td

    """
    if len(dd_list) != len(td_list):
        raise ValueError(
            f"td_list and dd_list must be the same length. Got {len(dd_list)} and {len(td_list)}")
    try:
        dd = next(dd for dd, tdd in zip(dd_list, td_list) if tdd > td)
        assert dd is not None
        return dd
    except AssertionError:
        raise ValueError(F"td: {td} is outside of data range")
    except StopIteration:
        raise ValueError(F"td: {td} is outside of data range")

    # except IndexError:
    #     raise ValueError(
    #         f"""{dd} is before or after day range {start_day} -> {start_day + len(td)/24}.
    #         {int((dd - start_day) * 24)} not in td data.""")


def calc_effective_temperature(
    t_acc: float,
    t_b: float,
    t_o: float,
    t_m: float,
) -> float:
    """Calculate the effective temperature using ambient temp and plant params.

    Parameters
    ----------
    t_acc : float
        Accumulated Ambient Temperature range for day
    t_b : float
        Base temperature [degrees]
    t_o : float
        Optimum temperature [degrees]
    t_m : float
        Maximum temperature [degrees]

    Returns
    -------
    t_eff: float
        Effective temperature [degrees]

    """
    t_mean = t_acc / 24
    t_eff = None

    if t_mean < t_b:
        t_eff = 0
    elif t_b <= t_mean <= t_o:
        t_eff = t_mean - t_b
    elif t_o < t_mean <= t_m:
        t_eff = (t_o - t_b) * (1 - ((t_mean - t_o) / (t_m - t_o)))
    elif t_mean >= t_m:
        t_eff = 0
    return t_eff


def get_td_dd(
    dd: int,
    td: float,
    # TODO: remove defaulte below
    start_day: int = None,
    start_day_td: float = None,
    # @Deprecated inputs
    sowing_day: int = None,
    season_Astart_temperature: float = None,
) -> float:
    """Get the difference between current td and td at start day.

    Parameters
    ----------
    dd : int
        Current Day (Julian)
    td : List[float]
        Day Thermal Time [deg day]
    start_day : int
        Day to find difference from
    start_day_td : float
        Thermal time at season start
    row_index: int
        Current Row index

    Returns
    -------
    td_dd: float
        Difference between current thermal time and td at season start
    """
    if sowing_day is not None:
        start_day = sowing_day
        start_day_td = season_Astart_temperature
        warnings.warn("Using sowing_day input is deprecated. Use start_day instead",
                      DeprecationWarning)
    if dd >= start_day:
        # < During the growing season we get the difference between the temperature ...
        td_dd = td - start_day_td
    else:
        td_dd = 0  # < before growing season td_dd is 0
    return td_dd

"""Model setup function.

Returns the arguments for the run model function.
"""
from math import ceil
from typing import Tuple, List
import datetime

from deprecated import deprecated
from pyDO3SE.External_State.External_State_Config import DayCalcMethods, InputMethod
from pyDO3SE.overrides import Main_Overrides
from pyDO3SE.util.error_handling import DayRangeError


def generate_dd(
    start_day: int,
    end_day: int,
    row_count: int,
    pad: bool,
) -> List[int]:
    assert start_day is not None, f"Must supply Location.start_day in config to calculate dd"
    row_count_a = row_count if end_day is None else (end_day - start_day) * 24
    row_count_b = row_count_a if not pad else ceil(row_count / 24) * 24
    return [start_day + dd for dd in range(ceil(row_count_b / 24)) for _ in range(24)][0:row_count_b]


def generate_hr(
    start_day: int,
    end_day: int,
    row_count: int,
    pad: bool,
) -> List[int]:
    assert start_day is not None, f"Must supply Location.start_day in config to calculate dd"
    row_count_a = row_count if end_day is None else (end_day - start_day) * 24
    row_count_b = row_count_a if not pad else ceil(row_count / 24) * 24
    return [hr for _ in range(ceil(row_count_b / 24)) for hr in range(24)][0:row_count_b]


def generate_row_index(
    start_day: int,
    end_day: int,
    row_count: int,
    pad: bool,
) -> List[int]:
    assert start_day is not None, f"Must supply Location.start_day in config to calculate dd"
    row_count_a = row_count if end_day is None else (end_day - start_day) * 24
    row_count_b = row_count_a if not pad else ceil(row_count / 24) * 24
    return [(24 * dd) + hr for dd in range(ceil(row_count_b / 24)) for hr in range(24)][0:row_count_b]


def generate_dd_and_hr(
    start_day: int,
    end_day: int,
    row_count: int,
    pad: bool = False,
):
    dd_out = generate_dd(start_day, end_day, row_count, pad)
    hr_out = generate_hr(start_day, end_day, row_count, pad)
    row_index_out = generate_row_index(start_day, end_day, row_count, pad)
    return dd_out, hr_out, row_index_out


def setup_hr(
    dds: List[int],
    hrs: List[int],
    row_indexes: List[int] = None,
    hr_offset: int = 0,
) -> Tuple[List[int], List[int], int, int]:
    """Setup hr data to manage timezone offsets and cropping.

    Parameters
    ----------
    dds : List[int]
        external state day data
    hrs : List[int]
        external state hr data
    row_indexes : List[int], optional
        external state row_index data, by default None
    hr_offset: float = 0
        Hour offset (for timezone correction)
    crop_to_day_start: bool = True
        Crop the data at the start to whole days only
    crop_to_day_end: bool = True
        Crop the data at the end to whole days only

    Returns
    -------
    Tuple[List[int], List[int], int, int]
        dd_out, hr_out, start_row, end_row

    """
    new_hrs = [(hr + hr_offset) % 24 for hr in hrs]
    day_offset = [1 if (hr + hr_offset) >= 24 else -1 if (hr + hr_offset) <
                  0 else 0 for hr in hrs]
    new_row_indexes = row_indexes if row_indexes is not None else list(
        range(len(new_hrs)))
    new_days = [dd + day_offset[i] for i, dd in enumerate(dds)]
    return new_days, new_hrs, new_row_indexes


def setup_dd(
    dds: List[int],
    hrs: List[int],
    row_indexes: List[int] = None,
    start_day: int = None,
    end_day: int = None,
    pad: bool = False,
) -> Tuple[List[int], List[int], int, int]:
    """Reconstruct dd data to manage wrapping around a year and cropping or padding.

    Parameters
    ----------
    dds : List[int]
        external state day data
    hrs : List[int]
        external state hr data
    row_indexes : List[int]
        external state row_index data
    start_day : _type_, optional
        start day index inclusive
    end_day : _type_, optional
        end day index inclusive
    pad: boolean, optional
        If True then will pad data to full day lengths

    Returns
    -------
    Tuple[List[int],List[int], int, int]
        dd_out, hr_out, start_row, end_row

    Raises
    ------
    DayRangeError
        _description_

    """
    dd_out = []
    assert dds is not None and dds[0] is not None, "Missing external state day data"
    d_prev = dds[0]
    wrap_d = 0
    start_row = None
    end_row = None
    zero_val = 1

    hr_start, hr_end = hrs[0], hrs[-1]

    if start_day is not None and start_day < dds[0]:
        raise DayRangeError(f"Start day: {start_day} is before start of data: {dds[0]}")

    if start_day is not None and start_day > dds[-1]:
        raise DayRangeError(f"Start day: {start_day} is after end of data: {dds[0]}")

    if end_day is not None and end_day < dds[0]:
        raise DayRangeError(f"End day: {end_day} is before start of data: {dds[0]}")

    if pad and (hr_start != 0 or hr_end != 23):
        # Repeats the first and last day to pad the data
        dd_padded = [dds[0] for _ in range(hr_start)] + dds + [dds[-1] for _ in range(23 - hr_end)]
        row_indexes_padded = list(range(-hr_start, len(hrs) + 24 - hr_end - 1)
                                  )
        hr_start = 0
        hr_end = 23
    else:
        dd_padded = dds
        row_indexes_padded = row_indexes or list(range(len(hrs)))

    dd_acc = []
    hr_out = []
    row_index_out = []
    curr_hour = hr_start
    for i, d in zip(row_indexes_padded, dd_padded):

        if d < d_prev:
            dd_acc.append(d_prev)
            zero_val = 1 - d
            wrap_d = sum(dd_acc) + zero_val
        new_dd = d + wrap_d

        if end_day is not None and new_dd > end_day:
            break
        if start_day is None or new_dd >= start_day:
            hr_out.append(curr_hour)
            dd_out.append(new_dd)
            row_index_out.append(i)
            if start_row is None:
                start_row = i
        end_row = i
        curr_hour = curr_hour + 1 if curr_hour < 23 else 0
        d_prev = d

    if start_day is not None and start_day > dd_out[-1]:
        raise DayRangeError(f"Start day: {start_day} is after end of data: {dds[-1]}")

    if end_day is not None and dd_out[-1] < end_day:
        raise DayRangeError(f"End day: {end_day} is after end of data: {dd_out[-1]}")

    return dd_out, hr_out, row_index_out, start_row, end_row


def convert_year_month_day_to_jd(
    year: List[int],
    month: List[int],
    dom: List[int],
    allow_wrap: bool = True,
    base_year: int = None,
) -> int:
    """Converts year month day of month input to julian day.

    Parameters
    ----------
    year : List[int]
        Year int
    month : List[int]
        month in 1-12
    dom : List[int]
        day of month 1-30
    allow_wrap : bool, optional
        If true allows dd to go back to 0 when it goes beyond end of year, by default False
    base_year : int, optional
        Used to determine if date is wrapped

    Returns
    -------
    int
        Current Julian day(dd)

    """
    d = datetime.datetime(year=year, month=month, day=dom).timetuple().tm_yday
    if not allow_wrap and base_year and year > base_year:
        d += sum([datetime.datetime(year=y, month=12,
                 day=31).timetuple().tm_yday for y in range(base_year, year)])
    return d


def get_day_from_input(
    years: List[int],
    months: List[int],
    doms: List[int],
    hrs: List[int],
    start_day: int = None,
    end_day: int = None,
    pad: bool = None,
) -> List[int]:
    dd = list(map(lambda args: convert_year_month_day_to_jd(
        *args, allow_wrap=False, base_year=years[0]), zip(years, months, doms)))

    year_start = years[0]
    hour_start = hrs[0]
    hour_end = hrs[-1]
    hour_end_pad = 24 - hour_end
    hrs_padded = list(range(0, hour_start)) + hrs + list(range(hour_end + 1, 24)) if pad else hrs
    dds_padded = [dd[0]for _ in range(0, hour_start)] \
        + dd \
        + [dd[-1] for _ in range(hour_end + 1, 24)] if pad else dd

    row_indexes = list(range(-hour_start, len(hrs) + hour_end_pad - 1)
                       ) if pad else list(range(len(hrs)))

    start_day_index = next((i for i, d in enumerate(dds_padded) if d >=
                           start_day)) if start_day is not None else 0
    end_day_index = next((i for i, d in sorted(enumerate(dds_padded), reverse=True)
                         if d <= end_day)) if end_day is not None else len(dds_padded) - 1

    dds_out = dds_padded[start_day_index: end_day_index + 1]
    hrs_out = hrs_padded[start_day_index: end_day_index + 1]
    row_indexes_out = row_indexes[start_day_index: end_day_index + 1]
    return dds_out, hrs_out, row_indexes_out


def get_day_crop_override(
    config_start_day: int,
    config_end_day: int,
    overrides: Main_Overrides = None,
) -> Tuple[int, int]:
    start_day_override = next(int(sd) for sd in [
        overrides and overrides.start_day,
        config_start_day,
        False,
    ] if sd is not None) or None

    end_day_override = next(int(sd) for sd in [
        overrides and overrides.end_day,
        config_end_day,
        False
    ] if sd is not None) or None

    if start_day_override is not None and end_day_override is not None and start_day_override > end_day_override:
        raise DayRangeError(f"Start day: {start_day_override} is after end day {end_day_override}")
    return start_day_override, end_day_override


def crop_data_to_full_days(
    hrs: List[int],
    row_indexes: List[int] = None,
    crop_to_day_start: bool = False,
    crop_to_day_end: bool = False,
):
    """Crop data to full days.

    We crop data by setting the row indexes that should be used.
    Later in the model we will crop the data to the row indexes.

    I.e. an output of [3,4,5] will later crop [11,22,33,44,55,66] to [33,44,55]

    """
    if not crop_to_day_start and not crop_to_day_end:
        return row_indexes
    new_hour_start_index = next(i for i, hr in enumerate(
        hrs) if hr == 0) if crop_to_day_start else 0
    new_hour_end_index = len(hrs) - next(i for i, hr in enumerate(reversed(hrs))
                                         if hr == 0) - 1 if crop_to_day_end else len(hrs)
    new_row_indexes = row_indexes[new_hour_start_index:new_hour_end_index] if row_indexes is not None else list(
        range(new_hour_start_index, new_hour_end_index))
    return new_row_indexes

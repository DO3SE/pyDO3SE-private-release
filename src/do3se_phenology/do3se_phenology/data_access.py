"""Some data prep tools."""

from typing import List
from pathlib import Path
import pandas as pd

from thermal_time.calcs import calc_thermal_time_range


def get_td_from_file(
    filepath: Path,
    t_b: float = 0,
    tsc_header: str = "Ts_C",
    strict: bool = True,
):
    """Load a file with temperature data and output thermal time."""
    df = pd.read_csv(filepath)
    if tsc_header not in df:
        if not strict:
            return None
        raise Exception(f"tsc_header: {tsc_header} not in \"{filepath}\" headers")
    td = calc_thermal_time_range(df[tsc_header].values, t_b)
    return td


def get_var_from_file(
    filepath: Path,
    var_header: str,
    strict: bool = True,
):
    """Load a var from file."""
    df = pd.read_csv(filepath)
    if var_header not in df:
        if not strict:
            return None
        raise Exception(f"{var_header} not in \"{filepath}\" headers")
    return df[var_header].values


def get_vars_from_file(
    filepath: Path,
    var_headers: List[str],
    strict: bool = True,
):
    """Load a var from file."""
    df = pd.read_csv(filepath)
    df.columns = [h.lower() for h in df.columns]
    var_headers = [h.lower() for h in var_headers]
    if strict:
        for var_header in var_headers:
            if var_header not in df:
                raise Exception(f"{var_header} not in \"{filepath}\" headers")

    return [df[var_header].values for var_header in var_headers if var_header in df]

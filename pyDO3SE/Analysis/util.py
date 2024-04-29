"""Analysis module utility functions."""

from typing import Any, List
from pyDO3SE.util.Objects import Field
import pandas as pd

MONTHS = [
    "jan",
    "feb",
    "mar",
    "apr",
    "may",
    "jun",
    "jul",
    "aug",
    "sep",
    "oct",
    "nov",
    "dec",
]


def output_log_to_field_data(
    output_logs: List[dict],
    field: Field,
) -> List[Any]:
    return [row.get(field.id, 'MISSING') for row in output_logs]


def day_of_year_to_month(day: int) -> int:
    # assume day start with 0
    # assume month start with 0
    # TODO: Handle leap year.
    day_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month = 0
    start = 0
    while month < len(day_in_month):
        end = start + day_in_month[month]
        if start <= day < end:
            return month
        start += day_in_month[month]
        if month == 11:
            month = 0
        else:
            month += 1

    return month


def output_log_to_data_with_mm(
    output_logs: List[dict],
    field: Field,
) -> List[Any]:
    # add mm to data
    return [
        {
            "hr": item["hr"],
            "mm": day_of_year_to_month(item["dd"]),
            field.id: item[field.id],
        }
        for item in output_logs
    ]

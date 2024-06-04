"""Analysis module utility functions."""

from typing import Any, List
from pyDO3SE.util.Objects import Field


def output_log_to_field_data(
    output_logs: List[dict],
    field: Field,
) -> List[Any]:
    return [row[field.id] for row in output_logs]

import pytest
from do3se_phenology.utils import wrap_day_of_year, generate_days_data


@pytest.mark.parametrize(
    "in_dd, expected_dd",
    [
        [-200 - 365,165],
        [-200,165],
        [-2,363],
        [-1,364],
        [0,0],
        [100,100],
        [364,364],
        [365,0],
        [365 + 1,1],
        [365 + 2,2],
        [365 + 365,0],
        [365 + (365*3),0],
        [365 + (365*3) + 1,1],
        [365 + (365*3) - 1,364],
    ]
)
def test_wrap_day_of_year(in_dd: int, expected_dd: int):
    assert wrap_day_of_year(in_dd) == expected_dd


@pytest.mark.parametrize(
    "day_count",
    [
        10,
        20,
        300,
        365,
        366,
        400
    ]
)
def test_generate_days_data(day_count: int):
    dd_data = generate_days_data(
        day_count=day_count,
    )
    assert len(dd_data) == day_count * 24
    assert dd_data[0] == 0
    assert dd_data[24] == 1
    assert dd_data[-1] == day_count - 1

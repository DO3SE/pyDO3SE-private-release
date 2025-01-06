import pytest
from pyDO3SE.Analysis.util import day_of_year_to_month


class TestDayOfYearToMonth:

    @pytest.mark.parametrize(["day","month"], [
        [0,0],
        [364,11],
        [365,0],
    ])
    def test_should_return_month(self, day, month):
        month_out = day_of_year_to_month(day)
        assert month == month_out

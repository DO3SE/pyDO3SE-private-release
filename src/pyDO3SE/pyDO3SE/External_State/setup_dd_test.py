from math import ceil
import pytest

# Internal Libraries
from pyDO3SE.util.error_handling import DayRangeError

# Libraries we are testing
from .setup_dd import (
    convert_year_month_day_to_jd,
    get_day_from_input,
    setup_dd,
    setup_hr,
)


class TestSetupDD:

    def test_some_days(self):
        day_count = 10
        hr_in = [hr for _ in range(day_count) for hr in range(24)]
        dd_in = [dd for dd in range(day_count) for _ in range(24)]
        row_indexes_in = list(range(len(hr_in)))
        start_day = 0
        end_day = 9
        dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
            dds=dd_in,
            hrs=hr_in,
            row_indexes=row_indexes_in,
            start_day=start_day,
            end_day=end_day,
        )
        assert len(dd_out) > 0
        assert len(dd_out) == end_row - start_row + 1
        assert dd_out[0] == start_day
        assert dd_out[-1] == end_day

    def test_with_cropped_hours(self):
        start_hour, end_hour = 3, 20
        day_count = 10
        hr_in = [hr for _ in range(day_count) for hr in range(24)]
        hr_in = hr_in[start_hour: len(hr_in) - (23 - end_hour)]
        dd_in = [dd for dd in range(day_count) for _ in range(24)]
        dd_in = dd_in[start_hour: len(dd_in) - (23 - end_hour)]
        row_indexes_in = list(range(len(hr_in)))
        start_day = dd_in[0]
        end_day = dd_in[-1]

        dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
            dds=dd_in,
            hrs=hr_in,
            row_indexes=row_indexes_in,
            start_day=start_day,
            end_day=end_day,
        )

        assert dd_out[0] == start_day

        assert len(dd_out) > 0
        assert len(dd_out) == end_row - start_row + 1
        print(len(dd_out))
        assert len(hr_in) % 24 == (24 - ((23 - end_hour) + start_hour)) % 24
        assert len(dd_out) % 24 == (24 - ((23 - end_hour) + start_hour)) % 24
        assert dd_out[-1] == end_day

    def test_with_padded_hours(self):
        start_hour, end_hour = 16, 0
        day_count = 10
        days_in_month = 30
        hr_in = [hr for _ in range(day_count) for hr in range(24)]
        hr_in = hr_in[start_hour: len(hr_in) - (23 - end_hour)]
        start_day = 1
        dd_in = [dd for dd in range(1, day_count + 1) for _ in range(24)]
        dd_in = dd_in[start_hour: len(dd_in) - (23 - end_hour)]
        dom = [start_day + ceil((i - hr) / 24) for i, hr in enumerate(hr_in)]
        dom = [1 + ((d - 1) % days_in_month) for d in dom]
        dd_in = dom
        row_indexes_in = list(range(len(hr_in)))

        dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
            dds=dd_in,
            hrs=hr_in,
            row_indexes=row_indexes_in,
            start_day=None,
            end_day=None,
            pad=True,
        )
        assert len(dd_out) == len(hr_out) == len(row_indexes_out)

        assert all([d is not None for d in dd_out])
        assert all([h is not None for h in hr_out])

        assert dd_out[0] == start_day
        assert dd_out[-1] == dd_in[-1]
        assert start_row == -start_hour
        assert end_row == len(dd_out) - start_hour - 1

        assert len(dd_out) > 0
        assert len(dd_out) == end_row - start_row + 1
        assert len(dd_out) % 24 == 0

    def test_some_days_with_start_offset(self):
        day_count = 10
        hr_in = [hr for _ in range(day_count) for hr in range(24)]
        dd_in = [dd for dd in range(day_count) for _ in range(24)]
        row_indexes_in = list(range(len(hr_in)))
        start_day = 2
        end_day = 9
        dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
            dds=dd_in,
            hrs=hr_in,
            row_indexes=row_indexes_in,
            start_day=start_day,
            end_day=end_day,
        )

        assert dd_out[0] == start_day
        assert start_row == 2 * 24
        assert end_row == 10 * 24 - 1
        assert row_indexes_out[0] == start_row

        assert len(dd_out) > 0
        assert len(dd_out) == end_row - start_row + 1
        assert dd_out[-1] == end_day

    def test_does_not_have_invalid_hours(self):
        hr_in = list(range(30))
        dd_in = [1 if i < 23 else 2 for i, _ in enumerate(hr_in)]
        dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
            dds=dd_in,
            hrs=hr_in,
        )
        assert hr_out[-1] < 23
        assert max(hr_out) == 23

    @pytest.mark.parametrize(['dd', 'last_val'], [
        [range(0, 2), 1],
        [range(1, 2), 1],
        [range(1, 365), 364],
        [range(2, 365), 364],
        [range(1, 366), 365],
        [range(1, 367), 366],
        [range(100, 367), 366],
        [range(100, 400), 399],
        [range(100, 800), 799],
        [range(0, 730), 729],
        [range(0, 733), 732],

        # index 0 to 364
        [list(range(0, 365)) + list(range(0, 3)), 367],
        [list(range(200, 365)) + list(range(0, 3)), 367],
        [list(range(364, 365)) + list(range(0, 1)), 365],
        [list(range(364, 365)) + list(range(0, 2)), 366],

        # index 1 to 365
        [list(range(1, 366)) + list(range(1, 3)), 367],
        [list(range(200, 366)) + list(range(1, 3)), 367],
        [list(range(365, 366)) + list(range(1, 3)), 367],
        [list(range(365, 366)) + list(range(1, 2)), 366],

        [list(range(1, 366)) + list(range(1, 366)) + list(range(1, 366)), 1095],
        [list(range(0, 365)) + list(range(0, 365)) + list(range(0, 365)), 1093],
    ])
    @pytest.mark.parametrize(['start_day', 'end_day'], [
        [None, None],
        [1, None],
        [2, None],
        [None, 1],
        [None, 10],
    ])
    @pytest.mark.parametrize(['start_hour', 'end_hour'], [
        [0, 23],
    ])
    def test_setup_dd(self, dd, last_val, start_day, end_day, start_hour, end_hour):
        hr_in = [hr for dd in dd for hr in range(24)]
        hr_in = hr_in[start_hour: len(hr_in) - (23 - end_hour)]
        dd_in = [dd for dd in dd for _ in range(24)]
        dd_in = dd_in[start_hour: len(dd_in) - (23 - end_hour)]

        row_indexes_in = list(range(len(hr_in)))

        if start_day is None and end_day is None:
            dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
                dd_in, hr_in, row_indexes_in)
            assert len(dd_out) == len(dd_in)
            assert dd_out[-1] == last_val
            assert start_row == 0
            assert end_row == len(dd_out) - 1
            assert all([b == a + 1 or b == a for i, (a, b)
                       in enumerate(zip(dd_out, dd_out[1:])) if (i + 1) % 24 == 0])
        elif start_day is not None and end_day is None:
            if dd_in[0] <= start_day <= dd_in[-1]:
                # Test that input is clipped by start day
                dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
                    dd_in, hr_in, row_indexes_in, start_day, end_day)
                assert len(dd_out) == len(dd_in) - (start_day - dd_in[0]) * 24
                assert start_row == (start_day - dd_in[0]) * 24
                assert end_row == len(dd_in) - 1

            elif start_day < dd_in[0]:
                with pytest.raises(DayRangeError):
                    dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
                        dd_in, hr_in, row_indexes_in, start_day, end_day)
            elif start_day > dd_in[-1]:
                with pytest.raises(DayRangeError):
                    dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
                        dd_in, hr_in, row_indexes_in, start_day, end_day)
            elif start_day == dd_in[0]:
                dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
                    dd_in, hr_in, row_indexes_in, start_day, end_day)
                assert len(dd_out) == len(dd_in)
                assert dd_out[-1] == last_val
                assert end_row == len(dd_out) - 1
            else:
                raise NotImplementedError(
                    f"This condition is not tested! start_day: {start_day} dd_in[0]: {dd_in[0]}, dd_in[-1]: {dd_in[-1]}")

        elif start_day is None and end_day is not None:
            if end_day < dd_in[0]:
                with pytest.raises(DayRangeError):
                    dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
                        dd_in, hr_in, row_indexes_in, start_day, end_day)
                return

            try:
                dd_out, hr_out, row_indexes_out, start_row, end_row = setup_dd(
                    dd_in, hr_in, row_indexes_in, start_day, end_day)
            except DayRangeError as e:
                # TODO: Assert should only be if end day outside data
                return
            if len(dd_out) > 0:
                assert dd_out[-1] == end_day
                assert start_row == 0
                assert end_row == len(dd_out) - 1

            assert len(dd_out) <= len(dd_in)

        else:
            raise NotImplementedError("This condition is not tested!")


class TestSetupHR:

    def test_should_offset_hours(self):
        day_count = 10
        hr_in = [hr for _ in range(day_count) for hr in range(24)]
        dd_in = [dd for dd in range(day_count) for _ in range(24)]
        row_indexes_in = list(range(len(hr_in)))
        offset = 3
        dd_out, hr_out, row_indexes_out = setup_hr(
            dds=dd_in,
            hrs=hr_in,
            row_indexes=row_indexes_in,
            hr_offset=offset,
        )
        assert len(dd_out) == len(hr_out) == len(row_indexes_out), "Lengths should be the same"
        assert len(hr_out) > 0
        assert hr_out[0] == hr_in[0] + offset
        assert all([0 <= hr < 24 for hr in hr_out])

    def test_should_offset_hours_with_negative_offset(self):
        day_count = 10
        hr_in = [hr for _ in range(day_count) for hr in range(24)]
        dd_in = [dd for dd in range(day_count) for _ in range(24)]
        row_indexes_in = list(range(len(hr_in)))
        offset = -3
        dd_out, hr_out, row_indexes_out = setup_hr(
            dds=dd_in,
            hrs=hr_in,
            row_indexes=row_indexes_in,
            hr_offset=offset,
        )
        assert len(dd_out) == len(hr_out) == len(row_indexes_out), "Lengths should be the same"
        assert len(hr_out) > 0
        assert hr_out[0] == (hr_in[0] + offset) % 24


class TestConvertYearMonthDayToJd:

    @pytest.mark.parametrize(['year', 'month', 'dom', 'expected'], [
        [2019, 1, 1, 1],
        [2019, 1, 2, 2],
        [2019, 1, 31, 31],
        [2019, 2, 1, 32],
        [2019, 3, 1, 60],
        [2019, 4, 1, 91],
        [2019, 5, 1, 121],
        [2019, 6, 1, 152],
        [2019, 7, 1, 182],
        [2019, 8, 1, 213],
        [2019, 9, 1, 244],
        [2019, 10, 1, 274],
        [2019, 11, 1, 305],
        [2019, 12, 1, 335],
        [2019, 12, 2, 336],
    ])
    def test_value(self, year, month, dom, expected):
        assert convert_year_month_day_to_jd(year, month, dom) == expected


class TestGetDayFromInput:

    def test_value(self):
        day_count = 100
        hr_start = 3
        days_in_month = 3  # Month is shorted to generate data only
        m_start = 12
        y_start = 2019
        dom_start = 2

        hours = [(hr + hr_start) % 24 for _ in range(day_count) for hr in range(24)]
        dom = [dom_start + ceil((i - hr) / 24) for i, hr in enumerate(hours)]
        months = [m_start + ceil(d / days_in_month) - 1 for d in dom]

        years = [y_start + ceil(m / 12) - 1 for m in months]
        months = [1 + ((m - 1) % 12) for m in months]
        dom = [1 + ((d - 1) % days_in_month) for d in dom]
        assert len(hours) == len(months) == len(dom)
        assert years[0] == y_start
        assert dom[0] == dom_start
        assert months[0] == m_start

        out_dd, out_hr, row_indexes = get_day_from_input(years, months, dom, hours)
        assert len(out_dd) == len(out_hr) == len(row_indexes) == len(hours)
        assert out_dd[0] == 336
        assert row_indexes[0] == 0
        assert row_indexes[-1] == len(row_indexes) - 1

    def test_value_pad(self):
        day_count = 100
        hr_start = 3
        days_in_month = 3  # Month is shorted to generate data only
        m_start = 12
        y_start = 2019
        dom_start = 2

        hours = [(hr + hr_start) % 24 for _ in range(day_count) for hr in range(24)]
        dom = [dom_start + ceil((i - hr) / 24) for i, hr in enumerate(hours)]
        months = [m_start + ceil(d / days_in_month) - 1 for d in dom]

        years = [y_start + ceil(m / 12) - 1 for m in months]
        months = [1 + ((m - 1) % 12) for m in months]
        dom = [1 + ((d - 1) % days_in_month) for d in dom]
        assert len(hours) == len(months) == len(dom)
        assert years[0] == y_start
        assert dom[0] == dom_start
        assert months[0] == m_start

        out_dd, out_hr, row_indexes = get_day_from_input(years, months, dom, hours, pad=True)
        assert len(out_dd) == len(out_hr) == len(row_indexes)
        assert len(out_dd) % 24 == 0
        assert row_indexes[0] == - hr_start
        assert out_hr[0] == 0
        assert out_hr[-1] == 23
        assert out_dd[0] == 336

    def test_value_crop(self):
        day_count = 300
        hr_start = 3
        days_in_month = 3  # Month is shorted to generate data only
        m_start = 12
        y_start = 2019
        dom_start = 2

        hours = [(hr + hr_start) % 24 for _ in range(day_count) for hr in range(24)]
        dom = [dom_start + ceil((i - hr) / 24) for i, hr in enumerate(hours)]
        months = [m_start + ceil(d / days_in_month) - 1 for d in dom]

        years = [y_start + ceil(m / 12) - 1 for m in months]
        months = [1 + ((m - 1) % 12) for m in months]
        dom = [1 + ((d - 1) % days_in_month) for d in dom]
        assert len(hours) == len(months) == len(dom)
        assert years[0] == y_start
        assert dom[0] == dom_start
        assert months[0] == m_start

        start_day_crop = 337
        end_day_crop = 339

        out_dd, out_hr, row_indexes = get_day_from_input(
            years, months, dom, hours, start_day=start_day_crop, end_day=end_day_crop, pad=True)
        assert len(out_dd) == len(out_hr) == len(row_indexes)
        assert len(out_dd) % 24 == 0
        assert out_hr[0] == 0
        assert out_hr[-1] == 23
        assert out_dd[0] >= start_day_crop
        assert out_dd[-1] <= end_day_crop

    def test_value_crop_overlap(self):
        day_count = 300
        hr_start = 3
        days_in_month = 3  # Month is shorted to generate data only
        m_start = 10
        y_start = 2019
        dom_start = 2

        hours = [(hr + hr_start) % 24 for _ in range(day_count) for hr in range(24)]
        dom = [dom_start + ceil((i - hr) / 24) for i, hr in enumerate(hours)]
        months = [m_start + ceil(d / days_in_month) - 1 for d in dom]

        years = [y_start + ceil(m / 12) - 1 for m in months]
        months = [1 + ((m - 1) % 12) for m in months]
        dom = [1 + ((d - 1) % days_in_month) for d in dom]
        assert len(hours) == len(months) == len(dom)
        assert years[0] == y_start
        assert dom[0] == dom_start
        assert months[0] == m_start

        start_day_crop = 20
        end_day_crop = 305

        out_dd, out_hr, row_indexes = get_day_from_input(
            years, months, dom, hours, start_day=start_day_crop, end_day=end_day_crop, pad=True)

        assert out_dd[0] < out_dd[-1]
        assert len(out_dd) == len(out_hr) == len(row_indexes)
        assert len(out_dd) % 24 == 0
        assert out_hr[0] == 0
        assert out_hr[-1] == 23
        assert out_dd[0] >= start_day_crop
        assert out_dd[-1] <= end_day_crop

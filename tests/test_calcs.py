"""Thermal time funciton tests."""

from math import floor
import pytest
from math import isclose
import numpy as np
from thermal_time.calcs import \
    calc_effective_temperature, \
    get_day_at_thermal_time, \
    get_td_dd, \
    calc_thermal_time_range, \
    get_thermal_time_at_day, \
    calc_thermal_time, \
    calc_effective_thermal_time_range


def test_calc_effective_temperature():
    """Test calc effective temperature func."""
    t_eff = calc_effective_temperature(sum([3.7925 for i in range(0, 24)]), 0, 20, 30)
    assert isclose(t_eff, 3.7925, abs_tol=1e-3)

    t_eff = calc_effective_temperature(sum([21.845 for i in range(0, 24)]), 0, 20, 30)
    assert isclose(t_eff, 16.31, abs_tol=1e-3)


def test_get_td_dd():
    td_dd = get_td_dd(dd=4, td=100, sowing_day=3,
                      season_Astart_temperature=30)
    assert isclose(td_dd, 70.0, abs_tol=1e-3)


def test_calc_thermal_time():
    TsC_data = np.full((8, 24), 3).reshape(8 * 24)
    td = calc_thermal_time(TsC_data[0:24], td_prev=10, t_b=0)
    assert isclose(td, 13)


def test_calc_thermal_time_random_vals():
    TsC_data = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
    td = calc_thermal_time(TsC_data[0:24], td_prev=10, t_b=0)
    assert isclose(td, 15.260, abs_tol=1e-3)


def test_calc_thermal_time_random_vals_base_1():
    TsC_data = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
    td = calc_thermal_time(TsC_data[0:24], td_prev=10, t_b=1)
    assert isclose(td, 14.260, abs_tol=1e-3)


class TestCalcThermalTimeRange:

    def test_calc_thermal_time_range(self):
        TsC_data = np.ones((8, 24)).reshape(8 * 24)
        dd_data = np.arange(floor(len(TsC_data) / 24))
        td = calc_thermal_time_range(
            TsC_data,
            dd_data,
        )
        first_hours = td.reshape((8, 24)).transpose()[0]
        print(first_hours)
        assert np.array_equal(first_hours, [0, 1, 2, 3, 4, 5, 6, 7])
        assert len(td) == len(TsC_data)


    def test_calc_thermal_time_range_offset(self):
        TsC_data = np.ones((8, 24)).reshape(8 * 24)
        dd_data = np.arange(floor(len(TsC_data) / 24))
        td = calc_thermal_time_range(
            TsC_data,
            dd_data,
            t_b=0,
            thermal_time_start=3,
            thermal_time_offset=0,
        )
        first_hours = td.reshape((8, 24)).transpose()[0]
        print(first_hours)
        assert np.array_equal(first_hours, [-3., -2., -1., 0., 1., 2., 3., 4.])
        assert len(td) == len(TsC_data)

    def test_calc_thermal_time_range_not_multiple_of_24(self):
        TsC_data = np.ones((8, 24)).reshape(8 * 24)
        dd_data = np.arange(floor(len(TsC_data) / 24))
        TsC_data = TsC_data[0:-2]
        dd_data = dd_data[0:-2]
        td = calc_thermal_time_range(
            TsC_data,
            dd_data,
        )
        assert len(td) == len(TsC_data)


def test_get_thermal_time_at_dd():
    TsC_data = np.array([i for i in range(10 * 24)])
    thermal_time_at_dd = get_thermal_time_at_day(3, TsC_data, start_day=0)
    assert thermal_time_at_dd == 3 * 24


def test_get_thermal_time_at_dd_offset_day():
    TsC_data = np.array([i for i in range(10 * 24)])
    start_day = 2
    dd = 3
    thermal_time_at_dd = get_thermal_time_at_day(dd, TsC_data, start_day=start_day)
    assert thermal_time_at_dd == (dd - start_day) * 24


def test_get_thermal_time_at_dd_out_of_range():
    TsC_data = np.array([i for i in range(10 * 24)])
    start_day = 2
    dd = 999
    with pytest.raises(ValueError) as e:
        get_thermal_time_at_day(dd, TsC_data, start_day=start_day)


def test_get_dd_at_thermal_time():
    TsC_data = np.array([i for i in range(10 * 24)])
    dd_data = np.arange(0, len(TsC_data))
    dd = get_day_at_thermal_time(TsC_data[120], TsC_data, dd_data)
    assert dd == 121


def test_get_dd_at_thermal_time_out_of_bounds():
    TsC_data = np.array([i for i in range(10 * 24)])
    dd_data = np.arange(0, len(TsC_data))
    with pytest.raises(ValueError) as e:
        get_day_at_thermal_time(99999, TsC_data, dd_data)


def test_get_dd_at_thermal_time_incorrect_input():
    TsC_data = np.arange(10)
    dd_data = np.arange(20)
    with pytest.raises(ValueError) as e:
        get_day_at_thermal_time(1, TsC_data, dd_data)


class TestCalcEffectiveThermalTimeRange:

    def test_correct_output(self):
        DAY_COUNT = 10
        TsC_data = ((np.sin(np.arange(0,3.14 * DAY_COUNT,3.14/24)) + 0.8) * 15)[0:DAY_COUNT*24]
        dd_data = np.array([dd for dd in range(DAY_COUNT) for _ in range(24) ])
        t_b=0
        t_o=20
        t_m=30
        thermal_time_start=0
        thermal_time_offset=0
        out = calc_effective_thermal_time_range(
            TsC_data,
            dd_data,
            t_b,
            t_o,
            t_m,
            thermal_time_start,
            thermal_time_offset,
        )
        assert out.shape == (DAY_COUNT * 24,)
        assert len(out) == DAY_COUNT * 24
        assert sum(out) > 1000
        assert isclose(sum(out),  9603.8735, abs_tol=1e-3)

    def test_correct_output_diff_shape(self):
        DAY_COUNT = 10
        TsC_data = ((np.sin(np.arange(0,3.14 * DAY_COUNT,3.14/24)) + 0.8) * 15)[0:DAY_COUNT*24]
        dd_data = np.array([dd for dd in range(DAY_COUNT) for _ in range(24) ])
        t_b=0
        t_o=20
        t_m=30
        thermal_time_start=0
        thermal_time_offset=0
        out = calc_effective_thermal_time_range(
            TsC_data[0:-2],
            dd_data[0:-2],
            t_b,
            t_o,
            t_m,
            thermal_time_start,
            thermal_time_offset,
        )
        assert out.shape == (DAY_COUNT * 24,)
        assert len(out) == DAY_COUNT * 24
        assert sum(out) > 1000
        assert isclose(sum(out),  9586.1776, abs_tol=1e-3)
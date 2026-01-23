"""Test the phenology functions."""

import pytest
import numpy as np
from math import isclose
from thermal_time.calcs import calc_thermal_time_range

from do3se_phenology.config import PhenologyKeyLengths
from do3se_phenology.state import LeafPhenologyStage, PhenologyStage
from do3se_phenology.phyllochron_dvi import (
    calc_dvi,
    calc_dvi_tt_PLF,
    calc_emergence_date_from_dvi,
    calc_life_stages_from_phyllochron,
    calc_phyllochron_from_dvi,
    calc_rpe,
    calc_t_l_from_phyllochron,
    calculate_t_l_from_leaf_f_phen,
    calculate_t_l_from_leaf_f_phen_b,
    estimate_t_lse_from_t_l,
    phyllochron_to_t_lem,
    phyllochron_to_t_lma,
    calc_t_l,
    calc_phyllochron,
    calc_emergence_rate,
    calc_emerged_leaf_count,
    calc_td_dd_per_leaf_pop,
    get_growing_populations,
)

PHYLLOCHRON = 82.1945
DL = 0.0679335


def test_phyllochron_to_t_lem():
    t_lem = phyllochron_to_t_lem(PHYLLOCHRON)
    assert isclose(t_lem, 147.95, abs_tol=1e-3)


def test_phyllochron_to_t_lma():
    t_lma = phyllochron_to_t_lma(PHYLLOCHRON)
    assert isclose(t_lma, 287.68, abs_tol=1e-3)


def test_calc_t_l_from_phyllochron():
    t_l = calc_t_l_from_phyllochron(PHYLLOCHRON)
    assert isclose(t_l, 435.63085, abs_tol=1e-3)
    assert calc_t_l_from_phyllochron(None) is None


def test_calc_life_stages_from_phyllochron():
    life_stages = calc_life_stages_from_phyllochron(PHYLLOCHRON, t_lse_constant=0.3)
    assert isclose(life_stages.t_l, 435.63085)
    # effectively ignored before phyllochron is set
    assert calc_life_stages_from_phyllochron(None, 0.3).t_l > 99998


def test_calc_t_l():
    t_l = calc_t_l(t_lem=147.95, t_lma=287.68)
    assert isclose(t_l, 435.63, abs_tol=1e-3)


def test_calc_phyllochron():
    phyllochron = calc_phyllochron(dl=DL)
    assert isclose(phyllochron, PHYLLOCHRON, abs_tol=1e-3)


def test_calc_phyllochron_from_dvi():
    phyllochron = calc_phyllochron_from_dvi(
        prev_dvi=-0.0590286,
        dvi=0.2030714,
        prev_pr=13.5754,
        pr=13.643524,
        prev_phyllochron=None,
    )
    assert isclose(phyllochron, 82.16100, abs_tol=1e-5)


def test_calc_phyllochron_from_dvi_b():
    t_emerg, d_emerg = calc_emergence_date_from_dvi(
        prev_dvi=-0.0590286,
        dvi=0.2030714,
        prev_t_emerg=99.9,
        prev_d_emerg=1,
        td=1000,
        dd=100,
    )
    assert d_emerg == 1
    assert isclose(t_emerg, 99.9, abs_tol=1e-5)

    t_emerg, d_emerg = calc_emergence_date_from_dvi(
        prev_dvi=-0.0590286,
        dvi=0.2030714,
        prev_t_emerg=None,
        prev_d_emerg=None,
        td=1000,
        dd=100,
    )
    assert d_emerg == 100
    assert isclose(t_emerg, 1000, abs_tol=1e-5)

    t_emerg, d_emerg = calc_emergence_date_from_dvi(
        prev_dvi=-1,
        dvi=-0.99,
        prev_t_emerg=None,
        prev_d_emerg=None,
        td=1000,
        dd=100,
    )
    assert d_emerg == None
    assert t_emerg == None


def test_calc_rpe():
    rpe = calc_rpe(24, 24, 0)
    assert isclose(rpe, 1, abs_tol=1e-3)


def test_calc_dvi():
    dvi = calc_dvi(
        prev_dvi=0.2210259,
        t_eff=3.281,
        tt_emr=35,
        tt_veg=1000,
        tt_rep=666.66667,
        rpe=1.0,
        dd=1,
        sowing_day=50,
    )

    assert isclose(dvi, -1, abs_tol=1e-5)

    dvi = calc_dvi(
        prev_dvi=0.2210259,
        t_eff=3.281,
        tt_emr=35,
        tt_veg=1000,
        tt_rep=666.66667,
        rpe=1.0,
        dd=100,
        sowing_day=50,
    )

    assert isclose(dvi, 0.2243069, abs_tol=1e-5)

    dvi = calc_dvi(
        prev_dvi=1.1234024,
        t_eff=14.63,
        tt_emr=35,
        tt_veg=1000,
        tt_rep=666.66667,
        rpe=1.0,
        dd=100,
        sowing_day=50,
    )

    assert isclose(dvi, 1.1453474, abs_tol=1e-5)


def test_calculate_t_l_from_leaf_f_phen():
    dd_start = 3
    dd_emerge = 2
    dd_mature = 3
    dd_sen = 2
    dd_end = 2
    dd_total = dd_start + dd_emerge + dd_mature + dd_sen + dd_end
    Ts_C_data = np.ones((dd_total, 24)).reshape(dd_total * 24)
    td_data = calc_thermal_time_range(Ts_C_data)
    dd_data = np.array([[d for i in range(24)] for d in range(dd_total)]).reshape(dd_total * 24)

    leaf_f_phen_data = np.concatenate((np.zeros(dd_start * 24), np.arange(0, 1, 1 / (24 * dd_emerge)), np.ones(
        24 * dd_mature), np.arange(1, 0, -1 / (24 * dd_sen)), np.zeros(24 * dd_end)))
    out = calculate_t_l_from_leaf_f_phen(
        td_data, leaf_f_phen_data, dd_data)

    assert isclose(out.t_lse, 2.31, abs_tol=1e-3)
    assert isclose(out.t_l, 15.217, abs_tol=1e-3)

    t_lse = out.t_lse
    t_l = out.t_l
    t_lma = out.t_lma  # t_lse / 0.33
    t_lem = out.t_lem  # 1.8 * (t_lma / 3.5)
    t_lep = t_lma - t_lse
    t_l_from_t_lse = t_lem + t_lma
    assert isclose(t_l_from_t_lse, t_l, abs_tol=1e-1)
    assert isclose(t_lse + t_lep, t_lma, abs_tol=1e-3)
    assert isclose(t_lse / 0.33, t_lma, abs_tol=1e-3)


def test_calculate_t_l_from_leaf_f_phen_b():
    dd_start = 3
    dd_emerge = 2
    dd_mature = 3
    dd_sen = 2
    dd_end = 2
    dd_total = dd_start + dd_emerge + dd_mature + dd_sen + dd_end
    Ts_C_data = np.ones((dd_total, 24)).reshape(dd_total * 24)
    td_data = calc_thermal_time_range(Ts_C_data)
    dd_data = np.array([[d for i in range(24)] for d in range(dd_total)]).reshape(dd_total * 24)

    leaf_f_phen_data = np.concatenate((np.zeros(dd_start * 24), np.arange(0, 1, 1 / (24 * dd_emerge)), np.ones(
        24 * dd_mature), np.arange(1, 0, -1 / (24 * dd_sen)), np.zeros(24 * dd_end)))
    out = calculate_t_l_from_leaf_f_phen_b(
        Ts_C_data, td_data, leaf_f_phen_data, dd_data, tt_emr=35, t_b=0, t_o=20, t_m=30)

    assert isclose(out.t_lse, 2, abs_tol=1e-3)
    assert isclose(out.t_l, 7.0, abs_tol=1e-3)
    assert out.sowing_day == 2

    t_lse = out.t_lse
    t_l = out.t_l

    t_lma = t_lse / 0.33
    t_lem = 1.8 * (t_lma / 3.5)

    t_lep = t_lma - t_lse

    # This method does not meet the below relationships
    # assert isclose(t_l, t_lem + t_lma, abs_tol=1e-3)
    # assert isclose(t_lse + t_lep, t_lma, abs_tol=1e-3)
    # assert isclose(t_lse / 0.33, t_lma, abs_tol=1e-3)


def test_calculate_t_l_from_leaf_f_phen_short():
    td_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    dd_data = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3]
    leaf_f_phen_data = [0, 1, 1, 0.8, 0.7, 0.3, 0.1, 0, 0, 0]
    out = calculate_t_l_from_leaf_f_phen(td_data, leaf_f_phen_data, dd_data)
    assert isclose(out.t_lse, 2.31, abs_tol=1e-3)
    assert isclose(out.t_l, 15.217, abs_tol=1e-3)


def test_estimate_t_lse_from_t_l():
    t_l = 1000
    t_lse = estimate_t_lse_from_t_l(t_l)
    assert isclose(t_lse, 219.373, abs_tol=1e-3)
    t_lma = t_lse / 0.33
    t_lem = 1.8 * (t_lma / 3.5)
    t_lep = t_lma - t_lse
    assert isclose(t_lem + t_lma, t_l, abs_tol=1e2)
    assert isclose(t_lse + t_lep, t_lma, abs_tol=1e-3)
    assert isclose(t_lse / 0.33, t_lma, abs_tol=1e-3)


@pytest.mark.parametrize(['nP', 't_flag_emerge', 'expected_emergence_rate'], [
    (1, 10, 0.00),
    (2, 10, 0.1),
    (3, 10, 0.2),
])
def test_calc_emergence_rate(nP, t_flag_emerge, expected_emergence_rate):
    out = calc_emergence_rate(nP, t_flag_emerge)
    assert isclose(out, expected_emergence_rate, abs_tol=1e-3)


def test_calc_emergence_rate_np_1():
    # Emergence rate should be 0 when nP is 1 so we only use flag leaf.
    out = calc_emergence_rate(1, 10)
    assert out == 0


class TestCalcEmergedLeafCount:

    @pytest.mark.parametrize(['td_dd', 'emergence_rate', 'expected_count'], [
        (0, 1, 0, ),  # at td_dd = 0
        (0, 0, 0, ),  # at td_dd = 0 and e_rate = 0
        (10, 0.1, 1, ),  # at td_dd = 10 and e_rate = 0.1
        (10, 0, 0, ),  # at td_dd = 10 and e_rate = 0
    ])
    def test_should_return_emerged_leaves_according_to_rate(self, td_dd, emergence_rate, expected_count):
        nP = 999
        t_emerge_flag = 999
        out = calc_emerged_leaf_count(nP, td_dd, emergence_rate, t_emerge_flag)
        assert out == expected_count

    @pytest.mark.parametrize('nP', [2, 4, 8])
    def test_should_limit_emerged_leaves_to_max_np(self, nP):
        td_dd = 999
        emergence_rate = 1
        t_emerge_flag = 999
        out = calc_emerged_leaf_count(nP, td_dd, emergence_rate, t_emerge_flag)
        assert out == nP

    @pytest.mark.parametrize('td_dd', [0, 10, 100])
    def test_should_use_flag_leaf_value_when_np_is_1(self, td_dd):
        nP = 1
        emergence_rate = 0
        t_emerge_flag = 9
        out = calc_emerged_leaf_count(nP, td_dd, emergence_rate, t_emerge_flag)
        if td_dd > t_emerge_flag:
            assert out == 1
        else:
            assert out == 0


class TestCalcTdddPerLeafPop:

    @pytest.mark.parametrize(['nP', 'emerged_leaf_count', 'td', 'td_prev', 'td_dd_prev', 'expected_td_dd'], [
        (1, 1, 100, 0, [10], [11]),
        (3, 2, 1, 0, np.zeros(3), [1, 1, 0]),
        (3, 2, 1, 0, np.zeros(5), [1, 1, 0]),
        (3, 0, 1, 0, np.zeros(5), [0, 0, 0]),
    ])
    def test_calc_td_dd_per_leaf_pop(self, nP, emerged_leaf_count, td, td_prev, td_dd_prev, expected_td_dd):
        t_emerg_to_flag_emerg = td - 11 if nP == 1 else None
        out = calc_td_dd_per_leaf_pop(nP, emerged_leaf_count, td, td_prev, td_dd_prev, t_emerg_to_flag_emerg)
        assert out == expected_td_dd

    def test_calc_td_dd_per_leaf_pop_flag_only(self):
        nP = 1
        emerged_leaf_count = 1
        td = 101
        td_prev = 0
        td_dd_prev = [0]
        t_emerg_to_flag_emerg = 100
        out = calc_td_dd_per_leaf_pop(nP, emerged_leaf_count, td, td_prev, td_dd_prev, t_emerg_to_flag_emerg)
        assert out == [1]


class TestGetGrowingPopulations:

    def test_get_growing_populations(self):
        td_dd_emerg = [0, 5, 20]
        leaf_population_t_lems = [10, 10]
        flag_leaf_t_lem = 10
        out = get_growing_populations(td_dd_emerg, leaf_population_t_lems, flag_leaf_t_lem)
        assert out == [False, True, False]


class TestCalcDviTTPLF:

    def test_can_get_dvi_from_thermal_time(self):
        dvi_interval=[
            (0, -1.000000001),
            (100, 0.0),
            (800.0, 1.0),
            (2000.0, 2.0),
        ]

        td = 0
        dvi = calc_dvi_tt_PLF(td, dvi_interval)
        assert isclose(dvi, -1.0, abs_tol=1e-3)

        td = 10
        dvi = calc_dvi_tt_PLF(td, dvi_interval)
        assert isclose(dvi, -0.90, abs_tol=1e-3)

        td = 900
        dvi = calc_dvi_tt_PLF(td, dvi_interval)
        assert isclose(dvi, 1.083, abs_tol=1e-3)

        td = 2000
        dvi = calc_dvi_tt_PLF(td, dvi_interval)
        assert isclose(dvi, 2.0, abs_tol=1e-3)


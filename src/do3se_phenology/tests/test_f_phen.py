from functools import partial
from itertools import accumulate
import pytest
from math import isclose

from do3se_phenology.f_phen import *


def test_f_phen_simple_PLF():
    """Tests the f_phen_simple_PLF function returns the correct value."""
    out = f_phen_simple_PLF(
        f_phen_1=30,
        f_phen_4=80,
        f_phen_a=20.20,
        f_phen_c=60.6,
        f_phen_e=30.4,
        SGS=1,
        EGS=200,
        dd=80,
    )
    assert isclose(out, 60.6, abs_tol=1e-3)


def test_tt_f_phen_simple_PLF():
    """Tests the tt_f_phen_simple_PLF function returns the correct value."""
    temperature_data = [2 for _ in range(1253)]
    accumulated_temperatures = list(accumulate(temperature_data))
    _get_f_phen = partial(tt_f_phen_simple_PLF,
                          t_f_phen_a=70,
                          t_f_phen_b=360,
                          t_f_phen_c=1145,
                          t_f_phen_d=1845,
                          f_phen_min=0.2,
                          td_at_sgs=100,
                          )
    f_phen = list(map(_get_f_phen, accumulated_temperatures))

    assert isclose(f_phen[0], 0.2, abs_tol=1e-3)
    assert isclose(f_phen[10], 0.2, abs_tol=1e-3)
    assert isclose(f_phen[100], 0.2883, abs_tol=1e-3)
    assert isclose(f_phen[500], 1, abs_tol=1e-3)
    assert isclose(f_phen[1000], 0, abs_tol=1e-3)
    assert isclose(f_phen[-1], 0, abs_tol=1e-3)


def test_f_phen_complex_PLF():
    """Tests the f_phen_complex_PLF function returns the correct value.

    The values to check have been matched against the DO3SE UI model
    """
    def calc_fphen_demo_p(dd):
        return f_phen_complex_PLF(
            f_phen_1=0,
            f_phen_2=1,
            f_phen_3=1,
            f_phen_4=45,
            f_phen_limA=0,
            f_phen_limB=0,
            f_phen_a=0.1,
            f_phen_b=1.0,
            f_phen_c=1.0,
            f_phen_d=1.0,
            f_phen_e=0.1,
            SGS=118,
            EGS=210,
            dd=dd,
        )
    assert isclose(calc_fphen_demo_p(dd=117), 0.0)  # ddadj = 364
    assert isclose(calc_fphen_demo_p(dd=118), 1.0)  # ddadj = 0 (365)
    assert isclose(calc_fphen_demo_p(dd=165), 1.0)
    assert isclose(calc_fphen_demo_p(dd=166), 0.98)
    assert isclose(calc_fphen_demo_p(dd=210), 0.1)
    assert isclose(calc_fphen_demo_p(dd=211), 0.0)


def test_leaf_f_phen_PLF():
    """Tests the leaf_f_phen_PLF function returns the correct value."""
    out = leaf_f_phen_PLF(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=12.0,
        leaf_f_phen_b=30.0,
        leaf_f_phen_c=8.0,
        Astart=1,
        Aend=100,
        dd=74,
    )
    assert isclose(out, 27.066, abs_tol=1e-3)


def test_leaf_f_phen_PLF_Astart():
    """Tests the leaf_f_phen_PLF function returns the correct value.

    Make sure that the fphen starts to increase at dd==Astart.
    """
    out = leaf_f_phen_PLF(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=152,
    )
    assert isclose(out, 0.0, abs_tol=1e-3)
    out = leaf_f_phen_PLF(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=153,
    )
    assert isclose(out, 0.8, abs_tol=1e-3)


def test_leaf_f_phen_PLF_Aend():
    """Tests the leaf_f_phen_PLF function returns the correct value.

    Make sure that the fphen hits 0 at dd==Aend + 1.
    """
    out = leaf_f_phen_PLF(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=208,
    )
    assert isclose(out, 0.2, abs_tol=1e-3)

    out = leaf_f_phen_PLF(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=209,
    )
    assert isclose(out, 0.0, abs_tol=1e-3)


def test_tt_leaf_f_phen_PLF():
    """Tests the tt_leaf_f_phen_PLF function returns the correct value."""
    temperature_data = [2 for _ in range(1253)]
    accumulated_temperatures = list(accumulate(temperature_data))
    _get_f_phen = partial(tt_leaf_f_phen_PLF,
                          t_leaf_f_phen_a=0.3,
                          t_leaf_f_phen_b=0.7,
                          t_leaf_f_phen_e=200,
                          t_leaf_f_phen_g=100,
                          t_leaf_f_phen_h=525,
                          t_leaf_f_phen_i=700,
                          t_astart=950,
                          td_at_sgs=100,
                          )
    leaf_f_phen = list(map(_get_f_phen, accumulated_temperatures))

    assert isclose(leaf_f_phen[0], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[10], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[600], 1, abs_tol=1e-3)
    assert isclose(leaf_f_phen[800], 0.8221, abs_tol=1e-3)
    assert isclose(leaf_f_phen[900], 0.592, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1000], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[-1], 0, abs_tol=1e-3)


def test_tt_leaf_f_phen_PLF_range():
    """Tests the tt_leaf_f_phen_PLF function returns the correct value."""
    temperature_data = [2 for _ in range(1253)]
    accumulated_temperatures = list(accumulate(temperature_data))

    leaf_f_phen = tt_leaf_f_phen_PLF_range(
        td_list=accumulated_temperatures,
        t_leaf_f_phen_a=0.3,
        t_leaf_f_phen_b=0.7,
        t_leaf_f_phen_e=200,
        t_leaf_f_phen_g=100,
        t_leaf_f_phen_h=525,
        t_leaf_f_phen_i=700,
        t_astart=950,
        td_at_sgs=100,
    )

    assert isclose(leaf_f_phen[0], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[10], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[600], 1, abs_tol=1e-3)
    assert isclose(leaf_f_phen[800], 0.8221, abs_tol=1e-3)
    assert isclose(leaf_f_phen[900], 0.592, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1000], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[-1], 0, abs_tol=1e-3)


def test_calc_leaf_f_phen_effect_on_V_cmax_25():
    out = calc_leaf_f_phen_effect_on_V_cmax_25(
        V_cmax_25_in=180,
        J_max_25_in=400,
        leaf_f_phen=1.23,
    )
    assert isclose(out.V_cmax_25, 221.4, abs_tol=1e-3)
    assert isclose(out.J_max_25, 492.0, abs_tol=1e-3)

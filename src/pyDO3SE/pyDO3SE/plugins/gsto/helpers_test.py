"""Gsto related helper functions.

Note: Values are random unless specified
"""

from math import isclose
from .helpers import (
    calc_gsto_leaf,
    calc_gsto_mean,
    temp_dep,
    temp_dep_inhibit,
)


def test_calc_gsto_leaf():
    """Tests the calc_gsto_leaf function returns the correct value."""
    out = calc_gsto_leaf(
        gp_gmax=99.9,
        gp_leaf_f_phen=0.3,
        gp_f_O3=0.3,
        gp_leaf_f_light=0.3,
        gp_fmin=0.3,
        gp_f_temp=0.3,
        gp_f_VPD=0.3,
        gp_f_SW=0.3,
    )
    assert isclose(out, 2.6973, abs_tol=1e-3)


def test_calc_gsto_mean():
    """Tests the calc_gsto_mean function returns the correct value."""
    out = calc_gsto_mean(
        gp_gmax=99.9,
        gp_gmorph=0.3,
        gp_f_phen=0.3,
        gp_f_light=0.3,
        gp_fmin=0.3,
        gp_f_temp=0.3,
        gp_f_VPD=0.3,
        gp_f_SW=0.3,
    )
    assert isclose(out, 0.80919, abs_tol=1e-3)


def test_temp_dep():
    """Tests the temp_dep function returns the correct value."""
    out = temp_dep(
        P_ref=99.9,
        T_ref=99.9,
        H_a=99.9,
        T=99.9,
        R=99.9,
    )
    assert isclose(out, 99.9, abs_tol=1e-3)


def test_temp_dep_inhibit():
    """Tests the temp_dep_inhibit function returns the correct value."""
    out = temp_dep_inhibit(
        P_ref=99.9,
        T_ref=99.9,
        H_a=99.9,
        H_d=99.9,
        S=99.9,
        T=99.9,
        R=99.9,
    )
    assert isclose(out, 99.9, abs_tol=1e-3)

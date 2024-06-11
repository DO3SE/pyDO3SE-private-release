"""Gsto related helper functions.

Note: Values are random unless specified
"""
from math import isclose
from .helpers import *


def test_calc_f_temp():
    """Tests the f_temp function returns the correct value."""
    out = calc_f_temp(
        Ts_C=22.3,
        T_min=1.3,
        T_opt=18,
        T_max=30.3,
        fmin=0.01,
    )
    assert isclose(out, 0.91602, abs_tol=1e-3)


def test_calc_f_temp_square_high():
    """Tests the calc_f_temp_square_high function returns the correct value."""
    out = calc_f_temp_square_high(
        Ts_C=22.3,
        T_min=1.3,
        T_opt=18,
        T_max=30.3,
        fmin=0.01,
    )
    assert isclose(out, 1.0, abs_tol=1e-3)


def test_calc_f_light():
    """Tests the calc_f_light function returns the correct value."""
    out = calc_f_light(
        f_lightfac=0.0105,
        PARsun=80.8,
        PARshade=30.7,
        LAIsunfrac=0.7,
    )
    assert isclose(out, 0.9167, abs_tol=1e-3)


def test_calc_leaf_f_light():
    """Tests the calc_leaf_f_light function returns the correct value."""
    out = calc_leaf_f_light(
        f_lightfac=0.0105,
        PAR=40.9,
    )
    assert isclose(out, 0.85950, abs_tol=1e-3)


def test_calc_f_light_method():
    out = calc_f_light_method(
        LAI=3.5,
        sinB=0.3,
        f_lightfac=0.0105,
        PARsun=80.3,
        PARshade=30.7,
        LAIsunfrac=0.7,
        PAR=40.9,
    )
    assert isclose(out.f_light, 0.9167, abs_tol=1e-3)
    assert isclose(out.leaf_f_light, 0.85950, abs_tol=1e-3)


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


def test_f_VPD_linear():
    """Tests the f_VPD_linear function returns the correct value."""
    out = f_VPD_linear(
        VPD=2,
        VPD_max=1.2,
        VPD_min=3.2,
        fmin=0.01,
    )
    assert isclose(out, 0.6040, abs_tol=1e-3)


def test_f_VPD_log():
    """Tests the f_VPD_log function returns the correct value."""
    out = f_VPD_log(
        VPD=2,
        fmin=0.01,
    )
    assert isclose(out, 0.58411, abs_tol=1e-3)


def test_inverse_f_VPD_linear():
    """Tests the inverse_f_VPD_linear function returns the correct value."""
    out = inverse_f_VPD_linear(
        f_VPD=12.3,
        VPD_max=99.9,
        VPD_min=1.0,
        fmin=1.5,
    )
    assert isclose(out, -2135.24000, abs_tol=1e-3)


def test_inverse_f_VPD_log():
    """Tests the inverse_f_VPD_log function returns the correct value."""
    out = inverse_f_VPD_log(
        f_VPD=1.3,
        fmin=12.3,
    )
    assert isclose(out, 0.606530, abs_tol=1e-3)


def test_f_SWP_exp():
    """Tests the f_SWP_exp function returns the correct value."""
    out = f_SWP_exp(
        a=2,
        b=3,
        fmin=99.9,
        SWP=99.9,
    )
    assert isclose(out, 1.0, abs_tol=1e-3)


def test_inverse_f_SWP_exp():
    """Tests the inverse_f_SWP_exp function returns the correct value."""
    out = inverse_f_SWP_exp(
        a=3,
        b=2,
        f_SWP=99.9,
    )
    assert isclose(out, -5.77061, abs_tol=1e-3)


def test_f_SWP_linear():
    """Tests the f_SWP_linear function returns the correct value."""
    out = f_SWP_linear(
        SWP_min=0,
        SWP_max=99.9,
        fmin=12.3,
        SWP=12.3,
    )
    assert isclose(out, 1.0, abs_tol=1e-3)


def test_f_PAW():
    """Tests the f_PAW function returns the correct value."""
    out = f_PAW(
        ASW_FC=99.9,
        fmin=12.3,
        ASW=99.9,
        ASW_max=50,
        ASW_min=0,
    )
    assert isclose(out, 1.0, abs_tol=1e-3)


def test_inverse_f_PAW():
    """Tests the inverse_f_PAW function returns the correct value."""
    out = inverse_f_PAW(
        ASW_FC=99.9,
        fmin=12.3,
        f_PAW=99.9,
        ASW_max=50,
        ASW_min=0,
    )
    assert isclose(out, -387.223008, abs_tol=1e-3)

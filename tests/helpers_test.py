"""Tests for met modules helpers."""

import numpy as np
import pytest

from do3se_met.helpers import (
    calc_humidity,
    calc_humidity_list,
    calc_vpd_daily_accumulation,
    calc_vpd_daily_accumulation_list,
)


def test_calc_humidity():
    """Test the output of calc_humidity."""
    with pytest.raises(Exception):
        calc_humidity(Ts_C_in=20, VPD_in=2, RH_in=3)
    with pytest.raises(Exception):
        calc_humidity(Ts_C_in=20, VPD_in=None, RH_in=None)

    out = calc_humidity(Ts_C_in=20, VPD_in=2, RH_in=None)
    assert out.esat == 2.3390469163992624
    assert out.RH == 0.14495088320895783
    assert out.eact == 0.33904691639926243
    assert out.VPD == 2

    out_b = calc_humidity(Ts_C_in=20, VPD_in=None, RH_in=2)
    assert out_b.esat == 2.3390469163992624
    assert out_b.RH == 2
    assert out_b.eact == 4.678093832798525
    assert out_b.VPD == -2.3390469163992624


def test_calc_humidity_list():
    """Test the output of calc_humidity_list."""
    Ts_C_list = [22, 23, 12, 33, 23, 14]
    VPD_list = [2, 3, 1, 2, 3, 3]
    RH_list = [0.2, 0.3, 0.1, 0.3, 0.2, 0.4]
    none_list = [None for i in range(6)]
    out = calc_humidity_list(Ts_C_list=Ts_C_list, VPD_list=VPD_list, RH_list=none_list)
    assert np.isclose(out.esat, [2.644796919516473, 2.8103575430330863, 1.4030231277532583, 5.031794864302203, 2.8103575430330863, 1.5991283056791965],
                      atol=1e-4).all()
    assert np.isclose(out.RH, [0.2437982722826054, -0.06747983274834189, 0.28725337436072235, 0.6025275167338613, -0.06747983274834189, -0.8760220736170462],
                      atol=1e-4).all()
    assert np.isclose(out.eact, [0.6447969195164731, -0.18964245696691373, 0.4030231277532583, 3.031794864302203, -0.18964245696691373, -1.4008716943208035],
                      atol=1e-4).all()
    assert out.VPD == VPD_list

    out = calc_humidity_list(Ts_C_list=Ts_C_list, VPD_list=none_list, RH_list=RH_list)
    assert np.isclose(out.esat, [2.644796919516473, 2.8103575430330863, 1.4030231277532583, 5.031794864302203, 2.8103575430330863, 1.5991283056791965],
                      atol=1e-4).all()
    assert out.RH == RH_list
    assert np.isclose(out.eact, [0.5289593839032947, 0.8431072629099259, 0.14030231277532584, 1.5095384592906609, 0.5620715086066173, 0.6396513222716786],
                      atol=1e-4).all()
    assert np.isclose(out.VPD, [2.1158375356131787, 1.9672502801231604, 1.2627208149779325, 3.522256405011542, 2.248286034426469, 0.9594769834075179],
                      atol=1e-4).all()


def test_calc_VPD_dd():
    """Test the output of calc_VPD_dd."""
    VPD_dd = calc_vpd_daily_accumulation(1.0, 0, False, 0)
    assert VPD_dd == 0
    VPD_dd = calc_vpd_daily_accumulation(1.0, 3.4, False, 3)
    assert VPD_dd == 3.4
    VPD_dd = calc_vpd_daily_accumulation(1.0, 3.4, True, 3)
    assert VPD_dd == 4.4


def test_calc_VPD_dd_list():
    """Test the output of calc_VPD_dd_list."""
    VPD_list = [1, 2, 3, 4, 5, 4]
    is_daylight_list = [False, False, True, True, True, False]
    hr_list = [0, 1, 2, 3, 0, 1, 2, 3]
    VPD_dd_list = calc_vpd_daily_accumulation_list(VPD_list, is_daylight_list, hr_list)
    assert VPD_dd_list == [0, 0, 3, 7, 5, 5]

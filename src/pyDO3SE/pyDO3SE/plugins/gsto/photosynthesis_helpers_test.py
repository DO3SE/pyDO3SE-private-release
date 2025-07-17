from math import isclose
from pyDO3SE.Config.ConfigEnums import FVPDMethods
from pyDO3SE.plugins.gsto.photosynthesis_helpers import calc_D_0, calc_D_0_f_VPD_Method


def test_f_VPD_method():
    D_0 = calc_D_0_f_VPD_Method(
        f_VPD_method=FVPDMethods.LINEAR,
        fmin=0.1,
        VPD_max=100,
        VPD_min=20,
    )
    assert isclose(D_0, 55.555, abs_tol=1e-3)
    D_0 = calc_D_0_f_VPD_Method(
        f_VPD_method=FVPDMethods.LOG,
        fmin=0.1,
        VPD_max=100,
        VPD_min=20,
    )
    assert isclose(D_0, 2.3009, abs_tol=1e-3)


def test_calc_D_0():
    D_0 = calc_D_0(
        D_0_method="constant",
        constant_D_0=12.3,
    )
    assert isclose(D_0, 12.3, abs_tol=1e-3)

    D_0 = calc_D_0(
        D_0_method="f_VPD",
        f_VPD_method=FVPDMethods.LINEAR,
        fmin=0.1,
        VPD_max=100,
        VPD_min=20,
    )

    assert isclose(D_0, 55.555, abs_tol=1e-3)

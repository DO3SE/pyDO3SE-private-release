from math import isclose

from pyDO3SE.plugins.soil_moisture.helpers import (
    SWC_to_SWP,
    SWP_to_SWC,
    get_soil_config, init_soilWater_state,
    soil_moisture_from_SWC,
    soil_moisture_from_SWP,
)
from pyDO3SE.plugins.soil_moisture.config import SOIL_LOAM, Soil_t


def test_check_soil():
    soil_config = get_soil_config('loam', None)
    assert soil_config == SOIL_LOAM
    soil_config = get_soil_config('custom', SOIL_LOAM)
    assert soil_config == SOIL_LOAM


def test_SWC_to_SWP():
    # TODO: Check realistic values
    SWP_AE = 0.3
    b = 0.2
    SWC = 0.3
    SWP = SWC_to_SWP(SWP_AE, b, SWC)
    assert isclose(SWP, 0.3177, abs_tol=1e-3)


def test_SWP_to_SWC():
    # TODO: Check realistic values
    SWP_AE = 0.3
    b = 0.2
    SWP = 0.3
    SWC = SWP_to_SWC(SWP_AE, b, SWP)
    assert isclose(SWC, 0.4, abs_tol=1e-3)


def test_soil_moisture_from_SWC():
    soil_config = Soil_t(0.2, 0.3, 0.4, 0.5)
    soil_moisture = soil_moisture_from_SWC(
        soil_config, PWP=0.3, root_depth=0.3, Sn_in=0.3)
    assert isclose(soil_moisture.Sn, 1.6855, abs_tol=1e-3)
    assert isclose(soil_moisture.SWP, 0.3, abs_tol=1e-3)
    assert isclose(soil_moisture.ASW, 0.0, abs_tol=1e-3)
    assert isclose(soil_moisture.SMD, -0.41567, abs_tol=1e-3)


def test_soil_moisture_from_SWP():
    soil_config = Soil_t(0.2, 0.3, 0.4, 0.5)
    soil_moisture = soil_moisture_from_SWP(
        soil_config, PWP=0.3, root_depth=0.3, SWP=0.9)
    assert isclose(soil_moisture.Sn, 1.6855, abs_tol=1e-3)
    assert isclose(soil_moisture.SWP, 0.3, abs_tol=1e-3)
    assert isclose(soil_moisture.ASW, 0.0, abs_tol=1e-3)
    assert isclose(soil_moisture.SMD, -0.41567, abs_tol=1e-3)


def test_soil_moisture_from_SWP_match_DO3SE_UI_row_26():
    soil_config = Soil_t(6.58, 0.29, -0.00188, 0.0002286)
    soil_moisture = soil_moisture_from_SWP(
        soil_config,
        PWP=-4,
        root_depth=0.75,
        SWP=-0.0156003293092588,
    )
    assert isclose(soil_moisture.Sn, 0.289270281792, abs_tol=1e-3)
    assert isclose(soil_moisture.SWP, -0.0158611088991, abs_tol=1e-3)
    assert isclose(soil_moisture.ASW, 0.123334348202, abs_tol=1e-3)
    assert isclose(soil_moisture.SMD, 0.000547282397747, abs_tol=1e-3)


def test_init_soilWater_state():
    out = init_soilWater_state(
        SWP_min=-1.25,
        SWP_AE=-0.00188,
        soil_b=6.58,
        Fc_m=0.29,
        root_depth=0.75,
        SWC_sat=0.4,
    )
    assert isclose(out.Sn, 0.29, abs_tol=1e-8)
    assert isclose(out.ASW, 0.12388163529178195, abs_tol=1e-8)
    assert isclose(out.SWP, -0.0156003293, abs_tol=1e-8)
    assert isclose(out.SMD, 0.0, abs_tol=1e-8)
    assert isclose(out.ASW_FC, 0.123881635291781, abs_tol=1e-8)

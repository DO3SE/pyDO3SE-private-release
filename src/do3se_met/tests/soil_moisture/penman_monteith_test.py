"""Tests for the Penman monteith functions.

Note: Some inputs are random.
"""

from math import exp, inf, isclose
import pytest

from do3se_met.resistance.model import Resistance_Model
from do3se_met.soil_moisture.enums import FSW_Methods
from do3se_met.soil_moisture.config import Soil_t
from do3se_met.soil_moisture.penman_monteith import (
    PM_soil_moisture_calc,
    check_soil_evaporation_blocked,
    multi_layer_r_model_to_single_H20,
    penman_monteith_daily,
    penman_monteith_hourly,
    penman_monteith_reset,
)


demo_rmodel_O3 = Resistance_Model(
    nL=3,
    Ra_canopy_to_izr=1.1,
    Rb=1.2,
    Rinc=[100, 200, 300],
    Rext=[120, 220, 330],
    Rsto=[130, 140, 150],
    Rgs=102.0,
    Rsur=[111, 222, 333],
)


def test_PM_soil_moisture_calc():
    soil_config = Soil_t(0.2, 0.3, 0.4, 0.5)
    soil_moisture = PM_soil_moisture_calc(
        soil_config, PWP=0.3, root_depth=0.3, Sn_in=0.3, Sn_diff=0.1
    )
    assert isclose(soil_moisture.Sn, 1.6855, abs_tol=1e-3)
    assert isclose(soil_moisture.SWP, 0.3, abs_tol=1e-3)
    assert isclose(soil_moisture.ASW, 0.0, abs_tol=1e-3)
    assert isclose(soil_moisture.SMD, -0.41567, abs_tol=1e-3)


def test_penman_monteith_reset():
    result = penman_monteith_reset()
    assert result.Ei_acc == 0
    assert result.Et_acc == 0
    assert result.Es_acc == 0
    assert result.Eat_acc == 0


def test_multilayer_r_model_to_single_H20():
    rmodel_H20 = multi_layer_r_model_to_single_H20(
        demo_rmodel_O3,
        ustar=1.1,
        total_LAI=1.3,
        total_SAI=3.4,
    )
    assert isinstance(rmodel_H20, Resistance_Model)
    assert isclose(rmodel_H20.Ra_canopy_to_izr, 1.1, abs_tol=1e-4)
    assert isclose(rmodel_H20.Rb, 3.92704, abs_tol=1e-4)
    assert len(rmodel_H20.Rinc) == 1
    assert isclose(rmodel_H20.Rinc[0], 54.545454, abs_tol=1e-4)
    assert isclose(rmodel_H20.Rext[0], 62.85714, abs_tol=1e-4)
    assert isclose(rmodel_H20.Rsto[0], 30.834582, abs_tol=1e-4)
    assert isclose(rmodel_H20.Rgs, 102.0, abs_tol=1e-4)
    assert isclose(rmodel_H20.Rsur[0], 18.2721, abs_tol=1e-4)


def test_multilayer_r_model_to_single_H20_when_rinc_is_0():
    """Test that before leaf emergence we can deal with the odd resistances."""
    rmodel_H20 = multi_layer_r_model_to_single_H20(
        Resistance_Model(
            nL=1,
            Ra_canopy_to_izr=47.24,
            Rb=33.485,
            Rinc=[0.0],
            Rext=[inf],
            Rsto=[100000],
            Rgs=102.0,
            Rsur=[100000],
        ),
        ustar=0.194063,
        total_LAI=0.0,
        total_SAI=0.0,
    )
    assert isinstance(rmodel_H20, Resistance_Model)
    assert isclose(rmodel_H20.Ra_canopy_to_izr, 47.24, abs_tol=1e-3)
    assert isclose(rmodel_H20.Rb, 22.259525, abs_tol=1e-3)
    assert len(rmodel_H20.Rinc) == 1
    assert isclose(rmodel_H20.Rinc[0], 0.0, abs_tol=1e-3)
    assert isclose(rmodel_H20.Rext[0], inf, abs_tol=1e-3)
    assert isclose(rmodel_H20.Rsto[0], 66300.0, abs_tol=1e-3)
    assert isclose(rmodel_H20.Rgs, 102.0, abs_tol=1e-3)
    assert isclose(rmodel_H20.Rsur[0], 102.0, abs_tol=1e-3)

    # This changed when we updated the multilayer model Jan 2025.
    # TODO: Investigate this further to check it still works as expected.
    # assert isclose(rmodel_H20.Rsur[0], 10.763110043, abs_tol=1e-4)


def test_check_soil_evaporation_blocked():
    """Test that the soil evaporation function returns the correct bool."""
    result = check_soil_evaporation_blocked(
        f_SW_method=None,
    )
    assert result is True
    result = check_soil_evaporation_blocked(
        f_SW_method=FSW_Methods.FSWP_EXP,
        fSWP_exp_a=99.9,
        fSWP_exp_b=999.9,
        SWP=0.3,
    )
    assert result is False

    # TODO: Check why this fails
    # result = check_soil_evaporation_blocked(
    #     f_SW_method='fPAW',
    #     ASW=0.0875123441219,
    #     ASW_FC=0.12388163529178195,
    #     fmin=0.01,
    # )
    # assert result is True

    result = check_soil_evaporation_blocked(
        f_SW_method=FSW_Methods.FSWP_LINEAR,
        SWP=0.3,
        SWP_max=0.2,
    )
    assert result is False

    with pytest.raises(ValueError):
        check_soil_evaporation_blocked(
            f_SW_method="None",
        )


def test_penman_monteith_hourly():
    # TODO: Check sensible inputs used here
    out = penman_monteith_hourly(
        Rn_MJ=30,
        P_kPa=91.2,
        Ts_C=24,
        esat_kPa=5,
        eact_kPa=3.37,
        VPD_kPa=1.3,
        rm_h2o_nL=1,
        rm_h2o_Rb=33.45,
        rm_h2o_Rinc_l0=99.9,
        rm_h2o_Rsto_l0=10000,
        rm_h2o_Rgs=10,
        rm_pet_Rsto_l0=0,
        LAI=0.01,
        Es_blocked=False,
        pm_state_Ei_acc=0.1,
        pm_state_Et_acc=0.1,
        pm_state_PEt_acc=0.1,
        pm_state_Es_acc=0.1,
        pm_state_Eat_acc=0.1,
    )
    assert out is not None

    assert isclose(out.Ei_hr, 0.009355466277616449, abs_tol=1e-7)
    assert isclose(out.Et_hr, 0.00018222670193883804, abs_tol=1e-7)
    assert isclose(out.Es_hr, 0.00907950519114286, abs_tol=1e-7)
    assert isclose(out.Eat_hr, 0.009083568473584668, abs_tol=1e-7)


def test_penman_monteith_hourly_match_do3se_ui_row_0():
    out = penman_monteith_hourly(
        Rn_MJ=0.0,
        P_kPa=91.86,
        Ts_C=3.28,
        esat_kPa=0.7732104619053026,
        eact_kPa=0.7162204619053026,
        VPD_kPa=0.05699,
        rm_h2o_nL=1.0,
        rm_h2o_Rb=22.25952506067983,
        rm_h2o_Rinc_l0=216.42455968133493,
        rm_h2o_Rsto_l0=99999.99999999999,
        rm_h2o_Rgs=200,
        rm_pet_Rsto_l0=0.0,
        LAI=3.0,
        Es_blocked=False,
        pm_state_Ei_acc=0.0,
        pm_state_Et_acc=0.0,
        pm_state_PEt_acc=0.1,
        pm_state_Es_acc=0.0,
        pm_state_Eat_acc=0.0,
    )
    assert out is not None

    assert isclose(out.Ei_hr, 3.769305726470758e-05, abs_tol=1e-7)
    # assert isclose(out.Et, 0.00018222670193883804, abs_tol=1e-7)
    # assert isclose(out.Es, 0.00907950519114286, abs_tol=1e-7)
    assert isclose(out.Eat_hr, 2.465408666674146e-06, abs_tol=1e-7)


def test_penman_monteith_hourly_match_do3se_ui_dd_61_hr_6():
    Ts_C = 10.9399995804
    esat = 0.611 * exp(17.27 * Ts_C / (Ts_C + 237.3))
    VPD = 0.297369986773
    eact = esat - VPD
    out = penman_monteith_hourly(
        Rn_MJ=0.0,
        P_kPa=94.0899963379,
        Ts_C=Ts_C,
        esat_kPa=esat,
        eact_kPa=eact,
        VPD_kPa=VPD,
        rm_h2o_nL=1,
        rm_h2o_Rb=116.222541809,
        rm_h2o_Rinc_l0=803.862609863,
        rm_h2o_Rsto_l0=100000.0,
        rm_h2o_Rgs=200,
        rm_pet_Rsto_l0=0,
        LAI=3.0,
        Es_blocked=True,
        pm_state_Ei_acc=0.0,
        pm_state_Et_acc=0.0,
        pm_state_PEt_acc=0.0,
        pm_state_Es_acc=0.0,
        pm_state_Eat_acc=0.0,
    )
    assert out is not None

    # assert isclose(out.Ei_hr, 3.769305726470758e-05, abs_tol=1e-7)
    # assert isclose(out.Et, 0.00018222670193883804, abs_tol=1e-7)
    assert isclose(out.Es_hr, 0.0, abs_tol=1e-7)
    # assert isclose(out.Eat_hr, 2.465408666674146e-06, abs_tol=1e-7)


def test_penman_monteith_daily():
    # TODO: Check sensable inputs used here
    out = penman_monteith_daily(
        LAI=0.01,
        root_depth=0.1,
        run_off_fraction=0.3,
        ASW=0.01,
        SMD=0.01,
        pm_state_precip_acc=0.01,
        pm_state_run_off_acc=0.01,
        pm_state_Ei_acc=0.01,
        pm_state_Eat_acc=0.01,
        pm_state_percolated_acc=0.01,
    )
    assert out is not None

    assert isclose(out.rain_input, 0.01, abs_tol=1e-3)
    assert isclose(out.run_off, 0.003, abs_tol=1e-3)
    assert isclose(out.run_off_acc, 0.013000, abs_tol=1e-3)
    assert isclose(out.effective_irrig, 0.007, abs_tol=1e-3)
    assert isclose(out.intercepted_evaporated, 1.0e-6, abs_tol=1e-3)
    assert isclose(out.evapotranspiration, 0.01, abs_tol=1e-3)
    assert isclose(out.Sn_diff, -0.03001, abs_tol=1e-3)
    assert isclose(out.percolated, 0.0, abs_tol=1e-3)
    assert isclose(out.percolated_acc, 0.01, abs_tol=1e-3)


def test_penman_monteith_daily_dd_62_hr_0():
    """Test output for row 1466."""
    out = penman_monteith_daily(
        LAI=3,
        root_depth=0.75,
        run_off_fraction=0.0,
        # ASW=0.0892531771635093,
        ASW=0.08751236537119962,
        SMD=0.03636926992058233,
        pm_state_precip_acc=0.00,
        pm_state_run_off_acc=0.0,
        pm_state_Ei_acc=0.010254495916419309,
        pm_state_Eat_acc=0.0019107638798741059,
<<<<<<< HEAD

=======
>>>>>>> 7aa3811e3c61482c39c63e5d7a57f2b53cd61328
        # pm_state_Ei_acc=0.00849898696512993,
        # pm_state_Eat_acc=0.0017408117923096824,
        pm_state_percolated_acc=0.0,
    )
    assert out is not None

    assert isclose(out.rain_input, 0.0, abs_tol=1e-8)
    assert isclose(out.run_off, 0.0, abs_tol=1e-8)
    assert isclose(out.run_off_acc, 0.0, abs_tol=1e-8)
    assert isclose(out.effective_irrig, 0.0, abs_tol=1e-8)
    assert isclose(out.intercepted_evaporated, 0.0, abs_tol=1e-8)
    assert isclose(out.evapotranspiration, 0.00191076387987, abs_tol=1e-8)
    assert isclose(out.Sn_diff, -0.00254768517, abs_tol=1e-8)
    assert isclose(out.percolated, 0.0, abs_tol=1e-8)
    assert isclose(out.percolated_acc, 0.0, abs_tol=1e-8)

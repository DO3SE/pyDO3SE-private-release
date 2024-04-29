"""Tests for Ewert helpers."""
import numpy as np
from math import isclose
from pyDO3SE.Config.ConfigEnums import FVPDMethods
from pyDO3SE.constants.physical_constants import DRATIO_O3_CO2
from pyDO3SE.plugins.O3.helpers import calc_O3up_accumulation

from .ewert_helpers import (
    Damage_Factors,
    calc_CO2_assimilation_rate,
    calc_CO2_supply,
    calc_humidity_defecit_fVPD,
    calc_input_factors,
    calc_mean_gsto,
    calc_ozone_damage_factors,
    calc_ozone_impact_on_lifespan,
    calc_senescence_factor,
    calc_stomatal_conductance,
)

t_lem_constant = 0.15
t_lse_constant = 0.33

t_l = 800
t_lem = t_l * t_lem_constant
t_lma = t_l - t_lem
t_lse = t_lma * t_lse_constant
t_lep = t_l - (t_lem + t_lse)  # check this


def test_calc_input_factors():
    """Test calc_input_factors output."""
    out = calc_input_factors(
        Tleaf_C=20.0,
        Q=800 * 4.57,
        V_cmax_25=180.0,
        J_max_25=400.0,
        R_d_coeff=0.015,
    )
    assert isclose(out.Tleaf_K, 293.15, abs_tol=1e-3)
    assert isclose(out.Gamma_star, 32.95310, abs_tol=1e-3)
    assert isclose(out.K_C, 234.425, abs_tol=1e-3)
    assert isclose(out.K_O, 216.752, abs_tol=1e-3)
    assert isclose(out.R_d, 1.796, abs_tol=1e-3)
    assert isclose(out.J_max, 311.687, abs_tol=1e-3)
    assert isclose(out.V_cmax, 119.785, abs_tol=1e-3)
    assert isclose(out.J, 300.359, abs_tol=1e-3)
    assert isclose(out.e_sat_i, 2339.047, abs_tol=1e-3)
    assert isclose(out.Gamma, 40.4835, abs_tol=1e-3)


def test_calc_ozone_damage_factors():
    """Test calc_ozone_damage_factors output."""
    out = calc_ozone_damage_factors(
        gamma_1=0.06,
        gamma_2=0.0045,
        gamma_3=0.5,
        cL3=0.2,
        O3up=88.8,
        O3up_acc=1250.0,
        fO3_d_prev=0.81,
        td_dd=770.0,
        t_lem=t_lem,
        t_lma=t_lma,
        hr=12,
        is_daylight=True,
    )

    assert isinstance(out, Damage_Factors)
    assert isclose(out.fO3_h, 0.660, abs_tol=1e-3)
    assert isclose(out.fO3_d, 0.534, abs_tol=1e-3)
    assert isclose(out.f_LA, 0.044, abs_tol=1e-3)
    assert isclose(out.fO3_l, 0.475, abs_tol=1e-3)
    assert isclose(out.rO3, 0.818, abs_tol=1e-3)


def test_calc_ozone_impact_on_lifespan():
    """Test calc_ozone_impact_on_lifespan output."""
    out = calc_ozone_impact_on_lifespan(
        t_lep=t_lep,
        t_lse=t_lse,
        t_lem=t_lem,
        fO3_l=0.375,
        gamma_4_senes=1,
        gamma_5_harvest=1,
    )
    assert isclose(out.t_lma_O3, 539.75, abs_tol=1e-3)
    assert isclose(out.t_lse_limited, 368.9, abs_tol=1e-3)
    assert isclose(out.t_lep_limited, 170.85, abs_tol=1e-3)


def test_calc_ozone_impact_on_lifespan_alt_senesence_method():
    """Test calc_ozone_impact_on_lifespan output.

    With this method t_lma should not change as ozone only effects onset of senesence
    t_lep should reduce bringing forward the date of onset of senesense.
    t_lse increase to make harvest date unchanged.
    """
    out = calc_ozone_impact_on_lifespan(
        t_lep=t_lep,
        t_lse=t_lse,
        t_lem=t_lem,
        fO3_l=0.375,
        gamma_4_senes=1,
        gamma_5_harvest=0,
    )
    assert isclose(out.t_lma_O3, t_lma, abs_tol=1e-3)
    assert isclose(out.t_lse_limited, 509.15, abs_tol=1e-3)
    assert isclose(out.t_lep_limited, 170.85, abs_tol=1e-3)


def test_calc_ozone_impact_on_lifespan_alt_senesence_method_2():
    """Test calc_ozone_impact_on_lifespan output.

    With this method t_lse should reduce but t_lep is uneffected
    """
    out = calc_ozone_impact_on_lifespan(
        t_lep=t_lep,
        t_lse=t_lse,
        t_lem=t_lem,
        fO3_l=0.375,
        gamma_4_senes=0,
        gamma_5_harvest=1,
    )
    assert isclose(out.t_lma_O3, 539.75, abs_tol=1e-3)
    assert isclose(out.t_lse_limited, 84.15, abs_tol=1e-3)
    assert isclose(out.t_lep_limited, t_lep, abs_tol=1e-3)

def test_calc_ozone_impact_on_lifespan_alt_senesence_method_3():
    """Test calc_ozone_impact_on_lifespan output.

    With this method fO3_l is disabled
    """
    out = calc_ozone_impact_on_lifespan(
        t_lep=t_lep,
        t_lse=t_lse,
        t_lem=t_lem,
        fO3_l=0.375,
        gamma_4_senes=0,
        gamma_5_harvest=0,
    )
    assert isclose(out.t_lma_O3, t_lma, abs_tol=1e-3)
    assert isclose(out.t_lse_limited, t_lse, abs_tol=1e-3)
    assert isclose(out.t_lep_limited, t_lep, abs_tol=1e-3)


def test_calc_ozone_impact_on_lifespan_high_ozone():
    """Test calc_ozone_impact_on_lifespan output with high ozone.

    This should shorten the t_lse and t_lep period
    """
    out = calc_ozone_impact_on_lifespan(
        t_lep=t_lep,
        t_lse=t_lse,
        t_lem=t_lem,
        fO3_l=0.8,
        gamma_4_senes=1,
        gamma_5_harvest=1,
    )
    assert isclose(out.t_lma_O3, 635.12, abs_tol=1e-3)
    assert isclose(out.t_lse_limited, 270.64, abs_tol=1e-3)
    assert isclose(out.t_lep_limited, 364.48, abs_tol=1e-3)


def test_calc_senescence_factor():
    """Test calc_senescence_factor output."""
    fO3_l = 0.9
    f_LS = calc_senescence_factor(
        td_dd=710,
        t_l_O3=t_lem + fO3_l * (t_lep + t_lse),
        # t_l_estimate=t_l,
        # fO3_l=0.375,
        t_lem=t_lem,
        t_lep_limited=t_lep * fO3_l,
        # t_lma=255.0,
        t_lse_limited=t_lse * fO3_l,
    )
    assert isclose(f_LS, 0.1089324, abs_tol=1e-3)


def test_full_hour():
    """Test full_hour output."""
    O3up = 88.8
    td_dd = 400
    O3up_acc = calc_O3up_accumulation(
        O3up=O3up,
        O3up_prev=O3up - 5,
        O3up_acc=232,
        td_dd=td_dd,
        td_dd_prev=td_dd - 10,
    )
    assert O3up_acc == 1095

    ozone_damage_factors = calc_ozone_damage_factors(
        gamma_1=0.06,
        gamma_2=0.0045,
        gamma_3=0.5,
        cL3=0.2,
        O3up=88.8,
        O3up_acc=O3up_acc,
        fO3_d_prev=0.81,
        td_dd=td_dd,
        t_lem=t_lem,
        t_lma=t_lma,
        is_daylight=True,
        hr=12,
    )
    new_lifespan = calc_ozone_impact_on_lifespan(
        t_lep=t_lep,
        t_lse=t_lse,
        t_lem=t_lem,
        fO3_l=ozone_damage_factors.fO3_l,
        gamma_4_senes=1,
        gamma_5_harvest=1,
    )

    f_LS = calc_senescence_factor(
        td_dd=td_dd,
        # t_l_estimate=t_l,
        # fO3_l=ozone_damage_factors.fO3_l,
        t_lem=t_lem,
        t_lep_limited=new_lifespan.t_lep_limited,
        # t_lma=new_lifespan.t_lma_O3,
        t_l_O3=new_lifespan.t_l_O3,
        t_lse_limited=new_lifespan.t_lse_limited,
    )
    assert isclose(f_LS, 0.913, abs_tol=1e-2)


def test_calc_CO2_assimilation_rate():
    """Test calc_CO2_assimilation_rate output."""
    out = calc_CO2_assimilation_rate(
        c_i_in=378.55,
        V_cmax=119.78,
        Gamma_star=32.95,
        K_C=234.42,
        K_O=216.75,
        fO3_d=0.534,
        f_LS=0.357,
        J=300.36,
        R_d=0.32,
    )

    assert isclose(out.A_c, 9.3938, abs_tol=1e-3)
    assert isclose(out.A_j, 58.3892, abs_tol=1e-3)
    assert isclose(out.A_p, 59.89, abs_tol=1e-3)
    assert isclose(out.A_n, 9.0738, abs_tol=1e-3)
    assert out.A_n_limit_factor == 'A_c'


def test_calc_stomatal_conductance():
    """Test calc_stomatal_conductance output."""
    out = calc_stomatal_conductance(
        g_sto_0=20000.0,
        m=8.12,
        Gamma=34.277,
        g_bl=1469999.9,
        c_a=391.0,
        A_n=12.865,
        f_SW=1,
        f_VPD=0.63,
    )

    assert isclose(out, 210907.634860, abs_tol=1e-3)


def test_calc_stomatal_conductance_An_0():
    """Test calc_stomatal_conductance output when A_n is 0."""
    out = calc_stomatal_conductance(
        g_sto_0=20000.0,
        m=8.12,
        Gamma=34.277,
        g_bl=1469999.9,
        c_a=391.0,
        A_n=0,
        f_SW=1,
        f_VPD=0.63,
    )

    assert isclose(out, 20000.0, abs_tol=1e-3)


def test_calc_stomatal_conductance_An_neg():
    """Test calc_stomatal_conductance output when A_n is 0."""
    out = calc_stomatal_conductance(
        g_sto_0=20000.0,
        m=8.12,
        Gamma=34.277,
        g_bl=1469999.9,
        c_a=391.0,
        A_n=-10,
        f_SW=1,
        f_VPD=0.63,
    )

    assert isclose(out, 20000, abs_tol=1e-3)


def test_calc_stomatal_conductance_f_SW():
    """Test calc_stomatal_conductance output when A_n is 0."""
    out = calc_stomatal_conductance(
        g_sto_0=20000.0,
        m=8.12,
        Gamma=34.277,
        g_bl=1469999.9,
        c_a=391.0,
        A_n=12.865,
        f_SW=0.5,
        f_VPD=0.63,
    )

    assert isclose(out, 115453.81743, abs_tol=1e-3)


def test_calc_humidity_defecit_fVPD():
    """Test calc_humidity_defecit_fVPD outputs correctly."""
    out = calc_humidity_defecit_fVPD(
        g_sto_in=20000.0,
        D_0=0.75,
        e_a=1000,
        g_bl=1469999.9,
        e_sat_i=2339.05,
        f_VPD_method=FVPDMethods.LEUNING,
    )
    assert isclose(out, 0.362130572, abs_tol=1e-3)


eact = 0.7


def test_calc_CO2_supply():
    """Test calc_CO2_supply output."""
    out = calc_CO2_supply(
        A_n=13.865,
        c_a=391.0,
        g_sto=145197.0,
        g_bl=1469999.9,
    )
    assert isclose(out, 282.5872427, abs_tol=1e-3)


def test_calc_CO2_supply_b():
    """Test calc_CO2_supply output."""
    out = calc_CO2_supply(
        A_n=13.865,
        c_a=391.0,
        g_sto=150 * (1 / DRATIO_O3_CO2) * 1000,
        g_bl=1469999.9,
    )
    assert isclose(out, 281.6702297332, abs_tol=1e-3)


# def test_calc_mean_gsto():
#     nL = 4
#     out = calc_mean_gsto(
#         [
#             [10, 10, 10, 10],  # pop 0
#             [20, 20, 20, 20],  # pop 1
#             [30, 30, 30, 30],  # pop2
#         ],
#         [
#             [1.0, 0.5, 0.2, 0.0],  # pop 0
#             [0.0, 0.5, 0.7, 0.0],  # pop 1
#             [0.0, 0.0, 0.1, 0.0],  # pop2
#         ],
#         nL,
#     )
#     assert len(out) == nL
#     assert isclose(out[0], 10.0, abs_tol=1e-3)
#     assert isclose(out[1], 15.0, abs_tol=1e-3)
#     assert isclose(out[2], 19.0, abs_tol=1e-3)
#     assert isclose(out[3], 0.0, abs_tol=1e-3)


def test_calc_mean_gsto():
    nL = 4
    leaf_lai_list = [
        # [pop0, pop1, pop2]
        [1.0, 0.0, 0.0],  # layer 0,
        [0.5, 0.5, 0.0],  # layer 1
        [0.2, 0.7, 0.1],  # layer 2
        [0.0, 0.0, 0.0],  # layer 3
    ]

    leaf_gsto = [
        # [pop0, pop1, pop2]
        [10.0, 20.0, 30.0],  # layer 0,
        [10.0, 20.0, 30.0],  # layer 1
        [10.0, 20.0, 30.0],  # layer 2
        [10.0, 20.0, 30.0],  # layer 3
    ]

    out = calc_mean_gsto(
        leaf_gsto,
        leaf_lai_list,
        nL,
    )
    assert len(out) == nL
    assert isclose(out[0], 10.0, abs_tol=1e-3)
    assert isclose(out[1], 15.0, abs_tol=1e-3)
    assert isclose(out[2], 19.0, abs_tol=1e-3)
    assert isclose(out[3], 0.0, abs_tol=1e-3)

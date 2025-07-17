"""Ewert model tests."""

import os
from dataclasses import asdict
from math import isclose
from pyDO3SE.Config.ConfigEnums import FVPDMethods
from do3se_met.enums import GAS

import pprint
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d  # noqa: F401
from tests.test_helpers import long_test
from thermal_time import calcs as td_functions


from pyDO3SE.constants.physical_constants import DRATIO_O3_CO2
from do3se_met.irradiance import (
    MLMC_sunlit_LAI,
    calc_Idrctt_Idfuse,
    calc_PAR_sun_shade,
)
from do3se_met.solar_position import calc_solar_elevation
from do3se_met.resistance import calc_g_bv


from .ewert import (
    CO2_Constant_Loop_Inputs,
    CO2_loop_State,
    CO2_Concentration_Args,
    ModelOptions,
    co2_concentration_in_stomata_iteration,
    co2_concentration_in_stomata_loop,
    ewert_leaf_pop,
    ewert_leaf_pop_cubic,
    co2_concentration_in_stomata_cubic,
)

matplotlib.use("Agg")

SAVE_OUTPUTS_LOCATION = os.path.dirname(os.path.realpath(__file__)) + "/output/"

t_lem_constant = 0.15
t_lse_constant = 0.33

t_l = 800
t_lem = t_l * t_lem_constant
t_lma = t_l - t_lem
t_lse = t_lma * t_lse_constant
t_lep = t_l - (t_lem + t_lse)  # check this


def test_co2_concentration_in_stomata_iteration(snapshot):
    """Test co2_concentration_in_stomata_iteration output."""
    state_out = co2_concentration_in_stomata_iteration(
        CO2_Constant_Loop_Inputs(
            c_a=391.0,
            e_a=1000.0,
            g_bl=1469999.0,
            g_sto_0=20000,
            m=8.12,
            D_0=2.27 * 1e3,
            Gamma=34.277,
            Gamma_star=32.95,
            V_cmax=119.0,
            K_C=234.42,
            K_O=216.75,
            J=300.36,
            R_d=0.32,
            e_sat_i=2339.05,
            f_SW=1.0,
            fO3_d=1.0,
            f_LS=1.0,
            f_VPD=1.0,
            fmin=0.1,
        ),
        CO2_loop_State(
            c_i=0.0,
            c_i_diff=0,
            g_sto=20000,
            A_n=0,
        ),
        ModelOptions(
            f_VPD_method=FVPDMethods.LEUNING,
            co2_concentration_max_iterations=5,
        ),
    )
    assert isinstance(state_out, CO2_loop_State)
    snapshot.assert_match(
        pprint.pformat(state_out, indent=4), "co2_concentration_in_stomata_iteration"
    )


def test_co2_concentration_in_stomata_iteration_b(snapshot):
    """Test co2_concentration_in_stomata_iteration_b output."""
    state_out = co2_concentration_in_stomata_loop(
        CO2_Constant_Loop_Inputs(
            c_a=391.0,
            e_a=1000.0,
            g_bl=1469999.0,
            g_sto_0=20000,
            m=8.12,
            D_0=2.27 * 1e3,
            Gamma=34.277,
            Gamma_star=32.95,
            V_cmax=119.0,
            K_C=234.42,
            K_O=216.75,
            J=300.36,
            R_d=0.32,
            e_sat_i=2339.05,
            f_SW=1,
            fO3_d=0.8,
            f_LS=0.91,
            f_VPD=1.0,
            fmin=0.1,
        ),
        ModelOptions(
            f_VPD_method=FVPDMethods.LEUNING,
            co2_concentration_max_iterations=5,
        ),
    )
    assert isinstance(state_out, CO2_loop_State)

    snapshot.assert_match(
        pprint.pformat(state_out, indent=4), "co2_concentration_in_stomata_iteration_b"
    )


# def test_ewert(snapshot):
#     """Test ewert output."""
#     LAI = 0.01

#     out = ewert(
#         nP=1,
#         Tleaf_C=20.0,
#         PARsun=800,
#         PARshade=800,
#         V_cmax_25=180.0,
#         J_max_25=400.0,
#         layer_LAI=LAI,
#         LAIsunfrac=0.9,

#         t_lse_constant=t_lse_constant,
#         t_l_estimate=[t_l],
#         t_lem=t_lem,
#         t_lep=t_lep,
#         t_lse=t_lse,
#         t_lma=t_lma,

#         c_a=391.0,
#         eact=1.0,
#         g_bv=1469999.0,
#         g_sto_0=20000,
#         m=8.12,
#         D_0=2.27,
#         O3up=0.2,
#         O3up_acc=20,
#         fO3_d_prev=[0.89],
#         td_dd=[24.1],
#         gamma_1=0.06,
#         gamma_2=0.0045,
#         gamma_3=0.5,
#         hr=0,
#         f_SW=1,
#         f_VPD=1.0,
#         leaf_pop_distribution=[LAI],
#         is_daylight=True,
#         use_O3_damage=True,
#         opt_full_night_recovery=True,
#         f_VPD_method=FVPDMethods.LEUNING,
#         R_d_coeff=0.015,
#     )

#     snapshot.assert_match(asdict(out))

#     g_sto = DRATIO_O3_CO2 * (max(0.0, out.g_sv) / 1000)
#     assert isclose(g_sto, 705.431, abs_tol=1e-3)

#     # Check outputs that should be per leaf population
#     assert type(out.g_sv) == float
#     assert type(out.f_LS) == type([])
#     assert type(out.f_LA) == type([])
#     assert type(out.A_n) == float
#     assert type(out.A_c) == float
#     assert type(out.A_j) == float
#     assert type(out.A_p) == float
#     assert type(out.A_n_limit_factor) == type([])
#     assert type(out.A_n_limit_factor[0]) == str
#     assert type(out.R_d) == float
#     assert type(out.fO3_h_out) == type([])
#     assert type(out.fO3_d_out) == type([])
#     assert type(out.fO3_l) == type([])
#     assert type(out.c_i) == float
#     assert type(out.t_l_O3) == type([])
#     assert type(out.f_VPD) == float
#     assert type(out.v_cmax) == float
#     assert type(out.j_max) == float


# def test_ewert_senescence(snapshot):
#     """Test ewert output."""
#     LAI = 0.01

#     out = ewert(
#         nP=1,
#         Tleaf_C=24.36,
#         PARsun=483.0,
#         PARshade=115.0,
#         V_cmax_25=180.0,
#         J_max_25=400.0,
#         layer_LAI=LAI,
#         LAIsunfrac=0.9,

#         t_lse_constant=t_lse_constant,
#         t_l_estimate=[t_l],
#         t_lem=t_lem,
#         t_lep=t_lep,
#         t_lse=t_lse,
#         t_lma=t_lma,

#         c_a=391.0,
#         eact=1.934,
#         g_bv=1927886.0,
#         g_sto_0=20000,
#         m=8.12,
#         D_0=2.27,
#         O3up=0.2,
#         O3up_acc=71,
#         fO3_d_prev=[0.89],
#         td_dd=[691.0],
#         gamma_1=0.06,
#         gamma_2=0.0045,
#         gamma_3=0.5,
#         hr=0,
#         f_SW=1,
#         f_VPD=1.0,
#         leaf_pop_distribution=[LAI],
#         is_daylight=True,
#         use_O3_damage=True,
#         opt_full_night_recovery=True,
#         f_VPD_method=FVPDMethods.LEUNING,
#         R_d_coeff=0.015,
#     )
#     snapshot.assert_match(asdict(out))
#     g_sto = DRATIO_O3_CO2 * (max(0.0, out.g_sv) / 1000)
#     assert isclose(g_sto, 283.526, abs_tol=1e-3)


# def test_ewert_f_sw(snapshot):
#     """Test ewert output."""
#     LAI = 0.01

#     out = ewert(
#         nP=1,
#         Tleaf_C=20.0,
#         PARsun=800,
#         PARshade=800,
#         V_cmax_25=180.0,
#         J_max_25=400.0,
#         layer_LAI=LAI,
#         LAIsunfrac=0.9,
#         leaf_pop_distribution=[LAI],
#         t_lse_constant=t_lse_constant,
#         t_l_estimate=[t_l],
#         t_lem=t_lem,
#         t_lep=t_lep,
#         t_lse=t_lse,
#         t_lma=t_lma,

#         c_a=391.0,
#         eact=1.0,
#         g_bv=1469999.0,
#         g_sto_0=20000,
#         m=8.12,
#         D_0=2.27,
#         O3up=0.2,
#         O3up_acc=20,
#         fO3_d_prev=[0.89],
#         td_dd=[24.1],
#         gamma_1=0.06,
#         gamma_2=0.0045,
#         gamma_3=0.5,
#         hr=0,
#         f_SW=0.5,
#         f_VPD=None,
#         is_daylight=True,
#         use_O3_damage=True,
#         opt_full_night_recovery=True,
#         f_VPD_method=FVPDMethods.LEUNING,
#         R_d_coeff=0.015,
#     )
#     snapshot.assert_match(asdict(out))
#     g_sto = DRATIO_O3_CO2 * (max(0.0, out.g_sv) / 1000)
#     assert isclose(g_sto, 293.295, abs_tol=1e-3)


# def test_ewert_at_lai_0(snapshot):
#     """Test ewert output before season start.

#     This is to ensure we don't get any NaN values.
#     """
#     LAI = 0.00

#     out = ewert(
#         nP=1,
#         Tleaf_C=20.0,
#         PARsun=800,
#         PARshade=800,
#         V_cmax_25=180.0,
#         J_max_25=400.0,
#         layer_LAI=LAI,
#         LAIsunfrac=0.9,
#         leaf_pop_distribution=[LAI],
#         t_lse_constant=t_lse_constant,
#         t_l_estimate=[t_l],
#         t_lem=t_lem,
#         t_lep=t_lep,
#         t_lse=t_lse,
#         t_lma=t_lma,

#         c_a=391.0,
#         eact=1.0,
#         g_bv=1469999.0,
#         g_sto_0=20000,
#         m=8.12,
#         D_0=2.27,
#         O3up=0.2,
#         O3up_acc=20,
#         fO3_d_prev=[0.89],
#         td_dd=[24.1],
#         gamma_1=0.06,
#         gamma_2=0.0045,
#         gamma_3=0.5,
#         hr=0,
#         f_SW=1,
#         f_VPD=1.0,
#         is_daylight=True,
#         use_O3_damage=True,
#         opt_full_night_recovery=True,
#         f_VPD_method=FVPDMethods.LEUNING,
#         R_d_coeff=0.015,
#     )

#     snapshot.assert_match(asdict(out))
#     g_sto = DRATIO_O3_CO2 * (max(0.0, out.g_sv) / 1000)
#     assert isclose(g_sto, 0, abs_tol=1e-3)


# @long_test
# def test_gsto_output_full_year(snapshot):
#     """Test full gsto_pn method for a full year."""
#     DEMO_DATA_LOCATION = './examples/spanish_wheat/data/spanish_wheat_data_test.csv'
#     SAVE_OUTPUTS_LOCATION = './pyDO3SE/plugins/gsto/ewert/output/'

#     # import data
#     df = pd.read_csv(DEMO_DATA_LOCATION)
#     data = {}
#     data["Ts_C"] = df['Ts_C'].values.reshape((365, 24))
#     data["u"] = df['u'].values.reshape((365, 24))
#     data["VPD"] = df['VPD'].values.reshape((365, 24))
#     data["PAR"] = df['PAR, W m-2'].values.reshape((365, 24))
#     data["P"] = df['P, kPa'].values.reshape((365, 24))
#     data["O3"] = df['O3, ppb'].values.reshape((365, 24))
#     data['fst'] = df['Fst (nmol/m^2/s)'].values.reshape((365, 24))
#     data['Lai'] = df['LAI'].values.reshape((365, 24))

#     td_data = np.array(td_functions.calc_thermal_time_range(df['Ts_C'].values))

#     # previous_day_output = None
#     output_data = []
#     input_data = []
#     D_0 = 2.27
#     day_count = 365
#     fO3_d_prev = [1]
#     O3up_prev = 0.0
#     O3up_acc = 0.0
#     td_dd_prev = 0.0
#     # config
#     season_Astart = 153
#     season_Astart_temp = td_functions.get_thermal_time_at_day(season_Astart, td_data, 0)
#     print(season_Astart_temp)
#     ewert_constant_inputs = {
#         "nP": 1,
#         "t_l_estimate": [t_l],
#         "t_lse_constant": t_lse_constant,
#         "t_lma": t_lma,
#         "t_lem": t_lem,
#         "t_lep": t_lep,
#         "t_lse": t_lse,
#         "gamma_1": 0.06,
#         "gamma_2": 0.0045,
#         "gamma_3": 0.5,
#         "g_sto_0": 20000,
#         "m": 8.12,
#         "V_cmax_25": 180.0,
#         "J_max_25": 400,
#         "use_O3_damage": True,
#         "opt_full_night_recovery": True,
#         "f_VPD_method": FVPDMethods.LEUNING,
#     }

#     for dd in range(day_count):
#         fo3_d_acc = 0
#         td = td_data[dd * 24]

#         for hr in range(24):
#             # O3_nmol = O3_ppb_to_nmol(data["Ts_C"][dd][hr], data["P"][dd][hr], data["O3"][dd][hr])
#             O3up = data['fst'][dd][hr]
#             td_dd = td_functions.get_td_dd(dd, td, season_Astart, season_Astart_temp)
#             g_bv = calc_g_bv(0.02, data["u"][dd][hr], GAS.H2O)
#             LAI = data['Lai'][dd][hr]
#             sinB = calc_solar_elevation(40.43, -3.7, dd, hr)
#             LAIsunfrac = MLMC_sunlit_LAI(1, 1, [[LAI]], sinB)[0][0]
#             PAR = data["PAR"][dd][hr]
#             P = data["P"][dd][hr]
#             Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
#             cosA = 0.5
#             PARsun, PARshade = calc_PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, LAI)
#             # PARsun = data["PAR"][dd][hr]
#             # PARshade = data["PAR"][dd][hr]
#             is_daylight = calc_is_daylight(PARsun)

#             ewert_inputs = {
#                 **ewert_constant_inputs,

#                 "PARsun": PARsun,
#                 "PARshade": PARshade,
#                 "layer_LAI": LAI,
#                 "LAIsunfrac": LAIsunfrac,
#                 "leaf_pop_distribution": [LAI],
#                 "D_0": D_0,
#                 "g_bv": g_bv,
#                 "td_dd": [td_dd],
#                 "fO3_d_prev": fO3_d_prev,
#                 "O3up_acc": O3up_acc,
#                 "Tleaf_C": data["Ts_C"][dd][hr],
#                 "O3up": O3up,
#                 "eact": data["VPD"][dd][hr],
#                 "c_a": 391.0,
#                 "is_daylight": is_daylight,
#                 "hr": hr,
#                 "f_SW": 1,
#                 "f_VPD": 1,
#                 "R_d_coeff": 0.015,
#             }
#             hr_output = ewert(**ewert_inputs)

#             O3up_acc = O3up_acc + ((O3up + O3up_prev) / 2 * (td_dd - td_dd_prev))
#             assert len(hr_output.f_LS) == 1
#             output_data.append(hr_output)
#             input_data.append(ewert_inputs)

#             fO3_d_prev = hr_output.fO3_d_out
#             td_dd_prev = td_dd
#             O3up_prev = O3up

#         fO3_d_prev = [fo3_d_acc]

#     A_n_full = np.array([i.A_n for i in output_data])
#     A_n_limit_factor_full = np.array([i.A_n_limit_factor for i in output_data])
#     A_c_full = np.array([i.A_c for i in output_data])
#     A_j_full = np.array([i.A_j for i in output_data])
#     A_p_full = np.array([i.A_p for i in output_data])
#     fO3_h_full = np.array([i.fO3_h_out for i in output_data])
#     fO3_d_full = np.array([i.fO3_d_out for i in output_data])
#     f_LS_full = np.array([i.f_LS for i in output_data])
#     R_d_full = np.array([i.R_d for i in output_data])
#     g_sto_full = np.array([DRATIO_O3_CO2 * (max(0.0, i.g_sv) / 1000) for i in output_data])
#     g_sv_full = np.array([i.g_sv for i in output_data])

#     df_out = pd.DataFrame([asdict(i) for i in output_data])
#     # .drop([
#     #     'Tleaf_C'
#     # ], axis=1)
#     df_in = pd.DataFrame(input_data)

#     # export data
#     df_out.to_csv(SAVE_OUTPUTS_LOCATION + 'gsto_output_full_year.csv')
#     df_in.to_csv(SAVE_OUTPUTS_LOCATION + 'gsto_output_full_year_inputs.csv')

#     # plots
#     x = np.array([[k for k in range(24)] for j in [i for i in range(day_count)]]).flatten()
#     y = np.array([[j for k in range(24)] for j in [i for i in range(day_count)]]).flatten()

#     fig = plt.figure()
#     fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
#     ax = fig.add_subplot(111, projection='3d')

#     z = g_sv_full
#     ax.scatter(x, y, z, c=z)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'g_sv.png')

#     ax.clear()
#     z = g_sto_full
#     ax.scatter(x, y, z, c=z)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'g_sto.png')

#     ax.clear()
#     z = [float(i[0]) for i in f_LS_full]
#     ax.scatter(x, y, z, c=z)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'f_LS.png')

#     ax.clear()
#     z = A_n_full
#     ax.scatter(x, y, z, c=z)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'A_n.png')

#     ax.clear()
#     z = A_n_full
#     c = ([0 if i[-1] == 'A_c' else 1 if i[-1] == 'A_j' else 2 if i[-1] ==  # noqaW504
#           'A_p' else 5 for i in A_n_limit_factor_full])
#     ax.scatter(x, y, z, c=c)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'A_n_limit.png')

#     ax.clear()
#     z = A_c_full
#     ax.scatter(x, y, z, c=z)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'A_c.png')

#     ax.clear()
#     z = fO3_d_full
#     ax.scatter(x, y, z, c=z)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'fO3_d.png')

#     ax.clear()
#     z = fO3_h_full
#     ax.scatter(x, y, z, c=z)
#     ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'fO3_h.png')

#     # print max and min
#     snapshot.assert_match('max: {0}, min: {1}'.format(A_n_full.max(), A_n_full.min()), 'A_n - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(A_c_full.max(), A_c_full.min()), 'A_c - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(A_j_full.max(), A_j_full.min()), 'A_j - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(A_p_full.max(), A_p_full.min()), 'A_p - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(fO3_h_full.max(), fO3_h_full.min()), 'fO3_h - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(fO3_d_full.max(), fO3_d_full.min()), 'fO3_d - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(R_d_full.max(), R_d_full.min()), 'R_d - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(g_sto_full.max(), g_sto_full.min()), 'g_sto - min/max')  # noqa: E501
#     snapshot.assert_match('max: {0}, min: {1}'.format(g_sv_full.max(), g_sv_full.min()), 'g_sv - min/max')  # noqa: E501

#     # print random sample
#     np.random.seed(0)

#     snapshot.assert_match(np.array2string(np.random.choice(A_n_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'A_n')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice(A_c_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'A_c')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice(A_j_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'A_j')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice(A_p_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'A_p')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice([i[0] for i in fO3_h_full], 30, replace=False), precision=5, separator=',', threshold=np.inf), 'fO3_h')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice([i[0] for i in fO3_d_full], 30, replace=False), precision=5, separator=',', threshold=np.inf), 'fO3_d')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice(R_d_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'R_d')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice(g_sto_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'g_sto')  # noqa: E501
#     snapshot.assert_match(np.array2string(np.random.choice(g_sv_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'g_sv')  # noqa: E501
#     # snapshot.assert_match(np.array2string(np.random.choice(g_bv_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'g_bv')  # noqa: E501


# class TestMultiplePopulations:

#     def test_should_run_with_multiple_populations(self, snapshot):
#         """Test ewert output."""
#         LAI = 0.01

#         out = ewert(
#             nP=3,
#             Tleaf_C=20.0,
#             PARsun=800,
#             PARshade=800,
#             V_cmax_25=180.0,
#             J_max_25=400.0,
#             layer_LAI=LAI,
#             LAIsunfrac=0.9,

#             c_a=391.0,
#             eact=1.0,
#             g_bv=1469999.0,
#             g_sto_0=20000,
#             m=8.12,
#             D_0=2.27,
#             f_SW=1,
#             f_VPD=1.0,
#             leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
#             f_VPD_method=FVPDMethods.LEUNING,
#             R_d_coeff=0.015,
#         )

#         snapshot.assert_match(asdict(out))

#         g_sto = DRATIO_O3_CO2 * (max(0.0, out.g_sv) / 1000)
#         assert isclose(g_sto, 505.388, abs_tol=1e-3)

#         # Check outputs that should be per leaf population
#         assert type(out.g_sv) == float
#         assert type(out.f_LS) == type([])
#         assert type(out.f_LA) == type([])
#         assert type(out.A_n) == float
#         assert type(out.A_c) == float
#         assert type(out.A_j) == float
#         assert type(out.A_p) == float
#         assert type(out.A_n_limit_factor) == type([])
#         assert type(out.A_n_limit_factor[0]) == str
#         assert type(out.R_d) == float
#         assert type(out.fO3_h_out) == type([])
#         assert type(out.fO3_d_out) == type([])
#         assert type(out.fO3_l) == type([])
#         assert type(out.c_i) == float
#         assert type(out.t_l_O3) == type([])
#         assert type(out.f_VPD) == float
#         assert type(out.v_cmax) == float
#         assert type(out.j_max) == float

#     def test_should_set_f_ls_correctly(self):
#         """Test ewert output."""
#         LAI = 0.01

#         out = [ewert(
#             nP=3,
#             Tleaf_C=20.0,
#             PARsun=800,
#             PARshade=800,
#             V_cmax_25=180.0,
#             J_max_25=400.0,
#             layer_LAI=LAI,
#             LAIsunfrac=0.9,


#             c_a=391.0,
#             eact=1.0,
#             g_bv=1469999.0,
#             g_sto_0=20000,
#             m=8.12,
#             D_0=2.27,
#             f_SW=1,
#             f_VPD=1.0,
#             leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
#             f_VPD_method=FVPDMethods.LEUNING,
#             R_d_coeff=0.015,
#         ) for i in range(2000)]


class TestEwertLeafPop:
    def test_should_run_with_multiple_layers(self, snapshot):
        """Test ewert output."""
        LAI = 0.01
        nL = 3
        out = ewert_leaf_pop(
            nL=nL,
            g_sto_0=20000,
            m=8.12,
            V_cmax_25=[180.0 for _ in range(nL)],
            J_max_25=[400.0 for _ in range(nL)],
            R_d_coeff=0.015,
            PARsun=[800 for _ in range(nL)],
            PARshade=[800 for _ in range(nL)],
            LAIsunfrac=[0.9, 0.8, 0.7],
            D_0=2.27,
            fmin=0.1,
            g_bv=[1469999.0 for _ in range(nL)],
            # leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
            layer_lai_frac=[0.33, 0.33, 0.33],
            layer_lai=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
            Tleaf_C=[20.0 for _ in range(nL)],
            f_SW=1,
            f_VPD=[1.0 for _ in range(nL)],
            f_LS=0.89,
            fO3_d=0.89,
            eact=1.0,
            c_a=391.0,
            f_VPD_method=FVPDMethods.LEUNING,
        )

        snapshot.assert_match(pprint.pformat(asdict(out), indent=4), "ewert_leaf_pop")

        g_sto = DRATIO_O3_CO2 * (max(0.0, out.g_sv_per_layer[0]) / 1000)
        assert isclose(g_sto, 616.307, abs_tol=1e-3)

        # Check outputs that should be per leaf population
        assert type(out.g_sv_per_layer) == list
        assert type(out.A_n) == float
        assert type(out.A_c) == float
        assert type(out.A_j) == float
        assert type(out.A_p) == float
        assert type(out.A_n_limit_factor) == type([])
        assert type(out.A_n_limit_factor[0]) == str
        assert type(out.R_d) == float
        assert type(out.c_i) == float
        assert type(out.f_VPD) == list
        assert type(out.v_cmax) == float
        assert type(out.j_max) == float

    @long_test
    def test_gsto_output_full_year(self, snapshot):
        """Test full gsto_pn method for a full year."""
        DEMO_DATA_LOCATION = "./examples/spanish_wheat/data/spanish_wheat_data_test.csv"
        SAVE_OUTPUTS_LOCATION = "./pyDO3SE/plugins/gsto/ewert/output/"
        nL = 3
        # import data
        df = pd.read_csv(DEMO_DATA_LOCATION)
        data = {}
        data["Ts_C"] = df["Ts_C"].values.reshape((365, 24))  # type: ignore
        data["u"] = df["u"].values.reshape((365, 24))  # type: ignore
        data["VPD"] = df["VPD"].values.reshape((365, 24))  # type: ignore
        data["PAR"] = df["PAR, W m-2"].values.reshape((365, 24))  # type: ignore
        data["P"] = df["P, kPa"].values.reshape((365, 24))  # type: ignore
        data["O3"] = df["O3, ppb"].values.reshape((365, 24))  # type: ignore
        data["fst"] = df["Fst (nmol/m^2/s)"].values.reshape((365, 24))  # type: ignore
        data["Lai"] = df["LAI"].values.reshape((365, 24))  # type: ignore

        td_data = np.array(td_functions.calc_thermal_time_range(df["Ts_C"].values))  # type: ignore

        # previous_day_output = None
        output_data = []
        input_data = []
        D_0 = 2.27
        day_count = 365
        # config
        season_Astart = 153

        ewert_constant_inputs = {
            "nL": nL,
            "g_sto_0": 20000,
            "m": 8.12,
            "V_cmax_25": [180.0 for _ in range(nL)],
            "J_max_25": [400 for _ in range(nL)],
            "R_d_coeff": 0.015,
            "f_VPD_method": FVPDMethods.LEUNING,
        }

        for dd in range(day_count):
            for hr in range(24):
                # O3_nmol = O3_ppb_to_nmol(data["Ts_C"][dd][hr], data["P"][dd][hr], data["O3"][dd][hr])
                fO3_d = 1.0  # Get from data
                f_LS = 1.0  # Get from data
                g_bv = [calc_g_bv(0.02, data["u"][dd][hr], GAS.H2O) for _ in range(nL)]
                LAI = [data["Lai"][dd][hr] * iL for iL in range(nL)]
                sinB = calc_solar_elevation(40.43, -3.7, dd, hr)
                LAIsunfrac = MLMC_sunlit_LAI(nL, 1, [[L] for L in LAI], sinB)
                LAIsunfrac = [f[0] for f in LAIsunfrac]
                PAR = data["PAR"][dd][hr]
                P = data["P"][dd][hr]
                Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
                cosA = 0.5
                PARsun, PARshade = list(
                    zip(*[calc_PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, L) for L in LAI])
                )
                gsto_prev = [
                    DRATIO_O3_CO2 * (max(0.0, output_data[-1].g_sv_per_layer[iL]) / 1000)
                    if len(output_data) > 0
                    else ewert_constant_inputs["g_sto_0"]
                    for iL in range(nL)
                ]
                ewert_inputs = {
                    **ewert_constant_inputs,
                    "PARsun": PARsun,
                    "PARshade": PARshade,
                    "LAIsunfrac": LAIsunfrac,
                    # "layer_lai_frac": np.mean(LAI, axis=0),
                    "layer_lai_frac": [lai/sum(LAI) for lai in LAI],
                    "layer_lai": LAI,
                    "D_0": D_0,
                    "g_bv": g_bv,
                    "P": data["P"][dd][hr],
                    "g_sto_prev": gsto_prev,
                    "Tleaf_C": [data["Ts_C"][dd][hr] for _ in range(nL)],
                    "eact": data["VPD"][dd][hr],
                    "c_a": 391.0,
                    "fO3_d": fO3_d,
                    "f_LS": f_LS,
                    "f_SW": 1,
                    "f_VPD": [1 for _ in range(nL)],
                    "fmin": 0.1,
                }
                hr_output = ewert_leaf_pop_cubic(**ewert_inputs)

                output_data.append(hr_output)
                input_data.append(ewert_inputs)

        A_n_full = np.array([i.A_n for i in output_data])
        A_n_limit_factor_full = np.array([i.A_n_limit_factor for i in output_data])
        A_c_full = np.array([i.A_c for i in output_data])
        A_j_full = np.array([i.A_j for i in output_data])
        A_p_full = np.array([i.A_p for i in output_data])
        R_d_full = np.array([i.R_d for i in output_data])
        g_sto_full = np.array(
            [DRATIO_O3_CO2 * (max(0.0, i.g_sv_per_layer[0]) / 1000) for i in output_data]
        )
        g_sv_full = np.array([i.g_sv_per_layer[0] for i in output_data])

        df_out = pd.DataFrame([asdict(i) for i in output_data])

        df_in = pd.DataFrame(input_data)

        # export data
        os.makedirs(SAVE_OUTPUTS_LOCATION, exist_ok=True)
        df_out.to_csv(SAVE_OUTPUTS_LOCATION + "gsto_output_full_year.csv")
        df_in.to_csv(SAVE_OUTPUTS_LOCATION + "gsto_output_full_year_inputs.csv")

        # plots
        x = np.array([[k for k in range(24)] for j in [i for i in range(day_count)]]).flatten()
        y = np.array([[j for k in range(24)] for j in [i for i in range(day_count)]]).flatten()

        fig = plt.figure()
        fig = plt.figure(figsize=(18, 16), dpi=80, facecolor="w", edgecolor="k")
        ax = fig.add_subplot(111, projection="3d")

        z = g_sv_full
        ax.scatter(x, y, z, c=z)
        ax.figure.savefig(SAVE_OUTPUTS_LOCATION + "g_sv.png")  # type: ignore

        ax.clear()
        z = g_sto_full
        ax.scatter(x, y, z, c=z)
        ax.figure.savefig(SAVE_OUTPUTS_LOCATION + "g_sto.png")  # type: ignore

        ax.clear()
        z = A_n_full
        ax.scatter(x, y, z, c=z)
        ax.figure.savefig(SAVE_OUTPUTS_LOCATION + "A_n.png")  # type: ignore

        ax.clear()
        z = A_n_full
        c = [
            0
            if i[-1] == "A_c"
            else 1
            if i[-1] == "A_j"
            else 2
            if i[-1]  # noqaW504
            == "A_p"
            else 5
            for i in A_n_limit_factor_full
        ]
        ax.scatter(x, y, z, c=c)
        ax.figure.savefig(SAVE_OUTPUTS_LOCATION + "A_n_limit.png")  # type: ignore

        ax.clear()
        z = A_c_full
        ax.scatter(x, y, z, c=z)
        ax.figure.savefig(SAVE_OUTPUTS_LOCATION + "A_c.png")  # type: ignore

        # print max and min
        snapshot.assert_match(
            "max: {0}, min: {1}".format(A_n_full.max(), A_n_full.min()), "A_n - min/max"
        )  # noqa: E501
        snapshot.assert_match(
            "max: {0}, min: {1}".format(A_c_full.max(), A_c_full.min()), "A_c - min/max"
        )  # noqa: E501
        snapshot.assert_match(
            "max: {0}, min: {1}".format(A_j_full.max(), A_j_full.min()), "A_j - min/max"
        )  # noqa: E501
        snapshot.assert_match(
            "max: {0}, min: {1}".format(A_p_full.max(), A_p_full.min()), "A_p - min/max"
        )  # noqa: E501
        snapshot.assert_match(
            "max: {0}, min: {1}".format(R_d_full.max(), R_d_full.min()), "R_d - min/max"
        )  # noqa: E501
        snapshot.assert_match(
            "max: {0}, min: {1}".format(g_sto_full.max(), g_sto_full.min()),
            "g_sto - min/max",
        )  # noqa: E501
        snapshot.assert_match(
            "max: {0}, min: {1}".format(g_sv_full.max(), g_sv_full.min()),
            "g_sv - min/max",
        )  # noqa: E501

        # print random sample
        np.random.seed(0)

        snapshot.assert_match(
            np.array2string(
                np.random.choice(A_n_full, 30, replace=False),
                precision=5,
                separator=",",
                threshold=np.inf,  # type: ignore
            ),
            "A_n",
        )  # noqa: E501
        snapshot.assert_match(
            np.array2string(
                np.random.choice(A_c_full, 30, replace=False),
                precision=5,
                separator=",",
                threshold=np.inf,  # type: ignore
            ),
            "A_c",
        )  # noqa: E501
        snapshot.assert_match(
            np.array2string(
                np.random.choice(A_j_full, 30, replace=False),
                precision=5,
                separator=",",
                threshold=np.inf,  # type: ignore
            ),
            "A_j",
        )  # noqa: E501
        snapshot.assert_match(
            np.array2string(
                np.random.choice(A_p_full, 30, replace=False),
                precision=5,
                separator=",",
                threshold=np.inf,  # type: ignore
            ),
            "A_p",
        )  # noqa: E501
        snapshot.assert_match(
            np.array2string(
                np.random.choice(R_d_full, 30, replace=False),
                precision=5,
                separator=",",
                threshold=np.inf,  # type: ignore
            ),
            "R_d",
        )  # noqa: E501
        snapshot.assert_match(
            np.array2string(
                np.random.choice(g_sto_full, 30, replace=False),
                precision=5,
                separator=",",
                threshold=np.inf,  # type: ignore
            ),
            "g_sto",
        )  # noqa: E501
        snapshot.assert_match(
            np.array2string(
                np.random.choice(g_sv_full, 30, replace=False),
                precision=5,
                separator=",",
                threshold=np.inf,  # type: ignore
            ),
            "g_sv",
        )  # noqa: E501
        # snapshot.assert_match(np.array2string(np.random.choice(g_bv_full, 30, replace=False), precision=5, separator=',', threshold=np.inf), 'g_bv')  # noqa: E501


def test_negative_A_n_values():
    R_d = 0.19323428021561032
    constant_inputs = CO2_Constant_Loop_Inputs(
        c_a=391.0,
        e_a=635.2238437279115,
        g_bl=3478654.912462574,
        g_sto_0=20000.0,
        m=6,
        D_0=2.7,
        Gamma=14.333665361703794,
        Gamma_star=12.52438053902872,
        V_cmax=12.882285347707354,
        K_C=30.750578109134757,
        K_O=85.49206630570801,
        J=40.30643977341539,
        R_d=R_d,
        e_sat_i=747.3221690916605,
        f_SW=1.0,
        f_LS=0,
        fO3_d=1.0,
        f_VPD=1.0,
        fmin=0.1,
    )
    model_options = ModelOptions(
        f_VPD_method=FVPDMethods.DISABLED,
        co2_concentration_balance_threshold=0.001,
        co2_concentration_max_iterations=50,
    )
    final_state = co2_concentration_in_stomata_loop(
        constant_inputs=constant_inputs,
        model_options=model_options,
    )
    assert final_state.A_n >= -R_d


def test_negative_A_n_values_using_cubic_method():
    R_d = 0.19323428021561032
    cubic_inputs = CO2_Concentration_Args(
        c_a=391.0,
        e_a=635.2238437279115,
        g_bl=3478654.912462574,
        g_sto_0=20000.0,
        g_sto_prev=20000.0,
        P=101.325,
        f_VPD=0.5,
        m=6,
        D_0=2.7,
        Gamma=14.333665361703794,
        Gamma_star=12.52438053902872,
        V_cmax=12.882285347707354,
        K_C=30.750578109134757,
        K_O=85.49206630570801,
        J=40.30643977341539,
        R_d=R_d,
        e_sat_i=747.3221690916605,
        f_SW=1.0,
        f_LS=0,
        fO3_d=1.0,
        fmin=0.1,
    )

    model_options = ModelOptions(
        f_VPD_method=FVPDMethods.DISABLED,
        co2_concentration_balance_threshold=0.001,
        co2_concentration_max_iterations=50,
    )
    final_state = co2_concentration_in_stomata_cubic(
        constant_inputs=cubic_inputs,
        model_options=model_options,
    )
    assert final_state.A_n >= -R_d


def test_negative_A_n_values_compare_methods():
    # OLD METHOD
    R_d = 0.19323428021561032
    constant_inputs = CO2_Constant_Loop_Inputs(
        c_a=391.0,
        e_a=3000.2238437279115,
        g_bl=1478654.912462574,
        g_sto_0=10000.0,
        m=6,
        D_0=2.7,
        Gamma=14.333665361703794,
        Gamma_star=12.52438053902872,
        V_cmax=82.882285347707354,
        K_C=30.750578109134757,
        K_O=85.49206630570801,
        J=40.30643977341539,
        R_d=R_d,
        e_sat_i=747.3221690916605,
        f_SW=1.0,
        f_LS=1.0,
        fO3_d=1.0,
        f_VPD=1.0,
        fmin=0.1,
    )

    cubic_inputs = CO2_Concentration_Args(
        c_a=391.0,
        e_a=3000.2238437279115,
        g_bl=1478654.912462574,
        g_sto_0=10000.0,
        g_sto_prev=40000.0,
        P=101.325,
        f_VPD=0.5,
        m=6,
        D_0=2.7,
        Gamma=14.333665361703794,
        Gamma_star=12.52438053902872,
        V_cmax=82.882285347707354,
        K_C=30.750578109134757,
        K_O=85.49206630570801,
        J=40.30643977341539,
        R_d=R_d,
        e_sat_i=747.3221690916605,
        f_SW=1.0,
        f_LS=1.0,
        fO3_d=1.0,
        fmin=0.1,
    )
    model_options = ModelOptions(
        f_VPD_method=FVPDMethods.DANIELSSON,
        co2_concentration_balance_threshold=0.001,
        co2_concentration_max_iterations=50,
    )
    final_state = co2_concentration_in_stomata_cubic(
        constant_inputs=cubic_inputs,
        model_options=model_options,
    )

    # NEW METHOD
    final_state_cubic = co2_concentration_in_stomata_loop(
        constant_inputs=constant_inputs,
        model_options=model_options,
    )

    assert final_state_cubic.A_n > -R_d
    print(final_state_cubic)
    print(final_state)
    # assert final_state_cubic.iterations != final_state.iterations
    assert isclose(final_state_cubic.A_p, final_state.A_p, rel_tol=0.1)
    assert isclose(final_state_cubic.A_c, final_state.A_c, rel_tol=0.1)
    assert isclose(final_state_cubic.A_j, final_state.A_j, rel_tol=0.1)
    assert isclose(final_state_cubic.A_n, final_state.A_n, rel_tol=0.1)
    assert final_state_cubic.A_n_limit_factor == final_state.A_n_limit_factor
    # assert isclose(final_state_cubic.A_p, final_state.A_p, abs_tol=1e-1)  False = isclose(6.441142673853677, 6.247908393638067, abs_tol=0.1)
    assert isclose(final_state_cubic.g_sto, final_state.g_sto, rel_tol=0.1)
    assert isclose(final_state_cubic.c_i, final_state.c_i, rel_tol=0.1)
    print("iterations", final_state.iterations, final_state_cubic.iterations)
    print("A_n_limiting_factor", final_state.A_n_limit_factor, final_state_cubic.A_n_limit_factor)
    print("A_n", final_state.A_n, final_state_cubic.A_n)
    print("A_j", final_state.A_j, final_state_cubic.A_j)
    print("A_p", final_state.A_p, final_state_cubic.A_p)
    print("A_c", final_state.A_c, final_state_cubic.A_c)
    print("cubic", final_state_cubic)
    print("iterative", final_state)


class TestEwertCubic:
    def test_outputs_f_vpd(self):
        nL = 3
        out = ewert_leaf_pop_cubic(
            nL=nL,
            g_sto_0=20000,
            m=8.12,
            V_cmax_25=[180.0 for _ in range(3)],
            J_max_25=[400.0 for _ in range(3)],
            R_d_coeff=0.015,
            PARsun=[800 for _ in range(3)],
            PARshade=[800 for _ in range(3)],
            LAIsunfrac=[0.9, 0.8, 0.7],
            layer_lai_frac=[0.33, 0.33, 0.33],
            layer_lai=[0.01 * 0.33, 0.01 * 0.33, 0.01 * 0.33],
            Tleaf_C=[20.0 for _ in range(3)],
            g_sto_prev=[20000.0 for _ in range(3)],
            f_VPD=[0.5 for _ in range(nL)],
            P=101.325,
            D_0=2.27,
            fmin=0.1,
            g_bv=[1469999.0 for _ in range(nL)],
            f_SW=1,
            f_LS=0.89,
            fO3_d=0.89,
            eact=1.0,
            c_a=391.0,
            f_VPD_method=FVPDMethods.LEUNING,
            co2_concentration_balance_threshold=0.001,
            co2_concentration_max_iterations=50,
        )
        assert len(out.f_VPD) == nL
        assert 0 < out.f_VPD[0] < 1
        assert 0 < out.f_VPD[1] < 1
        assert 0 < out.f_VPD[2] < 1

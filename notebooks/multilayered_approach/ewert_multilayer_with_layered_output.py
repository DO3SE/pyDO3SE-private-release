"""Experiment with ewert outputs per layer

TODO: We currently assume constant O3 across layers
We need to implement O3 deposition and also store the accumulated fst per layer

"""
# %%
import os
from dataclasses import asdict
from math import isclose

from numpy.lib.function_base import average
from pyDO3SE.util.Objects import Field
from pyDO3SE.Analysis.charts import multi_series_annual_graph

from helpers.list_helpers import flatten_list

import pytest
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d  # noqa: F401

from pyDO3SE.constants.physical_constants import DRATIO
from pyDO3SE.plugins.met.irradiance import (
    MLMC_sunlit_LAI, calc_Idrctt_Idfuse,
    calc_PAR_sun_shade,
    calc_is_daylight,
)
from pyDO3SE.plugins.met.solar_position import calc_solar_elevation
from pyDO3SE.plugins.met.wind import calc_windspeed_parameters

from pyDO3SE.plugins.met.thermal_time import get_season_astart_temperature, calc_thermal_time
from pyDO3SE.plugins.gsto.photosynthesis_helpers import calc_g_bv
from pyDO3SE.plugins.met.thermal_time import get_td_dd

from pyDO3SE.plugins.O3.helpers import O3_ppb_to_nmol

from pyDO3SE.plugins.gsto.ewert.ewert import (
    Output_Shape, ewert,
)

# %%

SAVE_OUTPUTS_LOCATION = './notebooks/multilayered_approach/output/'
DEMO_DATA_LOCATION = './examples/spanish_wheat/spanish_wheat_data_test.csv'

t_lem_constant = 0.15
t_lse_constant = 0.33

t_l = 800
t_lem = t_l * t_lem_constant
t_lma = t_l - t_lem
t_lse = t_l - (t_lma * t_lse_constant)
t_lep = t_l - (t_lem + t_lse)  # check this


# fractions = [0, 0.1, 0.6, 1.0, 2, 10, 100]
fractions = [0, 0.1, 0.6, 1.0, 2]
variation_count = len(fractions)


def get_ozone_fac(i, count):
    return fractions[i]


def test_gsto_output_full_year_multilayer(O3frac, nL):
    """Test full gsto_pn method for a full year."""

    # import data
    df = pd.read_csv(DEMO_DATA_LOCATION)
    data = {}
    data["Ts_C"] = df['Ts_C'].values.reshape((365, 24))
    data["u"] = df['u'].values.reshape((365, 24))
    data["VPD"] = df['VPD'].values.reshape((365, 24))
    data["PAR"] = df['PAR, W m-2'].values.reshape((365, 24))
    data["P"] = df['P, kPa'].values.reshape((365, 24))
    df['O3, ppb'] = df['O3, ppb'] * O3frac
    df['Fst (nmol/m^2/s)'] = df['Fst (nmol/m^2/s)'] * O3frac
    # data["O3"] = df['O3, ppb'].values.reshape((365, 24))
    data['fst'] = df['Fst (nmol/m^2/s)'].values.reshape((365, 24))
    data['Lai'] = df['LAI'].values.reshape((365, 24))
    td_data = np.array(calc_thermal_time(df['Ts_C'].values))

    # previous_day_output = None
    output_data = [[] for iL in range(nL)]
    input_data = []
    D_0 = 2.27
    day_count = 365
    fO3_d_prev = 1
    O3up_prev = 0.0
    O3up_acc = 0.0
    td_dd_prev = 0.0
    # config
    season_Astart = 153
    season_Astart_temp = get_season_astart_temperature(season_Astart, td_data)
    ewert_constant_inputs = {
        "t_l_estimate": t_l,
        "t_lse_constant": t_lse_constant,
        "t_lma": t_lma,
        "t_lem": t_lem,
        "t_lep": t_lep,
        "t_lse": t_lse,
        "gamma_1": 0.06,
        "gamma_2": 0.0045,
        "gamma_3": 0.5,
        "g_sto_0": 20000,
        "m": 8.12,
        "V_cmax_25": 180.0,
        "J_max_25": 400,
    }

    for dd in range(day_count):
        # if dd != 50: continue
        # fo3_d_acc = 0
        td = td_data[dd]

        for hr in range(24):
            # if hr != 11: continue
            O3up = data['fst'][dd][hr]

            td_dd = get_td_dd(dd, td, season_Astart, season_Astart_temp)
            u = data["u"][dd][hr] + 0.01  # TODO: Make per layer
            g_bv = calc_g_bv(0.02, u)
            LAI = [10 / nL for i in range(nL)]
            # LAI = [data['Lai'][dd][hr] / nL for i in range(nL)]
            LAI_c = [sum(LAI[0:i]) for i in range(nL)]
            # print(f'LAI: {LAI} LAI_c {LAI_c}')
            sinB = calc_solar_elevation(40.43, -3.7, dd, hr)
            LAIsunfrac = MLMC_sunlit_LAI(nL, 1, [[lai] for lai in LAI], sinB)
            PAR = data["PAR"][dd][hr]
            P = data["P"][dd][hr]
            Idrctt, Idfuse = calc_Idrctt_Idfuse(PAR, sinB, P)
            cosA = 0.5
            # PARsun shade should input LAI_c for layer
            # print(f'PAR: {PAR} Idrctt :{Idrctt}, Idfuse: {Idfuse}, sinB: {sinB}, cosA: {cosA}, LAI_c: {LAI_c}')
            PARSunShade = [calc_PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, lai_c) for lai_c in LAI_c]
            is_daylight = calc_is_daylight(PAR)
            # print(PARSunShade)
            # print(PARSunShade)

            # For each pass run for each layer.
            # LAI varies per layer
            # O3 varies per layer
            # PARsun shade varies per layer
            # Wind varies per layer
            hr_layer_outputs = [{} for iL in range(nL)]

            for iL in range(nL):
                PARsun, PARshade = PARSunShade[iL]
                ewert_inputs = {
                    **ewert_constant_inputs,
                    "PARsun": PARsun,
                    "PARshade": PARshade,
                    "LAI": LAI[iL],
                    "LAIsunfrac": LAIsunfrac[iL][0],
                    "D_0": D_0,
                    "g_bv": g_bv,
                    "td_dd": td_dd,
                    "fO3_d_prev": fO3_d_prev,  # This should be per layer
                    "O3up_acc": O3up_acc,  # This should be per layer
                    "Tleaf_C": data["Ts_C"][dd][hr],
                    "O3up": O3up,  # This should be per layer
                    "eact": data["VPD"][dd][hr],
                    "sinB": sinB,
                    "c_a": 391.0,
                    "is_daylight": is_daylight,
                }
                from pprint import pprint
                hr_layer_outputs[iL] = ewert(**ewert_inputs)

            # TODO: Upscale layer outputs to canopy output

            hr_output_canopy = Output_Shape(
                Tleaf_C=average([o.Tleaf_C for o in hr_layer_outputs]),
                g_sv=sum([o.g_sv for o in hr_layer_outputs]),
                f_LS=average([o.f_LS for o in hr_layer_outputs]),
                A_n=sum([o.A_n for o in hr_layer_outputs]),
                A_c=sum([o.A_c for o in hr_layer_outputs]),
                A_j=sum([o.A_j for o in hr_layer_outputs]),
                A_p=sum([o.A_p for o in hr_layer_outputs]),
                A_n_limit_factor=hr_layer_outputs[0].A_n_limit_factor,
                R_d=average([o.R_d for o in hr_layer_outputs]),
                fO3_h_out=average([o.fO3_h_out for o in hr_layer_outputs]),
                fO3_d_out=average([o.fO3_d_out for o in hr_layer_outputs]),
                fO3_l=average([o.fO3_l for o in hr_layer_outputs]),
                c_i=sum([o.c_i for o in hr_layer_outputs]),
            )
            # hr_output_layers = [Output_Shape(**hr_layer_outputs[iL]) for iL in range(nL)]

            O3up_acc = O3up_acc + ((O3up + O3up_prev) / 2 * (td_dd - td_dd_prev))

            [output_data[iL].append(hr_layer_outputs[iL]) for iL in range(nL)]
            # input_data.append(ewert_inputs)

            fO3_d_prev = hr_output_canopy.fO3_d_out
            td_dd_prev = td_dd
            O3up_prev = O3up

        # fO3_d_prev = fo3_d_acc
    return output_data


# %%
nL = 3
output_data_multilayer = test_gsto_output_full_year_multilayer(1, nL)
output_data_singlelayer = test_gsto_output_full_year_multilayer(1, 1)

# %%
output_data_multilayer = np.array(output_data_multilayer)
print(output_data_multilayer.shape)
output_data_multilayer[0][0]

# %%
plt.plot([i.c_i for i in output_data_multilayer[0]])
plt.plot([i.c_i for i in output_data_multilayer[1]])
plt.plot([i.c_i for i in output_data_multilayer[2]])

# %%
plt.plot([i.g_sv for i in output_data_multilayer[0]])
plt.plot([i.g_sv for i in output_data_multilayer[1]])
plt.plot([i.g_sv for i in output_data_multilayer[2]])

# %%
#


def get_canopy_data(field_data, nL):
    canopy = [sum([field_data[iL][i] for iL in range(nL)]) for i, d in enumerate(field_data[0])]
    return canopy


get_canopy_data([[1, 2, 3], [1, 2, 3], [1, 2, 3]], 3)
# %%
fields = ['g_sv', 'c_i', 'A_n']
for f in fields:
    field_data = [[getattr(i, f) for i in output_data_multilayer[iL]] for iL in range(nL)]
    canopy = [sum([field_data[iL][i] for iL in range(nL)]) for i, d in enumerate(field_data[0])]
    print(get_canopy_data(field_data, nL)[100])
    print(np.array(field_data)[:, 0])
    print(canopy[100])
    # field_data.append(canopy)
    field_data.append(getattr(i, f) for i in output_data_singlelayer[0])
    field_data = np.array([list(i) for i in field_data])

    print(field_data.shape)
    fig = multi_series_annual_graph(
        f,
        field_data,
        [i for i in range(nL)] + ['single'],
        # [i for i in range(nL)] + ['canopy'] + ['single'],
        # [i for i in range(nL)],
        None,
        label_x_days=30,
        average_step=24,
        figsize=(8, 3),
        dpi=300,
        chart_id=f,
        output_dir=SAVE_OUTPUTS_LOCATION,
    )

    print(np.max(field_data))


# %%
output_data_singlelayer = np.array(output_data_singlelayer)
print(output_data_singlelayer.shape)

# %%
plt.plot([i.c_i for i in output_data_singlelayer[0]])

# %%
plt.plot([i.g_sv for i in output_data_singlelayer[0]])


# %%
# for i in range(variation_count):
#     o3frac = get_ozone_fac(i, variation_count)
#     print(o3frac)
#     test_gsto_output_full_year_multilayer(1, 3)


# %%
# loaded_data = pd.read_csv(SAVE_OUTPUTS_LOCATION + f'gsto_output_full_year_{0}.csv')

# data_all = {f: [] for f in loaded_data.columns}
# for i in range(variation_count):
#     o3frac = get_ozone_fac(i, variation_count)
#     loaded_data = pd.read_csv(SAVE_OUTPUTS_LOCATION + f'gsto_output_full_year_{o3frac}.csv')
#     # print(loaded_data.head())
#     for f in loaded_data.columns:
#         data_all[f].append(loaded_data[f].values)

# fields = ['Tleaf_C', 'g_sv', 'f_LS', 'A_n', 'A_c', 'A_j', 'A_p', '',
#           'R_d', 'fO3_h_out', 'fO3_d_out', 'fO3_l', 'c_i', 'Canopy_A_n']

# start = 150
# end = 220
# for f in loaded_data.columns:
#     if f not in fields:
#         continue
#     series_titles = ['{:.0f}%'.format(100 * get_ozone_fac(i, variation_count))
#                      for i in range(variation_count)]

#     multi_series_annual_graph(
#         f,
#         [x[start * 24:end * 24] for x in data_all[f]],
#         series_titles,
#         Field(f, float),
#         output_dir=SAVE_OUTPUTS_LOCATION,
#         label_x_days=30,
#         average_step=24,
#         chart_id=f'{f}-comparison',
#         linestyles=[{
#             'style': 'dashed' if s == 'DO3SE_ui' else 'solid',
#             'width': 0.8,
#             'zorder': 100 if s == 'DO3SE_ui' else i,
#         } for i, s in enumerate(series_titles)],
#         figsize=(8, 3),
#         dpi=300,
#         output_format='svg',
#         start_day=start,
#         end_day=end,
#     )

# # %%
# print([x[0 * 24:1 * 24] for x in data_all['f_LS']])

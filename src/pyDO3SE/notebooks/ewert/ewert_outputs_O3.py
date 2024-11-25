"""Experiment with ewert outputs when we vary O3."""
# %%
import os
from dataclasses import asdict
from math import isclose
from pyDO3SE.util.Objects import Field
from pyDO3SE.Analysis.charts import multi_series_annual_graph

import pytest
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d  # noqa: F401


from pyDO3SE.constants.physical_constants import DRATIO
from pyDO3SE.plugins.met.irradiance import (
    calc_Idrctt_Idfuse,
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
    ewert,
)

# %%

SAVE_OUTPUTS_LOCATION = os.path.dirname(os.path.realpath(__file__)) + '/output/'
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


def test_gsto_output_full_year(O3frac):
    """Test full gsto_pn method for a full year."""
    if os.environ.get('TQUICK', 'False') == 'True':
        pytest.skip('Skip test in TQUICK mode as it takes a while to run')

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
    print(data['fst'][0])
    td_data = np.array(calc_thermal_time(df['Ts_C'].values))

    # previous_day_output = None
    output_data = []
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
        # fo3_d_acc = 0
        td = td_data[dd]

        for hr in range(24):
            # O3_nmol = O3_ppb_to_nmol(data["Ts_C"][dd][hr], data["P"][dd][hr], data["O3"][dd][hr])
            O3up = data['fst'][dd][hr]
            td_dd = get_td_dd(dd, td, season_Astart, season_Astart_temp)
            # u = data["u"][dd][hr] + 0.01
            g_bv = calc_g_bv(0.02, data["u"][dd][hr])
            LAI = data['Lai'][dd][hr]
            sinB = calc_solar_elevation(40.43, -3.7, dd, hr)
            PAR = data["PAR"][dd][hr]
            P = data["P"][dd][hr]
            Idrctt, Idfuse = calc_Idrctt_Idfuse(PAR, sinB, P)
            cosA = 0.5
            PARsun, PARshade = calc_PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, LAI)
            # PARsun = data["PAR"][dd][hr]
            # PARshade = data["PAR"][dd][hr]
            is_daylight = calc_is_daylight(PARsun)

            ewert_inputs = {
                **ewert_constant_inputs,
                "PARsun": PARsun,
                "PARshade": PARshade,
                "LAI": LAI,
                "D_0": D_0,
                "g_bv": g_bv,
                "td_dd": td_dd,
                "fO3_d_prev": fO3_d_prev,
                "O3up_acc": O3up_acc,
                "Tleaf_C": data["Ts_C"][dd][hr],
                "O3up": O3up,
                "eact": data["VPD"][dd][hr],
                "sinB": sinB,
                "c_a": 391.0,
                "is_daylight": is_daylight,
            }
            hr_output = ewert(**ewert_inputs)

            O3up_acc = O3up_acc + ((O3up + O3up_prev) / 2 * (td_dd - td_dd_prev))
            output_data.append(hr_output)
            input_data.append(ewert_inputs)

            fO3_d_prev = hr_output.fO3_d_out
            td_dd_prev = td_dd
            O3up_prev = O3up

        # fO3_d_prev = fo3_d_acc

    df_out = pd.DataFrame([asdict(i) for i in output_data])
    # .drop([
    #     'Tleaf_C'
    # ], axis=1)
    df_in = pd.DataFrame(input_data)

    # export data
    df_out.to_csv(SAVE_OUTPUTS_LOCATION + f'gsto_output_full_year_{O3frac}.csv')
    df_in.to_csv(SAVE_OUTPUTS_LOCATION + f'gsto_output_full_year_inputs.csv')

# %%


# ======== RUN MODEL ========== #
for i in range(variation_count):
    o3frac = get_ozone_fac(i, variation_count)
    print(o3frac)
    test_gsto_output_full_year(o3frac)


# %%
# ======== GRAPH DATA ========== #
loaded_data = pd.read_csv(SAVE_OUTPUTS_LOCATION + f'gsto_output_full_year_{0}.csv')

data_all = {f: [] for f in loaded_data.columns}
for i in range(variation_count):
    o3frac = get_ozone_fac(i, variation_count)
    loaded_data = pd.read_csv(SAVE_OUTPUTS_LOCATION + f'gsto_output_full_year_{o3frac}.csv')
    # print(loaded_data.head())
    for f in loaded_data.columns:
        data_all[f].append(loaded_data[f].values)

fields = ['Tleaf_C', 'g_sv', 'f_LS', 'A_n', 'A_c', 'A_j', 'A_p', '',
          'R_d', 'fO3_h_out', 'fO3_d_out', 'fO3_l', 'c_i', 'Canopy_A_n']

start = 150
end = 220
for f in loaded_data.columns:
    if f not in fields:
        continue
    series_titles = ['{:.0f}%'.format(100 * get_ozone_fac(i, variation_count))
                     for i in range(variation_count)]

    multi_series_annual_graph(
        f,
        [x[start * 24:end * 24] for x in data_all[f]],
        series_titles,
        Field(f, float),
        output_dir=SAVE_OUTPUTS_LOCATION,
        label_x_days=30,
        average_step=24,
        chart_id=f'{f}-comparison',
        linestyles=[{
            'style': 'dashed' if s == 'DO3SE_ui' else 'solid',
            'width': 0.8,
            'zorder': 100 if s == 'DO3SE_ui' else i,
        } for i, s in enumerate(series_titles)],
        figsize=(8, 3),
        dpi=300,
        output_format='svg',
        start_day=start,
        end_day=end,
    )

# %%
print([x[0 * 24:1 * 24] for x in data_all['f_LS']])

# %%
# ======== 3D GRAPH DATA ========== #
loaded_data = pd.read_csv(SAVE_OUTPUTS_LOCATION + f'gsto_output_full_year_{1.0}.csv')
day_count = 365

# plots
x = np.array([[k for k in range(24)] for j in [i for i in range(day_count)]]).flatten()
y = np.array([[j for k in range(24)] for j in [i for i in range(day_count)]]).flatten()

fig = plt.figure()
fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection='3d')

# z = loaded_data.g_sv
# ax.scatter(x, y, z, c=z)
# ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'g_sv.png')

ax.clear()
g_sto_full = np.array([DRATIO * (max(0.0, i) / 1000) for i in loaded_data['g_sv']])
z = g_sto_full
ax.scatter(x, y, z, c=z)
ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'g_sto.png')

# %%

fields = ['Tleaf_C', 'g_sv', 'f_LS', 'A_n', 'A_c', 'A_j', 'A_p', '',
          'R_d', 'fO3_h_out', 'fO3_d_out', 'fO3_l', 'c_i', 'Canopy_A_n']
for f in fields:
    try:
        ax.clear()
        z = loaded_data[f]
        ax.scatter(x, y, z, c=z)
        ax.figure.savefig(SAVE_OUTPUTS_LOCATION + f + '.png')
    except:
        pass

# %%

# A_C is purple, A_j is blue, A_p is yellow
ax.clear()
z = loaded_data['A_n']
c = ([0 if i == 'A_c' else 1 if i == 'A_j' else 2 if i ==  # noqaW504
      'A_p' else 5 for i in loaded_data['A_n_limit_factor']])
ax.scatter(x, y, z, c=c)
ax.figure.savefig(SAVE_OUTPUTS_LOCATION + 'A_n_limit.png')

# %%
from pyDO3SE.plugins.gsto.ewert import ewert
from math import isclose
from pyDO3SE.constants.physical_constants import DRATIO_O3_CO2
import numpy as np
import pandas as pd
import os
from dataclasses import asdict, is_dataclass
from math import isclose

import matplotlib
import matplotlib.pyplot as plt

from pyDO3SE.plugins.met.irradiance import (
    MLMC_sunlit_LAI, calc_Idrctt_Idfuse,
    calc_PAR_sun_shade, calc_PAR_sun_shade_farq_b,
    calc_is_daylight,
)
from pyDO3SE.plugins.met.solar_position import calc_solar_elevation

from pyDO3SE.plugins.met.thermal_time import get_season_astart_thermal_time, calc_thermal_time
from pyDO3SE.plugins.gsto.photosynthesis_helpers import calc_g_bv
from pyDO3SE.plugins.met.thermal_time import get_td_dd

%load_ext autoreload

# %%
%autoreload 2
t_lem_constant = 0.15
t_lse_constant = 0.33
season_Astart = 10

t_l = 800
t_lem = t_l * t_lem_constant
t_lma = t_l - t_lem
t_lse = t_lma * t_lse_constant
t_lep = t_l - (t_lem + t_lse)  # check this

assert t_l == t_lem + t_lep + t_lse
# %%
print(t_l, t_lem, t_lma, t_lse, t_lep)

day_count = 180
DEMO_DATA_LOCATION = './examples/spanish_wheat/spanish_wheat_data_test.csv'
SAVE_OUTPUTS_LOCATION = './pyDO3SE/plugins/gsto/ewert/output/'
df = pd.read_csv(DEMO_DATA_LOCATION)
data = {}

data["Ts_C"] = df['Ts_C'].values.reshape((365, 24))
data["u"] = df['u'].values.reshape((365, 24))
data["VPD"] = df['VPD'].values.reshape((365, 24))
data["PAR"] = df['PAR, W m-2'].values.reshape((365, 24))
data["P"] = df['P, kPa'].values.reshape((365, 24))
data["O3"] = df['O3, ppb'].values.reshape((365, 24))
data['fst'] = df['Fst (nmol/m^2/s)'].values.reshape((365, 24))
data['Lai'] = df['LAI'].values.reshape((365, 24))

td_data = np.array(calc_thermal_time(df['Ts_C'].values))
season_Astart_temp = get_season_astart_thermal_time(season_Astart, td_data, 0)
td_dd_full = [get_td_dd(dd, td_data, season_Astart, season_Astart_temp, dd * 24)
              for dd in range(day_count) for _ in range(24)]
t_l_day = next(x[0] for x in enumerate(td_dd_full) if x[1] > t_l)
t_lem_day = next(x[0] for x in enumerate(td_dd_full) if x[1] > t_lem)
t_lep_day = next(x[0] for x in enumerate(td_dd_full) if x[1] > t_lem + t_lep)

plt.plot(td_dd_full)
plt.axhline(y=t_l)
plt.axhline(y=t_lem)
plt.axhline(y=t_lem + t_lep)
plt.axvline(x=t_l_day)
plt.axvline(x=t_lem_day)
plt.axvline(x=t_lep_day)
# %%


def run_ewert(O3_override=None, O3_mult=None):
    """Test full gsto_pn method for a full year."""
    # previous_day_output = None
    output_data = []
    input_data = []
    D_0 = 2.27
    fO3_d_prev = 1
    O3up_prev = 0.0
    O3up_acc = 0.0
    td_dd_prev = 0.0
    # config
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
        fo3_d_acc = 0
        for hr in range(24):
            # O3_nmol = O3_ppb_to_nmol(data["Ts_C"][dd][hr], data["P"][dd][hr], data["O3"][dd][hr])
            O3up = O3_override if O3_override else data['fst'][dd][hr] * \
                O3_mult if O3_mult else data['fst'][dd][hr]
            td_dd = get_td_dd(dd, td_data, season_Astart, season_Astart_temp, dd * 24)
            g_bv = calc_g_bv(0.02, data["u"][dd][hr])
            LAI = 3
            # LAI = data['Lai'][dd][hr]
            sinB = calc_solar_elevation(40.43, -3.7, dd, hr)
            LAIsunfrac = MLMC_sunlit_LAI(1, 1, [[LAI]], sinB)[0][0]
            PAR = data["PAR"][dd][hr]
            P = data["P"][dd][hr]
            Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
            cosA = 0.5
            PARsun, PARshade = calc_PAR_sun_shade_farq_b(Idrctt, Idfuse, sinB, cosA, 1)
            # PARsun = data["PAR"][dd][hr]
            # PARshade = data["PAR"][dd][hr]
            is_daylight = calc_is_daylight(PARsun)

            ewert_inputs = {
                **ewert_constant_inputs,
                "PARsun": PARsun,
                "PARshade": PARshade,
                "LAI": LAI,
                "LAIsunfrac": LAIsunfrac,
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
                "use_O3_damage": True,
            }
            hr_output = ewert.ewert(**ewert_inputs)

            O3up_acc = O3up_acc + ((O3up + O3up_prev) / 2 * (td_dd - td_dd_prev))

            output_data.append(hr_output)
            input_data.append(ewert_inputs)

            fO3_d_prev = hr_output.fO3_d_out
            td_dd_prev = td_dd
            O3up_prev = O3up

        fO3_d_prev = fo3_d_acc

    return output_data, input_data


def get_attr_list(data, attr):
    return np.array([getattr(i, attr, None) if is_dataclass(i) else i.get(attr, None) for i in data])


output_data, input_data = run_ewert(0.5)
plt.plot(get_attr_list(output_data[24 * 0:24 * 180], 'f_LS'))
plt.axvline(x=t_l_day, color="red")
plt.axvline(x=t_lem_day, color="green")
plt.axvline(x=t_lep_day, color="blue")
# %%
fig, axs = plt.subplots(6, figsize=(8, 6), dpi=150)
axs[0].axvline(x=t_l_day, color="red")
axs[0].axvline(x=t_lem_day, color="red")
axs[0].axvline(x=t_lep_day, color="red")
for O3 in [0, 0.01, 0.1, 1]:
    output_data, input_data = run_ewert(O3)

    axs[0].plot(get_attr_list(output_data[24 * 0:24 * 180], 'f_LS'), label=O3)
    axs[0].set_title('f_LS')
    axs[1].plot(get_attr_list(output_data[24 * 0:24 * 180], 'f_LA'), label=O3)
    axs[1].set_title('f_LA')
    axs[2].plot(get_attr_list(output_data[24 * 0:24 * 180], 'fO3_d_out'), label=O3)
    axs[2].set_title('fO3_d')
    axs[3].plot(get_attr_list(output_data[24 * 0:24 * 180], 'g_sv'), label=O3)
    axs[3].set_title('g_sv')
    axs[4].plot(get_attr_list(output_data[24 * 0:24 * 180], 't_l_O3'), label=O3)
    axs[4].set_title('t_l_O3')
    axs[5].plot(get_attr_list(output_data[24 * 0:24 * 180], 'fO3_l'), label=O3)
    axs[5].set_title('fO3_l')

    axs[0].legend()

# %%
# fig, axs = plt.subplots(6, figsize=(8, 6), dpi=150)
# axs[0].axvline(x=t_l_day, color="red")
# axs[0].axvline(x=t_lem_day, color="red")
# axs[0].axvline(x=t_lep_day, color="red")
# for O3 in [None]:
#     # for O3 in [None, 0.00000001]:
#     output_data, input_data = run_ewert(O3)

#     axs[0].plot(get_attr_list(output_data[24 * 0:24 * 180], 'f_LS'), label=O3 or 'data')
#     axs[0].set_title('f_LS')
#     axs[1].plot(get_attr_list(output_data[24 * 0:24 * 180], 'f_LA'), label=O3 or 'data')
#     axs[1].set_title('f_LA')
#     axs[2].plot(get_attr_list(output_data[24 * 0:24 * 180], 'fO3_d_out'), label=O3 or 'data')
#     axs[2].set_title('fO3_d')
#     axs[3].plot(get_attr_list(output_data[24 * 0:24 * 180], 'g_sv'), label=O3 or 'data')
#     axs[3].set_title('g_sv')
#     axs[4].plot(get_attr_list(output_data[24 * 0:24 * 180], 't_l_O3'), label=O3 or 'data')
#     axs[4].set_title('t_l_O3')
#     axs[5].plot(get_attr_list(output_data[24 * 0:24 * 180], 'fO3_l'), label=O3 or 'data')
#     axs[5].set_title('fO3_l')
#     print(output_data[-1].t_l_O3)
#     final_t_l = next(x[0] for x in enumerate(td_dd_full) if x[1] > output_data[-1].t_l_O3)
#     axs[0].axvline(x=final_t_l, color="green")

#     axs[0].legend()

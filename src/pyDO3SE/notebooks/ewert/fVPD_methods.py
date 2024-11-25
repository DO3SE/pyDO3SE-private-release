# %%
from matplotlib import pyplot as plt
import os
from dataclasses import asdict
from math import isclose
from pyDO3SE.constants.enums import GAS

import pytest
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d  # noqa: F401
from tests.test_helpers import long_test


from pyDO3SE.constants.physical_constants import DRATIO_O3_CO2
from do3se_met.irradiance import (
    MLMC_sunlit_LAI, calc_Idrctt_Idfuse,
    calc_PAR_sun_shade,
    calc_is_daylight,
)
from do3se_met.solar_position import calc_solar_elevation

from pyDO3SE.plugins.gsto.photosynthesis_helpers import calc_g_bv
from thermal_time.calcs import get_td_dd


from pyDO3SE.plugins.gsto.ewert.ewert import (
    CO2_Constant_Loop_Inputs, CO2_loop_State, ModelOptions,
    co2_concentration_in_stomata_iteration,
    co2_concentration_in_stomata_loop,
    ewert,
)
f_VPD_method = [
    # "linear",
    "photosynthesis",
][0]


constant_inputs = CO2_Constant_Loop_Inputs(
    c_a=391.0,
    e_a=759.7060764490742,
    g_bl=1676057.8748957326,
    g_sto_0=10000.0,
    m=8.12,
    D_0=0.75,
    O3up=0.8705724925761334,
    O3up_acc=654.2909552766799,
    fO3_d_prev=0.9970964447727888,
    td_dd=300.5987499999997,
    gamma_1=0.06,
    gamma_2=0.0045,
    gamma_3=0.5,
    is_daylight=True,
    t_lse_constant=0.33,
    t_l_estimate=800.11214305013493,
    t_lem=155.5852561302345,
    t_lep=202.69301423633325,
    t_lse=99.83387268356714,
    t_lma=302.5268869199004,
    Gamma=53.5135179983076,
    Gamma_star=50.82724005285872,
    V_cmax=220.2995017637242,
    K_C=582.3161124585115,
    K_O=328.8128555882607,
    J=368.72268771159526,
    R_d=0.5872449303570749,
    e_sat_i=3874.6260764490744,
    hr=16,
    f_SW=1.0,
    f_VPD=0.2,
)


tolerance = 0.001
max_iterations = 100
c_i_start = 0
initial_c_i_diff = 1
model_options = ModelOptions(
    opt_full_night_recovery=False,
    f_VPD_method=f_VPD_method,
    use_O3_damage=True,
)

state = CO2_loop_State(
    c_i=c_i_start,
    c_i_diff=initial_c_i_diff,
    g_sto=constant_inputs.g_sto_0,
)

out_vals = []
while state.c_i_diff > tolerance and state.iterations < max_iterations:
    # print(state.iterations)
    state = co2_concentration_in_stomata_iteration(constant_inputs, state, model_options)
    out_vals.append(state)
fig, axs = plt.subplots(2, 1)
axs[0].plot([s.c_i for s in out_vals])
axs[1].plot([s.g_sto for s in out_vals])
state
# print(state_out)
# assert isinstance(state_out, CO2_loop_State)

# assert isclose(state_out.c_i, 303.181273378, abs_tol=1e-3)
# assert isclose(state_out.c_i_diff, 0.0005581661, abs_tol=1e-5)
# assert isclose(state_out.g_sto, 662369.31742, abs_tol=1e-3)
# assert isclose(state_out.f_LS, 1, abs_tol=1e-3)
# assert isclose(state_out.A_n, 35.966032219, abs_tol=1e-3)
# assert isclose(state_out.A_c, 36.28603221, abs_tol=1e-3)
# assert isclose(state_out.A_p, 59.5, abs_tol=1e-3)
# assert isclose(state_out.A_j, 54.9788411629, abs_tol=1e-3)
# assert isclose(state_out.fO3_d, 0.8628995, abs_tol=1e-3)
# assert isclose(state_out.fO3_h, 0.96955, abs_tol=1e-3)
# assert isclose(state_out.iterations, 19, abs_tol=1e-3)
# assert state_out.A_n_limit_factor == 'A_c'

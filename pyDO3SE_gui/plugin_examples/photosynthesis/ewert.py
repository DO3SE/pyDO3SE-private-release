"""Ewert model streamlit app."""

from dataclasses import asdict
import json
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d  # noqa: F401


from pyDO3SE.constants.physical_constants import DRATIO
from pyDO3SE.plugins.met.irradiance import (
    calc_Idrctt_Idfuse,
    calc_PAR_sun_shade,
    calc_is_daylight,
)
from pyDO3SE.plugins.met.solar_position import calc_solar_elevation

from pyDO3SE.plugins.met.thermal_time import get_season_astart_temperature, calc_thermal_time
from pyDO3SE.plugins.gsto.photosynthesis_helpers import calc_g_bv
from pyDO3SE.plugins.met.thermal_time import get_td_dd


from pyDO3SE.plugins.gsto.ewert.ewert import (
    CO2_Constant_Loop_Inputs, CO2_loop_State,
    co2_concentration_in_stomata_iteration,
    co2_concentration_in_stomata_loop,
    ewert,
)

t_lem_constant = 0.15
t_lse_constant = 0.33


@st.cache
def co2_concentration_in_stomata_iteration_fn(*args, **kwargs):
    co2_concentration_in_stomata_iteration(*args, **kwargs)


@st.cache
def get_demo_data(DEMO_DATA_LOCATION):
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
    return data, td_data


@st.cache
def ewert_cached(*args, **kwargs):
    return ewert(*args, **kwargs)


@st.cache
def ewert_full_year_cached(*args, **kwargs):
    pass


def gui_ewert_full_year():
    DEMO_DATA_LOCATION = './examples/spanish_wheat/data/spanish_wheat_data_test.csv'

    day_count = st.slider("Day Count", 0, 365, 10, 1)
    t_l = 800.0  # st.sidebar.slider("t_l_estimate", 0.0, 10000.0, 800.0, 0.1)
    t_lem = t_l * t_lem_constant
    t_lma = t_l - t_lem
    t_lse = t_l - (t_lma * t_lse_constant)
    t_lep = t_l - (t_lem + t_lse)  # check this

    # import data
    data, td_data = get_demo_data(DEMO_DATA_LOCATION)

    # previous_day_output = None
    output_data = []
    input_data = []
    D_0 = 2.27
    fO3_d_prev = 1
    O3up_acc = 0
    # config
    season_Astart = 180
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

    progress_bar = st.progress(0)
    for dd in range(day_count):
        fo3_d_acc = 0
        td = td_data[dd]

        for hr in range(24):
            progress_bar.progress(((dd * 24) + hr) / (day_count * 24))
            O3up = data['fst'][dd][hr]

            td_dd = get_td_dd(dd, td, season_Astart, season_Astart_temp)
            g_bv = calc_g_bv(0.02, data["u"][dd][hr])
            LAI = 0.01
            sinB = calc_solar_elevation(40.43, -3.7, dd, hr)
            PAR = data["PAR"][dd][hr]
            P = data["P"][dd][hr]
            Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
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
            hr_output = ewert_cached(**ewert_inputs)

            fO3_d_prev = hr_output.fO3_d_out
            O3up_acc = O3up_acc + O3up
            output_data.append(hr_output)
            input_data.append(ewert_inputs)
        fO3_d_prev = fo3_d_acc

    st.write('Complete')
    df_out = pd.DataFrame([asdict(i) for i in output_data])
    df_out['g_sto'] = [DRATIO * (max(0.0, i.g_sv) / 1000) for i in output_data]

    st.write('Input Data')
    st.dataframe(input_data)
    st.write('Output Data')
    st.dataframe([asdict(i) for i in output_data])
    # st.line_chart(df_out['A_n'])
    # st.line_chart(df_out['g_sto'])
    st.line_chart(df_out['g_sto'])

    # plots
    x = np.array([[k for k in range(24)] for j in [i for i in range(day_count)]]).flatten()
    y = np.array([[j for k in range(24)] for j in [i for i in range(day_count)]]).flatten()

    selection = st.selectbox("Choose variable to graph", df_out.columns)

    fig = plt.figure()
    fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111, projection='3d')

    z = df_out[selection].values
    ax.scatter(x, y, z, c=z)
    st.pyplot(fig)

    # fig = plt.figure()
    # fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
    # ax = fig.add_subplot(111, projection='3d')
    # z = df_out['A_n'].values
    # c = [0 if i == 'A_c' else 1 if i == 'A_j' else 2 if i ==
    #      'A_p' else 5 for i in df_out['A_n_limit_factor'].values]
    # ax.scatter(x, y, z, c=c)
    # st.pyplot(fig)


def gui_ewert():
    st.title("pyDO3SE - Photosynthesis - ewert model")

    t_l = 800.0  # st.sidebar.slider("t_l_estimate", 0.0, 10000.0, 800.0, 0.1)
    Tleaf_C = st.sidebar.slider("Tleaf_C", 0.0, 100.0, 20.0, 0.1)
    PARsun = st.sidebar.slider("PARsun", 0.0, 3000.0, 800.0, 0.1)
    PARshade = st.sidebar.slider("PARshade", 0.0, 3000.0, 400.0, 0.1)
    V_cmax_25 = st.sidebar.slider("V_cmax_25", 0.0, 1000.0, 180.0, 0.1)
    J_max_25 = st.sidebar.slider("J_max_25", 0.0, 1000.0, 400.0, 0.1)
    LAI = st.sidebar.slider("LAI", 0.0, 1.0, 0.01, 0.01)
    sinB = st.sidebar.slider("sinB", 0.0, 2.0, 0.3, 0.01)
    c_a = st.sidebar.slider("c_a", 0.0, 10000.0, 391.0, 0.1)
    eact = st.sidebar.slider("eact", 0.0, 2.0, 1.0, 0.1)
    g_bv = st.sidebar.slider("g_bv", 0.0, 3000000.0, 1469999.0, 0.1)
    g_sto_0 = st.sidebar.slider("g_sto_0", 0.0, 30000.0, 20000.0, 0.1)
    m = st.sidebar.slider("m", 0.0, 100.0, 8.12, 0.01)
    D_0 = st.sidebar.slider("D_0", 0.0, 10.0, 2.27, 0.01)
    O3up = st.sidebar.slider("O3up", 0.0, 100.0, 80.0, 0.1)
    O3up_acc = st.sidebar.slider("O3up_acc", 0.0, 100.0, 30000.0, 0.1)
    fO3_d_prev = st.sidebar.slider("fO3_d_prev", 0.0, 1.0, 0.89, 0.01)
    td_dd = st.sidebar.slider("td_dd", 0.0, 1000.0, 24.1, 0.01)
    gamma_1 = st.sidebar.slider("gamma_1", 0.0, 1.0, 0.06, 0.001)
    gamma_2 = st.sidebar.slider("gamma_2", 0.0, 1.0, 0.0045, 0.0001)
    gamma_3 = st.sidebar.slider("gamma_3", 0.0, 1.0, 0.5, 0.001)
    is_daylight = st.sidebar.checkbox("is_daylight", True)

    t_lem = t_l * t_lem_constant
    t_lma = t_l - t_lem
    t_lse = t_l - (t_lma * t_lse_constant)
    t_lep = t_l - (t_lem + t_lse)  # check this
    out = ewert(
        t_lse_constant=t_lse_constant,
        t_l_estimate=t_l,
        t_lem=t_lem,
        t_lep=t_lep,
        t_lse=t_lse,
        t_lma=t_lma,

        Tleaf_C=Tleaf_C,
        PARsun=PARsun,
        PARshade=PARshade,
        V_cmax_25=V_cmax_25,
        J_max_25=J_max_25,
        LAI=LAI,
        sinB=sinB,
        c_a=c_a,
        eact=eact,
        g_bv=g_bv,
        g_sto_0=g_sto_0,
        m=m,
        O3up=O3up,
        O3up_acc=O3up_acc,
        D_0=D_0,
        fO3_d_prev=fO3_d_prev,
        td_dd=td_dd,
        gamma_1=gamma_1,
        gamma_2=gamma_2,
        gamma_3=gamma_3,
        is_daylight=is_daylight,
    )
    st.json(json.dumps(asdict(out)))


def gui_ewert_co2_concentration_in_stomata_iteration():
    st.title("pyDO3SE - Photosynthesis - ewert - gui_ewert_co2_concentration_in_stomata_iteration")

    c_a = st.sidebar.slider("c_a", 0.0, 800.0, 391.0, 0.1)
    e_a = st.sidebar.slider("e_a", 0.0, 10000.0, 1000.0, 0.1)
    g_bl = 1469999.0  # st.sidebar.slider("g_bl", 0.0, 2000000.0, 1469999.0, 0.1)
    g_sto_0 = 20000.0  # st.sidebar.slider("g_sto_0", 0.0, 40000.0, 20000.0, 0.1)
    m = 8.12  # st.sidebar.slider("m", 0.0, 10.0, 8.12, 0.1)
    D_0 = 2.27 * 1e3  # st.sidebar.slider("D_0", 0.0, 10000.0, 2.27 * 1e3, 0.1)
    O3up = 80.0  # st.sidebar.slider("O3up_prev", 0.0, 1000000.0, 80.0, 0.1)
    O3up_acc = 300.0  # st.sidebar.slider("O3up_acc_prev", 0.0, 10000000.0, 300.0, 0.1)
    fO3_d_prev = 0.89  # st.sidebar.slider("fO3_d_prev", 0.0, 1.0, 0.89, 0.1)
    td_dd = 24.1  # st.sidebar.slider("td_dd", 0.0, 2000.0, 24.1, 0.1)
    gamma_1 = 0.06  # st.sidebar.slider("gamma_1", 0.0, 1.0, 0.06, 0.1)
    gamma_2 = 0.0045  # st.sidebar.slider("gamma_2", 0.0, 1.0, 0.0045, 0.1)
    gamma_3 = 0.5  # st.sidebar.slider("gamma_3", 0.0, 1.0, 0.5, 0.1)
    is_daylight = st.sidebar.checkbox("is_daylight", True)
    t_l = 800.0  # st.sidebar.slider("t_l_estimate", 0.0, 10000.0, 800.0, 0.1)
    Gamma = 34.277  # st.sidebar.slider("Gamma", 0.0, 100.0, 34.277, 0.1)
    Gamma_star = 32.95  # st.sidebar.slider("Gamma_star", 0.0, 100.0, 32.95, 0.1)
    V_cmax = 119.0  # st.sidebar.slider("V_cmax", 0.0, 1000.0, 119.0, 0.1)
    K_C = 234.42  # st.sidebar.slider("K_C", 0.0, 1000.0, 234.42, 0.1)
    K_O = 216.75  # st.sidebar.slider("K_O", 0.0, 1000.0, 216.75, 0.1)
    J = 300.36  # st.sidebar.slider("J", 0.0, 1000.0, 300.36, 0.1)
    R_d = 0.32  # st.sidebar.slider("R_d", 0.0, 1.0, 0.32, 0.1)
    e_sat_i = 2339.05  # st.sidebar.slider("e_sat_i", 0.0, 10000.0, 2339.05, 0.1)

    t_lem = t_l * t_lem_constant
    t_lma = t_l - t_lem
    t_lse = t_l - (t_lma * t_lse_constant)
    t_lep = t_l - (t_lem + t_lse)  # check this

    st.text('Loading')

    state_out = co2_concentration_in_stomata_iteration(
        CO2_Constant_Loop_Inputs(
            c_a=c_a,
            e_a=e_a,
            g_bl=g_bl,
            g_sto_0=g_sto_0,
            m=m,
            D_0=D_0,
            O3up=O3up,
            O3up_acc=O3up_acc,
            fO3_d_prev=fO3_d_prev,
            td_dd=td_dd,
            gamma_1=gamma_1,
            gamma_2=gamma_2,
            gamma_3=gamma_3,
            is_daylight=is_daylight,
            t_lse_constant=t_lse_constant,
            t_l_estimate=t_l,
            t_lem=t_lem,
            t_lep=t_lep,
            t_lse=t_lse,
            t_lma=t_lma,
            Gamma=Gamma,
            Gamma_star=Gamma_star,
            V_cmax=V_cmax,
            K_C=K_C,
            K_O=K_O,
            J=J,
            R_d=R_d,
            e_sat_i=e_sat_i,
        ),
        CO2_loop_State(
            c_i=0.0,
            c_i_diff=0,
            g_sto=20000,
            A_n=0,
            fO3_d=1.0,
            fO3_h=1.0,
        )
    )

    st.json(json.dumps(asdict(state_out)))


def gui_ewert_co2_concentration_in_stomata_loop():
    st.title("pyDO3SE - Photosynthesis - ewert - gui_ewert_co2_concentration_in_stomata_loop")

    c_a = st.sidebar.slider("c_a", 0.0, 800.0, 391.0, 0.1)
    e_a = st.sidebar.slider("e_a", 0.0, 10000.0, 1000.0, 0.1)
    g_bl = st.sidebar.slider("g_bl", 0.0, 2000000.0, 1469999.0, 0.1)
    g_sto_0 = st.sidebar.slider("g_sto_0", 0.0, 40000.0, 20000.0, 0.1)
    m = 8.12  # st.sidebar.slider("m", 0.0, 10.0, 8.12, 0.1)
    D_0 = 2.27 * 1e3  # st.sidebar.slider("D_0", 0.0, 10000.0, 2.27 * 1e3, 0.1)
    O3up = 80.0  # st.sidebar.slider("O3up", 0.0, 1000000.0, 80.0, 0.1)
    O3up_acc = 300.0  # st.sidebar.slider("O3up_acc", 0.0, 10000000.0, 300.0, 0.1)
    fO3_d_prev = 0.89  # st.sidebar.slider("fO3_d_prev", 0.0, 1.0, 0.89, 0.1)
    td_dd = st.sidebar.slider("td_dd", 0.0, 2000.0, 24.1, 0.1)
    gamma_1 = 0.06  # st.sidebar.slider("gamma_1", 0.0, 1.0, 0.06, 0.1)
    gamma_2 = 0.0045  # st.sidebar.slider("gamma_2", 0.0, 1.0, 0.0045, 0.1)
    gamma_3 = 0.5  # st.sidebar.slider("gamma_3", 0.0, 1.0, 0.5, 0.1)
    is_daylight = st.sidebar.checkbox("is_daylight", True)
    t_l = 800.0  # st.sidebar.slider("t_l_estimate", 0.0, 10000.0, 800.0, 0.1)
    Gamma = 34.277  # st.sidebar.slider("Gamma", 0.0, 100.0, 34.277, 0.1)
    Gamma_star = 32.95  # st.sidebar.slider("Gamma_star", 0.0, 100.0, 32.95, 0.1)
    V_cmax = 119.0  # st.sidebar.slider("V_cmax", 0.0, 1000.0, 119.0, 0.1)
    K_C = 234.42  # st.sidebar.slider("K_C", 0.0, 1000.0, 234.42, 0.1)
    K_O = 216.75  # st.sidebar.slider("K_O", 0.0, 1000.0, 216.75, 0.1)
    J = 300.36  # st.sidebar.slider("J", 0.0, 1000.0, 300.36, 0.1)
    R_d = 0.32  # st.sidebar.slider("R_d", 0.0, 1.0, 0.32, 0.1)
    e_sat_i = 2339.05  # st.sidebar.slider("e_sat_i", 0.0, 10000.0, 2339.05, 0.1)

    t_lem = t_l * t_lem_constant
    t_lma = t_l - t_lem
    t_lse = t_l - (t_lma * t_lse_constant)
    t_lep = t_l - (t_lem + t_lse)  # check this

    st.text('Loading')

    state_out = co2_concentration_in_stomata_loop(
        CO2_Constant_Loop_Inputs(
            c_a=c_a,
            e_a=e_a,
            g_bl=g_bl,
            g_sto_0=g_sto_0,
            m=m,
            D_0=D_0,
            O3up=O3up,
            O3up_acc=O3up_acc,
            fO3_d_prev=fO3_d_prev,
            td_dd=td_dd,
            gamma_1=gamma_1,
            gamma_2=gamma_2,
            gamma_3=gamma_3,
            is_daylight=is_daylight,
            t_lse_constant=t_lse_constant,
            t_l_estimate=t_l,
            t_lem=t_lem,
            t_lep=t_lep,
            t_lse=t_lse,
            t_lma=t_lma,
            Gamma=Gamma,
            Gamma_star=Gamma_star,
            V_cmax=V_cmax,
            K_C=K_C,
            K_O=K_O,
            J=J,
            R_d=R_d,
            e_sat_i=e_sat_i,
        )
    )

    st.json(json.dumps(asdict(state_out)))

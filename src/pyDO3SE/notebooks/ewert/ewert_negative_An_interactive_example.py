# %%
from ipywidgets import interact
import ipywidgets as widgets
import numpy as np
from matplotlib import pyplot as plt

from pyDO3SE.plugins.gsto.ewert.ewert import (
    CO2_Constant_Loop_Inputs,
    CO2_loop_State,
    ModelOptions,
    co2_concentration_in_stomata_iteration,
    FVPDMethods,
    ewert_leaf_pop,
    ewert_leaf_pop_alt,
)


# %%
# Interactive example of a failed convergence
@interact(
    max_iterations=widgets.IntSlider(min=1, max=100, step=1, value=50),
    threshold=widgets.FloatSlider(min=0.001, max=1, step=0.001, value=0.01),
    g_bl=widgets.FloatSlider(
        min=1000000, max=10000000, step=100000, value=4594662.190354731
    ),
    c_a=widgets.FloatSlider(min=300, max=500, step=1, value=391.0),
    e_a=widgets.FloatSlider(min=600, max=800, step=1, value=712.529686719428),
    g_sto_0=widgets.FloatSlider(min=5000, max=20000, step=100, value=10000.0),
    m=widgets.FloatSlider(min=1, max=10, step=0.1, value=4),
    D_0=widgets.FloatSlider(min=1, max=5, step=0.1, value=2.2),
    Gamma=widgets.FloatSlider(min=50, max=70, step=0.1, value=57.89952113732494),
    Gamma_star=widgets.FloatSlider(min=40, max=60, step=0.1, value=48.51315602386692),
    V_cmax=widgets.FloatSlider(min=10, max=200, step=1, value=159.98275062496984),
    K_C=widgets.FloatSlider(min=500, max=600, step=1, value=528.0415080745356),
    K_O=widgets.FloatSlider(min=300, max=400, step=1, value=314.4035643006643),
    J=widgets.FloatSlider(min=100, max=200, step=1, value=180.94373072965217),
    R_d=widgets.FloatSlider(min=1, max=2, step=0.01, value=1.5998275062496985),
    e_sat_i=widgets.FloatSlider(min=3000, max=4000, step=10, value=3670.6412787571467),
    f_SW=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
    f_LS=widgets.FloatSlider(min=0, max=1, step=0.1, value=1),
    fO3_d=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
    f_VPD=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
)
def run_loop_failed(
    max_iterations=50,
    threshold=0.01,
    g_bl=4594662.190354731,
    c_a=391.0,
    e_a=712.529686719428,
    g_sto_0=10000.0,
    m=4,
    D_0=2.2,
    Gamma=57.89952113732494,
    Gamma_star=48.51315602386692,
    V_cmax=159.98275062496984,
    K_C=528.0415080745356,
    K_O=314.4035643006643,
    J=180.94373072965217,
    R_d=1.5998275062496985,
    e_sat_i=3670.6412787571467,
    f_SW=1.0,
    f_LS=1,
    fO3_d=1.0,
    f_VPD=1.0,
):
    constant_inputs = CO2_Constant_Loop_Inputs(
        c_a=c_a,
        e_a=e_a,
        g_bl=g_bl,
        g_sto_0=g_sto_0,
        m=m,
        D_0=D_0,
        Gamma=Gamma,
        Gamma_star=Gamma_star,
        V_cmax=V_cmax,
        K_C=K_C,
        K_O=K_O,
        J=J,
        R_d=R_d,
        e_sat_i=e_sat_i,
        f_SW=f_SW,
        f_LS=f_LS,
        fO3_d=fO3_d,
        f_VPD=f_VPD,
    )
    state = CO2_loop_State(
        c_i=0.0,
        c_i_diff=1,
        g_sto=10000.0,
        A_n=0.0,
        A_c=0.0,
        A_p=0.0,
        A_j=0.0,
        A_n_limit_factor=None,
        f_VPD=None,
        iterations=0,
    )
    model_options = ModelOptions(
        f_VPD_method=FVPDMethods.DANIELSSON,
        co2_concentration_balance_threshold=threshold,
        co2_concentration_max_iterations=max_iterations,
    )
    iteration_state_s = []
    while (
        state.c_i_diff > model_options.co2_concentration_balance_threshold
        and state.iterations < model_options.co2_concentration_max_iterations
    ):
        state = co2_concentration_in_stomata_iteration(
            constant_inputs, state, model_options
        )
        iteration_state_s.append(state)
    len(iteration_state_s)
    plt.close()
    fig, axs = plt.subplots(1, 5, figsize=(20, 5))
    x = np.arange(len(iteration_state_s))
    c = [cc % 2 for cc in range(len(iteration_state_s))]
    axs[0].scatter(x, [state.c_i_diff for state in iteration_state_s], c=c)
    axs[0].plot(x, [state.c_i_diff for state in iteration_state_s])
    axs[0].scatter(
        np.arange(len(iteration_state_s)),
        [state.c_i_diff for state in iteration_state_s],
        c=c,
    )
    axs[0].set_title("CO2 concentration diff")
    axs[0].set_xlabel("Iteration")
    axs[0].set_ylabel("CO2 concentration difference")
    axs[0].hlines(
        model_options.co2_concentration_balance_threshold,
        0,
        len(iteration_state_s),
        colors="r",
        linestyles="dashed",
    )
    axs[0].text(0.5, 100, iteration_state_s[-1].A_n_limit_factor)

    axs[1].scatter(x, [state.A_n for state in iteration_state_s], c=c)
    axs[1].set_title("A_n")
    axs[2].scatter(x, [state.g_sto for state in iteration_state_s], c=c)
    axs[2].set_title("gsto")
    axs[3].scatter(x, [state.c_i for state in iteration_state_s], c=c)
    axs[3].set_title("CO2")
    if len(iteration_state_s) == max_iterations:
        axs[0].text(0.5, 50, "Did not converge", c="red")

    axs[4].scatter(x, [state.f_VPD for state in iteration_state_s])
    axs[4].set_title("f_VPD")


# %%
# Interactive example of a succesful convergence.
@interact(
    max_iterations=widgets.IntSlider(min=1, max=100, step=1, value=50),
    threshold=widgets.FloatSlider(min=0.001, max=1, step=0.001, value=0.01),
    c_a=widgets.FloatSlider(min=300, max=500, step=1, value=391.0),
    e_a=widgets.FloatSlider(min=600, max=800, step=1, value=1013.8537256970482),
    g_bl=widgets.FloatSlider(
        min=10, max=10000000, step=100000, value=2684466.893155262
    ),
    g_sto_0=widgets.FloatSlider(min=5000, max=20000, step=100, value=10000.0),
    m=widgets.FloatSlider(min=1, max=10, step=0.1, value=4),
    D_0=widgets.FloatSlider(min=1, max=5, step=0.1, value=2.2),
    Gamma=widgets.FloatSlider(min=0, max=70, step=0.1, value=22.573130096737337),
    Gamma_star=widgets.FloatSlider(min=0, max=60, step=0.1, value=20.206142428451084),
    V_cmax=widgets.FloatSlider(min=0, max=200, step=1, value=37.03248350125308),
    K_C=widgets.FloatSlider(min=0, max=600, step=1, value=83.94772418584698),
    K_O=widgets.FloatSlider(min=0, max=400, step=1, value=135.42252855524103),
    J=widgets.FloatSlider(min=0, max=300, step=1, value=36.197064918673746),
    R_d=widgets.FloatSlider(min=0, max=2, step=0.01, value=0.37032483501253083),
    e_sat_i=widgets.FloatSlider(min=0, max=4000, step=10, value=1317.014541860481),
    f_SW=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
    f_LS=widgets.FloatSlider(min=0, max=1, step=0.1, value=1),
    fO3_d=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
    f_VPD=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
)
def run_loop_converges(
    max_iterations=50,
    threshold=0.01,
    c_a=391.0,
    e_a=1013.8537256970482,
    g_bl=2684466.893155262,
    g_sto_0=10000.0,
    m=4,
    D_0=2.2,
    Gamma=22.573130096737337,
    Gamma_star=20.206142428451084,
    V_cmax=37.03248350125308,
    K_C=83.94772418584698,
    K_O=135.42252855524103,
    J=36.197064918673746,
    R_d=0.37032483501253083,
    e_sat_i=1317.014541860481,
    f_SW=1.0,
    f_LS=1,
    fO3_d=1.0,
    f_VPD=1.0,
):
    constant_inputs = CO2_Constant_Loop_Inputs(
        c_a=c_a,
        e_a=e_a,
        g_bl=g_bl,
        g_sto_0=g_sto_0,
        m=m,
        D_0=D_0,
        Gamma=Gamma,
        Gamma_star=Gamma_star,
        V_cmax=V_cmax,
        K_C=K_C,
        K_O=K_O,
        J=J,
        R_d=R_d,
        e_sat_i=e_sat_i,
        f_SW=f_SW,
        f_LS=f_LS,
        fO3_d=fO3_d,
        f_VPD=f_VPD,
    )
    state = CO2_loop_State(
        c_i=0.0,
        c_i_diff=1,
        g_sto=10000.0,
        A_n=0.0,
        A_c=0.0,
        A_p=0.0,
        A_j=0.0,
        A_n_limit_factor=None,
        f_VPD=None,
        iterations=0,
    )
    model_options = ModelOptions(
        f_VPD_method=FVPDMethods.DANIELSSON,
        co2_concentration_balance_threshold=threshold,
        co2_concentration_max_iterations=max_iterations,
    )
    iteration_state_s = []
    while (
        state.c_i_diff > model_options.co2_concentration_balance_threshold
        and state.iterations < model_options.co2_concentration_max_iterations
    ):
        state = co2_concentration_in_stomata_iteration(
            constant_inputs, state, model_options
        )
        iteration_state_s.append(state)
    len(iteration_state_s)
    plt.close()
    fig, axs = plt.subplots(1, 5, figsize=(20, 5))
    x = np.arange(len(iteration_state_s))
    c = [cc % 2 for cc in range(len(iteration_state_s))]
    axs[0].scatter(x, [state.c_i_diff for state in iteration_state_s], c=c)
    axs[0].plot(x, [state.c_i_diff for state in iteration_state_s])
    axs[0].scatter(
        np.arange(len(iteration_state_s)),
        [state.c_i_diff for state in iteration_state_s],
        c=c,
    )
    axs[0].set_title("CO2 concentration diff")
    axs[0].set_xlabel("Iteration")
    axs[0].set_ylabel("CO2 concentration difference")
    axs[0].hlines(
        model_options.co2_concentration_balance_threshold,
        0,
        len(iteration_state_s),
        colors="r",
        linestyles="dashed",
    )
    axs[0].text(0.5, 100, iteration_state_s[-1].A_n_limit_factor)

    axs[1].scatter(x, [state.A_n for state in iteration_state_s], c=c)
    axs[1].set_title("A_n")
    axs[2].scatter(x, [state.g_sto for state in iteration_state_s], c=c)
    axs[2].set_title("gsto")
    axs[3].scatter(x, [state.c_i for state in iteration_state_s], c=c)
    axs[3].set_title("CO2")
    if len(iteration_state_s) == max_iterations:
        axs[0].text(0.5, 50, "Did not converge", c="red")

    axs[4].scatter(x, [state.f_VPD for state in iteration_state_s])
    axs[4].set_title("f_VPD")


# %%

from notebooks.ewert.calculate_A_n_via_cubic_equations import (
    calc_A_n
)

# Interactive example of a failed convergence
@interact(
    max_iterations=widgets.IntSlider(min=1, max=100, step=1, value=50),
    threshold=widgets.FloatSlider(min=0.001, max=1, step=0.001, value=0.01),
    g_bl=widgets.FloatSlider(
        min=1000000, max=10000000, step=100000, value=4594662.190354731
    ),
    c_a=widgets.FloatSlider(min=300, max=500, step=1, value=391.0),
    e_a=widgets.FloatSlider(min=600, max=800, step=1, value=712.529686719428),
    g_sto_0=widgets.FloatSlider(min=5000, max=20000, step=100, value=10000.0),
    m=widgets.FloatSlider(min=1, max=10, step=0.1, value=4),
    D_0=widgets.FloatSlider(min=1, max=5, step=0.1, value=2.2),
    Gamma=widgets.FloatSlider(min=50, max=70, step=0.1, value=57.89952113732494),
    Gamma_star=widgets.FloatSlider(min=40, max=60, step=0.1, value=48.51315602386692),
    V_cmax=widgets.FloatSlider(min=10, max=200, step=1, value=159.98275062496984),
    K_C=widgets.FloatSlider(min=500, max=600, step=1, value=528.0415080745356),
    K_O=widgets.FloatSlider(min=300, max=400, step=1, value=314.4035643006643),
    J=widgets.FloatSlider(min=100, max=200, step=1, value=180.94373072965217),
    R_d=widgets.FloatSlider(min=1, max=2, step=0.01, value=1.5998275062496985),
    e_sat_i=widgets.FloatSlider(min=3000, max=4000, step=10, value=3670.6412787571467),
    f_SW=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
    f_LS=widgets.FloatSlider(min=0, max=1, step=0.1, value=1),
    fO3_d=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
    f_VPD=widgets.FloatSlider(min=0, max=1, step=0.1, value=1.0),
)
def run_loop_failed_alt(
    max_iterations=50,
    threshold=0.01,
    g_bl=4594662.190354731,
    c_a=391.0,
    e_a=712.529686719428,
    g_sto_0=10000.0,
    m=4,
    D_0=2.2,
    Gamma=57.89952113732494,
    Gamma_star=48.51315602386692,
    V_cmax=159.98275062496984,
    K_C=528.0415080745356,
    K_O=314.4035643006643,
    J=180.94373072965217,
    R_d=1.5998275062496985,
    e_sat_i=3670.6412787571467,
    f_SW=1.0,
    f_LS=1,
    fO3_d=1.0,
    f_VPD=1.0,
):
    A_n = calc_A_n(
        O2=20900,
        P=101,
        K_O=K_O,
        V_cmax=V_cmax,
        R_d=R_d,
        G_0c=0.006625,
        G_1c=5.625,
        g_bv=0.21,
        c_a=c_a,
        RH=1.0,
        Gamma=Gamma,
        K_C=K_C,
        C_s=361,
        J=J,
        negative_A_n=False
    )
    print(A_n)

#     O2 - Partial pressure of atmospheric O2 [Pa] - which was set to 20900 in the model
# G_0c - which is set to 0.006625 as a fixed model parameter
# G_1c - which is set to 5.625 as a fixed model parameter
# C_s - CO2 concentration at leaf surface, which we could link to the multilayer concentration calculations
# g_bv - Boundary layer conductance, this is in the model but it seems the values are different and some work is needed to work out which value to use but I've currently set it to 0.21
# %%


from itertools import product
from matplotlib import pyplot as plt

from pyDO3SE.plugins.gsto.ewert.ewert import (
    CO2_Constant_Loop_Inputs,
    CO2_loop_State,
    ModelOptions,
    co2_concentration_in_stomata_iteration,
    FVPDMethods,
    ewert_leaf_pop,
    ewert_leaf_pop_alt,
)

# TODO: Get inputs from negative run
g_sto_0_list = [20000, 40000]
m_list = [8.12]
V_cmax_25_list = [40, 80, 120, 160, 180.0]
J_max_25_list = [400]
R_d_coeff_list = [0.015]
PARsun_list = [300, 400,500,600,700,800,900,1000]
PARshade_list = [300, 400,500,600,700,800,900,1000]
D_0_list = [2.27]
g_bv_list = [1469999.0]
Tleaf_C_list = [0.0,2.0,4.0,6.0,15.0,20.0,25.0,30.0,35.0,40.0]
eact_list = [1.0]
LAI = 0.01
nL = 3

results = []

for [
    g_sto_0,
    m,
    V_cmax_25,
    J_max_25,
    R_d_coeff,
    PARsun,
    PARshade,
    D_0,
    g_bv,
    Tleaf_C,
    eact,
] in product(
    g_sto_0_list,
    m_list,
    V_cmax_25_list,
    J_max_25_list,
    R_d_coeff_list,
    PARsun_list,
    PARshade_list,
    D_0_list,
    g_bv_list,
    Tleaf_C_list,
    eact_list,
):
    out = ewert_leaf_pop(
        nL=nL,
        g_sto_0=g_sto_0,
        m=8.12,
        V_cmax_25=[V_cmax_25 for _ in range(nL)],
        J_max_25=[J_max_25 for _ in range(nL)],
        R_d_coeff=R_d_coeff,
        PARsun=[PARsun for _ in range(nL)],
        PARshade=[PARshade for _ in range(nL)],
        LAIsunfrac=[0.9, 0.8, 0.7],
        D_0=D_0,
        g_bv=g_bv,
        # leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        layer_lai_frac=[0.33, 0.33, 0.33],
        layer_lai=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        Tleaf_C=[Tleaf_C for _ in range(nL)],
        f_SW=1,
        f_VPD=1.0,
        f_LS=0.89,
        fO3_d=0.89,
        eact=eact,
        c_a=391.0,
        f_VPD_method=FVPDMethods.LEUNING,
    )

    out_alt = ewert_leaf_pop_alt(
       nL=nL,
        g_sto_0=g_sto_0,
        m=8.12,
        V_cmax_25=[V_cmax_25 for _ in range(nL)],
        J_max_25=[J_max_25 for _ in range(nL)],
        R_d_coeff=R_d_coeff,
        PARsun=[PARsun for _ in range(nL)],
        PARshade=[PARshade for _ in range(nL)],
        LAIsunfrac=[0.9, 0.8, 0.7],
        D_0=D_0,
        g_bv=g_bv,
        # leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        layer_lai_frac=[0.33, 0.33, 0.33],
        layer_lai=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        Tleaf_C=[Tleaf_C for _ in range(nL)],
        f_SW=1,
        f_VPD=1.0,
        f_LS=0.89,
        fO3_d=0.89,
        eact=eact,
        c_a=391.0,
        f_VPD_method=FVPDMethods.LEUNING,
    )
    results.append([out.A_n, out_alt.A_n])

x = [r[0] for r in results]
y = [r[1] for r in results]

min_val = min(min(x), min(y))
max_val = max(max(x), max(y))
print(min_val, max_val)
print(len(results))
# plt.plot([min_val, min_val], [max_val, max_val])
plt.scatter(x, y)

# %%

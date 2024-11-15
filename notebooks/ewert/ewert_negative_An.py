# %%

import pandas as pd
from do3se_met.enums import GAS
from do3se_met.resistance.resistance import calc_g_bv
from pyDO3SE.plugins.gsto.ewert.ewert import (
    CO2_Constant_Loop_Inputs, CO2_loop_State, ModelOptions,
    co2_concentration_in_stomata_iteration,
    co2_concentration_in_stomata_loop,
    ewert_leaf_pop,
    calc_CO2_supply,
    FVPDMethods,
)
from ipywidgets import interact, interact_manual
import ipywidgets as widgets
from math import exp
import numpy as np
from matplotlib import pyplot as plt
# %%



def temp_dep(
    P_ref: float,
    T_ref: float,
    H_a: float,
    T: float,
    R: float,
) -> float:
    """Calculate parameter from temperature dependence curve.

    Parameters
    ----------
    P_ref (float):
        Parameter value at T_ref [unitless]
    T_ref (float):
        Reference temperature [degrees K]
    H_a (float):
        Activation energy [J/mol]
    T (float):
        Temperature [degrees K]
    R (float):
        Universal gas constant [J K-1 mol-1] (CONSTANT)
    Returns
    -------
    P (float):
        [unitless]
    """
    P = P_ref * exp((H_a * (T - T_ref)) / (T_ref * R * T))
    return P


def deg_to_kel(v: float) -> float:
    return 273.15 + v


R: float = 8.314472  # real, parameter
E_Gamma_star = 37830.0  # activation energy for C-comp-point [J/mol]            Medlyn2002
Gamma_star_25 = 42.75  # CO2 compensation point at T= 25    [micro mol/mol]    Medlyn2002
x_values = np.arange(-20, 50, 1)
Tleaf_K_values = deg_to_kel(x_values)
# Tleaf_K_values = [deg_to_kel(5), deg_to_kel(15), deg_to_kel(25), deg_to_kel(35), deg_to_kel(45)]
y = [
    temp_dep(
        Gamma_star_25, deg_to_kel(25), E_Gamma_star, Tleaf_K, R
    ) for Tleaf_K in Tleaf_K_values
]

plt.scatter(x_values, y)

# %%
# import ipywidgets as widgets
# %%

@interact
def A_n_values_check(
        c_a=391.0,
        e_a=635.2238437279115,
        g_bl=3478654.912462574,
        g_sto_0=20000.0,
        m=6,
        D_0=2.7,
        Gamma=14.333665361703794,
        Gamma_star=12.52438053902872,
        V_cmax=141.70322298110568,
        K_C=30.750578109134757,
        K_O=85.49206630570801,
        J=40.30643977341539,
        R_d=0.19323428021561032,
        e_sat_i=747.3221690916605,
        f_SW=1.0,
        f_LS=0,
        fO3_d=1.0,
        f_VPD=1.0
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
    model_options = ModelOptions(f_VPD_method=FVPDMethods.DISABLED,
                                 co2_concentration_balance_threshold=0.001,
                                 co2_concentration_max_iterations=50)
    final_state = co2_concentration_in_stomata_loop(
        constant_inputs=constant_inputs,
        model_options=model_options,
    )

    return final_state
# %%


@interact
def A_n_values_check(
    g_sto_0=20000,
    m=8.12,
    V_cmax_25=117.0,
    J_max_25=228.0,
    PARsun=224.758,
    PARshade=93.0,
    LAIsunfrac=0.81909,
    D_0=2.2,
    g_bv=3720747.41893,
    Tleaf_C=26.649947,
    f_SW=1.0,
    f_VPD=1.0,
    f_LS=1.0,
    fO3_d=0.99696062,
    eact=0.525388732,
):
    nL = 3
    LAI = 2.25
    final_state = ewert_leaf_pop(
        nL=nL,
        g_sto_0=g_sto_0,
        m=m,
        V_cmax_25=[V_cmax_25 for _ in range(nL)],
        J_max_25=[J_max_25 for _ in range(nL)],
        R_d_coeff=0.015,
        PARsun=[PARsun for _ in range(nL)],
        PARshade=[PARshade for _ in range(nL)],
        LAIsunfrac=[0.9, 0.8, 0.7],
        D_0=D_0,
        g_bv=g_bv,
        layer_lai_frac=[0.33, 0.33, 0.33],
        layer_lai=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        Tleaf_C=[Tleaf_C for _ in range(nL)],
        f_SW=f_SW,
        f_VPD=f_VPD,
        f_LS=f_LS,
        fO3_d=fO3_d,
        eact=eact,
        c_a=391.0,
        f_VPD_method=FVPDMethods.LEUNING,
    )

    return final_state.A_n


# %%


@interact
def A_n_values_check(
    g_sto_0=10.5300173629,
    m=8.12,
    V_cmax_25=117.91,
    J_max_25=228.0,
    PARsun=229.23417752214127,
    PARshade=58.4001157700821,
    LAIsunfrac=0.81909,
    D_0=2.2,
    Tleaf_C=27.11390,
    f_SW=1.0,
    f_VPD=1.0,
    f_LS=1.0,
    fO3_d=0.99696062,
    eact=0.521354,
    # Ext
    u=3.28309301,
):
    nL = 3
    LAI = 2.25
    g_bv = calc_g_bv(
        Lm=0.02,
        u=u,
        gas=GAS.H2O,
    )
    final_state = ewert_leaf_pop(
        nL=nL,
        g_sto_0=g_sto_0,
        m=m,
        V_cmax_25=[V_cmax_25 for _ in range(nL)],
        J_max_25=[J_max_25 for _ in range(nL)],
        R_d_coeff=0.015,
        PARsun=[PARsun for _ in range(nL)],
        PARshade=[PARshade for _ in range(nL)],
        LAIsunfrac=[0.9, 0.8, 0.7],
        D_0=D_0,
        g_bv=g_bv,
        layer_lai_frac=[0.33, 0.33, 0.33],
        layer_lai=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        Tleaf_C=[Tleaf_C for _ in range(nL)],
        f_SW=f_SW,
        f_VPD=f_VPD,
        f_LS=f_LS,
        fO3_d=fO3_d,
        eact=eact,
        c_a=391.0,
        f_VPD_method=FVPDMethods.LEUNING,
    )

    return final_state.A_n


# %%
df = pd.read_csv('tests/key_processes/photosynthesis/tmp/0.csv')
df.A_n_canopy.plot()

# %%
# @interact
# def func(
#         A_n=50,
#         c_a=391.0,
#         g_sto=20000,
#         g_bl=4594662.190354731,
# ):
#     return calc_CO2_supply(
#             A_n=A_n,
#             c_a=c_a,
#             g_sto=g_sto,
#             g_bl=g_bl,
#         )
A_n_min = -60
A_n_max = 50
g_sto_min = 10000
g_sto_max = 70000
x = np.linspace(A_n_min, A_n_max, 20)
y = np.linspace(g_sto_min, g_sto_max, 20)
xx, yy = np.meshgrid(x, y)
zz = np.zeros_like(xx)
for i in range(len(x)):
    for j in range(len(y)):
        zz[i, j] = calc_CO2_supply(
            A_n=x[i],
            c_a=391.0,
            g_sto=y[j],
            g_bl=4594662.190354731,
        )

# ax.contourf(xx, yy, zz, levels=20, cmap='viridis')
fig, ax = plt.subplots(figsize=(20, 10))
mat = ax.matshow(zz)
ax.set_xticks(ticks=np.arange(len(x)), labels=np.round(x, 0))
ax.set_yticks(ticks=np.arange(len(y)), labels=np.round(y, 0))
ax.set_xlabel('A_n')
ax.set_ylabel('g_sto')
ax.set_title('CO2 supply')

# add colorbar
cbar = plt.colorbar(mat)

plt.show()
# %%
np.diag(range(15)).shape


# %%
constant_inputs = CO2_Constant_Loop_Inputs(c_a=391.0, e_a=712.529686719428, g_bl=4594662.190354731, g_sto_0=10000.0, m=4, D_0=2.2, Gamma=57.89952113732494, Gamma_star=48.51315602386692,
                                           V_cmax=159.98275062496984, K_C=528.0415080745356, K_O=314.4035643006643, J=180.94373072965217, R_d=1.5998275062496985, e_sat_i=3670.6412787571467, f_SW=1.0, f_LS=1, fO3_d=1.0, f_VPD=1.0)
state = CO2_loop_State(c_i=0.0, c_i_diff=1, g_sto=10000.0, A_n=0.0, A_c=0.0,
                       A_p=0.0, A_j=0.0, A_n_limit_factor=None, f_VPD=None, iterations=0)
model_options = ModelOptions(f_VPD_method=FVPDMethods.DANIELSSON,
                             co2_concentration_balance_threshold=0.001, co2_concentration_max_iterations=5000)
iteration_state_s = []
while state.c_i_diff > model_options.co2_concentration_balance_threshold and state.iterations < model_options.co2_concentration_max_iterations:
    state = co2_concentration_in_stomata_iteration(
        constant_inputs, state, model_options)
    iteration_state_s.append(state)
len(iteration_state_s)

# %%
plt.plot([state.c_i_diff for state in iteration_state_s])
plt.scatter(np.arange(len(iteration_state_s)), [state.c_i_diff for state in iteration_state_s])
plt.title('CO2 concentration difference between iterations')
plt.xlabel('Iteration')
plt.ylabel('CO2 concentration difference')
plt.hlines(model_options.co2_concentration_balance_threshold, 0, len(iteration_state_s), colors='r', linestyles='dashed')
print(min([state.c_i_diff for state in iteration_state_s]), max([state.c_i_diff for state in iteration_state_s]))
# %%
plt.plot([state.A_n for state in iteration_state_s])
plt.scatter(np.arange(len(iteration_state_s)), [state.A_n for state in iteration_state_s])
plt.title('Net assimilation')
plt.xlabel('Iteration')
plt.ylabel('Net assimilation')


# %%
# c_i
plt.plot(np.arange(len(iteration_state_s)), [state.c_i for state in iteration_state_s])
plt.scatter(np.arange(len(iteration_state_s)), [state.c_i for state in iteration_state_s])
plt.xlabel('Iteration')
plt.ylabel('c_i')
plt.title('CO2 concentration in stomata')

# %%
# g_sto
plt.plot(np.arange(len(iteration_state_s)), [state.g_sto for state in iteration_state_s])
plt.scatter(np.arange(len(iteration_state_s)), [state.g_sto for state in iteration_state_s])
plt.xlabel('Iteration')
plt.ylabel('g_sto')
plt.title('Stomatal conductance')

# %%


# ax.contourf(xx, yy, zz, levels=20, cmap='viridis')
fig, ax = plt.subplots(figsize=(20, 10))
mat = ax.matshow(zz)
ax.set_xticks(ticks=np.arange(len(x)), labels=np.round(x, 0))
ax.set_yticks(ticks=np.arange(len(y)), labels=np.round(y, 0))
ax.set_xlabel('A_n')
ax.set_ylabel('g_sto')
ax.set_title('CO2 supply')

# # add colorbar
cbar = plt.colorbar(mat)

xb = [state.A_n for state in iteration_state_s]
yb = [state.g_sto for state in iteration_state_s]
cb = [state.c_i_diff for state in iteration_state_s]
print(min(xb), max(xb), min(yb), max(yb))
# Scale x and y to between A_n_min and A_n_max and g_sto_min and g_sto_max
xb_scaled = len(x) * (np.array(xb) - A_n_min) / (A_n_max - A_n_min)
yb_scaled = len(y) * (np.array(yb) - g_sto_min) / (g_sto_max - g_sto_min)

sctc = ax.scatter(xb_scaled, yb_scaled, c=cb, cmap='Reds')
cbar = plt.colorbar(sctc)

# %%
xb
# %%


plt.show()

# %%

@interact(
    max_iterations=widgets.IntSlider(min=1, max=100, step=1, value=50),
    threshold=widgets.FloatSlider(min=0.001, max=1, step=0.001, value=0.01),
    g_bl=widgets.FloatSlider(min=1000000, max=10000000, step=100000, value=4594662.190354731),
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
        c_a=c_a, e_a=e_a, g_bl=g_bl, g_sto_0=g_sto_0, m=m, D_0=D_0, Gamma=Gamma, Gamma_star=Gamma_star,
        V_cmax=V_cmax, K_C=K_C, K_O=K_O, J=J, R_d=R_d, e_sat_i=e_sat_i,
        f_SW=f_SW, f_LS=f_LS, fO3_d=fO3_d, f_VPD=f_VPD)
    state = CO2_loop_State(c_i=0.0, c_i_diff=1, g_sto=10000.0, A_n=0.0, A_c=0.0,
                        A_p=0.0, A_j=0.0, A_n_limit_factor=None, f_VPD=None, iterations=0)
    model_options = ModelOptions(f_VPD_method=FVPDMethods.DANIELSSON,
                                co2_concentration_balance_threshold=threshold, co2_concentration_max_iterations=max_iterations)
    iteration_state_s = []
    while state.c_i_diff > model_options.co2_concentration_balance_threshold and state.iterations < model_options.co2_concentration_max_iterations:
        state = co2_concentration_in_stomata_iteration(
            constant_inputs, state, model_options)
        iteration_state_s.append(state)
    len(iteration_state_s)
    plt.close()
    plt.plot([state.c_i_diff for state in iteration_state_s])
    plt.scatter(np.arange(len(iteration_state_s)), [state.c_i_diff for state in iteration_state_s])
    plt.title('CO2 concentration difference between iterations')
    plt.xlabel('Iteration')
    plt.ylabel('CO2 concentration difference')
    plt.hlines(model_options.co2_concentration_balance_threshold, 0, len(iteration_state_s), colors='r', linestyles='dashed')
    # print(min([state.c_i_diff for state in iteration_state_s]), max([state.c_i_diff for state in iteration_state_s]))

# %%

@interact(
    max_iterations=widgets.IntSlider(min=1, max=100, step=1, value=50),
    threshold=widgets.FloatSlider(min=0.001, max=1, step=0.001, value=0.01),
    c_a=widgets.FloatSlider(min=300, max=500, step=1, value=391.0),
    e_a=widgets.FloatSlider(min=600, max=800, step=1, value=1013.8537256970482),
    g_bl=widgets.FloatSlider(min=10, max=10000000, step=100000, value=2684466.893155262),
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
        c_a=c_a, e_a=e_a, g_bl=g_bl, g_sto_0=g_sto_0, m=m, D_0=D_0, Gamma=Gamma, Gamma_star=Gamma_star,
        V_cmax=V_cmax, K_C=K_C, K_O=K_O, J=J, R_d=R_d, e_sat_i=e_sat_i,
        f_SW=f_SW, f_LS=f_LS, fO3_d=fO3_d, f_VPD=f_VPD)
    state = CO2_loop_State(c_i=0.0, c_i_diff=1, g_sto=10000.0, A_n=0.0, A_c=0.0,
                        A_p=0.0, A_j=0.0, A_n_limit_factor=None, f_VPD=None, iterations=0)
    model_options = ModelOptions(f_VPD_method=FVPDMethods.DANIELSSON,
                                co2_concentration_balance_threshold=threshold, co2_concentration_max_iterations=max_iterations)
    iteration_state_s = []
    while state.c_i_diff > model_options.co2_concentration_balance_threshold and state.iterations < model_options.co2_concentration_max_iterations:
        state = co2_concentration_in_stomata_iteration(
            constant_inputs, state, model_options)
        iteration_state_s.append(state)
    len(iteration_state_s)
    plt.close()
    fig, axs = plt.subplots(1,5, figsize=(20, 5))
    x = np.arange(len(iteration_state_s))
    c = [cc%2 for cc in range(len(iteration_state_s))]
    axs[0].scatter(x,[state.c_i_diff for state in iteration_state_s], c=c)
    axs[0].plot(x,[state.c_i_diff for state in iteration_state_s])
    axs[0].scatter(np.arange(len(iteration_state_s)), [state.c_i_diff for state in iteration_state_s], c=c)
    axs[0].set_title('CO2 concentration diff')
    axs[0].set_xlabel('Iteration')
    axs[0].set_ylabel('CO2 concentration difference')
    axs[0].hlines(model_options.co2_concentration_balance_threshold, 0, len(iteration_state_s), colors='r', linestyles='dashed')
    axs[0].text(0.5,100, iteration_state_s[-1].A_n_limit_factor)


    axs[1].scatter(x,[state.A_n for state in iteration_state_s], c=c)
    axs[1].set_title('A_n')
    axs[2].scatter(x,[state.g_sto for state in iteration_state_s], c=c)
    axs[2].set_title('gsto')
    axs[3].scatter(x,[state.c_i for state in iteration_state_s], c=c)
    axs[3].set_title('CO2')
    if len(iteration_state_s) == max_iterations:
        axs[0].text(0.5, 50, 'Did not converge', c="red")

    axs[4].scatter(x,[state.f_VPD for state in iteration_state_s])
    axs[4].set_title('f_VPD')
# f_VPD

    # print(min([state.c_i_diff for state in iteration_state_s]), max([state.c_i_diff for state in iteration_state_s]))
# %%
iteration_state_s[-1]

# %%
c_a = 391
ind = -2
A_n = iteration_state_s[ind].A_n
g_sto = iteration_state_s[ind].g_sto
g_bl = 2684466.893155262
co2_supply = c_a - ((A_n * (1 / g_sto + 1.37 / g_bl)) * 1e6)
c_i_in = iteration_state_s[ind].c_i
c_i_diff = iteration_state_s[ind].c_i_diff
co2_supply, c_i_in, c_i_in - (c_i_in - co2_supply) / 2, c_i_diff
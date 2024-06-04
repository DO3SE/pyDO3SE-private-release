# %%
from matplotlib import pyplot as plt
from pyDO3SE.plugins.gsto.ewert.ewert import *
from pyDO3SE.plugins.gsto.ewert.ewert_helpers import *


t_lem_constant = 0.15
t_lse_constant = 0.33

t_l = 800
t_lem = t_l * t_lem_constant
t_lma = t_l - t_lem
t_lse = t_l - (t_lma * t_lse_constant)
t_lep = t_l - (t_lem + t_lse)  # check this


# %%
def get_vals(O3):
    return [co2_concentration_in_stomata_iteration(
        CO2_Constant_Loop_Inputs(
            c_a=391.0,
            e_a=1000.0,
            g_bl=1469999.0,
            g_sto_0=20000,
            m=8.12,
            D_0=2.27 * 1e3,
            O3up_acc=i * O3,
            O3up=21.0,
            fO3_d_prev=0.89,
            td_dd=i,
            gamma_1=0.06,
            gamma_2=0.0045,
            gamma_3=0.5,
            is_daylight=True,
            t_lse_constant=t_lse_constant,
            t_l_estimate=t_l,
            t_lem=t_lem,
            t_lep=t_lep,
            t_lse=t_lse,
            t_lma=t_lma,
            Gamma=34.277,
            Gamma_star=32.95,
            V_cmax=119.0,
            K_C=234.42,
            K_O=216.75,
            J=300.36,
            R_d=0.32,
            e_sat_i=2339.05,
        ),
        CO2_loop_State(
            c_i=0.0,
            c_i_diff=0,
            g_sto=20000,
            A_n=0,
            fO3_d=1.0,
            fO3_h=1.0,
        )

    ) for i in range(800)]


# %%
for i in range(0, 100):
    plt.plot([o.f_LS for o in get_vals(i / 100)])
# plt.plot([o.f_LS for o in get_vals(0.1)])
# %%

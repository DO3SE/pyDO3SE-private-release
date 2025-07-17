# %%

from itertools import product
from matplotlib import pyplot as plt
import numpy as np

from pyDO3SE.plugins.gsto.ewert.ewert import (
    FVPDMethods,
    ewert_leaf_pop,
    ewert_leaf_pop_alt,
)

negative_result_inputs = {
    "nL": 1,
    "g_sto_0": 10000.0,
    "m": 4,
    "V_cmax_25": [137.0],
    "J_max_25": [228.0],
    "R_d_coeff": 0.01,
    "PARsun": [
        186.5,
    ],
    "PARshade": [
        15.70,
    ],
    "LAIsunfrac": [
        1.0,
    ],
    "D_0": 2.2,
    "g_bv": 3766811.8338,
    "layer_lai_frac": [
        1.0,
    ],
    "layer_lai": [
        0.564,
    ],
    "Tleaf_C": [
        27.113,
    ],
    "f_SW": 1.0,
    "f_VPD": 1.0,
    "f_LS": 1,
    "fO3_d": 0.996,
    "eact": 0.521,
    "c_a": 391.0,
    "f_VPD_method": FVPDMethods.DANIELSSON,
}

# TODO: Get inputs from negative run
# g_sto_0_list = [10000.0]
# m_list = [4]
# V_cmax_25_list = [ 137.0]
# J_max_25_list = [228.0]
# R_d_coeff_list = [0.01]
# PARsun_list = [186.0]
# PARshade_list = [15.7]
# D_0_list = [2.2]
# g_bv_list = [3766811.956]
# Tleaf_C_list = [27.0]
# eact_list = [0.5, 1.0]
# f_LS_list = [1.0]
# f_O3d_list = [1.0]

g_sto_0_list = [10000]
m_list = [4]
V_cmax_25_list = [40, 137, 180.0]
J_max_25_list = [228, 400]
R_d_coeff_list = [0.015]
PARsun_list = [200, 1000]
PARshade_list = [15, 120]
D_0_list = [2.2]
g_bv_list = [1469999.0, 3766811.956]
Tleaf_C_list = [15.0, 27.0]
eact_list = [0.5, 1.0]
f_LS_list = [0.4, 1.0]
f_O3d_list = [0.4, 1.0]

# g_sto_0_list = [10000]
# m_list = [4]
# V_cmax_25_list = [137]
# J_max_25_list = [228]
# R_d_coeff_list = [0.015]
# PARsun_list = [200]
# PARshade_list = [10]
# D_0_list = [2.2]
# g_bv_list = [2758244.956]
# Tleaf_C_list = [15.0,20.0,25.0]
# eact_list = [0.6]

# g_sto_0_list = [10000]
# m_list = [4]
# V_cmax_25_list = [40,137, 180.0]
# J_max_25_list = [228, 400]
# R_d_coeff_list = [0.015]
# PARsun_list = [200,600,1000]
# PARshade_list = [10, 80, 120]
# D_0_list = [2.2]
# g_bv_list = [1469999.0, 2758244.956]
# Tleaf_C_list = [0.0,15.0,25.0]
# eact_list = [0.6,1.0]

LAI = 0.5
nL = 1


print(
    len(
        list(
            product(
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
                f_LS_list,
                f_O3d_list,
            )
        )
    )
)

# %%
print(nL)
results = []
inputs = []
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
    f_LS,
    f_O3d,
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
    f_LS_list,
    f_O3d_list,
):
    kwargs_in = dict(
        nL=nL,
        g_sto_0=g_sto_0,
        m=m,
        V_cmax_25=[V_cmax_25 for _ in range(nL)],
        J_max_25=[J_max_25 for _ in range(nL)],
        R_d_coeff=R_d_coeff,
        PARsun=[PARsun for _ in range(nL)],
        PARshade=[PARshade for _ in range(nL)],
        LAIsunfrac=[1.0],
        D_0=D_0,
        g_bv=g_bv,
        # leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        layer_lai_frac=[1.0],
        layer_lai=[LAI],
        Tleaf_C=[Tleaf_C for _ in range(nL)],
        f_SW=1,
        f_VPD=1.0,
        f_LS=f_LS,
        fO3_d=f_O3d,
        eact=eact,
        c_a=391.0,
        f_VPD_method=FVPDMethods.LEUNING,
    )
    inputs.append(kwargs_in)
    out = ewert_leaf_pop(
        **kwargs_in,
        # nL=nL,
        # g_sto_0=g_sto_0,
        # m=m,
        # V_cmax_25=[V_cmax_25 for _ in range(nL)],
        # J_max_25=[J_max_25 for _ in range(nL)],
        # R_d_coeff=R_d_coeff,
        # PARsun=[PARsun for _ in range(nL)],
        # PARshade=[PARshade for _ in range(nL)],
        # LAIsunfrac=[1.0],
        # D_0=D_0,
        # g_bv=g_bv,
        # # leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        # layer_lai_frac=[1.0],
        # layer_lai=[LAI],
        # Tleaf_C=[Tleaf_C for _ in range(nL)],
        # f_SW=1,
        # f_VPD=1.0,
        # f_LS=f_LS,
        # fO3_d=f_O3d,
        # eact=eact,
        # c_a=391.0,
        # f_VPD_method=FVPDMethods.LEUNING,
    )

    out_alt = ewert_leaf_pop_alt(
        **kwargs_in,
        # nL=nL,
        # g_sto_0=g_sto_0,
        # m=m,
        # V_cmax_25=[V_cmax_25 for _ in range(nL)],
        # J_max_25=[J_max_25 for _ in range(nL)],
        # R_d_coeff=R_d_coeff,
        # PARsun=[PARsun for _ in range(nL)],
        # PARshade=[PARshade for _ in range(nL)],
        # LAIsunfrac=[1.0],
        # D_0=D_0,
        # g_bv=g_bv,
        # # leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
        # layer_lai_frac=[1.0],
        # layer_lai=[LAI],
        # Tleaf_C=[Tleaf_C for _ in range(nL)],
        # f_SW=1,
        # f_VPD=1.0,
        # f_LS=f_LS,
        # fO3_d=f_O3d,
        # eact=eact,
        # c_a=391.0,
        # f_VPD_method=FVPDMethods.LEUNING,
    )
    results.append([out, out_alt])
print(len(results))

print(kwargs_in)
# %%
# Run with negative

# out = ewert_leaf_pop(
#     **negative_result_inputs
# )
# out_alt = ewert_leaf_pop_alt(
#     **negative_result_inputs
# )
# results.append([out, out_alt])

# %%
print(kwargs_in)
print(negative_result_inputs)
# %%
#  Plot compare A_n
ckey = "eact"
x = [r[0].A_n for r in results]
y = [r[1].A_n for r in results]
c = [k[ckey] for k in inputs]

min_val = min(min(x), min(y))
max_val = max(max(x), max(y))
print(min_val, max_val)
# plt.plot([min_val, min_val], [max_val, max_val])
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
cx = ax.scatter(x, y, c=c)
ax.set_xlabel("Original")
ax.set_ylabel("Alternative")
ax.set_title("A_n")
ax.plot([0, max_val], [0, max_val], color="k", linestyle="-", linewidth=2)
ax.set_aspect("equal")
cbar = plt.colorbar(cx)
cbar.set_label(ckey)

# %%
#  Plot compare gsv

x = [np.mean(r[0].g_sv_per_layer) for r in results]
y = [np.mean(r[1].g_sv_per_layer) for r in results]
c = [k[ckey] for k in inputs]

min_val = min(min(x), min(y))
max_val = max(max(x), max(y))
print(min_val, max_val)

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
cx = ax.scatter(x, y, c=c)
ax.set_xlabel("Original")
ax.set_ylabel("Alternative")
ax.set_title("mean g_sv_per_layer")
ax.plot([0, max_val], [0, max_val], color="k", linestyle="-", linewidth=2)
ax.set_aspect("equal")

cbar = plt.colorbar(cx)
cbar.set_label(ckey)

# %%



out_alt = ewert_leaf_pop(
    **inputs
    # nL=nL,
    # g_sto_0=g_sto_0,
    # m=4,
    # V_cmax_25=[V_cmax_25 for _ in range(nL)],
    # J_max_25=[J_max_25 for _ in range(nL)],
    # R_d_coeff=R_d_coeff,
    # PARsun=[PARsun for _ in range(nL)],
    # PARshade=[PARshade for _ in range(nL)],
    # LAIsunfrac=[0.9, 0.8, 0.7],
    # D_0=D_0,
    # g_bv=g_bv,
    # # leaf_pop_distribution=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
    # layer_lai_frac=[0.33, 0.33, 0.33],
    # layer_lai=[LAI * 0.33, LAI * 0.33, LAI * 0.33],
    # Tleaf_C=[Tleaf_C for _ in range(nL)],
    # f_SW=1,
    # f_VPD=1.0,
    # f_LS=0.0,
    # fO3_d=1.0,
    # eact=eact,
    # c_a=391.0,
    # f_VPD_method=FVPDMethods.LEUNING,
)
from pprint import pprint

pprint(out_alt.__dict__, indent=4)

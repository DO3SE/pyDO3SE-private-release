# %%
from helpers.list_helpers import offset
from helpers.util import PLF_value
from matplotlib import pyplot as plt
# %%
LAI = 117
LAI_1 = 21
SGS = 118
EGS = 200


# %%
# FORTRAN
fphen_limA = 0
fphen_limB = 0
fphen_a = 0.1
fphen_b = 0
fphen_c = 1.0
fphen_d = 0
fphen_e = 0.1
fphen_1 = 0
fphen_2 = 0
fphen_3 = 0
fphen_4 = 45
leaf_fphen_method = "day PLF"
leaf_fphen_a = 0.8
leaf_fphen_b = 1.0
leaf_fphen_c = 0.2
leaf_fphen_1 = 15
leaf_fphen_2 = 40


# %%
Astart = 115
Aend = 200

# PYTHON


def calc_leaf_fphen(dd):
    gs_values = [
        Astart,
        Astart,
        (Astart + leaf_fphen_1),
        (Aend - leaf_fphen_2),
        Aend,
        Aend,
    ]
    fphen_values = [
        0.0,
        leaf_fphen_a,
        leaf_fphen_b,
        leaf_fphen_b,
        leaf_fphen_c,
        0.0,
    ]

    # Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    gs_offset = offset(gs_values, float(Astart), 365.0)

    gs_values_in_size_order = all([a <= b for a, b in zip(gs_offset[0: 5], gs_offset[1: 6])])
    if not gs_values_in_size_order:
        raise ValueError("fphen_simple_PLF: points not in order")

    dd_adj = dd - Astart if Astart - dd < 0 else dd - Astart + 365

    # Lookup value in PLF
    func = [gs_offset, fphen_values]
    leaf_fphen = PLF_value(func, float(dd_adj))
    return leaf_fphen


results = [calc_leaf_fphen(dd) for dd in range(365)]
plt.plot(results)
# %%


def calc_leaf_fphen(dd):
    if (dd < Astart or dd > Aend):
        leaf_fphen = 0
    elif (dd < (Astart + leaf_fphen_1)):
        leaf_fphen = leaf_fphen_a + (leaf_fphen_b - leaf_fphen_a) * (dd - Astart) / leaf_fphen_1
    elif (dd < (Aend - leaf_fphen_2)):
        leaf_fphen = leaf_fphen_b
    elif (dd <= Aend):
        leaf_fphen = leaf_fphen_c + (leaf_fphen_b - leaf_fphen_c) * (Aend - dd) / leaf_fphen_2
    return leaf_fphen


results = [calc_leaf_fphen(dd) for dd in range(365)]
plt.plot(results)

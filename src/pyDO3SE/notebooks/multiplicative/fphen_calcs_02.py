# %%
from helpers.util import PLF_value
from helpers.list_helpers import offset
from matplotlib import pyplot as plt

# %%
# FORTRAN ===============
# TODO: Check this matches Fortran output


def Calc_fphen(
    f_phen_1: int,
    f_phen_2: int,
    f_phen_3: int,
    f_phen_4: int,
    f_phen_limA: int,
    f_phen_limB: int,
    f_phen_a: float,
    f_phen_b: float,
    f_phen_c: float,
    f_phen_d: float,
    f_phen_e: float,
    SGS: int,
    EGS: int,
    dd: int,
):
    dd = dd - 117
    SGS = SGS - 117
    EGS = EGS - 117
    if dd < SGS or dd > EGS:
        f_phen = 0.0
    # elif dd < (SGS + f_phen_1):
    #     f_phen = f_phen_a + (f_phen_b - f_phen_a) * (dd - SGS) / f_phen_1
    # elif dd < f_phen_limA:
    #     f_phen = f_phen_b
    # elif dd < (f_phen_limA + f_phen_2):
    #     f_phen = f_phen_b - (f_phen_b - f_phen_c) * (dd - f_phen_limA) / f_phen_2
    # elif dd < (f_phen_limB - f_phen_3):
    #     f_phen = f_phen_c
    # elif dd < f_phen_limB:
    #     f_phen = f_phen_d - (f_phen_d - f_phen_c) * (f_phen_limB - dd) / f_phen_3
    elif dd < (EGS - f_phen_4):
        f_phen = f_phen_d
    elif dd <= EGS:
        f_phen = f_phen_e + (f_phen_d - f_phen_e) * (EGS - dd) / f_phen_4
    return f_phen


def calc_fphen_demo(dd): return Calc_fphen(
    f_phen_1=0,
    f_phen_2=1,
    f_phen_3=1,
    f_phen_4=45,
    f_phen_limA=0,
    f_phen_limB=0,
    f_phen_a=0.1,
    f_phen_b=1.0,
    f_phen_c=1.0,
    f_phen_d=1.0,
    f_phen_e=0.1,
    SGS=118,
    EGS=210,
    dd=dd,
)


# %%
# Calc for day 117
out = calc_fphen_demo(
    dd=117,
)
assert calc_fphen_demo(dd=117) == 0.0
assert calc_fphen_demo(dd=118) == 1.0
assert calc_fphen_demo(dd=165) == 1.0
assert calc_fphen_demo(dd=166) == 0.98
assert calc_fphen_demo(dd=210) == 0.1
assert calc_fphen_demo(dd=211) == 0.0

# %%
plt.plot([calc_fphen_demo(dd=dd) for dd in range(0, 365)])


# %%

# ======== PYTHON COMPLEX ======================================================

def f_phen_complex_PLF(
    f_phen_1: int,
    f_phen_2: int,
    f_phen_3: int,
    f_phen_4: int,
    f_phen_limA: int,
    f_phen_limB: int,
    f_phen_a: float,
    f_phen_b: float,
    f_phen_c: float,
    f_phen_d: float,
    f_phen_e: float,
    SGS: int,
    EGS: int,
    dd: int,
):
    # TODO: Python  list length is short 1
    # TODO: Check these are correct below
    use_complex_components = f_phen_limA and f_phen_limB

    complex_middle_gs_values = [
        f_phen_limA,
        f_phen_limA + f_phen_2,
        f_phen_limB - f_phen_3,
        f_phen_limB
    ] if use_complex_components else [
        SGS + f_phen_1,
        SGS + f_phen_1,
        EGS - f_phen_4,
        EGS - f_phen_4,
    ]

    complex_middle_fphen_values = [
        f_phen_b,
        f_phen_c,
        f_phen_c,
        f_phen_d,
    ] if use_complex_components else [
        f_phen_b, f_phen_b, f_phen_d, f_phen_d
    ]

    gs_values = [
        SGS,
        SGS,
        SGS + f_phen_1,
        *complex_middle_gs_values,
        EGS - f_phen_4,
        EGS,
        EGS,
    ]
    fphen_values = [
        0.0,
        f_phen_a,
        f_phen_b,
        *complex_middle_fphen_values,
        f_phen_d,
        f_phen_e,
        0.0,
    ]
    print(gs_values)
    print(fphen_values)
    # Correct to here!!!

    # Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    gs_offset = offset(gs_values, float(SGS), 365.0)
    print(gs_offset)
    gs_values_in_size_order: bool = all([a <= b for a, b in zip(gs_offset[0: 9], gs_offset[1: 10])])

    print(gs_values_in_size_order)

    if not gs_values_in_size_order:
        print(gs_values)
        raise ValueError("f_phen_simple_PLF: points not in order")

    dd_adj = dd - SGS if SGS - dd <= 0 else dd - SGS + 365

    # Lookup value in PLF
    func = [gs_offset, fphen_values]

    # fphen
#   F 0.0, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.0
#   P 0.0, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.0
# gs
#   F 0.0, 0.0, 0.0, 247.0, 248.0, 246.0, 247.0, 47.0, 92.0, 92.0
  # P 0.0, 0.0, 0.0, 0.0,   0.0,  47.0,  47.0,  47.0, 92.0, 92.0
    # dd_adj = 364 # is correct
    # dd = 118
    print(dd_adj)
    print(func)
    f_phen = PLF_value(func, float(dd_adj))
    print(f_phen)
    print("end -------------")
    return f_phen


def calc_fphen_demo_p(dd): return f_phen_complex_PLF(
    f_phen_1=0,
    f_phen_2=1,
    f_phen_3=1,
    f_phen_4=45,
    f_phen_limA=0,
    f_phen_limB=0,
    f_phen_a=0.1,
    f_phen_b=1.0,
    f_phen_c=1.0,
    f_phen_d=1.0,
    f_phen_e=0.1,
    SGS=118,
    EGS=210,
    dd=dd,
)


assert calc_fphen_demo_p(dd=118) == 1.0  # ddadj = 0 (365)

# %%
plt.plot([calc_fphen_demo(dd=dd) - calc_fphen_demo_p(dd=dd) for dd in range(0, 365)])
# %%
plt.plot([calc_fphen_demo_p(dd=dd) for dd in range(0, 365)])

# %%
# print([calc_fphen_demo_p(dd=dd) for dd in range(117, 120)])
# assert calc_fphen_demo_p(dd=117) == 0.0 # ddadj = 364
# assert calc_fphen_demo_p(dd=118) == 1.0 # ddadj = 0 (365)
# assert calc_fphen_demo_p(dd=165) == 1.0
# assert calc_fphen_demo_p(dd=166) == 0.98
assert calc_fphen_demo_p(dd=210) == 0.1
# assert calc_fphen_demo_p(dd=211) == 0.0

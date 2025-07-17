# %%
%load_ext autoreload
%autoreload 2
from do3se_met import irradiance as met_irrad_helpers
met_irrad_helpers.calc_PAR_sun_shade_farq_b
# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# %%

# Idrctt 0-400
# Idfuse 0-175
# sinb - 0-0.8

met_irrad_helpers.calc_PAR_sun_shade_farq_b(
    100,
    100,
    0.8,
    0.5,
    LAI_c=5,
)
# %%

IdrcttList = [

0,
0,
0,
0,
0,
0,
0,
0,
0.053079056061155,
0.319355870292217,
0.728195873884611,
7.77938525349878,
21.7941410192125,
24.8608681806006,
48.436546934961,
69.3580481898954,
88.2661527697272,
0.458873074581877,
0,
0,
0,
0,
0,
0,

]
IdfuseList = [
    0,0,0,0,5,10,10,25,50,75,100,100,100,100,100,100,100,75,50,25,0,0,0,0
]
sinBList = [
    0,
    0,
    0,
    0,
    0,
    0,
    0.132977743690365,
    0.349302674899292,
    0.539259748045431,
    0.689903702512779,
    0.790968401755291,
    0.835566453536904,
    0.820658574330465,
    0.747260711467457,
    0.620374808006292,
    0.448647928578978,
    0.2437829762264,
    0.019741158896044,
    0,
    0,
    0,
    0,
    0,
    0,
]

assert len(IdrcttList) == 24, f"len: {len(IdrcttList)}"
assert len(IdfuseList) == 24, f"len: {len(IdfuseList)}"
assert len(sinBList) == 24, f"len: {len(sinBList)}"

# %%
par_sun_vals = []
par_shade_vals = []
for Idrctt, Idfuse, sinB in zip(IdrcttList, IdfuseList, sinBList):
    parsun, parshade = (
        met_irrad_helpers.calc_PAR_sun_shade_farq_b(
            Idrctt,
            (Idfuse *0) + 1,
            sinB=sinB,
            # sinB=1.5708 - sinB,
            cosA=0.5,
            LAI_c=99
        )
    )
    par_sun_vals.append(parsun)
    par_shade_vals.append(parshade)

plt.plot(par_sun_vals, label="par_sun")
# plt.plot(par_shade_vals, label="par_shade")
plt.plot(IdrcttList, label="Idrctt")
# plt.plot((Idfuse *0) + 1, label="Idfuse")
plt.legend()
# %%


from math import cos
cos(0.1)
# %%
fig, ax = plt.subplots(1, figsize=(15, 5))
for LAI_c in range(1, 10, 1):
    par_sun_vals = []
    par_shade_vals = []
    for Idrctt, Idfuse, sinB in zip(IdrcttList, IdfuseList, sinBList):
        parsun, parshade = (
            met_irrad_helpers.calc_PAR_sun_shade_farq_b(
                Idrctt,
                (Idfuse *0) + 1,
                sinB=sinB,
                # sinB=1.5708 - sinB,
                cosA=0.5,
                LAI_c=LAI_c
            )
        )
        par_sun_vals.append(parsun)
        par_shade_vals.append(parshade)

    ax.plot(par_sun_vals, label="par_sun")
    ax.plot(par_shade_vals, label="par_shade")
    # ax.plot(IdrcttList, label="Idrctt")
    # ax.plot((Idfuse *0) + 1, label="Idfuse")
    ax.legend()
# %%
fig, ax = plt.subplots(1, figsize=(15, 5))
par_sun_vals = []
par_shade_vals = []
LAI_values = np.array(
    [
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                1e-07
            ],
            [
                0.01
            ],
            [
                0.02
            ],
            [
                0.05
            ],
            [
                0.12
            ],
            [
                0.18
            ],
            [
                0.26
            ],
            [
                0.23
            ],
            [
                0.11
            ],
            [
                0.02
            ]
        ]
) * 5
# LAI_c equals cumulative LAI
LAI_c_vals = np.cumsum(list(reversed(LAI_values)))


Idrctt, Idfuse, sinB = 50, 100, 0.8
for LAI_c in LAI_c_vals:
    parsun, parshade = (
        met_irrad_helpers.calc_PAR_sun_shade_farq_b(
            Idrctt,
            Idfuse,
            sinB=sinB,
            # sinB=1.5708 - sinB,
            cosA=0.5,
            LAI_c=LAI_c
        )
    )
    par_sun_vals.append(parsun)
    par_shade_vals.append(parshade)

ax.plot(par_sun_vals, label="par_sun")
ax.plot(par_shade_vals, label="par_shade")#
ax.plot(LAI_c_vals, label="LAI_c")
# ax.plot(IdrcttList, label="Idrctt")
# ax.plot((Idfuse *0) + 1, label="Idfuse")
ax.legend()

# %%
iL = 0
nL = 10
list(range(iL + 1, nL))
# %%

"""Demonstrates an overview of the multilayered approach.

This an experiment with using simplified versions of all the processes to see how
we combine them all in a multilayer approach.
"""

# %%
from matplotlib import pyplot as plt
from pyDO3SE.util.Objects import Field
from pyDO3SE.Analysis.charts import multi_series_annual_graph


nL = 5
LAI = [0.01 for i in range(nL)]
LAI_c = [sum(LAI[0:i]) for i in range(nL)]

LAI_c
# %%

# 1. Calculate ozone deposition to each layer


def calc_ozone_deposition_per_layer(O3_in, LAI, LAI_c):
    O3_per_layer = [la * (O3_in - (O3_in * lac)) for la, lac in zip(LAI, LAI_c)]
    return O3_per_layer


O3_measured = 10
O3_per_layer = calc_ozone_deposition_per_layer(O3_measured, LAI, LAI_c)
O3_per_layer
# %%

# 2. Calculate the PARsun and PARshade per layer


def calc_PAR_sun_per_layer(PAR, LAI):
    PARsun_per_layer = [lai * PAR for lai in LAI]
    return PARsun_per_layer


def calc_PAR_shade_per_layer(PAR, LAI):
    sun_shade_fract = 0.6  # fraction of par that is shaded par
    PARsun_per_layer = [sun_shade_fract * lai * PAR for lai in LAI]
    return PARsun_per_layer


# %%
LAI_sunlit = [lai - lai_c / nL for lai, lai_c in zip(LAI, LAI_c)]
LAI_shaded = [lai_c / nL for lai, lai_c in zip(LAI, LAI_c)]
LAI, LAI_sunlit
# %%
PAR_measured = 600
PARsun_per_layer = calc_PAR_sun_per_layer(PAR_measured, LAI_sunlit)
PARshade_per_layer = calc_PAR_sun_per_layer(PAR_measured, LAI_shaded)
PARsun_per_layer, PARshade_per_layer

# %%

PARsun_per_layer, PARshade_per_layer, O3_per_layer
# %%

# 3. Calculate Photosynthesis


def calc_photosynthesis_per_layer(PAR, O3, gsto_in):
    # print('----')
    # print(PAR, O3)
    fLS = 1 - O3 + 0.04 if O3 > 0.04 else 1
    A_c = 1 * fLS
    A_j = PAR / 6.0 if PAR > 3 else 0
    A_n = min(A_c, A_j)
    # print(A_c, A_j, A_n)
    gsto = gsto_in * A_n
    # print('----')
    return A_n, gsto


gsto_per_layer = [0 for iL in range(nL)]
for iL in range(nL):
    gsto_in = 1000
    sun_values = calc_photosynthesis_per_layer(PARsun_per_layer[iL], O3_per_layer[iL], gsto_in)
    shade_values = calc_photosynthesis_per_layer(PARshade_per_layer[iL], O3_per_layer[iL], gsto_in)
    gsto_per_layer[iL] = sun_values[1] + shade_values[1]
    print(sun_values, shade_values)
gsto_per_layer


# %%
def multi_layer_simple(nL, LAI, O3_measured, PAR_measured):
    LAI_c = [sum(LAI[0:i]) for i in range(nL)]
    LAI_sunlit = [lai - lai_c / nL for lai, lai_c in zip(LAI, LAI_c)]
    LAI_shaded = [lai_c / nL for lai, lai_c in zip(LAI, LAI_c)]
    LAI, LAI_sunlit
    PARsun_per_layer = calc_PAR_sun_per_layer(PAR_measured, LAI_sunlit)
    PARshade_per_layer = calc_PAR_sun_per_layer(PAR_measured, LAI_shaded)

    O3_per_layer = calc_ozone_deposition_per_layer(O3_measured, LAI, LAI_c)

    gsto_per_layer = [0 for iL in range(nL)]
    for iL in range(nL):
        gsto_in = 1000
        sun_values = calc_photosynthesis_per_layer(PARsun_per_layer[iL], O3_per_layer[iL], gsto_in)
        shade_values = calc_photosynthesis_per_layer(
            PARshade_per_layer[iL], O3_per_layer[iL], gsto_in)
        gsto_per_layer[iL] = sun_values[1] + shade_values[1]
        # print(sun_values, shade_values)
    return gsto_per_layer


# %%
LAI = [0.01 for i in range(40)]
data = []

data.append(multi_layer_simple(24, LAI, 100, 800))
data.append(multi_layer_simple(24, LAI, 100, 600))
data.append(multi_layer_simple(24, LAI, 100, 400))
data.append(multi_layer_simple(24, LAI, 100, 350))
data.append(multi_layer_simple(24, LAI, 100, 300))

multi_series_annual_graph(
    'gsto',
    data,
    [i for i in range(len(data))],
    Field(1, float),
    0, 1,
)
plt.show()

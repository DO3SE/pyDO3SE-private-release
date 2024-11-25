# %% [markdown]
"""Compare the old and new multi-layer ozone deposition setups.

See email chain on 19/06/24
"""

# %%
import pandas as pd
import numpy as np
from scipy.linalg.lapack import sgesv
from matplotlib import pyplot as plt
from math import log, sqrt, pi, exp
from typing import  List

from do3se_met.physical_constants import VON_KAR as K
from do3se_met.model_constants import izR
from do3se_met.resistance.model import Resistance_Model
from do3se_met.resistance import *
from do3se_met.wind import ustar_from_velocity_simple, calc_layer_windspeed
from collections import namedtuple

OzoneConcentrationOutput = namedtuple('OzoneConcentrationOutput', 'O3_i micro_O3 Vd_i Vd')
# %% [markdown]
"""# Setup Constants"""
# %%

nL: int = 5
nLC: int = 1
u_at_canopy_top: float = 2.5
canopy_height: float = 24
# canopy_height: float = 1
SAI_values: List[List[float]] = [[0.1] for i in range(nL)]
LAI_values: List[List[float]] = [[0.1] for i in range(nL)]
mean_gsto_values: List[List[float]] = [[50] for i in range(nL)]
Rsoil: float = 200
Ts_C: float = 20
# O3_ppb_zR = 1
O3_ppb_zR = 30
top_layer_index = nL - 1
measured_height: float = 30
layer_heights = np.linspace(0, canopy_height, nL + 1)[1:]
print(layer_heights)
Lm = 0.1
SAI = np.sum(SAI_values)
W = [Lm * lai[0] for lai in LAI_values]
u_values = [calc_layer_windspeed(canopy_height, w, SAI, u_at_canopy_top, z) for z, w in zip(layer_heights, W)]
ustar_values = [ustar_from_velocity_simple(u, z=canopy_height, z0=0.1) for u in u_values]
# ustar_canopy_top: float = 0.8 # taken from the average of treescapes evergreen tree runs
ustar_canopy_top= ustar_values[top_layer_index]

# %%
print("ustar_values", ustar_values)
# %% [markdown]
"""First we setup an example using the existing setup"""
# %%

def calc_resistance_model(
    nL: int,
    nLC: int,
    ustar: float,
    canopy_height: float,
    SAI_values: List[List[float]],
    LAI_values: List[List[float]],
    mean_gsto_values: List[List[float]],
    Rsoil: float,
    Ts_C: float,
    snow_depth: float = None,
    Rext_base: float = 2500.0,
    Rb_diff: float = 0.000015,
    izr: float = izR,
    measured_height: float = 10,
    MIN_CANOPY_HEIGHT: float = 0.01,
    CANOPY_D: float = 0.7,
    CANOPY_Z0: float = 0.1,
) -> Resistance_Model:
    """This is a stripped back version of the function in resistance.py"""

    TOP_LAYER = -1
    LAI_sum_per_layer = [sum(LAI_values[iL]) for iL in range(nL)]
    SAI_sum_per_layer = [sum(SAI_values[iL]) for iL in range(nL)]
    mean_gsto_per_layer = [
        sum(
            [mean_gsto_values[iL][iLC] * LAI_values[iL][iLC] / LAI_sum_per_layer[iL] if LAI_sum_per_layer[iL] > 0 else 0
             for iLC in range(nLC)]
        )
        for iL in range(nL)
    ]

    bulk_gsto_per_layer = [
        sum(
            [LAI_values[iL][iLC] * mean_gsto_values[iL][iLC] * LAI_values[iL][iLC] / LAI_sum_per_layer[iL] if LAI_sum_per_layer[iL] > 0 else 0
             for iLC in range(nLC)]
        )
        for iL in range(nL)
    ]

    canopy_height_lim = max(MIN_CANOPY_HEIGHT, canopy_height)
    canopy_height_d = canopy_height_lim * CANOPY_D
    canopy_height_zo = canopy_height_lim * CANOPY_Z0

    Ra_canopy_to_izr = calc_Ra_simple(
        ustar,
        z1=canopy_height_zo,
        z2=izr,
        d=canopy_height_d,
    )

    Ra_canopy_top_to_izr = calc_Ra_simple(
        ustar,
        z1=canopy_height_lim - canopy_height_d,
        z2=izr,
        d=canopy_height_d,
    )

    Rb = calc_Rb(ustar, Rb_diff)

    Ra_measured_to_izr = calc_Ra_simple(
        ustar,
        z1=measured_height - canopy_height_d,
        z2=izr,
        d=canopy_height_d,
    )
    Rinc: List[float] = [calc_Rinc(sum(SAI_sum_per_layer), canopy_height_lim, ustar)
                         for iL in range(nL)]

    Rext: List[float] = [calc_Rext(Rext_base, SAI_sum_per_layer[iL]) for iL in range(nL)]

    Rsto: List[float] = [calc_Rsto(mean_gsto_per_layer[iL]) for iL in range(nL)]
    Rsto_c: List[float] = [calc_Rsto(bulk_gsto_per_layer[iL]) for iL in range(nL)]
    Rgs = calc_Rgs(Rsoil, snow_depth, Ts_C)

    # TODO: Should allow input LAI here

    Rsur: List = calc_Rsur_multilayer(
        nL,
        [Rb for _ in range(nL)],
        Rsto,
        Rext,
        LAI_sum_per_layer,
        SAI_sum_per_layer
    )
    Rsur_c: List = calc_Rsur_multilayer(
        nL,
        [Rb for _ in range(nL)],
        Rsto_c,
        Rext,
        LAI_sum_per_layer,
        SAI_sum_per_layer
    )


    # This is experimental and not used in standard runs.
    Rtotal: List[float] = calc_Rtotal(nL, Rsur, Rinc, Rgs)

    return Resistance_Model(
        nL,
        Ra_measured_to_izr=Ra_measured_to_izr,
        Ra_canopy_to_izr=Ra_canopy_to_izr,
        Ra_canopy_top_to_izr=Ra_canopy_top_to_izr,
        Rb=Rb,
        Rinc=Rinc,
        Rext=Rext,
        Rsto=Rsto,
        Rsto_c=Rsto_c,
        Rgs=Rgs,
        Rsur=Rsur,
        Rsur_c=Rsur_c,
        Rtotal=Rtotal,
    )

def calc_canopy_ozone_concentration(
    O3_ppb_zR: float,
    Ra_ref_canopy: float,
    Ra_ref_measured: float,
    Ra_tar_canopy: float,
    Ra_tar_canopy_top: float,
    Rsur_ref: float,
    Rsur_top_layer: float,
    Rb_ref: float,
    Rb_top_layer: float,
) -> OzoneConcentrationOutput:
    r"""This is a stripped back version of the function in deposition.py"""
    Vd_i = 1.0 / (Ra_ref_canopy + Rb_ref + Rsur_ref)
    O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_ref_measured * Vd_i))
    Vd = 1.0 / (Ra_tar_canopy + Rb_top_layer + Rsur_top_layer)
    O3_ppb = O3_ppb_i * (1.0 - (Ra_tar_canopy_top * Vd))

    return OzoneConcentrationOutput(
        O3_i=O3_ppb_i,
        micro_O3=O3_ppb,
        Vd_i=Vd_i,
        Vd=Vd,
    )

def calc_multi_layer_O3_ozone_concentration(
    nL: int,
    O3_in: float,
    rm_Ra: float,
    rm_Rinc: List[float],
    rm_Rsur: List[float],
    rm_Rgs: float,
) -> List[float]:
    """Stripped back version of the function in deposition.py"""
    C = np.full((nL + 1), 0, dtype=float)
    X = np.full((nL + 1, nL + 1), 0, dtype=float)

    bigR = np.array([rm_Ra] + rm_Rinc)
    assert bigR.shape == (nL + 1,)
    smallR = np.array(rm_Rsur + [rm_Rgs])
    assert smallR.shape == (nL + 1,)


    # for j in range(0, nL + 1):
    #     X[0:j, j] = bigR[0:j]
    #     X[j, j] = X[j, j] + smallR[j]
    #     if j < nL:  # TODO: Check this
    #         X[j + 1, j] = -smallR[j]
    # Nathans fix
    for j in range(0, nL + 1):
        X[0:j+1, j] = bigR[0:j+1]
        if j < nL:  # TODO: Check this
            X[j, j] = X[j, j] + smallR[j]
            X[j + 1, j] = -smallR[j]


    C[0] = O3_in
    lu, IPIV, C_out, info = sgesv(X, C)

    if info != 0:
        raise Exception('SGESV Failed')
    C_final = smallR * C_out
    O3_out = C_final[0: nL]
    return C_final
    # return O3_out


# %%

resistance_model = calc_resistance_model(
    nL=nL,
    nLC=nLC,
    ustar=ustar_canopy_top,
    canopy_height=canopy_height,
    SAI_values=SAI_values,
    LAI_values=LAI_values,
    mean_gsto_values=mean_gsto_values,
    Rsoil=Rsoil,
    Ts_C=Ts_C,
    measured_height=measured_height,
)
# Assume reference and modelled canopy are the same
resistance_model_ref = resistance_model



ozone_deposition = calc_canopy_ozone_concentration(
    O3_ppb_zR=O3_ppb_zR,
    Ra_ref_canopy=resistance_model.Ra_canopy_to_izr,
    Ra_ref_measured=resistance_model.Ra_measured_to_izr,
    Ra_tar_canopy=resistance_model.Ra_canopy_to_izr,
    Ra_tar_canopy_top=resistance_model.Ra_canopy_top_to_izr,
    Rsur_ref=resistance_model.Rsur[top_layer_index],
    Rsur_top_layer=resistance_model.Rsur[top_layer_index],
    Rb_ref=resistance_model.Rb,
    Rb_top_layer=resistance_model.Rb,
)

print(ozone_deposition._asdict())

multi_layer_ozone = calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=O3_ppb_zR,
    # rm_Ra=resistance_model.Ra_canopy_to_izr,
    rm_Ra = 0,
    rm_Rinc=resistance_model.Rinc,
    rm_Rsur=resistance_model.Rsur,
    rm_Rgs=resistance_model.Rgs,
)
plt.plot(multi_layer_ozone)
# %% [markdown]
"""Now implementing the new setup"""
# TODO: In updated version use bulk Rsur

# %%

def calc_Rtotal(
    nL: int,
    Rsur: List[float],
    Rinc: List[float],
    Rgs: float,
) -> List[float]:
    """Stripped back version of the function in resistance.py"""
    # TODO: Check list order is correct
    tmp = [None for _ in range(nL + 1)]
    tmp[0] = Rgs
    for i in range(1, nL + 1):
        tmp[i] = 1 / (1 / Rsur[i - 1] + 1 / (Rinc[i - 1] + tmp[i - 1]))
    Rtotal = tmp[1:nL + 1]
    return Rtotal



def calc_resistance_model_alt(
    nL: int,
    nLC: int,
    ustar_values: List[float],
    canopy_height: float,
    SAI_values: List[List[float]],
    LAI_values: List[List[float]],
    mean_gsto_values: List[List[float]],
    Rsoil: float,
    Ts_C: float,
    snow_depth: float = None,
    Rext_base: float = 2500.0,
    Rb_diff: float = 0.000015,
    izr: float = izR,
    measured_height: float = 10,
    MIN_CANOPY_HEIGHT: float = 0.01,
    CANOPY_D: float = 0.7,
    CANOPY_Z0: float = 0.1,
) -> Resistance_Model:
    """This is a stripped back version of the function in resistance.py"""

    TOP_LAYER = -1
    LAI_sum_per_layer = [sum(LAI_values[iL]) for iL in range(nL)]
    SAI_sum_per_layer = [sum(SAI_values[iL]) for iL in range(nL)]
    mean_gsto_per_layer = [
        sum(
            [mean_gsto_values[iL][iLC] * LAI_values[iL][iLC] / LAI_sum_per_layer[iL] if LAI_sum_per_layer[iL] > 0 else 0
             for iLC in range(nLC)]
        )
        for iL in range(nL)
    ]

    bulk_gsto_per_layer = [
        sum(
            [LAI_values[iL][iLC] * mean_gsto_values[iL][iLC] * LAI_values[iL][iLC] / LAI_sum_per_layer[iL] if LAI_sum_per_layer[iL] > 0 else 0
             for iLC in range(nLC)]
        )
        for iL in range(nL)
    ]
    print(MIN_CANOPY_HEIGHT, canopy_height)
    canopy_height_lim = max(MIN_CANOPY_HEIGHT, canopy_height)
    canopy_height_d = canopy_height_lim * CANOPY_D
    canopy_height_zo = canopy_height_lim * CANOPY_Z0

    Ra_canopy_to_izr = calc_Ra_simple(
        ustar_values[TOP_LAYER],
        z1=canopy_height_zo,
        z2=izr,
        d=canopy_height_d,
    )

    Ra_canopy_top_to_izr = calc_Ra_simple(
        ustar_values[TOP_LAYER],
        z1=canopy_height_lim - canopy_height_d,
        z2=izr,
        d=canopy_height_d,
    )

    Rb = [calc_Rb(ustar, Rb_diff) for ustar in ustar_values]

    Ra_measured_to_izr = calc_Ra_simple(
        ustar_values[TOP_LAYER],
        z1=measured_height - canopy_height_d,
        z2=izr,
        d=canopy_height_d,
    )
    Rinc: List[float] = [calc_Rinc(sum(SAI_sum_per_layer), canopy_height_lim, ustar_values[iL])
                         for iL in range(nL)]

    Rext: List[float] = [calc_Rext(Rext_base, SAI_sum_per_layer[iL]) for iL in range(nL)]

    Rsto: List[float] = [calc_Rsto(mean_gsto_per_layer[iL]) for iL in range(nL)]
    Rsto_c: List[float] = [calc_Rsto(bulk_gsto_per_layer[iL]) for iL in range(nL)]
    Rgs = calc_Rgs(Rsoil, snow_depth, Ts_C)

    # TODO: Check if we should be summing up layers below
    Rsur: List = calc_Rsur_multilayer(
        nL,
        Rb,
        Rsto,
        Rext,
        LAI_sum_per_layer,
        SAI_sum_per_layer
    )
    Rsur_c: List = calc_Rsur_multilayer(
        nL,
        Rb,
        Rsto_c,
        Rext,
        LAI_sum_per_layer,
        SAI_sum_per_layer
    )


    # This is experimental and not used in standard runs.
    # TODO: Check RTotal calcs
    # TODO: Check Rsur or Rsur bulk
    Rtotal: List[float] = calc_Rtotal(nL, Rsur_c, Rinc, Rgs)

    return Resistance_Model(
        nL,
        Ra_measured_to_izr=Ra_measured_to_izr,
        Ra_canopy_to_izr=Ra_canopy_to_izr,
        Ra_canopy_top_to_izr=Ra_canopy_top_to_izr,
        Rb=Rb,
        Rinc=Rinc,
        Rext=Rext,
        Rsto=Rsto,
        Rsto_c=Rsto_c,
        Rgs=Rgs,
        Rsur=Rsur,
        Rsur_c=Rsur_c,
        Rtotal=Rtotal,
    )

def calc_canopy_ozone_concentration_alt(
    O3_ppb_zR: float,
    # Ra_ref_canopy: float,
    Ra_ref_measured: float,
    Ra_tar_canopy: float,
    RTotal_ref: float,
    RTotal_tar: float,
) -> OzoneConcentrationOutput:
    r"""This is a stripped back version of the function in deposition.py"""
    Vd_i = 1.0 / (Ra_ref_measured + RTotal_ref)
    O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_ref_measured * Vd_i))
    Vd = 1.0 / (Ra_tar_canopy + RTotal_tar)
    O3_ppb = O3_ppb_i * (1.0 - (Ra_tar_canopy * Vd))

    return OzoneConcentrationOutput(
        O3_i=O3_ppb_i,
        micro_O3=O3_ppb,
        Vd_i=Vd_i,
        Vd=Vd,
    )

def calc_multi_layer_O3_ozone_concentration_alt(
    nL: int,
    O3_in: float,
    rm_Ra: float,
    rm_Rinc: List[float],
    rm_Rsur: List[float],
    rm_Rgs: float,
) -> List[float]:
    """Stripped back version of the function in deposition.py"""
    C = np.full((nL + 1), 0, dtype=float)
    X = np.full((nL + 1, nL + 1), 0, dtype=float)

    # bigR = np.array([rm_Ra] + rm_Rinc)
    # assert bigR.shape == (nL + 1,)
    # smallR = np.array(rm_Rsur + [rm_Rgs])
    # assert smallR.shape == (nL + 1,)
    bigR = np.array([rm_Ra] + rm_Rinc)
    assert bigR.shape == (nL + 1,)
    smallR = np.array(rm_Rsur + [rm_Rgs])
    assert smallR.shape == (nL + 1,)

    # for j in range(0, nL + 1):
    #     X[0:j, j] = bigR[0:j]
    #     X[j, j] = X[j, j] + smallR[j]
    #     if j < nL:  # TODO: Check this
    #         X[jI + 1, j] = -smallR[j]
    # Nathans fix
    for j in range(0, nL + 1):
        X[0:j+1, j] = bigR[0:j+1]
        if j < nL:  # TODO: Check this
            X[j, j] = X[j, j] + smallR[j]
            X[j + 1, j] = -smallR[j]
    C[0] = O3_in

    lu, IPIV, C_out, info = sgesv(X, C)

    if info != 0:
        raise Exception('SGESV Failed')
    C_final = smallR * C_out
    O3_out = C_final[0: nL]
    print(len(O3_out), len(C_final))
    # return O3_out
    return C_final

# %%

print("u_values", u_values)
print("total sai", SAI)
print("Lm_LAI", W)
resistance_model_alt = calc_resistance_model_alt(
    nL=nL,
    nLC=nLC,
    ustar_values=ustar_values,
    canopy_height=canopy_height,
    SAI_values=SAI_values,
    LAI_values=LAI_values,
    mean_gsto_values=mean_gsto_values,
    Rsoil=Rsoil,
    Ts_C=Ts_C,
    measured_height=measured_height,
)
# %%

# Assume reference and modelled canopy are the same
# resistance_model_ref_alt = resistance_model_alt

ozone_deposition = calc_canopy_ozone_concentration_alt(
    O3_ppb_zR=O3_ppb_zR,
    # TODO: CHeck the difference between measured and canopy
    # Ra_ref_canopy=resistance_model_alt.Ra_canopy_to_izr,
    Ra_ref_measured=resistance_model_alt.Ra_measured_to_izr,
    # TODO: Check diff between top layer and canopy
    Ra_tar_canopy=resistance_model_alt.Ra_canopy_to_izr,
    # Ra_tar_canopy_top=resistance_model_alt.Ra_canopy_top_to_izr,
    # Rsur_ref=resistance_model_alt.Rsur[top_layer_index],
    # Rsur_top_layer=resistance_model_alt.Rsur[top_layer_index],
    # Rb_ref=resistance_model_alt.Rb,
    # Rb_top_layer=resistance_model_alt.Rb,
    RTotal_ref=resistance_model_alt.Rtotal[top_layer_index],
    RTotal_tar=resistance_model_alt.Rtotal[top_layer_index],
)

print(ozone_deposition._asdict())

multi_layer_ozone_alt = calc_multi_layer_O3_ozone_concentration_alt(
    nL=nL,
    O3_in=O3_ppb_zR,
    rm_Ra=resistance_model_alt.Ra_canopy_to_izr,
    rm_Rinc=list(reversed(resistance_model_alt.Rinc)),
    rm_Rsur=list(reversed(resistance_model_alt.Rsur_c)),
    rm_Rgs=resistance_model_alt.Rgs,
)
plt.plot(multi_layer_ozone_alt, label="new multi layer")
plt.plot(multi_layer_ozone, label="pyDO3SE")
plt.title("Ozone")
plt.xlabel("Layer")
plt.ylabel("Ozone concentration")
print(multi_layer_ozone, multi_layer_ozone_alt)
plt.legend()
# %%

# Assume reference and modelled canopy are the same
# resistance_model_ref_alt = resistance_model_alt

ozone_deposition = calc_canopy_ozone_concentration_alt(
    O3_ppb_zR=O3_ppb_zR,
    # TODO: CHeck the difference between measured and canopy
    # Ra_ref_canopy=resistance_model_alt.Ra_canopy_to_izr,
    Ra_ref_measured=resistance_model_alt.Ra_measured_to_izr,
    # TODO: Check diff between top layer and canopy
    Ra_tar_canopy=resistance_model_alt.Ra_canopy_to_izr,
    # Ra_tar_canopy_top=resistance_model_alt.Ra_canopy_top_to_izr,
    # Rsur_ref=resistance_model_alt.Rsur[top_layer_index],
    # Rsur_top_layer=resistance_model_alt.Rsur[top_layer_index],
    # Rb_ref=resistance_model_alt.Rb,
    # Rb_top_layer=resistance_model_alt.Rb,
    RTotal_ref=resistance_model_alt.Rtotal[top_layer_index],
    RTotal_tar=resistance_model_alt.Rtotal[top_layer_index],
)

print(ozone_deposition._asdict())

multi_layer_ozone_alt = calc_multi_layer_O3_ozone_concentration_alt(
    nL=nL,
    O3_in=ozone_deposition.O3_i,
    rm_Ra=resistance_model_alt.Ra_canopy_to_izr,
    rm_Rinc=list(reversed(resistance_model_alt.Rinc)),
    rm_Rsur=list(reversed(resistance_model_alt.Rsur_c)),
    rm_Rgs=resistance_model_alt.Rgs,
)
plt.scatter(list(reversed([0, *layer_heights])), multi_layer_ozone_alt, label="multilayer model")
plt.scatter(list(reversed([0, *layer_heights])), multi_layer_ozone, label="pyDO3SE")
print(multi_layer_ozone, multi_layer_ozone_alt)
plt.title("Ozone")
plt.xlabel("Layer height")
plt.ylabel("Ozone concentration")
plt.legend()
# %%
layer_heights

# """Add new version with per layer Rb"""

# def calc_Rtotal_new(
#     nL: int,
#     Rsur: List[float],
#     Rinc: List[float],
#     Rgs: float,
# ) -> List[float]:
#     """Stripped back version of the function in resistance.py"""
#     # TODO: Check list order is correct
#     tmp = [None for _ in range(nL + 1)]
#     tmp[0] = Rgs
#     for i in range(1, nL + 1):
#         tmp[i] = 1 / (1 / Rsur[i - 1] + 1 / (Rinc[i - 1] + tmp[i - 1]))
#     Rtotal = tmp[1:nL + 1]
#     return Rtotal


# def calc_resistance_model_new(
#     nL: int,
#     nLC: int,
#     ustar: List[float],
#     canopy_height: float,
#     SAI_values: List[List[float]],
#     LAI_values: List[List[float]],
#     mean_gsto_values: List[List[float]],
#     Rsoil: float,
#     Ts_C: float,
#     snow_depth: float = None,
#     Rext_base: float = 2500.0,
#     Rb_diff: float = 0.000015,
#     izr: float = izR,
#     measured_height: float = 10,
#     MIN_CANOPY_HEIGHT: float = 0.01,
#     CANOPY_D: float = 0.7,
#     CANOPY_Z0: float = 0.1,
# ) -> Resistance_Model:
#     """This is a stripped back version of the function in resistance.py"""

#     TOP_LAYER = -1
#     LAI_sum_per_layer = [sum(LAI_values[iL]) for iL in range(nL)]
#     SAI_sum_per_layer = [sum(SAI_values[iL]) for iL in range(nL)]
#     mean_gsto_per_layer = [
#         sum(
#             [mean_gsto_values[iL][iLC] * LAI_values[iL][iLC] / LAI_sum_per_layer[iL] if LAI_sum_per_layer[iL] > 0 else 0
#              for iLC in range(nLC)]
#         )
#         for iL in range(nL)
#     ]

#     bulk_gsto_per_layer = [
#         sum(
#             [LAI_values[iL][iLC] * mean_gsto_values[iL][iLC] * LAI_values[iL][iLC] / LAI_sum_per_layer[iL] if LAI_sum_per_layer[iL] > 0 else 0
#              for iLC in range(nLC)]
#         )
#         for iL in range(nL)
#     ]
#     print(MIN_CANOPY_HEIGHT, canopy_height)
#     canopy_height_lim = max(MIN_CANOPY_HEIGHT, canopy_height)
#     canopy_height_d = canopy_height_lim * CANOPY_D
#     canopy_height_zo = canopy_height_lim * CANOPY_Z0

#     Ra_canopy_to_izr = calc_Ra_simple(
#         ustar[top_layer_index],
#         z1=canopy_height_zo,
#         z2=izr,
#         d=canopy_height_d,
#     )

#     Ra_canopy_top_to_izr = calc_Ra_simple(
#         ustar[top_layer_index],
#         z1=canopy_height_lim - canopy_height_d,
#         z2=izr,
#         d=canopy_height_d,
#     )

#     # TODO: Multilayer Rb
#     Rb = [calc_Rb(us, Rb_diff) for us in ustar]

#     Ra_measured_to_izr = calc_Ra_simple(
#         ustar[top_layer_index],
#         z1=measured_height - canopy_height_d,
#         z2=izr,
#         d=canopy_height_d,
#     )
#     Rinc: List[float] = [calc_Rinc(sum(SAI_sum_per_layer), canopy_height_lim, ustar)
#                          for iL in range(nL)]

#     Rext: List[float] = [calc_Rext(Rext_base, SAI_sum_per_layer[iL]) for iL in range(nL)]

#     Rsto: List[float] = [calc_Rsto(mean_gsto_per_layer[iL]) for iL in range(nL)]
#     Rsto_c: List[float] = [calc_Rsto(bulk_gsto_per_layer[iL]) for iL in range(nL)]
#     Rgs = calc_Rgs(Rsoil, snow_depth, Ts_C)

#     # TODO: Check if we should be summing up layers below
#     Rsur: List = calc_Rsur_multilayer(
#         nL,
#         Rb,
#         Rsto,
#         Rext,
#         LAI_sum_per_layer,
#         SAI_sum_per_layer
#     )
#     Rsur_c: List = calc_Rsur_multilayer(
#         nL,
#         Rb,
#         Rsto_c,
#         Rext,
#         LAI_sum_per_layer,
#         SAI_sum_per_layer
#     )


#     # This is experimental and not used in standard runs.
#     # TODO: Check RTotal calcs
#     # TODO: Check Rsur or Rsur bulk
#     Rtotal: List[float] = calc_Rtotal_new(nL, Rsur_c, Rinc, Rgs)

#     return Resistance_Model(
#         nL,
#         Ra_measured_to_izr=Ra_measured_to_izr,
#         Ra_canopy_to_izr=Ra_canopy_to_izr,
#         Ra_canopy_top_to_izr=Ra_canopy_top_to_izr,
#         Rb=Rb,
#         Rinc=Rinc,
#         Rext=Rext,
#         Rsto=Rsto,
#         Rsto_c=Rsto_c,
#         Rgs=Rgs,
#         Rsur=Rsur,
#         Rsur_c=Rsur_c,
#         Rtotal=Rtotal,
#     )

# def calc_canopy_ozone_concentration_new(
#     O3_ppb_zR: float,
#     # Ra_ref_canopy: float,
#     Ra_ref_measured: float,
#     Ra_tar_canopy: float,
#     RTotal_ref: float,
#     RTotal_tar: float,
# ) -> OzoneConcentrationOutput:
#     r"""This is a stripped back version of the function in deposition.py"""
#     Vd_i = 1.0 / (Ra_ref_measured + RTotal_ref)
#     O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_ref_measured * Vd_i))
#     Vd = 1.0 / (Ra_tar_canopy + RTotal_tar)
#     O3_ppb = O3_ppb_i * (1.0 - (Ra_tar_canopy * Vd))

#     return OzoneConcentrationOutput(
#         O3_i=O3_ppb_i,
#         micro_O3=O3_ppb,
#         Vd_i=Vd_i,
#         Vd=Vd,
#     )

# def calc_multi_layer_O3_ozone_concentration_new(
#     nL: int,
#     O3_in: float,
#     rm_Ra: float,
#     rm_Rinc: List[float],
#     rm_Rsur: List[float],
#     rm_Rgs: float,
# ) -> List[float]:
#     """Stripped back version of the function in deposition.py"""
#     C = np.full((nL + 1), 0, dtype=float)
#     X = np.full((nL + 1, nL + 1), 0, dtype=float)

#     bigR = np.array([rm_Ra] + rm_Rinc)
#     assert bigR.shape == (nL + 1,)
#     smallR = np.array(rm_Rsur + [rm_Rgs])
#     assert smallR.shape == (nL + 1,)


#     # for j in range(0, nL + 1):
#     #     X[0:j, j] = bigR[0:j]
#     #     X[j, j] = X[j, j] + smallR[j]
#     #     if j < nL:  # TODO: Check this
#     #         X[j + 1, j] = -smallR[j]
#     # Nathans fix
#     for j in range(0, nL + 1):
#         X[0:j+1, j] = bigR[0:j+1]
#         if j < nL:  # TODO: Check this
#             X[j, j] = X[j, j] + smallR[j]
#             X[j + 1, j] = -smallR[j]
#     C[0] = O3_in

#     lu, IPIV, C_out, info = sgesv(X, C)

#     if info != 0:
#         raise Exception('SGESV Failed')
#     C_final = smallR * C_out
#     O3_out = C_final[0: nL]
#     return O3_out

# # %%

# resistance_model_new = calc_resistance_model_new(
#     nL=nL,
#     nLC=nLC,
#     ustar=ustar,
#     canopy_height=canopy_height,
#     SAI_values=SAI_values,
#     LAI_values=LAI_values,
#     mean_gsto_values=mean_gsto_values,
#     Rsoil=Rsoil,
#     Ts_C=Ts_C,
#     measured_height=measured_height,
# )
# # %%

# # Assume reference and modelled canopy are the same
# # resistance_model_ref_new = resistance_model_new

# ozone_deposition = calc_canopy_ozone_concentration_new(
#     O3_ppb_zR=O3_ppb_zR,
#     # TODO: CHeck the difference between measured and canopy
#     # Ra_ref_canopy=resistance_model_new.Ra_canopy_to_izr,
#     Ra_ref_measured=resistance_model_new.Ra_measured_to_izr,
#     # TODO: Check diff between top layer and canopy
#     Ra_tar_canopy=resistance_model_new.Ra_canopy_to_izr,
#     # Ra_tar_canopy_top=resistance_model_new.Ra_canopy_top_to_izr,
#     # Rsur_ref=resistance_model_new.Rsur[top_layer_index],
#     # Rsur_top_layer=resistance_model_new.Rsur[top_layer_index],
#     # Rb_ref=resistance_model_new.Rb,
#     # Rb_top_layer=resistance_model_new.Rb,
#     RTotal_ref=resistance_model_new.Rtotal[top_layer_index],
#     RTotal_tar=resistance_model_new.Rtotal[top_layer_index],
# )

# print(ozone_deposition._asdict())

# multi_layer_ozone_new = calc_multi_layer_O3_ozone_concentration_new(
#     nL=nL,
#     O3_in=O3_ppb_zR,
#     rm_Ra=resistance_model_alt.Ra_canopy_to_izr,
#     rm_Rinc=resistance_model_alt.Rinc,
#     rm_Rsur=resistance_model_alt.Rsur_c,
#     rm_Rgs=resistance_model_alt.Rgs,
# )
# plt.plot(multi_layer_ozone_new, label="New per layer Rb")
# plt.plot(multi_layer_ozone_alt, label="DO3SE model")
# plt.plot(multi_layer_ozone, label="pyDO3SE")
# plt.legend()
# # %%

# %%

resistance_model_alt.Rsur
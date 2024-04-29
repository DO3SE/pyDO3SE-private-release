# %%
"""Helper functinos associated with the met module."""

from math import exp, isclose
from pyDO3SE.constants.physical_constants import DIFF_O3
from typing import List, NamedTuple, Tuple
from collections import namedtuple
import numpy as np
from scipy.linalg.lapack import sgesv

from pyDO3SE.plugins.met.wind import ustar_from_velocity
from pyDO3SE.plugins.resistance.helpers import calc_Ra_simple, calc_Rb, calc_deposition_velocity
from pyDO3SE.constants.model_constants import CANOPY_D, CANOPY_Z0


estimate_ustar = ustar_from_velocity
ra_simple = calc_Ra_simple
rb_func = calc_Rb
# %%


def calc_O3(
    O3: float,
    canopy_height: float,
    u_50: float,
    h_O3_in: float,
    z_O3: float,
    Rsur_top_layer: float,
    Rb_top_layer: float,
    Ra_top_layer: float,
    ustar: float,
    # Rtotal_top_layer: float,
):
    # use Constants, only: k, izR, v, DO3, Pr, Ts_K
    # use Inputs, only: O3_ppb_zR, uh_i, Ts_C, P, ustar
    # use Inputs, only: estimate_ustar
    # use Variables, only: Ra, Rb, Rsur
    # use Variables, only: Vd, O3_ppb, O3_nmol_m3, Vd_i, O3_ppb_i, Ra_ref_i, &
    #                         Ra_ref, Ra_O3zR_i, Ra_tar_i
    # use Parameters, only: O3zR, O3_d, O3_zo, h, d, zo
    # use R, only: ra_simple, rb_func => rb

    uh_i = u_50
    izR = 50
    Rsur = Rsur_top_layer
    Rb = Rb_top_layer
    Ra = Ra_top_layer
    h_O3 = canopy_height

    O3_d = h_O3 * CANOPY_D
    O3_z0 = h_O3 * CANOPY_Z0
    h = canopy_height
    z0 = h * CANOPY_Z0
    zo = z0
    O3_zo = O3_z0

    DO3 = DIFF_O3
    O3zR = z_O3
    O3_ppb_zR = O3
    d = h * CANOPY_D

    # FORTRAN COPY ==========

    M_O3 = 48.0  # Molecular weight of O3 (g)

    # real :: ustar_ref, Rb_ref, Vn

    # ustar over reference canopy
    ustar_ref = estimate_ustar(uh_i, izR - O3_d, O3_zo)
    # Ra between reference canopy and izR
    Ra_ref_i = ra_simple(ustar_ref, O3_zo + O3_d, izR, O3_d)
    # Rb for reference canopy
    Rb_ref = rb_func(ustar_ref, DO3)
    # Deposition velocity at izR over reference canopy
    # (assuming that Rsur_ref = Rsur)
    Vd_i = 1.0 / (Ra_ref_i + Rb_ref + Rsur)
    # Ra between measurement height and izR
    Ra_O3zR_i = ra_simple(ustar_ref, O3zR, izR, O3_d)
    # O3 concentration at izR
    O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_O3zR_i * Vd_i))
    # Ra between target canopy and izR
    # (ustar already calculated for target canopy)
    Ra_tar_i = ra_simple(ustar, zo + d, izR, d)
    # Deposition velocity at izR over target canopy
    Vd = 1.0 / (Ra_tar_i + Rb + Rsur)
    # O3 concentration at target canopy
    # (Ra already calculated between canopy height and izR)
    O3_ppb = O3_ppb_i * (1.0 - (Ra * Vd))

    # Specific molar volume of an ideal gas at current temp + pressure
    # Vn = 8.314510 * ((Ts_C + Ts_K) / P)
    # # Convert to nmol/m^3
    # O3_nmol_m3 = (1.0 / Vn) * O3_ppb * M_O3 * 20.833  # 1 microgram O3 = 20.833 nmol/m^3

    return O3_ppb, O3_ppb_i


O3_ppb, O3_50 = calc_O3(
    O3=33.904,
    canopy_height=1.0,
    u_50=2.93485212326,
    h_O3_in=None,
    z_O3=50.0,
    Rsur_top_layer=276.901275635,  # TODO: should come from args
    Rb_top_layer=31.2906856537,  # TODO: should come from args
    Ra_top_layer=64.121711731,  # TODO: should come from args
    ustar=0.19406299293,
)
print(O3_50, O3_ppb)

assert O3_50 == 33.904
assert isclose(O3_ppb, 28.27369, abs_tol=1e-5)

# %%
# 4208
O3_ppb, O3_50 = calc_O3(
    O3=35.1069984436,
    canopy_height=1.0,
    u_50=1.35454714298,
    h_O3_in=None,
    z_O3=50.0,
    Rsur_top_layer=62.4227142334,  # TODO: should come from args
    Rb_top_layer=67.7964859009,  # TODO: should come from args
    Ra_top_layer=138.930374146,  # TODO: should come from args
    ustar=0.0895675346255,
)
print(O3_50, O3_ppb)

assert O3_50 == 35.1069984436
assert isclose(O3_ppb, 18.7981281281, abs_tol=1e-5)

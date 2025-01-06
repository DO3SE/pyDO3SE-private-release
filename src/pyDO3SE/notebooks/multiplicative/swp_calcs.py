# %%
# CHECK SWP IS WORKING
from math import exp
from collections import namedtuple
from pyDO3SE.constants.physical_constants import DRATIO
# PYTHON

# %%


def SWC_to_SWP(SWP_AE, b, SWC) -> float:
    SWC_sat = 0.4  # Saturated soil water content for soil water release curve

    SWP = SWP_AE * ((SWC_sat / SWC)**b)
    return SWP


def SWP_to_SWC(SWP_AE, b, SWP):
    SWC_sat = 0.4  # Saturated soil water content for soil water release curve

    SWC = 1.0 / (((SWP / SWP_AE)**(1.0 / b)) / SWC_sat)
    return SWC


# %%
PWP = -4.0
root_depth = 0.75
# soil config
b = 6.58
FC = 0.29
SWP_AE = -0.00188
Ksat = 0.0002286
SWC_sat = 0.4

# SWP_AE=None
# b=None
# FC=None
# FC=None

# %%
# soil_moisture_from_SWP

T0 = 273.15


def get_swp(
    ASW,

    precip_acc,
    VPD_kPa,
    Ts_C,
    P_kPa,
    Rn_MJ,
    esat_kPa,
    eact_kPa,
    rm_Rb,
    LAI,
    rm_Rsto_l0,
    Sn_in_prev,
    rm_Rinc_l0,
    Es_blocked,

    rm_Rgs,
    Ei_acc,
    Et_acc,
    Es_acc,
    AEt_acc,
):
    # === penman_monteith_hourly

    Rb_H2O = rm_Rb
    Rinc = rm_Rinc_l0
    Rsto_H2O = rm_Rsto_l0
    Rsoil = rm_Rgs

    VPD = VPD_kPa * 1000

    Rn = Rn_MJ * 1000000
    P = P_kPa * 1000
    esat = esat_kPa * 1000
    eact = eact_kPa * 1000

    # ============ CHECKED TO HERE
    Tvir = (Ts_C + T0) / (1 - (0.378 * (eact / P)))
    delta = ((4098 * esat) / ((Ts_C + 237.3)**2))
    lmd = (2501000 - (2361 * Ts_C))
    psychro = 1628.6 * (P / lmd)
    Pair = (0.003486 * (P / Tvir))
    Cair = (0.622 * ((lmd * psychro) / P))

    G = 0.1 * Rn

    Et_1 = (delta * (Rn - G)) / lmd
    Et_2 = 3600 * Pair * Cair * VPD / Rb_H2O / lmd

    Ei_3 = delta + psychro
    Ei_hr = (Et_1 + Et_2) / Ei_3 / 1000

    Et_3 = delta + psychro * (1 + Rsto_H2O / Rb_H2O)
    Et_hr = (Et_1 + Et_2) / Et_3 / 1000

    if Es_blocked:
        Es_hr = 0
    else:
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lmd
        Es_2 = ((3600 * Pair * Cair * VPD) - (delta * Rinc *  # noqa: W504
                                              ((Rn - G) - (Es_Rn - Es_G)))) / (Rinc + Rb_H2O) / lmd
        Es_3 = delta + (psychro * (1.0 + (Rsoil / (Rb_H2O + Rinc))))
        Es_hr = (Es_1 + Es_2) / Es_3 / 1000

    # Calculate Eat from Et and Es (after Shuttleworth and Wallace, 1985)
    SW_a = (delta + psychro) * Rb_H2O
    SW_s = (delta + psychro) * Rinc + (psychro * Rsoil)
    SW_c = psychro * Rsto_H2O  # Boundary layer resistance = 0
    C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
    C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
    if Es_hr <= 0:
        AEt_hr = Et_hr
    else:
        AEt_hr = (C_canopy * Et_hr) + (C_soil * Es_hr)

    # Accumulate values # OK
    Ei_acc = Ei_acc + Ei_hr
    Et_acc = Et_acc + Et_hr
    Es_acc = Es_acc + Es_hr
    AEt_acc = AEt_acc + AEt_hr

    # ===
    # =====Calc Penmonteith Daily
    # THIS PART in PEN calc daily

    intercepted_evaporated = min(precip_acc, 0.0001 * LAI, Ei_acc)  # ok

    # TODO: In ui we do not limit evapotranspiration
    max_ET = ASW + precip_acc - intercepted_evaporated
    evapotranspiration = min(max_ET, AEt_acc)

    P_input = precip_acc - intercepted_evaporated  # OK

    delta_SM = P_input - evapotranspiration  # OK

    Sn_diff = delta_SM / root_depth  # ok

    # =====Calc SWP soil_moisture_from_SWC
    PWP_vol = 1.0 / (((PWP / SWP_AE)**(1.0 / b)) / SWC_sat)  # OK

    # Sn_in = 1.0 / (((SWP / SWP_AE)**(1.0 / b)) / SWC_sat)  # ?? Check if done elsewhere in ui

    # Constrain soil water content to be between field capacity and PWP
    Sn_in = Sn_in_prev + Sn_diff  # from input args
    Sn = max(PWP_vol, min(FC, Sn_in))  # OK

    # Calculate soil water potential (SWP)
    SWP = SWP_AE * ((SWC_sat / Sn)**b)  # OK

    # Calculate available soil water (ASW)
    ASW = (Sn - PWP_vol) * root_depth  # OK

    # Calculate soil moisture deficit (SMD)
    SMD = (FC - Sn) * root_depth  # OK

    return SWP


get_swp(-0.01)

# %%

# FORTRAN


def calc_SWP(
    precip_acc,
    VPD,
    Ts_C,
    P,
    Rn_MJ,
    esat_kPa,
    eact_kPa,
    Rb_H2O,
    LAI,
    Rsto_c,
    Rsto_PEt,
    Sn_in,
    Rinc,
    Es_blocked,
    Rsoil,
    Ei_acc=0,
    PEt_acc=0,
    Et_acc=0,
    Es_acc=0,
    AEt_acc=0
):
    Rsto_H2O = Rsto_c * DRATIO
    # Calc_Penman_Monteith (HOURLY)
    Ts_K = 273.15

    VPD_Pa = VPD * 1000

    Rn = Rn_MJ * 1000000.0
    P_Pa = P * 1000
    esat = 1000 * esat_kPa
    eact = 1000 * eact_kPa

    # ============ CHECKED TO HERE
    Tvir = (Ts_C + Ts_K) / (1 - (0.378 * (eact / P_Pa)))
    delta = ((4098 * esat) / ((Ts_C + 237.3)**2))
    lamb = (2501000 - (2361 * Ts_C))
    psychro = 1628.6 * (P_Pa / lamb)
    Pair = (0.003486 * (P_Pa / Tvir))
    Cair = (0.622 * ((lamb * psychro) / P_Pa))

    G = 0.1 * Rn

    Et_1 = (delta * (Rn - G)) / lamb
    Et_2 = 3600 * Pair * Cair * VPD_Pa / Rb_H2O / lamb

    Ei_3 = delta + psychro
    Ei_hr = (Et_1 + Et_2) / Ei_3 / 1000

    # PEt not in pyDO3SE
    PEt_3 = delta + psychro * (1 + (Rsto_PEt * DRATIO) / Rb_H2O)
    PEt_hr = (Et_1 + Et_2) / PEt_3 / 1000

    Et_3 = delta + psychro * (1 + Rsto_H2O / Rb_H2O)
    Et_hr = (Et_1 + Et_2) / Et_3 / 1000

    if (Es_blocked):
        Es_hr = 0
    else:
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lamb
        Es_2 = ((3600 * Pair * Cair * VPD_Pa) - (delta * Rinc *
                                                 ((Rn - G) - (Es_Rn - Es_G)))) / (Rinc + Rb_H2O) / lamb
        Es_3 = delta + (psychro * (1.0 + (Rsoil / (Rb_H2O + Rinc))))
        Es_hr = (Es_1 + Es_2) / Es_3 / 1000

    SW_a = (delta + psychro) * Rb_H2O
    SW_s = (delta + psychro) * Rinc + (psychro * Rsoil)
    SW_c = psychro * Rsto_H2O  # Boundary layer resistance = 0
    C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
    C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
    if (Es_hr <= 0):
        AEt_hr = Et_hr
    else:
        AEt_hr = (C_canopy * Et_hr) + (C_soil * Es_hr)

    Ei_acc = Ei_acc + Ei_hr
    PEt_acc = PEt_acc + PEt_hr
    Et_acc = Et_acc + Et_hr
    Es_acc = Es_acc + Es_hr
    AEt_acc = AEt_acc + AEt_hr

    # ==================== Calc swp
    evapotranspiration = AEt = AEt_acc
    # TODO: Where is PWP_vol set
    PWP_vol = 1.0 / (((PWP / SWP_AE)**(1.0 / b)) / SWC_sat)

    if (precip_acc > 0):
        intercepted_evaporated = min(Ei, 0.0001 * LAI)
        P_input = precip_acc - intercepted_evaporated
    else:
        P_input = 0
    # Can't lose water through Ei
    P_input = max(0.0, P_input)

    delta_SM = (P_input - evapotranspiration)
    Sn_diff = delta_SM / root_depth

    # Calculate new Sn, with field capacity as a maximum
    Sn = max(PWP_vol, min(FC, Sn_in + Sn_diff))
    # per_vol = Sn * 100

    # Calculate ASW and SWP for new water content
    SWP = SWP_AE * ((SWC_sat / Sn)**b)
    ASW = (Sn - PWP_vol) * root_depth

    # Calculate SMD for new water content
    SMD = (FC - Sn) * root_depth

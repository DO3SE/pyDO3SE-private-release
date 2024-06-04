"""O3 Helpers.

Functions that calculate O3 state.

Taken from the DO3SE-UI model
"""

from math import isclose
from do3se_met.conversion import deg_to_kel
from typing import List, NamedTuple, Tuple
from collections import namedtuple

from do3se_met.resistance import (
    calc_Ra_simple,
    calc_Ra_with_heat_flux,
    calc_Rb,
)

from pyDO3SE.constants.model_constants import izR

from pyDO3SE.plugins.resistance.helpers import (
    calc_Rext,
    calc_Rgs,
    calc_Rinc,
    calc_Rsto,
    calc_Rsur, calc_Rsur_multilayer,
    calc_Rtotal,
    calc_leaf_gb,
    calc_leaf_rb,
)
from pyDO3SE.plugins.resistance.model import Leaf_Resistance_Model, Resistance_Model

from pyDO3SE.constants.model_constants import CANOPY_D, CANOPY_Z0, DT, LEAF_G_O3
from pyDO3SE.constants.physical_constants import DIFF_O3


def calc_resistance_model(
    nL: int,
    nLC: int,
    ustar: float,
    canopy_height: float,
    SAI_values: List[List[float]],
    LAI_values: List[List[float]],
    mean_gsto_values: List[List[float]],
    Rsoil: float,
    L: float = None,
    ra_calc_method: str = "simple",
    rsur_calc_method: str = "single_layer",
    rext_calc_method: str = "const",
    Ts_C: float = None,
    snow_depth: float = None,
    Rext_base: float = 2500.0,
    izr: float = izR
) -> Resistance_Model:
    """Calculate the resistance model for O3 over the target canopy.

    Parameters
    ----------
    nL: int
        Number of model layers
    nLC: int
        Number of model components
    ustar: float
        Friction velocity over target canopy [m/s]
    canopy_height: float
        The total height of the canopy [m]
    SAI_values: List[float]
        full model SAI values with shape(nL, nLC)
    LAI_values: List[float]
        full model LAI values with shape(nL, nLC)
    mean_gsto_values: float
        full model mean_gsto values with shape(nL, nLC) [mmol O3 m-2 PLA s-1]
    Rsoil: float
        Soil resistance?
    L: float
        Monin-Obukhov length [m]
    ra_calc_method: str
        Method of calculating Ra. Options:
        - "simple" -
        - "heat_flux" -
    rsur_calc_method: str
        Method of calculating Rsur. Options:
        - "single_layer" -
        - "multi_layer" -
    rext_calc_method: str
        Method of calculating Rext. Options:
        - "const" - use Rext_const value
        - "calculate_SAI" - calculate using Rext base * SAI
    Ts_C: float
        Air Temperature [Degrees]
    Hd: float
        Sensible heat flux [W/m^2] only needed for Ra heat flux calc
    P: float
        Air Pressure [kPa] only needed for Ra heat flux calc
    snow_depth: float
        Measured snow depth. Default None
    Rext_const: float
        Rext base value
    izr: float
        Decoupled height [m]


    Returns
    -------
    resistance_model: Resistance_Model
        A new instance of the Resistance Model

    """
    if nL > 1 and rsur_calc_method == "single_layer":
        raise ValueError("Cannot use single layer rsur calc for multilayer model")

    TOP_LAYER = -1
    LAI_sum_per_layer = [sum(LAI_values[iL]) for iL in range(nL)]
    SAI_sum_per_layer = [sum(SAI_values[iL]) for iL in range(nL)]
    bulk_gsto_sum_per_layer = [sum([mean_gsto_values[iL][iLC]*LAI_values[iL][iLC] for iLC in range(nLC)]) for iL in range(nL)]
    # TODO: Check that mean gsto correct here
    mean_gsto_per_layer = [
        sum([mean_gsto_values[iL][iLC]*LAI_values[iL][iLC]/LAI_sum_per_layer[iL] if LAI_sum_per_layer[iL]>0 else 0 for iLC in range(nLC)]) for iL in range(nL)]

    # TODO: Check difference between Ra and Ra_c
    # TODO: Should we use different Ra calc for Ra_c
    # Ra_c uses (canopy_height * (CANOPY_D + CANOPY_Z0) instead of canopy_height
    # for z0 in simple calc

    Ra_c = calc_Ra_simple(ustar, (canopy_height * (CANOPY_D + CANOPY_Z0)),
                          izr, canopy_height * CANOPY_D)

    # TODO: Explain which Ra this is!
    Ra = calc_Ra_with_heat_flux(
        ustar,
        z1=canopy_height,
        z2=izr,
        L=L,
    ) if ra_calc_method == "heat_flux" else calc_Ra_simple(
        ustar,
        z1=canopy_height,
        z2=izr,
        d=canopy_height * CANOPY_D,
    )

    Rb = calc_Rb(ustar, DIFF_O3)

    # TODO: Rinc should use layer depth instead of canopy height
    # TODO: to use in a multilayer model, what does h represent?  Height above
    # ground, or thickness of layer?
    # Rinc: List[float] = [calc_Rinc(SAI_sum_per_layer[iL], layer depth, ustar) for iL in range(nL)]
    Rinc: List[float] = [calc_Rinc(sum(SAI_sum_per_layer), canopy_height, ustar)
                         for iL in range(nL)]

    if rext_calc_method == "const":
        Rext: List[float] = [Rext_base for _ in range(nL)]
    elif rext_calc_method == "calculate_SAI":
        Rext: List[float] = [calc_Rext(SAI_sum_per_layer[iL]) for iL in range(nL)]
    else:
        raise ValueError("Invalid rext calc method")

    Rsto: List[float] = [calc_Rsto(mean_gsto_per_layer[iL]) for iL in range(nL)]
    Rsto_c: List[float] = [calc_Rsto(bulk_gsto_sum_per_layer[iL]) for iL in range(nL)]
    Rgs = calc_Rgs(Rsoil, snow_depth, Ts_C)

    # Rsur not calculated for layers
    # if rsur_calc_method == "single_layer":
    #     Rsur: list = [calc_Rsur(
    #         Rb,
    #         Rsto[TOP_LAYER],
    #         Rext[TOP_LAYER],
    #         Rinc[TOP_LAYER],
    #         Rgs,
    #         LAI_sum_per_layer[TOP_LAYER],
    #         SAI_sum_per_layer[TOP_LAYER],
    #     ) for _ in range(nL)]
    # elif rsur_calc_method == "multi_layer":
    #     Rsur: List = calc_Rsur_multilayer(
    #         nL,
    #         Rb,
    #         Rsto,
    #         Rext,
    #         LAI_sum_per_layer,
    #         SAI_sum_per_layer
    #     )
    # else:
    #     raise ValueError("Invalid Rsur calc method")

    # TODO: Check that this is the correct way to deal with Rsur for Rtotal
    if rsur_calc_method == "single_layer":
        Rsur_c: list = [calc_Rsur(
            Rb,
            Rsto_c[TOP_LAYER],
            Rext[TOP_LAYER],
            Rinc[TOP_LAYER],
            Rgs,
            LAI_sum_per_layer[TOP_LAYER],
            SAI_sum_per_layer[TOP_LAYER],
        ) for _ in range(nL)]
    elif rsur_calc_method == "multi_layer":
        Rsur_c: List = calc_Rsur_multilayer(
            nL,
            Rb,
            Rsto_c,
            Rext,
            LAI_sum_per_layer,
            SAI_sum_per_layer
        )
    else:
        raise ValueError("Invalid Rsur calc method")

    Rtotal: List[float] = calc_Rtotal(nL, Rsur_c, Rinc, Rgs)

    return Resistance_Model(
        nL,
        Ra_c=Ra_c,
        Ra=Ra,
        Rb=Rb,
        Rinc=Rinc,
        Rext=Rext,
        Rsto=Rsto,
        Rgs=Rgs,
        Rsur=Rsur_c,
        Rtotal=Rtotal,
    )


def O3_ppb_to_nmol(
    Ts_C: float,
    P: float,
    O3_ppb: float,
):
    """Convert O3 ppb to nmol.

    Parameters
    ----------
    Ts_C : float
        Air Temperature [degressC]
    P : float
        Pressure [kPa]
    O3_ppb : float
        Ozone concentration in ppb [ppb]

    Returns
    -------
    float
        Ozone concentration [umol/m\u00b2/s]

    """
    M_O3 = 48.0      # Molecular weight of O3 (g)
    R = 8.314510  # Ideal gas constant
    T = deg_to_kel(Ts_C)

    Vn = R * (T / P)
    O3_nmol_m3 = (1.0 / Vn) * O3_ppb * M_O3 * 20.833
    return O3_nmol_m3


def calc_leaf_resistance_model(
    nL: int,
    Lm: float,
    u_per_layer: List[float],
    leaf_gsto_per_layer: List[float],
) -> Leaf_Resistance_Model:
    """Setup the leaf O3 resistance model.

    Parameters
    ----------
    nL: int
        Number of model layers
    Lm: float
        Leaf dimension
    u_per_layer: float
        wind speed for each layer [m/s]
    leaf_gsto_per_layer: List[float]
        leaf_gsto [mmol m-2 s-1]

    Returns
    -------
    Leaf_Resistance_Model
    """
    Rb_per_layer: List[float] = [calc_leaf_rb(calc_leaf_gb(
        LEAF_G_O3, Lm, u_per_layer[iL])) for iL in range(nL)]

    # TODO: Implement advanced Rext equation
    # Calculating Rext currently disabled to align with DO3SE ui
    # Rext_per_layer: List[List[float]] = [[calc_Rext(1.0)
    #                                       for iLC in range(nLC)]
    #                                      for iL in range(nL)]

    # TODO: Get Rext_const from config
    Rext_per_layer: List[float] = [2500 for iL in range(nL)]

    Rsto_per_layer: List[float] = [calc_Rsto(leaf_gsto_per_layer[iL]) for iL in range(nL)]

    leaf_r_model = Leaf_Resistance_Model(
        Rb=Rb_per_layer,
        Rext=Rext_per_layer,
        Rsto=Rsto_per_layer,
    )
    return leaf_r_model


def calc_leaf_resistance_model_old(
    nL: int,
    nLC: int,
    Lm_per_component: List[float],  # LC%Lm
    u_per_layer: List[float],  # met%u
    leaf_gsto_values: List[List[float]],
) -> List[List[Leaf_Resistance_Model]]:
    """Setup the leaf O3 resistance model.

    Parameters
    ----------
    nL: int
        Number of model layers
    nLC: int
        Number of model components
    Lm_per_component: float
        Lm for each component
    u_per_layer: float
        wind speed for each layer [m/s]
    leaf_gsto_values: float
        leaf_gsto for each layer and component [mmol m-2 s-1]

    Returns
    -------
    Leaf_Resistance_Model with shape(nL, nLC)
    """
    Rb_per_layer: List[List[float]] = [
        [calc_leaf_rb(calc_leaf_gb(LEAF_G_O3, Lm_per_component[iLC], u_per_layer[iL]))
            for iLC in range(nLC)]
        for iL in range(nL)
    ]

    # Calculating Rext currently disabled to align with DO3SE ui
    # Rext_per_layer: List[List[float]] = [[calc_Rext(1.0)
    #                                       for iLC in range(nLC)]
    #                                      for iL in range(nL)]

    Rext_per_layer: List[List[float]] = [[2500
                                          for iLC in range(nLC)]
                                         for iL in range(nL)]

    Rsto_per_layer: List[List[float]] = [[calc_Rsto(leaf_gsto_values[iL][iLC])
                                          for iLC in range(nLC)]
                                         for iL in range(nL)]

    leaf_r_models = [[Leaf_Resistance_Model(
        Rb=Rb_per_layer[iL][iLC],
        Rext=Rext_per_layer[iL][iLC],
        Rsto=Rsto_per_layer[iL][iLC],
    ) for iLC in range(nLC)]
        for iL in range(nL)]
    return leaf_r_models


def stomatal_flux_rate(
    leaf_r_Rsto: float,
    leaf_r_Rext: float,
    leaf_r_Rb: float,
) -> float:
    """Calculate the rate of stomatal flux for a leaf resistance model.

    Parameters
    ----------
    leaf_r_Rsto : float
        Leaf resistance model Rsto
    leaf_r_Rext : float
        Leaf resistance model Rext
    leaf_r_Rb : float
        Leaf resistance model Rb

    Returns
    -------
    fst: float
        Leaf resistance stomatal flux
    """
    leaf_r = 1.0 / ((1.0 / leaf_r_Rsto) + (1.0 / leaf_r_Rext))
    flux_rate = (1.0 / leaf_r_Rsto) * (leaf_r / (leaf_r_Rb + leaf_r))
    return flux_rate


PODOutput = namedtuple('PODOutput', 'POD_0 POD_Y')


def calc_POD(
    POD_0_prev: float,
    POD_Y_prev: float,
    has_emerged: bool,
    is_growing: bool,
    Y: float,
    fst: float,
) -> PODOutput:
    """Calculate the Phytotoxic Ozone Dose values.
    Calculate the accumulated stomatal flux above threshold Y

    Calc_AFstY in DO3SE UI model

    Parameters
    ----------
    POD_0_prev : float
        POD_0 from previous hour
    POD_Y_prev : float
        POD_Y from previous hour
    has_emerged: bool
        leaf has emerged
    Y : float
        Phytotoxic Ozone Dose threshold
    fst : float
        Stomatal ozone flux (nmol O3 m-2 PLA s-1)

    Returns (As namedtuple)
    -------
    POD_0: float
        Phytotoxic Ozone Dose, no threshold [mmol m-2 PLA]
    POD_Y: float
        Phytotoxic Ozone Dose above threshold Y [mmol m-2 PLA]
    """
    if not has_emerged or is_growing:
        return PODOutput(0, 0)
    POD_0 = POD_0_prev + ((fst * DT) / 1000000)
    POD_Y = POD_Y_prev + ((max(0.0, fst - Y) * DT) / 1000000)
    return PODOutput(
        POD_0=POD_0,
        POD_Y=POD_Y,
    )


def calc_OT(
    is_daylight: bool,
    f_phen: float,
    leaf_f_phen: float,
    micro_O3: float,
) -> Tuple[float, float]:
    OT_0 = micro_O3 / 1000 if is_daylight and leaf_f_phen > 0 else 0
    OT_40 = max(0.0, (micro_O3 - 40) / 1000) if is_daylight and f_phen > 0 else 0
    return OT_0, OT_40


OT_out = namedtuple('OT_out', 'OT_0 OT_40 AOT_0 AOT_40')


def calc_OT_acc(
    is_daylight: bool,
    f_phen: float,
    leaf_f_phen: float,
    micro_O3: float,
    AOT_0_prev: float,
    AOT_40_prev: float,
) -> NamedTuple:
    """Calculate the accumulated OT when is daylight and leaf f_phen > 0.

    Parameters
    ----------
    is_daylight : bool
        True if is daylight
    f_phen : float
        f_phen ?
    leaf_f_phen : float
        leaf_f_phen ?
    micro_O3: float
        O3 at layer [ppb]
    AOt_0_prev: float
        Accumulated OT_0 at previous hour
    AOt_40_prev: float
        Accumulated OT_40 at previous hour

    Returns (As namedtuple)
    -------
    OT_0: float
        layer O3 / 1000
    OT_40: float
        layer (O3 - 40) / 1000
    AOt_0: float
        Accumulated OT_0
    AOt_40: float
        Accumulated OT_40

    """
    OT_0, OT_40 = calc_OT(
        is_daylight,
        f_phen,
        leaf_f_phen,
        micro_O3,
    )
    AOT_0 = AOT_0_prev + OT_0
    AOT_40 = AOT_40_prev + OT_40

    return OT_out(
        OT_0=OT_0,
        OT_40=OT_40,
        AOT_0=AOT_0,
        AOT_40=AOT_40,
    )


def calc_OT_leaf(
    is_daylight: bool,
    f_phen: float,
    leaf_f_phen: float,
    micro_O3: List[float],
    AOT_0_prev: float,
    AOT_40_prev: float,
    nL: int,
) -> float:
    OT_0_per_layer, OT_40_per_layer = list(zip(*[
        calc_OT(
            is_daylight,
            f_phen,
            leaf_f_phen,
            micro_O3[i],
        ) for i in range(nL)
    ]))
    OT_0 = sum(OT_0_per_layer)
    OT_40 = sum(OT_40_per_layer)
    AOT_0 = AOT_0_prev + OT_0
    AOT_40 = AOT_40_prev + OT_40

    return OT_out(
        OT_0=OT_0,
        OT_40=OT_40,
        AOT_0=AOT_0,
        AOT_40=AOT_40,
    )


def calc_FO3_eff(
    fst: float,
    FO3_eff_prev: float,
    F_0: float,
) -> float:
    """Calculate FO3_eff if using martin2000 O3 method

    Parameters
    ----------
    fst : float
        Stomatal ozone flux [nmol O3 m-2 PLA s-1]
    FO3_eff_prev : float
        FO3_eff from previous hour [nmol O3 m-2 PLA]
    F_0 : float
        Threshold flux for ozone damage [nmol O3 m-2 PLA s-1]

    Returns
    -------
    FO3_eff: float
        (Accumulated) effective ozone dose [nmol O3 m-2 PLA]
    """
    FO3_eff = max(0.0, FO3_eff_prev + (fst - F_0) * DT)
    return FO3_eff


def calc_o3_dose(
    nL,
    nLC,
    O3_per_layer: List[float],
    leaf_rmodel_nLnLC: List[List[Leaf_Resistance_Model]],
    leaf_gsto_nLnLC: List[List[float]],
    f_phen_nLnLC: List[List[float]],
    leaf_f_phen_nLnLC: List[List[float]],
    Y_nLnLC: List[List[float]],
    POD_0_prev_nLnLC: List[List[float]],
    POD_Y_prev_nLnLC: List[List[float]],
    AOT_0_prev_nLnLC: List[List[float]],
    AOT_40_prev_nLnLC: List[List[float]],
    FO3_eff_prev_nLnLC: List[List[float]],
    F_0_nLnLC: List[List[float]],
    O3_method: str,
    is_daylight: bool,
):
    """Calculate the O3 dose.


    # TODO: fix reliance on MAX_RSTO...
    """
    raise NotImplementedError()
    # O3_nmol_per_layer = [O3_ppb_to_nmol(O3_per_layer[iL]) for iL in range(nLC)]
    # fst_nLnLC = [[calc_fst_alt(O3_nmol_per_layer[iL],
    #                            leaf_rmodel_nLnLC[iL][iLC],
    #                            leaf_gsto_nLnLC[iL][iLC])
    #               for iLC in range(nLC)]
    #              for iL in range(nL)]

    # POD_0_nLnLC = [[POD_0_prev_nLnLC[iL][iLC] + ((fst_nLnLC[iL][iLC] * DT) / 1000000)
    #                 for iLC in range(nLC)]
    #                for iL in range(nL)]
    # POD_Y_nLnLC = [[POD_Y_prev_nLnLC[iL][iLC] + ((max(0.0, fst_nLnLC[iL][iLC]
    # - Y_nLnLC[iL][iLC]) * DT) / 1000000)
    #                 for iLC in range(nLC)]
    #                for iL in range(nL)]

    # # Default OT0/OT40 to 0
    # OT_0_default = 0
    # OT_40_default = 0
    # OT_0_nLnLC = [[
    #     O3_per_layer[iL] / 1000
    #     if is_daylight and leaf_f_phen_nLnLC[iL][iLC] > 0 else OT_0_default
    #     for iLC in range(nLC)]
    #     for iL in range(nL)]

    # OT_40_nLnLC = [[
    #     max(0.0, O3_per_layer[iL] - 40) / 1000
    #     if is_daylight and f_phen_nLnLC[iL][iLC] > 0 else OT_40_default
    #     for iLC in range(nLC)]
    #     for iL in range(nL)]

    # AOT_0_nLnLC = [[AOT_0_prev_nLnLC[iL][iLC] + OT_0_nLnLC[iL][iLC]
    #                 for iLC in range(nLC)]
    #                for iL in range(nL)]

    # AOT_40_nLnLC = [[AOT_40_prev_nLnLC[iL][iLC] + OT_40_nLnLC[iL][iLC]
    #                  for iLC in range(nLC)]
    #                 for iL in range(nL)]

    # FO3_eff_nLnLC = [[max(0.0, FO3_eff_prev_nLnLC[iL][iLC]
    # + (fst_nLnLC[iL][iLC] - F_0_nLnLC) * DT)
    #                   if O3_method == 'martin2000' else None
    #                   for iLC in range(nLC)]
    #                  for iL in range(nL)]
    # return AOT_0_nLnLC, POD_Y_nLnLC, POD_0_nLnLC, AOT_40_nLnLC, FO3_eff_nLnLC


def calc_fst(
    Gsto_l: float,
    Rb_l: float,
    Rsto_l: float,
    Rext: float,
    O3_nmol_m3: float,
) -> float:
    """Calculate the upper leaf stomatal ozone flux (Hourly)

    This is used for Ewert Photosynthesis.

    Fst = O3(h) * gsto * (rc/(rb + rc))

    Parameters
    ----------
    # Config Parameters
    Rext: float
        External plant cuticle resistance [s/m]


    # State Variables
    Rb_l: float
        Leaf boundary resistance for this layer [s/m]
    Gsto_l: float
        (leaf_gsto) Leaf stomatal conductance [mmol/m^2/s]
    Rsto_l: float
        Leaf stomatal ozone resistance [s/m]
    O3_nmol_m3: float
        Ozone concentration at canopy [nmol/m^3]

    Returns
    -------
    Fst: float
        Upper leaf stomatal O3 flux (Fst, nmol/m\u00b2/s)

    """

    if (Gsto_l > 0):
        leaf_r = 1.0 / ((1.0 / Rsto_l) + (1.0 / Rext))  # leaf resistance in s/m
        Fst = O3_nmol_m3 * (1 / Rsto_l) * (leaf_r / (Rb_l + leaf_r))
    else:
        Fst = 0
    return Fst


def calc_fst_leaf(
    Gsto_l: List[float],
    Rb_l: List[float],
    Rsto_l: List[float],
    Rext: List[float],
    O3_nmol_m3: List[float],
    fLAI: List[float],
    nL: int,
) -> float:
    """Calculate the fst for a specific leaf population

    Parameters
    ----------
    Gsto_l : float
        The leaf stomatal conductance for this leaf population [mmol/m\u00b2/s]
    Rb_l: float
        Leaf boundary resistance for this layer [s/m]
    Rsto_l : List[float]
        stomatal resistance per layer[s/m]
    Rext : List[float]
        external plant cuticle resistance per layer[s/m]
    O3_nmol_m3 : List[float]
        Ozone concentration at each layer [nmol/m\u00b3]
    fLAI : float
        fraction of leaf population at each layer
    nL: int
        Number of layers

    Returns
    -------
    float
        leaf stomatal O3 flux (Fst) [nmol/m\u00b2/s]

    """
    total_fLAI = sum(fLAI)
    if total_fLAI == 0:
        return 0
    assert isclose(total_fLAI, 1.0, abs_tol=1e-3)

    fst_per_layer = [
        calc_fst(
            Gsto_l[i],
            Rb_l[i],
            Rsto_l[i],
            Rext[i],
            O3_nmol_m3[i],
        ) * fLAI[i] for i in range(nL)
    ]
    fst = sum(fst_per_layer) / sum(fLAI)
    return fst


def calc_O3up_accumulation(
    O3up: float,
    O3up_prev: float,
    O3up_acc: float,
    td_dd: float,
    td_dd_prev: float,
) -> float:
    """Calculate the accumulated Ozone uptake.

    As per Ewerr & Porter Equation 17

    Parameters
    ----------
    O3up_out: float
        Ozone uptake [nmol O3 PLA m-2 s-1]
    O3up_prev: float
        Ozone uptake from previous hour [nmol O3 PLA m-2 s-1]
    O3up_acc: float
        Accumulated Ozone uptake ([nmol O3 PLA m-2 s-1])
    td_dd: float
        difference between current td and td at season_Astart [Thermal Time]
    td_dd_prev: float
        difference between previous hour td and td at season_Astart [Thermal Time]

    Returns
    -------
    O3up_acc: float
        Accumulated Ozone uptake ([nmol O3 PLA m-2 s-1])
    """
    O3up_integration = ((O3up + O3up_prev) / 2) * (td_dd - td_dd_prev)
    O3up_acc_out = O3up_acc + O3up_integration
    return O3up_acc_out


# # TODO: Should these come from config?
# #  Properties of vegetation over which O3 concentration is measured
# O3_h = 25     # Canopy height[m]


# def calc_O3_concentration(
#         O3_ppb_zR: float, uh_i: float, TsC: float, P: float, ustar: float) -> float:
#     """Calculate the ozone concentration at the canopy (Hourly)

#     This procedure results in the calculation of deposition velocity (Vd) and
#     the ozone concentration at the canopy in both parts-per-billion and
#     nmol/m^3

#     Constants: izR, v, DIFF_DO3, Pr, Ts_K

#     External_State
#     O3_ppb_zR(O3),
#     uh_i,
#     Ts_C,
#     P,
#     ustar

#     Variables: Ra, Rb, Rsur
#     Variables: Vd, O3_ppb, O3_nmol_m3, Vd_i, O3_ppb_i, Ra_ref_i, &
#                         Ra_ref, Ra_O3zR_i, Ra_tar_i
#     Parameters: O3_h, O3zR, h, d, zo
#     R: ra_simple, rb_func => rb

#     Taken from DO3SE-UI model
# """
#     M_O3: float = 48.0   # Molecular weight of O3 (g)

#     O3_d = calc_veg_displacement_height(O3_h)
#     O3_zo = calc_veg_roughness_length(O3_h)

#     # ustar over reference canopy
#     ustar_ref = estimate_ustar(uh_i, izR - O3_d, O3_zo)
#     # Ra between reference canopy and izR
#     Ra_ref_i = ra_simple(ustar_ref, O3_zo + O3_d, izR, O3_d)
#     # Rb for reference canopy
#     Rb_ref = rb_func(ustar_ref, DIFF_O3)
#     # Deposition velocity at izR over reference canopy
#     # (assuming that Rsur_ref = Rsur)
#     Vd_i = 1.0 / (Ra_ref_i + Rb_ref + Rsur)
#     # Ra between measurement height and izR
#     Ra_O3zR_i = ra_simple(ustar_ref, O3zR, izR, O3_d)
#     # O3 concentration at izR
#     O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_O3zR_i * Vd_i))
#     # Ra between target canopy and izR
#     # (ustar already calculated for target canopy)
#     Ra_tar_i = ra_simple(ustar, zo + d, izR, d)
#     # Deposition velocity at izR over target canopy
#     Vd = 1.0 / (Ra_tar_i + Rb + Rsur)
#     # O3 concentration at target canopy
#     # (Ra already calculated between canopy height and izR)
#     O3_ppb = O3_ppb_i * (1.0 - (Ra * Vd))

#     # Specific molar volume of an ideal gas at current temp + pressure
#     Vn = 8.314510 * ((Ts_C + Ts_K) / P)
#     # Convert to nmol/m^3
#     O3_nmol_m3 = (1.0 / Vn) * O3_ppb * M_O3 * 20.833  # 1 microgram O3 = 20.833 nmol / m ^ 3

#     !==========================================================================
#     ! Calculate the ozone concentration at the canopy
#     !
#     ! This procedure results in the calculation of deposition velocity (Vd) and
#     ! the ozone concentration at the canopy in both parts-per-billion and
#     ! nmol/m^3
#     !==========================================================================
#     subroutine Calc_O3_Concentration()
#         use Constants, only: izR, v, DO3, Pr, Ts_K
#         use Inputs, only: O3_ppb_zR, uh_i, Ts_C, P, ustar
#         use Inputs, only: estimate_ustar
#         use Variables, only: Ra, Rb, Rsur
#         use Variables, only: Vd, O3_ppb, O3_nmol_m3, Vd_i, O3_ppb_i, Ra_ref_i, &
#                              Ra_ref, Ra_O3zR_i, Ra_tar_i
#         use Parameters, only: O3zR, O3_d, O3_zo, h, d, zo
#         use R, only: ra_simple, rb_func => rb

#         real, parameter :: M_O3 = 48.0      ! Molecular weight of O3 (g)

#         real :: ustar_ref, Rb_ref, Vn

#         ! ustar over reference canopy
#         ustar_ref = estimate_ustar(uh_i, izR - O3_d, O3_zo)
#         ! Ra between reference canopy and izR
#         Ra_ref_i = ra_simple(ustar_ref, O3_zo + O3_d, izR, O3_d)
#         ! Rb for reference canopy
#         Rb_ref = rb_func(ustar_ref, DO3)
#         ! Deposition velocity at izR over reference canopy
#         ! (assuming that Rsur_ref = Rsur)
#         Vd_i = 1.0 / (Ra_ref_i + Rb_ref + Rsur)
#         ! Ra between measurement height and izR
#         Ra_O3zR_i = ra_simple(ustar_ref, O3zR, izR, O3_d)
#         ! O3 concentration at izR
#         O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_O3zR_i * Vd_i))
#         ! Ra between target canopy and izR
#         ! (ustar already calculated for target canopy)
#         Ra_tar_i = ra_simple(ustar, zo + d, izR, d)
#         ! Deposition velocity at izR over target canopy
#         Vd = 1.0 / (Ra_tar_i + Rb + Rsur)
#         ! O3 concentration at target canopy
#         ! (Ra already calculated between canopy height and izR)
#         O3_ppb = O3_ppb_i * (1.0 - (Ra * Vd))

#         ! Specific molar volume of an ideal gas at current temp + pressure
#         Vn = 8.314510 * ((Ts_C + Ts_K) / P)
#         ! Convert to nmol/m^3
#         O3_nmol_m3 = (1.0/Vn) * O3_ppb * M_O3 * 20.833  ! 1 microgram O3 = 20.833 nmol/m^3
#     end subroutine Calc_O3_Concentration

#     !==========================================================================
#     ! Calculate the total ozone flux to the vegetated surface
#     !==========================================================================
#     subroutine Calc_Ftot()
#         use Variables, only: Ftot, O3_nmol_m3, Vd

#         Ftot = O3_nmol_m3 * Vd
#     end subroutine Calc_Ftot


#     !==========================================================================
#     ! Calculate the accumulated stomatal flux above threshold Y
#     !==========================================================================
#     subroutine Calc_AFstY()
#         use Variables, only: Fst, AFstY, AFst0
#         use Parameters, only: Y

#         ! Fst == 0 if Gsto_l == 0 (and Gsto_l == 0 if leaf_fphen == 0), so no
#         ! need to checeaf_fphen
#         AFst0 = AFst0 + ((Fst*60*60)/1000000)
#         AFstY = AFstY + ((max(0.0, Fst - Y)*60*60)/1000000)
#     end subroutine Calc_AFstY

#     !==========================================================================
#     ! Calculate the accumulated OT40
#     !==========================================================================
#     subroutine Calc_AOT40()
#         use Inputs, only: R
#         use Variables, only: OT0, OT40, AOT0, AOT40, O3_ppb, fphen, leaf_fphen

#         ! Default OT0 and OT40 to 0
#         OT0 = 0
#         OT40 = 0

#         ! Only accumulate OT when global radiation > 50 W/m^2
#         if (R > 50.0) then
#             ! Only accumulate OT0 when leaf_fphen > 0
#             if (leaf_fphen > 0) then
#                 OT0 = O3_ppb / 1000
#             end if
#             ! Only accumulate OT40 when fphen > 0
#             if (fphen > 0) then
#                 OT40 = max(0.0, O3_ppb - 40) / 1000
#             end if
#         end if

#         ! Accumulate
#         AOT0 = AOT0 + OT0
#         AOT40 = AOT40 + OT40
#     end subroutine Calc_AOT40

#     !==========================================================================
#     ! Set fO3 to 1.0, so it is ignored by Gsto calculation
#     !==========================================================================
#     subroutine Calc_fO3_Ignore()
#         use Variables, only: fO3

#         fO3 = 1.0
#     end subroutine Calc_fO3_Ignore

#     !==========================================================================
#     ! Calculate fO3 for wheat
#     !==========================================================================
#     subroutine Calc_fO3_Wheat()
#         use Variables, only: AFst0, fO3

#         fO3 = ((1+(AFst0/11.5)**10)**(-1))
#     end subroutine Calc_fO3_Wheat

#     !==========================================================================
#     ! Calculate fO3 for potato
#     !==========================================================================
#     subroutine Calc_fO3_Potato()
#         use Variables, only: AOT0, fO3

#         fO3 = ((1+(AOT0/40)**5)**(-1))
#     end subroutine Calc_fO3_Potato

# end module O3

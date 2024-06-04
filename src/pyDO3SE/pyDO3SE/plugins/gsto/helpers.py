"""Helper functions for gsto processes."""

from collections import namedtuple
from math import exp, log
from typing import NamedTuple

from pyDO3SE.constants.physical_constants import PAR_Wm2_to_photons


def calc_f_temp(
    Ts_C: float,
    T_min: float,
    T_opt: float,
    T_max: float,
    fmin: float,
) -> float:
    """Temperature effect on stomatal conductance.

    Describes a sort of bell-shaped curve with minima at T_min and T_max and a
    maximum at T_opt.


    Parameters
    ----------
    Ts_C : float
        Ambient Air temperature [deg C]
    T_min : float
        minimum temperature
    T_opt : float
        optimum temperature
    T_max : float
        maximum temperature
    fmin : float
        Minimum f_temp

    Returns
    -------
    float
        Temperature effect on stomatal conductance
    """
    bt = (T_max - T_opt) / (T_opt - T_min)

    if Ts_C < T_max:
        f_temp = ((Ts_C - T_min) / (T_opt - T_min)) * ((T_max - Ts_C) / (T_max - T_opt))**bt
        f_temp = max(fmin, min(1.0, f_temp))
    else:
        f_temp = 0
    return f_temp


def calc_f_temp_square_high(
    Ts_C: float,
    T_min: float,
    T_opt: float,
    T_max: float,
    fmin: float,
) -> float:
    """Temperature effect on stomatal conductance without reduction after T_opt.

    Follows f_temp exactly before T_opt, but the high end of the function is
    square - it continues at 1.0 until dropping to 0.0 at T_max.


    Parameters
    ----------
    Ts_C : float
        Ambient Air temperature [deg C]
    T_min : float
        minimum temperature
    T_opt : float
        optimum temperature
    T_max : float
        maximum temperature
    fmin : float
        Minimum f_temp

    Returns
    -------
    float
        Temperature effect on stomatal conductance
    """

    if (Ts_C >= T_max):
        f_temp_square_high = 0.0
    elif (Ts_C >= T_opt):
        f_temp_square_high = 1.0
    else:
        f_temp_square_high = calc_f_temp(Ts_C, T_min, T_opt, T_max, fmin)

    return f_temp_square_high


def calc_leaf_f_light(
    f_lightfac: float,
    PAR: float,
) -> float:
    """Calculate solar irradiance effect on leaf stomatal conductance.

    Parameters
    ----------
    f_lightfac : float
        Single leaf f_light coefficient
    PAR : float
        Photosynthetically active radiation [W m-2]

    Returns
    -------
    float
        solar irradiance effect on leaf stomatal conductance
    """
    # Note: Conversion of PAR to ppfd
    leaf_f_light = 1.0 - exp(-f_lightfac * (PAR * PAR_Wm2_to_photons))
    return leaf_f_light


def calc_f_light(
    f_lightfac: float,
    PARsun: float,
    PARshade: float,
    LAIsunfrac: float,
) -> float:
    """Calculate solar irradiance effect on canopy stomatal conductance.

    Leaf stomatal conductances for sun and shade leaves are combined according
    to give a canopy mean.


    Parameters
    ----------
    f_lightfac : float
        Single leaf f_light coefficient
    PARsun : float
        PAR received by sunlit leaves [W m-2]
    PARshade : float
        PAR received by shaded leaves [W m-2]
    LAIsunfrac : float
        Fraction of canopy component that is sunlit

    Returns
    -------
    float
        solar irradiance effect on canopy stomatal conductance
    """
    # TODO: does this need albedo?
    Flightsun = calc_leaf_f_light(f_lightfac, PARsun)
    Flightshade = calc_leaf_f_light(f_lightfac, PARshade)

    f_light = LAIsunfrac * Flightsun + (1.0 - LAIsunfrac) * Flightshade
    return f_light


def calc_f_light_method(
    LAI: float,
    sinB: float,
    f_lightfac: float,
    PARsun: float,
    PARshade: float,
    LAIsunfrac: float,
    PAR: float,
) -> NamedTuple:
    """Calculates the f_light for multiplicative gsto calculations.

    Calc_Flight in DO3SE UI

    Parameters
    ----------
    LAI : float
        Leaf area index
    sinB : float
        sin of solar elevation
    f_lightfac : float
        single leaf f_light coefficient
    PARsun : float
        Sunlit PAR
    PARshade : float
        Shaded PAR
    LAIsunfrac : float
        fraction of LAI that is sunlit
    PAR : float
        Photosynthetically active radiation

    Returns
    -------
    NamedTuple
        f_light, leaf_f_light
    """
    Output = namedtuple('Output', 'f_light leaf_f_light')
    if LAI > 0 and sinB > 0:
        # Calculate f_light and leaf_f_light
        # TODO: attenuate PAR properly through the canopy
        f_light = calc_f_light(f_lightfac, PARsun, PARshade, LAIsunfrac)
        # TODO: "grassland multilayer" model used leaf_flight = Flightsun, i.e.
        #       leaf_f_light(gc%f_lightfac, ll%met%PARsun) - which version is right?
        leaf_f_light = calc_leaf_f_light(f_lightfac, PAR)
    else:
        f_light = 0.0
        leaf_f_light = 0.0

    return Output(
        f_light=f_light,
        leaf_f_light=leaf_f_light,
    )


def calc_gsto_leaf(
    gp_gmax: float,
    gp_leaf_f_phen: float,
    gp_f_O3: float,
    gp_leaf_f_light: float,
    gp_fmin: float,
    gp_f_temp: float,
    gp_f_VPD: float,
    gp_f_SW: float,
) -> float:
    """Mean stomatal conductance of top leaf [mmol O3 m-2 PLA s-1].

    Parameters
    ----------
    gp_gmax : float
        Maximum stomatal conductance [mmol O3 m-2 PLA s-1]
    gp_leaf_f_phen : float
        Phenology-related effect on leaf gsto [fraction]
    gp_f_O3 : float
        O3 effect on gsto [fraction]
    gp_leaf_f_light : float
        Irradiance effect on leaf gsto [fraction]
    gp_fmin : float
        Minimum stomatal conductance [fraction]
    gp_f_temp : float
        Temperature effect on gsto [fraction]
    gp_f_VPD : float
        VPD effect on gsto [fraction]
    gp_f_SW : float
        Soil water effect on gsto [fraction]

    Returns
    -------
    float
        Mean stomatal conductance of top leaf [mmol O3 m-2 PLA s-1]
    """
    gsto_leaf = gp_gmax * min(gp_leaf_f_phen, gp_f_O3) * gp_leaf_f_light * \
        max(gp_fmin, gp_f_temp * gp_f_VPD * gp_f_SW)
    return gsto_leaf


def calc_gsto_mean(
    gp_gmax: float,
    gp_gmorph: float,
    gp_f_phen: float,
    gp_f_light: float,
    gp_fmin: float,
    gp_f_temp: float,
    gp_f_VPD: float,
    gp_f_SW: float,
) -> float:
    """Mean canopy stomatal conductance (mmol O3 m-2 PLA s-1).

    TODO: can SMD canopy gsto just use this with f_SW=1.0 ?

    Parameters
    ----------
    gp_gmax : float
        Maximum stomatal conductance [mmol O3 m-2 PLA s-1]
    gp_gmorph : float
        Sun/shade morphology factor [fraction]
    gp_f_phen : float
        Phenology-related effect on gsto [fraction]
    gp_f_light : float
        Irradiance effect on gsto [fraction]
    gp_fmin : float
        Minimum stomatal conductance [fraction]
    gp_f_temp : float
        Temperature effect on gsto [fraction]
    gp_f_VPD : float
        VPD effect on gsto [fraction]
    gp_f_SW : float
        Soil water effect on gsto [fraction]

    Returns
    -------
    float
        Mean canopy stomatal conductance (mmol O3 m-2 PLA s-1)
    """
    gsto_mean = gp_gmax * gp_gmorph * gp_f_phen * gp_f_light * \
        max(gp_fmin, gp_f_temp * gp_f_VPD * gp_f_SW)
    return gsto_mean


def temp_dep(
    P_ref: float,
    T_ref: float,
    H_a: float,
    T: float,
    R: float,
) -> float:
    """Calculate parameter from temperature dependence curve.

    Parameters
    ----------
    P_ref (float):
        Parameter value at T_ref [unitless]
    T_ref (float):
        Reference temperature [degrees K]
    H_a (float):
        Activation energy [J/mol]
    T (float):
        Temperature [degrees K]
    R (float):
        Universal gas constant [J K-1 mol-1] (CONSTANT)
    Returns
    -------
    P (float):
        [unitless]
    """
    P = P_ref * exp((H_a * (T - T_ref)) / (T_ref * R * T))
    return P


def temp_dep_inhibit(
    P_ref: float,
    T_ref: float,
    H_a: float,
    H_d: float,
    S: float,
    T: float,
    R: float,
) -> float:
    """Calculate parameter from temperature dependence curve.

    with high temperature inhibition

    Parameters
    ----------
    P_ref: float
        Parameter value at T_ref
    T_ref: float
        Reference temperature [degrees K]
    H_a: float
        Activation energy [J/mol]
    H_d: float
        Deactivation energy [J/mol]
    S: float
        Entropy term [J/(mol*K)]
    T: float
        Temperature [degrees K]
    R (float):
        Universal gas constant (J K-1 mol-1) (CONSTANT)
    Returns
    -------
    P (float):
        [Description]
    """
    P = P_ref * exp((H_a * (T - T_ref)) / (T_ref * R * T)) * \
        ((1 + exp((T_ref * S - H_d) / (T_ref * R))) / (1 + exp((T * S - H_d) / (T * R))))
    return P


def f_VPD_linear(
    VPD: float,
    VPD_max: float,
    VPD_min: float,
    fmin: float,
) -> float:
    """Calculate the VPD effect on stomatal conductance.

    from a linear function
    from (VPD_max, 1.0) to (VPD_min, fmin) - VPD_max should be smaller than
    VPD_min.

    Parameters
    ----------
    VPD : float
        Vapour pressure defecit [kPa]
    VPD_max : float
        VPD for maximum gsto [kPa]
    VPD_min : float
        VPD for minimum gsto [kPa]
    fmin : float
        Minimum f_VPD

    Returns
    -------
    float
        VPD effect on stomatal conductance
    """
    slope = (1.0 - fmin) / (VPD_min - VPD_max)
    f_VPD_linear_out = max(fmin, min(1.0, (fmin + (VPD_min - VPD) * slope)))
    return f_VPD_linear_out


def f_VPD_log(
    VPD: float,
    fmin: float,
) -> float:
    """Calculate VPD effect on gsto using simple logarithmic relationship.

    Parameters
    ----------
    VPD : float
        Vapour pressure defecit [kPa]
    fmin : float
        Minimum f_VPD

    Returns
    -------
    float
        VPD effect on stomatal conductance
    """
    f_VPD_log_out = max(fmin, min(1.0, (1.0 - 0.6 * log(VPD))))
    return f_VPD_log_out


def inverse_f_VPD_linear(
        f_VPD: float,
        VPD_max: float,
        VPD_min: float,
        fmin: float,
) -> float:
    """Calculate the VPD for a given f_VPD based on the linear VPD model.

    Parameters
    ----------
    f_VPD: float
        ___
    VPD_max: float
        VPD for maximum gsto [kPa]
    VPD_min: float
        VPD for minimum gsto [kPa]
    fmin: float
        Minimum f_VPD

    Returns
    -------
    D_0: float
        __

    """
    slope: float = (1.0 - fmin) / (VPD_min - VPD_max)
    D_0 = VPD_min - ((f_VPD - fmin) / slope)
    return D_0


def inverse_f_VPD_log(
    f_VPD: float,
    fmin: float,
) -> float:
    """Calculate the VPD for a given f_VPD.

    based on the simple logarithmic VPD model.

    Parameters
    ----------
    f_VPD: float
        [Description]
    fmin: float
        Minimum f_VPD

    Returns
    -------
    D_0: float
        [Description]

    """
    # TODO: fmin not having any effect here!
    D_0 = exp((1 - f_VPD) / 0.6)
    return D_0


def f_SWP_exp(
    fmin: float,
    SWP: float,
    a: float = 0.355,
    b: float = -0.706,
) -> float:
    r"""Calculate the soil water effect on stomatal conductance.

    using an exponential relationship.

    \f[
        f_{SW} = \min(1, \max(f_{min}, a \times (-SWP)^b))
    \f]

    Parameters
    ----------
    a : float
        Parameterization constant set in config
    b : float
        Parameterization constant set in config
    fmin : float
        Minimum f_SWP
    SWP : float
        Soil water potential [MPa]

    Returns
    -------
    float
        f_SWP
    """
    f_SWP_exp_out = min(1.0, max(fmin, a * (-SWP)**b))
    return f_SWP_exp_out


def inverse_f_SWP_exp(
    a: float,
    b: float,
    f_SWP: float,
) -> float:
    """Calculate the SWP (MPa) for a given f_SWP.

    based on the exponential relationship.

    Useful for determining SWP_max when using exponential f_SWP instead of
    linear f_SWP.

    Parameters
    ----------
    a : float
        [description]
    b : float
        [description]
    f_SWP : float
        Target f_SWP

    Returns
    -------
    inverse_f_SWP_exp: float
        [description]
    """

    inverse_f_SWP_exp_out = -(f_SWP / a)**(1.0 / b)
    return inverse_f_SWP_exp_out


def f_SWP_linear(
    SWP_min: float,
    SWP_max: float,
    fmin: float,
    SWP: float,
) -> float:
    """Calculate the soil water effect on stomatal conductance with a linear.

    function between (SWP_min, fmin) and (SWP_max, 1.0).

    Parameters
    ----------
    SWP_min : float
        SWP for minimum gsto [MPa]
    SWP_max : float
        SWP for maximum gsto [MPa]
    fmin : float
        Minimum f_SWP
    SWP : float
        Soil water potential [MPa]

    Returns
    -------
    f_SWP_linear: float
        [description]
    """
    f_SWP_linear_val = fmin + (1 - fmin) * ((SWP_min - SWP) / (SWP_min - SWP_max))
    f_SWP_linear_out = min(1.0, max(fmin, f_SWP_linear_val))
    return f_SWP_linear_out


def f_PAW(
    ASW_FC: float,
    fmin: float,
    ASW: float,
    ASW_min: float,
    ASW_max: float,
) -> float:
    """Calculate the soil water effect on stomatal conductance.

    with a linear relationship based on available soil water (ASW).


    (Bueker et al 2012)
    2.5.2 Plant available water method (PAW)
    In this approach fSW is assumed to be related to PAW (where PAW=FC-PAWmin)...

    Parameters
    ----------
    ASW_FC : float
        Available soil water at field capacity [m]
    fmin : float
        Minimum f_PAW
    ASW : float
        Available soil water [m]

    Returns
    -------
    f_PAW: float
        soil water effect on stomatal conductance


    References
    ----------
    Bueker et al (2012)


    """
    f_PAW_val = fmin + (1 - fmin) * ((100 * (ASW / ASW_FC)) - ASW_min) / (ASW_max - ASW_min)
    f_PAW_out = min(1.0, max(fmin, f_PAW_val))
    return f_PAW_out


def inverse_f_PAW(
    ASW_FC: float,
    fmin: float,
    f_PAW: float,
    ASW_max: float = 50,
    ASW_min: float = 0,
) -> float:
    """Calculate the ASW (m) for a given f_PAW based on the f_PAW relationship.

    Parameters
    ----------
    ASW_FC : float
        Available soil water at field capacity [m]
    fmin : float
        Minimum f_PAW
    f_PAW : float
        Target f_PAW
    ASW_min: float = 0.0
        ASW (available soil water) for minimum gsto (percent of ASW at field capacity)
    ASW_max: float = 50.0
        ASW (available soil water) for maximum gsto (percent of ASW at field capacity)


    Returns
    -------
    inverse_f_PAW: float
        [description]
    """
    inverse_f_PAW_out = (ASW_min + ((f_PAW - fmin) / (1 - fmin))
                         * (ASW_max - ASW_min)) * (ASW_FC / 100)  # noqa:W503
    return inverse_f_PAW_out

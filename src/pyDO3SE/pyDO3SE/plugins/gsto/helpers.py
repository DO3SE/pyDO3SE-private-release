"""Helper functions for gsto processes."""

from math import exp


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

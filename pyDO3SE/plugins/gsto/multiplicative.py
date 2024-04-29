from typing import NamedTuple


class Output_Shape(NamedTuple):
    new_leaf_gsto: float
    new_mean_gsto: float


def calc_gsto_leaf(
    gmax: float,
    leaf_f_phen: float,
    f_O3: float,
    leaf_f_light: float,
    fmin: float,
    f_temp: float,
    f_VPD: float,
    f_SW: float,
) -> float:
    """Mean stomatal conductance of top leaf (mmol O3 m-2 PLA s-1).

    type(GstoParams_t), intent(in) :: gp    !< Multiplicative gsto parameters
    """
    gsto_leaf = gmax * min(leaf_f_phen, f_O3) * leaf_f_light * \
        max(fmin, f_temp * f_VPD * f_SW)
    return gsto_leaf


def calc_gsto_mean(
    gmax: float,
    gmorph: float,
    f_phen: float,
    f_light: float,
    fmin: float,
    f_temp: float,
    f_VPD: float,
    f_SW: float,
) -> float:
    """Mean canopy stomatal conductance (mmol O3 m-2 PLA s-1).

    TODO: can SMD canopy gsto just use this with f_SW=1.0 ?

    type(GstoParams_t), intent(in) :: gp    !< Multiplicative gsto parameters

    """
    gsto_mean = gmax * gmorph * f_phen * f_light * \
        max(fmin, f_temp * f_VPD * f_SW)
    return gsto_mean


def apply_VPD_crit(
    VPD_crit: float,
    VPD_dd: float,
    old_gsto: float,
    new_gsto: float,
) -> float:
    """Limit stomatal conductance if the daily accumulated VPD threshold has been exceeded.

    Parameters
    ----------
    VPD_crit : float
        Accumulated VPD threshold [kPa]
    VPD_dd : float
        Accumulated VPD for the day [kPa]
    old_gsto : float
        Previous stomatal conductance
    new_gsto : float
        New stomatal conductance

    Returns
    -------
    float
        gsto
    """
    if (VPD_dd >= VPD_crit):
        return min(old_gsto, new_gsto)
    else:
        return new_gsto


def multiplicative(
    # From Config
    VPD_crit: float,
    # From State
    VPD_dd: float,
    initial_leaf_gsto: float,
    initial_mean_gsto: float,
    gmax: float,
    gmorph: float,
    f_phen: float,
    f_light: float,
    fmin: float,
    f_temp: float,
    f_VPD: float,
    f_SW: float,
    leaf_f_phen: float,
    f_O3: float,
    leaf_f_light: float,
) -> Output_Shape:
    """Multiplicative model for calculating gsto.

    References
    ----------
    Simpson et al (2012)

    Parameters
    ----------
    VPD_crit : float
        Critical daily VPD threshold above which stomatal conductance will stop
        increasing [kPa].
    VPD_dd : float
        Daily VPD sum during daylight hours [kPa]
    initial_leaf_gsto : float
         Leaf stomatal conductance [mmol O3 m-2 PLA s-1]
    initial_mean_gsto : float
        Canopy mean stomatal conductance [mmol O3 m-2 PLA s-1]
    fmin : float
        Minimum stomatal conductance [fraction]
    gmax : float
        Maximum stomatal conductance [mmol O3 m-2 PLA s-1]
    gmorph : float
        Sun/shade morphology factor [fraction]
    f_phen : float
        Phenology-related effect on gsto [fraction]
    leaf_f_phen : float
        Phenology-related effect on leaf gsto [fraction]
    f_light : float
        Irradiance effect on gsto [fraction]
    leaf_f_light : float
        Irradiance effect on leaf gsto [fraction]
    f_temp : float
        Temperature effect on gsto [fraction]
    f_VPD : float
        VPD effect on gsto [fraction]
    f_SW : float
        Soil water effect on gsto [fraction]
    f_O3 : float
        O3 effect on gsto [fraction]

    Returns
    -------
    Output_Shape
        leaf_gsto, mean_gsto
    """
    gsto_leaf = calc_gsto_leaf(
        gmax,
        leaf_f_phen,
        f_O3,
        leaf_f_light,
        fmin,
        f_temp,
        f_VPD,
        f_SW,
    )
    gsto_mean = calc_gsto_mean(
        gmax,
        gmorph,
        f_phen,
        f_light,
        fmin,
        f_temp,
        f_VPD,
        f_SW,
    )
    # Estimate Canopy gsto from mean gsto
    new_leaf_gsto = apply_VPD_crit(VPD_crit, VPD_dd, initial_leaf_gsto, gsto_leaf)
    new_mean_gsto = apply_VPD_crit(VPD_crit, VPD_dd, initial_mean_gsto, gsto_mean)

    return Output_Shape(
        new_leaf_gsto=new_leaf_gsto,
        new_mean_gsto=new_mean_gsto,
    )

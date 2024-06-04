"""Penman Monteith Soil Moisture Model.

To calculate Ei, Es and Et we calculate their hourly rate and accumulate it for the day.
At the end of the day the final accumulated value is used as Ei for the next day.
This means that the Ei used in following calculation is actually the Ei from the previous day.
"""

from collections import namedtuple
from typing import NamedTuple
from math import exp, inf

from do3se_met.resistance import calc_Rb

from pyDO3SE.plugins.gsto.helpers import inverse_f_PAW, inverse_f_SWP_exp
from pyDO3SE.constants.physical_constants import DIFF_H2O, DRATIO_O3_H20, T0
from pyDO3SE.plugins.resistance.helpers import calc_Rsto, calc_Rsur
from pyDO3SE.plugins.resistance.model import Resistance_Model
from pyDO3SE.plugins.soil_moisture.helpers import soil_moisture_from_SWC
from pyDO3SE.plugins.soil_moisture.config import Soil_t


def get_initial_SWC(SWC_in: float = None, FC: float = None) -> float:
    """Set up soil water data source parameters"""
    SWC = FC if SWC_in is None else SWC_in
    assert SWC is not None
    return SWC


# TODO: Specify args instead of using soil_config
def PM_soil_moisture_calc(
    soil_config: Soil_t, PWP: float, root_depth: float, Sn_in: float, Sn_diff: float
) -> NamedTuple:
    """Apply change in soil moisture from previous day's P-M summary."""
    Output = namedtuple('Output', 'Sn SWP ASW SMD')
    output: Output = soil_moisture_from_SWC(soil_config, PWP, root_depth, Sn_in + Sn_diff)
    return output

# TODO: remove below as doesn't help clean up PM hourly calc
# def calc_evaporation_from_canopy(
#     Rn: float,
#     delta: float,
#     psychro: float,
#     Pair: float,
#     lmd: float,
#     Cair: float,
#     VPD: float,
#     Rb_H2O: float,
# ) -> float:
#     """Calculate the evaporation from the canopy using resistances.

#     Parameters
#     ----------
#     Rn : float
#         Net radiation (J)
#     delta : float
#         [description]
#     psychro : float
#         [description]
#     Pair : float
#         [description]
#     lmd : float
#         [description]
#     Cair : float
#         [description]
#     VPD : float
#         Vapour pressure deficit (Pa)
#     Rb_H2O : float
#         resistance model - Quasi-laminar boundary layer resistance [s m-1]

#     Returns
#     -------
#     Ei: float
#         Evaporation from canopy [m]
#     """
#     G = 0.1 * Rn

#     Et_1 = (delta * (Rn - G)) / lmd
#     Et_2 = 3600 * Pair * Cair * VPD / Rb_H2O / lmd

#     Ei_3 = delta + psychro
#     Ei = (Et_1 + Et_2) / Ei_3 / 1000
#     return Ei


def calc_PET_gsto(
    gmax: float,
    f_phen: float,
    f_light: float,
    f_temp: float,
    f_VPD: float,
    LAI: float,
) -> float:
    """Calculate the potential gsto limited by multiplicative fractions.

    This is only valid for the multiplicative model.

    Parameters
    ----------
    gmax: float
        Maximum stomatal conductance [mmol O3 m-2 PLA s-1]
    f_phen: float
        Phenology-related effect on gsto [fraction]
    f_light: float
        Irradiance effect on gsto [fraction]
    f_temp: float
        Temperature effect on gsto [fraction]
    f_VPD: float
        VPD effect on gsto [fraction]
    LAI: float
        Leaf area index [m2/m2]
    """
    Gsto_PEt = gmax * f_phen * f_light * f_temp * f_VPD * LAI
    Rsto_PEt = calc_Rsto(Gsto_PEt)
    return Rsto_PEt


def penman_monteith_hourly(
    Rn_MJ: float,
    P_kPa: float,
    Ts_C: float,
    esat_kPa: float,
    eact_kPa: float,
    VPD_kPa: float,
    rm_h2o_nL: int,
    rm_h2o_Rb: float,
    rm_h2o_Rinc_l0: float,
    rm_h2o_Rsto_l0: float,
    rm_h2o_Rgs: float,
    rm_pet_Rsto_l0: float,
    LAI: float,
    Es_blocked: bool,
    pm_state_Ei_acc: float,
    pm_state_Et_acc: float,
    pm_state_PEt_acc: float,
    pm_state_Es_acc: float,
    pm_state_Eat_acc: float,
) -> NamedTuple:
    """Hourly Penman-Monteith calculations for evaporation and transpiration.

    This model (probably) makes some one-layer assumptions, so don't allow
    multi-layer resistance model.

    TODO: should we still be using Rsoil, instead of rm.Rgs?
    NOTE: PEt only calculated in multiplicative model

    Parameters
    ----------
    Rn_MJ: float
        Net radiation (MJ)
    P_kPa: float
        Atmospheric pressure (kPa)
    Ts_C: float
        Air temperature [degrees C]
    esat_kPa: float
        Saturated vapour pressure (kPa)
    eact_kPa: float
        Actual vapour pressure (kPa)
    VPD_kPa: float
        Vapour pressure deficit (kPa)
    rm_h2o_nL: int
        number of layers in resistance model
    rm_h2o_Rb: float
        resistance model - Quasi-laminar boundary layer resistance [s m-1]
    rm_h2o_Rinc_l0: float
        resistance model (layer 0) - In-canopy aerodynamic resistance [s m-1]
    rm_h2o_Rsto_l0: float
        resistance model (layer 0) - Stomatal resistance [s m-1] (H20)
    rm_h2o_Rgs: float
        resistance model - Ground surface resistance [s m-1]
    rm_pet_Rsto_l0: float
        Potential transpiration
    LAI: float
        Leaf area index [m2 m-2]
    Es_blocked: bool
        Is soil evaporation blocked?
    pm_state_Ei_acc: float
        Penman montein Evaporation from canopy accumulated[m]
    pm_state_Et_acc: float
        Penman montein Plant transpiration accumulated[m]
    pm_state_PEt_acc: float
        Penman montein potential transpiration accumulated[m]
    pm_state_Es_acc: float
        Penman montein Soil surface evaporation accumulated[m]
    pm_state_Eat_acc: float
        Penman montein Evapotranspiration accumulated[m]

    Returns (As namedtuple)
    -------
    Ei: float
        Penman montein Evaporation from canopy[m]
    Et: float
        Penman montein Plant transpiration[m]
    Es: float
        Penman montein Soil surface evaporation[m]
    Eat: float
        Penman montein Evapotranspiration[m]
    Ei_acc: float
        Penman montein Evaporation from canopy accumulated[m]
    Et_acc: float
        Penman montein Plant transpiration accumulated[m]
    Es_acc: float
        Penman montein Soil surface evaporation accumulated[m]
    Eat_acc: float
        Penman montein Evapotranspiration accumulated[m]
    """
    Output = namedtuple(
        'Output', 'Ei_hr Et_hr PEt_hr Es_hr Eat_hr, Ei_acc Et_acc PEt_acc Es_acc Eat_acc')
    # real :: Tvir, delta, lambda, psychro, Pair, Cair, G

    # real :: Et_1, Et_2, Ei_3, Et_3
    # real :: t, Es_Rn, Es_G, Es_1, Es_2, Es_3
    # real :: SW_a, SW_s, SW_c, C_canopy, C_soil

    # This model (probably) makes some one-layer assumptions, so don't allow
    # multi-layer resistance model.
    assert rm_h2o_nL == 1

    # Convert units and associate values
    Rn = Rn_MJ * 1000000
    P = P_kPa * 1000
    esat = esat_kPa * 1000
    eact = eact_kPa * 1000
    VPD = VPD_kPa * 1000

    Rb_H2O = rm_h2o_Rb
    Rinc = rm_h2o_Rinc_l0
    Rsto_H2O = rm_h2o_Rsto_l0
    Rsto_PEt = rm_pet_Rsto_l0
    Rsoil = rm_h2o_Rgs
    Ei_acc = pm_state_Ei_acc
    Et_acc = pm_state_Et_acc
    # Es_acc = pm_state_Es_acc
    Eat_acc = pm_state_Eat_acc

    should_calculate_PEt = rm_pet_Rsto_l0 is not None

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

    if should_calculate_PEt:
        PEt_3 = delta + psychro * (1 + Rsto_PEt / Rb_H2O)
        PEt_hr = (Et_1 + Et_2) / PEt_3 / 1000

    if Es_blocked:
        Es_hr = 0
    else:
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lmd
        Es_2 = (
            (3600 * Pair * Cair * VPD) - (delta * Rinc * ((Rn - G) - (Es_Rn - Es_G)))
        ) / (Rinc + Rb_H2O) / lmd
        Es_3 = delta + (psychro * (1.0 + (Rsoil / (Rb_H2O + Rinc))))
        Es_hr = (Es_1 + Es_2) / Es_3 / 1000

    # Calculate Eat from Et and Es (after Shuttleworth and Wallace, 1985)
    SW_a = (delta + psychro) * Rb_H2O
    SW_s = (delta + psychro) * Rinc + (psychro * Rsoil)
    SW_c = psychro * Rsto_H2O  # Boundary layer resistance = 0
    C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
    C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
    if Es_hr <= 0:
        Eat_hr = Et_hr
    else:
        Eat_hr = (C_canopy * Et_hr) + (C_soil * Es_hr)

    # Accumulate values
    Ei_acc = pm_state_Ei_acc + Ei_hr
    Et_acc = pm_state_Et_acc + Et_hr
    PEt_acc = should_calculate_PEt and pm_state_PEt_acc + PEt_hr
    Es_acc = pm_state_Es_acc + Es_hr
    Eat_acc = pm_state_Eat_acc + Eat_hr

    return Output(
        Ei_hr=Ei_hr,
        Et_hr=Et_hr,
        PEt_hr=should_calculate_PEt and PEt_hr,
        Es_hr=Es_hr,
        Eat_hr=Eat_hr,
        Ei_acc=Ei_acc,
        Et_acc=Et_acc,
        PEt_acc=should_calculate_PEt and PEt_acc,
        Es_acc=Es_acc,
        Eat_acc=Eat_acc,
    )


def penman_monteith_daily(
    LAI: float,
    root_depth: float,
    run_off_fraction: float,
    ASW: float,
    SMD: float,
    pm_state_precip_acc: float,
    pm_state_run_off_acc: float,
    pm_state_Ei_acc: float,
    pm_state_Eat_acc: float,
    pm_state_percolated_acc: float,
):
    """Daily soil water content update from accumulated Penman-Monteith values.

    Sets state.Sn_diff (change in soil water content) and clears daily
    accumulators.  state.Sn_diff isn't constrained by field capacity here, so
    that consideration must be handled elsewhere.

    Also processes run-off and deep percolation according to SMD configuration.
    Both "most recent" and "accumulated" values for these are stored in state.

    Parameters
    ----------
    LAI: float
        Leaf area index [m2 m-2]
    root_depth: float
        Root depth [m]
    run_off_fraction: float
        Amount of precipitation that is lost as run-off [fraction]
    ASW: float
        Current available soil water, above PWP [m]
    SMD: float
        Current soil moisture defecit, from FC [m]
    pm_state_precip_acc: float
        Accumulated precipitation [m]
    pm_state_run_off_acc: float
        Accumulated run-off [m]
    pm_state_Ei_acc: float
        Accumulated evaporation from canopy [m]
    pm_state_Eat_acc: float
        Accumulated evapotranspiration [m]
    pm_state_percolated_acc: float
        Accumulated deep percolation [m]

    Returns (as named tuple)
    -------
    rain_input: float
        Input from rainfall + irrigation [m]
    run_off: float
        Loss to run-off [m]
    run_off_acc: float
        Accumlated run-off [m]
    effective_irrig: float
        Effective irrigation [m]
    intercepted_evaporated: float
        Loss to evaporation of intercepted [m]
    evapotranspiration: float
        Loss to evapotranspiration
    Sn_diff: float
        Latest change in Sn [m3 m-3]
    percolated: float
        Loss to deep percolation [m]
    percolated_acc: float
        Accumulated deep percolation [m]
    """
    Output = namedtuple(
        'Output',
        [
            'rain_input',
            'run_off',
            'run_off_acc',
            'effective_irrig',
            'intercepted_evaporated',
            'evapotranspiration',
            'Sn_diff',
            'percolated',
            'percolated_acc',
        ],
    )
    # real :: max_ET, delta_SM

    # Start with full amount of precipitation
    # Estimate loss to run-off
    rain_input = pm_state_precip_acc
    run_off = run_off_fraction * rain_input
    run_off_acc = pm_state_run_off_acc + run_off

    # Calculate "effective irrigation"
    effective_irrig = rain_input - run_off

    # Estimate loss of intercepted precipitation to evaporation.  Intercepted
    # precipitation is estimated as 0.0001*LAI, which is therefore a limit on
    # how much can be evaporated.
    intercepted_evaporated = min(effective_irrig, 0.0001 * LAI, pm_state_Ei_acc)

    # Can't lose water below PWP, so constrain evapotranspiration to ASW by
    # restricting evapotranspiration.
    max_ET = ASW + effective_irrig - intercepted_evaporated
    evapotranspiration = min(max_ET, pm_state_Eat_acc)

    # Total balance = input - run_off - evaporated - evapotranspiration
    delta_SM = effective_irrig - intercepted_evaporated - evapotranspiration
    # Converted to volumetric change using root depth.
    Sn_diff = delta_SM / root_depth

    # Amount that will go to deep percolation = remainder if water balance
    # refills soil water, i.e. if it is greater than SMD.
    percolated = max(0.0, delta_SM - SMD)
    percolated_acc = pm_state_percolated_acc + percolated
    return Output(
        rain_input=rain_input,
        run_off=run_off,
        run_off_acc=run_off_acc,
        effective_irrig=effective_irrig,
        intercepted_evaporated=intercepted_evaporated,
        evapotranspiration=evapotranspiration,
        Sn_diff=Sn_diff,
        percolated=percolated,
        percolated_acc=percolated_acc,
    )


def penman_monteith_reset() -> NamedTuple:
    """Reset daily accumulated Penman-Monteith values.

    # Outputs
    - Ei_acc
    - Et_acc
    - Es_acc
    - Eat_acc
    """
    Output = namedtuple('Output', 'Ei_acc Et_acc Es_acc Eat_acc')
    # Clear accumulated variables
    return Output(
        Ei_acc=0.0,
        Et_acc=0.0,
        Es_acc=0.0,
        Eat_acc=0.0,
    )


def multi_layer_r_model_to_single_H20(
    rmodel: Resistance_Model,
    # nL: int,
    ustar: float,
    total_LAI: float,
    total_SAI: float,
) -> Resistance_Model:
    """Adapt multi-layer O3 resistance model to single-layer H2O resistance."""
    # Below only needed to be called once in model?
    Ra: float = rmodel.Ra
    Rb: float = calc_Rb(ustar, DIFF_H2O)
    # Combine multi-layer into single-layer
    # NOTE: We have to manage division by zero here
    Rinc_sum = sum([(1.0 / r) if r > 0 else inf for r in rmodel.Rinc[0:rmodel.nL]])
    Rinc: float = 1.0 / Rinc_sum if Rinc_sum > 0 else inf
    Rext_sum = sum([(1.0 / r) if r > 0 else inf for r in rmodel.Rext[0:rmodel.nL]])
    Rext: float = 1.0 / Rext_sum if Rext_sum > 0 else inf

    # Our gsto is in O3nmol so we must convert it to H2Onmol for H2O resistance model
    # TODO: Not sure this conversion is correct
    Rsto: float = DRATIO_O3_H20 * 1.0 / sum([(1.0 / r) for r in rmodel.Rsto[0:rmodel.nL]])
    # Rsto: float = 1.0 / sum([(1.0 / r) for r in rmodel.Rsto[0:rmodel.nL]])

    Rgs: float = rmodel.Rgs
    Rsur: float = calc_Rsur(Rb,
                            Rsto,
                            Rext,
                            Rinc,
                            Rgs,  # Rsoil = Rgs
                            total_LAI,
                            total_SAI)

    rmodel_H20 = Resistance_Model(1, Ra=Ra, Rb=Rb, Rinc=[Rinc], Rext=[
                                  Rext], Rsto=[Rsto], Rgs=Rgs, Rsur=[Rsur])
    return rmodel_H20


def check_soil_evaporation_blocked(
    f_SW_method: str,
    SWP: float = None,
    SWP_max: float = None,
    ASW: float = None,
    ASW_FC: float = None,
    fmin: float = None,
    fSWP_exp_a: float = None,
    fSWP_exp_b: float = None,
) -> bool:
    """Is soil evaporation blocked?"""
    # TODO: this assumes the first land cover is the only one that matters
    # TODO: Implement pn_gsto methods
    Es_blocked: bool = True

    if f_SW_method in ['fSWP exp', 'fLWP exp']:
        assert SWP is not None
        assert fSWP_exp_a is not None
        assert fSWP_exp_b is not None
        Es_blocked = SWP < inverse_f_SWP_exp(fSWP_exp_a, fSWP_exp_b, 1.0)
    if f_SW_method == 'fSWP linear':
        assert SWP is not None
        assert SWP_max is not None
        Es_blocked = SWP < SWP_max
    if f_SW_method == 'fPAW':
        assert ASW is not None
        assert ASW_FC is not None
        assert fmin is not None
        # ASW_MAX = 50.0
        # Es_blocked = (ASW < (ASW_FC * (ASW_MAX / 100.0)))
        # TODO: Below has different output for DO3SE UI
        Es_blocked = ASW < inverse_f_PAW(ASW_FC, fmin, 1.0)
    if f_SW_method == "disabled":
        # TODO: This is set to calculate to match DO3SE UI
        # Es_blocked = SWP < SWP_max
        Es_blocked = True

    # TODO: Move these checks to config validation
    if f_SW_method is not None and \
            f_SW_method not in ['fSWP exp', 'fLWP exp', 'fSWP linear', 'fPAW', 'disabled']:
        raise ValueError(f'{f_SW_method} is an invald fSWP linear method')

    # Ensure output is bool (Issue caused by numpy bool_)
    return bool(Es_blocked)

"""Functions specific to the plant Phenology i.e plant life cycle using the JULES/Ewert Phyllochron method.



    ------------------------t_l---------------|
    -----t_lem---|----------t_lma-------------|
                 |-------t_lep--------|-t_lse-|
                 |                    |       |
    = = = = = = = - - - - - - - - - - |       |
                \\                            |
                  \\                   \      |
                    \\                        |
                      \\                \     |
    f_LA = = =          \\                    |
    f_LS - - -            \\             \    |
                            \\                |
                              \\          \   |
                                \\            |
                                  \\       \  |
                                    \\        |
                                      \\    \ |
                                        \\    |
                                          \\ \|
                                            \\|
"""


import numpy as np
from math import ceil
import warnings
from enum import Enum
from typing import List, NamedTuple, Tuple
from deprecated import deprecated
from collections import namedtuple
from thermal_time.calcs import calc_effective_temperature
from .units import *

from do3se_phenology.state import LeafPhenologyStage, PhenologyStage
from do3se_phenology.config import PhenologyKeyLengths

def phyllochron_to_t_lem(phyllochron: float) -> float:
    """Get the leaf emerging leaf time interval from the phyllochron.

    .. math:: t_lem = 1.8 * phyllochron

    Parameters
    ----------
    phyllochron : float
        :data:`phyllochron`

    Returns
    -------
    t_lem: float
        Emerging leaf thermal interval [deg days]

    """
    t_lem = 1.8 * phyllochron
    return t_lem


def phyllochron_to_t_lma(phyllochron: float) -> float:
    """Get the leaf mature life time inteval from the phyllochron.

    Parameters
    ----------
    phyllochron : float
        See docs

    Returns
    -------
    t_lma: float
        Mature leaf thermal interval [deg days]
    """
    t_lma = 3.5 * phyllochron
    return t_lma


def calc_t_l(t_lem: float, t_lma: float) -> float:
    """Total lifespan of the leaf from the thermal interval of the emerging leaf and mature leaf.

    Parameters
    ----------
    t_lem: float
        Emerging leaf thermal interval [deg days]
    t_lma: float
        Mature leaf thermal interval [deg days]

    Returns
    -------
    t_l: float
        Life span of leaf in thermal time [deg days]
    """
    t_l = t_lem + t_lma
    return t_l


def calc_t_l_from_phyllochron(phyllochron: float) -> float:
    """Calculate the life span of the leaf from the phyllochron.

    Parameters
    ----------
    phyllochron : float
        the Phyllochron

    Returns
    -------
    t_l: float
        Life span of leaf in thermal time [deg days]
    """
    if phyllochron is None:
        # phyllochron not set yet
        return None
    t_lem = phyllochron_to_t_lem(phyllochron)
    t_lma = phyllochron_to_t_lma(phyllochron)
    t_l = calc_t_l(t_lem, t_lma)
    return t_l


def calc_life_stages_from_phyllochron(
        phyllochron: float, t_lse_constant: float) -> float:
    """Calculate the life span of the leaf from the phyllochron.

    Parameters
    ----------
    phyllochron : float
        the Phyllochron
    t_lse_constant: float
        t_lse as fraction of t_l

    Returns (As named tuple)
    -------
    t_l: float
        Life span of leaf in thermal time [deg days]
    t_lem: float
        Emerging leaf thermal interval [deg days]
    t_lep: float
        Time during which leaf is expanding [thermal time]
    t_lse: float
        Time that the leaf is senescing [thermal time]
    t_lma: float
        Mature leaf thermal interval [deg days]

    """
    Output = namedtuple('output', 't_l t_lem t_lep t_lse t_lma')
    if phyllochron is None:
        # phyllochron not set yet
        return Output(99999, 99999, 99999, 99999, 99999)
    t_lem = phyllochron_to_t_lem(phyllochron)
    t_lma = phyllochron_to_t_lma(phyllochron)
    t_l = calc_t_l(t_lem, t_lma)
    t_lse = t_lma * t_lse_constant
    t_lep = t_lma - t_lse

    return Output(
        t_l=t_l or None,
        t_lem=t_lem or None,
        t_lep=t_lep or None,
        t_lse=t_lse or None,
        t_lma=t_lma or None,
    )


def calc_phyllochron(dl: float) -> float:
    """Calculate the phyllochron from the change in day length.

    Parameters
    ----------
    dl : float
        Change in day length from photoperiod calculation

    Returns
    -------
    phyllochron: float
        See docs
    """
    phyllochron = 1 / (0.0104 + 0.026 * dl)
    return phyllochron


def calc_phyllochron_from_dvi(
        prev_dvi: float,
        dvi: float,
        prev_pr: float,
        pr: float,
        prev_phyllochron: float
) -> float:
    """Check if phyllochron needs setting.

    Will be set when DVI goes from < 0 to > 0

    Parameters
    ----------
    prev_dvi : float
        The previous day development index
    dvi : float
        The development index
    prev_pr : float
        The previous day photoperiod value
    pr : float
        The photoperiod
    prev_phyllochron : [type]
        The current set phyllochron value

    Returns
    -------
    Phyllochron: float
        The phyllochron value
    """
    if prev_dvi > 0 and not prev_phyllochron:
        raise ValueError('Phyllachron should have been set!')
    if prev_phyllochron:  # or null
        return prev_phyllochron
    if prev_dvi < 0 and dvi > 0:
        return calc_phyllochron(pr - prev_pr)


def calc_emergence_date_from_dvi(
        prev_dvi: float,
        dvi: float,
        prev_t_emerg: float,
        prev_d_emerg: int,
        td: float,
        dd: int,
) -> Tuple:
    """Check if t_emerg needs setting.

    Will be set when DVI goes from < 0 to > 0
    Sets both the julian day and thermal time

    Parameters
    ----------
    prev_dvi : float
        The previous day development index
    dvi : float
        The development index
    prev_t_emerg : float
        Previous t_emerg value
    td: float
        Current td value

    Returns
    -------
    t_emerg: float
        emergence thermal time
    """
    if prev_dvi > 0 and not prev_t_emerg:
        raise ValueError('emergence date should have been set!')
    if prev_t_emerg:  # or null
        return prev_t_emerg, prev_d_emerg
    if prev_dvi < 0 and dvi > 0:
        return td, dd
    return None, None


def calc_rpe(p: float, p_crit: float, p_sens: float) -> float:
    """Calculate the relative photoperiod effect(RPE).

    Parameters
    ----------
    p : float
        photoperiod (day length)
    p_crit : float
        critical photoperiod
    p_sens : float
        sensitivity of development rate to photoperiod

    Returns
    -------
    rpe: float
        relative photoperiod effect
    """
    rpe = 1 - (p - p_crit) * p_sens
    return rpe


def calc_dvi(
        prev_dvi: float,
        t_eff: float,
        tt_emr: float,
        tt_veg: float,
        tt_rep: float,
        rpe: float,
        dd: int,
        sowing_day: int,
) -> float:
    """Calculate the development index (DVI).

    Parameters
    ----------
    prev_dvi : float
        Previous day development index [dimensionless]
    t_eff : float
        effective temperature [degrees]
    tt_emr : float
        thermal time interval between sowing and emergence [deg days]
    tt_veg : float
        thermal time interval between emergence and flowering [deg days]
    tt_rep : float
        thermal time interval between flowering and maturity/harvest [deg days]
    rpe : float
        relative photoperiod effect [dimensionless]
    dd: int
        Current julian day
    sowing_day: int
        Sowing day

    Returns
    -------
    dvi: float
        Development index [dimensionless]

    """
    dvi_change = 0

    assert tt_emr is not None
    assert tt_veg is not None
    assert tt_rep is not None

    if dd <= sowing_day:
        # Before sowing date DVI is -1
        return -1

    # TODO: t_eff should be 0 at SGS
    if -1 <= prev_dvi < 0:
        dvi_change = t_eff / tt_emr
    elif 0 <= prev_dvi < 1:
        dvi_change = (t_eff / tt_veg) * rpe
    elif 1 <= prev_dvi < 2:
        dvi_change = t_eff / tt_rep
    else:
        return 2
        # TODO: Handle dvi above 2
        # raise ValueError('DVI should not be more than 2')

    dvi = prev_dvi + dvi_change
    return dvi if dvi <= 2 else 2


def estimate_t_l_from_t_lse(t_lse: float):
    """Estimate t_l from t_lse by reversing Ewert, Porter and Jules equations.

    t_lma = t_lse / 0.33
    t_lma = 3.5 * phyllochron
    phyllochron = t_lma / 3.5
    t_lem = 1.8 * phyllochron
    t_lem = 1.8 *  t_lma / 3.5
    t_l = t_lma + t_lem
    t_l = (t_lse/0.33) + 1.8 * ( t_lse / 0.33)/3.5
    t_l = t_lse(3 +1.5584)
    t_l = 4.55844 * t_lse

    Parameters
    ----------
    t_lse : float
        [description]

    Returns
    -------
    t_l: float
        [description]

    """
    t_l = 4.55844 * t_lse
    return t_l


def estimate_t_lse_from_t_l(
    t_l: float,
):
    """Estimate the value of t_lse from t_l using  Ewert, Porter and Jules equations.

    t_l   = t_lma + t_lem

    t_lma = 3.5 * phyllochron
    t_lem = 1.8 * phyllochron
    t_l   = 5.3 * phyllochron

    t_lse = 0.33 * t_lma
    t_lse = 0.33 * 3.5 * phyllochron
    t_lse = (t_l / 5.3) * 0.33 * 3.5


    t_lma = t_lse / 0.33
    t_lma = 3.5 * phyllochron
    phyllochron = t_lma / 3.5
    t_lem = 1.8 * phyllochron
    t_lem = 1.8 *  t_lma / 3.5
    t_l = t_lma + t_lem
    t_l = (t_lse/0.33) + 1.8 * ( t_lse / 0.33)/3.5
    t_l = t_lse(3 +1.5584)
    t_l = 4.55844 * t_lse
    t_lse = t_l / 4.55844

    """
    t_lse = t_l / 4.55844
    return t_lse


def calculate_day_from_td(
    td: List[float],
    dd: List[int],
    td_dd: float,
):
    d = None

    for(
        dd,
        td,
    ) in zip(
            dd, td,
    ):
        if td >= td_dd:
            d = dd
            break
    return d


class Leaf_f_Phen_Stages(Enum):
    START = 0
    EMERGING = 1
    MATURE = 2
    SENESCENCE = 3
    END = 4


def calculate_t_l_from_leaf_f_phen(
    td_full: List[float],
    leaf_f_phen: List[float],
    dd_full: List[int],
    f_t_lse: Fraction = 0.33,
    f_t_lem: Fraction = 0.54,  # TODO: Should be 0.54 as td_percentage
) -> NamedTuple:
    """Estimate the t_l when given leaf_f_phen input data.

    This method assumes t_lma matches length of leaf_f_phen.

    Parameters
    ----------
    td_full : List[float]
        thermal time for full dataset
    leaf_f_phen : List[float]
        leaf_f_phen for full dataset
    dd : List[int]
        dd for full dataset
    f_t_lse: Fraction
        fraction of t_lma that is t_lse

    Returns
    -------
    namedtuple
        Key phenology stages

    """
    Output = namedtuple('Output', [
        "t_lse",
        "t_l",
        "t_lem",
        "t_lep",
        "t_lma",
        "AStart",
        "t_emerg",
        "d_emerg",
    ])
    # Get the gradient of leaf_f_phen at each time point
    leaf_f_phen_diff = [0] + [a - b for a, b in zip(leaf_f_phen[1:], leaf_f_phen[:-1])]
    # row_indexes = [i for i in range(len(leaf_f_phen))]
    t_lma_start = None
    t_lma_end = None
    AStart = None
    f_phen_start = None

    for (
        leaf_f_phen_dd,
        leaf_f_phen_diff,
        td,
        dd,
    ) in zip(
        leaf_f_phen,
        leaf_f_phen_diff,
        td_full,
        dd_full,
    ):
        if leaf_f_phen_diff > 0 and t_lma_start is None:
            f_phen_start = dd
            AStart = dd - 1

            # When leaf_f_phen first increases above 0 we assume this is the t_l start
            t_lma_start = td
        if t_lma_end is None and f_phen_start is not None and leaf_f_phen_diff == 0 and leaf_f_phen_dd == 0:
            t_lma_end = td

    if t_lma_start is None:
        warnings.warn("Leaf f phen data does not start with 0")
        t_lma_end = td_full[0]
    if t_lma_end is None:
        warnings.warn("Leaf f phen data does not end with 0")
        t_lma_end = td_full[-1]
    assert t_lma_end is not None
    t_lma = t_lma_end - t_lma_start
    t_lse = t_lma * f_t_lse
    t_l = t_lma / (1 - f_t_lem)
    t_lem = t_l - t_lma
    t_lep = t_lma - t_lse
    t_emerg = t_lma_end - t_l
    d_emerg = calculate_day_from_td(td_full, dd_full, t_emerg)

    return Output(
        t_lse=t_lse,
        t_l=t_l,
        t_lem=t_lem,
        t_lep=t_lep,
        t_lma=t_lma,
        AStart=AStart,
        t_emerg=t_emerg,
        d_emerg=d_emerg,
    )


@deprecated(version="0.0.6", reason="WIP and not checked")
def calculate_t_l_from_leaf_f_phen_b(
    Ts_C_full: List[float],
    td_full: List[float],
    leaf_f_phen: List[float],
    dd_full: List[int],
    tt_emr: float,
    t_b: float,
    t_o: float,
    t_m: float,
):
    """Estimate the t_l when given leaf_f_phen input data.

    WIP METHOD Using Tte, Tts and Ttg from Porter 1984
    This method does not result in the correct relationship between t_l and t_lse. See tests.


    Parameters
    ----------
    td_full : List[float]
        thermal time for full dataset
    leaf_f_phen : List[float]
        leaf_f_phen for full dataset
    dd : List[int]
        dd for full dataset
    t_lse_gradient : float, optional
        the gradient at which senesence(t_lse) is occuring, by default -0.01
    method: str, optional
        the method to use, 't_l' or 't_lse', default = t_l
        t_l => assumes t_l is the length between start and end of leaf_f_phen
        t_lse => assumes t_lse is the length of leaf_f_phen when gradient is -ve

    Returns
    -------
    namedtuple('Output', 't_lse t_l t_lem t_lma t_lep sowing_day SGS')
        Key phenology stages

    """
    Output = namedtuple('Output', 't_lse t_l t_lem t_lma t_lep sowing_day SGS')
    # Get the gradient of leaf_f_phen at each time point
    leaf_f_phen_diff = [0] + [a - b for a, b in zip(leaf_f_phen[1:], leaf_f_phen[:-1])]
    row_indexes = [i for i in range(len(leaf_f_phen))]

    sowing_day = None
    t_l_start = None
    t_l_end = None
    SGS = None

    fphen_start_row_index = None
    td_start = None
    td_end_emergence = None
    Ttg = 0
    Tts = 0
    Tte = 0

    stage = Leaf_f_Phen_Stages.START

    for (
        row_index,
        leaf_f_phen_row,
        leaf_f_phen_diff,
        td,
        dd,
    ) in zip(
        row_indexes,
        leaf_f_phen,
        leaf_f_phen_diff,
        td_full,
        dd_full,
    ):
        if stage == Leaf_f_Phen_Stages.START and leaf_f_phen_diff > 0:
            stage = Leaf_f_Phen_Stages.EMERGING
            t_l_start = td
            fphen_start_row_index = row_index
            SGS = dd - 1  # TODO: This is an arbitary number to make LAI work
            td_start = td

        if stage == Leaf_f_Phen_Stages.EMERGING and leaf_f_phen_row == 1:
            stage = Leaf_f_Phen_Stages.MATURE
            Ttg = td - td_start
            td_end_emergence = td

        if stage == Leaf_f_Phen_Stages.MATURE and leaf_f_phen_diff < 0:
            stage = Leaf_f_Phen_Stages.SENESCENCE
            Tte = td - td_end_emergence

        if stage == Leaf_f_Phen_Stages.SENESCENCE and leaf_f_phen_diff == 0:
            stage = Leaf_f_Phen_Stages.END
            Tts = td - td_end_emergence
            t_l_end = td

    assert stage == Leaf_f_Phen_Stages.END

    t_l = t_l_end - t_l_start
    t_lse = Tts - Tte  # estimate_t_lse_from_t_l(t_l)
    t_lma = Tts  # t_lse / 0.33
    t_lem = Ttg  # 1.8 * (t_lma / 3.5)
    t_lep = Tte  # t_lma - t_lse

    phyllochron = t_lma / 3.5

    # phyllochron = 1 / (0.0104 + 0.026 * dl)
    # d1 = ((1 / phyllochron) - 0.0104) / 0.026

    dvi = 0
    # Reverse dvi equation
    da = dd_full[0]
    db = dd_full[-1]
    t_eff_full = [calc_effective_temperature(
        sum(Ts_C_full[row_index:row_index + 24]), t_b, t_o, t_m) for _ in range(24) for row_index in range(dd_full[0] * 24, dd_full[-1] * 24, 24)]

    for t_eff, dd in reversed(list(zip(t_eff_full[0:fphen_start_row_index], dd_full[0:fphen_start_row_index]))):
        if dvi < -1:
            sowing_day = dd + 1
            break
        dvi_diff = t_eff / tt_emr
        dvi -= dvi_diff

    if sowing_day is None:
        raise ValueError("Sowing day is before start day of data so must be set in config")

    return Output(
        t_lse=t_lse,
        t_l=t_l,
        t_lem=t_lem,
        t_lep=t_lep,
        t_lma=t_lma,
        sowing_day=sowing_day,
        SGS=SGS,
    )


def calc_dvi_range(teff_data, rpe_data, dd_data, hrs_data, tt_emr, tt_veg, tt_rep, sowing_day):
    dvi = [-1]
    for t, rp, d, hr in zip(teff_data, rpe_data, dd_data, hrs_data):
        if hr == 0:
            new_dvi = calc_dvi(dvi[-1], t, tt_emr, tt_veg, tt_rep, rp, d, sowing_day)
            dvi.append(new_dvi)
    return dvi


def get_dvi_PLF(
    tt_emr: TimeUnit,
    tt_veg: TimeUnit,
    tt_rep: TimeUnit,
    x_offset: TimeUnit = 0,
) -> PiecewiseFunction[TimeUnit, float]:
    """Get a piecewise function for DVI."""
    x_values = [
        x_offset,
        tt_emr + x_offset,
        tt_emr + tt_veg + x_offset,
        tt_emr + tt_veg + tt_rep + x_offset,
    ]
    y_values = [
        -1.000000001,
        0.0,
        1.0,
        2.0,
    ]
    return x_values, y_values


def calc_emergence_rate(
    nP: int,
    t_emerg_to_flag_emerg: float
) -> float:
    """Calculate the rate of emergence as a function of td.

    Sub populations emergence is spread evenly between plant emergence and flag leaf emergence.

    If nP == 1 then this will be 0. In this case we only use flag leaf definition.

    If nP > 1 the nP-1 populations are spread between plant emergence and flag leaf emergence


    Parameters
    ----------

    nP: int
        Total number of leaf populations inc. flag leaf
    t_emerg_to_flag_emerg: float
        thermal time between plant emergence and flag leaf emergence


    Returns
    -------
    float
        Rate at which leaf populations emerge as a fn(td)

    """
    emergence_rate = (nP - 1) / t_emerg_to_flag_emerg
    return emergence_rate


def calc_emerged_leaf_count(
    nP: int,
    td_dd: float,
    emergence_rate: float,
    t_emerge_flag: float,
) -> int:
    """Calculate the number of emerged leaf populations.

    If nP is 1 then we use the flag leaf only.

    Parameters
    ----------
    nP: int
        Number of leaf populations
    td_dd : float
        Thermal time since plant emergence
    emergence_rate : float
        Rate at which leaf populations emerge as a fn(td)
    t_emerge_flag: float
        thermal time between plant emergence and flag leaf emergence

    Returns
    -------
    int
        Number of emerged leaf populations

    """
    emerged_leaf_count = max(0, min(nP, ceil(td_dd * emergence_rate))
                             ) if nP > 1 else 1 if td_dd > t_emerge_flag else 0
    return emerged_leaf_count


def calc_td_dd_per_leaf_pop(
    nP: int,
    emerged_leaf_count: int,
    td: float,
    td_prev: float,
    td_dd_prev: List[float],
    t_emerg_to_flag_emerg: float,
) -> List[float]:
    """Calculate the thermal time difference between leaf pop emergence and current.

    Parameters
    ----------
    nP : int
        Number of leaf populations
    emerged_leaf_count : int
        Number of leaf populations that have emerged
    td : float
        Thermal time since plant emergence
    td_dd_prev : float
        The previous leaf populations thermal time diff since emergence

    Returns
    -------
    List[float]
        Leaf populations thermal time difference since emergence

    """
    if nP == 1:
        return [max(0, td - t_emerg_to_flag_emerg)]
    td_diff = td - td_prev
    td_per_leaf_pop = [td_dd_prev[i] + td_diff if i < emerged_leaf_count else 0 for i in range(nP)]
    return td_per_leaf_pop


def get_growing_populations(
    td_dd_emerg: List[float],
    leaf_population_t_lems: List[float],
    flag_leaf_t_lem: float,
) -> List[bool]:
    """Calculate which leaf populations are still growing.

    True when thermal time is between leaf emergence and end of leaf t_lem.
    t_lem is the the length of time a leaf population grows.

    Parameters
    ----------
    td_dd_emerg : List[float]
        List of thermal time since leaf pop emergence for each leaf pop
    leaf_population_t_lems : List[float]
        List of the t_lem for each leaf population
    flag_leaf_t_lem: float
        t_lem for flag leaf

    Returns
    -------
    List[bool]
        Boolean list of which populations are growing

    """
    all_tlems = [*leaf_population_t_lems, flag_leaf_t_lem]
    growing_populations = [0 < td < t_lem for td, t_lem in zip(td_dd_emerg, all_tlems)]
    return growing_populations


def get_plant_phenology_stage(
    td_dd: float,
    key_lengths: PhenologyKeyLengths,
    phenology_stage: PhenologyStage,
) -> PhenologyStage:
    """Get the plant phenology stage.

    Parameters
    ----------
    td_dd : float
        thermal time since sowing
    key_lengths : PhenologyKeyLengths
        key phenology lengths
    phenology_stage: PhenologyStage
        Current phenology stage

    Returns
    -------
    PhenologyStage
        Output phenology stage

    """
    if phenology_stage < PhenologyStage.EMERGED:
        return phenology_stage
    elif td_dd <= key_lengths.emerg_to_astart:
        return PhenologyStage.EMERGED
    elif td_dd <= key_lengths.emerg_to_end:
        return PhenologyStage.ASTART
    else:
        return PhenologyStage.HARVEST

def get_leaf_phenology_stage(
    td_dd: float,
    t_lem: float,
    t_lep: float,
    t_lse: float,
) -> LeafPhenologyStage:
    if td_dd <= 0:
        return LeafPhenologyStage.NOT_EMERGED
    elif td_dd <= t_lem:
        return LeafPhenologyStage.GROWING
    elif td_dd <= t_lem + t_lep:
        return LeafPhenologyStage.MATURE
    elif td_dd <= t_lem + t_lep + t_lse:
        return LeafPhenologyStage.SENESCENCE
    else:
        return LeafPhenologyStage.FULLY_SENESED



def calc_dvi_tt_PLF(
    td: float,
    dvi_interval: List[float],
) -> float:
    dvi_x = [i[0] for i in dvi_interval]
    dvi_y = [i[1] for i in dvi_interval]

    dvi = np.interp(td, dvi_x, dvi_y)
    return dvi

def get_dvi_range_from_species_config(species_config, td):
    td_at_sgs = species_config.key_dates_td.sowing
    tt_emr = species_config.key_lengths_td.sowing_to_emerge
    tt_veg = species_config.key_lengths_td.emerg_to_veg
    tt_rep = species_config.key_lengths_td.veg_to_harvest

    dvi_x = [td_at_sgs, td_at_sgs + tt_emr, td_at_sgs +
             tt_emr + tt_veg, td_at_sgs + tt_emr + tt_veg + tt_rep]
    dvi_y = [-1, 0, 1, 2]
    dvi = np.interp(td, dvi_x, dvi_y)
    return dvi

""" Phenology definitions as a fraction of total season thermal time.


    f_Astart : Fraction, optional
        Fraction of total growing season between sowing and Astart, by default Wheat.f_Astart [UNIT]
    f_mid_anthesis : Fraction, optional
        Fraction of total growing season between sowing and mid anthesis, by default Wheat.f_mid_anthesis [UNIT]
    f_fphen_a : Fraction, optional
        Fraction of total growing season between sowing and plant emergence, by default Wheat.f_fphen_a [UNIT]
    f_fphen_b : Fraction, optional
        Fraction of total growing season between sowing and f_phen_b(Double Ridge), by default Wheat.f_fphen_b [UNIT]
    f_fphen_c : Fraction, optional
        Fraction of total growing season between sowing and Mid Anthesis, by default Wheat.f_fphen_c [UNIT]
    f_fphen_d : Fraction, optional
        Fraction of thermal time between sowing and harvest(Should always be 1.0), by default Wheat.f_fphen_d [UNIT]
    f_tt_emr : Fraction, optional
        Fraction of total growing season between sowing and Emergence, by default Wheat.f_tt_emr [UNIT]
    f_tt_veg : Fraction, optional
        Fraction of total growing season between sowing and the Vegatitive stage, by default Wheat.f_tt_veg [UNIT]
    f_tt_rep : Fraction, optional
        thermal time interval between flowering and maturity/harvest [deg days]
        Fraction of thermal time interval between flowering and maturity/harvest, by default Wheat.f_tt_rep [UNIT]


"""

from math import isclose
from typing import List, NamedTuple
from collections import namedtuple
import warnings
from .units import *
from .presets.wheat import Wheat

from do3se_phenology.f_phen import tt_f_phen_simple_PLF_value, tt_leaf_f_phen_PLF_value


def get_dvi_from_season_length(
    sgs: ThermalTime,
    egs: ThermalTime,
    f_tt_emr: Fraction = Wheat.f_tt_emr,
    f_tt_veg: Fraction = Wheat.f_tt_veg,
    f_tt_rep: Fraction = Wheat.f_tt_rep,
):
    season_length = egs - sgs

    tt_emr = season_length * f_tt_emr
    tt_veg = season_length * f_tt_veg
    tt_rep = season_length * f_tt_rep
    return tt_emr, tt_veg, tt_rep

CanopyTDIntervals = namedtuple('CanopyTDIntervals', [
        "t_Astart",
        "t_mid_anthesis",
        "t_fphen_a",
        "t_fphen_b",
        "t_fphen_c",
        "t_fphen_d",
        "tt_emr",
        "tt_veg",
        "tt_rep",
    ])

def get_canopy_td_intervals_f(
    t_sgs: ThermalTime,
    t_egs: ThermalTime,
    f_Astart: Fraction = Wheat.f_Astart,
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
    f_fphen_a: Fraction = Wheat.f_fphen_a,
    f_fphen_b: Fraction = Wheat.f_fphen_b,
    f_fphen_c: Fraction = Wheat.f_fphen_c,
    f_fphen_d: Fraction = Wheat.f_fphen_d,
    f_tt_emr: Fraction = Wheat.f_tt_emr,
    f_tt_veg: Fraction = Wheat.f_tt_veg,
    f_tt_rep: Fraction = Wheat.f_tt_rep,
) -> CanopyTDIntervals:
    """Get the canopy thermal time intervals using fractions of growing season.

    NOTE: t_sgs is thermal time at sowing date and t_egs is thermal time at end of canopy senesence.

    Parameters
    ----------
    t_sgs : ThermalTime
        Thermal time at sowing date [Thermal Time]
    t_egs : ThermalTime
        Thermal time at end of growing season(Harvest) [Thermal Time]
    f_Astart : Fraction, optional
        Fraction of total growing season between sowing and Astart, by default Wheat.f_Astart [UNIT]
    f_mid_anthesis : Fraction, optional
        Fraction of total growing season between sowing and mid anthesis, by default Wheat.f_mid_anthesis [UNIT]
    f_fphen_a : Fraction, optional
        Fraction of total growing season between sowing and plant emergence, by default Wheat.f_fphen_a [UNIT]
    f_fphen_b : Fraction, optional
        Fraction of total growing season between sowing and f_phen_b(Double Ridge), by default Wheat.f_fphen_b [UNIT]
    f_fphen_c : Fraction, optional
        Fraction of total growing season between sowing and Mid Anthesis, by default Wheat.f_fphen_c [UNIT]
    f_fphen_d : Fraction, optional
        Fraction of thermal time between sowing and harvest(Should always be 1.0), by default Wheat.f_fphen_d [UNIT]
    f_tt_emr : Fraction, optional
        Fraction of total growing season between sowing and Emergence, by default Wheat.f_tt_emr [UNIT]
    f_tt_veg : Fraction, optional
        Fraction of total growing season between sowing and the Vegatitive stage, by default Wheat.f_tt_veg [UNIT]
    f_tt_rep : Fraction, optional
        thermal time interval between flowering and maturity/harvest [deg days]
        Fraction of thermal time interval between flowering and maturity/harvest, by default Wheat.f_tt_rep [UNIT]

    Returns
    -------
    Named tuple
        Thermal time equivalents of fractions.
    """

    assert t_sgs is not None
    assert t_egs is not None
    season_td = t_egs - t_sgs
    t_Astart = season_td * f_Astart
    t_mid_anthesis = season_td * f_mid_anthesis
    t_fphen_a = season_td * f_fphen_a
    t_fphen_b = season_td * f_fphen_b
    t_fphen_c = season_td * f_fphen_c
    t_fphen_d = season_td * f_fphen_d
    tt_emr = season_td * f_tt_emr
    tt_veg = season_td * f_tt_veg
    tt_rep = season_td * f_tt_rep

    return CanopyTDIntervals(
        t_Astart=t_Astart,
        t_mid_anthesis=t_mid_anthesis,
        t_fphen_a=t_fphen_a,
        t_fphen_b=t_fphen_b,
        t_fphen_c=t_fphen_c,
        t_fphen_d=t_fphen_d,
        tt_emr=tt_emr,
        tt_veg=tt_veg,
        tt_rep=tt_rep,
    )


# NOT USED YET! #
def get_current_f_phen_from_t_sgs_t_egs(
    td: ThermalTime,
    t_sgs: ThermalTime,
    t_egs: ThermalTime,
    f_phen_min: Fraction = Wheat.f_phen_min,
    f_Astart: Fraction = Wheat.f_Astart,
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
    f_fphen_a: Fraction = Wheat.f_fphen_a,
    f_fphen_b: Fraction = Wheat.f_fphen_b,
    f_fphen_c: Fraction = Wheat.f_fphen_c,
    f_fphen_d: Fraction = Wheat.f_fphen_d,
) -> float:
    """Get the current phenology from thermal time and fractions.

    Parameters
    ----------
    td: ThermalTime,
        Current thermal time since model start
    t_sgs: ThermalTime,
        thermal time at season start(Assumed sowing date)
    t_egs: ThermalTime,
        thermal time at season end
    f_phen_min: Fraction,
        species specific parameter(See phenology fraction docs)
    f_Astart: Fraction = Wheat.f_Astart,
        species specific parameter(See phenology fraction docs)
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
        species specific parameter(See phenology fraction docs)
    f_fphen_a: Fraction = Wheat.f_fphen_a,
        species specific parameter(See phenology fraction docs)
    f_fphen_b: Fraction = Wheat.f_fphen_b,
        species specific parameter(See phenology fraction docs)
    f_fphen_c: Fraction = Wheat.f_fphen_c,
        species specific parameter(See phenology fraction docs)
    f_fphen_d: Fraction = Wheat.f_fphen_d,
        species specific parameter(See phenology fraction docs)

    Returns
    -------
    float
        Current f_phen

    """
    out = get_canopy_td_intervals_f(
        t_sgs=t_sgs,
        t_egs=t_egs,
        f_Astart=f_Astart,
        f_mid_anthesis=f_mid_anthesis,
        f_fphen_a=f_fphen_a,
        f_fphen_b=f_fphen_b,
        f_fphen_c=f_fphen_c,
        f_fphen_d=f_fphen_d,
    )
    f_phen = tt_f_phen_simple_PLF_value(
        td=td,
        t_f_phen_a=out.t_fphen_a,
        t_f_phen_b=out.t_fphen_b,
        t_f_phen_c=out.t_fphen_c,
        t_f_phen_d=out.t_fphen_d,
        f_phen_min=f_phen_min,
        td_at_sgs=t_sgs,
    )
    return f_phen


def get_current_leaf_f_phen_from_t_sgs_t_egs(
    td: ThermalTime,
    t_sgs: ThermalTime,
    t_egs: ThermalTime,
    f_Astart: Fraction = Wheat.f_Astart,
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
    f_fphen_1_ets: Fraction = Wheat.f_fphen_a,
    f_fphen_3_ets: Fraction = Wheat.f_fphen_b,
    f_fphen_4_ets: Fraction = Wheat.f_fphen_c,
    f_fphen_5_ets: Fraction = Wheat.f_fphen_d,
    f_leaf_f_phen_a: Fraction = Wheat.leaf_f_phen_a,
    f_leaf_f_phen_b: Fraction = Wheat.leaf_f_phen_b,
) -> float:
    """Get the current phenology from thermal time and fractions.

    Parameters
    ----------
    td: ThermalTime,
        Current thermal time since model start
    t_sgs: ThermalTime,
        thermal time at season start(Assumed sowing date)
    t_egs: ThermalTime,
        thermal time at season end
    f_Astart: Fraction = Wheat.f_Astart,
        species specific parameter(See phenology fraction docs)
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
        species specific parameter(See phenology fraction docs)
    f_fphen_1_ets: float
        species specific parameter(See phenology fraction docs)
    f_fphen_3_ets: float
        species specific parameter(See phenology fraction docs)
    f_fphen_4_ets: float
        species specific parameter(See phenology fraction docs)
    f_fphen_5_ets: float
        species specific parameter(See phenology fraction docs)
    f_leaf_f_phen_a: float
        species specific parameter(See phenology fraction docs)
    f_leaf_f_phen_b: float
        species specific parameter(See phenology fraction docs)


    Returns
    -------
    float
        Current flag leaf_f_phen

    """

    out = get_leaf_td_intervals_f(
        t_sgs=t_sgs,
        t_egs=t_egs,
        f_Astart=f_Astart,
        f_mid_anthesis=f_mid_anthesis,
        f_fphen_1_ets=f_fphen_1_ets,
        f_fphen_3_ets=f_fphen_3_ets,
        f_fphen_4_ets=f_fphen_4_ets,
        f_fphen_5_ets=f_fphen_5_ets,
    )
    leaf_f_phen = tt_leaf_f_phen_PLF_value(
        td=td,
        t_leaf_f_phen_a=f_leaf_f_phen_a,
        t_leaf_f_phen_b=f_leaf_f_phen_b,
        t_leaf_f_phen_e=out.t_fphen_1_ets,
        t_leaf_f_phen_g=out.t_fphen_3_ets,
        t_leaf_f_phen_h=out.t_fphen_4_ets,
        t_leaf_f_phen_i=out.t_fphen_5_ets,
        t_astart=out.t_Astart,
        td_at_sgs=t_sgs,
    )
    return leaf_f_phen


def get_canopy_td_intervals(
    td_list: List[float],
    dd_list: List[int],
    t_sgs: ThermalTime,
    t_egs: ThermalTime,
    f_Astart: Fraction = Wheat.f_Astart,
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
    f_fphen_a: Fraction = Wheat.f_fphen_a,
    f_fphen_b: Fraction = Wheat.f_fphen_b,
    f_fphen_c: Fraction = Wheat.f_fphen_c,
    f_fphen_d: Fraction = Wheat.f_fphen_d,
    f_tt_emr: Fraction = Wheat.f_tt_emr,
    f_tt_veg: Fraction = Wheat.f_tt_veg,
    f_tt_rep: Fraction = Wheat.f_tt_rep,
):
    """Get the canopy phenology stages in thermal time intervals using percentages of provided growing season.

    See tt_f_phen_simple_PLF for explaination of parameters.

    Note that tt_emr, tt_veg and tt_rep assume no photoperiod impact.

    Parameters
    ----------
    td_list : List[float]
        Daily thermal time values
    dd_list : List[int]
        julian day indexes
    t_sgs : ThermalTime
        thermal time at start of growing season
    t_egs : ThermalTime
        thermal time at end of growing season
    f_Astart : Fraction, optional
        fraction of thermal time at fully developed leaf(Astart), by default 0.54
    f_mid_anthesis : Fraction, optional
        fraction of thermal time at mid anthesis , by default 0.60
    f_fphen_a : Fraction, optional
        fraction of thermal time from SGS to plant emergence, by default 0.05
    f_fphen_b : Fraction, optional
        fraction of thermal time from SGS to Double ridge stage, by default 0.2
    f_fphen_c : Fraction, optional
        fraction of thermal time from SGS to Mid Anthesis, half way through flowering, by default 0.6
    f_fphen_d : Fraction, optional
        fraction of thermal time from SGS to Harvest, by default 1.0

    Returns
    -------
    NamedTuple
        phenology parameters in thermal time and julian days

    """
    Output = namedtuple('Output', [
        "SGS",
        "EGS",
        "Astart",
        "mid_anthesis",
        "fphen_a",
        "fphen_b",
        "fphen_c",
        "fphen_d",
        "dd_emr",
        "dd_veg",
        "dd_rep",
        "t_Astart",
        "t_mid_anthesis",
        "t_fphen_a",
        "t_fphen_b",
        "t_fphen_c",
        "t_fphen_d",
        "tt_emr",
        "tt_veg",
        "tt_rep",
    ])
    assert t_sgs is not None
    assert t_egs is not None
    season_td = t_egs - t_sgs

    t_Astart = season_td * f_Astart
    t_mid_anthesis = season_td * f_mid_anthesis
    t_fphen_a = season_td * f_fphen_a
    t_fphen_b = season_td * f_fphen_b
    t_fphen_c = season_td * f_fphen_c
    t_phen_d = season_td * f_fphen_d
    tt_emr = season_td * f_tt_emr
    tt_veg = season_td * f_tt_veg
    tt_rep = season_td * f_tt_rep

    sgs = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs))
    egs = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_egs))
    Astart = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + t_Astart))
    mid_anthesis = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + t_mid_anthesis))
    fphen_a = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + t_fphen_a))
    fphen_b = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + t_fphen_b))
    fphen_c = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + t_fphen_c))
    fphen_d = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + t_phen_d))

    dd_emr = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + tt_emr))
    dd_veg = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + tt_veg))
    dd_rep = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >= t_sgs + tt_rep))

    return Output(
        SGS=sgs,
        EGS=egs,
        Astart=Astart,
        mid_anthesis=mid_anthesis,
        fphen_a=fphen_a,
        fphen_b=fphen_b,
        fphen_c=fphen_c,
        fphen_d=fphen_d,
        dd_emr=dd_emr,
        dd_veg=dd_veg,
        dd_rep=dd_rep,
        t_Astart=t_Astart,
        t_mid_anthesis=t_mid_anthesis,
        t_fphen_a=t_fphen_a,
        t_fphen_b=t_fphen_b,
        t_fphen_c=t_fphen_c,
        t_fphen_d=t_phen_d,
        tt_emr=tt_emr,
        tt_veg=tt_veg,
        tt_rep=tt_rep,
    )


def get_leaf_td_intervals_f(
    t_sgs: ThermalTime,
    t_egs: ThermalTime,
    f_Astart: Fraction = Wheat.f_Astart,
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
    f_fphen_1_ets: Fraction = Wheat.f_fphen_1_ets,
    f_fphen_3_ets: Fraction = Wheat.f_fphen_3_ets,
    f_fphen_4_ets: Fraction = Wheat.f_fphen_4_ets,
    f_fphen_5_ets: Fraction = Wheat.f_fphen_5_ets,
    f_tt_fst_acc: Fraction = Wheat.f_tt_fst_acc,
    f_t_lem: Fraction = Wheat.f_t_lem,
    f_t_lma: Fraction = Wheat.f_t_lma,
    f_t_lep: Fraction = Wheat.f_t_lep,
    f_t_lse: Fraction = Wheat.f_t_lse,

):
    """Get the leaf phenology stages in thermal time intervals using percentages of provided growing season.

    Parameters
    ----------
    t_sgs : ThermalTime
        The start of the growing season in thermal time.
    t_egs : ThermalTime
        The end of the growing season in thermal time.
    f_Astart : Fraction, optional
        The fraction of the growing season at which Astart occurs, default is Wheat.f_Astart.
    f_mid_anthesis : Fraction, optional
        The fraction of the growing season at which mid-anthesis occurs, default is Wheat.f_mid_anthesis.
    f_fphen_1_ets : Fraction, optional
        The fraction of the growing season at which the first phenological stage occurs, default is Wheat.f_fphen_1_ets.
    f_fphen_3_ets : Fraction, optional
        The fraction of the growing season at which the third phenological stage occurs, default is Wheat.f_fphen_3_ets.
    f_fphen_4_ets : Fraction, optional
        The fraction of the growing season at which the fourth phenological stage occurs, default is Wheat.f_fphen_4_ets.
    f_fphen_5_ets : Fraction, optional
        The fraction of the growing season at which the fifth phenological stage occurs, default is Wheat.f_fphen_5_ets.
    f_tt_fst_acc: Fraction, optional
        The fraction of the t_lem at which fst_acc begins, default is Wheat.f_tt_fst_acc.
    f_t_lem : Fraction, optional
        The fraction of the growing season at which the leaf emergence occurs, default is Wheat.f_t_lem.
    f_t_lma : Fraction, optional
        The fraction of the growing season at which the leaf maturation occurs, default is Wheat.f_t_lma.
    f_t_lep : Fraction, optional
        The fraction of the growing season at which the leaf expansion occurs, default is Wheat.f_t_lep.
    f_t_lse : Fraction, optional
        The fraction of the growing season at which the leaf senescence occurs, default is Wheat.f_t_lse.

    Returns
    -------
    _type_
        _description_

    """
    Output = namedtuple('Output', [
        "t_Astart",
        "t_mid_anthesis",
        "t_fphen_1_ets",
        "t_fphen_3_ets",
        "t_fphen_4_ets",
        "t_fphen_5_ets",
        "t_fst_acc",
        "t_lem",
        "t_lse",
        "t_lma",
        "t_lep",
        "t_l",
        "t_emerg",
    ])
    assert t_sgs is not None
    assert t_egs is not None
    season_td = t_egs - t_sgs

    t_Astart = season_td * f_Astart
    t_mid_anthesis = season_td * f_mid_anthesis
    t_fphen_1_ets = season_td * f_fphen_1_ets
    t_fphen_3_ets = season_td * f_fphen_3_ets
    t_fphen_4_ets = season_td * f_fphen_4_ets
    t_fphen_5_ets = season_td * f_fphen_5_ets
    t_lem = season_td * f_t_lem
    t_lse = season_td * f_t_lse
    t_lma = season_td * f_t_lma
    t_lep = season_td * f_t_lep
    t_l = t_lem + t_lma
    t_fst_acc = t_lem * f_tt_fst_acc if f_tt_fst_acc is not None else None
    t_emerg = t_egs - t_l
    return Output(
        t_Astart=t_Astart,
        t_mid_anthesis=t_mid_anthesis,
        t_fphen_1_ets=t_fphen_1_ets,
        t_fphen_3_ets=t_fphen_3_ets,
        t_fphen_4_ets=t_fphen_4_ets,
        t_fphen_5_ets=t_fphen_5_ets,
        t_fst_acc=t_fst_acc,
        t_lem=t_lem,
        t_lse=t_lse,
        t_lma=t_lma,
        t_lep=t_lep,
        t_l=t_l,
        t_emerg=t_emerg,
    )


def get_leaf_td_intervals(
    td_list: List[float],
    dd_list: List[int],
    t_sgs: ThermalTime,
    t_egs: ThermalTime,
    f_Astart: Fraction = Wheat.f_Astart,
    f_mid_anthesis: Fraction = Wheat.f_mid_anthesis,
    f_fphen_1_ets: Fraction = Wheat.f_fphen_1_ets,
    f_fphen_3_ets: Fraction = Wheat.f_fphen_3_ets,
    f_fphen_4_ets: Fraction = Wheat.f_fphen_4_ets,
    f_fphen_5_ets: Fraction = Wheat.f_fphen_5_ets,
    f_t_lem: Fraction = Wheat.f_t_lem,
    f_t_lma: Fraction = Wheat.f_t_lma,
    f_t_lep: Fraction = Wheat.f_t_lep,
    f_t_lse: Fraction = Wheat.f_t_lse,
):
    """Get the leaf phenology stages in thermal time intervals using percentages of provided growing season.

    Parameters
    ----------
    td_list : List[float]
        Daily thermal time values
    dd_list : List[int]
        julian day indexes
    t_sgs : ThermalTime
        thermal time at start of growing season
    t_egs : ThermalTime
        thermal time at end of growing season
    f_Astart : Fraction, optional
        fraction of thermal time at fully developed leaf(Astart), by default 0.53
    f_mid_anthesis : Fraction, optional
        fraction of thermal time at mid anthesis , by default 0.63
    f_fphen_1_ets : Fraction, optional
        fraction of thermal time between Flag leaf fully developed(Astart) to mid anthesis, by default 0.12
    f_fphen_3_ets : Fraction, optional
        fraction of thermal time between Mid anthesis to start of flag leaf ageing, by default 0.05
    f_fphen_4_ets : Fraction, optional
        fraction of thermal time between Mid anthesis to start of flag leaf senesence, by default 0.22
    f_fphen_5_ets : Fraction, optional
        fraction of thermal time between Start of flag leaf senescence to harvest, by default 0.42
    f_t_lem : Fraction, optional
        fraction of thermal time to Astart used Ewert model, by default 0.53
    f_t_lma : Fraction, optional
        fraction of thermal time from Astart to harvest used Ewert model, by default 0.47
    f_t_lep : Fraction, optional
        fraction of thermal time from Astart to start of senesence used Ewert model, by default 0.35
    f_t_lse : Fraction, optional
        fraction of thermal time for senesence to harvest used Ewert model, by default 0.53

    Returns
    -------
    NamedTuple
        phenology parameters in thermal time and julian days
    """
    Output = namedtuple('Output', [
        "SGS",
        "EGS",
        "Astart",
        "mid_anthesis",
        "fphen_1_ets",
        "fphen_3_ets",
        "fphen_4_ets",
        "fphen_5_ets",
        "t_Astart",
        "t_mid_anthesis",
        "t_fphen_1_ets",
        "t_fphen_3_ets",
        "t_fphen_4_ets",
        "t_fphen_5_ets",
        "t_lem",
        "t_lse",
        "t_lma",
        "t_lep",
        "lem",
        "lse",
        "lma",
        "lep",
    ])
    assert t_sgs is not None
    assert t_egs is not None
    season_td = t_egs - t_sgs

    t_Astart = season_td * f_Astart
    t_mid_anthesis = season_td * f_mid_anthesis
    t_fphen_1_ets = season_td * f_fphen_1_ets
    t_fphen_3_ets = season_td * f_fphen_3_ets
    t_fphen_4_ets = season_td * f_fphen_4_ets
    t_fphen_5_ets = season_td * f_fphen_5_ets
    t_lem = season_td * f_t_lem
    t_lse = season_td * f_t_lse
    t_lma = season_td * f_t_lma
    t_lep = season_td * f_t_lep

    sgs = next((dd for dd, tdd in zip(dd_list, td_list) if tdd > t_sgs))
    egs = next((dd for dd, tdd in zip(dd_list, td_list) if tdd > t_egs))
    Astart = next((dd for dd, tdd in zip(dd_list, td_list) if tdd > t_sgs + t_Astart)) - sgs
    mid_anthesis = next((dd for dd, tdd in zip(dd_list, td_list)
                         if tdd > t_sgs + t_mid_anthesis)) - sgs
    fphen_1_ets = next((dd for dd, tdd in zip(dd_list, td_list)
                        if tdd > t_sgs + t_Astart + t_fphen_1_ets)) - sgs - Astart
    fphen_3_ets = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >
                        t_sgs + t_Astart + t_fphen_1_ets + t_fphen_3_ets)) - sgs - Astart - fphen_1_ets
    fphen_4_ets = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >
                        t_sgs + t_Astart + t_fphen_1_ets + t_fphen_4_ets)) - sgs - Astart - fphen_1_ets
    fphen_5_ets = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >
                        t_sgs + t_Astart + t_fphen_1_ets + t_fphen_5_ets)) - sgs - Astart - fphen_1_ets

    lem = next((dd for dd, tdd in zip(dd_list, td_list) if tdd > t_sgs + t_lem)) - sgs
    lep = next((dd for dd, tdd in zip(dd_list, td_list) if tdd > t_sgs + t_lem + t_lep)) - sgs - lem
    lse = next((dd for dd, tdd in zip(dd_list, td_list) if tdd >
                t_sgs + t_lem + t_lep + t_lse)) - sgs - lem - lep
    lma = next((dd for dd, tdd in zip(dd_list, td_list) if tdd > t_sgs + t_lem + t_lma)) - sgs - lem

    return Output(
        SGS=sgs,
        EGS=egs,
        Astart=Astart,
        mid_anthesis=mid_anthesis,
        fphen_1_ets=fphen_1_ets,
        fphen_3_ets=fphen_3_ets,
        fphen_4_ets=fphen_4_ets,
        fphen_5_ets=fphen_5_ets,
        t_Astart=t_Astart,
        t_mid_anthesis=t_mid_anthesis,
        t_fphen_1_ets=t_fphen_1_ets,
        t_fphen_3_ets=t_fphen_3_ets,
        t_fphen_4_ets=t_fphen_4_ets,
        t_fphen_5_ets=t_fphen_5_ets,
        t_lem=t_lem,
        t_lse=t_lse,
        t_lma=t_lma,
        t_lep=t_lep,
        lem=lem,
        lse=lse,
        lma=lma,
        lep=lep,
    )


def estimate_growing_season_from_flag_leaf_period(
    Astart: TimeUnit,
    Aend: TimeUnit,
    f_leaf_f_fphen: Fraction = Wheat.f_fphen_1_ets + Wheat.f_fphen_5_ets,
) -> TimeUnit:
    growing_season_length = (Aend - Astart) / f_leaf_f_fphen
    return growing_season_length


def calculate_growing_season_from_leaf_f_phen_data(
    leaf_f_phen_data: List[float],
    time_data: List[TimeUnit],
    f_leaf_f_fphen: Fraction = Wheat.f_fphen_1_ets + Wheat.f_fphen_5_ets,
) -> Tuple[TimeUnit, TimeUnit, TimeUnit, TimeUnit, TimeUnit]:
    """Estimate the growing season start and end dates from leaf f phen data."""

    # Get the gradient of leaf_f_phen at each time point
    leaf_f_phen_diff_data = [0] + [a - b for a,
                                   b in zip(leaf_f_phen_data[1:], leaf_f_phen_data[:-1])]

    Astart = None
    Aend = None
    mature_leaf_start = None
    leaf_f_phen_mature_gradient = None
    senes_start = None
    for (
        leaf_f_phen,
        leaf_f_phen_diff,
        td,
        td_prev,
    ) in zip(
        leaf_f_phen_data,
        leaf_f_phen_diff_data,
        time_data,
        [0] + list(time_data),
    ):
        if leaf_f_phen_diff > 0 and Astart is None:
            # When leaf_f_phen first increases above 0 we assume this is the Astart
            Astart = td_prev
        if mature_leaf_start is None and Astart is not None and leaf_f_phen_diff < 0:
            mature_leaf_start = td_prev
            leaf_f_phen_mature_gradient = leaf_f_phen_diff
        if senes_start is None and mature_leaf_start is not None and not isclose(leaf_f_phen_diff - leaf_f_phen_mature_gradient, 0, abs_tol=1e-3):
            senes_start = td_prev
        if Aend is None and Astart is not None and leaf_f_phen_diff == 0 and leaf_f_phen == 0:
            Aend = td_prev
    if (Aend is None):
        warnings.warn("Leaf f phen data does not end with 0")
        Aend = time_data[-1]
    growing_season_length = estimate_growing_season_from_flag_leaf_period(
        Astart,
        Aend,
        f_leaf_f_fphen,
    )
    senes_start = senes_start or mature_leaf_start + (Aend - mature_leaf_start) / 2
    SGS = Aend - growing_season_length
    return SGS, Astart, Aend, mature_leaf_start, senes_start

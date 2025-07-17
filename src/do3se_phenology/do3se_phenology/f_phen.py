"""
F_phen functions.
=================

This module contains a number of functions for calculating f_phen and leaf_f_phen
from thermal time or julian day intervals.


"""

from math import isclose
from typing import List, NamedTuple
from collections import namedtuple
import numpy as np

from .utils import PLF_value
from .utils import offset
from .units import *


def get_fphen_PLF_ss(
    f_phen_1: int,
    f_phen_4: int,
    f_phen_a: float,
    f_phen_c: float,
    f_phen_e: float,
    SGS: int,
    EGS: int,
) -> PiecewiseFunction[TimeUnit, float]:
    """Get the piecewise function for fphen.

    Uses old UI method.
    """
    gs_values = [SGS, SGS, (SGS + f_phen_1), (EGS - f_phen_4), EGS, EGS]
    fphen_values = [0.0, f_phen_a, f_phen_c, f_phen_c, f_phen_e, 0.0]
    return gs_values, fphen_values


def get_fphen_PLF(
    f_phen_a: float,
    f_phen_b: float,
    f_phen_c: float,
    f_phen_d: float,
    f_phen_min: float,
    x_offset: float = 0,
) -> PiecewiseFunction[TimeUnit, float]:
    """Get the piecewise function for fphen.

    Assumes x = 0 is SGS.

    """
    x_values = [
        x_offset,
        f_phen_a + x_offset,
        f_phen_b + x_offset,
        f_phen_c + x_offset,
        f_phen_d + x_offset,
    ]
    y_values = [f_phen_min, f_phen_min, 1, 1, 0]
    return x_values, y_values


def get_leaf_fphen_PLF_ss(
    leaf_f_phen_1: int,
    leaf_f_phen_2: int,
    leaf_f_phen_a: float,
    leaf_f_phen_b: float,
    leaf_f_phen_c: float,
    Astart: int,
    Aend: int,
) -> PiecewiseFunction[TimeUnit, float]:
    """Get the piecewise function for leaf fphen.

    Uses old UI method.

    """

    # x values
    gs_values = [
        Astart,
        Astart,
        (Astart + leaf_f_phen_1),
        (Aend - leaf_f_phen_2),
        Aend,
        Aend,
    ]
    # y values
    fphen_values = [
        0.0,
        leaf_f_phen_a,
        leaf_f_phen_b,
        leaf_f_phen_b,
        leaf_f_phen_c,
        0.0,
    ]
    return gs_values, fphen_values


def get_leaf_fphen_PLF(
    leaf_f_phen_a: float,
    leaf_f_phen_b: float,
    leaf_f_phen_e: float,
    leaf_f_phen_g: float,
    leaf_f_phen_h: float,
    leaf_f_phen_i: float,
    Astart: float,
) -> PiecewiseFunction[TimeUnit, float]:
    """Get the piecewise function for leaf fphen.

    Assumes x = 0 is SGS.
    """
    assert isclose(leaf_f_phen_a + leaf_f_phen_b, 1)
    x_values = [
        0,  # SGS
        Astart,
        Astart,
        Astart + leaf_f_phen_e,
        Astart + leaf_f_phen_e + leaf_f_phen_g,
        Astart + leaf_f_phen_e + leaf_f_phen_h,
        Astart + leaf_f_phen_e + leaf_f_phen_i,
        Astart + leaf_f_phen_e + leaf_f_phen_i,  # EGS
    ]
    y_values = [
        0,
        0,
        1,
        1,
        1,
        leaf_f_phen_b,
        0,
        0,
    ]

    return x_values, y_values


def f_phen_simple_PLF(
    f_phen_1: int,
    f_phen_4: int,
    f_phen_a: float,
    f_phen_c: float,
    f_phen_e: float,
    SGS: int,
    EGS: int,
    dd: int,
) -> float:
    r"""Phenology effect on stomatal conductance.

    according to a simple piecewise linear function.

    To handle all situations, including winter growing seasons, everything is
    re-indexed to be relative to SGS=0.


    .. code-block:: python


                                    anthesis
                               __________|___ <f_phen_c>
                              /              \
                             /                \
             <f_phen_a>  _  /                  \
                           |                    \  _<f_phen_e>
                           |                     |
        f_phen_min  |______|                     |________ 0

                    |Sowing day or DAY 0
                    |------| SGS
                           |---| f_phen_1
                                           |-----| f_phen_4
                    |----------------------------| EGS

    Parameters
    ----------
    f_phen_1: int
        Time from f_phen_a to f_phen_b [days]
    f_phen_4: int
        Time from f_phen_d to f_phen_e [days]
    f_phen_a: float
        f_phen at SGS [fraction]
    f_phen_c: float
        max f_phen [fraction]
    f_phen_e: float
        f_phen at EGS [fraction]
    SGS : int
        Start of growing season [julian day]
    EGS : int
        End of growing season [julian day]
    dd : int
        days since Jan 1st [julian day]

    Returns
    -------
    f_phen: float
        [description]

    """
    gs_values = [SGS, SGS, (SGS + f_phen_1), (EGS - f_phen_4), EGS, EGS]
    fphen_values = [0.0, f_phen_a, f_phen_c, f_phen_c, f_phen_e, 0.0]

    # Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    gs_offset = offset(gs_values, float(SGS), 365.0)

    gs_values_in_size_order = all([a <= b for a, b in zip(gs_offset[0: 5], gs_offset[1: 6])])
    if not gs_values_in_size_order:
        raise ValueError("f_phen_simple_PLF: points not in order")

    dd_adj = dd - SGS if SGS - dd <= 0 else dd - SGS + 365

    # Lookup value in PLF
    func = [gs_offset, fphen_values]
    f_phen = PLF_value(func, float(dd_adj))
    return f_phen


def tt_f_phen_simple_PLF(
    td: float,
    t_f_phen_a: float,
    t_f_phen_b: float,
    t_f_phen_c: float,
    t_f_phen_d: float,
    f_phen_min: float,
    td_at_sgs: float,
) -> float:
    r"""Calculate phenology effect fraction on stomatal conductance using thermal time intervals.

    .. code-block:: python

                               ______________ 1
                              /              \
                             /                \
                            /                  \
                           /                    \
        f_phen_min   _____/                      \
                                                  \
                                                   \________  0

                    |Sowing day
                    |------| t_f_phen_a
                    |----------| t_f_phen_b
                    |-----------------------| t_f_phen_c
                    |-------------------------------| t_f_phen_d


    Parameters
    ----------
    td: float
        Thermal time [degC days]
    t_f_phen_a: float
        thermal time between sowing day and f_phen_a [degC days]
    t_f_phen_b: float
        thermal time between sowing day and f_phen_b [degC days]
    t_f_phen_c: float
        thermal time between sowing day and f_phen_c [degC days]
    t_f_phen_d: float
        thermal time between sowing day and f_phen_d [degC days]
    f_phen_min: float
        thermal time between sowing day and phen_min [fraction]

    Returns
    -------
    f_phen: float
        [description]

    """
    td_s = td - td_at_sgs  # td between sgs and current td
    if td_s < t_f_phen_a:
        return f_phen_min
    if td_s < t_f_phen_b:
        return (1 - f_phen_min) * ((td_s - t_f_phen_a) / (t_f_phen_b - t_f_phen_a)) + f_phen_min
    if td_s < t_f_phen_c:
        return 1
    if td_s < t_f_phen_d:
        return (1 - f_phen_min) * (t_f_phen_d - td_s) / (t_f_phen_d - t_f_phen_c) + f_phen_min
    if td_s >= t_f_phen_d:
        return 0
    return


def tt_f_phen_simple_PLF_range(
    td_list: List[float],
    t_f_phen_a: float,
    t_f_phen_b: float,
    t_f_phen_c: float,
    t_f_phen_d: float,
    f_phen_min: float,
    td_at_sgs: float,
) -> List[float]:
    return [
        tt_f_phen_simple_PLF(
            td,
            t_f_phen_a,
            t_f_phen_b,
            t_f_phen_c,
            t_f_phen_d,
            f_phen_min,
            td_at_sgs,
        ) for td in td_list
    ]


def f_phen_complex_PLF(
    f_phen_1: int,
    f_phen_2: int,
    f_phen_3: int,
    f_phen_4: int,
    f_phen_limA: int,
    f_phen_limB: int,
    f_phen_a: float,
    f_phen_b: float,
    f_phen_c: float,
    f_phen_d: float,
    f_phen_e: float,
    SGS: int,
    EGS: int,
    dd: int,
) -> float:
    r"""Phenology effect on stomatal conductance.

    According to a more complicated piecewise linear function.

    To handle all situations, including winter growing seasons, everything is
    re-indexed to be relative to SGS=0.

    If not f_phen_limA and f_phen_limB then ignore between limA and limB

    .. code-block:: python

                        f_phen_b                  f_phen_d
                       __________                __________
                      /         :\              /:         \\
                     /          : \  f_phen_c  / :          \\
                    /           :  \__________/  :           \\
                   /            :<2>          <3>:            \\
        f_phen_a  /             :                :             \\
                 | <1>          :                :          <4> \  f_phen_e
                 |              :                :               |
          _______|              :                :               |________
                SGS            limA            limB             EGS

    Parameters
    ----------
    f_phen_1: int
        Time from f_phen_a to f_phen_b [days]
    f_phen_2: int
        Time from f_phen_b to f_phen_c [days]
    f_phen_3: int
        Time from f_phen_c to f_phen_d [days]
    f_phen_4: int
        Time from f_phen_d to f_phen_e [days]
    f_phen_limA: int
        Start of soil water limitation [days]
    f_phen_limB: int
        End of soil water limitation [days]
    f_phen_a: float
        f_phen at SGS [?]
    f_phen_b: float
        f_phen at SGS [?]
    f_phen_c: float
        f_phen at SGS [?]
    f_phen_d: float
        f_phen at SGS [?]
    f_phen_e: float
        f_phen at SGS [?]
    SGS : int
        Start of growing season [day of year]
    EGS : int
        End of growing season [day of year]
    dd : int
        Day of year [day of year]

    Returns
    -------
    f_phen: float
        [description]

    """
    use_complex_middle_components = f_phen_limA and f_phen_limB

    complex_middle_gs_values = [
        f_phen_limA,
        f_phen_limA + f_phen_2,
        f_phen_limB - f_phen_3,
        f_phen_limB
    ] if use_complex_middle_components else [
        SGS + f_phen_1,
        SGS + f_phen_1,
        EGS - f_phen_4,
        EGS - f_phen_4,
    ]

    complex_middle_fphen_values = [
        f_phen_b,
        f_phen_c,
        f_phen_c,
        f_phen_d,
    ] if use_complex_middle_components else [
        f_phen_b, f_phen_b, f_phen_d, f_phen_d
    ]

    gs_values = [
        SGS,
        SGS,
        SGS + f_phen_1,
        *complex_middle_gs_values,
        EGS - f_phen_4,
        EGS,
        EGS,
    ]
    fphen_values = [
        0.0,
        f_phen_a,
        f_phen_b,
        *complex_middle_fphen_values,
        f_phen_d,
        f_phen_e,
        0.0,
    ]

    # Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    gs_offset = offset(gs_values, float(SGS), 365.0)
    gs_values_are_in_size_order: bool = all(
        [a <= b for a, b in zip(gs_offset[0: 9], gs_offset[1: 10])])

    if not gs_values_are_in_size_order:
        print(gs_values)
        raise ValueError("f_phen_simple_PLF: points not in order")

    dd_adj = dd - SGS if SGS - dd <= 0 else dd - SGS + 365

    # Lookup value in PLF
    func = [gs_offset, fphen_values]

    f_phen = PLF_value(func, float(dd_adj))
    return f_phen


def leaf_f_phen_PLF(
    leaf_f_phen_1: int,
    leaf_f_phen_2: int,
    leaf_f_phen_a: float,
    leaf_f_phen_b: float,
    leaf_f_phen_c: float,
    Astart: int,
    Aend: int,
    dd: int
) -> float:
    r"""Phenology effect on leaf stomatal conductance.

    According to a simple piecewise linear function.

    To handle all situations, including winter growing seasons, everything is
    re-indexed to be relative to Astart=0.

    .. code-block:: python

                        b
                    ______________
                  /              \\
                 /                \\
             a  /                  \\
                |<1>              <2>\  c
                |                     |
         _______|                     |________
            Astart                  Aend
        ~~~


    Parameters
    ----------
    leaf_f_phen_1: int
        Time from _a to _b [days]
    leaf_f_phen_2: int
        Time from _b to _c [days]
    leaf_f_phen_a: float
        f_phen at Astart [?]
    leaf_f_phen_b: float
        f_phen at mid-season peak [?]
    leaf_f_phen_c: float
        f_phen at Aend [?]
    Astart : int
        Start of accumulation period
    Aend : int
        End of accumulation period
    dd : int
        Day of year [day of year]

    Returns
    -------
    leaf_f_phen: float
        [description]

    """
    # x values
    gs_values = [
        Astart,
        Astart,
        (Astart + leaf_f_phen_1),
        (Aend - leaf_f_phen_2),
        Aend,
        Aend,
    ]
    # y values
    fphen_values = [
        0.0,
        leaf_f_phen_a,
        leaf_f_phen_b,
        leaf_f_phen_b,
        leaf_f_phen_c,
        0.0,
    ]

    # Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    gs_offset = offset(gs_values, float(Astart), 365.0)

    gs_values_in_size_order = all([a <= b for a, b in zip(gs_offset[0: 5], gs_offset[1: 6])])
    if not gs_values_in_size_order:
        raise ValueError("f_phen_simple_PLF: points not in order")

    dd_adj = dd - Astart if Astart - dd <= 0 else dd - Astart + 365

    # Lookup value in PLF
    func = [gs_offset, fphen_values]
    leaf_f_phen = PLF_value(func, float(dd_adj))
    return leaf_f_phen


def tt_leaf_f_phen_PLF(
    td: float,
    t_leaf_f_phen_a: float,
    t_leaf_f_phen_b: float,
    t_leaf_f_phen_e: float,
    t_leaf_f_phen_g: float,
    t_leaf_f_phen_h: float,
    t_leaf_f_phen_i: float,
    t_astart: float,
    td_at_sgs: float,
) -> float:
    r"""Calculate tt_leaf_f_phen_plf.

    .. code-block:: python

                Mid Anthesis    f
                 ____|__________/
                |              \\
                |               \\
             a  |                \\
                |<1>              <2>\  c
                |                     |
         _______|                     |________
            Astart                  Aend


        |-------| t_astart
                |----| t_leaf_f_phen_e (fphen_1_ETS)
                     |----------| t_leaf_f_phen_g (fphen_3_ETS)
                     |-------------| t_leaf_f_phen_h (fphen_4_ETS)
                     |----------------| t_leaf_f_phen_i (fphen_5_ETS)

    Parameters
    ----------
    td: float
        Thermal time [degC days]
    t_leaf_f_phen_a : float
        Gradient of descent during seed setting [fraction]
    t_leaf_f_phen_b : float
        Gradient of descent during senescence [fraction]
    t_leaf_f_phen_e : float
        Thermal time between Astart and leaf_f_phen_e(Mid anthesis) [DegC Days]
    t_leaf_f_phen_g : float
        Thermal time between Mid Anthesis and start of seed setting [DegC Days]
    t_leaf_f_phen_h : float
        Thermal time between Mid anthesis and start of senesence [DegC Days]
    t_leaf_f_phen_i : float
        Thermal time between Mid anthesis and end of senescence [DegC Days]
    t_astart : float
        Thermal time between Sowing day and Astart [DegC Days]

    Returns
    -------
    [type]
        [description]

    """
    td_a = td - t_leaf_f_phen_e - t_astart - td_at_sgs  # thermal time at Mid anthesis
    if td_a < -t_leaf_f_phen_e:
        return 0
    if td_a <= t_leaf_f_phen_g:
        return 1
    if td_a > t_leaf_f_phen_g and td_a <= t_leaf_f_phen_h:
        return 1 - t_leaf_f_phen_a / (t_leaf_f_phen_h - t_leaf_f_phen_g) * (td_a - t_leaf_f_phen_g)
    if td_a > t_leaf_f_phen_h and td_a <= t_leaf_f_phen_i:
        return t_leaf_f_phen_b - t_leaf_f_phen_b / (t_leaf_f_phen_i - t_leaf_f_phen_h) * (td_a - t_leaf_f_phen_h)
    else:
        return 0


def calc_leaf_f_phen_effect_on_V_cmax_25(
    V_cmax_25_in: float,
    J_max_25_in: float,
    leaf_f_phen: float,
) -> NamedTuple:
    """Calculate the effect of leaf_f_phen on V_cmax_25 and J_max_25.

    Parameters
    ----------
    V_cmax_25_in : float
        Maximum catalytic rate at 25 degrees [umol m-2 s-1]
    J_max_25_in : float
        Maximum rate of electron transport at 25 degrees [umol m-2 s-1]
    leaf_f_phen : float
        Phenology-related effect on leaf gsto [fraction]

    Returns (As namedtuple)
    -------
    V_cmax_25
        Maximum catalytic rate at 25 degrees [umol m-2 s-1]
    J_max_25
        Maximum rate of electron transport at 25 degrees [umol m-2 s-1]

    """
    Output = namedtuple('Output', 'V_cmax_25 J_max_25')
    # Added min values here to avoid 0 division in photosynthesis module
    MIN_V_cmax_25 = 0.000001
    MIN_J_max_25 = 0.000001
    V_cmax_25 = max(MIN_V_cmax_25, V_cmax_25_in * leaf_f_phen)
    J_max_25 = max(MIN_J_max_25, J_max_25_in * leaf_f_phen)
    return Output(
        V_cmax_25=V_cmax_25,
        J_max_25=J_max_25,
    )


def tt_leaf_f_phen_PLF_range(
    td_list: List[float],
    t_leaf_f_phen_a: float,
    t_leaf_f_phen_b: float,
    t_leaf_f_phen_e: float,
    t_leaf_f_phen_g: float,
    t_leaf_f_phen_h: float,
    t_leaf_f_phen_i: float,
    t_astart: float,
    td_at_sgs: float,
) -> List[float]:
    return [tt_leaf_f_phen_PLF(
        td,
        t_leaf_f_phen_a,
        t_leaf_f_phen_b,
        t_leaf_f_phen_e,
        t_leaf_f_phen_g,
        t_leaf_f_phen_h,
        t_leaf_f_phen_i,
        t_astart,
        td_at_sgs,
    ) for td in td_list]

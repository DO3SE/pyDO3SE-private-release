"""Helper functions associated with the met module."""

from math import exp
from typing import List, NamedTuple
from collections import namedtuple
import warnings

def saturated_vapour_pressure(Ts_C: float) -> float:
    """Calculate saturated vapour pressure [kPa].

    Parameters
    ----------
    Ts_C: float
        Surface air temperature [degrees C]

    """
    return 0.611 * exp(17.27 * Ts_C / (Ts_C + 237.3))


def calc_humidity(
    Ts_C_in: float,
    VPD_in: float,
    RH_in: float,
) -> float:
    """Calculate humidity.

    met_humidity in DO3SE fortran model

    Parameters
    ----------
    Ts_C_in: float
        [Description][Unit]
    VPD_in: float
        [Description][Unit]
    RH_in: float
        [Description][Unit]
    """
    Output = namedtuple('output', 'esat RH eact VPD')

    esat = saturated_vapour_pressure(Ts_C_in)
    if VPD_in is None and RH_in is None:
        raise Exception('Must supply VPD or RH')
    elif RH_in is None:
        # Calculate relative humidity from VPD
        eact = esat - VPD_in
        RH = eact / esat
        VPD = VPD_in
    elif VPD_in is None:
        # Calculate VPD from relative humidity
        eact = esat * RH_in
        VPD = esat - eact
        RH = RH_in
    else:
        # TODO: Check if we should raise error here
        raise Exception("Supplied both VPD and RH")

    return Output(
        esat=esat,
        RH=RH,
        eact=eact,
        VPD=VPD,
    )


def calc_humidity_list(
    Ts_C_list: List[float],
    VPD_list: List[float],
    RH_list: List[float],
) -> NamedTuple:
    """Run Calc_humidity on a list of data.

    Parameters
    ----------
    Ts_C_list : List[float]
        [description]
    VPD_list : List[float]
        [description]
    RH_list : List[float]
        Relative humidity between 0 and 1

    Returns
    -------
    esat: List[float]
        [DESCRIPTION][Unit]
    RH: List[float]
        [DESCRIPTION][Unit]
    eact: List[float]
        [DESCRIPTION][Unit]
    VPD: List[float]
        [DESCRIPTION][Unit]
    """
    Output = namedtuple('Output', 'esat, RH, eact, VPD')
    esat_out = []
    RH_out = []
    eact_out = []
    VPD_out = []

    if RH_list is not None and len(RH_list) > 0 and RH_list[0] is not None and max(RH_list) > 1:
        warnings.warn("Rh should be between 0 and 1")
        # raise ValueError("Rh should be between 0 and 1")

    for Ts_C, VPD, RH in zip(Ts_C_list, VPD_list, RH_list):
        out = calc_humidity(Ts_C, VPD, RH)
        esat_out.append(out.esat)
        RH_out.append(out.RH)
        eact_out.append(out.eact)
        VPD_out.append(out.VPD)
    return Output(
        esat=esat_out,
        RH=RH_out,
        eact=eact_out,
        VPD=VPD_out,
    )


def calc_vpd_daily_accumulation(
    VPD: float,
    VPD_dd_in: float,
    is_daylight: bool,
    hr: int,
) -> float:
    """Calculate the daylight accumulated Vapour pressure deficit.

    during daylight hours per day

    Parameters
    ----------
    VPD: float
        [description][unit]
    VPD_dd_in: float
        [description][unit]
    is_daylight: bool
        [description][unit]
    hr: int
        [description][unit]

    Returns
    -------
    VPD_dd: float
        [description][unit]
    """
    # ! Reset VPD sum at start of day
    VPD_dd_prev = VPD_dd_in if hr > 0 else 0
    # ! Only accumulate VPD sum during daylight hours
    VPD_dd = VPD_dd_prev + VPD if is_daylight else VPD_dd_prev
    return VPD_dd


def calc_vpd_daily_accumulation_list(
    VPD_list: List[float],
    is_daylight_list: List[bool],
    hr_list: List[int],
) -> List[float]:
    """Calculate the accumulated Vapour pressure deficit for full data.

    Parameters
    ----------
    VPD_list : List[float]
        [description]
    is_daylight_list : List[bool]
        [description]
    hr_list : List[int]
        [description]

    Returns
    -------
    List[float]
        [description]
    """
    vpd_dd = 0
    vpd_dd_out = []
    for VPD, is_daylight, hr in zip(VPD_list, is_daylight_list, hr_list):
        vpd_dd = calc_vpd_daily_accumulation(VPD, vpd_dd, is_daylight, hr)
        vpd_dd_out.append(vpd_dd)
    return vpd_dd_out


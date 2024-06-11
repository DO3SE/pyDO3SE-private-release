""" Vernalisation functions.

The vernalisation factor reduces the accumulated thermal time between emergence and Anthesis.
It is effected by the daily min and maximum temperature.



References
----------
 - The APSIM-Wheat Module (7.5 R3008)

"""
from typing import List, Tuple

from do3se_phenology.state import PhenologyStage
from do3se_met.irradiance import calc_photoperiod, calc_photoperiod_factor

def calculate_vernalisation_factor(
    # state
    max_ambient_temp: float,
    min_ambient_temp: float,
    V_acc_prev: float,
    phenology_stage: PhenologyStage,

    # config
    v_T_max: float,
    v_T_min: float,
    PIV: float,
) -> Tuple[float, float, float, float]:
    """Calculate the vernalisation factor for single data point.

    Vernalisation only applies between Emergence and Floral initiation.

    Parameters
    ----------
    max_ambient_temp : float
        Maximum ambient temperataure [DegC]
    min_ambient_temp : float
        Minimum ambient temperataure [DegC]
    V_acc_prev: float,
        Previous accumulated Vernalisation
    v_T_max : float
        cultivar specific parameter
    v_T_min : float
        cultivar specific parameter
    PIV : float
        cultivar specific parameter

    Returns
    -------
    Tuple[float, float]
        Vernalisation Factor, Accumulated Vernalisation

    """

    # TODO: Astart should be flowering.
    if not(PhenologyStage.EMERGED <= phenology_stage < PhenologyStage.ASTART):
        return 1, 0, 0, 0

    max_crown_temp = 2 + (max_ambient_temp) * (0.4 + 0.0018 * (0 - 15) **
                                               2) if max_ambient_temp < 0 else max_ambient_temp
    min_crown_temp = 2 + (min_ambient_temp) * (0.4 + 0.0018 * (0 - 15) **
                                               2) if min_ambient_temp < 0 else min_ambient_temp

    t_leaf = (max_crown_temp + min_crown_temp) / 2  # Crown temp

    a = 1.4 - 0.0778 * t_leaf
    b = 0.5 + 13.44 * t_leaf / (max_ambient_temp - min_ambient_temp + 3)**2.0
    Vc = min(a, b)

    V_pos = Vc if min_ambient_temp < v_T_min and max_ambient_temp < v_T_max else 0

    # TODO: Check what the min value should be here
    Vd = min(0.5 * (max_ambient_temp - v_T_max), V_acc_prev)

    V_neg = Vd if max_ambient_temp > v_T_max and V_acc_prev > 10 else 0
    V_tot = V_pos + V_neg
    V_acc = V_acc_prev + V_tot

    Vf = min(1, max(0, 1 - (0.0054545 * PIV + 0.0003) * (50 - V_acc)))
    return Vf, V_acc, V_pos, V_neg


def calc_vernalised_thermal_time_range(
    td_data: List[float],
    hr_data: List[int],
    dd_data: List[int],
    max_ambient_temperatures: List[float],
    min_ambient_temperatures: List[float],
    t_emerg: float,
    t_flower: float,
    v_T_max: float=30,
    v_T_min: float=15,
    PIV: float=1.5,
    PID: float=None,
    lat: float=None,
) -> List[float]:
    """Calculate the vernalised thermal time for a range of data.

    Parameters
    ----------
    td_data : List[float]
        Thermal time data
    hr_data : List[int]
        Hour data
    dd_data : List[int]
        Julian day data
    max_ambient_temperatures : List[float]
        Max ambient daily temperatures
    min_ambient_temperatures : List[float]
        Min ambient daily temperature
    t_emerg : float
        thermal time at emergence
    t_flower : float
        thermal time at flowering
    v_T_max : int, optional
        Max vernalisation temp, by default 30
    v_T_min : int, optional
         min vernalisation temp, by default 15
    PIV : float, optional
        vernalisation sensativity factor, by default 1.5
    PID : float, optional
        photoperiod sensativity factor, by default None
    lat: float
        latitude, only required for photoperiod calcs

    Returns
    -------
    List[float]
        Vernalised thermal time values

    """
    Vf = 0
    td_f_values = [0]
    prev_td = 0

    for td, dd, hr, max_ambient_temp, min_ambient_temp in zip(
        td_data, dd_data, hr_data, max_ambient_temperatures, min_ambient_temperatures,
    ):
        if hr == 0:
            emerged = td_f_values[-1] > t_emerg
            flowering = td_f_values[-1] > t_flower
            phenology_stage = PhenologyStage.SOWN if not emerged else PhenologyStage.EMERGED if not flowering else PhenologyStage.ASTART
            photoperiod = calc_photoperiod(dd, lat) if PID is not None else None
            photoperiod_factor = calc_photoperiod_factor(photoperiod, PID) if PID is not None else 1
            Vf, V_acc, V_pos, V_neg = calculate_vernalisation_factor(
                max_ambient_temp=max_ambient_temp,
                min_ambient_temp=min_ambient_temp,
                V_acc_prev=0,
                phenology_stage=phenology_stage,
                v_T_max=v_T_max,
                v_T_min=v_T_min,
                PIV=PIV,
            )
        td_diff = td - prev_td

        td_f = td_f_values[-1] + td_diff * Vf * photoperiod_factor
        td_f_values.append(td_f)
        prev_td = td

    td_f_values = td_f_values[1:]
    return td_f_values

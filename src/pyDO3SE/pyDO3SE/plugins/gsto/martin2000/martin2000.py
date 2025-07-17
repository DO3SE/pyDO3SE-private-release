"""Functions derived from Martin 2000."""
from collections import namedtuple
from typing import NamedTuple
from deprecated import deprecated


@deprecated(version='4.0.0', reason="No longer supported")
def calc_O3_effect_on_V_xmax_25(
    K_z: float,
    FO3_eff: float,
    V_cmax_25_in: float,
    J_max_25_in: float,
) -> NamedTuple:
    """Calculate the effect of O3 on V_cmax_25 and J_max_25.

    Parameters
    ----------
    K_z : float
        Coefficient for ozone damage [dimensionless]
    FO3_eff : float
        (Accumulated) effective ozone dose [nmol O3 m-2 PLA]
    V_cmax_25_in : float
        Maximum catalytic rate at 25 degrees [umol m-2 s-1]
    J_max_25_in : float
        Maximum rate of electron transport at 25 degrees [umol m-2 s-1]

    Returns
    -------
    V_cmax_25
        Maximum catalytic rate at 25 degrees [umol m-2 s-1]
    J_max_25
        Maximum rate of electron transport at 25 degrees [umol m-2 s-1]
    """
    Output = namedtuple('Output', 'V_cmax_25 J_max_25')
    # Percentage reduction in V_cmax (converting FO3_eff from nmol to mmol)
    delta_V_cmax = K_z * (FO3_eff / 1e6)
    # Convert to multiplier
    mult_V_cmax = 1.0 - (delta_V_cmax / 100)
    # Reduce V_cmax and J_max
    V_cmax_25 = max(0.1 * V_cmax_25_in, min(V_cmax_25_in * mult_V_cmax, V_cmax_25_in))
    J_max_25 = max(0.1 * J_max_25_in, min(J_max_25_in * mult_V_cmax, J_max_25_in))
    return Output(
        V_cmax_25=V_cmax_25,
        J_max_25=J_max_25,
    )

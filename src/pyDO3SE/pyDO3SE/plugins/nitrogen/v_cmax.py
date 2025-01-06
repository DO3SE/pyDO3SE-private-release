from typing import List
from math import exp


def multilayer_vcmax25(
    layer_LAI: List[float],
    kN: float = 0.78,
    V_cmax_25: float = 90,
) -> List[float]:
    """Get the v_cmax_25 value at each layer using nitrogen relationship coefficients.

    Parameters
    ----------
    layer_LAI : List[float]
        per layer lai
    kN: float
        Nitrogen profile coefficient [-]
    Vcmax25: float
        Max carboxylation capacity [umol m2/s]

    Returns
    -------
    List[float]
        v_cmax_25i per layer v_cmax_25 [umol m2/s]


    References
    ----------
    - Williams et al 2017 https://gmd.copernicus.org/articles/10/1291/2017/gmd-10-1291-2017.pdf
    - Clark et al (2011) https://gmd.copernicus.org/articles/4/701/2011/gmd-4-701-2011.pdf

    """

    lai_c = [0.0] + [sum(layer_LAI[0:i + 1]) for i, _ in enumerate(layer_LAI)]
    total_lai = sum(layer_LAI)
    if total_lai == 0:
        return [V_cmax_25 for _ in layer_LAI]
    v_cmax_25i = [V_cmax_25 * exp(-kN * lai / sum(layer_LAI)) for lai in lai_c[:-1]]
    return v_cmax_25i

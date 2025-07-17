from typing import List
from math import exp


def calc_multilayer_nitrogen(
    layer_LAI: List[float],
    Nb: float = 0.05,
    k: float = 0.5,
    yp: float = 1.0,
) -> List[float]:
    """Get the nitrogen at each level using linear coefficients

    Parameters
    ----------
    layer_LAI : List[float]
        per layer lai
    Nb: float
        Leaf N concentration not associated with Anet [mmol/m3]
    k: float
        Canopy extinction co-efficient [m2 ground/m2 leaf]
    yp: float
        dimensionless co-efficient [-]
    Vcmax25: float
        Max carboxylation capacity [umol m2/s]

    Returns
    -------
    List[float]
        Nitrogen per layer [umol m2/s]


    References
    ----------
    - Williams et al 2017 https://gmd.copernicus.org/articles/10/1291/2017/gmd-10-1291-2017.pdf
    - Clark et al (2011) https://gmd.copernicus.org/articles/4/701/2011/gmd-4-701-2011.pdf

    """
    # TODO: Calculate N0
    # Leaf N concentration at the top of the canopy [mmol/m2]
    N0: float = 0.3
    lai_c = [0.0] + [sum(layer_LAI[0:i + 1]) for i, _ in enumerate(layer_LAI)]
    Ni = [N0 - (N0 - Nb) * (1 - exp(-k * lai)) ** yp for lai in lai_c]
    return Ni[:-1]

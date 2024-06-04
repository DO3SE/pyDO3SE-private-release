"""Unit conversions used by carbon allocation module.

TODO: Move these to met package.
"""


def umol_c_to_kg_c(val: float) -> float:
    """Convert from umols per s of carbon to kg per hour of carbon.

    Parameters
    ----------
    val : float
        input

    Returns
    -------
    float
        output
    """
    UMOL_TO_MOL = 1e-6
    MOL_TO_G = 12
    S_TO_H = 3600
    G_TO_KG = 1e-3
    return val * UMOL_TO_MOL * MOL_TO_G * S_TO_H * G_TO_KG

def deg_to_kel(v: float) -> float:
    return 273.15 + v


def convert_ppb_to_nmol(
    value_ppb: float,
    Ts_C: float,
    P: float,
    gas_molecular_weight: float,
    C: float = 20.833,  # TODO: What is this??
) -> float:
    """Convert a gas from ppb to nmol.

    # TODO: Check this is correct

    Uses ideal gas law https://en.wikipedia.org/wiki/Ideal_gas_law

    ## Ideal gas law
    PV = n * R * T

    Parameters
    ----------
    Ts_c : float
        air temperature [degrees]
    P : float
        Pressure [kPa]
    value_ppb : float
        value to convert [ppb]
    gas_molecular_weight: float,
        molecular weight of the gas to convert [g]
    C: float
        Constant???

    Returns
    -------
    float
        gas [umol/m^3] # Check output
    """
    R = 8.314510  # Ideal gas constant
    T = deg_to_kel(Ts_C)

    # Check below volume equation
    V = value_ppb * gas_molecular_weight * C  # Conversion of ppb to volume [m^2]

    # Check this is correct output
    n = (V * P) / (R * T)
    return n

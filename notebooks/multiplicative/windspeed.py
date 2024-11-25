
# %%
def calc_windspeed_parameters(
        h: float,
        u: float,
        o_top_chamber: bool = True,
        loc_h_u: float = None,
        loc_z_u: float = None,
        MIN_WINDSPEED: float = 0.01,
        MIN_USTAR: float = 0.0001,
) -> NamedTuple:
    """Calculate wind speed parameters.

    aka met_windspeed

    Parameters
    ----------
    h: float
        Canopy height [m]
    u: float
        wind speed from met data [m/s]
    o_top_chamber: bool
        OTC: Is open top chamber experiment
    loc_h_u: float
        canopy height for windspeed measurement, default = target canopy[m]
    loc_z_u: float
        Measurement height for windspeed[m]
    MIN_WINDSPEED: float
        Minimum windspeed value [m/s]
    MIN_USTAR: float
        Minimum friction velocity value [m/s]

    Returns
    -------
    Named tuple containing:
    u_50: float
        [description][Unit]
    ustar: float
        [description][Unit]
    micro_u: float
        [description][Unit]
    """
    Output = namedtuple('Output', 'u_50 ustar micro_u')

    h_u = h if o_top_chamber else loc_h_u if loc_h_u is not None else h
    z_u = h if o_top_chamber else loc_z_u if loc_z_u is not None else h

    u_d = h_u * CANOPY_D
    u_z0 = h_u * CANOPY_Z0
    d = h * CANOPY_D
    z0 = h * CANOPY_Z0

    # ! Find ustar over reference canopy
    ustar_ref = ustar_from_velocity(max(MIN_WINDSPEED, u), (z_u - u_d), u_z0)
    # Find windspeed at izR, over reference canopy
    u_50 = max(MIN_WINDSPEED, velocity_from_ustar(ustar_ref, (izR - u_d), u_z0))
    # Find ustar over target canopy, assuming that at izR windspeed will
    # be equal over both vegitations
    ustar = max(MIN_USTAR, ustar_from_velocity(u_50, (izR - d), z0))
    # Calculates windspeed at top layer canopy (uh)
    micro_u = max(MIN_WINDSPEED, velocity_from_ustar(ustar, (h - d), z0))

    return Output(
        u_50=u_50,
        ustar=ustar,
        micro_u=micro_u,
    )

"""Met module windspeed helpers.

Many of these are cross checked against EMEP calculations through discussions with Dave Simpson (2021).

## Ustar

To correctly calculate ustar we need to perform a number of iterations between the
values for ustar_ref and invL till they converge. There are 3 method to achieve this:

 - 1) Iteratively calculate ustar_ref and invL until they converge on the correct value.
 - 2) Input the grid ustar value and calculate invL from this.
 - 3) Use the previous hours ustar_ref value to calculate invL for the current hour.

In this model we use method 3.

"""

from math import exp, log, pi, sqrt, isclose
from typing import List, NamedTuple, Tuple
from collections import namedtuple
from decimal import Decimal as D
from deprecated import deprecated

from do3se_met.model_constants import  MIN_WINDSPEED, izR
from do3se_met.physical_constants import VON_KAR as K, Rmass, cp, g
from do3se_met.resistance import calc_PsiM

def calc_monin_obukhov_length(
    Tk: float,
    ustar: float,
    Hd: float,
    P: float,
) -> float:
    """Calculate the Monin-Obukhov Length.

    NOTE: Sometimes refered to as invL which is 1/L.

    # rho_surf = psurf/(RGAS_KG * Sub(iL)%t2 )
    # invL =  -KARMAN * GRAV * Hd / ( CP * rho_surf * ustar**3 * t2)

    Parameters
    ----------
    Tk : float
        Temperature [Kelvin]
    ustar : float
        Friction velocity [m/s]
    Hd : float
        Sensible heat flux [W/m^2] Note should be +ve in middle of summer day
    P : float
        Atmospheric Pressure [kPa]

    Returns
    -------
    float
        Monin-Obukhov Length (Used in ustar calculations)

    """

    # TODO: We should check if this should be negative on input
    Hd_f = Hd # max(Hd, 0.000000000001) # Limit turned off to match EMEP
    # Surface density of dry air (including conversion from to kPa to hPa??)
    rho = (P * 1000) / (Rmass * Tk)
    # Monin-Obukhov Length
    L = -(Tk * ustar**3 * rho * cp) / (K * g * (-Hd_f))

    return L


def calc_ustar_and_L(
    u: float,
    Hd: float,
    P: float,
    Tk: float,
    z_u: float,
    u_d: float,
    u_z0: float,
    initial_L: float = 1,
    ustar_ref_in: float = None,
    min_windspeed: float = MIN_WINDSPEED,
    min_ustar: float = 0.0001,
    MAX_ITERATIONS: int = 1,
    TOLERANCE: float = 1e-3,
) -> Tuple[float, float]:
    """Calculate ustar at reference height.

    Parameters
    ----------
    u : float
        Wind speed [m/s]
    Hd : float
        Sensible heat flux [W/m^2]
    P : float
        Air pressure [kPa]
    Tk : float
        Temperature [Kelvin]
    z_u : float
        Wind speed measured height [m]
    u_d : float
        Measured canopy displacement height [m]
    u_z0 : float
        Measured canopy roughness length [m]
    initial_L : float, optional
        Initial Monin-obukhov length, by default 1
    ustar_ref_in : float, optional
        Override ustar_ref, by default None
    min_windspeed : float, optional
        Max windspeed [m/s], by default 0.01
    MAX_ITERATIONS : int, optional
        Maximum iterations, by default 1
    TOLERANCE : float, optional
        Difference in ustar_ref to stop iterations at, by default 1e-3

    Returns
    -------
    Tuple[float, float]
        Ustar_ref, L

    """
    if ustar_ref_in is not None:
        L = calc_monin_obukhov_length(Tk, ustar_ref_in, Hd, P)
        return ustar_ref_in, L

    # Set initial state
    ustar_ref = min_ustar
    L = initial_L

    for _ in range(MAX_ITERATIONS):
        # ! Find ustar over reference canopy
        ustar_ref_next = max(min_ustar, ustar_from_velocity(max(min_windspeed, u), (z_u - u_d), u_z0, L))
        L = calc_monin_obukhov_length(Tk, ustar_ref_next, Hd, P)
        if isclose(ustar_ref_next, ustar_ref, abs_tol=TOLERANCE):
            ustar_ref = ustar_ref_next
            L = calc_monin_obukhov_length(Tk, ustar_ref_next, Hd, P)
            break
        ustar_ref = ustar_ref_next
    return ustar_ref, L



def calc_windspeed(
    h: float,
    u: float,
    z_u: float,
    u_z0: float,
    d: float,
    z0: float,
    L: float = None,
    u_d: float = None,
    ustar_ref: float = None,
    o_top_chamber: bool = False,
    min_windspeed: float = MIN_WINDSPEED,
    MIN_USTAR: float = 0.0001,
    izr: float = izR,
) -> NamedTuple:
    """Calculate wind speed parameters.

    aka met_windspeed

    Must supply ustar_ref or L

    Parameters
    ----------
    h: float
        Canopy height [m]
    u: float
        wind speed from met data [m/s]
    z_u: float
        Measured windspeed height[m]
    u_z0: float
        Measured canopy roughness length[m]
    d: float
        Canopy displacement height [m]
    z0: float
        Canopy roughness length [m]
    L: float
        Monin_Obukhov Length [m], optional
    u_d: float
        Measured Canopy displacement height [m], only needed for ustar_ref method
    ustar_ref: float
        Friction velocity [m/s], optional
    o_top_chamber: bool
        Is open top chamber [bool], default = False
    min_windspeed: float
        Minimum windspeed value [m/s], default = MIN_WINDSPEED
    MIN_USTAR: float
        Minimum friction velocity value [m/s], default = 0.00001
    izr: float
        Reference height [m], default = 50m

    Returns
    -------
    Named tuple containing:
    u_i: float
        windspeed at reference height[m/s]
    ustar: float
        Friction velocity[m/s]
    micro_u: float
        Windspeed at top of canopy[m/s]

    """
    # TODO: Check that ustar_ref calc is correct as inputed in UI
    Output = namedtuple('Output', 'u_i ustar micro_u')
    # TODO: Move these parameters to model setup

    # TODO: move these to parameter setup
    # # Measured canopy height
    # h_u = h if o_top_chamber else loc_h_u if loc_h_u is not None else h
    # # Measured windspeed height
    # z_u = h if o_top_chamber else loc_z_u if loc_z_u is not None else h

    # # Heights adjusted to match
    # # u_d = h_u * canopy_d
    # u_z0 = h_u * canopy_z0
    # d = h * canopy_d
    # z0 = h * canopy_z0

    z_u_a = z_u if not o_top_chamber else h
    u_z0_a = u_z0 if not o_top_chamber else z0
    u_d_a = u_d if not o_top_chamber else d

    # Find windspeed at izR, over reference canopy
    u_i = estimate_windspeed(u, z_u_a, izr, u_z0_a, d, 1/L) if L is not None \
        else velocity_from_ustar(ustar_ref, (izr - u_d_a), u_z0_a)
    u_i_lim =  max(min_windspeed,u_i)

    # Find ustar over target canopy, assuming that at izR windspeed will
    # be equal over both vegitations
    ustar =  ustar_from_velocity(u_i_lim, (izr - d), z0, L) if L is not None \
        else ustar_from_velocity_simple(u_i_lim, (izr - d), z0)
    ustar_lim  = max(MIN_USTAR,ustar)

    # Calculates windspeed at top layer canopy (uh)
    micro_u = estimate_windspeed(u_i_lim, izr, h, z0, d, 1/L) if L is not None \
        else velocity_from_ustar(ustar_lim, (h - d), z0)
    micro_u_lim =max(min_windspeed, micro_u)

    return Output(
        u_i=u_i_lim,
        ustar=ustar_lim,
        micro_u=micro_u_lim,
    )

@deprecated
def ustar_from_velocity_simple(
    u: float,
    z: float,
    z0: float,
    MIN_USTAR: float = 0.0001,
) -> float:
    """Estimate friction velocity from windspeed velocity.

    Parameters
    ----------
    u: float
        Velocity at height above boundary [m/s]
    z: float
        Height above boundary, e.g. z - d[m]
    z0: float
        Roughness length, height at which u=0[m]
    MIN_USTAR: float
        Minimum friction velocity value [m/s], default = 0.00001

    Returns
    -------
    ustar: float
        Output: friction velocity, ustar [m/s]

    """
    ustar = (u * K) / log(z / z0)
    ustar_lim = max(MIN_USTAR, ustar)
    return ustar_lim

def ustar_from_velocity(
    u: float,
    z: float,
    z0: float,
    L: float,
) -> float:
    """Estimate friction velocity from windspeed velocity.

    Parameters
    ----------
    u: float
        Velocity at height above boundary [m/s]
    z: float
        Height above boundary, e.g. z - d[m]
    z0: float
        Roughness length, height at which u=0[m]
    L: float
        Monin-Obukhov Length [m]

    Returns
    -------
    ustar: float
        Output: friction velocity, ustar [m/s]

    """

    psim_a = calc_PsiM(D(z) / D(L))
    psim_b = calc_PsiM(D(z0) / D(L))
    a = (D(z) / D(z0)).ln()
    ustar = D(u * K) / D(a - psim_a + psim_b)
    return float(ustar)

@deprecated
def velocity_from_ustar(
    ustar: float,
    z: float,
    z0: float
) -> float:
    """Estimate windspeed velocity from friction velocity.

    Parameters
    ----------
    ustar: float
        Friction velocity [m/s]
    z: float
        Height above boundary, e.g. z - d[m]
    z0: float
        Roughness length, height at which u=0[m]

    Returns
    -------
    u: float
        Output: velocity [m/s]

    """
    u = (ustar / K) * log(z / z0)
    return u


def estimate_windspeed(
    u_ref: float,
    z_ref: float,
    z: float,
    z0: float,
    d: float,
    invL: float,
) -> float:
    """Estimate the velocity by transfering measured velocity up or down.

    Parameters
    ----------
    u_ref : float
        Velocity at reference height[m/s]
    z_ref : float
        Reference height[m]
    z : float
        Target height [m]
    z0 : float
        Roughness length, height at which u=0 [m]
    d : float
        displacement height [m]
    invL : float
        Inverse Monin-Obukhov length [1/m]

    Returns
    -------
    float
        Windspeed [m/s]

    """
    # Updated to match EMEP
    psima = calc_PsiM(D((z-d)*invL))
    psimb = calc_PsiM(D(z0*invL))
    psimc = calc_PsiM(D((z_ref-d)*invL))
    a = ((D(z-d))/D(z0)).ln()
    b = ((D(z_ref-d))/D(z0)).ln()
    u = D(u_ref) * (a -psima + psimb)/(b - psimc + psimb)
    # u = u_ref * (log((z-d)/z0)  -psima + psimb)/(log((z_ref-d)/z0) -psimc + psimb)
    return float(u)


def calc_layer_windspeed(
    h: float,
    w: float,
    SAI: float,
    u_at_canopy_top: float,
    z: float,
    o_top_chamber: bool = False,
    layer_depth: float = 0,
) -> List[float]:
    """Calculate wind speeds at each layer in canopy.

    multi_layer_windspeed in DO3SE fortran model

    Parameters
    ----------
    h: float
        Canopy height[m]
    w: float
        Leaf width[m]
    SAI: float
        Stand area index [m2 m-2]
    u_at_canopy_top: float
        Wind speed at canopy height [m/s]
    z: float
        canopy layer height
    o_top_chamber: bool
        true if based on open top chamber model

    Returns
    -------
    u_z: List[float](iL,)
        Wind speed at each z [m/s]

    """
    if o_top_chamber:
        return u_at_canopy_top
    else:
        if z - layer_depth > h:
            # layer is still higher than top of canopy
            # TODO: Check this still works
            return u_at_canopy_top
        if w <= 0:
            return u_at_canopy_top
        assert w
        assert h
        # TODO: different definitions for lm
        # TODO: Check what we do when SAI is 0
        max_lm = 99999999
        lm = sqrt((4 * w * h) / (pi * SAI)) if SAI > 0 else max_lm

        a = sqrt((0.2 * SAI * h) / lm)

        # TODO: why do we `-1` here?
        u_z = u_at_canopy_top * exp(a * (z / h - 1))
        return u_z

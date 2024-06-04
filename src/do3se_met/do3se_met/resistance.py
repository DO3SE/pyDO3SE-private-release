"""Resistance module calculations (External to plant canopy)."""

from math import log, sqrt, pi, atan
from decimal import Decimal as D
from typing import Tuple
from deprecated import deprecated

from do3se_met.physical_constants import VON_KAR as K
from do3se_met.physical_constants import g, T0, cp, Rmass
from do3se_met.model_constants import izR


R_INF = 1000000


def calc_displacement_and_roughness_parameters(
    h: float,
    is_forest: bool,
) -> Tuple[float, float]:
    """Derive the displacement height (d) and roughness length (zo) of the

    vegetation under the windspeed measurement based on its height


    Parameters
    ----------
    h : float
        Height [m]
    is_forest : bool
        If true use forest vegitation parameters

    Returns
    -------
    Tuple[float,float]
        d, z0

    """
    if h is None:
        return None, None
    if is_forest:
        d = h * 0.78
        z0 = min(1.0, h * 0.07)
    else:
        d = h * 0.7
        z0 = max(0.001, h * 0.1)
    return d, z0


def calc_PsiH(zL: float) -> float:
    """Estimate the integral flux-gradient stability function for heat.

    Ref: Garratt, 1994, pp52-54
    VDHH modified - use van der Hurk + Holtslag?


    Parameters
    ----------
    zL : float
        Surface layer stability factor

    Returns
    -------
    float
        integral flux-gradient stability for heat

    """
    if zL < 0:  # unstable
        x = sqrt(1.0 - 16.0 * zL)
        stab_h = 2.0 * log((1.0 + x) / 2.0)
    else:             # stable
        # ESX if ( FluxPROFILE == "Ln95" ) then
        # ESX    stab_h = -( (1+2*a/3.0*zL)**1.5 + b*(zL-c/d)* exp(-d*zL) + (b*c/d-1) )
        # ESX else
        stab_h = -5.0 * zL
        # ESX end if
    return stab_h


def calc_PsiM(zL: D) -> D:
    """Estimate integral flux-gradient stability function for momentum.

    Out:
    PsiM = integral flux-gradient stability function for momentum
    Ref: Garratt, 1994, pp52-54 && EMEP discussions 2021

    NOTE: Use of Decimal

    Parameters
    ----------
    zL : float
        surface layer stability parameter, (z-d)/L
                                    ! notation must be preserved

    Returns
    -------
    float
        integral flux-gradient stability

    """
    if zL < 0:  # unstable
        x = sqrt(sqrt(D(1.0) - D(16.0) * zL))
        b = D(0.125 * (1.0 + x) * (1.0 + x) * (1.0 + x * x))
        stab_m = b.ln() + D(pi) / D(2.0) - D(2.0 * atan(x))
    else:             # stable
        # ESX if ( FluxPROFILE == "Ln95" ) then
        # ESX    stab_m = -( a*zL + b*(zl-c/d)*exp(-d*zL) + b*c/d)
        # ESX else
        stab_m = D(-5.0) * zL
        # ESX end if
    return stab_m


def calc_Ra_with_heat_flux(
    ustar: float,
    z1: float,
    z2: float,
    L: float,
) -> float:
    """Calculate Ra, Atmospheric resistance, taking into account heat flux data.

    Taken from DO3SE UI r.F90

    References
    ----------
    Garratt, 1994, pp.55-58
    Simpson et al., xxxx
    EMEP MSC-W Chemical transport Model

    Parameters
    ----------
    ustar : float
        Friction velocity [m/s]
    z1 : float
        lower height [m]
    z2 : float
        upper height [m]
    L : float
        Monin-Obukhov length [m]

    Returns
    -------
    float
        Air resistance

    """
    invL = 1 / L
    if z1 > z2:
        Ra = -999.0
    else:
        Ra = log(z2 / z1) - calc_PsiH(z2 * invL) + calc_PsiH(z1 * invL)
        Ra = Ra / (K * ustar)
    # TODO: Check if this should be allowed to go below 0
    Ra_lim = max(0, Ra)
    return Ra_lim


@deprecated
def calc_Ra_with_heat_flux_old(
    Ts_C: float,  # TODO: use kelvin as input instead of converting for each pass
    Hd: float,
    P: float,
    ustar: float,
    d: float,
    zo: float,
) -> float:
    """Calculate Ra, Atmospheric resistance, taking into account heat flux data.

    Taken from DO3SE UI r.F90

    References
    ----------
    Garratt, 1994, pp.55-58
    EMEP MSC-W Chemical transport Model

    Parameters
    ----------
    Ts_C : float
        Air Temperature [Degrees]
    Hd : float
        Sensible heat flux [W/m^2]
    P : float
        Air Pressure [kPa]
    ustar : float
        Friction velocity [m/s]
    d : float
        Displacement height (m)
    zo : float
        Aerodynamic roughness length (z0) [m]

    Returns
    -------
    Ra: float
        Atmospheric resistance
    """
    z2 = izR - d
    z1 = zo

    Tk = Ts_C + T0
    if Hd == 0:
        Hd = 0.000000000001

    # Surface density of dry air (including conversion from to Pa to kPa)
    rho = (P * 1000) / (Rmass * Tk)

    # Monin-Obukhov Length
    # TODO: Tk should be virtual temperature
    L = -ustar**3 * Tk * rho * cp / (K * g * Hd)

    Ezd = z2 / L
    Ezo = z1 / L

    # TODO: Check we don't need _m calculations. These came from DO3SE ui model
    if Ezd >= 0:
        Psi_h_zd = -5 * Ezd
        # Psi_m_zd = -5 * Ezd
    else:
        # Xzd_m = (1 - 16 * Ezd)**(1.0 / 4.0)
        Xzd_h = (1 - 16 * Ezd)**(1.0 / 2.0)
        # Psi_m_zd = log(((1 + Xzd_m**2) / 2) * ((1 + Xzd_m) / 2)**2) - 2 * atan(Xzd_m) + (pi / 2)
        Psi_h_zd = 2 * log((1 + Xzd_h**2) / 2)

    if Ezo >= 0:
        Psi_h_zo = -5 * Ezo
        # Psi_m_zo = -5 * Ezo
    else:
        # Xzo_m = (1 - 15 * Ezo)**(1 / 4.0)
        Xzo_h = (1 - 15 * Ezo)**(1 / 2.0)
        # Psi_m_zo = log(((1 + Xzo_m**2) / 2) * ((1 + Xzo_m) / 2)**2) - 2 * atan(Xzo_m) + (pi / 2)
        Psi_h_zo = 2 * log((1 + Xzo_h**2) / 2)

    Ra = (1 / (K * ustar)) * (log(z2 / z1) - Psi_h_zd + Psi_h_zo)
    return Ra


def calc_Ra_simple(ustar: float, z1: float, z2: float, d: float) -> float:
    r"""Calculate aerodynamic resistance (Ra, s m-1).

    between two heights using a simple, neutral stability model.
    Must satisfy z_2 <= z_1, z_2 > d and z_1 > d.

    Parameters
    ----------
    ustar: float
        Friction velocity [m/s]
    z1: float
        Lower height[m]
    z2: float
        Upper height[m]
    d: float
        Zero displacement height[m]
    """
    if z2 < z1:
        # TODO: Check that this is valid
        return 0
    if d > z2 or d > z1:
        # TODO: Check that this is valid
        return 0
    Ra = (1.0 / (ustar * K)) * log((z2 - d) / (z1 - d))
    return Ra


def calc_Rb(ustar: float, diff: float) -> float:
    """Calculate quasi-laminar boundary layer resistance (Rb, s m-1).

    based on a given friction velocity and diffusivity.

    References
    ----------
    Hicks et al (1987)

    Parameters
    ----------
    ustar: float
        Friction velocity [m/s]
    diff: float
        Molecular diffusivity in air (m2 s-1)

    Returns
    -------
    Rb: float
        quasi-laminar boundary layer resistance [s/m]
    """
    PR = 0.72    # Prandtl number
    V = 0.000015  # Kinematic viscosity of air at 20 C (m2 s-1)
    # K = K  # von Karman's constant

    Rb = (2.0 / (K * ustar)) * (((V / diff) / PR)**(2.0 / 3.0))
    return Rb


def calc_deposition_velocity(
    rmodel_Ra_c: float,
    rmodel_Rtotal_top_layer: float,
) -> float:
    """Calculate deposition velocity (\f$V_d\f$) from a canopy resistance model.

    Vd = 1/(Ra + Rb + Rsur)

    Parameters
    ----------
    rmodel_Ra_c : float
        Ra at canopy height
    rmodel_Rtotal_top_layer : float
        Should equal Rb + Rsur

    Returns
    -------
    float
        Deposition velocity

    """
    deposition_velocity = 1.0 / (rmodel_Ra_c + rmodel_Rtotal_top_layer)
    return deposition_velocity

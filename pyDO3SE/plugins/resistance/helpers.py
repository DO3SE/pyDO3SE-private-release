"""Resistance module calculations."""

from math import exp, inf, log, sqrt
from typing import List
from pyDO3SE.constants.physical_constants import VON_KAR as K
from pyDO3SE.constants.physical_constants import g, T0, cp, Rmass
from pyDO3SE.constants.model_constants import izR

from .constants import R_INF


def calc_Ra_with_heat_flux(
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
    raise Exception("Calc Ra moved to do3se_met package")
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
    raise Exception("Calc Ra moved to do3se_met package")
    if z2 > z1:
        # TODO: Check that this is valid
        return 0
    if d > z_2 or d > z_1:
        return 0
    Ra = (1.0 / (ustar * K)) * log((z2 - d) / (z1 - d))
    return Ra


def calc_Rb(ustar: float, diff: float) -> float:
    """Calculate quasi-laminar boundary layer resistance (Rb, s m-1).

    based on a given friction velocity and diffusivity.

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
    raise Exception("Calc Rb moved to do3se_met package")
    PR = 0.72    # Prandtl number
    V = 0.000015  # Kinematic viscosity of air at 20 C (m2 s-1)
    # K = K  # von Karman's constant

    Rb = (2.0 / (K * ustar)) * (((V / diff) / PR)**(2.0 / 3.0))
    return Rb


#   TODO: Check elemental
def calc_leaf_gb(
    G: float,
    Lm: float,
    u: float,
) -> float:
    """Calculate leaf-level quasi-laminar boundary layer conductance.

    (gb, mol m-2 s-1), for a particular kind of quantity specified by the base
    conductance of a single leaf surface, G.

    Parameters
    ----------
    G: float
        Leaf surface conductance [mol m-2 s-1]
    Lm: float
        Cross-wind leaf dimension[m]
    u: float
        Wind speed [m/s]
    """
    # G * 2 : from single surface to PLA (both sides of leaf)
    leaf_gb = (G * 2) * sqrt(u / Lm)
    return leaf_gb


def calc_leaf_rb(gb: float) -> float:
    """Calculate leaf-level quasi-laminar boundary layer resistance (rb, s m-1).

    from conductance (gb, mol m-2 s-1).

    Parameters
    ----------
    gb: float
        Leaf boundary layer conductance [mol m-2 s-1]

    Returns
    -------
    leaf_rb: float
        [Description]
    """
    # gb / 41 : 'mol m-2 s-1' to 'm s-1'
    # 1 / gb  : 'm s-1' to 's m-1'
    leaf_rb = 41 / gb
    return leaf_rb


def calc_Rinc(SAI: float, h: float, ustar: float) -> float:
    """Estimate in-canopy aerodynamic resistance (Rinc, s m-1).



    The in-canopy resistance (Rinc) determines the resistance to ozone transfer
    within the canopy and hence the amount of ozone available for deposition to
    the surface underlying the vegetation

    Rinc = b * SAI * (h/ustar)

    References
    ----------



    Parameters
    ----------
    SAI: float
        Stand area index [m2 m-2]
    h: float
        Vegetation height[m]
    ustar: float
        Friction velocity [m/s]

    Returns
    -------
    Rinc: float
        Canopy aerodynamic resistance [m/s]

    """
    MAX_RINC = inf  # How do we deal with infinity here
    Rinc_b: float = 14    # Rinc coefficient
    Rinc = Rinc_b * SAI * h / ustar if ustar > 0 else MAX_RINC
    return Rinc

#   !> Estimate in-canopy aerodynamic resistance (Rinc, s m-1).
#   !!
#   !! This is the experimental method developed for the Keenley grassland
#   !! multilayer model.
#   !!
#   !! TODO: decide if we should keep this method
#   pure real function Rinc_prototype(SAI, ustar)
#     real, intent(in) :: SAI   !< Stand area index [m2 m-2]
#     real, intent(in) :: ustar !< Friction velocity [m/s]

#     real, parameter :: Rinc_b = 14    ! Rinc coefficient

#     Rinc_prototype = Rinc_b * SAI * Rinc_b/ustar
#   end function Rinc_prototype


def calc_Rext(SAI: float) -> float:
    """Estimate external plant cuticle resistance (Rext, s m-1).

    Parameters
    ----------
    SAI: float
        Stand area index [m2 m-2]
    Returns
    -------
    Rext: float
        external plant cuticle resistance [s/m]
    """
    MAX_REXT = inf
    Rext_base: float = 2500
    Rext = Rext_base / SAI if SAI > 0 else MAX_REXT
    return Rext


def calc_Rsto(Gsto: float) -> float:
    """Convert stomatal conductance to stomatal resistance (Rsto, s m-1).

    The maximum stomatal resistance is capped to prevent infinite values
    when the conductance is 0.

    Parameters
    ----------
    Gsto: float
        Stomatal conductance [mmol m-2 s-1]

    Returns
    -------
    Rsto: float
        stomatal resistance [s/m]
    """
    MAX_RSTO: float = 100000
    # TODO: Conversion here should be temperature dependend. See UI changes

    # (gsto in m s-1) = 41000 * (gsto in mmol m-2 s-1)
    # (rsto in s m-1) = 1 / (gsto in m s-1)
    Rsto = min(MAX_RSTO, 41000.0 / Gsto) if Gsto != 0 else MAX_RSTO
    return Rsto


def calc_Rgs(
    Rsoil: float,
    snow_depth: float = None,
    Ts_C: float = None,
) -> float:
    """Calculate ground resistance.

    Snow input is optional.

    References
    ----------


    Parameters
    ----------
    Rsoil : float
        The rsoil base which is specific to landcover
    snow_depth : float, optional
        measured snow depth[m], by default None
    Ts_C : float, optional
        measured air temperature[degrees C], by default None

    Returns
    -------
    float
        Rgs - ground resistance
    """
    if snow_depth is None:
        return Rsoil
    else:
        Rgs_base = Rsoil  # TODO: Check this is correct
        F_t = exp(-0.2 * (1 + Ts_C))
        # TODO: Should calculate T2
        T2 = Ts_C  # Near surface air temperature (2m)
        # TODO: Check what S_d_max should be
        S_d_max = 10  # Snow depth when snow fraction assumed to be 1
        S_d = snow_depth
        f_snow = max(0, min(1, S_d / S_d_max))
        Rx_snow = 70 if Ts_C >= 1 else (70 * (2 - T2)) if (-1 >= Ts_C < 1) else 700
        Rgs_inv = ((1 - 2 * f_snow) / (F_t * Rgs_base)) + ((2 * f_snow) / Rx_snow)
        return 1 / Rgs_inv


def calc_Rsur_multilayer(
    nL: int,
    Rb: float,
    Rsto: List[float],
    Rext: List[float],
    LAI: List[float],
    SAI: List[float],
) -> List[float]:
    """Calculate per-layer surface resistance - combined Rb, Rsto and Rext.

    TODO: This outputs a different value to calc_Rsur when using single layer.
    Parameters

    ----------
    nL: int
        Number of layers
    Rb: float
        [Description]
    Rsto: List[float]
        [Description]
    Rext: List[float]
        [Description]
    LAI: List[float]
        [Description]
    SAI: List[float]
        [Description]

    TODO: per-layer Rb

    """
    Rsur = [None for _ in range(nL)]
    for i in range(nL):
        # TODO: let infinities happen and propagate through this? 1/Inf = 0?
        if (LAI[i] > 0):
            # LAI (and SAI) > 0, include Rsto and Rext components
            # TODO: Rsur should be whole canopy rather than multilayer
            Rsur[i] = Rb + 1 / (1 / Rsto[i] + 1 / Rext[i])
        elif (SAI[i] > 0):
            # Only SAI, omit the Rsto component
            Rsur[i] = Rb + Rext[i]
        else:
            # No foliage, very high resistance!
            # TODO: find a justification for this, probably based on Rsto
            # TODO: have an "R_INF" constant?
            Rsur[i] = R_INF
    return Rsur


def calc_Rsur(
    Rb: float,
    Rsto_c: float,
    Rext: float,
    Rinc: float,
    Rsoil: float,
    LAI: float,
    SAI: float,
) -> List[float]:
    """Calculate single layer surface resistance - combined Rb, Rsto and Rext.

    Parameters
    ----------
    Rb: float
        Quasi-laminar boundary layer resistance [s m-1]
    Rsto: List[float]
        Bulk Stomatal resistance [s m-1]) per layer
    Rext: List[float]
        External plant cuticle resistance [s m-1] per layer
    LAI: List[float]
        Leaf area index [m^2 m^2]
    SAI: List[float]
        Stand area index [m^2 m^2]

    TODO: per-layer Rb
    """
    if LAI > 0:
        Rsur = 1 / ((1 / Rsto_c) + (SAI / Rext) + (1 / (Rinc + Rsoil)))
    elif SAI > 0:
        Rsur = 1 / ((SAI / Rext) + (1 / (Rinc + Rsoil)))
    elif SAI == 0:
        # surely this is Rsur = Rsoil ?
        Rsur = 1 / (1 / Rsoil)
    return Rsur


def calc_Rtotal_reversed(
    nL: int,
    Rsur: List[float],
    Rinc: List[float],
    Rgs: float,
) -> List[float]:
    """Calculate multi-layer Rtotal.

    The total resistance for each layer and everything below that layer.

    The equation below will tend towards the lowest of Rsur or Rinc+prev_layer_total.
    If the previous layer is 0 then the resistance is equal to Rsur + Rinc

    NOTE: Top layer is layer 0
    This matches original fortran version.

    Parameters
    ----------
    nL: int
        Number of model layers
    Rsur: List[float]
        Combined surface resistance [s m-1] per layer
    Rinc: List[float]
        In-canopy aerodynamic resistance [s m-1] per layer
    Rgs: float
        Ground surface resistance [s m-1]

    """
    tmp = [None for _ in range(nL + 1)]
    tmp[-1] = Rgs
    for i in range(nL-1, -1, -1):
        tmp[i] = 1 / (1 / Rsur[i] + 1 / (Rinc[i] + tmp[i+1]))
    Rtotal = tmp[0:nL]
    return Rtotal

def calc_Rtotal(
    nL: int,
    Rsur: List[float],
    Rinc: List[float],
    Rgs: float,
) -> List[float]:
    """Calculate multi-layer Rtotal.

    The total resistance for each layer and everything below that layer.

    The equation below will tend towards the lowest of Rsur or Rinc+prev_layer_total.
    If the previous layer is 0 then the resistance is equal to Rsur + Rinc

    NOTE: Bottom layer is layer 0.

    Parameters
    ----------
    nL: int
        Number of model layers
    Rsur: List[float]
        Combined surface resistance [s m-1] per layer
    Rinc: List[float]
        In-canopy aerodynamic resistance [s m-1] per layer
    Rgs: float
        Ground surface resistance [s m-1]

    """
    tmp = [None for _ in range(nL + 1)]
    tmp[0] = Rgs
    for i in range(1, nL+1):
        tmp[i] = 1 / (1 / Rsur[i-1] + 1 / (Rinc[i-1] + tmp[i-1]))
    Rtotal = tmp[1:nL+1]
    return Rtotal

def calc_deposition_velocity(
    rmodel_Ra_c: float,
    rmodel_Rtotal_top_layer: float,
) -> float:
    """Calculate deposition velocity (\f$V_d\f$) from a canopy resistance model.

    Parameters
    ----------
    rmodel_Ra_c : float
        Aerodynamic resistance [s m-1] between 50m and inside the canopy.
    rmodel_Rtotal_top_layer : float
        Total resistance for each layer downwards [s m-1]

    Returns
    -------
    float
        Deposition velocity
    """
    deposition_velocity = 1.0 / (rmodel_Ra_c + rmodel_Rtotal_top_layer)
    return deposition_velocity

# TODO: Implement below if needed
#   !> Calculate the rate of stomatal flux for a leaf resistance model.
#   elemental real function stomatal_flux_rate(leaf_rmodel)
#     type(LeafResistanceModel_t), intent(in) :: leaf_rmodel

#     real :: leaf_r

#     leaf_r = 1.0 / ((1.0 / leaf_rmodel%Rsto) + (1.0 / leaf_rmodel%Rext))
#     stomatal_flux_rate = (1.0/leaf_rmodel%Rsto) * (leaf_r / (leaf_rmodel%Rb + leaf_r))
#   end function stomatal_flux_rate

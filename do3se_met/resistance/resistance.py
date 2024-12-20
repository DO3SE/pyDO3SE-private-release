"""Resistance module calculations (External to plant canopy)."""

from math import log, sqrt, pi, atan, inf, exp
from decimal import Decimal as D
from typing import Tuple, List
from deprecated import deprecated

from do3se_met.physical_constants import g, T0, cp, Rmass, VON_KAR as K
from do3se_met.model_constants import izR
from do3se_met.resistance.model import Resistance_Model, Leaf_Resistance_Model

R_INF = 1000000
NOT_SET = -99999999.0


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
    else:  # stable
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
    else:  # stable
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
    L = -(ustar**3) * Tk * rho * cp / (K * g * Hd)

    Ezd = z2 / L
    Ezo = z1 / L

    # TODO: Check we don't need _m calculations. These came from DO3SE ui model
    if Ezd >= 0:
        Psi_h_zd = -5 * Ezd
        # Psi_m_zd = -5 * Ezd
    else:
        # Xzd_m = (1 - 16 * Ezd)**(1.0 / 4.0)
        Xzd_h = (1 - 16 * Ezd) ** (1.0 / 2.0)
        # Psi_m_zd = log(((1 + Xzd_m**2) / 2) * ((1 + Xzd_m) / 2)**2) - 2 * atan(Xzd_m) + (pi / 2)
        Psi_h_zd = 2 * log((1 + Xzd_h**2) / 2)

    if Ezo >= 0:
        Psi_h_zo = -5 * Ezo
        # Psi_m_zo = -5 * Ezo
    else:
        # Xzo_m = (1 - 15 * Ezo)**(1 / 4.0)
        Xzo_h = (1 - 15 * Ezo) ** (1 / 2.0)
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
        return 0
    z1_lim = max(0.000001, z1)
    Ra = (1.0 / (ustar * K)) * log((z2 - d) / (z1_lim))
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
    PR = 0.72  # Prandtl number
    V = 0.000015  # Kinematic viscosity of air at 20 C (m2 s-1)
    # K = K  # von Karman's constant

    Rb = (2.0 / (K * ustar)) * (((V / diff) / PR) ** (2.0 / 3.0))
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
    if gb <= 0:
        return inf
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
    Rinc_b: float = 14  # Rinc coefficient
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


def calc_Rext(Rext_base: float, SAI: float) -> float:
    """Estimate external plant cuticle resistance (Rext, s m-1).

    Parameters
    ----------
    Rext_base: float
        base rext value
    SAI: float
        Stand area index [m2 m-2]
    Returns
    -------
    Rext: float
        external plant cuticle resistance [s/m]

    """
    MAX_REXT = inf
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
    Rsto = min(MAX_RSTO, 41000.0 / Gsto) if Gsto > 0 else MAX_RSTO
    return Rsto


def calc_Rgs(
    Rsoil: float,
    snow_depth: float | None = None,
    Ts_C: float | None = None,
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
        if Ts_C is None:
            raise ValueError("Ts_C must be provided if snow depth is provided")
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


    """
    Rsur = [NOT_SET for _ in range(nL)]
    for i in range(nL):
        # TODO: let infinities happen and propagate through this? 1/Inf = 0?
        if LAI[i] > 0:
            # LAI (and SAI) > 0, include Rsto and Rext components
            # TODO: Rsur should be whole canopy rather than multilayer
            Rsur[i] = Rb + 1 / (1 / Rsto[i] + 1 / Rext[i])
        elif SAI[i] > 0:
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
) -> float:
    """Calculate single layer surface resistance - combined Rb, Rsto and Rext.

    NOTE: Rext term = Rext_base/SAI, Rext_base = 2500

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


    """
    if LAI > 0:
        Rsur = 1 / ((1 / Rsto_c) + 1 / Rext + (1 / (Rinc + Rsoil)))
    elif SAI > 0:
        Rsur = 1 / ((1 / Rext) + (1 / (Rinc + Rsoil)))
    elif SAI == 0:
        # surely this is Rsur = Rsoil ?
        Rsur = 1 / (1 / Rsoil)
    else:
        raise ValueError(f"Could not calculate Rsur as LAI={LAI} and SAI={SAI}")
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
    tmp = [NOT_SET for _ in range(nL + 1)]
    tmp[-1] = Rgs
    for i in range(nL - 1, -1, -1):
        tmp[i] = 1 / (1 / Rsur[i] + 1 / (Rinc[i] + tmp[i + 1]))
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
    tmp = [NOT_SET for _ in range(nL + 1)]
    tmp[0] = Rgs
    for i in range(1, nL + 1):
        tmp[i] = 1 / (1 / Rsur[i - 1] + 1 / (Rinc[i - 1] + tmp[i - 1]))
    Rtotal = tmp[1 : nL + 1]
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


def calc_resistance_model(
    nL: int,
    nLC: int,
    ustar: float,
    canopy_height: float,
    SAI_values: List[List[float]],
    LAI_values: List[List[float]],
    mean_gsto_values: List[List[float]],
    Rsoil: float,
    L: float|None = None,
    ra_calc_method: str = "simple",
    rsur_calc_method: str = "single_layer",
    Ts_C: float|None = None,
    snow_depth: float|None = None,
    Rext_base: float = 2500.0,
    Rb_diff: float = 0.000015,
    izr: float = izR,
    measured_height: float = 10,
    MIN_CANOPY_HEIGHT: float = 0.01,
    CANOPY_D: float = 0.7,
    CANOPY_Z0: float = 0.1,
) -> Resistance_Model:
    """Calculate the resistance model for O3 over the target canopy.

    Parameters
    ----------
    nL: int
        Number of model layers
    nLC: int
        Number of model components
    ustar: float
        Friction velocity over target canopy [m/s]
    ustar_ref: float
        Friction velocity over reference canopy [m/s]
    canopy_height: float
        The total height of the canopy [m]
    SAI_values: List[float]
        full model SAI values with shape(nL, nLC)
    LAI_values: List[float]
        full model LAI values with shape(nL, nLC)
    mean_gsto_values: float
        full model mean_gsto values with shape(nL, nLC) [mmol O3 m-2 PLA s-1]
    Rsoil: float
        Soil resistance?
    L: float
        Monin-Obukhov length [m]
    ra_calc_method: str
        Method of calculating Ra. Options:
        - "simple" -
        - "heat_flux" -
    rsur_calc_method: str
        Method of calculating Rsur. Options:
        - "single_layer" -
        - "multi_layer" -
    Ts_C: float
        Air Temperature [Degrees]
    P: float
        Air Pressure [kPa] only needed for Ra heat flux calc
    snow_depth: float
        Measured snow depth. Default None
    Rext_const: float
        Rext base value
    Rb_diff: float
        Molecular diffusivity in air (m2 s-1)
    measured_height: float
        The height that the measurement was taken at [m]
    izr: float
        Decoupled height [m]
    MIN_CANOPY_HEIGHT: float = 0.01
        Use min canopy height to avoid div_zero error
    CANOPY_D: float = 0.7
        Canopy displacement (fraction of canopy height)[frac]
    CANOPY_Z0: float = 0.1
        Canopy roughness length (fraction of canopy height)[frac]

    Returns
    -------
    resistance_model: Resistance_Model
        A new instance of the Resistance Model

    """
    if nL > 1 and rsur_calc_method == "single_layer":
        raise ValueError("Cannot use single layer rsur calc for multilayer model")

    if ra_calc_method == "heat_flux" and L is None:
        raise ValueError("Monin-Obukhov length must be provided for heat flux Ra calc")

    TOP_LAYER = -1
    LAI_sum_per_layer = [sum(LAI_values[iL]) for iL in range(nL)]
    SAI_sum_per_layer = [sum(SAI_values[iL]) for iL in range(nL)]
    # bulk_gsto_sum_per_layer = [sum([mean_gsto_values[iL][iLC] * LAI_values[iL][iLC]
    #                                for iLC in range(nLC)]) for iL in range(nL)]
    # TODO: Check that mean gsto correct here
    mean_gsto_per_layer = [
        sum(
            [
                mean_gsto_values[iL][iLC] * LAI_values[iL][iLC] / LAI_sum_per_layer[iL]
                if LAI_sum_per_layer[iL] > 0
                else 0
                for iLC in range(nLC)
            ]
        )
        for iL in range(nL)
    ]

    # Ra between reference canopy and izR
    canopy_height_lim = max(MIN_CANOPY_HEIGHT, canopy_height)
    canopy_height_d = canopy_height_lim * CANOPY_D
    canopy_height_zo = canopy_height_lim * CANOPY_Z0

    Ra_canopy_to_izr = (
        calc_Ra_with_heat_flux(
            ustar,
            z1=canopy_height_zo,  # TODO: Check z0 - d
            z2=izr,
            L=L,
        )
        if ra_calc_method == "heat_flux" and L is not None
        else calc_Ra_simple(
            ustar,
            z1=canopy_height_zo,
            z2=izr,
            d=canopy_height_d,
        )
    )

    Ra_canopy_top_to_izr = (
        calc_Ra_with_heat_flux(
            ustar,
            z1=canopy_height_zo,  # TODO: Check z0 - d
            z2=izr,
            L=L,
        )
        if ra_calc_method == "heat_flux" and L is not None
        else calc_Ra_simple(
            ustar,
            z1=canopy_height_lim - canopy_height_d,
            z2=izr,
            d=canopy_height_d,
        )
    )

    Rb = calc_Rb(ustar, Rb_diff)

    Ra_measured_to_izr = (
        calc_Ra_with_heat_flux(
            ustar,
            z1=measured_height,
            z2=izr,
            L=L,
        )
        if ra_calc_method == "heat_flux" and L is not None
        else calc_Ra_simple(
            ustar,
            z1=measured_height - canopy_height_d,
            z2=izr,
            d=canopy_height_d,
        )
    )

    # TODO: Rinc should use layer depth instead of canopy height
    # TODO: to use in a multilayer model, what does h represent?  Height above
    # ground, or thickness of layer?
    # Rinc: List[float] = [calc_Rinc(SAI_sum_per_layer[iL], layer depth, ustar) for iL in range(nL)]
    Rinc: List[float] = [
        calc_Rinc(sum(SAI_sum_per_layer), canopy_height_lim, ustar) for iL in range(nL)
    ]

    Rext: List[float] = [calc_Rext(Rext_base, SAI_sum_per_layer[iL]) for iL in range(nL)]

    Rsto: List[float] = [calc_Rsto(mean_gsto_per_layer[iL]) for iL in range(nL)]
    # Rsto_c: List[float] = [calc_Rsto(bulk_gsto_sum_per_layer[iL]) for iL in range(nL)]
    Rgs = calc_Rgs(Rsoil, snow_depth, Ts_C)

    # TODO: Should allow input LAI here
    if rsur_calc_method == "single_layer":
        Rsur: list = [
            calc_Rsur(
                Rb,
                Rsto[TOP_LAYER],
                Rext[TOP_LAYER],
                Rinc[TOP_LAYER],
                Rgs,
                LAI_sum_per_layer[TOP_LAYER],
                SAI_sum_per_layer[TOP_LAYER],
            )
            for _ in range(nL)
        ]
    elif rsur_calc_method == "multi_layer":
        Rsur: List = calc_Rsur_multilayer(nL, Rb, Rsto, Rext, LAI_sum_per_layer, SAI_sum_per_layer)
    else:
        raise ValueError("Invalid Rsur calc method")

    Rtotal: List[float] = calc_Rtotal(nL, Rsur, Rinc, Rgs)

    return Resistance_Model(
        nL,
        Ra_measured_to_izr=Ra_measured_to_izr,
        Ra_canopy_to_izr=Ra_canopy_to_izr,
        Ra_canopy_top_to_izr=Ra_canopy_top_to_izr,
        Rb=Rb,
        Rinc=Rinc,
        Rext=Rext,
        Rsto=Rsto,
        Rgs=Rgs,
        Rsur=Rsur,
        Rtotal=Rtotal,
    )


def calc_leaf_resistance_model(
    nL: int,
    Lm: float,
    u_per_layer: List[float],
    leaf_gsto_per_layer: List[float],
    LEAF_G: float = 0.105,
) -> Leaf_Resistance_Model:
    """Setup the leaf O3 resistance model.

    Parameters
    ----------
    nL: int
        Number of model layers
    Lm: float
        Leaf dimension
    u_per_layer: float
        wind speed for each layer [m/s]
    leaf_gsto_per_layer: List[float]
        leaf_gsto [mmol m-2 s-1]
    LEAF_G: float
        ?? 0.105 for Ozone

    Returns
    -------
    Leaf_Resistance_Model
    """
    Rb_per_layer: List[float] = [
        calc_leaf_rb(calc_leaf_gb(LEAF_G, Lm, u_per_layer[iL])) for iL in range(nL)
    ]

    # TODO: Implement advanced Rext equation
    # Calculating Rext currently disabled to align with DO3SE ui
    # Rext_per_layer: List[List[float]] = [[calc_Rext(1.0)
    #                                       for iLC in range(nLC)]
    #                                      for iL in range(nL)]

    # TODO: Get Rext_const from config
    Rext_per_layer: List[float] = [2500 for iL in range(nL)]

    Rsto_per_layer: List[float] = [calc_Rsto(leaf_gsto_per_layer[iL]) for iL in range(nL)]

    leaf_r_model = Leaf_Resistance_Model(
        Rb=Rb_per_layer,
        Rext=Rext_per_layer,
        Rsto=Rsto_per_layer,
    )
    return leaf_r_model

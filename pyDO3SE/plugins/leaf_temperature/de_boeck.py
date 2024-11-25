"""Functions to calculate leaf temperature using the de Boeck method.


References
----------
- de Boeck, 2012 - De Boeck, H.J., De Groote, T. and Nijs, I., 2012. Leaf temperatures in glasshouses and open‐top chambers. New Phytologist, 194(4), pp.1155-1164.

"""
# %%
from math import sqrt, exp
from do3se_met.helpers import saturated_vapour_pressure
from pyDO3SE.constants.physical_constants import T0


def combined_leaf_conductance(
    gs_ab: float,
    gs_ad: float,
    ga: float,
) -> float:
    """Calculate the combined leaf conductance.

    Parameters
    ----------
    gs_ab : float
        Abaxial stomatal conductance (mol m-2 s-1)
    gs_ad : float
        Adaxial stomatal conductance (mol m-2 s-1)
    ga : float
        Boundary layer conductance (mol m-2 s-1)

    Returns
    -------
    float
        combined leaf conductance (mol m-2 s-1)

    """
    return (0.5 * gs_ab * ga / (gs_ab + ga)) + (0.5 * gs_ad * ga / (gs_ad + ga))


def get_leaf_temp_de_boeck(
    R: float,
    eact: float,
    T_air: float,
    initial_T_leaf: float,
    P: float,
    u_speed: float,
    g_vs: float,
    hypostomatous: bool,
    d: float,
    albedo: float,
    cloud_cover: float = 1.0,
    balance_threshold: float = 0.0010000000474974513,
    adjustment_factor: float = 0.019999999552965164,
    max_iterations: int = 50,
) -> float:
    """Calculate the leaf temperature using the de_Boeck method

    Parameters
    ----------
    R : float
        Global radiation (W m-2)
    eact : float
        Vapour pressure (kPa)
    T_air : float
        Air temperature (degrees C)
    initial_T_leaf : float
        T_leaf starting value (degrees C)
    P : float
        Air pressure (kPa)
    u_speed : float
        Wind speed (m s-1)
    g_vs : float
        Stomatal conductance to water vapour (umol m-2 s-1)
    hypostomatous : bool
        Are leaves hypostomatous?
    d : float
        Leaf characteristic dimension (m)
    albedo : float
        Surface albedo (fraction)
    cloud_cover : float
        Cloud cover (fraction)
    balance_threshold : float
        Threshold within which to accept leaf energy balance
    adjustment_factor : float
        Multiplier for energy balance based T_leaf adjustments
    max_iterations : int = 50
        Maximum iteration count

    Returns
    -------
    float
        Leaf temperature [DegreesC]

    """
    e_a = eact * 1e3
    pres = P * 1e3
    g_vs_mol = g_vs * 1e-6
    # Copied from DO3SE-UI pure function leaf_temp_de_Boeck
    # TODO: Add reference
    # Stefan-Boltzmann constant, W m-2 K-4 (Campbell & Norman (1998), p281, table A5)
    SBC: float = 5.670373e-8
    # Specific heat capacity of dry air, J mol-1 C-1 (Campbell & Norman (1998), p279, table A1)
    c_p: float = 29.3

    # Leaf short-wave absorptivity (de Boeck (2012); from Campbell & Norman (1998), p153, in text)
    alpha_s_leaf: float = 0.5
    # Leaf long-wave absorptivity (de Boeck, 2012)
    alpha_l_leaf: float = 0.97
    # Soil long-wave emissivity (de Boeck, 2012)
    eps_l_soil: float = 0.945
    # Leaf long-wave emissivity (de Boeck, 2012)
    eps_l_leaf: float = 0.97

    # Soil short-wave reflectivity
    rho_s_soil = albedo
    # Water vapour path length (de Boeck 2012)
    VPL = 46.5 * ((0.01 * e_a) / (T_air + T0))
    # Clear sky emissivity (de Boeck 2012)
    eps_ac = 1 - (1 + VPL) * exp(-(1.2 + 3.0 * VPL)**0.5)
    # Sky long-wave emmisivity (de Boeck (2012); from Campbell & Norman (1998), p164, eq. 10.12)
    eps_l_sky = (1 - 0.84 * cloud_cover) * eps_ac + 0.84 * cloud_cover
    # Latent heat of vapourisation, J mol-1 (de Boeck (2012))
    lam = -42.575 * T_air + 44994

    # Assume soil temperature is equal to air temperature
    # TODO: use the banded estimate from de Boeck (2012)?
    T_soil = T_air

    # Starting point
    T_leaf = initial_T_leaf
    T_leaf_adjust = 0.0

    for _ in range(max_iterations):
        # Apply adjustment (from previous iteration)
        T_leaf = T_leaf + T_leaf_adjust

        # Boundary layer conductances (de Boeck (2012))
        # XXX: Assume forced convection:
        g_Ha = 1.4 * 0.135 * sqrt(u_speed / d)
        g_vb = 1.4 * 0.147 * sqrt(u_speed / d)

        # Total conductivity to water vapour
        if (hypostomatous):
            g_v = combined_leaf_conductance(g_vs_mol, 0.0, g_vb)
        else:
            g_v = combined_leaf_conductance(g_vs_mol, g_vs_mol, g_vb)

        # Saturated vapour pressure for leaf temperature, Pa
        e_s_Tleaf = saturated_vapour_pressure(T_leaf) * 1000  # Converted from kPa to Pa

        # Absorbed short-wave (direct on top of leaf, reflected from soil on underside) (de Boeck (2012))
        R_s_in = (0.5 * R * alpha_s_leaf) + (0.5 * R * rho_s_soil * alpha_s_leaf)

        # Absorbed long-wave (emitted from soil on underside, emitted from sky on top) (de Boeck (2012))
        R_l_in = (0.5 * alpha_l_leaf * eps_l_soil * SBC * (T_soil + T0)**4) + \
            (0.5 * alpha_l_leaf * eps_l_sky * SBC * (T_air + T0)**4)

        # Emitted long-wave (de Boeck (2012))
        R_l_out = eps_l_leaf * SBC * (T_leaf + T0)**4

        # Sensible heat flux (de Boeck (2012))
        H = g_Ha * c_p * (T_leaf - T_air)

        # Latent heat flux (de Boeck (2012))
        lam_E = lam * g_v * (e_s_Tleaf - e_a) / pres

        # Energy balance equation (de Boeck (2012))
        energy_balance = R_s_in + R_l_in - R_l_out - H - lam_E

        if (abs(energy_balance) <= balance_threshold):
            break
        else:
            T_leaf_adjust = energy_balance * adjustment_factor
    T_leaf_adjusted = T_air + 0.2 * (T_leaf - T_air)
    return T_leaf_adjusted

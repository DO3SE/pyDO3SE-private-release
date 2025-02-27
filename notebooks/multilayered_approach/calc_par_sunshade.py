# %%
"""Calculations for PAR sun and shade."""
from functools import partial
from typing import Tuple
from math import cos, pi, exp
from scipy.integrate import quad as integrate

# %%
from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

# %%
# Model coefficients
sigma = 0.15  # leaf scattering coefficient of PAR (p_i + T_i)
k_d_alt = 0.719  # diffuse and scattered diffuse PAR extinction coefficient, 0.719
seaP = 101.325  # real, parameter


# Existing equations=======

# %%
# calc_Idrctt_Idfuse
@interact(
    PAR=widgets.FloatSlider(min=0.001, max=800, step=100, value=300),
    sinB=widgets.FloatSlider(min=0, max=2 * pi, step=0.1, value=0.5),
    P=fixed(90),
)
def calc_Idrctt_Idfuse(
    PAR: float,
    sinB: float,
    P: float,
) -> float:
    """Estimate diffuse and direct PAR components.

    Args:
    PAR: float       !< Photosynthetically active radiation [W m-2]
    sinB: float      !< sin() of solar elevation angle
    P: float         !< Atmospheric pressure [kPa]

    Returns:
    Idrctt: float   !< Direct PAR irradiance [W m-2]
    Idfuse: float   !< Diffuse PAR irradiance [W m-2]

    real :: m, pPARdir, pPARdif, pPARtotal, ST, fPARdir, fPARdif
    """

    if sinB > 0.0:
        m = 1.0 / sinB

        # Potential direct PAR
        pPARdir = 600 * exp(-0.185 * (P / seaP) * m) * sinB
        # Potential diffuse PAR
        pPARdif = 0.4 * (600 - pPARdir) * sinB
        # Potential total PAR
        pPARtotal = pPARdir + pPARdif

        # Sky transmissivity
        ST = max(0.21, min(0.9, PAR / pPARtotal))

        # Direct and diffuse fractions
        fPARdir = (pPARdir / pPARtotal) * (1.0 - ((0.9 - ST) / 0.7)**(2.0/3.0))
        fPARdif = 1 - fPARdir

        # Apply calculated direct and diffuse fractions to PARtotal
        Idrctt = fPARdir * PAR
        Idfuse = fPARdif * PAR
    else:
        Idrctt = 0.0
        Idfuse = 0.0

    return Idrctt, Idfuse


# ================

# %%
# calc_beam_irradiance_horiz
@interact(
    sigma=sigma,
)
def calc_beam_irradiance_horiz(
    sigma: float = sigma,
) -> float:
    """Calculate the beam irradiance for horiontal leaves."""
    P_h = (1 - (1 - sigma) ** 0.5) / (1 + (1 - sigma) ** 0.5)
    return P_h


# %%
# calc_beam_irradiance_uad
@interact(
    P_h=widgets.FloatSlider(min=0.1, max=1.0, step=0.1, value=0.2),
    k_b=widgets.FloatSlider(min=0.1, max=10.0, step=0.1, value=0.2),
)
def calc_beam_irradiance_uad(
    P_h: float,  # reflection coefficient of a canopy with horizontal leaves
    k_b: float,  # Beam radiation extinction coefficient (0.5 / sinB)
) -> float:
    """Calculate beam irradiance for uniform leaf angle distribution."""
    # TODO: Check this as square brackets inside of exp
    P_cb = 1 - exp(2 * P_h * k_b / (1 + k_b))
    return P_cb


# %%
# calc_diffuse_irradiance_refl
@interact(
    sinB=widgets.FloatSlider(min=0.01, max=2*pi, step=0.1, value=0.5),
    P_h=widgets.FloatSlider(min=0.01, max=2*pi, step=0.1, value=0.5),
    Ir_dfuse_0=widgets.FloatSlider(min=0.01, max=500, step=50, value=200),
)
def calc_diffuse_irradiance_refl(
    # P_cb: float,  # beam irradiance for uniform leaf angle distribution
    sinB: float,
    P_h: float,
    # N_d: float,  #
    Ir_dfuse_0: float,  # Diffuse PAR per unit ground area at top of canopy[umol m^-2 s^-2]
) -> float:
    """Calculate diffuse irradiance reflection coefficient EQ A21."""
    # f
    def f(alpha):
        """Calculate integral per rotation unit of sun[radian]"""
        # TODO: Check this is correct
        # We should be
        k_b = 0.5 / sinB
        # TODO: This probably varies based on leaf normal
        # k_b = 0.5 / (sinB * cos(alpha))
        # Diffuse photon radiance of the sky (per radian??)[umol m^-2 sr ^-1]
        N_d = Ir_dfuse_0 / (2 * pi)
        P_cb = calc_beam_irradiance_uad(P_h, k_b)
        return N_d * P_cb

    integ, err = integrate(f, 0, pi/2)
    P_cd = (1 / Ir_dfuse_0) * integ
    return P_cd


# %%
# calc_Ir_scattered_b
def calc_Ir_scattered_b(
    P_cb: float,  # beam irradiance for uniform leaf angle distribution
    Ir_beam_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    LAI_c: float,  # cumulative leaf-area index from top of canopy (L = 0 at top)
    k_b: float,  # Beam radiation extinction coefficient
    k_b_alt: float,  # beam and scattered beam PAR extinction coefficient, 0-46/sinj3
    sigma: float = sigma,  # leaf scattering coefficient of PAR (pi + Ti)
) -> float:
    """Eq ??."""
    # TODO: Check this (Large square brackets)
    Ir_bs = Ir_beam_0 * (
        ((1 - P_cb)*k_b_alt * exp(-k_b * LAI_c)) /
        -(1 - sigma) * k_b * exp(-k_b * LAI_c)
    )
    return Ir_bs


# %%
# calc_Ir_beam_sun
def calc_Ir_beam_sun(
    sinB: float,  # sine of solar elevation angle [radians]
    cosA: float,  # cosine of angle of beam irradiance to the leaf normal [radians]
    Ir_beam_0: float,  # beam PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    sigma: float = sigma,  # leaf scattering coefficient of PAR
) -> float:
    """Eq A11."""
    Ir_b = (1 - sigma) * Ir_beam_0 * cosA / sinB
    return Ir_b


# %%
# calc_Ir_diffuse
def calc_Ir_diffuse(
    P_cd: float,  # diffuse irradiance reflection coefficient
    Ir_dfuse_0: float,  # diffuse PAR per unit ground area at top of canopy
    LAI_c: float,  # Cumulative leaf area index from top of canopy
    k_d_alt: float = k_d_alt,  # Diffuse and scattered diffuse PAR extinction coefficient
) -> float:
    """Diffuse PAR per unit ground area Eq A5."""
    ir_diffuse = (1 - P_cd) * k_d_alt * Ir_dfuse_0 * exp(-k_d_alt * LAI_c)
    return ir_diffuse


# %%
# calc_PAR_shade
def calc_PAR_shade(
    Ir_diffuse: float,  # diffuse PAR per unit ground area
    Ir_scattered_b: float,  # absorbed scattered beam PAR per unit leaf area
) -> float:
    """Eq A7."""
    return Ir_diffuse + Ir_scattered_b


# %%
# calc_PAR_sun
def calc_PAR_sun(
    PAR_shade: float,  # Irradiance absorbed by shaded leaves
    Ir_beam_sun: float,  # Beam irradiance absorbed by sunlit leaves
) -> float:
    """Calculate total irradiance on sunlit leaves.

    Farquhar 1997 - Eq A12

    Parameters
    ----------
    PAR_shade : float
        Irradiance received by shaded leaves

    Returns
    -------
    PARsun: float
        Total irradiance on sunlit leaves.
    """
    return PAR_shade + Ir_beam_sun


# %%
# calc_PAR_sun_shade
def calc_PAR_sun_shade(
    Ir_beam_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    Ir_dfuse_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    sinB: float,
    cosA: float,
    LAI_c: float,  # cumulative leaf-area index from top of canopy (L = 0 at top)
) -> Tuple[float, float]:
    """Calculate the sun and shade PAR values.

    Uses equations from Farquhar 1997 (Simple scaling of photosynthesis from leaves to canopies
    without the errors of big-leaf models)
    """
    k_b = 0.5 / sinB  # Beam radiation extinction coefficient
    k_b_alt = 0.46 / sinB  # beam and scattered beam PAR extinction coefficient
    P_h = calc_beam_irradiance_horiz()
    P_cb = calc_beam_irradiance_uad(P_h, k_b)
    P_cd = calc_diffuse_irradiance_refl(sinB, P_h, Ir_dfuse_0)
    Ir_diffuse = calc_Ir_diffuse(k_d_alt, P_cd, Ir_dfuse_0, LAI_c)
    Ir_beam_sun = calc_Ir_beam_sun(sinB, cosA, Ir_beam_0)
    Ir_scattered_b = calc_Ir_scattered_b(P_cb, Ir_beam_0, LAI_c, k_b, k_b_alt)
    PAR_shade = calc_PAR_shade(Ir_diffuse, Ir_scattered_b)
    PAR_sun = calc_PAR_sun(PAR_shade, Ir_beam_sun)
    return PAR_sun, PAR_shade


# %%
# calc_PAR_sun_shade_b
def calc_PAR_sun_shade_b(
    Ir_beam_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    Ir_dfuse_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    sinB: float,
    cosA: float,
    LAI_c: float,  # cumulative leaf-area index from top of canopy (L = 0 at top)
) -> Tuple[float, float]:
    """Calculate the sun and shade PAR values.

    Uses equations from Farquhar 1997 (Simple scaling of photosynthesis from leaves to canopies
    without the errors of big-leaf models)
    """
    k_b = 0.5 / sinB  # Beam radiation extinction coefficient
    k_b_alt = 0.46 / sinB  # beam and scattered beam PAR extinction coefficient
    P_h = (1 - (1 - sigma) ** 0.5) / (1 + (1 - sigma) ** 0.5)
    P_cb = 1 - exp(2 * P_h * k_b / (1 + k_b))
    P_cd = calc_diffuse_irradiance_refl(sinB, P_h, Ir_dfuse_0)
    Ir_diffuse = (1 - P_cd) * k_d_alt * Ir_dfuse_0 * exp(-k_d_alt * LAI_c)
    Ir_beam_sun = (1 - sigma) * Ir_beam_0 * cosA / sinB
    Ir_scattered_b = Ir_beam_0 * (
        ((1 - P_cb)*k_b_alt * exp(-k_b * LAI_c)) /
        -(1 - sigma) * k_b * exp(-k_b * LAI_c)
    )
    PAR_shade = Ir_diffuse + Ir_scattered_b
    PAR_sun = PAR_shade + Ir_beam_sun
    return PAR_sun, PAR_shade

#
# %%

# %%
# full_par_sun_shade_calc


@interact(
    sinB=widgets.FloatSlider(min=0.1, max=2.5, step=0.1, value=0.5),
    PAR=widgets.FloatSlider(min=0.01, max=700, step=100, value=500),
    P=fixed(91),
    cosA=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.5),
    LAI_c=widgets.FloatSlider(min=0.001, max=1, step=0.001, value=0.1),
)
def full_par_sun_shade_calc(
    sinB=0.5,
    PAR=700,
    P=91,
    cosA=0.5,
    LAI_c=0.1,
):
    Ir_beam_0, Ir_dfuse_0 = calc_Idrctt_Idfuse(PAR, sinB, P)
    sun_shade = calc_PAR_sun_shade(Ir_beam_0, Ir_dfuse_0, sinB, cosA, LAI_c)
    # sun_shade = [calc_PAR_sun_shade(Ir_beam_0, Ir_dfuse_0, sinB, cosA, LAI_c) for i in range(100)]
    return {'PARsun': sun_shade[0], 'PARshade': sun_shade[1]}

# %%

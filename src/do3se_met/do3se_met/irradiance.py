"""Met helpers that calcualte irradiance.

References
----------
 - DE PURY, D.G.G. and FARQUHAR, G.D. (1997), Simple scaling of photosynthesis from leaves
 to canopies without the errors of big‐leaf models. Plant, Cell & Environment
 - Tao et al (2012) - Spatiotemporal changes of wheat phenology in China under the effects of
temperature, day length and cultivar thermal characteristics

"""

from math import acos, cos, exp, pi, radians, sin, sqrt, asin, tan
from itertools import zip_longest
from typing import List, NamedTuple, Tuple
from collections import namedtuple
from deprecated import deprecated
from scipy.integrate import quad as integrate

from do3se_met.error_handling import InputError
from do3se_met.solar_position import solar_declination, solar_noon
from do3se_met.physical_constants import (
    GSC,
    PAR_Wm2_to_photons,
    PARfrac,
    # SBC,
    T0,
    Rn_MJ_to_W,
    seaP
)

from do3se_met.model_constants import MIN_DAYLIGHT_R


def calc_is_daylight(global_radiation: float) -> bool:
    """Calculate if global radiation is above specified amount.

    Is it currently daylight?  Uses the accepted criteria for accumulating
    OT40, i.e. when global radiation > 50.0 W m-2.

    Parameters
    ----------
    global_radiation: float
        Global radiation [W m^2]

    Returns
    -------
    is_daylight: bool
        True if radiation is above 50 W m-2 [bool]

    """
    return global_radiation > MIN_DAYLIGHT_R


def calc_photoperiod(dd: int,
                     lat: float) -> float:
    """Calculate the photoperiod using site latitude and Julian day number.

    Parameters
    ----------
    dd : int
        Julian day number(DOY)
    lat : float
        Site Latitude (degrees)

    Returns
    -------
    photoperiod: float
        Day length [???]

    """
    dec = asin(0.3978 * sin(((2 * pi) * (dd - 80) / 365) +  # noqa:W504
                    (0.0335 * (sin(2 * pi * dd) - sin(2 * pi * 80)) / 365)))
    perc =((-0.10453) / (cos(radians(lat))) * (cos(dec)))- tan(radians(lat)) * tan(dec)

    if perc < -1:
        return 24
    if perc > 1:
        return 0
    pr = 24 * (acos((perc))) / pi  # noqa:W503
    return pr

def calc_photoperiod_factor(photoperiod: float, PID: float) -> float:
    """Calculate the photoperiod factor.

    Between 0->1

    References
    ----------
    Tao et al(2012) - Eq 5.

    Parameters
    ----------
    photoperiod: float
        THe photoperiod from calc_photoperiod
    PID: float
        A cultivar specific parameter controling the photoperiod sensativity

    Returns
    -------
    photoperiod: float
        Between 0->1

    """
    return 1-((PID/10000)*(20-photoperiod)**2)


def MLMC_sunlit_LAI(
    nL: int,
    nLC: int,
    LAI: List[List[float]],
    sinB: float,
) -> float:
    """Multi-layer multi-component sunlit LAI model.

    LAI and LAIsunfrac must be the same dimension(nL,nLC).


    Parameters
    ----------
    nL: int
        Number of layers
    nLC: int
        Number of components
    LAI : List[List[float]]
        Leaf area index [m^2/m^2] with shape (nL, nLC)
    sinB : float
        sin() of solar elevation angle

    Returns
    -------
    List[List[float]]
        sunlit LAI  with shape (nL, nLC)

    """
    LAIsunfrac = [[None for i in range(nLC)] for iLC in range(nL)]
    sunLAI_acc = [None for i in range(nL + 1)]

    sunLAI_acc[0] = 0.0
    for iL in range(nL):
        # How much of "canopy so far" is sunlit?
        # TODO: should top layer be 1.0
        lai_c = sum([sum(LAI_iL) for LAI_iL in LAI[0:iL + 1]])
        sunLAI_acc[iL + 1] = sunlit_LAI(lai_c, sinB)
        # How much of that is in this layer?
        sunLAI_layer = sunLAI_acc[iL + 1] - sunLAI_acc[iL]
        # Fraction of LAI which is sunlit
        LAI_layer = sum(LAI[iL])
        if (LAI_layer > 0.0):
            LAIsunfrac[iL] = [sunLAI_layer / LAI_layer for i in range(nLC)]
        else:
            LAIsunfrac[iL] = [0 for i in range(nLC)]

    return LAIsunfrac


def sunlit_LAI(LAI: float, sinB: float) -> float:
    """Calculate the sunlit LAI from total LAI and solar elevation angle.

    # Farquhar 1997 Eq 18

    Parameters
    ----------
    LAI: float
        Leaf area index [m^2/m^2]
    sinB: float
        sin() of solar elevation angle

    """
    # TODO: (2 * sinB) should be sinB/cosA
    if (LAI > 0.0 and sinB > 0.0):
        return ((1 - exp(-0.5 * LAI / sinB)) * (2 * sinB))
    else:
        return 0.0


def calc_Idrctt_Idfuse(
    sinB: float,
    P: float,
    PAR: float = None,
    cloudFrac: float = None,
) -> float:
    """Estimate above canopy diffuse and direct PAR components.

    Must provide either cloud cover or PAR.

    Idrctt == PARdir in DO3SE_UI
    Idfuse == PARdif in DO3SE_UI

    Parameters
    ----------
    sinB: float
        sin() of solar elevation angle
    P: float
        Atmospheric pressure [kPa]
    PAR: float
        Photosynthetically active radiation [W m-2] default None
    cloudFrac: float
        cloud cover fraction. default none

    Returns
    -------
    Idrctt: float
        Direct PAR irradiance [W m-2]
    Idfuse: float
        Diffuse PAR irradiance [W m-2]

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

        # Sky transmissivity from PAR or cloud frac
        if PAR is not None:
            ST = max(0.21, min(0.9, PAR / pPARtotal))
        elif cloudFrac is not None:
            ST = 1.0 - 0.75 * (cloudFrac ** 3.4)
        else:
            raise ValueError("Must supply PAR or cloudFrac")

        # Direct and diffuse fractions
        # A = 0.9
        # B = 0.7
        if ST < 0.9:
            fPARdir = (pPARdir / pPARtotal) * (1 - ((0.9 - ST) / 0.7)**(2.0 / 3.0))
        else:
            fPARdir = (pPARdir / pPARtotal)
        fPARdif = 1 - fPARdir

        # Apply calculated direct and diffuse fractions to PARtotal
        PAR_out = PAR if PAR is not None else ST * pPARtotal
        Idrctt = fPARdir * PAR_out
        Idfuse = fPARdif * PAR_out

    else:
        Idrctt = 0.0
        Idfuse = 0.0
        PAR_out = PAR or 0
    assert PAR_out is not None, "Could not define PAR"
    return Idrctt, Idfuse, PAR_out


def check_required_inputs(
    PAR_in: float = None,
    Idrctt_in: float = None,
    Idfuse_in: float = None,
    PPFD_in: float = None,
    R_in: float = None,
    Rn_in: float = None,
    sinB: float = None,
    P: float = None,
    cloudfrac: float = None,
) -> str:
    if PAR is None:
        Idrctt, Idfuse, PAR = calc_Idrctt_Idfuse(sinB, P, PAR=PAR_in, cloudFrac=cloudfrac) if (cloudfrac is not None and sinB is not None and P is not None) \
            else [Idrctt_in, Idfuse_in, Idrctt_in + Idfuse_in] if (Idrctt_in is not None and Idfuse_in is not None) \
            else [None, None, PPFD_in / PAR_Wm2_to_photons] if PPFD_in is not None \
            else [None, None, R_in * PARfrac] if R_in is not None \
            else [None, None, Rn_in * Rn_MJ_to_W * PARfrac] if Rn_in is not None \
            else [None, None, None]
        if PAR is None:
            raise InputError("PAR", "Could not define PAR check inputs")
    if Idfuse is None or Idrctt is None:
        Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
        if Idrctt is None:
            raise InputError("Idrctt", "Could not define Idrctt check inputs")
        if Idfuse is None:
            raise InputError("Idfuse", "Could not define Idfuse check inputs")



def calc_radiation(
    PAR_in: float = None,
    Idrctt_in: float = None,
    Idfuse_in: float = None,
    PPFD_in: float = None,
    R_in: float = None,
    Rn_in: float = None,
    sinB: float = None,
    P: float = None,
    cloudfrac: float = None,
) -> NamedTuple:
    """Calculate missing radiation data from provided data.

    pure subroutine met_radiation from Fortran model

    Parameters
    ----------
    PAR_in : float, optional
        PAR, by default None [umol/m^2/s]
    Idrctt_in : float, optional
        [description], by default None [UNIT]
    Idfuse_in : float, optional
        [description], by default None [UNIT]
    PPFD_in : float, optional
        [description], by default None [UNIT]
    R_in : float, optional
        [description], by default None [UNIT]
    Rn_in : float, optional
        Net radiation, by default None [Mj/m^2/h^1]
    sinB : float, optional
        [description], by default None [UNIT]
    P : float, optional
        [description], by default None [UNIT]
    cloudfrac: float, optional
        Cloud fraction [Fraction]

    Returns
    -------
    NamedTuple
        PAR PPFD Idrctt Idfuse R Rn

    """
    Output = namedtuple('Output', 'PAR PPFD Idrctt Idfuse R Rn')

    Idrctt = Idrctt_in
    Idfuse = Idfuse_in
    PAR = PAR_in
    PPFD = PPFD_in
    R = R_in
    Rn = Rn_in

    if PAR is None:
        Idrctt, Idfuse, PAR = calc_Idrctt_Idfuse(sinB, P, PAR=PAR_in, cloudFrac=cloudfrac) if (cloudfrac is not None and sinB is not None and P is not None) \
            else [Idrctt_in, Idfuse_in, Idrctt_in + Idfuse_in] if (Idrctt_in is not None and Idfuse_in is not None) \
            else [None, None, PPFD_in / PAR_Wm2_to_photons] if PPFD_in is not None \
            else [None, None, R_in * PARfrac] if R_in is not None \
            else [None, None, Rn_in * Rn_MJ_to_W * PARfrac] if Rn_in is not None \
            else [None, None, PAR_in]
        if PAR is None:
            raise InputError("PAR", "Could not define PAR check inputs. Must supply ")
    if Idfuse is None or Idrctt is None:
        Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
        if Idrctt is None:
            raise InputError("Idrctt", "Could not define Idrctt check inputs")
        if Idfuse is None:
            raise InputError("Idfuse", "Could not define Idfuse check inputs")

    if PPFD is None:
        PPFD = PAR * PAR_Wm2_to_photons

    if R is None:
        R = PAR / PARfrac

    if Rn is None:
        Rn = R / Rn_MJ_to_W

    return Output(
        Idrctt = Idrctt,
        Idfuse = Idfuse,
        PAR = PAR,
        PPFD = PPFD,
        R = R,
        Rn = Rn,
    )


def calc_radiation_list(
    PAR_list: List[float] = [],
    Idrctt_list: List[float] = [],
    Idfuse_list: List[float] = [],
    PPFD_list: List[float] = [],
    R_list: List[float] = [],
    Rn_list: List[float] = [],
    sinB_list: List[float] = [],
    P_list: List[float] = [],
    cloudfrac_list: List[float] = [],
) -> NamedTuple:
    """Run calc_radiation on list of data.

    Parameters
    ----------
    PAR_list : List[float], optional
        [description], by default None [UNIT]
    Idrctt_list : List[float], optional
        [description], by default None [UNIT]
    Idfuse_list : List[float], optional
        [description], by default None [UNIT]
    PPFD_list : List[float], optional
        [description], by default None [UNIT]
    R_list : List[float], optional
        [description], by default None [UNIT]
    R_list : List[float], optional
        [description], by default None [UNIT]
    sinB_list : List[float], optional
        [description], by default None [UNIT]
    P_list : List[float], optional
        [description], by default None [UNIT]
    cloudfrac_list: List[float], optional
        Cloud fraction [Fraction]

    Returns
    -------
    NamedTuple
        [description]

    """
    Output = namedtuple('Output', 'PAR Idrctt Idfuse PPFD R Rn cloudfrac')
    PAR_out = []
    Idrctt_out = []
    Idfuse_out = []
    PPFD_out = []
    R_out = []
    Rn_out = []
    _PAR_list = PAR_list if PAR_list is not None else []
    _Idrctt_list = Idrctt_list if Idrctt_list is not None else []
    _Idfuse_list = Idfuse_list if Idfuse_list is not None else []
    _PPFD_list = PPFD_list if PPFD_list is not None else []
    _R_list = R_list if R_list is not None else []
    _Rn_list = Rn_list if Rn_list is not None else []
    _sinB_list = sinB_list if sinB_list is not None else []
    _P_list = P_list if P_list is not None else []
    _cloudfrac_list = cloudfrac_list if cloudfrac_list is not None else []


    for PAR_in, Idfuse_in, Idfuse_in, PPFD_in, R, Rn, sinB, P, cloudfrac \
            in zip_longest(_PAR_list, _Idrctt_list, _Idfuse_list, _PPFD_list, _R_list, _Rn_list, _sinB_list, _P_list, _cloudfrac_list):
        out = calc_radiation(PAR_in, Idfuse_in, Idfuse_in, PPFD_in, R, Rn, sinB, P, cloudfrac)
        PAR_out.append(out.PAR)
        Idrctt_out.append(out.Idrctt)
        Idfuse_out.append(out.Idfuse)
        PPFD_out.append(out.PPFD)
        R_out.append(out.R)
        Rn_out.append(out.Rn)

    return Output(
        PAR=PAR_out,
        Idrctt=Idrctt_out,
        Idfuse=Idfuse_out,
        PPFD=PPFD_out,
        R=R_out,
        Rn=Rn_out,
        cloudfrac=cloudfrac,
    )


def calc_net_radiation(
    lat: float,
    lon: float,
    elev: float,
    albedo: float,
    dd: int,
    hr: int,
    sinB: float,
    R: float,
    Ts_C: float,
    eact: float,
) -> float:
    """Estimate net radiation [MJ m-2 h-1].

    Parameters
    ----------
    lat: float
        Latitude [degrees North]
    lon: float
        Longitude [degrees East]
    elev: float
        Elevation [m above sea level]
    albedo: float
        Surface albedo [fraction]
    dd: int
        Day of year [1--365]
    hr: int
        Hour of day [0--23]
    sinB: float
        sin() of solar elevation angle
    R: float
        Global radiation [W m-2]
    Ts_C: float
        Surface air temperature [degrees C]
    eact: float
        Actual vapour pressure [kPa]

    Returns
    -------
    net_radiation: float
        [description][Unit]

    """
    # Where s is the Stefan-Boltzman constant (4.903 10-9 MJ K-4 m-2 h-1, Tk
    # is air temperature (°K), ea is the actual vapour pressure (kPa), R is the
    # actual global radiation (MJ m-2 h-1) and pR is the potential global
    # radiation (MJ m-2 h-1).
    SBC = 4.903e-9 / 24  # ! Stephan Boltzman constant

    if sinB <= 0:
        net_radiation = 0.0
    else:
        # Latitude in radians
        lat_rad = radians(lat)

        # Convert global radiation W m-2 to MJ m-2 s-1
        R_MJ = R * 0.0036

        # Hour-angle of the sun
        t0_ = solar_noon(lon, dd)
        h = radians(15 * (hr - t0_))
        h1 = h - (pi / 24)
        h2 = h + (pi / 24)

        dr = 1 + (0.033 * cos(((2 * pi) / 365) * dd))
        dec = solar_declination(dd)
        # External radiation (with fix to stop div by zero)
        # TODO: fix this to be less hackish
        Re = max(0.00000000001,
                 ((12 * 60) / pi) * GSC * dr * ((h2 - h1) * sin(lat_rad) * sin(dec) + cos(lat_rad) * cos(dec) * (sin(h2) - sin(h1))))  # noqa:E501
        # TODO: what was this for?
        # Re = max(0.0, ((12*60)/pi)*Gsc*dr*sinB)

        # Calculate net longwave radiation
        assert elev is not None
        assert Re is not None
        pR = (0.75 + (2e-5 * elev)) * Re

        Rnl = max(0.0, (SBC * ((Ts_C + T0)**4)) * (0.34 - (0.14 * sqrt(eact)))
                  * ((1.35 * (min(1.0, R_MJ / pR))) - 0.35))  # noqa:W503
        Rns = (1 - albedo) * R_MJ

        net_radiation = max(0.0, Rns - Rnl)
    return net_radiation


def get_net_radiation(
    Rn_in: float = None,
    lat: float = None,
    lon: float = None,
    elev: float = None,
    albedo: float = None,
    dd: int = None,
    hr: int = None,
    sinB: float = None,
    R: float = None,
    Ts_C: float = None,
    eact: float = None,
) -> float:
    """Get the net radiation (Rn) using calc_net_radiation if Rn is None.

    Parameters
    ----------
    Rn_in : float, optional
        [description], by default None [UNIT]
    lat : float, optional
        [description], by default None [UNIT]
    lon : float, optional
        [description], by default None [UNIT]
    elev : float, optional
        [description], by default None [UNIT]
    albedo : float, optional
        [description], by default None [UNIT]
    dd : int, optional
        [description], by default None [UNIT]
    hr : int, optional
        [description], by default None [UNIT]
    sinB : float, optional
        [description], by default None [UNIT]
    R : float, optional
        [description], by default None [UNIT]
    Ts_C : float, optional
        [description], by default None [UNIT]
    eact : float, optional
        [description], by default None [UNIT]

    Returns
    -------
    float
        [description]

    """
    Rn = Rn_in if Rn_in is not None \
        else calc_net_radiation(lat, lon, elev, albedo,
                                dd, hr, sinB, R, Ts_C, eact)
    return Rn


def calc_PAR_sun_shade_multilayer(
    Idrctt: float,
    Idfuse: float,
    sinB: float,
    cosA: List[float],
    LAI: List[float],
) -> Tuple[float, float]:
    """Calculate PAR sun and shade for multilayer model.

    Calculates average values.

    # TODO: Check if this needs replacing for actual multilayer model.

    Parameters
    ----------
    Idrctt: float
        Direct PAR irradiance [W m-2]
    Idfuse: float
        Diffuse PAR irradiance [W m-2]
    sinB: float
        sin() of solar elevation angle
    cosA: float
        cos(A), A = mean leaf inclination (0.5 = 60 degrees) [radians]
    LAI: float
        Leaf area index [m2 m-2]

    Returns
    -------
    Tuple[float, float]
        PARsun, PARshade

    """
    total_layer_LAI = sum(LAI)
    average_cosA = sum(cosA) / len(cosA)
    return calc_PAR_sun_shade(Idrctt, Idfuse, sinB, average_cosA, total_layer_LAI)


def calc_PAR_sun_shade(
    Idrctt: float,
    Idfuse: float,
    sinB: float,
    cosA: float,
    LAI: float,
) -> Tuple[float, float]:
    """Estimate PAR received by sun and shade leaves within the canopy.

    ! TODO: multi-layer PAR
    Taken from Fortran DO3SE-model

    Parameters
    ----------
    Idrctt: float
        Direct PAR irradiance [W m-2]
    Idfuse: float
        Diffuse PAR irradiance [W m-2]
    sinB: float
        sin() of solar elevation angle
    cosA: float
        cos(A), A = mean leaf inclination (0.5 = 60 degrees) [radians]
    LAI: float
        Leaf area index [m2 m-2]

    Returns
    -------
    NamedTuple containing:
    PARsun: float
        PAR received by sunlit leaves [W m-2]
    PARshade: float
        PAR received by shaded leaves [W m-2]

    """
    Output = namedtuple('Output', 'PARsun PARshade')
    if sinB > 0.0:
        # PAR flux densities evaluated using method of Norman(1982, p.79):
        # "conceptually, 0.07 represents a scattering coefficient"

        PARshade = Idfuse * exp(-0.5 * LAI**0.8) + \
            0.07 * Idrctt * (1.1 - (0.1 * LAI)) * exp(-sinB)
        PARsun = Idrctt * 0.8 * (cosA / sinB) + PARshade
    else:
        PARshade = 0.0
        PARsun = 0.0
    out = Output(PARsun, PARshade)
    return out


# Equations below taken from Farquhar 1997
# Model coefficients
# TODO: Move these to config or constants
# sigma = 0.15  # leaf scattering coefficient of PAR (p_i + T_i)
k_d_alt = 0.719  # diffuse and scattered diffuse PAR extinction coefficient, 0.719
# seaP = 101.325  # real, parameter


def calc_beam_irradiance_horiz(
    sigma: float = 0.15,  # leaf scattering coefficient of PAR (p_i + T_i)
) -> float:
    """Calculate the beam irradiance for horiontal leaves.

    Farqhuar 1997 Eq A19

    """
    P_h = (1 - (1 - sigma) ** 0.5) / (1 + (1 - sigma) ** 0.5)
    return P_h


def calc_beam_irradiance_uad(
    P_h: float,  # reflection coefficient of a canopy with horizontal leaves
    k_b: float,  # Beam radiation extinction coefficient (0.5 / sinB)
) -> float:
    """Calculate beam irradiance for uniform leaf angle distribution.
    As a function of solar elevation

    Farquhar 1997 Eq A19.

    """
    P_cb = 1 - exp(-2 * P_h * k_b / (1 + k_b))  # Note missing '-' in paper after exp
    return P_cb


def calc_diffuse_irradiance_refl(
    # P_cb: float,  # beam irradiance for uniform leaf angle distribution
    sinB: float,
    P_h: float,
    # N_d: float,  #
    Ir_dfuse_0: float,  # Diffuse PAR per unit ground area at top of canopy[umol m^-2 s^-2]
) -> float:
    """Calculate diffuse irradiance reflection coefficient EQ A21.

    'In the case of a uniform leaf angle distribution, with leaf scattering
    coeff(sigma) = 0.15, P_h = 0.041 so that with a uniform sky radiance the canopy
    reflection coefficient for diffuse irradiance(P_cd) is calculated as 0.036

    # TODO: This needs checking

    """
    # f
    def f(alpha):
        """Calculate integral per rotation unit of sun[radian]"""
        # TODO: Check this is correct
        # We should be
        k_b = 0.5 / sinB
        # TODO: Ir_dfuse_0 and P_cb should be recalculated here as a function of alpha
        # Diffuse photon radiance of the sky (per radian??)[umol m^-2 sr ^-1]
        N_d = Ir_dfuse_0 / (2 * pi)
        P_cb = calc_beam_irradiance_uad(P_h, k_b)
        return N_d * P_cb

    integ, err = integrate(f, 0, pi / 2)
    P_cd = (1 / Ir_dfuse_0) * integ
    return P_cd


def calc_Ir_scattered_b(
    P_cb: float,  # beam irradiance for uniform leaf angle distribution
    Ir_beam_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    LAI_c: float,  # cumulative leaf-area index from top of canopy (L = 0 at top)
    k_b: float,  # Beam radiation extinction coefficient
    k_b_alt: float,  # beam and scattered beam PAR extinction coefficient, 0-46/sinj3
    sigma: float,  # leaf scattering coefficient of PAR (pi + Ti)
) -> float:
    """Calculate the scattered beam irradiance.

    Farquhar 1997 EqA8

    Parameters
    ----------
    P_cb: float
        beam irradiance for uniform leaf angle distribution
    Ir_beam_0: float
        PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    LAI_c: float
        cumulative leaf-area index from top of canopy (L = 0 at top)
    k_b: float
        Beam radiation extinction coefficient
    k_b_alt: float
        beam and scattered beam PAR extinction coefficient, 0-46/sinj3
    sigma: float
        leaf scattering coefficient of PAR (pi + Ti)

    Returns
    -------
    Ir_bs: float
        [description][unit]

    """
    Ir_bs = Ir_beam_0 * (
        ((1 - P_cb) * k_b_alt * exp(-k_b * LAI_c))
        - (1 - sigma) * k_b * exp(-k_b * LAI_c)  # noqa:W503
    )
    return Ir_bs


def calc_Ir_beam_sun(
    sinB: float,  # sine of solar elevation angle [radians]
    cosA: float,  # cosine of angle of beam irradiance to the leaf normal [radians]
    Ir_beam_0: float,  # beam PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    sigma: float = 0.15,  # leaf scattering coefficient of PAR
) -> float:
    """Eq A11."""
    Ir_b = (1 - sigma) * Ir_beam_0 * cosA / sinB
    return Ir_b


def calc_Ir_diffuse(
    P_cd: float,  # diffuse irradiance reflection coefficient
    Ir_dfuse_0: float,  # diffuse PAR per unit ground area at top of canopy
    LAI_c: float,  # Cumulative leaf area index from top of canopy
    k_d_alt: float = k_d_alt,  # Diffuse and scattered diffuse PAR extinction coefficient
) -> float:
    """Diffuse PAR per unit ground area Eq A5."""
    ir_diffuse = (1 - P_cd) * k_d_alt * Ir_dfuse_0 * exp(-k_d_alt * LAI_c)
    return ir_diffuse


def calc_PAR_shade(
    Ir_diffuse: float,  # diffuse PAR per unit ground area
    Ir_scattered_b: float,  # absorbed scattered beam PAR per unit leaf area
) -> float:
    """Eq A7."""
    return Ir_diffuse + Ir_scattered_b


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


def calc_PAR_sun_shade_farq(
    Ir_beam_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    Ir_dfuse_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    sinB: float,
    cosA: float,
    LAI_c: float,  # cumulative leaf-area index from top of canopy (L = 0 at top)
    sigma: float = 0.15,
) -> Tuple[float, float]:
    """Calculate the sun and shade PAR values.

    Uses equations from Farquhar 1997 (Simple scaling of photosynthesis from leaves to canopies
    without the errors of big-leaf models)

    """
    k_b = 0.5 / sinB  # Beam radiation extinction coefficient
    k_b_alt = 0.46 / sinB  # beam and scattered beam PAR extinction coefficient
    P_h = calc_beam_irradiance_horiz(sigma)
    P_cb = calc_beam_irradiance_uad(P_h, k_b)
    P_cd = calc_diffuse_irradiance_refl(sinB, P_h, Ir_dfuse_0)
    Ir_diffuse = calc_Ir_diffuse(P_cd, Ir_dfuse_0, LAI_c, k_d_alt)
    Ir_beam_sun = calc_Ir_beam_sun(sinB, cosA, Ir_beam_0, sigma)
    Ir_scattered_b = calc_Ir_scattered_b(P_cb, Ir_beam_0, LAI_c, k_b, k_b_alt, sigma)
    PAR_shade = calc_PAR_shade(Ir_diffuse, Ir_scattered_b)
    PAR_sun = calc_PAR_sun(PAR_shade, Ir_beam_sun)

    return PAR_sun, PAR_shade


def calc_PAR_sun_shade_farq_b(
    Ir_beam_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    Ir_dfuse_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    sinB: float,
    cosA: float,
    LAI_c: float,  # cumulative leaf-area index from top of canopy (LAI = 0 at top)
    sigma: float = 0.15,
) -> Tuple[float, float]:
    """Calculate the sun and shade PAR values.

    This is the same as calc_PAR_sun_shade but combined all the equations into a single function.
    This makes it more efficient to run

    Uses equations from Farquhar 1997 (Simple scaling of photosynthesis from leaves to canopies
    without the errors of big-leaf models)

    """
    Output = namedtuple('Output', 'PARsun, PARshade')
    if sinB <= 0:
        return Output(0, 0)
    if Ir_dfuse_0 + Ir_beam_0 < 0.001:
        return Output(0, 0)

    k_b = 0.5 / sinB  # Beam radiation extinction coefficient
    k_b_alt = 0.46 / sinB  # beam and scattered beam PAR extinction coefficient
    P_h = (1 - (1 - sigma) ** 0.5) / (1 + (1 - sigma) ** 0.5)
    P_cb = 1 - exp(-2 * P_h * k_b / (1 + k_b))
    P_cd = calc_diffuse_irradiance_refl(sinB, P_h, Ir_dfuse_0)
    Ir_diffuse = (1 - P_cd) * k_d_alt * Ir_dfuse_0 * exp(-k_d_alt * LAI_c)
    Ir_beam_sun = (1 - sigma) * Ir_beam_0 * cosA / sinB
    Ir_scattered_b = Ir_beam_0 * (
        ((1 - P_cb) * k_b_alt * exp(-k_b * LAI_c))
        - (1 - sigma) * k_b * exp(-k_b * LAI_c)  # noqa:W503
    )
    PAR_shade = Ir_diffuse + Ir_scattered_b
    PAR_sun = PAR_shade + Ir_beam_sun
    return Output(PAR_sun, PAR_shade)


@deprecated(version="0", reason="Only required to match old UI output. Use get_parSunShade_from_par instead.")
def calc_PAR_sun_shade_UI(
    PAR: float,
    sinB: float,
    P: float,
    cosA: float,
    LAI: float,
) -> Tuple[float, float]:
    """Estimate PAR received by sun and shade leaves within the canopy using UI method.

    Method taken from DO3SE UI fortran model.

    Parameters
    ----------
    PAR: float
        Input PAR data [W m-2]
    sinb: float
        Solar elevation [radians]
    P: float
        Air Pressure [kPa]
    cosA: float
        Leaf angle[radians]
    LAI: float
        Leaf area

    Returns
    -------
    NamedTuple containing:
    PARsun: float
        PAR received by sunlit leaves [W m-2]
    PARshade: float
        PAR received by shaded leaves [W m-2]

    """
    Output = namedtuple('Output', 'PARsun PARshade')
    if sinB <= 0:
        return Output(0, 0)
    Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
    C = 4.57
    PARshade = Idfuse * C * exp(-0.5 * (LAI**0.8)) + 0.07 * Idrctt * C * \
        (1.1 - (0.1 * LAI)) * exp(-sinB)
    PARsun = Idrctt * C * 0.8 * (cosA / sinB) + PARshade
    return Output(PARsun / C, PARshade / C)


def get_parSunShade_from_par(PAR, sinB, P, cosA, LAI):
    """Calculate PAR sun shade to match DO3SE MODEL output.

        Parameters
    ----------
    PAR: float
        Input PAR data [W m-2]
    sinb: float
        Solar elevation [radians]
    P: float
        Air Pressure [kPa]
    cosA: float
        Leaf angle[radians]
    LAI: float
        Leaf area

    Returns
    -------
    NamedTuple containing:
    PARsun: float
        PAR received by sunlit leaves [W m-2]
    PARshade: float
        PAR received by shaded leaves [W m-2]

    """
    Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
    return calc_PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, LAI)

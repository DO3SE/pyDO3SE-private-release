"""Functions and data models used by the Ewert model.

Key interactive plots:
f_LS, f_LA and fO3_l:
https://www.desmos.com/calculator/h66l70w4wq

Note: There are some variations from the original Ewert paper that are documented
in the individual function docstrings below.

References
----------

- Ewert, F., Porter, J.R., 2000. Ozone effects on wheat in relation to CO2: modelling short-term and long-term responses of leaf photosynthesis and leaf duration. Global Change Biology 6, 735–750. https://doi.org/10.1046/j.1365-2486.2000.00351.x
- Nikolov, N.T., Massman, W.J., Schoettle, A.W., 1995. Coupling biochemical and biophysical processes at the leaf level: an equilibrium photosynthesis model for leaves of C3 plants. Ecological Modelling 80, 205–235. https://doi.org/10.1016/0304-3800(94)00072-P
- LEUNING, R., 1995. A critical appraisal of a combined stomatal-photosynthesis model for C3 plants. Plant, Cell & Environment 18, 339–355. https://doi.org/10.1111/j.1365-3040.1995.tb00370.x
- Wang, Y.-P., Leuning, R., 1998. A two-leaf model for canopy conductance, photosynthesis and partitioning of available energy I:: Model description and comparison with a multi-layered model. Agricultural and Forest Meteorology 91, 89–111. https://doi.org/10.1016/S0168-1923(98)00061-6
- von Caemmerer, S., Farquhar, G.D., 1981. Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta 153, 376–387. https://doi.org/10.1007/BF00384257
- D. B. Clark et al.: JULES: carbon fluxes and vegetation dynamics
- Feng, Y. et al. (2022) ‘Identifying and modelling key physiological traits that confer tolerance or sensitivity to ozone in winter wheat’, Environmental Pollution. Elsevier Ltd, 304(April), p. 119251. doi: 10.1016/j.envpol.2022.119251.

"""

from typing import NamedTuple, List
from math import sqrt
import numpy as np
from pyDO3SE.constants.physical_constants import R
from pyDO3SE.util.error_handling import ConfigError
from pyDO3SE.Config.ConfigEnums import FVPDMethods
from do3se_met.conversion import deg_to_kel
from do3se_met.helpers import saturated_vapour_pressure
from pyDO3SE.plugins.gsto.helpers import temp_dep_inhibit, temp_dep
from pyDO3SE.plugins.gsto.constants import (
    A_j_a,
    A_j_b,
    Gamma_star_25,
    E_Gamma_star,
    E_K_C,
    K_C_25,
    K_O_25,
    E_K_O,
    H_a_jmax,
    H_d_jmax,
    S_V_jmax,
    H_a_vcmax,
    H_d_vcmax,
    S_V_vcmax,
    alpha,
    Teta,
    O_i,
)


def calc_A_n_product_limited(V_cmax: float, R_d: float) -> float:
    """Calculate the product-limited photosynthesis rate

    Parameters
    ----------
    V_cmax: float
        Maximum carboxylation rate [umol/m^2/s]
    R_d: float
        Dark respiration rate [umol/m^2/s CO2]

    Returns
    -------
    A_p
        Product-limited Photosythesis Rate [umol/m^2/s]

    References
    ----------
    - Ewert, F., Porter, J.R., 2000 - Eq8


    """
    return (0.5 * V_cmax) - R_d


def calc_A_n_rubisco_limited(
    O2: float,
    P: float,
    K_O: float,
    V_cmax: float,
    R_d: float,
    g_sto_0: float,
    G_1c: float,
    g_bl: float,
    c_a: float,
    RH: float,
    Gamma: float,
    Gamma_star: float,
    K_C: float,
    fO3_d: float,  # we need to include the ozone factors
    f_LS: float,  # need the senescene factor for Ac aswell
    f_SW: float,  # need the soil water factor to calculate Ac
    negative_A_n: bool,
):
    """Calculate the rubisco-limited photosynthesis rate

    Based on the following papers:
    (Masutomi, Y., 2023)

    With the following modifications:
    - The equation has been modified to include the CO2 compensation point using


    Parameters
    ----------
    O2: float
        Partial pressure of atmospheric O2 [Pa]
    P: float
        Surface Air Pressure [kPa]
    K_O: float
        Michaelis-Menten constant for O2 [mmol/mol]
    V_cmax: float
        Maximum carboxylation rate [umol/m^2/s]
    R_d: float
        Dark respiration rate [umol/m^2/s CO2]
    g_sto_0: float
        See gsto_0 [umol/m^2/s CO2]
    G_1c: 5.625
        the species specific sensitivity of stomatal conductance to VPD and stomatal conductance [dimensionless]
    g_bl: float
        Boundary Layer Conductance [mmol/m^2/s H2O]
    c_a: float
        Atmospheric CO2 concentration [ppm]
    RH: float
        Relative Humidity []
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    Gamma_star: float
        CO2 compensation point in the absence of respiration [umol/mol CO2]
    K_C: float
        Michaelis constant CO2 [ppm]
    C_s: float
        CO2 concentration at leaf surface [ppm]
    J: float
        Rate of electron transport [umol/m^2/s]
    fO3_d: float
        Hourly accumulated ozone impact factor [dimensionless][0-1]
        Cumulative ozone effect [dimensionless]
    f_LS: float
        factor effect of leaf senescence on A_c [Dimensionless][0-1]
    f_SW: float
        soil water effect on gsto [dimenstionless] [0-1]
    negative_A_n: boolean
        Boolean variable to select a negative A_n

    Returns
    -------
    A_c
        Rubisco-Limited Photosythesis Rate [umol/m^2/s]


    References
    ----------

    - Masutomi, Y., 2023. The appropriate analytical solution for coupled leaf photosynthesis
    - Baldocchi D. An analytical solution for coupled leaf photosynthesis and stomatal conductance models. Tree Physiol. 1994;14(7_9):1069-1079. doi:10.1093/treephys/14.7-8-9.1069
    and stomatal conductance models for C3 plants. Ecological Modelling, 481, p.110306.

    """
    # Converts Pressure from kPa (in DO3SE code) to Pa (in this code)
    P = 1000 * P

    #Convert gsto0 from micro moles m-2 PLA s-1 of H2O to moles m-2 PLA s-1 of CO2
    g_bl_mol_co2=g_bl*1e-6/1.37

    #Convert gsto0 from micro moles m-2 PLA s-1 of CO2 to moles m-2 PLA s-1 of CO2
    g_sto_0_mol=g_sto_0*1e-6

    # Define Factors for cubic Equations
    alpha = g_bl_mol_co2 * c_a
    beta = (g_bl_mol_co2 * G_1c * RH) - g_sto_0_mol
    gamma_c = (
        (V_cmax * fO3_d * f_LS * f_SW * Gamma_star)
        + (K_C * (1 + ((1000 * O2 / P) / K_O)) * R_d)
    )  # in our Ac and Aj we use compensation point in the absence of respiration so need Gamma_star
    zeta_c = V_cmax * fO3_d * f_LS * f_SW - R_d
    eta_c = (alpha * zeta_c) - (g_bl_mol_co2 * gamma_c)
    xi_c = alpha + zeta_c + (K_C * (1 + (1000 * O2 / P) / K_O) * g_bl_mol_co2)

    lambda_c = f_LS * f_SW * fO3_d * V_cmax
    b_with_comp_point = -g_bl_mol_co2 * (g_sto_0_mol + g_bl_mol_co2) * Gamma
    c_with_comp_point = (
        (1 / K_O)
        * g_bl_mol_co2
        * Gamma
        * (
            c_a * g_sto_0_mol * g_bl_mol_co2 * K_O
            + g_sto_0_mol * K_O * (-R_d + lambda_c)
            + g_bl_mol_co2 * (g_sto_0_mol * K_C * (K_O + (1000 * O2 / P)) - K_O * R_d + K_O * lambda_c)
        )
    )
    d_with_comp_point = (
        g_sto_0_mol
        * (g_bl_mol_co2**2)
        * Gamma
        * (
            K_O * K_C * R_d
            + K_C * (1000 * O2 / P) * R_d
            + lambda_c * K_O * Gamma_star
            + c_a * K_O * (R_d - lambda_c)
        )
    ) / K_O

    # Define Cubic Equation Coefficients
    if negative_A_n == False:
        a = beta - g_bl_mol_co2
        b = (
            (g_sto_0_mol * alpha)
            - (beta * xi_c)
            + (g_bl_mol_co2 * alpha)
            + (g_bl_mol_co2 * zeta_c)
            + b_with_comp_point
        )
        c = (
            -(g_sto_0_mol * alpha * xi_c)
            + (beta * eta_c)
            - (g_bl_mol_co2 * alpha * zeta_c)
            + c_with_comp_point
        )
        d = g_sto_0_mol * alpha * eta_c + d_with_comp_point
    elif negative_A_n == True:
        a = 0
        b = g_sto_0_mol + g_bl_mol_co2
        c = (g_sto_0_mol * xi_c) + (g_bl_mol_co2 * zeta_c)
        d = g_sto_0_mol * eta_c

    # Solve Cubic Equation
    # roots = np.roots([a, b, c, d])[np.isreal(np.roots([a, b, c, d]))]
    # roots = [np.real(x) for x in roots]
    roots = np.roots([a, b, c, d])
    roots = [np.real(x) for x in roots if (abs(x.imag)<1e-10)]
    return roots


def calc_A_n_rubp_limited(
    P: float,
    R_d: float,
    g_sto_0: float,
    G_1c: float,
    g_bl: float,
    c_a: float,
    RH: float,
    Gamma: float,
    Gamma_star: float,
    J: float,
    negative_A_n: bool,
):
    """Calculate the RuBP-limited photosynthesis rate


    Based on the following papers:
    (Masutomi, Y., 2023)

    With the following modifications:
    - The equation has been modified to include the CO2 compensation point using


    Parameters
    ----------
    O2: float
        Partial pressure of atmospheric O2 [Pa]
    P: float
        Surface Air Pressure [kPa]
    K_O: float
        Michaelis-Menten constant for O2 [mmol/mol]
    V_cmax: float
        Maximum carboxylation rate [umol/m^2/s]
    R_d: float
        Dark respiration rate [umol/m^2/s CO2]
    g_sto_0: float
        the minimum stomatal conductance [umol m-2 s-1 CO2]
    G_1c: float
       the species specific sensitivity of stomatal conductance to VPD and stomatal conductance [dimensionless]
    g_bl: float
        Boundary Layer Conductance [mmol/m^2/s H2O]
    c_a: float
        Atmospheric CO2 concentration [ppm]
    RH: float
        Relative Humidity []
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    K_C: float
        Michaelis constant CO2 [ppm]
    C_s: float
        CO2 concentration at leaf surface [ppm]
    J: float
        Rate of electron transport [umol/m^2/s]
    negative_A_n: boolean
        Boolean variable to select a negative A_n

    Returns
    -------
    A_j
        RuBP-Limited Photosythesis Rate [umol/m^2/s]

    References
    ----------
    - Masutomi, Y., 2023. The appropriate analytical solution for coupled leaf photosynthesis
    - Baldocchi D. An analytical solution for coupled leaf photosynthesis and stomatal conductance models. Tree Physiol. 1994;14(7_9):1069-1079. doi:10.1093/treephys/14.7-8-9.1069
    and stomatal conductance models for C3 plants. Ecological Modelling, 481, p.110306.

    """
    # Converts Pressure from kPa (in DO3SE code) to Pa (in this code)
    P = 1000 * P


    #Convert gsto0 from micro moles m-2 PLA s-1 of H2O to moles m-2 PLA s-1 of CO2
    g_bl_mol_co2=g_bl*1e-6/1.37

    #Convert gsto0 from micro moles m-2 PLA s-1 of CO2 to moles m-2 PLA s-1 of CO2
    g_sto_0_mol=g_sto_0*1e-6

    # Define Factors for cubic Equations
    alpha = g_bl_mol_co2 * c_a
    beta = (g_bl_mol_co2 * G_1c * RH) - g_sto_0_mol
    gamma_j = (J * Gamma_star) + (8 * Gamma_star * R_d)
    zeta_j = J - (A_j_a * R_d)
    eta_j = (alpha * zeta_j) - (g_bl_mol_co2 * gamma_j)
    xi_j = (A_j_a * alpha) + zeta_j + (A_j_b * Gamma_star * g_bl_mol_co2)

    # define extra terms for cubic coefficients
    b_with_comp_point = -A_j_a * g_bl_mol_co2 * (g_sto_0_mol + g_bl_mol_co2) * Gamma
    c_with_comp_point = (
        g_bl_mol_co2
        * Gamma
        * (
            A_j_a * c_a * g_sto_0_mol * g_bl_mol_co2
            + g_sto_0_mol * (J - A_j_a * R_d)
            + g_bl_mol_co2 * (J - A_j_a * R_d + A_j_b * g_sto_0_mol * Gamma_star)
        )
    )
    d_with_comp_point = (
        g_sto_0_mol
        * (g_bl_mol_co2**2)
        * Gamma
        * (-c_a * J + A_j_a * c_a * R_d + J * Gamma_star + A_j_b * R_d * Gamma_star)
    )

    # Define Cubic Equation Coefficients
    if negative_A_n == False:
        a = A_j_a * (beta - g_bl_mol_co2)
        b = (
            (A_j_a * g_sto_0_mol * alpha)
            - (beta * xi_j)
            + (A_j_a * g_bl_mol_co2 * alpha)
            + (g_bl_mol_co2 * zeta_j)
            + b_with_comp_point
        )
        c = (
            -(g_sto_0_mol * alpha * xi_j)
            + (beta * eta_j)
            - (g_bl_mol_co2 * alpha * zeta_j)
            + c_with_comp_point
        )
        d = g_sto_0_mol * alpha * eta_j + d_with_comp_point
    elif negative_A_n == True:
        a = 0
        b = A_j_a * (g_sto_0_mol + g_bl_mol_co2)
        c = (g_sto_0_mol * xi_j) + (g_bl_mol_co2 * zeta_j)
        d = g_sto_0_mol * eta_j

    # Solve Cubic Equation
    # roots = np.roots([a, b, c, d])[np.isreal(np.roots([a, b, c, d]))]
    # roots = [np.real(x) for x in roots]
    roots = np.roots([a, b, c, d])
    roots = [np.real(x) for x in roots if (abs(x.imag)<1e-10)]
    return roots


class Ewert_Input_Factors(NamedTuple):
    """Ewert inputs factor output shape.

    These are the outputs to the pre-iteration calculations.
    """

    Tleaf_K: float
    Gamma_star: float
    K_C: float
    K_O: float
    R_d: float
    J_max: float
    V_cmax: float
    J: float
    e_sat_i: float
    Gamma: float


def calc_input_factors(
    Tleaf_C: float,  # < Leaf temperature [degrees C]
    Q: float,  # < PPFD [umol/m^2/s] e.g. PARsun data converted to umol
    V_cmax_25: float,  # < Maximum catalytic rate at 25 degrees [umol/m^2/s]
    J_max_25: float,  # < Maximum rate of electron transport at 25 degrees [umol/m^2/s]
    R_d_coeff: float,  # < Dark respiration coefficient
) -> Ewert_Input_Factors:
    """Calculate input parameters to ewert pass.

    These are the values calculated before running the CO2 iteration.

    Parameters
    ----------
    Tleaf_C:float
        Leaf temperature [degrees C]
    Q:float
        PPFD [umol/m^2/s] e.g. PARsun or PARshade data converted to umol [umol /m^2/s (Photons?)]
    V_cmax_25:float
        Maximum catalytic rate at 25 degrees [umol/m^2/s]
    J_max_25:float
        Maximum rate of electron transport at 25 degrees [umol/m^2/s]
    R_d_coeff: float
        Dark respiration coefficient [Fraction] (Clark et al 2011)


    Returns
    -------
        Ewert_Input_Factors named tuple:

    Gamma_star: float
        CO2 comp. point without day resp.  [micro mol/mol]
    K_C: float
        Michaelis constant CO2 [micro mol/mol]
    K_O: float
        Michaelis constant O2 [mmol/mol]
    R_d: float
        day respiration rate [micro mol/(m^2*s)]
    J_max: float
        Max rate of electron transport     [micro mol/(m^2*s)]
    V_cmax: float
        Max catalytic rate of Rubisco      [micro mol/(m^2*s)]
    J: float
        Rate of electron transport         [micro mol/(m^2*s)]
    e_sat_i: float
        internal saturation vapour pressure[Pa]
    Gamma: float
        CO2 compensation point             [micro mol/mol]

    """
    #  get Leaf temperature in Kelvin
    Tleaf_K = deg_to_kel(Tleaf_C)

    #  Calculation of the model variables which are only
    #  dependent on environmental conditions:
    #  Using Photosynthesis helper functions

    # [micro mol/mol]
    Gamma_star = temp_dep(Gamma_star_25, deg_to_kel(25), E_Gamma_star, Tleaf_K, R)

    # [micro mol/mol]
    K_C = temp_dep(K_C_25, deg_to_kel(25), E_K_C, Tleaf_K, R)

    # [micro mol/mol]
    K_O = temp_dep(K_O_25, deg_to_kel(25), E_K_O, Tleaf_K, R)

    # [umol/m^2/s]
    J_max = temp_dep_inhibit(
        J_max_25, deg_to_kel(25), H_a_jmax, H_d_jmax, S_V_jmax, Tleaf_K, R
    )

    # [umol/m^2/s]
    V_cmax = temp_dep_inhibit(
        V_cmax_25, deg_to_kel(25), H_a_vcmax, H_d_vcmax, S_V_vcmax, Tleaf_K, R
    )

    # [micro mol/(m^2*s)]
    # R_d = temp_dep(R_d_20, deg_to_kel(20), E_R_d, Tleaf_K, R)
    # Dark respiration rate (Clark et al 2011)
    R_d = V_cmax * R_d_coeff

    #  Electron transport rate
    # Equation 4 in Ewert Paper
    # [mol electron?]
    J = (
        J_max
        + alpha * Q
        - sqrt((J_max + alpha * Q) ** 2 - 4 * alpha * Q * Teta * J_max)
    ) / (2 * Teta)  # noqa: 501

    # [Pa]
    e_sat_i = 1000 * saturated_vapour_pressure(Tleaf_C)

    # TODO: Handle V_cmax being 0
    # [mmol/mol]
    Gamma = (Gamma_star + (K_C * R_d * (1 + (O_i / K_O)) / V_cmax)) / (
        1 - (R_d / V_cmax)
    )

    return Ewert_Input_Factors(
        Tleaf_K=Tleaf_K,
        Gamma_star=Gamma_star,
        K_C=K_C,
        K_O=K_O,
        R_d=R_d,
        J_max=J_max,
        V_cmax=V_cmax,
        J=J,
        e_sat_i=e_sat_i,
        Gamma=Gamma,
    )


def calc_fO3_h(
    O3up: float,
    gamma_1: float = 0.060,
    gamma_2: float = 0.0045,
) -> float:
    """Calculate f03 as per Ewert & Porter(2000) equation 10.

    # fO3_h = min(max(1 + gamma_1 - gamma_2 * O3up, 0), 1)

    Parameters
    ----------
    O3up_out: float
        Ozone uptake [nmol O3 PLA m-2 s-1 O3]
    gamma_1: float
        short term damage coefficient [dimensionless]
    gamma_2: float
        short term damage coefficient [nmol m-2 O3]

    Returns
    -------
    fO3_h: float
        Hourly ozone impace factor [dimensionless][0-1] see ewert

    References
    ----------
    - Ewert, F., Porter, J.R., 2000

    """
    lower_bound = gamma_1 / gamma_2
    upper_bound = (1 + gamma_1) / gamma_2

    if O3up <= lower_bound:
        fO3_h = 1
    elif O3up >= upper_bound:
        fO3_h = 0
    else:
        fO3_h = 1 + gamma_1 - gamma_2 * O3up

    return fO3_h


def calc_f_LA(t_lem: float, t_lma: float, td_dd: float) -> float:
    """Calculate f_LA as per Ewert(2000) equations 14-15.

    f_LA = min(max(1 - (td_dd-t_lem)/(t_lma), 0), 1)
    eq 14/15 between t_lem and t_lem+t_lma

    https://www.desmos.com/calculator/h66l70w4wq

    Parameters
    ----------
    t_lem: float
        time from seed to end of emerging leaf/start of mature leaf [Thermal Time]
    t_lma: float
        Full time of mature leaf (t_lep + t_lse) [Thermal Time]
    td_dd: float
        difference between current td and td at season_Astart [Thermal Time]

    Returns
    -------
    f_LA: float
        Leaf repair capacity age factor. Between 0-1 over t_lma [dimensionless]


    References
    ----------
    - Ewert, F., Porter, J.R., 2000

    """
    lower_bound = t_lem
    upper_bound = t_lem + t_lma
    if td_dd <= lower_bound:
        f_LA = 1
    elif td_dd >= upper_bound:
        f_LA = 0
    else:
        f_LA = 1 - (td_dd - t_lem) / (t_lma)

    return f_LA


def calc_fo3_l(
    O3up_acc: float,
    gamma_3: float,
    cL3: float,
) -> float:
    """Calculate the effect of cumulative ozone on lifespan.

    NOTE: fO3_l equation has been modified from original Ewert method to match the Feng et al (2022) method

    Parameters
    ----------
    O3up_acc : float
        Accumulated Ozone uptake [nmol PLA m-2 s-1 O3]
    gamma_3 : float
        long term damage coefficient [umol m-2 O3]
    cL3: float
        critical cumulative stomatal ozone flux [umol m-2 O3]. I.e the point at which
        ozone accumulation affects the onset of senescence


    Returns
    -------
    float
        Ozone exposure factor [umol/m^2 O3]


    References
    ----------
    - Ewert, F., Porter, J.R., 2000
    - Feng et al (2022)

    """
    # eq 17 # note conversion to umol
    return min(1, max(0, 1 - (gamma_3 * ((O3up_acc / 1000) - cL3))))


class Damage_Factors(NamedTuple):
    """Ozone damage factors."""

    fO3_h: float
    fO3_d: float
    fO3_l: float
    f_LA: float
    f_LS: float
    rO3: float
    t_lep_limited: float
    t_lse_limited: float


def calc_ozone_damage_factors(
    gamma_1: float,
    gamma_2: float,
    gamma_3: float,
    cL3: float,
    O3up: float,
    O3up_acc: float,
    fO3_d_prev: float,
    td_dd: float,
    t_lem: float,
    t_lma: float,
    is_daylight: bool,
    hr: int,
    opt_full_night_recovery: bool = False,
) -> Damage_Factors:
    """Ewert calculate fO3_h, fO3_d and fO3_l.

    NOTE: leaf life span values are estimated without Ozone impact from the phyllochron


    Parameters
    -----------
    gamma_1: float
        short term damage coefficient [dimensionless]
    gamma_2: float
        short term damage coefficient [nmol m-2 O3]
    gamma_3: float
        long term damage coefficient [umol m-2 O3]
    cL3: float
        critical cumulative stomatal ozone flux [umol m-2 O3]. I.e the point at which
        ozone accumulation affects the onset of senescence
    O3up: float
        Ozone uptake [nmol PLA m-2 s-1 O3]
    O3up_acc: float
        Accumulated Ozone uptake [nmol PLA m-2 s-1 O3]
    fO3_d_prev: float
        Cumulative ozone effect from previous hour [dimensionless][0-1]
    td_dd: float
        difference between current td and td at season_Astart [Thermal Time]
    t_lem: float
        time from seed to end of emerging leaf/start of mature leaf [Thermal Time]
    t_lma: float
        Full time of mature leaf (t_lep + t_lse) [Thermal Time]
    is_daylight: bool
        True when PAR is greater than 50 W m^2 [bool]
    hr: int
        Hour of day [0-23]
    opt_full_night_recovery: bool
        If True leaf recovers every none daylight hour, else only at hr==0

    Returns
    -------
    fO3_h: float
        Hourly ozone impace factor [dimensionless][0-1]
    fO3_d: float
        Hourly accumulated ozone impace factor [dimensionless][0-1]
    fO3_l: float
        Ozone exposure factor [umol/m^2 O3]
    f_LA: float
        Leaf repair capacity age factor. Between 0-1 over t_lma [dimensionless]
    rO3: float
        Leaf capacity to recover from Ozone damage [dimensionless][0-1]

    References
    ----------
    - Ewert, F., Porter, J.R., 2000

    """
    fO3_h = calc_fO3_h(O3up, gamma_1, gamma_2)
    f_LA = calc_f_LA(t_lem, t_lma, td_dd)
    rO3 = fO3_d_prev + (1 - fO3_d_prev) * f_LA  # eq 13

    fO3_d_base = fO3_d_prev * fO3_h  # eq 11 fO3_d degrades during the day
    fO3_d_night = rO3 * fO3_h  # eq 12  fO3_d recovers during the night
    if opt_full_night_recovery:
        if is_daylight:
            fO3_d = fO3_d_base
        else:
            fO3_d = fO3_d_night
    else:
        if hr == 0:
            fO3_d = fO3_d_night
        else:
            fO3_d = fO3_d_base

    fO3_l = calc_fo3_l(O3up_acc, gamma_3, cL3)

    return Damage_Factors(
        fO3_h=fO3_h,
        fO3_d=fO3_d,
        fO3_l=fO3_l,
        f_LA=f_LA,
        rO3=rO3,
        # NOTE: These are calculated elsewhere
        f_LS=None,  # TODO: Remove this # type: ignore
        t_lep_limited=None, # type: ignore
        t_lse_limited=None, # type: ignore
    )


def calc_all_ozone_damage_factors(
    gamma_1: float,
    gamma_2: float,
    gamma_3: float,
    gamma_4_senes: float,
    gamma_5_harvest: float,
    cL3: float,
    O3up: float,
    O3up_acc: float,
    fO3_d_prev: float,
    td_dd: float,
    t_lem: float,
    t_lep: float,
    t_lse: float,
    is_daylight: bool,
    hr: int,
    opt_full_night_recovery: bool = False,
    use_O3_damage: bool = True,
) -> Damage_Factors:
    """Ewert calculate fO3_h, fO3_d and fO3_l.

    Parameters
    -----------
    gamma_1: float
        short term damage coefficient [dimensionless]
    gamma_2: float
        short term damage coefficient [nmol m-2 O3]
    gamma_3: float
        long term damage coefficient [umol m-2 O3]
    gamma_4_senes: float
        Fraction impact of fO3_l on start of senesense [fraction]
    gamma_5_harvest: float
        Fraction impact of fO3_l on harvest date [fraction]
    cL3: float
        critical cumulative stomatal ozone flux [umol m-2 O3]. I.e the point at which
        ozone accumulation affects the onset of senescence
    O3up: float
        Ozone uptake [nmol PLA m-2 s-1 O3]
    O3up_acc: float
        Accumulated Ozone uptake [nmol PLA m-2 s-1 O3]
    fO3_d_prev: float
        Cumulative ozone effect from previous hour [dimensionless][0-1]
    td_dd: float
        difference between current td and td at season_Astart [Thermal Time]
    t_lem: float
        time from seed to end of emerging leaf/start of mature leaf [Thermal Time]
    t_lma: float
        Full time of mature leaf (t_lep + t_lse) [Thermal Time]
    t_lep: float
        Time during which leaf is expanding [Thermal Time]
    t_lse: float
        Time that the leaf is senescing [Thermal Time]
    is_daylight: bool
        True when PAR is greater than 50 W m^2 [bool]
    hr: int
        Hour of day [0-23]
    opt_full_night_recovery: bool
        If True leaf recovers every none daylight hour, else only at hr==0
    use_O3_damage: bool
        If True then leaf life span is reduced by ozone


    Returns
    -------
    fO3_h: float
        Hourly ozone impace factor [dimensionless][0-1]
    fO3_d: float
        Hourly accumulated ozone impace factor [dimensionless][0-1]
    fO3_l: float
        Ozone exposure factor [umol/m^2 O3]
    f_LA: float
        Leaf repair capacity age factor. Between 0-1 over t_lma [dimensionless]
    rO3: float
        Leaf capacity to recover from Ozone damage [dimensionless][0-1]

    References
    ----------
    - Ewert, F., Porter, J.R., 2000

    """
    fO3_l = calc_fo3_l(O3up_acc, gamma_3, cL3) if use_O3_damage else 1
    fO3_h = calc_fO3_h(O3up, gamma_1, gamma_2) if use_O3_damage else 1

    lifespan_with_ozone = calc_ozone_impact_on_lifespan(
        t_lep=t_lep,
        t_lse=t_lse,
        t_lem=t_lem,
        fO3_l=fO3_l,
        gamma_4_senes=gamma_4_senes,
        gamma_5_harvest=gamma_5_harvest,
    )

    f_LA = calc_f_LA(t_lem, lifespan_with_ozone.t_lma_O3, td_dd)
    rO3 = fO3_d_prev + (1 - fO3_d_prev) * f_LA if use_O3_damage else 1  # eq 13

    fO3_d_base = fO3_d_prev * fO3_h  # eq 11 fO3_d degrades during the day
    fO3_d_night = rO3 * fO3_h  # eq 12  fO3_d recovers during the night

    if opt_full_night_recovery:
        if is_daylight:
            fO3_d = fO3_d_base
        else:
            fO3_d = fO3_d_night
    else:
        if hr == 0:
            fO3_d = fO3_d_night
        else:
            fO3_d = fO3_d_base

    f_LS = calc_senescence_factor(
        td_dd=td_dd,
        t_l_O3=lifespan_with_ozone.t_l_O3,
        t_lem=t_lem,
        t_lep_limited=lifespan_with_ozone.t_lep_limited,
        t_lse_limited=lifespan_with_ozone.t_lse_limited,
    )

    return Damage_Factors(
        fO3_h=fO3_h,
        fO3_d=fO3_d,
        fO3_l=fO3_l,
        f_LA=f_LA,
        f_LS=f_LS,
        rO3=rO3,
        t_lep_limited=lifespan_with_ozone.t_lep_limited,
        t_lse_limited=lifespan_with_ozone.t_lse_limited,
    )


class Leaf_Life_Span_Values(NamedTuple):
    """Values associated with leaf lifespan."""

    t_lma_O3: float
    t_lse_limited: float
    t_lep_limited: float
    t_l_O3: float


def calc_ozone_impact_on_lifespan(
    t_lep: float,
    t_lse: float,
    t_lem: float,
    fO3_l: float,
    gamma_4_senes: float = 1,
    gamma_5_harvest: float = 1,
) -> Leaf_Life_Span_Values:
    """Ewert hour pass.

    Note: leaf life span values are estimated without Ozone impact from
    cultivar parameters and thermal time.

    eq 16 # reduces estimated leaf age by accumulated fst

    NOTE: We have modified this equation from original Ewert Eq 16
    so that senescence can start earlier but still finish at the same date
    as none ozone impacted leaf.

    To modify this effect we adjust the parameter gamma_4 and gamma_5 parameters.
    A gamma_4 value of 0 means that fO3_l has no effect on the start of senesence.
    A gamma_5 value of 0 means that fO3_l has no effect on the end of senesence.


    A interactive graph demostrating this can be found here:
    https://www.desmos.com/calculator/dydaugsddd


    Original Ewert method located here:
    https://www.desmos.com/calculator/odbnt9a6bv


    Parameters
    -----------
    t_lep: float
        Time during which leaf is expanding [Thermal Time]
    t_lse: float
        Time that the leaf is senescing [Thermal Time]
    t_lem: float
        time from seed to end of emerging leaf/start of mature leaf [Thermal Time]
    fO3_l: float
        Ozone exposure factor [umol/m^2 CHECK GAS]
    gamma_4_senes: float
        Fraction impact of fO3_l on start of senesense [fraction]
    gamma_5_harvest: float
        Fraction impact of fO3_l on harvest date [fraction]

    Returns
    -------
    t_lma_O3: float
        t_lma asjusted for ozone impact [Thermal Time]
    t_lse_limited: float
        t_lse asjusted for ozone impact [Thermal Time]
    t_lep_limited: float
        t_lep asjusted for ozone impact [Thermal Time]

    References
    ----------
    - Ewert, F., Porter, J.R., 2000

    """
    # t_lma_O3 = (t_lep + t_lse) * fO3_l  #
    # t_lep_limited = t_lep * fO3_l
    # t_lse_limited =t_lse * (fO3_l + fO3_l_f*(1 - fO3_l)) + fO3_l_f*(t_lep - t_lep_limited) # See note on modified t_lse_limited eq
    g4 = gamma_4_senes
    g5 = gamma_5_harvest

    t_lep_limited = t_lep * (1 - g4 * (1 - fO3_l))
    t_lse_limited = t_lse * (1 - g5 * (1 - fO3_l)) + (t_lep - t_lep_limited)
    t_lma_O3 = t_lep_limited + t_lse_limited

    # non senescing period assumed to be remaining time
    t_l_O3 = t_lem + t_lma_O3
    return Leaf_Life_Span_Values(
        t_lma_O3=t_lma_O3,
        t_lse_limited=t_lse_limited,
        t_lep_limited=t_lep_limited,
        t_l_O3=t_l_O3,
    )


def calc_senescence_factor(
    td_dd: float,
    t_l_O3: float,
    t_lem: float,
    t_lep_limited: float,
    t_lse_limited: float,
) -> float:
    """Calculate senescence factor (f_LS).

    This is the factor that accounts for effect of leaf senescence on A_c
    Ewert & Porter(2000) Eq 18

    Note: we have replaced the bottom line of eq 18 with t_lse_o3. This was to ensure that end
    of senescence matched the end of t_l.

    https://www.desmos.com/calculator/h66l70w4wq

    Parameters
    ----------
    td_dd: float
        difference between current td and td at leaf emergence [Thermal Time]
    t_l_O3: float
        Ozone exposed estimated life span of leaf [Thermal Time]
    t_lem: float
        time from seed to end of emerging leaf/start of mature leaf [Thermal Time]
    t_lep_limited: float
        Time during which leaf is expanding [Thermal Time]
    t_lse_limited: float
        t_lse asjusted for ozone impact [Thermal Time]

    Returns
    -------
    f_LS: float
        factor effect of leaf senescence on A_c [Dimensionless][0-1]

    References
    ----------
    - Ewert, F., Porter, J.R., 2000

    """
    s_signal = t_lem + t_lep_limited  # new senescense signal
    upper_bound = t_l_O3
    # assert isclose(t_l_O3, t_lem + t_lep_limited + t_lse_limited, abs_tol=1e-3)

    if td_dd <= s_signal:
        f_LS = 1
    elif td_dd >= upper_bound:
        f_LS = 0
    else:
        # NOTE: Below function uses O3 vars instead of original vars(See docs)
        f_LS = 1 - ((td_dd - t_lem - t_lep_limited) / t_lse_limited)
        # f_LS = 1 - ((td_dd - t_lem - t_lep_limited) / (t_lma / fO3_l - t_lep_limited))

    return f_LS


class CO2_assimilation_rate_factors(NamedTuple):
    """Variables associated with CO2 assimilation rate."""

    A_c: float
    A_j: float
    A_p: float
    A_n: float
    A_n_limit_factor: str


def calc_CO2_assimilation_rate(
    c_i_in: float,
    V_cmax: float,
    Gamma_star: float,
    K_C: float,
    K_O: float,
    fO3_d: float,
    f_LS: float,
    J: float,
    R_d: float,
) -> CO2_assimilation_rate_factors:
    """Calculate the assimilation rates.

    Ewert & Porter(2000) Eq 8

    Parameters
    ----------
    c_i_in: float
        See ewert
    V_cmax: float
        Max catalytic rate of Rubisco      [micro mol/(m^2*s)]
    Gamma_star: float
        See ewert
    K_C: float
        Michaelis constant CO2             [micro mol/mol]
    K_O: float
        Michaelis constant O2              [mmol/mol]
    fO3_d: float
        Hourly accumulated ozone impact factor [dimensionless][0-1]
        Cumulative ozone effect [dimensionless]
    f_LS: float
        factor effect of leaf senescence on A_c [Dimensionless][0-1]
    J: float
        Rate of electron transport         [micro mol/(m^2*s)]
    R_d: float
        day respiration rate [micro mol/(m^2*s)]

    Returns
    -------
    A_c: float
        eq 8 - Rubisco activity limited rate of photosynthesis
    A_j: float
        eq 3 - RuBP regeneration (electron transport) limited assimilation rate (A_q in paper)
    A_p: float
        eq ? - Triose phosphate utilisation limited assimilation rate
    A_n: float
        eq 1 - CO2 Assimilation Rate [umol m-2 s-1 CO2]
    A_n_limit_factor: str
        The factor that has limited CO2 assimilation (A_c/A_j/A_p)

    References
    ----------
    - Ewert, F., Porter, J.R., 2000 - Eq8

    """
    A_c = (
        V_cmax
        * ((c_i_in - Gamma_star) / (c_i_in + (K_C * (1 + (O_i / K_O)))))
        * fO3_d
        * f_LS
    )
    A_j = J * ((c_i_in - Gamma_star) / ((A_j_a * c_i_in) + (A_j_b * Gamma_star)))

    A_p = 0.5 * V_cmax

    A_n = min(A_c, A_j, A_p) - R_d

    A_n_limit_factor = sorted(
        zip(["A_c", "A_j", "A_p"], [A_c, A_j, A_p]), key=lambda tup: tup[1]
    )[0][0]

    return CO2_assimilation_rate_factors(
        A_c=A_c,
        A_j=A_j,
        A_p=A_p,
        A_n=A_n,
        A_n_limit_factor=A_n_limit_factor,
    )


def does_A_n_solve_requirements(
    A_n: float,
    g_sto_0: float,
    g_sto_prev: float,
    e_a: float,
    g_bl: float,
    e_sat_i: float,
    D_0: float,
    f_VPD: float,
    m: float,
    Gamma: float,
    c_a: float,
    f_SW: float,
    f_VPD_method: FVPDMethods,
    negative_A_n: bool,
):
    """Determines if a suitable solution to the cubic equation satisfies the requirements to be the A_n value

    g_sto > 0
    c_i > 0

    Based on the methods by (Masutomi, Y., 2023)

    Parameters
    ----------
    A_n: float
        CO2 Assimilation Rate [umol m-2 s-1 CO2]
    g_sto_0: float
        Closed stomata conductance [umol/m^2/s CO2]
    g_sto_prev: float
        g_sto from previous hour [umol/m^2/s CO2]
    e_a: float
        Ambient vapour pressure [Pa]
    g_bl: float
        Boundary layer conductance to H2O vapour [umol m-2 PLA s-1 H2O]
    e_sat_i: float
        Internal saturation vapour pressure [Pa]
    D_0: float
        The VPD at which g_sto is reduced by a factor of 2 [kPa]
    f_VPD: float
        Humidity Function (Leuning 1995)
    m: float
        Species-specific sensitivity to An [dimensionless]
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    c_a: float
        Atmospheric CO2 concentration [ppm]
    f_SW: float
        Soil water stress factor [dimensionless]
    negative_A_n: bool
        Boolean variable to select a negative A_n
    f_VPD_method: FVPDMethods
        Method to calculate f_VPD

    Returns
    -------
    boolean
        True if the root is a suitable A_n value, False if not


    References
    ----------
    - Masutomi, Y., 2023. The appropriate analytical solution for coupled leaf photosynthesis
    and stomatal conductance models for C3 plants. Ecological Modelling, 481, p.110306.

    """

    # TODO: We are calculating these multiple times. There may be a more efficient way to do this
    f_VPD = (
        calc_humidity_defecit_fVPD(
            # TODO: We should use prev hour gsto here instead of gsto_0
            g_sto_in=g_sto_prev,
            e_a=e_a,
            g_bl=g_bl,
            e_sat_i=e_sat_i,
            D_0=D_0,
            f_VPD_method=f_VPD_method,
        )
        if f_VPD_method in [FVPDMethods.LEUNING, FVPDMethods.DANIELSSON]
        else f_VPD
    )

    g_sto = calc_stomatal_conductance(
        g_sto_0=g_sto_0,
        m=m,
        Gamma=Gamma,
        g_bl=g_bl,
        c_a=c_a,
        A_n=A_n,
        f_SW=f_SW,
        f_VPD=f_VPD,
    )
    c_i = calc_CO2_supply(
        A_n=A_n,
        c_a=c_a,
        g_sto=g_sto,
        g_bl=g_bl,
    )

    # This first if statement is to avoid overly large roots
    # Not a proper solution but a workaround
    if abs(A_n) > 100:
        return False
    if not negative_A_n:
        if A_n < 0:
            return False

    if (g_sto >= 0) and (c_i >= 0):
        return True
    else:
        return False


def calc_CO2_assimilation_rate_cubic(
    V_cmax: float,
    Gamma_star: float,
    K_C: float,
    K_O: float,
    fO3_d: float,
    f_LS: float,
    J: float,
    R_d: float,
    G_1c: float,
    c_a: float,
    RH: float,
    Gamma: float,
    P: float,
    g_sto_0: float,
    g_sto_prev: float,
    e_a: float,
    g_bl: float,
    e_sat_i: float,
    D_0: float,
    f_VPD: float,
    m: float,
    f_SW: float,
    f_VPD_method: FVPDMethods,
) -> CO2_assimilation_rate_factors:
    """Calculate the assimilation rates.

    The original method is based on :
    Ewert & Porter(2000) Eq 8

    In Dec 2024 we revised the method of calculating the assimilation rates to use
    a cubic equation as in (Masutomi, Y., 2023)

    Parameters
    ----------
    V_cmax: float
        Max catalytic rate of Rubisco [micro mol/(m^2*s)]
    Gamma_star: float
        See ewert
    K_C: float
        Michaelis constant CO2 [micro mol/mol]
    K_O: float
        Michaelis constant O2 [mmol/mol]
    fO3_d: float
        Hourly accumulated ozone impact factor [dimensionless][0-1]
    f_LS: float
        Factor effect of leaf senescence on A_c [Dimensionless][0-1]
    J: float
        Rate of electron transport [micro mol/(m^2*s)]
    R_d: float
        Day respiration rate [micro mol/(m^2*s)]
    G_1c: float
        Model parameter [mol/m^2/s]
    c_a: float
        Atmospheric CO2 concentration [ppm]
    RH: float
        Relative Humidity []
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    P: float
        Surface Air Pressure [kPa]
    g_sto_0: float
        Closed stomata conductance [umol/m^2/s CO2]
    g_sto_prev: float
        g_sto from previous hour [umol/m^2/s CO2]
    e_a: float
        Ambient vapour pressure [Pa]
    g_bl: float
        Boundary layer conductance to H2O vapour [umol m-2 PLA s-1 H2O]
    e_sat_i: float
        Internal saturation vapour pressure [Pa]
    D_0: float
        The VPD at which g_sto is reduced by a factor of 2 [kPa]
    f_VPD: float
        Humidity Function (Leuning 1995)
    m: float
        Species-specific sensitivity to An [dimensionless]
    f_SW: float
        Soil water stress factor [dimensionless]
    f_VPD_method: FVPDMethods
        Method to calculate f_VPD

    Returns
    -------
    A_c: float
        Eq 8 - Rubisco activity limited rate of photosynthesis
    A_j: float
        Eq 3 - RuBP regeneration (electron transport) limited assimilation rate (A_q in paper)
    A_p: float
        Eq ? - Triose phosphate utilisation limited assimilation rate
    A_n: float
        Eq 1 - CO2 Assimilation Rate [umol m-2 s-1 CO2]
    A_n_limit_factor: str
        The factor that has limited CO2 assimilation (A_c/A_j/A_p)

    References
    ----------
    - Ewert, F., Porter, J.R., 2000 - Eq8
    - Masutomi, Y., 2023. The appropriate analytical solution for coupled leaf photosynthesis
    and stomatal conductance models for C3 plants. Ecological Modelling, 481, p.110306.

    """

    O2 = 20900  # Partial pressure of atmospheric Oxygen

    negative_A_n = False
    A_p = calc_A_n_product_limited(V_cmax, R_d)
    A_c_roots = calc_A_n_rubisco_limited(
        O2,
        P,
        K_O,
        V_cmax,
        R_d,
        g_sto_0,
        G_1c,
        g_bl,
        c_a,
        RH,
        Gamma,
        Gamma_star,
        K_C,
        fO3_d,
        f_LS,
        f_SW,
        negative_A_n,
    )
    A_c_root = next(
        (
            root
            for root in A_c_roots
            if does_A_n_solve_requirements(
                root,
                g_sto_0,
                g_sto_prev,
                e_a,
                g_bl,
                e_sat_i,
                D_0,
                f_VPD,
                m,
                Gamma,
                c_a,
                f_SW,
                f_VPD_method,
                negative_A_n,
            )
        ),
        np.nan,
    )
    # if np.isnan(A_c_root):
    #     A_c_root = next(
    #     (
    #         root
    #         for root in A_c_roots
    #         if does_A_n_solve_requirements(
    #             root,
    #             g_sto_0,
    #             g_sto_prev,
    #             e_a,
    #             g_bl,
    #             e_sat_i,
    #             D_0,
    #             f_VPD,
    #             m,
    #             Gamma,
    #             c_a,
    #             f_SW,
    #             f_VPD_method,
    #             True,
    #         )
    #     ),
    #     np.nan,
    #     )
    A_c = ((A_c_root + R_d) * fO3_d * f_LS) - R_d
    A_j_roots = calc_A_n_rubp_limited(
        P,
        R_d,
        g_sto_0,
        G_1c,
        g_bl,
        c_a,
        RH,
        Gamma,
        Gamma_star,
        J,
        negative_A_n,
    )
    A_j = next(
        (
            root
            for root in A_j_roots
            if does_A_n_solve_requirements(
                root,
                g_sto_0,
                g_sto_prev,
                e_a,
                g_bl,
                e_sat_i,
                D_0,
                f_VPD,
                m,
                Gamma,
                c_a,
                f_SW,
                f_VPD_method,
                negative_A_n,
            )
        ),
        np.nan,
    )
    if any(a is None for a in [A_c, A_j, A_p]):
        raise ValueError(
            f"One of the assimilation rates is None A_c={A_c}, A_j={A_j}, A_p={A_p}"
        )
    A_n = np.nanmin([A_p, A_c, A_j])

    A_n_limit_factor = sorted(
        zip(["A_c", "A_j", "A_p"], [A_c, A_j, A_p]),
        key=lambda tup: np.isnan(tup[1]) and 10000000 or tup[1],
    )[0][0]

    return CO2_assimilation_rate_factors(
        A_c=A_c,
        A_j=A_j,
        A_p=A_p,
        A_n=A_n,
        A_n_limit_factor=A_n_limit_factor,
    )


def calc_humidity_defecit_fVPD(
    g_sto_in: float,
    e_a: float,
    g_bl: float,
    e_sat_i: float,
    D_0: float,
    f_VPD_method,
) -> float:
    """Calculate the humidity defecit from the gsto.

    Only calculated if VPD not calculated externally.

    Parameters
    ----------
    g_sto_in: float
        g_sto from previous iteration [umol/m^2/s CO2]
    e_a: float
        Ambient vapour pressure (Pa)
    g_bl: float
        Boundary layer conductance to H2O vapour ([umol m-2 PLA s-1 H2O])
    e_sat_i: float
        internal saturation vapour pressure[Pa]
    D_0: float
        "The VPD at which g_sto is reduced by a factor of 2" [kPa] (Leuning et al. 1998)

    Returns
    -------
    fVPD: float
        Humidity Function (Leuning 1995)

    References
    ----------
    - Nikolov, N.T., Massman, W.J., Schoettle, A.W., 1995
    - LEUNING, R., 1995
    - Wang, Y.-P., Leuning, R., 1998

    """

    # Unit conversions
    # Equation 5 is in mol
    # Convert from umol to mol and Pa to kPa
    g_sto_in_mol = g_sto_in * 1e-6
    g_bl_mol = g_bl * 1e-6
    e_sat_i_kPa = e_sat_i * 1e-3
    e_a_kPa = e_a * 1e-3
    # Surface humidity
    h_s = (g_sto_in_mol * e_sat_i_kPa + g_bl_mol * e_a_kPa) / (
        e_sat_i_kPa * (g_sto_in_mol + g_bl_mol)
    )  # Nikolov et al 1995
    # Convert relative humidity to humidity defecit
    d_s = e_sat_i_kPa - (e_sat_i_kPa * h_s)

    if f_VPD_method == FVPDMethods.LEUNING:
        f_VPD = 1 / (1 + (d_s / D_0))  # Leuning 1995
    elif f_VPD_method == FVPDMethods.DANIELSSON:
        f_VPD = 1 / (1 + (d_s / D_0) ** 8)  # Leuning 1995
    else:
        raise ConfigError(f"Invalid f_VPD_method {f_VPD_method}")

    return f_VPD


def calc_stomatal_conductance(
    g_sto_0: float,
    m: float,
    Gamma: float,
    g_bl: float,
    c_a: float,
    A_n: float,
    f_SW: float,
    f_VPD: float,
) -> float:
    """Calculate the stomatal conductance using photosynthesis method.

    Note: Ewert Equation 5 is in mol

    Parameters
    ----------
    g_sto_in: float
        g_sto from previous iteration [umol/m^2/s CO2]
    g_sto_0: float
        Closed stomata conductance [umol/m^2/s CO2]
    m: float
        Species-specific sensitivity to An [dimensionless]
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    g_bl: float
        Boundary layer conductance to H2O vapour ([umol m-2 PLA s-1 H2O])
    c_a: float
        CO2 concentration [ppm]
    A_n: float
        net CO2 assimilation [umol m-2 s-1 CO2]
    f_VPD: float
        Humidity Function (Leuning 1995)

    Returns
    -------
    g_sto: float
        Stomatal Conductance for CO2 [umol/m^2/s CO2]

    References
    ----------
    - Ewert, F., Porter, J.R., 2000
    - LEUNING, R., 1995

    """
    # Unit conversions
    # Convert from umol to mol
    g_sto_0_mol = g_sto_0 * 1e-6
    g_bl_mol = g_bl * 1e-6
    # NOTE: micro mol/(m^2*s) values are not converted

    # Surface CO2
    # g_bl converted from H2O to CO2
    c_s = c_a - (A_n * (1.37 / g_bl_mol))

    # Stomatal conductance

    # Equation 5 from Ewert paper and Leuning 1995
    g_eq_top = m * A_n * f_SW * f_VPD  # micro mol/(m^2*s)
    g_eq_bottom = c_s - Gamma

    # Restricted to min g_sto_out of g_sto_0 to improve loop convergence
    g_sto_out_mol = g_sto_0_mol + max(0, (g_eq_top / g_eq_bottom))

    # Convert from mol to umol
    g_sto_out_umol = g_sto_out_mol * 1e6
    return g_sto_out_umol


def calc_CO2_supply(
    A_n: float,
    c_a: float,
    g_sto: float,
    g_bl: float,
) -> float:
    """Calculate the CO2 supply at end of iteration.

    The original equation:
    c_i_sup = c_a - ((A_n * (1.6 / g_sto + 1.37 / g_bl)) * 1e6)
    We replaced 1.6/g_sto with 1/g_sto because g_sto is assumed to already
    have been converted from H2O to CO2

    Parameters
    ----------
    A_n: float
        net CO2 assimilation [umol m-2 s-1 CO2]
    c_a: float
        CO2 concentration [ppm]
    g_sto: float
        Stomatal Conductance [umol/m^2/s CO2]
    g_bl: float
        Boundary layer conductance to CO2 vapour ([umol m-2 PLA s-1 CO2])

    Returns
    -------
    c_i_sup: float
        CO2 concentration inside stomata possible through supply [ppm]

    References
    ----------
    -  von Caemmerer & Farquhar (1981)

    """
    # CO2 supply
    # NOTE: umol to mol conversions taken into account for g_sto and g_bl

    # gsto does not need converting here as already in CO2 umol
    # g_bl is converted from H2O to CO2 using Dratio
    c_i_sup = c_a - ((A_n * (1 / g_sto + 1.37 / g_bl)) * 1e6)

    return c_i_sup


def calc_mean_gsto(
    leaf_gsto_list: List[List[float]],  # shape(nL, nP)
    leaf_fLAI_list: List[List[float]],  # shape(nL, nP)
    nL: int,
) -> List[float]:
    """Calculate the mean gsto per layer given percent of each leaf population at each layer.


    Parameters
    ----------
    leaf_gsto_list : List[List[float]]
        [description]

    Returns
    -------
    List[float]
        mean_gsto per layer
    """
    return [
        sum((g * L) for g, L in zip(leaf_gsto_list[iL], leaf_fLAI_list[iL]))
        for iL in range(nL)
    ]

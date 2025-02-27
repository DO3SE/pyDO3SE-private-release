from dataclasses import dataclass
from typing import NamedTuple, List
from pyDO3SE.Config.ConfigEnums import FVPDMethods
from .enums import AdjustNegativeAnMethods

@dataclass
class Output_Shape:
    """Ewert Output Shape.

    Note list types are per leaf population.

    Parameters
    ----------
    g_sv_per_layer: float
        Stomatal conductance to water vapour
        Raw photosynthesis-based stomatal conductance [umol m-2 s-1 CO2]
    g_sv_sunlit: float
        Stomatal conductance to water vapour at top sunlit part of canopy
        Raw photosynthesis-based stomatal conductance [umol m-2 s-1 CO2]
    A_n: float
        eq 1 - net CO2 Assimilation Rate [umol m-2 s-1 CO2]
    A_n_sunlit: float
        eq 1 - net CO2 Assimilation Rate at the top sunlit part of canopy [umol m-2 s-1 CO2]
    A_c: float
        eq 8 - Rubisco activity limited rate of photosynthesis [umol m-2 s-1 CO2]
    A_j: float
        eq 3 - RuBP regeneration (electron transport)
        limited assimilation rate (A_q in paper) [umol m-2 s-1 CO2]
    A_p: float
        eq ? - Triose phosphate utilisation limited assimilation rate [umol m-2 s-1 CO2]
    A_n_limit_factor: str
        The factor that has limited CO2 assimilation [str][A_c/A_j/A_p]
    R_d: float
        day respiration rate [micro mol/(m^2*s) CO2]
    c_i: float
        CO2 concentration inside stomata [ppm]
    c_i_sunlit: float
        CO2 concentration inside stomata for top sunlit part of the canopy [ppm]
    fVPD: float
        Humidity Function (Leuning 1995)
    v_cmax: float
        catalytic rate [umol m-2 s-1]
    J_max: float
        rate of electron transport [umol m-2 s-1]
    loop_iterations: int
        Number of iterations in convergence loop. Useful for debugging

    """

    g_sv_per_layer: List[float]
    g_sv_sunlit: float
    A_n: float
    A_n_sunlit: float
    A_c: float
    A_j: float
    A_p: float
    A_n_limit_factor: List[str | None]
    R_d: float
    c_i: float
    c_i_sunlit: float
    f_VPD: list[float]
    v_cmax: float
    j_max: float
    # TODO: Delete this when negative A_n resolved
    f_VPD_alt: list[float]=None
    loop_iterations: int=None


@dataclass
class ModelOptions:
    """Ewert model options.

    f_VPD_method: FVPDMethods
        If "leuning" or "Danielsson" then f_VPD is calculated internally else should be provided as inputk
    co2_concentration_balance_threshold: float = 0.001
        Threshold (from 0) to consider co2 concentration equation as "balanced"
    co2_concentration_max_iterations: int = 50
        Maximum number of iterations to find co2 concentration solution
    adjust_negative_A_n: AdjustNegativeAnMethods = AdjustNegativeAnMethods.FALSE,
        If True then allow negative A_n values, else return NaN
        If "last_resort" then allow negative A_n values if no other solution is found
        If "clip" then clip negative A_n values to -R_d


    """

    f_VPD_method: FVPDMethods = FVPDMethods.LEUNING
    co2_concentration_balance_threshold: float = 0.001
    co2_concentration_max_iterations: int = 50
    adjust_negative_A_n: AdjustNegativeAnMethods = AdjustNegativeAnMethods.FALSE


@dataclass
class CO2_loop_State:
    """Variables that change over the CO2 convergence loop.

    Parameters
    ----------
    c_i: float
        CO2 concentration inside stomata [umol CO2]
    c_i_diff: float
        difference in c_i between iterations [ppm]
    g_sto: float
        stomatal conductance within the loop [umol/ m^2/s gsto]
    A_n: float
        eq 1 - CO2 Assimilation Rate  [umol m-2 s-1 CO2]
    A_c: float
        eq 8 - Rubisco activity limited rate of photosynthesis [umol m-2 s-1 CO2]
    A_p: float
        eq ? - Triose phosphate utilisation limited assimilation rate [umol m-2 s-1 CO2]
    A_j: float
        eq 3 - RuBP regeneration (electron transport) limited assimilation rate (A_q in paper)
         [umol m-2 s-1 CO2]
    A_n_limit_factor: str
        The factor that has limited CO2 assimilation (A_c/A_j/A_p)
    iterations: float
        Number of iterations in convergence loop
    f_VPD: list[float]
        Humidity Function (Leuning 1995)

    """

    c_i: float
    c_i_diff: float
    g_sto: float
    A_n: float = 0.0
    A_c: float = 0.0
    A_p: float = 0.0
    A_j: float = 0.0
    A_n_limit_factor: str | None = None
    f_VPD: float | None = None
    # TODO: Delete this when negative A_n resolved
    f_VPD_alt: float | None = None
    iterations: int = 0


@dataclass
class CO2_Constant_Loop_Inputs:
    """Inputs to the CO2_convergence loop that remain constant throughout loop.

    Parameters
    ------
    D_0: float
        "The VPD at which g_sto is reduced by a factor of 2" [Pa] (Leuning et al. 1998)
    fmin: float
        Minimum fVPD [fraction]
    c_a: float
        CO2 concentration [ppm]
    e_a: float
        Ambient vapour pressure [Pa]
    g_bl: float
        Boundary layer conductance to H2O vapour [umol m-2 PLA s-1 H2O?]
    g_sto_0: float
        Closed stomata conductance [umol/m^2/s CO2]
    m: float
        Species-specific sensitivity to An [dimensionless]
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    Gamma_star: float
        CO2 comp. point without day resp.  [umol/mol CO2]
    V_cmax: float
        Max catalytic rate of Rubisco [umol/(m^2*s) CO2]
    K_C: float
        Michaelis constant CO2 [umol/mol CHECK GAS]
    K_O: float
        Michaelis constant O2 [mmol/mol CHECK GAS]
    J: float
        Rate of electron transport [umol/(m^2*s) CHECK GAS]
    R_d: float
        day respiration rate [umol/(m^2*s) CHECK GAS]
    e_sat_i: float
        internal saturation vapour pressure[Pa]
    f_SW: float
        Soil water influence on photosynthesis [0-1]
    f_LS: float
        factor effect of leaf senescence on A_c [Dimensionless][0-1]
    fO3_d: float
        Hourly accumulated ozone impace factor [dimensionless][0-1]
    f_VPD: float
        VPD effect on gsto [fraction] (Optional if pre-calculated)
    RH: float
        Relative humidity [%] Only needed for cubic method
    P: float
        Air pressure [kPa]

    """

    c_a: float
    e_a: float
    g_bl: float
    g_sto_0: float
    m: float
    D_0: float
    fmin: float
    Gamma: float
    Gamma_star: float
    V_cmax: float
    K_C: float
    K_O: float
    J: float
    R_d: float
    e_sat_i: float
    f_SW: float
    f_LS: float
    fO3_d: float
    f_VPD: float




@dataclass
class CO2_Concentration_Args:
    """Inputs to the CO2 concentration cubic calculation

    Parameters
    ------
    D_0: float
        "The VPD at which g_sto is reduced by a factor of 2" [Pa] (Leuning et al. 1998)
    fmin: float
        Minimum fVPD [fraction]
    c_a: float
        CO2 concentration [ppm]
    e_a: float
        Ambient vapour pressure [Pa]
    g_bl: float
        Boundary layer conductance to H2O vapour [umol m-2 PLA s-1 H2O?]
    g_sto_0: float
        Closed stomata conductance [umol/m^2/s CO2]
    g_sto_prev: float
        g_sto from previous hour [umol/m^2/s CO2]
    m: float
        Species-specific sensitivity to An [dimensionless]
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    Gamma_star: float
        CO2 comp. point without day resp.  [umol/mol CO2]
    V_cmax: float
        Max catalytic rate of Rubisco [umol/(m^2*s) CO2]
    K_C: float
        Michaelis constant CO2 [umol/mol CHECK GAS]
    K_O: float
        Michaelis constant O2 [mmol/mol CHECK GAS]
    J: float
        Rate of electron transport [umol/(m^2*s) CHECK GAS]
    R_d: float
        day respiration rate [umol/(m^2*s) CHECK GAS]
    e_sat_i: float
        internal saturation vapour pressure[Pa]
    f_SW: float
        Soil water influence on photosynthesis [0-1]
    f_LS: float
        factor effect of leaf senescence on A_c [Dimensionless][0-1]
    fO3_d: float
        Hourly accumulated ozone impace factor [dimensionless][0-1]
    f_VPD: float
        VPD effect on gsto [fraction] (Optional if pre-calculated)
    P: float
        Air pressure [kPa]

    """

    c_a: float
    e_a: float
    g_bl: float
    g_sto_0: float
    g_sto_prev: float
    m: float
    D_0: float
    fmin: float
    Gamma: float
    Gamma_star: float
    V_cmax: float
    K_C: float
    K_O: float
    J: float
    R_d: float
    e_sat_i: float
    f_SW: float
    f_LS: float
    fO3_d: float
    f_VPD: float
    P: float


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


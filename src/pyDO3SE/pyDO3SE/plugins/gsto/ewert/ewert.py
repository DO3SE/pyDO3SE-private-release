r"""The Ewert photosynthesis module.

This module uses the Ewert & Porter (2000) method for calculating stomatal conductance and
leaf CO2 absorbtion.
It is based on a converging loop that finds CO2 absorption and photosynthesis(A_n) based stomatal
conductance(gsto). It converges until the change in CO2 between iterations is below a threshold.

Gsto is in CO2 for all equations.

.. code-block:: python

    ------------------------t_l---------------|
    -----t_lem---|----------t_lma-------------|
                 |-------t_lep--------|-t_lse-|
                 |                    |       |
    = = = = = = = - - - - - - - - - - |       |
                \\                            |
                  \\                   \      |
                    \\                        |
                      \\                \     |
    f_LA = = =          \\                    |
    f_LS - - -            \\             \    |
                            \\                |
                              \\          \   |
                                \\            |
                                  \\       \  |
                                    \\        |
                                      \\    \ |
                                        \\    |
                                          \\ \|
                                            \\|


"""
from dataclasses import dataclass, replace
from typing import List
from pyDO3SE.Config.ConfigEnums import FVPDMethods

from data_helpers.functional_helpers import fp_while

from .ewert_helpers import (
    Ewert_Input_Factors,
    calc_CO2_supply,
    calc_humidity_defecit_fVPD,
    calc_input_factors,
    calc_stomatal_conductance,
    calc_CO2_assimilation_rate,
)


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
    fVPD: float
        Humidity Function (Leuning 1995)
    v_cmax: float
        catalytic rate [umol m-2 s-1]
    J_max: float
        rate of electron transport [umol m-2 s-1]


    """

    g_sv_per_layer: List[float] = None
    g_sv_sunlit: float = None
    A_n: float = None
    A_n_sunlit: float = None
    A_c: float = None
    A_j: float = None
    A_p: float = None
    A_n_limit_factor: float = None
    R_d: float = None
    c_i: float = None
    f_VPD: float = None
    v_cmax: float = None
    j_max: float = None


@dataclass
class ModelOptions:
    """Ewert model options.

    f_VPD_method: FVPDMethods
        If "leuning" or "Danielsson" then f_VPD is calculated internally else should be provided as inputk
    co2_concentration_balance_threshold: float = 0.001
        Threshold (from 0) to consider co2 concentration equation as "balanced"
    co2_concentration_max_iterations: int = 50
        Maximum number of iterations to find co2 concentration solution

    """
    f_VPD_method: FVPDMethods = FVPDMethods.LEUNING
    co2_concentration_balance_threshold: float = 0.001
    co2_concentration_max_iterations: int = 50


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
    f_VPD: float
        Humidity Function (Leuning 1995)

    """

    c_i: float
    c_i_diff: float
    g_sto: float
    A_n: float = 0.0
    A_c: float = 0.0
    A_p: float = 0.0
    A_j: float = 0.0
    A_n_limit_factor: str = None
    f_VPD: float = None
    iterations: int = 0


@dataclass
class CO2_Constant_Loop_Inputs():
    """Inputs to the CO2_convergence loop that remain constant throughout loop.

    Parameters
    ------
    D_0: float
        "The VPD at which g_sto is reduced by a factor of 2" [Pa] (Leuning et al. 1998)
    c_a: float
        CO2 concentration [ppm]
    e_a: float
        Ambient vapour pressure [Pa]
    g_bl: float
        Boundary layer conductance to H2O vapour [umol m-2 PLA s-1 H2O?]
    g_sto_0: float
        Closed stomata conductance [umol/m^2/s CHECK GAS]
    m: float
        Species-specific sensitivity to An [dimensionless]
    Gamma: float
        CO2 compensation point [umol/mol CO2]
    Gamma_star: float
        # CO2 comp. point without day resp.  [umol/mol CO2]
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


    """

    c_a: float
    e_a: float
    g_bl: float
    g_sto_0: float
    m: float
    D_0: float
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
    f_VPD: float = None


def co2_concentration_in_stomata_iteration(
    constant_inputs: CO2_Constant_Loop_Inputs,
    state_current: CO2_loop_State,
    model_options: ModelOptions,
) -> CO2_loop_State:
    """Single iteration to find CO2 concentration (c_i) and stomatal conductance (g_sto).

    by allowing CO2 to converge.



    Parameters
    ----------
    constant_inputs: CO2_Constant_Loop_Inputs
        Inputs that are constant throughout the loop
    state_current: CO2_loop_State
        state that changes through the loop
    model_options: ModelOptions
        Ewert photosynthesis model options

    Returns
    -------
    state: CO2_loop_State
        state after iteration

    """
    # 1. get args from input and state

    # args from previous loop
    c_i_in = state_current.c_i
    g_sto_in = state_current.g_sto

    # 2. Run calculations

    co2_assimilation_rate_values = calc_CO2_assimilation_rate(
        c_i_in=c_i_in,
        V_cmax=constant_inputs.V_cmax,
        Gamma_star=constant_inputs.Gamma_star,
        K_C=constant_inputs.K_C,
        K_O=constant_inputs.K_O,
        fO3_d=constant_inputs.fO3_d,
        f_LS=constant_inputs.f_LS,
        J=constant_inputs.J,
        R_d=constant_inputs.R_d,
    )

    f_VPD = calc_humidity_defecit_fVPD(
        g_sto_in=g_sto_in,
        e_a=constant_inputs.e_a,
        g_bl=constant_inputs.g_bl,
        e_sat_i=constant_inputs.e_sat_i,
        D_0=constant_inputs.D_0,
        f_VPD_method=model_options.f_VPD_method,
    ) if model_options.f_VPD_method in [FVPDMethods.LEUNING, FVPDMethods.DANIELSSON] \
        else constant_inputs.f_VPD

    g_sto = calc_stomatal_conductance(
        g_sto_0=constant_inputs.g_sto_0,
        m=constant_inputs.m,
        Gamma=constant_inputs.Gamma,
        g_bl=constant_inputs.g_bl,
        c_a=constant_inputs.c_a,
        A_n=co2_assimilation_rate_values.A_n,
        f_SW=constant_inputs.f_SW,
        f_VPD=f_VPD,
    )

    co2_supply = calc_CO2_supply(
        A_n=co2_assimilation_rate_values.A_n,
        c_a=constant_inputs.c_a,
        g_sto=g_sto,
        g_bl=constant_inputs.g_bl,
    )

    c_i = c_i_in - (c_i_in - co2_supply) / 2

    new_state = CO2_loop_State(
        c_i=c_i,
        c_i_diff=abs(c_i_in - co2_supply),
        g_sto=g_sto,
        A_n=co2_assimilation_rate_values.A_n,
        A_c=co2_assimilation_rate_values.A_c,
        A_j=co2_assimilation_rate_values.A_j,
        A_p=co2_assimilation_rate_values.A_p,
        A_n_limit_factor=co2_assimilation_rate_values.A_n_limit_factor,
        iterations=state_current.iterations + 1,
        f_VPD=f_VPD,
    )
    return new_state


def co2_concentration_in_stomata_loop(
    constant_inputs: CO2_Constant_Loop_Inputs,
    model_options: ModelOptions,
) -> CO2_loop_State:
    """Run an iterative loop to find cO2 concentration by allowing CO2 to converge.

    This should be ran iteratively until the c_i from both methods converge

    Parameters
    ----------
    constant_inputs: CO2_Constant_Loop_Inputs
        Inputs that do not change in the CO2 loop
    model_options: ModelOptions
        Ewert photosynthesis model options

    Returns
    -------
        A CO2_loop_state named tuple that represents the state at the end of the loop
        See CO2_loop_State class for list of outputs

    """
    tolerance = model_options.co2_concentration_balance_threshold
    max_iterations = model_options.co2_concentration_max_iterations
    c_i_start = 0.0
    initial_c_i_diff = 1

    initial_state = CO2_loop_State(
        c_i=c_i_start,
        c_i_diff=initial_c_i_diff,
        g_sto=constant_inputs.g_sto_0,
    )
    final_state: CO2_loop_State = fp_while(
        lambda state: state.c_i_diff > tolerance,
        lambda state: co2_concentration_in_stomata_iteration(constant_inputs, state, model_options),
        initial_state,
        max_iter=max_iterations,
    )
    return final_state


# @deprecated(reason="use ewert_leaf_pop")
# def ewert(
#     # Config Inputs
#     nP: int,
#     g_sto_0: float,
#     m: float,
#     V_cmax_25: float,
#     J_max_25: float,

#     # State Inputs
#     PARsun: float,
#     PARshade: float,
#     layer_LAI: float,
#     LAIsunfrac: float,
#     D_0: float,
#     g_bv: float,

#     Tleaf_C: float,
#     eact: float,
#     c_a: float,
#     R_d_coeff: float,
#     f_SW: float,
#     f_VPD: float,
#     leaf_pop_distribution: List[float],

#     # TODO: Check ok as lists
#     f_LA: List[float],
#     f_LS: List[float],
#     fO3_d: List[float],

#     # model_options
#     f_VPD_method: FVPDMethods,
# ) -> Output_Shape:
#     """Run the Ewert Photosynthesis model.

#     NOTE: Be careful of units if modifying

#     Parameters
#     ----------
#     nP: int
#         Number of leaf populations
#     g_sto_0: float
#         Closed stomata conductance [umol/m^2/s CO2]
#         The stomatal conductance when A_n -> 0 and incident irradiance, I -> 0
#     m: float
#         Species-specific sensitivity to An [dimensionless]
#     V_cmax_25: float
#         Maximum catalytic rate at 25 degrees [umol/m^2/s]
#     J_max_25: float
#         Maximum rate of electron transport at 25 degrees [umol/m^2/s]
#     PARsun: float
#         PAR received by sunlit leaves [W m-2]
#     PARshade: float
#         PAR received by shaded leaves [W m-2]
#     layer_LAI: float
#         Leaf area index for current layer [m2 m-2]
#     D_0: float
#         "The VPD at which g_sto is reduced by a factor of 2" [kPa] (Leuning et al. 1998)
#     g_bv: float
#         boundary layer conductance for forced convection [umol m-2 s-1 H2O]
#     Tleaf_C: float
#         Leaf Temperature [degrees celsius]
#     eact: float
#         Ambient vapour pressure [kPa]
#     sinB: float
#         sin() of solar elevation angle [degrees]
#     c_a: float
#         CO2 concentration [ppm CO2]
#     R_d_coeff: float
#         Dark respiration coefficient [Fraction] (Clark et al 2011)
#     f_SW: float
#         Soil water influence on photosynthesis [0-1]
#     f_VPD: float
#         VPD effect on gsto [fraction] (Use if pre-calculated)
#     leaf_pop_distribution: float
#         distribution of lai between leaf populations
#     f_LA: float
#         Leaf repair capacity age factor. Between 0-1 over t_lma [dimensionless]
#     f_LS: float
#         factor effect of leaf senescence on A_c [Dimensionless][0-1]
#     fO3_d: float
#         Hourly accumulated ozone impace factor [dimensionless][0-1]
#     f_VPD_method: FVPDMethods
#         If "photosynthesis" then f_VPD is calculated internally else should be provided as input

#     Returns
#     -------
#         a named tuple of shape Output_Shape. See Output_Shape class for list of outputs

#     """
#     # TODO: Do we need to run when total PAR == 0
#     if layer_LAI == 0 or sum(leaf_pop_distribution) == 0:
#         # TODO: Check if any of these should be a min values
#         return Output_Shape(
#             g_sv=0,
#             A_n=0,
#             A_c=0,
#             A_j=0,
#             A_p=0,
#             A_n_limit_factor='NA',
#             R_d=0,
#             c_i=0,
#             f_VPD=f_VPD,
#         )
#     populations_to_run = [i > 0 for i in leaf_pop_distribution]
#     zero_state = CO2_loop_State(
#         c_i=0,
#         c_i_diff=0,
#         g_sto=0,
#         A_n=0,
#         A_c=0,
#         A_p=0,
#         A_j=0,
#         A_n_limit_factor=0,
#         f_VPD=0,
#         iterations=0,
#     )

#     model_options = ModelOptions(
#         f_VPD_method,
#     )

#     e_a = eact * 1e3

#     # ========== Run for PARsun
#     inputs_sun: Ewert_Input_Factors = calc_input_factors(
#         Tleaf_C=Tleaf_C,
#         Q=PARsun * 4.57,
#         V_cmax_25=V_cmax_25,
#         J_max_25=J_max_25,
#         R_d_coeff=R_d_coeff,
#     )

#     loop_inputs_sun: List[CO2_Constant_Loop_Inputs] = [CO2_Constant_Loop_Inputs(
#         c_a=c_a,
#         e_a=e_a,
#         g_bl=g_bv,
#         g_sto_0=g_sto_0,
#         m=m,
#         D_0=D_0,
#         Gamma=inputs_sun.Gamma,
#         Gamma_star=inputs_sun.Gamma_star,
#         V_cmax=inputs_sun.V_cmax,
#         K_C=inputs_sun.K_C,
#         K_O=inputs_sun.K_O,
#         J=inputs_sun.J,
#         R_d=inputs_sun.R_d,
#         e_sat_i=inputs_sun.e_sat_i,
#         f_SW=f_SW,
#         f_LA=f_LA,
#         f_LS=f_LS,
#         fO3_d=fO3_d,
#         f_VPD=f_VPD,
#     ) for iP in range(nP)]

#     state_out_sun: List[CO2_loop_State] = [
#         co2_concentration_in_stomata_loop(loop_inputs_sun[iP], model_options)
#         if should_run else zero_state
#         for iP, should_run in enumerate(populations_to_run)
#     ]

#     # ============= Run for PARshade

#     inputs_shade: Ewert_Input_Factors = calc_input_factors(
#         Tleaf_C=Tleaf_C,
#         Q=PARshade * 4.57,
#         V_cmax_25=V_cmax_25,
#         J_max_25=J_max_25,
#         R_d_coeff=R_d_coeff,
#     )

#     loop_inputs_shade: List[CO2_Constant_Loop_Inputs] = [replace(
#         loop_inputs_sun[iP],
#         Gamma=inputs_shade.Gamma,
#         Gamma_star=inputs_shade.Gamma_star,
#         V_cmax=inputs_shade.V_cmax,
#         K_C=inputs_shade.K_C,
#         K_O=inputs_shade.K_O,
#         J=inputs_shade.J,
#         R_d=inputs_shade.R_d,
#         e_sat_i=inputs_shade.e_sat_i,
#     ) for iP in range(nP)]

#     state_out_shade: List[CO2_loop_State] = [
#         co2_concentration_in_stomata_loop(loop_inputs_shade[iP], model_options)
#         if should_run else zero_state
#         for iP, should_run in enumerate(populations_to_run)
#     ]

#     # Get sun shade parts
#     # TODO: total_lai is a temp fix as layer_lai and leaf pop do not match
#     total_lai = sum(leaf_pop_distribution)
#     sun_fracs = [LAIsunfrac * (pop_lai / total_lai) for pop_lai in leaf_pop_distribution]
#     shade_fracs = [(1 - LAIsunfrac) * (pop_lai / total_lai) for pop_lai in leaf_pop_distribution]
#     try:
#         assert isclose(sum(sun_fracs) + sum(shade_fracs), 1.0)
#     except AssertionError as e:
#         print(sun_fracs)
#         print(shade_fracs)
#         raise e
#     shade_frac = 1 - LAIsunfrac

#     # Get Total sun shade values
#     # TODO: Check g_sv upscaling
#     g_sv_canopy = sum([state_out_sun[iP].g_sto * sun_fracs[iP] +
#                        state_out_shade[iP].g_sto * shade_fracs[iP] for iP in range(nP)])
#     g_sv_per_pop = [state_out_sun[iP].g_sto * LAIsunfrac +
#                     state_out_shade[iP].g_sto * shade_frac for iP in range(nP)]

#     A_n = sum([state_out_sun[iP].A_n * sun_fracs[iP] + state_out_shade[iP].A_n *
#                shade_fracs[iP] for iP in range(nP)])
#     A_c = sum([state_out_sun[iP].A_c * sun_fracs[iP] + state_out_shade[iP].A_c *
#                shade_fracs[iP] for iP in range(nP)])
#     A_j = sum([state_out_sun[iP].A_j * sun_fracs[iP] + state_out_shade[iP].A_j *
#                shade_fracs[iP] for iP in range(nP)])
#     A_p = sum([state_out_sun[iP].A_p * sun_fracs[iP] + state_out_shade[iP].A_p *
#                shade_fracs[iP] for iP in range(nP)])
#     c_i = sum([state_out_sun[iP].c_i * sun_fracs[iP] + state_out_shade[iP].c_i *
#                shade_fracs[iP] for iP in range(nP)])
#     R_d = sum([inputs_sun.R_d * sun_fracs[iP] + inputs_shade.R_d *
#                shade_fracs[iP] for iP in range(nP)])

#     # TODO: How do we deal with multi pop f_VPD
#     f_VPD_out = [state_out_sun[iP].f_VPD * LAIsunfrac +
#                  state_out_shade[iP].f_VPD * shade_frac for iP in range(nP)][0]

#     v_cmax_out = inputs_sun.V_cmax * LAIsunfrac + inputs_shade.V_cmax * shade_frac
#     j_max_out = inputs_sun.J_max * LAIsunfrac + inputs_shade.J_max * shade_frac

#     # Other
#     A_n_limit_factor = [state_out_sun[i].A_n_limit_factor for i in range(nP)]

#     return Output_Shape(
#         g_sv=g_sv_canopy,
#         A_n=A_n,
#         A_c=A_c,
#         A_j=A_j,
#         A_p=A_p,
#         A_n_limit_factor=A_n_limit_factor,
#         R_d=R_d,
#         c_i=c_i,
#         f_VPD=f_VPD_out,
#         v_cmax=v_cmax_out,
#         j_max=j_max_out,
#     )


def ewert_leaf_pop(
    # Config Inputs
    nL: int,
    g_sto_0: float,
    m: float,
    V_cmax_25: List[float],
    J_max_25: List[float],
    R_d_coeff: float,

    # State Inputs per layer
    PARsun: List[float],
    PARshade: List[float],
    LAIsunfrac: List[float],
    layer_lai_frac: List[float],
    layer_lai: List[float],
    Tleaf_C: List[float],

    # State Inputs full leaf pop
    D_0: float,
    g_bv: float,
    f_SW: float,
    f_VPD: float,

    f_LS: float,
    fO3_d: float,

    # external state
    eact: float,
    c_a: float,


    # model_options
    f_VPD_method: FVPDMethods,
    co2_concentration_balance_threshold: float = 0.001,
    co2_concentration_max_iterations: int = 50,
) -> Output_Shape:
    """Run the Ewert Photosynthesis model.

    On a per leaf population basis.

    NOTE: Be careful of units if modifying

    Parameters
    ----------
    nL: int
        Number of layers
    g_sto_0: float
        Closed stomata conductance [umol/m^2/s CO2]
        The stomatal conductance when A_n -> 0 and incident irradiance, I -> 0
    m: float
        Species-specific sensitivity to An [dimensionless]
    V_cmax_25: List[float]
        Maximum catalytic rate at 25 degrees [umol/m^2/s]
    J_max_25: List[float]
        Maximum rate of electron transport at 25 degrees [umol/m^2/s]
    R_d_coeff: float
        Dark respiration coefficient [Fraction] (Clark et al 2011)

    PARsun: List[float]
        PAR received by sunlit leaves [W m-2]
    PARshade: List[float]
        PAR received by shaded leaves [W m-2]
    D_0: float
        "The VPD at which g_sto is reduced by a factor of 2" [kPa] (Leuning et al. 1998)
    g_bv: float
        boundary layer conductance for forced convection [umol m-2 s-1 H2O]
    layer_lai_frac: List[Fraction]
        distribution of lai between layers as fraction of total layer LAI
    layer_lai: List[float]
        distribution of lai between layers
    Tleaf_C: List[float]
        Leaf Temperature [degrees celsius]
    f_SW: float
        Soil water influence on photosynthesis [0-1]
    f_VPD: float
        VPD effect on gsto [fraction] (Use if pre-calculated)
    f_LS: float
        factor effect of leaf senescence on A_c [Dimensionless][0-1]
    fO3_d: float
        Hourly accumulated ozone impace factor [dimensionless][0-1]

    eact: float
        Ambient vapour pressure [kPa]
    c_a: float
        CO2 concentration [ppm CO2]
    f_VPD_method: FVPDMethods
        If "photosynthesis" then f_VPD is calculated internally else should be provided as input
    co2_concentration_balance_threshold: float = 0.001
        Threshold (from 0) to consider co2 concentration equation as "balanced"
    co2_concentration_max_iterations: int = 50
        Maximum number of iterations to find co2 concentration solution

    Returns
    -------
        a named tuple of shape Output_Shape. See Output_Shape class for list of outputs

    """
    # TODO: Check how we can skip runs when lai is 0
    if sum(layer_lai_frac) == 0:
        return Output_Shape(
            g_sv_sunlit=0.0,
            g_sv_per_layer=[0.0 for _ in range(nL)],
            A_n=0.0,
            A_n_sunlit=0.0,
            A_c=0.0,
            A_j=0.0,
            A_p=0.0,
            A_n_limit_factor='NA',
            R_d=0.0,
            c_i=0.0,
            f_VPD=f_VPD,
            v_cmax=0.0,
            j_max=0.0,
        )
    layers_to_run = [i > 0 and (psun + pshade) > 0 for psun, pshade,
                     i in zip(PARsun, PARshade, layer_lai_frac)]

    # layers_to_run = [True for psun, pshade,
    #                  i in zip(PARsun, PARshade, layer_lai_frac)]

    # Output state of layers that have no sun or lai
    zero_state = CO2_loop_State(
        c_i=0,
        c_i_diff=0,
        g_sto=g_sto_0,
        A_n=0,
        A_c=0,
        A_p=0,
        A_j=0,
        A_n_limit_factor='NA',
        f_VPD=f_VPD,
        iterations=0,
    )

    model_options = ModelOptions(
        f_VPD_method,
        co2_concentration_balance_threshold=co2_concentration_balance_threshold,
        co2_concentration_max_iterations=co2_concentration_max_iterations,
    )

    e_a = eact * 1e3

    # ========== Run for PARsun
    inputs_sun: List[Ewert_Input_Factors] = [calc_input_factors(
        Tleaf_C=Tleaf_C[iL],
        Q=PARsun[iL] * 4.57,
        V_cmax_25=V_cmax_25[iL],
        J_max_25=J_max_25[iL],
        R_d_coeff=R_d_coeff,
    ) for iL in range(nL)]

    loop_inputs_sun: List[CO2_Constant_Loop_Inputs] = [CO2_Constant_Loop_Inputs(
        c_a=c_a,
        e_a=e_a,
        g_bl=g_bv,
        g_sto_0=g_sto_0,
        m=m,
        D_0=D_0,
        Gamma=inputs_sun[iL].Gamma,
        Gamma_star=inputs_sun[iL].Gamma_star,
        V_cmax=inputs_sun[iL].V_cmax,
        K_C=inputs_sun[iL].K_C,
        K_O=inputs_sun[iL].K_O,
        J=inputs_sun[iL].J,
        R_d=inputs_sun[iL].R_d,
        e_sat_i=inputs_sun[iL].e_sat_i,
        f_SW=f_SW,
        f_LS=f_LS,
        fO3_d=fO3_d,
        f_VPD=f_VPD,
    ) for iL in range(nL)]

    state_out_sun: List[CO2_loop_State] = [
        co2_concentration_in_stomata_loop(loop_inputs_sun[iP], model_options)
        if should_run else zero_state
        for iP, should_run in enumerate(layers_to_run)
    ]

    # ============= Run for PARshade

    inputs_shade: List[Ewert_Input_Factors] = [calc_input_factors(
        Tleaf_C=Tleaf_C[iL],
        Q=PARshade[iL] * 4.57,
        V_cmax_25=V_cmax_25[iL],
        J_max_25=J_max_25[iL],
        R_d_coeff=R_d_coeff,
    ) for iL in range(nL)]

    loop_inputs_shade: List[CO2_Constant_Loop_Inputs] = [replace(
        loop_inputs_sun[iL],
        Gamma=inputs_shade[iL].Gamma,
        Gamma_star=inputs_shade[iL].Gamma_star,
        V_cmax=inputs_shade[iL].V_cmax,
        K_C=inputs_shade[iL].K_C,
        K_O=inputs_shade[iL].K_O,
        J=inputs_shade[iL].J,
        R_d=inputs_shade[iL].R_d,
        e_sat_i=inputs_shade[iL].e_sat_i,
    ) for iL in range(nL)]

    state_out_shade: List[CO2_loop_State] = [
        co2_concentration_in_stomata_loop(loop_inputs_shade[iP], model_options)
        if should_run else zero_state
        for iP, should_run in enumerate(layers_to_run)
    ]

    # Get sun shade parts
    sun_fracs = LAIsunfrac
    shade_fracs = [(1 - LAIsunfrac[iL]) for iL in range(nL)]

    # Get Total sun shade values
    # gsv is outputed as a mean value for each layer.
    g_sv_per_layer = [state_out_sun[iL].g_sto * sun_fracs[iL] +
                      state_out_shade[iL].g_sto * shade_fracs[iL] for iL in range(nL)]

    # NOTE: below upscaled to total leaf values.
    # We get the weighted average values for each layer using sun fractions as weights
    # We then multiply the average at each layer by the LAI of this leaf pop at that layer
    # We then sum the layers to get the canopy scaled actual values for this leaf pop
    A_n = sum([layer_lai[iL] * (state_out_sun[iL].A_n * sun_fracs[iL] +
                                state_out_shade[iL].A_n * shade_fracs[iL]) for iL in range(nL)])
    A_c = sum([layer_lai[iL] * (state_out_sun[iL].A_c * sun_fracs[iL] +
                                state_out_shade[iL].A_c * shade_fracs[iL]) for iL in range(nL)])
    A_j = sum([layer_lai[iL] * (state_out_sun[iL].A_j * sun_fracs[iL] +
                                state_out_shade[iL].A_j * shade_fracs[iL]) for iL in range(nL)])
    A_p = sum([layer_lai[iL] * (state_out_sun[iL].A_p * sun_fracs[iL] +
                                state_out_shade[iL].A_p * shade_fracs[iL]) for iL in range(nL)])
    c_i = sum([layer_lai[iL] * (state_out_sun[iL].c_i * sun_fracs[iL] +
                                state_out_shade[iL].c_i * shade_fracs[iL]) for iL in range(nL)])
    R_d = sum([layer_lai[iL] * (inputs_sun[iL].R_d * sun_fracs[iL] + inputs_shade[iL].R_d *
                                shade_fracs[iL]) for iL in range(nL)])

    # Log only outputs
    # TODO: Check how we upscale f_VPD(Note only used for logging)
    f_VPD_out = [state_out_sun[iL].f_VPD * LAIsunfrac[iL] +
                 state_out_shade[iL].f_VPD * shade_fracs[iL] for iL in range(nL)][0]

    v_cmax_out = sum([layer_lai_frac[iL] * (inputs_sun[iL].V_cmax * LAIsunfrac[iL] +
                                            inputs_shade[iL].V_cmax * shade_fracs[iL]) for iL in range(nL)])
    j_max_out = sum([layer_lai_frac[iL] * (inputs_sun[iL].J_max * LAIsunfrac[iL] +
                                           inputs_shade[iL].J_max * shade_fracs[iL]) for iL in range(nL)])
    A_n_sunlit = state_out_sun[nL - 1].A_n
    g_sv_sunlit = state_out_sun[nL - 1].g_sto
    # Other
    A_n_limit_factor = [state_out_sun[iL].A_n_limit_factor for iL in range(nL)]

    return Output_Shape(
        g_sv_per_layer=g_sv_per_layer,
        g_sv_sunlit=g_sv_sunlit,
        A_n=A_n,
        A_n_sunlit=A_n_sunlit,
        A_c=A_c,
        A_j=A_j,
        A_p=A_p,
        A_n_limit_factor=A_n_limit_factor,
        R_d=R_d,
        c_i=c_i,
        f_VPD=f_VPD_out,
        v_cmax=v_cmax_out,
        j_max=j_max_out,
    )

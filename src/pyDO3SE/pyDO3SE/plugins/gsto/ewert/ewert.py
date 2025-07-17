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

from dataclasses import replace
from typing import List
from deprecated import deprecated
from pyDO3SE.Config.ConfigEnums import FVPDMethods
from .enums import AdjustNegativeAnMethods

from data_helpers.functional_helpers import fp_while

from .ewert_helpers import (
    Ewert_Input_Factors,
    calc_CO2_supply,
    calc_humidity_defecit_fVPD,
    calc_input_factors,
    calc_stomatal_conductance,
    calc_CO2_assimilation_rate,
    calc_CO2_assimilation_rate_cubic,
)

from .types import (
    CO2_Constant_Loop_Inputs,
    CO2_Concentration_Args,
    CO2_loop_State,
    ModelOptions,
    Output_Shape,
)


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
    if model_options.adjust_negative_A_n == AdjustNegativeAnMethods.CLIP_IN_LOOP and co2_assimilation_rate_values.A_n < 0:
        co2_assimilation_rate_values = co2_assimilation_rate_values._replace(
            A_n=0,
            A_n_limit_factor="clip",
        )

    f_VPD = (
        calc_humidity_defecit_fVPD(
            g_sto_in=g_sto_in,
            e_a=constant_inputs.e_a,
            g_bl=constant_inputs.g_bl,
            e_sat_i=constant_inputs.e_sat_i,
            D_0=constant_inputs.D_0,
            fmin=constant_inputs.fmin,
            f_VPD_method=model_options.f_VPD_method,
        )
        if model_options.f_VPD_method in [FVPDMethods.LEUNING, FVPDMethods.DANIELSSON]
        else constant_inputs.f_VPD
    )

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


def co2_concentration_in_stomata_cubic(
    constant_inputs: CO2_Concentration_Args,
    model_options: ModelOptions,
) -> CO2_loop_State:
    """Single iteration to find CO2 concentration (c_i) and stomatal conductance (g_sto).

    by allowing CO2 to converge.


    Parameters
    ----------
    constant_inputs: CO2_Concentration_Args
        Inputs provided to CO2 concentration cubic equations
    model_options: ModelOptions
        Ewert photosynthesis model options

    Returns
    -------
    state: CO2_loop_State
        state after iteration

    """
    co2_assimilation_rate_values = calc_CO2_assimilation_rate_cubic(
        V_cmax=constant_inputs.V_cmax,
        Gamma_star=constant_inputs.Gamma_star,
        K_C=constant_inputs.K_C,
        K_O=constant_inputs.K_O,
        fO3_d=constant_inputs.fO3_d,
        f_LS=constant_inputs.f_LS,
        J=constant_inputs.J,
        R_d=constant_inputs.R_d,
        g_bl=constant_inputs.g_bl,
        g_sto_0=constant_inputs.g_sto_0,
        g_sto_prev=constant_inputs.g_sto_prev,
        c_a=constant_inputs.c_a,
        Gamma=constant_inputs.Gamma,
        P=constant_inputs.P,
        e_a=constant_inputs.e_a,
        e_sat_i=constant_inputs.e_sat_i,
        D_0=constant_inputs.D_0,
        fmin=constant_inputs.fmin,
        f_VPD=constant_inputs.f_VPD,
        m=constant_inputs.m,
        f_SW=constant_inputs.f_SW,
        f_VPD_method=model_options.f_VPD_method,
        # TODO: Remove this once resolved neg a_n issues
        adjust_negative_A_n=model_options.adjust_negative_A_n,
    )

    f_VPD = (
        calc_humidity_defecit_fVPD(
            g_sto_in=constant_inputs.g_sto_prev,
            e_a=constant_inputs.e_a,
            g_bl=constant_inputs.g_bl,
            e_sat_i=constant_inputs.e_sat_i,
            D_0=constant_inputs.D_0,
            fmin=constant_inputs.fmin,
            f_VPD_method=model_options.f_VPD_method,
        )
        if model_options.f_VPD_method in [FVPDMethods.LEUNING, FVPDMethods.DANIELSSON]
        else constant_inputs.f_VPD
    )
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

    # TODO: Delete this when negative A_n resolved
    f_VPD_alt = (
        calc_humidity_defecit_fVPD(
            g_sto_in=g_sto,
            e_a=constant_inputs.e_a,
            g_bl=constant_inputs.g_bl,
            e_sat_i=constant_inputs.e_sat_i,
            D_0=constant_inputs.D_0,
            fmin=constant_inputs.fmin,
            f_VPD_method=model_options.f_VPD_method,
        )
        if model_options.f_VPD_method in [FVPDMethods.LEUNING, FVPDMethods.DANIELSSON]
        else constant_inputs.f_VPD
    )

    co2_supply = calc_CO2_supply(
        A_n=co2_assimilation_rate_values.A_n,
        c_a=constant_inputs.c_a,
        g_sto=g_sto,
        g_bl=constant_inputs.g_bl,
    )

    new_state = CO2_loop_State(
        c_i=co2_supply,
        c_i_diff=0,
        g_sto=g_sto,
        A_n=co2_assimilation_rate_values.A_n,
        # NOTE: We add R_d back in for consistency with legacy outputs
        A_c=co2_assimilation_rate_values.A_c + constant_inputs.R_d,
        A_j=co2_assimilation_rate_values.A_j + constant_inputs.R_d,
        A_p=co2_assimilation_rate_values.A_p + constant_inputs.R_d,
        A_n_limit_factor=co2_assimilation_rate_values.A_n_limit_factor,
        iterations=1,
        f_VPD=f_VPD,
        f_VPD_alt=f_VPD_alt,
    )
    return new_state


def co2_concentration_in_stomata_loop(
    constant_inputs: CO2_Constant_Loop_Inputs,
    model_options: ModelOptions,
) -> CO2_loop_State:
    """Run an iterative loop to find cO2 concentration by allowing CO2 to converge.

    This works by calculating the stomatal conductance (g_sto) and CO2 assimilation rate (A_n)
    using the current c_i and g_sto values. It then calculates the new c_i based on the
    CO2 assimilation rate and g_sto. The new c_i is then used to calculate the new g_sto and
    CO2 assimilation rate. This process is repeated until the difference between the new and
    old c_i is less than the tolerance. The final state is returned.

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
        lambda state: co2_concentration_in_stomata_iteration(
            constant_inputs, state, model_options
        ),
        initial_state,
        max_iter=max_iterations,
    )

    return final_state


@deprecated(reason="use ewert_leaf_pop_cubic instead")
def ewert_leaf_pop(
    # Config Inputs
    nL: int,
    g_sto_0: float,
    m: float,
    V_cmax_25: List[float],
    J_max_25: List[float],
    R_d_coeff: float,
    fmin: float,
    # State Inputs per layer
    PARsun: List[float],
    PARshade: List[float],
    LAIsunfrac: List[float],
    layer_lai_frac: List[float],
    layer_lai: List[float],
    Tleaf_C: List[float],
    # State Inputs full leaf pop
    D_0: float,
    g_bv: list[float],
    f_SW: float,
    f_VPD: list[float],
    f_LS: float,
    fO3_d: float,
    # external state
    eact: float,
    c_a: float,
    # model_options
    f_VPD_method: FVPDMethods,
    co2_concentration_balance_threshold: float = 0.001,
    co2_concentration_max_iterations: int = 50,
    adjust_negative_A_n: AdjustNegativeAnMethods = AdjustNegativeAnMethods.FALSE,
    **kwargs, # Added to allow cubic args to pass through
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
    fmin: float
        Minimum value of f_vpd [0-1]
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
    adjust_negative_A_n: AdjustNegativeAnMethods = True,
        If True then allow negative A_n values, else return NaN
        If "clip" then clip negative A_n values to -R_d

    Returns
    -------
        a named tuple of shape Output_Shape. See Output_Shape class for list of outputs

    """
    # TODO: Check how we can skip runs when lai is 0
    if sum(layer_lai_frac) == 0 or f_LS == 0:
        return Output_Shape(
            g_sv_sunlit=0.0,
            g_sv_per_layer=[0.0 for _ in range(nL)],
            A_n=0.0,
            A_n_sunlit=0.0,
            A_c=0.0,
            A_j=0.0,
            A_p=0.0,
            A_n_limit_factor=["NA" for _ in range(nL)],
            R_d=0.0,
            c_i=0.0,
            c_i_sunlit=0.0,
            f_VPD=f_VPD or [1.0 for _ in range(nL)],
            # TODO: Remove when we resolve neg anet issue
            f_VPD_alt=f_VPD or [1.0 for _ in range(nL)],
            v_cmax=0.0,
            j_max=0.0,
            loop_iterations=-1,
        )
    layers_to_run = [
        i > 0 and (psun + pshade) > 0 for psun, pshade, i in zip(PARsun, PARshade, layer_lai_frac)
    ]

    # Output state of layers that have no sun or lai
    zero_state = CO2_loop_State(
        c_i=0,
        c_i_diff=0,
        g_sto=g_sto_0,
        A_n=0,
        A_c=0,
        A_p=0,
        A_j=0,
        A_n_limit_factor="NA",
        f_VPD=1.0,
        iterations=0,
    )

    model_options = ModelOptions(
        f_VPD_method,
        co2_concentration_balance_threshold=co2_concentration_balance_threshold,
        co2_concentration_max_iterations=co2_concentration_max_iterations,
        adjust_negative_A_n=adjust_negative_A_n,
    )

    e_a = eact * 1e3

    # ========== Run for PARsun
    inputs_sun: List[Ewert_Input_Factors] = [
        calc_input_factors(
            Tleaf_C=Tleaf_C[iL],
            Q=PARsun[iL] * 4.57,
            V_cmax_25=V_cmax_25[iL],
            J_max_25=J_max_25[iL],
            R_d_coeff=R_d_coeff,
        )
        for iL in range(nL)
    ]

    loop_inputs_sun: List[CO2_Constant_Loop_Inputs] = [
        CO2_Constant_Loop_Inputs(
            c_a=c_a,
            e_a=e_a,
            g_bl=g_bv[iL],
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
            fmin=fmin,
            f_SW=f_SW,
            f_LS=f_LS,
            fO3_d=fO3_d,
            f_VPD=f_VPD[iL],
        )
        for iL in range(nL)
    ]

    state_out_sun: List[CO2_loop_State] = [
        co2_concentration_in_stomata_loop(loop_inputs_sun[iP], model_options)
        if should_run
        else zero_state
        for iP, should_run in enumerate(layers_to_run)
    ]

    # ============= Run for PARshade

    inputs_shade: List[Ewert_Input_Factors] = [
        calc_input_factors(
            Tleaf_C=Tleaf_C[iL],
            Q=PARshade[iL] * 4.57,
            V_cmax_25=V_cmax_25[iL],
            J_max_25=J_max_25[iL],
            R_d_coeff=R_d_coeff,
        )
        for iL in range(nL)
    ]

    loop_inputs_shade: List[CO2_Constant_Loop_Inputs] = [
        replace(
            loop_inputs_sun[iL],
            Gamma=inputs_shade[iL].Gamma,
            Gamma_star=inputs_shade[iL].Gamma_star,
            V_cmax=inputs_shade[iL].V_cmax,
            K_C=inputs_shade[iL].K_C,
            K_O=inputs_shade[iL].K_O,
            J=inputs_shade[iL].J,
            R_d=inputs_shade[iL].R_d,
            e_sat_i=inputs_shade[iL].e_sat_i,
        )
        for iL in range(nL)
    ]

    state_out_shade: List[CO2_loop_State] = [
        co2_concentration_in_stomata_loop(loop_inputs_shade[iP], model_options)
        if should_run
        else zero_state
        for iP, should_run in enumerate(layers_to_run)
    ]

    # Get sun shade parts
    sun_fracs = LAIsunfrac
    shade_fracs = [(1 - LAIsunfrac[iL]) for iL in range(nL)]

    # Get Total sun shade values
    # gsv is outputed as a mean value for each layer.
    g_sv_per_layer = [
        state_out_sun[iL].g_sto * sun_fracs[iL] + state_out_shade[iL].g_sto * shade_fracs[iL]
        for iL in range(nL)
    ]

    # NOTE: below upscaled to total leaf values.
    # We get the weighted average values for each layer using sun fractions as weights
    # We then multiply the average at each layer by the LAI of this leaf pop at that layer
    # We then sum the layers to get the canopy scaled actual values for this leaf pop
    A_n = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_n * sun_fracs[iL] + state_out_shade[iL].A_n * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_c = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_c * sun_fracs[iL] + state_out_shade[iL].A_c * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_j = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_j * sun_fracs[iL] + state_out_shade[iL].A_j * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_p = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_p * sun_fracs[iL] + state_out_shade[iL].A_p * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    c_i = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].c_i * sun_fracs[iL] + state_out_shade[iL].c_i * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    R_d = sum(
        [
            layer_lai[iL]
            * (inputs_sun[iL].R_d * sun_fracs[iL] + inputs_shade[iL].R_d * shade_fracs[iL])
            for iL in range(nL)
        ]
    )

    # Log only outputs
    # TODO: Check how we upscale f_VPD(Note only used for logging)
    f_VPD_out = [
        (state_out_sun[iL].f_VPD or 0) * LAIsunfrac[iL]
        + (state_out_shade[iL].f_VPD or 0) * shade_fracs[iL]
        for iL in range(nL)
    ]

    # TODO: Delete this when negative A_n resolved
    f_VPD_alt_out = [
        (state_out_sun[iL].f_VPD_alt or 0) * LAIsunfrac[iL]
        + (state_out_shade[iL].f_VPD_alt or 0) * shade_fracs[iL]
        for iL in range(nL)
    ]

    v_cmax_out = sum(
        [
            layer_lai_frac[iL]
            * (inputs_sun[iL].V_cmax * LAIsunfrac[iL] + inputs_shade[iL].V_cmax * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    j_max_out = sum(
        [
            layer_lai_frac[iL]
            * (inputs_sun[iL].J_max * LAIsunfrac[iL] + inputs_shade[iL].J_max * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_n_sunlit = state_out_sun[nL - 1].A_n
    g_sv_sunlit = state_out_sun[nL - 1].g_sto
    c_i_sunlit = state_out_sun[nL - 1].c_i
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
        c_i_sunlit=c_i_sunlit,
        f_VPD=f_VPD_out,
        # TODO: Delete this when negative A_n resolved
        f_VPD_alt=f_VPD_alt_out,
        v_cmax=v_cmax_out,
        j_max=j_max_out,
        loop_iterations=max(s.iterations for s in state_out_sun + state_out_shade),
    )


def ewert_leaf_pop_cubic(
    # Config Inputs
    nL: int,
    g_sto_0: float,
    m: float,
    V_cmax_25: List[float],
    J_max_25: List[float],
    R_d_coeff: float,
    fmin: float,
    # State Inputs per layer
    PARsun: List[float],
    PARshade: List[float],
    LAIsunfrac: List[float],
    layer_lai_frac: List[float],
    layer_lai: List[float],
    Tleaf_C: List[float],
    g_sto_prev: List[float],
    P: float,
    # State Inputs full leaf pop
    D_0: float,
    g_bv: list[float],
    f_SW: float,
    f_VPD: list[float],
    f_LS: float,
    fO3_d: float,
    # external state
    eact: float,
    c_a: float,
    # model_options
    f_VPD_method: FVPDMethods,
    co2_concentration_balance_threshold: float = 0.001,
    co2_concentration_max_iterations: int = 1,
    adjust_negative_A_n: AdjustNegativeAnMethods = AdjustNegativeAnMethods.FALSE,
) -> Output_Shape:
    """Run the Ewert Photosynthesis model.

    On a per leaf population basis.

    NOTE: Be careful of units if m
        f_VPD_alt=f_VPD_alt_out,
        # TODO: Delete this when negative A_n resolved

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
    fmin: float
        Minimum value of f_vpd [0-1]

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
    RH: float
        Relative humidity [%]
    P: float
        Air pressure [kPa]
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
    adjust_negative_A_n: AdjustNegativeAnMethods = AdjustNegativeAnMethods.FALSE,
        If True then allow negative A_n values, else return NaN
        If "last_resort" then allow negative A_n values if no other solution is found
        If "clip" then clip negative A_n values to -R_d


    Returns
    -------
        a named tuple of shape Output_Shape. See Output_Shape class for list of outputs

    """
    if sum(layer_lai_frac) == 0 or f_LS == 0:
        return Output_Shape(
            g_sv_sunlit=g_sto_0,
            g_sv_per_layer=[g_sto_0 for _ in range(nL)],
            A_n=0.0,
            A_n_sunlit=0.0,
            A_c=0.0,
            A_j=0.0,
            A_p=0.0,
            A_n_limit_factor=["NA" for _ in range(nL)],
            R_d=0.0,
            c_i=0.0,
            c_i_sunlit=0.0,
            f_VPD=f_VPD or [1.0 for _ in range(nL)],
            # TODO: Remove when we resolve neg anet issue
            f_VPD_alt=f_VPD or [1.0 for _ in range(nL)],
            v_cmax=0.0,
            j_max=0.0,
            loop_iterations=-1,
        )
    layers_to_run = [
        i > 0 and (psun + pshade) > 0 for psun, pshade, i in zip(PARsun, PARshade, layer_lai_frac)
    ]

    # Output state of layers that have no sun or lai
    zero_state = CO2_loop_State(
        c_i=0,
        c_i_diff=0,
        g_sto=g_sto_0,
        A_n=0,
        A_c=0,
        A_p=0,
        A_j=0,
        A_n_limit_factor="NA",
        f_VPD=1.0,
        f_VPD_alt=1.0,
        iterations=0,
    )

    model_options = ModelOptions(
        f_VPD_method,
        co2_concentration_balance_threshold=co2_concentration_balance_threshold,
        co2_concentration_max_iterations=co2_concentration_max_iterations,
        adjust_negative_A_n=adjust_negative_A_n,
    )

    e_a = eact * 1e3

    # ========== Run for PARsun
    inputs_sun: List[Ewert_Input_Factors] = [
        calc_input_factors(
            Tleaf_C=Tleaf_C[iL],
            Q=PARsun[iL] * 4.57,
            V_cmax_25=V_cmax_25[iL],
            J_max_25=J_max_25[iL],
            R_d_coeff=R_d_coeff,
        )
        for iL in range(nL)
    ]

    co2_concentration_inputs_sun: List[CO2_Concentration_Args] = [
        CO2_Concentration_Args(
            c_a=c_a,
            e_a=e_a,
            g_bl=g_bv[iL],
            g_sto_0=g_sto_0,
            g_sto_prev=g_sto_prev[iL],
            m=m,
            D_0=D_0,
            fmin=fmin,
            Gamma=inputs_sun[iL].Gamma,
            # TODO: We no longer use gamma_star is this correct?
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
            f_VPD=f_VPD[iL],
            P=P,
        )
        for iL in range(nL)
    ]

    state_out_sun: List[CO2_loop_State] = [
        co2_concentration_in_stomata_cubic(
            co2_concentration_inputs_sun[iL],
            model_options,
        )
        if should_run
        else zero_state
        for iL, should_run in enumerate(layers_to_run)
    ]

    # ============= Run for PARshade

    inputs_shade: List[Ewert_Input_Factors] = [
        calc_input_factors(
            Tleaf_C=Tleaf_C[iL],
            Q=PARshade[iL] * 4.57,
            V_cmax_25=V_cmax_25[iL],
            J_max_25=J_max_25[iL],
            R_d_coeff=R_d_coeff,
        )
        for iL in range(nL)
    ]

    co2_concentration_inputs_shade: List[CO2_Concentration_Args] = [
        replace(
            co2_concentration_inputs_sun[iL],
            Gamma=inputs_shade[iL].Gamma,
            # TODO: We no longer use gamma_star is this correct?
            # Gamma_star=inputs_shade[iL].Gamma_star,
            V_cmax=inputs_shade[iL].V_cmax,
            K_C=inputs_shade[iL].K_C,
            K_O=inputs_shade[iL].K_O,
            J=inputs_shade[iL].J,
            R_d=inputs_shade[iL].R_d,
            e_sat_i=inputs_shade[iL].e_sat_i,
        )
        for iL in range(nL)
    ]

    state_out_shade: List[CO2_loop_State] = [
        co2_concentration_in_stomata_cubic(co2_concentration_inputs_shade[iL], model_options)
        if should_run
        else zero_state
        for iL, should_run in enumerate(layers_to_run)
    ]

    # Get sun shade parts
    sun_fracs = LAIsunfrac
    shade_fracs = [(1 - LAIsunfrac[iL]) for iL in range(nL)]

    # Get Total sun shade values
    # gsv is outputed as a mean value for each layer.
    g_sv_per_layer = [
        state_out_sun[iL].g_sto * sun_fracs[iL] + state_out_shade[iL].g_sto * shade_fracs[iL]
        for iL in range(nL)
    ]

    # NOTE: below upscaled to total leaf values.
    # We get the weighted average values for each layer using sun fractions as weights
    # We then multiply the average at each layer by the LAI of this leaf pop at that layer
    # We then sum the layers to get the canopy scaled actual values for this leaf pop
    A_n = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_n * sun_fracs[iL] + state_out_shade[iL].A_n * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_c = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_c * sun_fracs[iL] + state_out_shade[iL].A_c * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_j = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_j * sun_fracs[iL] + state_out_shade[iL].A_j * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_p = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].A_p * sun_fracs[iL] + state_out_shade[iL].A_p * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    c_i = sum(
        [
            layer_lai[iL]
            * (state_out_sun[iL].c_i * sun_fracs[iL] + state_out_shade[iL].c_i * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    R_d = sum(
        [
            layer_lai[iL]
            * (inputs_sun[iL].R_d * sun_fracs[iL] + inputs_shade[iL].R_d * shade_fracs[iL])
            for iL in range(nL)
        ]
    )

    # Log only outputs
    # TODO: Check how we upscale f_VPD(Note only used for logging)
    f_VPD_out = [
        (state_out_sun[iL].f_VPD or 0) * LAIsunfrac[iL]
        + (state_out_shade[iL].f_VPD or 0) * shade_fracs[iL]
        for iL in range(nL)
    ]

    # TODO: Delete this when negative A_n resolved
    f_VPD_alt_out = [
        (state_out_sun[iL].f_VPD_alt or 0) * LAIsunfrac[iL]
        + (state_out_shade[iL].f_VPD_alt or 0) * shade_fracs[iL]
        for iL in range(nL)
    ]

    v_cmax_out = sum(
        [
            layer_lai_frac[iL]
            * (inputs_sun[iL].V_cmax * LAIsunfrac[iL] + inputs_shade[iL].V_cmax * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    j_max_out = sum(
        [
            layer_lai_frac[iL]
            * (inputs_sun[iL].J_max * LAIsunfrac[iL] + inputs_shade[iL].J_max * shade_fracs[iL])
            for iL in range(nL)
        ]
    )
    A_n_sunlit = state_out_sun[nL - 1].A_n
    g_sv_sunlit = state_out_sun[nL - 1].g_sto
    c_i_sunlit = state_out_sun[nL - 1].c_i

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
        c_i_sunlit=c_i_sunlit,
        f_VPD=f_VPD_out,
        # TODO: Delete this when negative A_n resolved
        f_VPD_alt=f_VPD_alt_out,
        v_cmax=v_cmax_out,
        j_max=j_max_out,
        loop_iterations=max(s.iterations for s in state_out_sun + state_out_shade),
    )

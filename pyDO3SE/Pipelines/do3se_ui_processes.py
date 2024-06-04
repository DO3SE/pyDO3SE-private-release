r"""
A modified version of default processes that matches more closely to DO3SE UI.
NOTE: The output of this now matches default_processes for multiplicative models

input calcs

## HOURLY PROCESSES ==============
# TODO: Check these are implemented correctly in external state setup
accumulate precip
calc_ustar
calc_r_PAR
Calc_sinB
SB_Calc_Rn
Calc_humidity

## GSTO SETUP
Calc_Flight - calc_f_light_process + calc_PARsun shade + MLMC_sunlit_LAI_process # OK
Calc_ftemp - calc_f_temp_process
Calc_fVPD - calc_f_VPD_process
SB_Calc_fO3 - f_O3_process
SB_Calc_fXWP - f_SW_process

## Resistance Calcs
SB_Calc_Ra - \
Calc_Rb -     \ - calc_resistance_model_process
Calc_Rgs -    /- & calc_leaf_resistance_model_process # TODO: Check
Calc_Rinc -  / - & multi_layer_r_model_to_single_H20 # TODO: Check

# GSTO CALCS
SB_Calc_Tleaf - ? Should input estimate method? (Uses de_boeck)
VPDcrit_prepare - Done in estate setup
SB_Calc_gsto - gsto_multiplicative_process
    calc_gsto_pn

VPDcrit_apply - gsto_multiplicative_process

## RESISTANCE CALCS
Calc_Rsto - \ - As above
Calc_Rsur - /

# CALC_SOIL_MOISTURE_CHANGES
SB_Calc_Es_blocked - check_soil_evaporation_blocked (Only has linear and fpaw method)
Calc_Penman_Monteith - penman_monteith_hourly


SB_Calc_LWP - ?
Calc_fLWP - ?

Calc_O3_Concentration - calc_canopy_ozone_concentration_process
Calc_Ftot - ?
Calc_Fst - calc_fst_process
Calc_AFstY - calc_POD_process
Calc_AOT40 - calc_POD_process

## DAILY PROCESSES ==============
Calc_LAI() - calc_canopy_LAI
SB_Calc_SAI() - calc_canopy_SAI
Calc_fphen() - f_phen_method_process
SB_Calc_leaf_fphen() - leaf_f_phen_process ?
 # ??
Calc_precip_acc() - accumulate_precipitation_process
C # OKalc_Penman_Monteith_daily() - penman_monteith_daily ?
Calc_SWP() - ??
SB_Calc_fSWP() - ??
Calc_SWP_meas() - ??
Calc_fPAW() - ?? f_SW_process??
"""

from copy import deepcopy
from typing import List
from data_helpers.list_helpers import flatten_list

from proflow.ProcessRunnerCls import advance_time_step_process
from proflow.Objects.Process import Process
from proflow.Objects.Interface import I
from proflow.Switch import switch
from proflow.helpers import set_value, skip
from proflow.logger import log_values

from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.Config.ConfigEnums import CanopyHeightMethods, FVPDMethods, LAIMethods
from do3se_phenology.config import FPhenMethods, LeafFPhenMethods

# Processes
from pyDO3SE.plugins.gsto import helpers as gsto_helpers
from do3se_phenology import f_phen
from pyDO3SE.plugins.gsto.multiplicative import multiplicative
from do3se_met import helpers as met_helpers
from do3se_met import irradiance as met_irrad_helpers
from do3se_met import wind as met_wind_helpers
from do3se_met.irradiance import calc_Idrctt_Idfuse, calc_PAR_sun_shade
from do3se_phenology import canopy_structure
# from pyDO3SE.plugins.soil_moisture import helpers as SMD_helpers
from pyDO3SE.plugins.soil_moisture import penman_monteith as SMD_PM_helpers
from pyDO3SE.plugins.O3 import helpers as O3_helpers
from pyDO3SE.plugins.O3.helpers import calc_fst

from pyDO3SE.Pipelines.validation_processes import setup_validation


def O3_ppb_to_nmol_process(iL: int) -> Process:
    """Convert O3 ppb to O3 nmol."""
    return Process(
        func=O3_helpers.O3_ppb_to_nmol,
        comment="Convert O3 ppb to O3 nmol",
        state_inputs=lambda state: [
            I(state.canopy_layers[iL].micro_met.micro_O3, as_='O3_ppb'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.Ts_C[row_index], as_='Ts_C'),
            I(e_state.P[row_index], as_='P'),
        ],
        state_outputs=lambda result: [
            (result, f'canopy_layers.{iL}.micro_met.micro_O3_nmol_m3')
        ]
    )


def calc_leaf_resistance_model_process(nLC: int, nL: int) -> Process:
    """Reset and setup Leaf Resistance Model."""
    return Process(
        func=O3_helpers.calc_leaf_resistance_model,
        comment="Reset and setup Leaf Resistance Model",
        config_inputs=lambda config: [
            I(config.Land_Cover.nL, as_='nL'),
            I(config.Land_Cover.nLC, as_='nLC'),
            I([config.Land_Cover.parameters[iLC].Lm for iLC in range(nLC)],
                as_='Lm_per_component')
        ],
        state_inputs=lambda state: [
            I([state.canopy_layers[iL].micro_met.micro_u for iL in range(nLC)],
                as_='u_per_layer'),
            I([[state.canopy_layer_component[iL][iLC].leaf_gsto
                for iLC in range(nLC)] for iL in range(nL)], as_='leaf_gsto_values'),
        ],
        state_outputs=lambda result: [
            (result[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.leaf_rmodel_O3')
            for iL in range(nL) for iLC in range(nLC)
        ]
    )


def MLMC_sunlit_LAI_process(nL: int, nLC: int) -> Process:
    """Calculate the sunlit and shaded LAI fractions."""
    return Process(
        # TODO: Check this matches Do3SE UI in f_light calc
        func=met_irrad_helpers.MLMC_sunlit_LAI,
        comment="Estimate sunlit LAI fractions",
        config_inputs=lambda config: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.nLC, as_='nLC'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.sinB[row_index], as_="sinB"),
        ],
        state_inputs=lambda state: [
            I([[state.canopy_layer_component[iL][iLC].LAI for iLC in range(nLC)]
               for iL in range(nL)], as_='LAI'),
        ],
        state_outputs=lambda result: [
            (result[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.LAIsunfrac')
            for iL in range(nL)
            for iLC in range(nLC)
        ],
    )


def calc_LAI_total_process(nL, nLC) -> Process:
    return Process(
        func=lambda LAI_values: sum([lai for LC in LAI_values for lai in LC]),
        comment="get LAI_total as sum of all layers and components",
        state_inputs=lambda state: [
                I([[state.canopy_layer_component[iL][iLC].LAI for iL in range(nL)]
                    for iLC in range(nLC)], as_='LAI_values')
        ],
        state_outputs=lambda result: [(result, 'canopy.LAI_total')],
    )


def calc_rmodel_h20_process() -> Process:
    return Process(
        func=SMD_PM_helpers.multi_layer_r_model_to_single_H20,
        comment="Adapt multi-layer O3 resistance model to single-layer H2O resistance",
        state_inputs=lambda state: [
            I(state.canopy.rmodel_O3, as_='rmodel'),
            I(state.met.ustar, as_='ustar'),
            I(state.canopy.LAI_total, as_='total_LAI'),
            I(state.canopy.SAI_total, as_='total_SAI'),
        ],
        state_outputs=lambda result: [(result, 'canopy.rmodel_H2O')],
    )


def calc_windspeed_parameters_process() -> Process:
    """Calculate estimated windspeed at canopy."""
    return Process(
        func=met_wind_helpers.calc_windspeed_parameters,
        comment="Calculate estimated windspeed at canopy",
        config_inputs=lambda config: [
                I(config.Location.OTC, as_='o_top_chamber'),
                I(config.Location.h_u, as_='loc_h_u'),
                I(config.Location.z_u, as_='loc_z_u'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.canopy_height, as_='h'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.u[row_index], as_="u"),
        ],
        additional_inputs=lambda: [
            I(0.1, as_='MIN_WINDSPEED'),
        ],
        state_outputs=lambda result: [
            (result.u_i, 'met.u_i'),
            (result.ustar, 'met.ustar'),
            (result.micro_u, 'canopy_layers.0.micro_met.micro_u'),  # micro_u at top layer
        ]
    )


def calc_canopy_height(nL: int, height_method: str) -> List[Process]:
    """Calculate the height of each layer.

    Canopy height is either an input or taken from the "primary land cover".
    """
    return [
        switch(
            gate=height_method,
            comment="Get canopy height",
            options={
                CanopyHeightMethods.CONSTANT: Process(
                    func=canopy_structure.height_method_constant,
                    comment="Use height of primary land cover if height method is constant",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.parameters[config.Land_Cover.primary_LC].height,
                          as_='primary_land_cover_height')
                    ],
                    state_outputs=lambda result: [(result, 'canopy.canopy_height')],
                ),
                CanopyHeightMethods.INPUT: Process(
                    func=canopy_structure.height_method_input,
                    comment="Do nothing if using input",
                ),
            }

        ),
        [
            Process(
                func=lambda canopy_height, layer_height: canopy_height * layer_height,
                comment=' '.join(["Calculate height of top of each canopy as a fraction of total",
                                  "canopy height"]),
                config_inputs=lambda config, iL=iL: [
                    I(config.Land_Cover.layer_height_frac[iL], as_='layer_height')],
                state_inputs=lambda state: [
                    I(state.canopy.canopy_height, as_='canopy_height'),
                ],
                state_outputs=lambda result, iL=iL: [
                    (result, f'canopy_layers.{iL}.layer_height'),
                ],
            ) for iL in range(nL)
        ],
    ]


# def soil_moisture_from_SWP_process():
#     return Process(
#         # TODO: SWP out is incorrect
#         func=SMD_helpers.soil_moisture_from_SWP,
#         comment="SWP",
#         config_inputs=lambda config: [
#             I(config.soil_moisture.soil_config, as_='soil_config'),
#             I(config.soil_moisture.PWP, as_='PWP'),
#             I(config.soil_moisture.root, as_='root_depth'),
#         ],
#         state_inputs=lambda state: [
#             I(state.canopy.SMD.SWP, as_='SWP'),
#         ],
#         state_outputs=lambda result: [
#             (result.Sn, 'canopy.SMD.Sn'),
#             (result.SWP, 'canopy.SMD.SWP'),
#             (result.ASW, 'canopy.SMD.ASW'),
#             (result.SMD, 'canopy.SMD.SMD'),
#         ]
#     )


def soil_moisture_PM_process():
    return Process(
        func=SMD_PM_helpers.PM_soil_moisture_calc,
        comment="P-M - soil moisture calc",
        config_inputs=lambda config: [
            I(config.soil_moisture.soil_config, as_='soil_config'),
            I(config.soil_moisture.PWP, as_='PWP'),
            I(config.soil_moisture.root, as_='root_depth'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.SMD.Sn, as_='Sn_in'),
            I(state.canopy.PM.Sn_diff, as_='Sn_diff'),
        ],
        state_outputs=lambda result: [
            # ASW Output here is incorrect
            (result.Sn, 'canopy.SMD.Sn'),
            (result.SWP, 'canopy.SMD.SWP'),
            (result.ASW, 'canopy.SMD.ASW'),
            (result.SMD, 'canopy.SMD.SMD'),
        ]
    )


def penman_monteith_daily_process():
    return Process(
        func=SMD_PM_helpers.penman_monteith_daily,
        comment="Daily soil water content update from accumulated \
            Penman-Monteith values",
        config_inputs=lambda config: [
            I(config.soil_moisture.root, as_='root_depth'),
            I(config.soil_moisture.run_off_fraction, as_='run_off_fraction'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.LAI_total, as_='LAI'),
            I(state.canopy.SMD.ASW, as_='ASW'),
            I(state.canopy.SMD.SMD, as_='SMD'),
            I(state.canopy.PM.precip_acc_prev_day, as_='pm_state_precip_acc'),
            I(state.canopy.PM.run_off_acc, as_='pm_state_run_off_acc'),
            I(state.canopy.PM.Ei_acc, as_='pm_state_Ei_acc'),
            I(state.canopy.PM.Eat_acc, as_='pm_state_Eat_acc'),
            I(state.canopy.PM.percolated_acc, as_='pm_state_percolated_acc'),
        ],
        state_outputs=lambda result: [
            (result.rain_input, 'canopy.PM.rain_input'),
            (result.run_off, 'canopy.PM.run_off'),
            (result.run_off_acc, 'canopy.PM.run_off_acc'),
            (result.effective_irrig, 'canopy.PM.effective_irrig'),
            (result.intercepted_evaporated, 'canopy.PM.intercepted_evaporated'),
            (result.evapotranspiration, 'canopy.PM.evapotranspiration'),
            (result.Sn_diff, 'canopy.PM.Sn_diff'),
            (result.percolated, 'canopy.PM.percolated'),
            (result.percolated_acc, 'canopy.PM.percolated_acc'),
        ],
    )


def penman_monteith_reset_process():
    return Process(
        func=SMD_PM_helpers.penman_monteith_reset,
        comment="reset penman monteith daily accumulators",
        state_outputs=lambda result: [
            (result.Ei_acc, 'canopy.PM.Ei_acc'),
            (result.Et_acc, 'canopy.PM.Et_acc'),
            (result.Es_acc, 'canopy.PM.Es_acc'),
            (result.Eat_acc, 'canopy.PM.Eat_acc'),
        ],
    )


def accumulate_precipitation_process() -> Process:
    """Stores the accumulated precip over a day."""
    return Process(
        func=lambda precip_acc_dd, precip_mm: precip_acc_dd + (precip_mm / 1000),
        comment="Accumulate precipitation",
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.precip[row_index], as_='precip_mm'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.PM.precip_acc_dd, as_='precip_acc_dd'),
        ],
        state_outputs=lambda result: [
            (result, 'canopy.PM.precip_acc_dd'),
        ]
    )


def store_accumulate_precipitation_process() -> Process:
    return Process(
        func=set_value,
        comment="Accumulate precipitation",
        state_inputs=lambda state: [
            I(state.canopy.PM.precip_acc_dd, as_='precip_acc_prev_day'),
        ],
        additional_inputs=lambda: [
            I(0, as_='precip_acc_dd')
        ],
        state_outputs=lambda result: [
            (result['precip_acc_dd'], 'canopy.PM.precip_acc_dd'),
            (result['precip_acc_prev_day'], 'canopy.PM.precip_acc_prev_day'),
        ]
    )


def leaf_f_phen_process(iL: int, iLC: int, leaf_f_phen_method) -> Process:
    """Calculate leaf_f_phen using specified method.

    Method options:
    - "disabled" - Not set
    - "f_phen" -
    - "day PLF" -
    - "leaf tt day PLF" -
    """
    return switch(
        gate=leaf_f_phen_method,
        comment="Choose leaf_f_phen_method",
        options={
            LeafFPhenMethods.DISABLED: Process(
                func=skip,
                comment="leaf_f_phen_method - disabled",
            ),
            LeafFPhenMethods.F_PHEN: Process(
                func=set_value,
                comment="leaf_f_phen_method - f_phen",
                state_inputs=lambda state: [
                    I(state.canopy_layer_component[iL][iLC].gsto_params.f_phen,
                        as_='leaf_f_phen'),
                ],
                state_outputs=lambda result: [
                    (result['leaf_f_phen'],
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
                ]
            ),
            LeafFPhenMethods.DAY_PLF: Process(
                func=f_phen.leaf_f_phen_PLF,
                comment="leaf_f_phen_method - day PLF",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_1,
                        as_='leaf_f_phen_1'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_2,
                        as_='leaf_f_phen_2'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_a,
                        as_='leaf_f_phen_a'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_b,
                        as_='leaf_f_phen_b'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_c,
                        as_='leaf_f_phen_c'),
                    I(config.Land_Cover.parameters[iLC].season.Astart, as_='Astart'),
                    I(config.Land_Cover.parameters[iLC].season.Aend, as_='Aend'),
                ],
                state_inputs=lambda state: [
                    I(state.temporal.dd, as_='dd'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
                ]
            ),
            LeafFPhenMethods.TT_DAY_PLF: Process(
                func=f_phen.tt_leaf_f_phen_PLF,
                comment="NOT IMPLEMENTED leaf_f_phen_method - leaf tt day PLF",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_1,
                        as_='leaf_f_phen_1'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_2,
                        as_='leaf_f_phen_2'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_a,
                        as_='leaf_f_phen_a'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_b,
                        as_='leaf_f_phen_b'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_c,
                        as_='leaf_f_phen_c'),
                    I(config.Land_Cover.parameters[iLC].season.Astart, as_='Astart'),
                    I(config.Land_Cover.parameters[iLC].season.Aend, as_='Aend'),
                ],
                state_inputs=lambda state: [
                    I(state.temporal.dd, as_='dd'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
                ]
            ),
        }
    )


def f_phen_method_process(iL: int, iLC: int, f_phen_method: str) -> Process:
    """Calculate f_phen using specified method.

    Method options:
    - "simple day PLF" -
    - "complex day PLF" -
    - "tt day PLF" -
    """
    return switch(
        gate=f_phen_method,
        comment="select a f_phen method",
        options={
            FPhenMethods.SIMPLE_DAY_PLF: Process(
                func=f_phen.f_phen_simple_PLF,
                comment="f_phen_simple_PLF",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].season.SGS, as_='SGS'),
                    I(config.Land_Cover.parameters[iLC].season.EGS, as_='EGS'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_1, as_='f_phen_1'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_4, as_='f_phen_4'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_a, as_='f_phen_a'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_c, as_='f_phen_c'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_e, as_='f_phen_e'),
                ],
                state_inputs=lambda state: [
                    I(state.temporal.dd, as_='dd'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
                ]
            ),
            # TODO: Below not checked against UI
            # "complex day PLF": Process(
            #     func=f_phen.f_phen_complex_PLF,
            #     comment="f_phen_complex_PLF",
            #     config_inputs=lambda config: [
            #         I(config.Land_Cover.parameters[iLC].season.SGS, as_='SGS'),
            #         I(config.Land_Cover.parameters[iLC].season.EGS, as_='EGS'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_1, as_='f_phen_1'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_2, as_='f_phen_2'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_3, as_='f_phen_3'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_4, as_='f_phen_4'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_a, as_='f_phen_a'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_b, as_='f_phen_b'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_c, as_='f_phen_c'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_d, as_='f_phen_d'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_e, as_='f_phen_e'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_limA, as_='f_phen_limA'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_limB, as_='f_phen_limB'),
            #     ],
            #     state_inputs=lambda state: [
            #         I(state.temporal.dd, as_='dd'),
            #     ],
            #     state_outputs=lambda result: [
            #         (result,
            #             f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
            #     ]
            # ),
            # "tt day PLF": Process(
            #     func=f_phen.tt_f_phen_simple_PLF,
            #     comment="tt_f_phen_simple_PLF",
            #     config_inputs=lambda config: [
            #         I(config.Land_Cover.parameters[iLC].season.SGS, as_='SGS'),
            #         I(config.Land_Cover.parameters[iLC].season.EGS, as_='EGS'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_1, as_='f_phen_1'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_4, as_='f_phen_4'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_a, as_='f_phen_a'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_c, as_='f_phen_c'),
            #         I(config.Land_Cover.parameters[iLC].multip_gsto.f_phen_e, as_='f_phen_e'),
            #     ],
            #     state_inputs=lambda state: [
            #         I(state.temporal.dd, as_='dd'),
            #     ],
            #     state_outputs=lambda result: [
            #         (result,
            #             f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
            #     ],
            # ),
        }
    )


def calc_canopy_SAI(nL: int, nLC: int, SAI_method: str, primary_SAI_method: str) -> List[Process]:
    """Calculate the SAI throughout the canopy.

    Supplying LAI and SAI for multi-layer multi-component models can be
    awkward.  To generate a nL*nLC array for LAI/SAI, three schemes are
    possible:

      - "input" - All LAI/SAI values are known and supplied as input.
      - "input total" - Total LAI/SAI is known and supplied as input in the top-left cell
        (e.g. LAI(1,1)).  Values are divided among layers and land covers
        according to a (normalised) nL*nLC fLAI array.
      - "estimate total" - Total LAI/SAI is estimated from the properties of the "primary land
        cover" and divided according to fLAI.
    """
    return [
        switch(
            gate=SAI_method,
            comment="calc_canopy_SAI",
            options={
                "input": Process(
                    func=skip,
                    comment="Do nothing if using input",
                ),
                "input total": [Process(
                    func=lambda SAI, fLAI: SAI * fLAI,
                    comment="Spread single SAI value to layers and LCs",
                    config_inputs=lambda config, iL=iL, iLC=iLC: [
                        I(config.Land_Cover.fLAI[iL][iLC], as_='fLAI'),
                    ],
                    state_inputs=lambda state: [
                        # TODO: Use primary land cover instead of 0.0
                        I(state.canopy_layer_component[0][0].SAI, as_='SAI'),
                    ],
                    state_outputs=lambda result, iL=iL, iLC=iLC: [
                        (result, f'canopy_layer_component.{iL}.{iLC}.SAI')
                    ]
                ) for iL in range(nL) for iLC in range(nLC)],
                "estimate total": [
                    switch(
                        gate=primary_SAI_method,
                        comment="Use primary land cover's estimate of total SAI",
                        options={
                            "LAI": Process(
                                func=lambda LAI_values: sum([LAI for LAI in LAI_values]),
                                comment="Calc_SAI Estimate_total - LAI",
                                state_inputs=lambda state: [
                                    I([[state.canopy_layer_component[jL][jLC].LAI
                                        for jL in range(nL)]
                                        for jLC in range(nLC)], as_='LAI_values'),
                                ],
                                state_outputs=lambda result: [
                                    (result, 'canopy_layer_component.0.0.SAI')
                                ]
                            ),
                            "forest": Process(
                                func=lambda clc: sum([c.LAI for c in clc]) + 1.0,
                                comment="Calc_SAI_Estimate_total - forest",
                                state_inputs=lambda state: [
                                    I([[state.canopy_layer_component[jL][jLC].LAI
                                        for jL in range(nL)]
                                        for jLC in range(nLC)], as_='LAI_values'),
                                ],
                                state_outputs=lambda result: [
                                    (result, 'canopy_layer_component.0.0.SAI')
                                ]
                            ),
                            "wheat": Process(
                                func=canopy_structure.SAI_wheat_and_LAI,
                                comment="Calc_SAI_Estimate_total - wheat",
                                config_inputs=lambda config: [
                                    I(config.Land_Cover.parameters[0].season.SGS, as_='SGS'),
                                    I(config.Land_Cover.parameters[0].season.EGS, as_='EGS'),
                                    I(config.Land_Cover.parameters[0].season.LAI_1, as_='LAI_1'),
                                    I(config.Land_Cover.parameters[0].season.LAI_2, as_='LAI_2'),
                                    I(config.Land_Cover.parameters[0].season.LAI_a, as_='LAI_a'),
                                    I(config.Land_Cover.parameters[0].season.LAI_b, as_='LAI_b'),
                                    I(config.Land_Cover.parameters[0].season.LAI_c, as_='LAI_c'),
                                    I(config.Land_Cover.parameters[0].season.LAI_1, as_='LAI_1'),
                                    I(config.Land_Cover.parameters[0].season.LAI_d, as_='LAI_d'),
                                ],
                                state_inputs=lambda state: [
                                    I(state.temporal.dd, as_='dd'),
                                ],
                                state_outputs=lambda result: [
                                    (result, 'canopy_layer_component.0.0.SAI')
                                ]
                            ),
                        }
                    ),
                    Process(
                        func=lambda fLAI, SAI, nL, nLC: [
                            [SAI * fLAI[iL][iLC] for iLC in range(nLC)] for iL in range(nL)],
                        comment="Spread single SAI value to layers and LCs",
                        config_inputs=lambda config: [
                            I(config.Land_Cover.fLAI, as_='fLAI'),
                            I(config.Land_Cover.nL, as_='nL'),
                            I(config.Land_Cover.nLC, as_='nLC'),
                        ],
                        state_inputs=lambda state: [
                            I(state.canopy_layer_component[0][0].SAI, as_='SAI')
                        ],
                        state_outputs=lambda result: [
                            (result[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.SAI')
                            for iL in range(nL) \
                            for iLC in range(nLC)
                        ]
                    ),
                ],
            },
        ),
        Process(
            func=lambda SAI_values: sum([SAI for LC in SAI_values for SAI in LC]),
            comment="get layer SAI_total",
            state_inputs=lambda state: [
                    I([[state.canopy_layer_component[iL][iLC].SAI for iL in range(nL)]
                       for iLC in range(nLC)], as_='SAI_values')
            ],
            state_outputs=lambda result: [(result, 'canopy.SAI_total')],
        ),
    ]


def calc_canopy_LAI(nL: int, nLC: int, LAI_method: str) -> List[Process]:
    """Calculate LAI change in canopy.

    Supplying LAI and SAI for multi-layer multi-component models can be
    awkward.  To generate a nL*nLC array for LAI/SAI, three schemes are
    possible:

      - All LAI/SAI values are known and supplied as input.
      - Total LAI/SAI is known and supplied as input in the top-left cell
        (e.g. LAI(1,1)).  Values are divided among layers and land covers
        according to a (normalised) nL*nLC fLAI array.
      - Total LAI/SAI is estimated from the properties of the "primary land
        cover" and divided according to fLAI.
    """
    return [
        switch(
            gate=LAI_method,
            comment="Get LAI",
            options={
                LAIMethods.INPUT: Process(
                    func=skip,
                    comment="Do nothing if using input",
                ),
                LAIMethods.CONSTANT: [
                    Process(
                        func=lambda LAI, fLAI: LAI * fLAI,
                        comment="Set LAI to constant value",
                        config_inputs=lambda config, iL=iL, iLC=iLC: [
                            I(config.Land_Cover.fLAI[iL][iLC], as_='fLAI'),
                            I(config.Land_Cover.LAI, as_='LAI'),
                        ],
                        state_outputs=lambda result, iL=iL, iLC=iLC: [
                            (result, f'canopy_layer_component.{iL}.{iLC}.LAI'),
                        ],
                    ) for iL in range(nL) for iLC in range(nLC)
                ],
                LAIMethods.INPUT_TOTAL: [
                    Process(
                        func=lambda LAI, fLAI: LAI * fLAI,
                        comment="Spread single LAI value to layers and LCs",
                        config_inputs=lambda config, iL=iL, iLC=iLC: [
                            I(config.Land_Cover.fLAI[iL][iLC], as_='fLAI'),
                        ],
                        state_inputs=lambda state: [
                            # TODO: Use primary land cover instead of 0.0
                            I(state.canopy_layer_component[0][0].LAI, as_='LAI'),
                        ],
                        state_outputs=lambda result, iL=iL, iLC=iLC: [
                            (result, f'canopy_layer_component.{iL}.{iLC}.LAI'),
                        ]
                    ) for iL in range(nL) for iLC in range(nLC)
                ],
                LAIMethods.DVI_LIMITED_CONSTANT: [
                    Process(
                        func=lambda LAI, fLAI, dvi: LAI * fLAI if 0 < dvi < 2 else 0,
                        comment="Set LAI to constant value",
                        config_inputs=lambda config, iL=iL, iLC=iLC: [
                            I(config.Land_Cover.fLAI[iL][iLC], as_='fLAI'),
                            I(config.Land_Cover.LAI, as_='LAI'),
                        ],
                        state_inputs=lambda state: [
                            I(state.canopy_component[0].dvi, as_='dvi'),
                        ],
                        state_outputs=lambda result, iL=iL, iLC=iLC: [
                            (result, f'canopy_layer_component.{iL}.{iLC}.LAI'),
                        ],
                    ) for iL in range(nL) for iLC in range(nLC)
                ],
                LAIMethods.ESTIMATE_TOTAL: [
                    Process(
                        func=canopy_structure.LAI_method_estimate_total,
                        comment="""Use primary land cover's estimate of total LAI and spread over
                        layers and LCs""",
                        config_inputs=lambda config: [
                            I(config.Land_Cover.parameters[0].season.EGS, as_='EGS'),
                            I(config.Land_Cover.parameters[0].season.LAI_1, as_='LAI_1'),
                            I(config.Land_Cover.parameters[0].season.LAI_2, as_='LAI_2'),
                            I(config.Land_Cover.parameters[0].season.LAI_a, as_='LAI_a'),
                            I(config.Land_Cover.parameters[0].season.LAI_b, as_='LAI_b'),
                            I(config.Land_Cover.parameters[0].season.LAI_c, as_='LAI_c'),
                            I(config.Land_Cover.parameters[0].season.LAI_1, as_='LAI_1'),
                            I(config.Land_Cover.parameters[0].season.LAI_d, as_='LAI_d'),
                            I(config.Land_Cover.fLAI, as_='fLAI'),
                            I(config.Land_Cover.nL, as_='nL'),
                            I(config.Land_Cover.nLC, as_='nLC'),
                        ],
                        state_inputs=lambda state: [
                            I(state.temporal.dd, as_='dd'),
                            I(state.canopy_component[0].SGS, as_='SGS'),
                        ],
                        state_outputs=lambda result: [
                            (result[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.LAI')
                            for iL in range(nL) \
                            for iLC in range(nLC)
                        ]
                    ),
                ],
            }
        ),
        Process(
            func=lambda LAI_values: sum([lai for LC in LAI_values for lai in LC]),
            comment="get LAI_total as sum of all layers and components",
            state_inputs=lambda state: [
                    I([[state.canopy_layer_component[iL][iLC].LAI for iL in range(nL)]
                       for iLC in range(nLC)], as_='LAI_values')
            ],
            state_outputs=lambda result: [(result, 'canopy.LAI_total')],
        ),
        [Process(
            func=canopy_structure.calc_distribution_of_LAI_between_lcs,
            comment="Calculate the distribution of LAI between land covers",
            config_inputs=lambda config: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.nLC, as_='nLC'),
            ],
            state_inputs=lambda state: [
                I([[state.canopy_layer_component[jL][jLC].LAI for jL in range(nL)]
                   for jLC in range(nLC)], as_='LAI_values')
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result[iLC], f'canopy_component.{iLC}.LC_dist')
            ]
        ) for iLC in range(nLC)],
        Process(
            func=lambda MCs, LCPs: sum([(LC_dist * Lm) for Lm in LCPs for LC_dist in MCs]),
            comment="get layer LAI-weighted mean leaf width(Lm_LAI)",
            config_inputs=lambda config: [
                I([config.Land_Cover.parameters[iLC].Lm for iLC in range(nLC)], as_='LCPs'),
            ],
            state_inputs=lambda state: [
                I([state.canopy_component[iLC].LC_dist for iLC in range(nLC)], as_='MCs'),
            ],
            state_outputs=lambda result: [(result, 'canopy.Lm_LAI')],
        ),
    ]


def calc_PAR_sun_shade_process(iL, iLC) -> Process:
    def get_parSunShade_from_par(PAR, sinB, P, cosA, LAI):
        Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(sinB, P, PAR=PAR)
        return calc_PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, LAI)

    return Process(
        func=get_parSunShade_from_par,
        comment="Calc par sun shade",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].cosA, as_='cosA'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.sinB[row_index], as_='sinB'),
            I(e_state.PAR[row_index], as_='PAR'),
            I(e_state.P[row_index], as_='P'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.LAI_total, as_='LAI'),
        ],

        state_outputs=lambda result: [
            (result.PARsun, f'canopy_layers.{iL}.micro_met.PARsun'),
            (result.PARshade, f'canopy_layers.{iL}.micro_met.PARshade'),
        ],

    )


def calc_f_light_process(iL: int, iLC: int) -> Process:
    """Calculate f_light."""
    return Process(
        func=gsto_helpers.calc_f_light_method,
        comment="Calculate f_light",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].multip_gsto.f_lightfac, as_='f_lightfac'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.sinB[row_index], as_='sinB'),
            I(e_state.PAR[row_index], as_='PAR'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.LAI_total, as_='LAI'),
            I(state.canopy_layers[iL].micro_met.PARsun, as_='PARsun'),  # Check this
            I(state.canopy_layers[iL].micro_met.PARshade, as_='PARshade'),  # Check this
            I(state.canopy_layer_component[iL][iLC].LAIsunfrac, as_='LAIsunfrac'),
        ],
        state_outputs=lambda result: [
            (result.f_light, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_light'),
            (result.leaf_f_light,
             f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_light'),
        ],
    )


def calc_f_temp_process(iL: int, iLC: int, f_temp_method: str) -> Process:
    """Calculate f_temp using specified method.

    Method options:
    - "disabled" - Not set
    - "default" -
    - "square high" -
    """
    return switch(
        gate=f_temp_method,
        comment="Choose f_temp Method",
        options={
            "disabled": Process(
                func=skip,
                comment="f_temp - disabled",
            ),
            "default": Process(
                func=gsto_helpers.calc_f_temp,
                comment="f_temp - default",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_min, as_='T_min'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_opt, as_='T_opt'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_max, as_='T_max'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy_layer_component[iL][iLC].gsto_params.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.Ts_C[row_index], as_='Ts_C'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_temp'),
                ],
            ),
            "square high": Process(
                func=gsto_helpers.calc_f_temp_square_high,
                comment="f_temp - square high",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_min, as_='T_min'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_opt, as_='T_opt'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_max, as_='T_max'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy_layer_component[iL][iLC].gsto_params.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.Ts_C[row_index], as_='Ts_C'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_temp'),
                ],
            ),
        },
    )


def calc_f_VPD_process(iL: int, iLC: int, f_VPD_method: FVPDMethods) -> Process:
    """Calculate f_VPD using specified method.

    Method options:
     - "disabled" - Not set
     - "linear" -
     - "log" -
    """
    return switch(
        gate=f_VPD_method,
        comment="Choose f_VPD method",
        options={
            FVPDMethods.DISABLED: Process(
                func=skip,
                comment="f_VPD method - disabled",
            ),
            FVPDMethods.LINEAR: Process(
                func=gsto_helpers.f_VPD_linear,
                comment="f_VPD_method - linear",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.VPD_max, as_='VPD_max'),
                    I(config.Land_Cover.parameters[iLC].gsto.VPD_min, as_='VPD_min'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.VPD[row_index], as_='VPD'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_VPD'),

                ]
            ),
            FVPDMethods.LOG: Process(
                func=gsto_helpers.f_VPD_log,
                comment="f_VPD_method - log",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.VPD[row_index], as_='VPD'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_VPD'),

                ]
            ),
        },
    )


def f_O3_process(iL: int, iLC: int, f_O3_method: str) -> Process:
    """Calculate f_O3 using specified method.

    Method options:
     - "disabled" - Not set
     - "wheat" - Using wheat calc
     - "potato" - Using Potato calc
     - "Steph wheat" - Using Steph Wheat calc
    """
    return switch(
        gate=f_O3_method,
        comment="Choose f_O3 method",
        options={
            "disabled": Process(
                func=skip,
                comment="f_O3 method - disabled",
            ),
            "wheat": Process(
                func=lambda POD_0: ((1 + (POD_0 / 11.5)**10)**(-1)),
                comment="f_O3_method - wheat",
                state_inputs=lambda state: [
                    I(state.canopy_layer_component[iL][iLC].POD_0, as_='POD_0'),
                ],
                state_outputs=lambda result: [
                    (result, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_O3'),
                ],
            ),
            "potato": Process(
                func=lambda AOT_0: ((1 + (AOT_0 / 40)**5)**(-1)),
                comment="f_O3_method - potato",
                state_inputs=lambda state: [
                    I(state.canopy_layer_component[iL][iLC].AOT_0, as_='AOT_0'),
                ],
                state_outputs=lambda result: [
                    (result, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_O3'),
                ],
            ),
            "Steph wheat": Process(
                func=lambda POD_0: ((1 + (POD_0 / 27)**10)**(-1)),
                comment="f_O3_method - Steph wheat",
                state_inputs=lambda state: [
                    I(state.canopy_layer_component[iL][iLC].POD_0, as_='POD_0'),
                ],
                state_outputs=lambda result: [
                    (result, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_O3'),
                ],
            ),
        },
    )


def f_SW_process(iL: int, iLC: int, f_SW_method: str) -> Process:
    """Calculate f_sw using specified method.

    Method options:
     - "disabled" - Not set
     - "fSWP exp" -
     - "fSWP linear" -
     - "fPAW" -
    """
    return switch(
        gate=f_SW_method,
        comment="Choose f_SW method",
        options={
            "disabled": Process(
                func=set_value,
                comment="f_SWP_method - Disabled",
                state_outputs=lambda result: [
                    (1.0, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
            "fSWP exp": Process(
                func=gsto_helpers.f_SWP_exp,
                comment="f_SWP_method - fSWP exp",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].multip_gsto.fSWP_exp_a, as_='a'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.fSWP_exp_b, as_='b'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.SWP, as_='SWP'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
            "fSWP linear": Process(
                func=gsto_helpers.f_SWP_linear,
                comment="f_SWP_method - linear",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.SWP_min, as_='SWP_min'),
                    I(config.Land_Cover.parameters[iLC].gsto.SWP_max, as_='SWP_max'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.SWP, as_='SWP'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
            "fPAW": Process(
                func=gsto_helpers.f_PAW,
                comment="f_SWP_method - fPAW",
                config_inputs=lambda config: [
                    I(config.soil_moisture.ASW_FC, as_='ASW_FC'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.ASW, as_='ASW'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
        },
    )


def calc_resistance_model_process(
    nL: int,
    nLC: int,
    ra_calc_method: str,
    rsur_calc_method: str,
) -> Process:
    """Reset and Calculate the resistance model for O3 over the target canopy.

    Parameters
    ----------
    nL : int
        number of model layers
    nLC : int
        number of model components
    ra_calc_method : str
        ra calc method. Options:
        - "simple" -
        - "heat_flux" -

    Returns
    -------
    Process
        [description]
    """
    return Process(
        func=O3_helpers.calc_resistance_model,
        comment="Reset and Calculate the resistance model for O3 over the target canopy",
        config_inputs=lambda config: [
            I(config.Land_Cover.nL, as_='nL'),
            I(config.Location.Rsoil, as_='Rsoil'),
            I(ra_calc_method, as_='ra_calc_method'),
            I(rsur_calc_method, as_='rsur_calc_method'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            # These are only required for Ra heat flux
            I(e_state.Ts_C[row_index], as_='Ts_C'),
            I(e_state.Hd[row_index], as_='Hd'),
            I(e_state.P[row_index], as_='P'),
        ] if ra_calc_method == 'heat_flux' else [],
        state_inputs=lambda state: [
            I(state.met.ustar, as_='ustar'),
            I(state.canopy.canopy_height, as_='canopy_height'),
            I([[state.canopy_layer_component[iL][iLC].SAI for iL in range(nL)]
               for iLC in range(nLC)], as_='SAI_values'),
            I([[state.canopy_layer_component[iL][iLC].LAI for iL in range(nL)]
               for iLC in range(nLC)], as_='LAI_values'),
            I([[state.canopy_layer_component[iL][iLC].bulk_gsto for iL in range(nL)]
               for iLC in range(nLC)], as_='bulk_gsto_values'),
        ],
        state_outputs=lambda result: [
            (result, 'canopy.rmodel_O3'),
        ]
    )


def gsto_multiplicative_process(iL: int, iLC: int) -> Process:
    """Calculate Gsto using the multiplicative method."""
    return [
        Process(
            func=multiplicative,
            comment="Calculate gsto - multiplicative",
            config_inputs=lambda config: [
                I(config.Land_Cover.parameters[iLC].gsto.VPD_crit, as_='VPD_crit'),
                I(config.Land_Cover.parameters[iLC].multip_gsto.gmax, as_='gmax'),
                I(config.Land_Cover.parameters[iLC].multip_gsto.gmorph, as_='gmorph'),
                I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.VPD_dd[row_index], as_='VPD_dd')
            ],
            state_inputs=lambda state: [
                I(state.canopy_layer_component[iL][iLC].leaf_gsto, as_='initial_leaf_gsto'),
                I(state.canopy_layer_component[iL][iLC].mean_gsto, as_='initial_mean_gsto'),
                I(state.canopy_layer_component[iL][iLC].gsto_params.f_phen, as_='f_phen'),
                I(state.canopy_layer_component[iL][iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
                I(state.canopy_layer_component[iL][iLC].gsto_params.f_light, as_='f_light'),
                I(state.canopy_layer_component[iL]
                  [iLC].gsto_params.leaf_f_light, as_='leaf_f_light'),
                I(state.canopy_layer_component[iL][iLC].gsto_params.f_temp, as_='f_temp'),
                I(state.canopy_layer_component[iL][iLC].gsto_params.f_VPD, as_='f_VPD'),
                I(state.canopy_layer_component[iL][iLC].gsto_params.f_SW, as_='f_SW'),
                I(state.canopy_layer_component[iL][iLC].gsto_params.f_O3, as_='f_O3'),
            ],
            state_outputs=lambda result: [
                (result.new_leaf_gsto, f'canopy_layer_component.{iL}.{iLC}.leaf_gsto'),
                (result.new_mean_gsto, f'canopy_layer_component.{iL}.{iLC}.mean_gsto')
            ]
        ),
        Process(
            func=lambda mean_gsto, LAI: mean_gsto * LAI,
            comment="Scale mean gsto up to bulk gsto",
            state_inputs=lambda state: [
                I(state.canopy_layer_component[iL][iLC].mean_gsto, as_='mean_gsto'),
                I(state.canopy_layer_component[iL][iLC].LAI, as_='LAI'),
            ],
            state_outputs=lambda result: [
                (result, f'canopy_layer_component.{iL}.{iLC}.bulk_gsto')
            ],
        )]


def check_soil_evaporation_blocked_process() -> Process:
    return Process(
        # TODO: this assumes the first land cover is the only one that matters
        func=SMD_PM_helpers.check_soil_evaporation_blocked,
        comment="Is soil evaporation blocked?",
        # TODO: Split this into function per method
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[0].gsto.f_SW_method, as_='f_SW_method'),
            I(config.Land_Cover.parameters[0].gsto.SWP_max, as_='SWP_max'),
            I(config.Land_Cover.parameters[0].gsto.fmin, as_='fmin'),
            I(config.Land_Cover.parameters[0].gsto.fSWP_exp_a, as_='fSWP_exp_a'),
            I(config.Land_Cover.parameters[0].gsto.fSWP_exp_b, as_='fSWP_exp_b'),
            I(config.soil_moisture.ASW_FC, as_='ASW_FC'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.SMD.SWP, as_='SWP'),
            I(state.canopy.SMD.ASW, as_='ASW'),
        ],
        state_outputs=lambda result: [
            (result, 'canopy.Es_blocked'),
        ],
    )


def penman_monteith_hourly_process() -> Process:
    return Process(
        func=SMD_PM_helpers.penman_monteith_hourly,
        comment="Hourly Penman-Monteith calculations for evaporation and transpiration",
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.Rn[row_index], as_='Rn_MJ'),
            I(e_state.P[row_index], as_='P_kPa'),
            I(e_state.Ts_C[row_index], as_='Ts_C'),
            I(e_state.esat[row_index], as_='esat_kPa'),
            I(e_state.eact[row_index], as_='eact_kPa'),
            I(e_state.VPD[row_index], as_='VPD_kPa'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.rmodel_H2O.nL, as_='rm_h2o_nL'),
            I(state.canopy.rmodel_H2O.Rb, as_='rm_h2o_Rb'),
            I(state.canopy.rmodel_H2O.Rinc[0], as_='rm_h2o_Rinc_l0'),
            I(state.canopy.rmodel_H2O.Rsto[0], as_='rm_h2o_Rsto_l0'),
            I(state.canopy.rmodel_H2O.Rgs, as_='rm_h2o_Rgs'),
            I(state.canopy.LAI_total, as_='LAI'),
            I(state.canopy.Es_blocked, as_='Es_blocked'),
            I(state.canopy.PM.Ei_acc, as_='pm_state_Ei_acc'),
            I(state.canopy.PM.Et_acc, as_='pm_state_Et_acc'),
            I(state.canopy.PM.Es_acc, as_='pm_state_Es_acc'),
            I(state.canopy.PM.Eat_acc, as_='pm_state_Eat_acc'),
        ],
        state_outputs=lambda result: [
            (result.Ei_acc, 'canopy.PM.Ei_acc'),
            (result.Et_acc, 'canopy.PM.Et_acc'),
            (result.Es_acc, 'canopy.PM.Es_acc'),
            (result.Eat_acc, 'canopy.PM.Eat_acc'),
            (result.Ei_hr, 'canopy.PM.Ei_hr'),
            (result.Et_hr, 'canopy.PM.Et_hr'),
            (result.Es_hr, 'canopy.PM.Es_hr'),
            (result.Eat_hr, 'canopy.PM.Eat_hr'),
        ],
    )


def calc_canopy_ozone_concentration_process(
    canopy_ozone_method: str,
) -> Process:
    """Calculate Top Layer Canopy ozone."""
    return switch(
        gate=canopy_ozone_method,
        comment="Calculate canopy_ozone_method",
        options={
            "multiplicative": Process(
                func=met_helpers.calc_canopy_ozone_concentration_legacy,
                comment="Calculate Top Layer Canopy ozone",
                config_inputs=lambda config: [
                    I(config.Location.z_O3, as_='z_O3'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.O3[row_index], as_='O3'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.canopy_height, as_='canopy_height'),
                    I(state.canopy.rmodel_O3.Rsur[0], as_='Rsur_top_layer'),
                    I(state.canopy.rmodel_O3.Rb, as_='Rb_top_layer'),
                    I(state.canopy.rmodel_O3.Ra, as_='Ra_top_layer'),
                    I(state.met.u_i, as_='u_i'),
                    I(state.met.ustar, as_='ustar'),
                ],
                state_outputs=lambda result: [
                    (result.O3_i, 'met.O3_i'),
                    (result.micro_O3, 'canopy_layers.0.micro_met.micro_O3'),
                ]),
        },
        default_option=Process(
            func=met_helpers.calc_canopy_ozone_concentration,
            comment="Calculate Top Layer Canopy ozone",
            config_inputs=lambda config: [
                I(config.Location.h_O3, as_='h_O3_in'),
                I(config.Location.z_O3, as_='z_O3'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.O3[row_index], as_='O3'),
            ],
            state_inputs=lambda state: [
                I(state.canopy.canopy_height, as_='canopy_height'),
                I(state.canopy.rmodel_O3.Rtotal[0], as_='Rtotal_top_layer'),
                I(state.met.u_i, as_='u_i'),
            ],
            state_outputs=lambda result: [
                (result.O3_i, 'met.O3_i'),
                (result.micro_O3, 'canopy_layers.0.micro_met.micro_O3'),
            ]
        )

    )


def calc_fst_process(iL: int, iLC: int) -> Process:
    """Calculate the O3up (fst)."""
    return Process(
        func=calc_fst,
        comment="Calculate the O3up (fst)",
        gate=True,
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].Lm, as_='Lm'),
        ],
        state_inputs=lambda state: [
            I(state.canopy_layers[iL].micro_met.micro_O3_nmol_m3, as_='O3_nmol_m3'),
            I(state.canopy_layers[iL].micro_met.micro_u, as_='uh'),
            I(state.canopy_layer_component[iL]
                [iLC].leaf_rmodel_O3.Rsto, as_='Rsto_l'),
            I(state.canopy_layer_component[iL]
                [iLC].leaf_rmodel_O3.Rext, as_='Rext'),
            I(state.canopy_layer_component[iL][iLC].leaf_gsto or 0, as_='Gsto_l'),
        ],
        state_outputs=lambda result: [
            (result, f'canopy_layer_component.{iL}.{iLC}.O3up'),
        ],
    )


def calc_POD_process(iL, iLC) -> Process:
    """Calculate the Phytotoxic Ozone Dose."""
    return Process(
        func=O3_helpers.calc_POD,
        comment="Calculate the Phytotoxic Ozone Dose",
        config_inputs=lambda config, iL=iL, iLC=iLC: [
            I(config.Land_Cover.parameters[iLC].Y, as_='Y', required=True),
        ],
        state_inputs=lambda state, iL=iL, iLC=iLC: [
            I(state.canopy_layer_component[iL][iLC].O3up, as_='fst', required=True),
            I(state.prev_hour.canopy_layer_component[iL][iLC].POD_0,
                as_='POD_0_prev', required=True),
            I(state.prev_hour.canopy_layer_component[iL][iLC].POD_Y,
                as_='POD_Y_prev', required=True),
        ],
        state_outputs=lambda result, iL=iL, iLC=iLC: [
            (result.POD_0, f'canopy_layer_component.{iL}.{iLC}.POD_0'),
            (result.POD_Y, f'canopy_layer_component.{iL}.{iLC}.POD_Y'),
        ]
    )


def set_hour(hr: int) -> Process:
    """Set the state hour."""
    return Process(
        func=set_value,
        comment="Set Hour",
        additional_inputs=lambda: [I(hr, as_='hr')],
        state_outputs=lambda result: [(result['hr'], 'temporal.hr')]
    )


def sync_row_index() -> Process:
    """Sync the internal state row index with the external state row index."""
    return Process(
        func=set_value,
        comment="Set data row index",
        external_state_inputs=lambda e_state, row_index: [
            I(row_index, as_="row_index"),
        ],
        state_outputs=lambda result: [
            (result['row_index'], 'temporal.row_index')
        ]
    )


def store_prev_state() -> List[Process]:
    """Store a copy of the previous State."""
    return [
        Process(
            func=lambda: Model_State_Shape(),
            comment="Reset previous state",
            state_outputs=lambda result: [
                (result, 'prev_hour'),
            ]
        ),
        Process(
            func=set_value,
            comment="Copy State",
            state_inputs=lambda state: [
                I(state.temporal, as_='temporal'),
                I(state.external_met, as_='external_met'),
                I(state.met, as_='met'),
                I(state.canopy, as_='canopy'),
                I(state.canopy_layers, as_='canopy_layers'),
                I(state.canopy_component, as_='canopy_component'),
                I(state.canopy_layer_component, as_='canopy_layer_component'),
            ],
            # Note we need to deep copy here to ensure we don't modify prev state
            # when updating current state
            state_outputs=lambda result: [
                (deepcopy(result['temporal']), 'prev_hour.temporal'),
                (deepcopy(result['external_met']), 'prev_hour.external_met'),
                (deepcopy(result['met']), 'prev_hour.met'),
                (deepcopy(result['canopy']), 'prev_hour.canopy'),
                (deepcopy(result['canopy_layers']), 'prev_hour.canopy_layers'),
                (deepcopy(result['canopy_component']), 'prev_hour.canopy_component'),
                (deepcopy(result['canopy_layer_component']), 'prev_hour.canopy_layer_component'),
            ]
        ),
    ]


def hourly_processes(config: Config_Shape, hr: int) -> List[Process]:
    """Take the current hour and returns a list of processes.

    to be ran in this hour
    """
    # nL = config.Land_Cover.nL
    nLC = config.Land_Cover.nLC

    # Options
    ra_calc_method = config.resistance.ra_calc_method
    rsur_calc_method = "single_layer"
    # Per component options
    f_temp_method = [
        config.Land_Cover.parameters[iLC].multip_gsto.f_temp_method for iLC in range(nLC)]
    f_VPD_method = [
        config.Land_Cover.parameters[iLC].gsto.f_VPD_method for iLC in range(nLC)]
    # f_SW_method = [config.Land_Cover.parameters[iLC].multip_gsto.f_SW_method for iLC in range(nLC)]
    f_O3_method = [config.Land_Cover.parameters[iLC].multip_gsto.f_O3_method for iLC in range(nLC)]

    return [
        set_hour(hr),
        sync_row_index(),
        daily_start_processes(config) if hr == 0 else [],

        # TODO: Add hourly todos
        accumulate_precipitation_process(),
        #  HOURLY PROCESSES
        calc_windspeed_parameters_process(),

        # Added calc sun shade
        # calc_Idrctt_Idfuse_process(),
        calc_PAR_sun_shade_process(0, 0),
        MLMC_sunlit_LAI_process(1, 1),

        calc_f_light_process(0, 0),  # TODO: Check par and lai frac inputs
        calc_f_temp_process(0, 0, f_temp_method[0]),  # OK

        calc_f_VPD_process(0, 0, f_VPD_method[0]),  # OK
        f_O3_process(0, 0, f_O3_method[0]),  # OK


        # SB_Calc_Ra # OK
        # Calc_Rb # OK
        # Calc_Rgs # OK # Check rsoil set correctly
        # Calc_Rinc # OK
        # \/\/\/
        # TODO: We have moved this after gsto calc
        # calc_resistance_model_process(1, 1, ra_calc_method),  # ok

        # SB_Calc_Tleaf  # TODO: Check this

        gsto_multiplicative_process(0, 0),  # OK TODO: Double check this
        calc_resistance_model_process(1, 1, ra_calc_method, rsur_calc_method),  # ok
        calc_rmodel_h20_process(),  # ADDED to get Rb_h20
        calc_leaf_resistance_model_process(1, 1),
        check_soil_evaporation_blocked_process(),  # OK f_SW_method
        penman_monteith_hourly_process(),  # OK


        # SB_Calc_LWP - ?? TODO: Find this
        # Calc_fLWP - ?? TODO: Find this

        calc_canopy_ozone_concentration_process("multiplicative"),

        # Calc_Ftot ?? TODO: Find this
        O3_ppb_to_nmol_process(0),
        calc_fst_process(0, 0),  # OK
        calc_POD_process(0, 0),  # OK

        # Store the previous hours state
        store_prev_state(),
        # TODO: This should come from config or output_values_map
        log_values(lambda state: [
            I(state.temporal.hr, as_='hr'),
            I(state.temporal.dd, as_='dd'),
            I(state.temporal.row_index, as_='row_index'),
            I(state.canopy.LAI_total, as_='lai'),
            I(state.canopy.SAI_total, as_='sai'),
            # Note rsto that is stored is Rsto bulk
            I(state.canopy.rmodel_O3.Rsto[0], as_='rsto_c'),
            I(state.canopy.rmodel_O3.Ra, as_='ra'),
            I(state.canopy.rmodel_O3.Rb, as_='rb'),
            I(state.canopy.rmodel_O3.Rsur[0], as_='rsur'),
            I(state.canopy.rmodel_O3.Rinc[0], as_='rinc'),
            I(state.canopy.rmodel_H2O.Rsto[0], as_='rsto_h2o'),
            I(state.canopy_component[0].dvi, as_='dvi'),
            I(state.external_met.photoperiod, as_='photoperiod'),
            I(state.canopy_layers[0].micro_met.PARsun, as_='PARsun'),
            I(state.canopy_layers[0].micro_met.PARshade, as_='PARshade'),
            I(state.canopy_layer_component[0][0].A_n, as_='A_n'),
            I(state.canopy_layer_component[0][0].O3up, as_='fst'),
            I(state.canopy_layer_component[0][0].leaf_rmodel_O3.Rsto, as_='rsto_l'),
            I(state.canopy_layer_component[0][0].leaf_gsto, as_='gsto_l'),
            I(state.canopy_layer_component[0][0].mean_gsto,
              as_='gsto'),
            I(state.canopy_layer_component[0][0].bulk_gsto,
              as_='gsto_bulk'),
            I(state.canopy_layer_component[0][0].fO3_d, as_='fO3_d'),
            I(state.canopy_layer_component[0][0].c_i, as_='c_i'),
            I(state.canopy_layer_component[0][0].td_dd, as_='td_dd'),
            I(state.canopy_layers[0].micro_met.micro_O3_nmol_m3, as_='o3_nmol_m3'),
            I(state.canopy_layers[0].micro_met.micro_O3, as_='o3_ppb'),
            I(state.met.O3_i, as_='o3_ppb_i'),
            I(state.canopy.canopy_height, as_='canopy_height'),
            I(state.canopy_layer_component[0][0].POD_Y, as_='pody'),
            I(state.canopy_layer_component[0][0].POD_0, as_='pod0'),
            I(state.canopy.PM.precip_acc_prev_day, as_='precip_acc'),
            I(state.canopy.SMD.SWP, as_='swp'),
            I(state.canopy.SMD.ASW, as_='asw'),
            I(state.canopy.SMD.SMD, as_='smd'),
            I(state.canopy.SMD.Sn, as_='sn'),

            I(state.canopy.PM.Ei_hr, as_='ei'),
            I(state.canopy.PM.Et_hr, as_='et'),
            I(state.canopy.PM.Es_hr, as_='es'),

            I(state.canopy.PM.Ei_acc, as_='ei_acc'),
            I(state.canopy.PM.Et_acc, as_='et_acc'),
            I(state.canopy.PM.Es_acc, as_='es_acc'),

            I(state.canopy_layer_component[0][0].gsto_params.f_phen, as_='f_phen'),
            I(state.canopy_layer_component[0][0].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
            I(state.canopy_layer_component[0][0].gsto_params.f_light, as_='f_light'),
            I(state.canopy_layer_component[0]
              [0].gsto_params.leaf_f_light, as_='leaf_f_light'),

            I(state.canopy_layer_component[0][0].gsto_params.f_temp, as_='f_temp'),
            I(state.canopy_layer_component[0][0].gsto_params.f_VPD, as_='f_VPD'),
            I(state.canopy_layer_component[0][0].gsto_params.f_SW, as_='f_SW'),
            I(state.canopy_layer_component[0][0].gsto_params.f_O3, as_='f_O3'),
            I(state.met.ustar, as_='ustar'),
            I(state.met.u_i, as_='u50'),
        ]),
        advance_time_step_process(),
        # Process(
        #     func=print_row_index('Row index end hour'),
        #     ptype=ProcessType.TIME,
        # ),
    ]


# def set_row_index() -> Process:
#     """Set the current data row index."""
#     return Process(
#         func=get_row_index,
#         comment="Set data row index",
#         state_inputs=lambda state: [
#             I(state.temporal.dd, as_='day'),
#             I(state.temporal.hr, as_='hour'),
#         ],
#         state_outputs=lambda result: [
#             (result, 'temporal.row_index')
#         ]
#     )


def set_day(dd: int = None) -> Process:
    """Set the state day."""
    return Process(
        func=set_value,
        additional_inputs=lambda: [I(dd, as_='dd')],
        state_outputs=lambda result: [(result['dd'], 'temporal.dd')]
    ) if dd else Process(
        func=set_value,
        external_state_inputs=lambda e_state, row_index: [I(e_state.dd[row_index], as_='dd')],
        state_outputs=lambda result: [(result['dd'], 'temporal.dd')]
    )


def daily_start_processes(config: Config_Shape) -> List[Process]:
    """Get processes to be ran at start of day."""

    # nL = config.Land_Cover.nL
    nLC = config.Land_Cover.nLC

    # Options
    height_method = config.Land_Cover.height_method
    LAI_method = config.Land_Cover.LAI_method
    SAI_method = config.Land_Cover.SAI_method
    primary_SAI_method = config.Land_Cover.parameters[0].season.SAI_method

    # Per component options
    f_phen_method = [
        config.Land_Cover.parameters[iLC].multip_gsto.f_phen_method for iLC in range(nLC)]
    leaf_f_phen_method = [
        config.Land_Cover.parameters[iLC].multip_gsto.leaf_f_phen_method for iLC in range(nLC)]
    f_SW_method = [config.Land_Cover.parameters[iLC].gsto.f_SW_method for iLC in range(nLC)]

    return [
        set_day(),
        calc_canopy_height(1, height_method),  # ADDED
        # calc_windspeed_parameters_process(),
        calc_canopy_LAI(1, 1, LAI_method),  # OK
        calc_canopy_SAI(1, 1, SAI_method, primary_SAI_method),  # OK SAI_wheat_and_LAI
        calc_LAI_total_process(1, 1),

        # simple day PLF OK (Matches DO3SE UI Interface)
        f_phen_method_process(0, 0, f_phen_method[0]),
        leaf_f_phen_process(0, 0, leaf_f_phen_method[0]),  # simple day PLF OK

        # Moved to hourly calcs
        # accumulate_precipitation_process(),  # OK
        # Store accumulated Precip and reset daily count
        store_accumulate_precipitation_process(),
        penman_monteith_daily_process(),
        # soil_moisture_from_SWP_process(),  # OK
        soil_moisture_PM_process(),
        # TODO: Check where we reset PM state
        penman_monteith_reset_process(),

        f_SW_process(0, 0, f_SW_method[0]),  # TODO: Why do we call this again after hourly process
        # Calc_SWP_meas() - TODO: Check if this is needed
        # Calc_fPAW() - ?? f_SW_process??
    ]


def daily_end_processes(config: Config_Shape) -> List[Process]:
    """Take the current day and returns a list of processes to.

    # NOTE: Will not currently work because we are setting timestep in hour_processes
    # so this will be hour ahead

    be ran at the end of the day
    """
    return [

    ]


def daily_process_list(config: Config_Shape) -> List[Process]:
    """Take the current day and returns a list of processes to be ran for that day."""
    return [
        [hourly_processes(config, hr) for hr in range(24)],
    ]


def full_model_processes(
    config: Config_Shape,
    start_day: int = 1,
    end_day: int = 365,
) -> List[Process]:
    """Get a flattened list of all processes to be ran between dates."""
    return flatten_list([
        # init_state_processes(config),
        setup_validation(config),
        [daily_process_list(config) for _ in range(start_day, end_day + 1)],
    ])

"""Default state initialization processes
This takes the config and external state data and returns the
initialized state

"""
from dataclasses import replace
from typing import List
from proflow.Objects.Process import Process
from proflow.Objects.Interface import I
from proflow.helpers import set_value
from data_helpers.list_helpers import flatten_list
from do3se_phenology.config import ZeroDayOptions
from do3se_met.resistance.model import Resistance_Model
from do3se_met.soil_moisture.enums import SoilMoistureSource
from do3se_met.soil_moisture import penman_monteith as SMD_PM_helpers
from pyDO3SE import settings
from pyDO3SE.External_State.External_State_Config import ThermalTimeMethods
from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.Config import Config_Shape


def calculate_phenology_stages_process(config: Config_Shape, iLC: int):
    """Get the phenology processes defined by config.

    This is where we estimate the life span of the leaf and key stages.
    This can be calculated either using preset values from the config or estimating based
    on thermal time and thermal interval parameters.


    If using full season we first find the thermal time at emergence

    """

    return [
        # Process(
        #     func=set_value,
        #     comment="Set initial phenology state",
        #     config_inputs=lambda config, iLC=iLC: [
        #         I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl, as_='t_l'),
        #         I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_em, as_='t_lem'),
        #         I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_ma, as_='t_lma'),
        #         I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_ep, as_='t_lep'),
        #         I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_se, as_='t_lse'),
        #     ],
        #     state_outputs=lambda result, iLC=iLC: [
        #         (result['t_l'], f'canopy_component.{iLC}.t_l_estimate'),
        #         (result['t_lem'], f'canopy_component.{iLC}.t_lem'),
        #         (result['t_lma'], f'canopy_component.{iLC}.t_lma'),
        #         (result['t_lep'], f'canopy_component.{iLC}.t_lep'),
        #         (result['t_lse'], f'canopy_component.{iLC}.t_lse'),
        #     ],
        # ),
        Process(
            # As thermal time is 0 at Astart thermal time at 0 is negative
            # td_dd should == t_Astart at the astart row.
            func=lambda t_Astart, td_0: t_Astart + td_0,
            comment="Set initial td_dd if zero day is after emergence",
            gate=config.Land_Cover.phenology_options.zero_day == ZeroDayOptions.ASTART,
            config_inputs=lambda config: [
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.emerg_to_astart, as_='t_Astart'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.td[0], as_='td_0'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'prev_hour.canopy_component.{iLC}.td_dd'),
                (result, f'canopy_component.{iLC}.td_dd'),
            ],
        ),
        Process(
            func=set_value,
            comment="Set initial td_dd to 0",
            gate=config.Land_Cover.phenology_options.zero_day in [
                ZeroDayOptions.SOWING, ZeroDayOptions.EMERGENCE],
            additional_inputs=lambda: [
                I(0, as_='td_0'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result['td_0'], f'prev_hour.canopy_component.{iLC}.td_dd'),
                (result['td_0'], f'canopy_component.{iLC}.td_dd'),
            ],
        ),
        Process(
            func=set_value,
            comment="Set initial vernalised_thermal_time",
            gate=config.Met.thermal_time_method == ThermalTimeMethods.EXTERNAL,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.td[0], as_='td_0'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result['td_0'], f'prev_hour.canopy_component.{iLC}.td_v'),
                (result['td_0'], f'canopy_component.{iLC}.td_v'),
            ],
        ),
        Process(
            func=set_value,
            comment="Set initial thermal_time",
            gate=config.Met.thermal_time_method == ThermalTimeMethods.HOURLY,
            config_inputs=lambda config: [
                I(config.Met.thermal_time_offset, as_='td_offset')
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result['td_offset'], f'prev_hour.canopy_component.{iLC}.td'),
                (result['td_offset'], f'canopy_component.{iLC}.td'),
                (result['td_offset'], f'prev_hour.canopy_component.{iLC}.td_v'),
                (result['td_offset'], f'canopy_component.{iLC}.td_v'),
            ],
        ),
    ]


def setup_initial_previous_state(config: Config_Shape) -> List[Process]:
    """Setup the previous state for first pass."""
    return [
        Process(
            func=lambda: Model_State_Shape(),
            comment="Reset previous state",
            state_outputs=lambda result: [
                (result, 'prev_hour'),
            ],
        ),
    ]


def soil_moisture_setup(config: Config_Shape) -> List[Process]:
    return [
        Process(
            func=set_value,
            comment="Set initial Sn to 0",
            additional_inputs=lambda: [I(0, as_='Sn')],
            state_outputs=lambda result: [(result['Sn'], 'canopy.SMD.Sn')]
        ),
        Process(
            func=SMD_PM_helpers.soil_moisture_from_SWC,
            gate=config.soil_moisture.source == SoilMoistureSource.P_M,
            comment="P-M - Set initial SMD state",
            config_inputs=lambda config: [
                I(config.soil_moisture.soil_config, as_='soil_config'),
                I(config.soil_moisture.PWP, as_='PWP'),
                I(config.soil_moisture.root, as_='root_depth'),
                I(config.soil_moisture.initial_SWC, as_='Sn_in'),
            ],
            state_outputs=lambda result: [
                (result.Sn, 'canopy.SMD.Sn'),
                (result.SWP, 'canopy.SMD.SWP'),
                (result.ASW, 'canopy.SMD.ASW'),
                (result.SMD, 'canopy.SMD.SMD'),
            ],
        ),
    ]


def variable_setup(config: Config_Shape) -> List[Process]:
    return [
        Process(
            func=set_value,
            config_inputs=lambda config: [
                I(config.Location.start_day, as_='dd'),
            ],
            state_outputs=lambda result: [
                (result['dd'], 'temporal.dd'),
            ],
        ),
        Process(
            func=lambda nL: Resistance_Model(nL),
            comment="Init O3 Resistance Model",
            config_inputs=lambda config: [I(config.Land_Cover.nL, as_='nL')],
            state_outputs=lambda result: [(result, 'canopy.rmodel_O3')]
        ),
        [
            Process(
                func=set_value,
                comment="Initialise gsto params",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.gmax, as_='gmax'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.gmorph, as_='gmorph'),
                ],
                state_outputs=lambda result: [
                    (result['fmin'], f'canopy_layer_component.{iL}.{iLC}.gsto_params.fmin'),
                    (result['gmax'], f'canopy_layer_component.{iL}.{iLC}.gsto_params.gmax'),
                    (result['gmorph'], f'canopy_layer_component.{iL}.{iLC}.gsto_params.gmorph'),
                ],
            )
            for iL in range(config.Land_Cover.nL)
            for iLC in range(config.Land_Cover.nLC)
        ],
    ]


def layer_pop_trim(config: Config_Shape) -> List[Process]:
    """Trim model state to correct number of populations and layers."""

    def trim_state(state, nLC, nL, nP):
        assert nP <= settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS, \
            f"Must set environment variable MAX_NUM_OF_LEAF_POPULATIONS={nP}. Currently: {settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS}"
        assert nL <= settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS, \
            f"Must set environment variable MAX_NUM_OF_CANOPY_LAYERS={nL}. Currently: {settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS}"
        assert nLC <= settings.global_settings.MAX_NUM_OF_CANOPY_COMPONENTS, \
            f"Must set environment variable MAX_NUM_OF_CANOPY_COMPONENTS={nLC}. Currently: {settings.global_settings.MAX_NUM_OF_CANOPY_COMPONENTS}"

        return replace(
            state,
            canopy_component_population=[
                [state.canopy_component_population[iLC][iP] for iP in range(nP)] for iLC in range(nLC)],
            canopy_layers=[state.canopy_layers[iL] for iL in range(nL)],
            canopy_component=[state.canopy_component[iLC] for iLC in range(nLC)],
            canopy_layer_component=[[state.canopy_layer_component[iL][iLC]
                                     for iLC in range(nLC)] for iL in range(nL)],
        )

    return [
        Process(
            func=trim_state,
            comment="Trim state",
            config_inputs=lambda config: [
                I(config.Land_Cover.nLC, as_='nLC'),
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.nP, as_='nP'),
            ],
            state_inputs=lambda state: [
                I(state, as_='state'),
            ],
            state_outputs=lambda result: [
                (result.temporal, 'temporal'),
                (result.external_met, 'external_met'),
                (result.met, 'met'),
                (result.canopy, 'canopy'),
                (result.canopy_layers, 'canopy_layers'),
                (result.canopy_component, 'canopy_component'),
                (result.canopy_layer_component, 'canopy_layer_component'),
                (result.canopy_component_population, 'canopy_component_population'),
            ],
        ),
    ]


def state_init_processes(
    config: Config_Shape,
) -> List[Process]:
    ''' Returns a flattened list of all processes to be ran over a annual cycle'''
    # rows = [get_row(dd, hr) for dd in range(start_day, end_day) for hr in range(24)]
    # row_count = len(rows)
    return flatten_list([
        layer_pop_trim(config),
        setup_initial_previous_state(config),
        [calculate_phenology_stages_process(config, iLC) for iLC in range(config.Land_Cover.nLC)],
        soil_moisture_setup(config),
        variable_setup(config),
    ])


# TODO: Add method to calculate from phyllochron
    # Process(
    #     func=calc_life_stages_from_phyllochron,
    #     comment="Calculate the life stages from the phyllochron",
    #     # TODO: Should only be set at emergence
    #     config_inputs=lambda config, iLC=iLC: [
    #         I(config.Land_Cover.parameters[iLC].pn_gsto.t_lse_constant,
    #             as_='t_lse_constant'),
    #     ],
    #     state_inputs=lambda state, iLC=iLC: [
    #         I(state.canopy_component[iLC].phyllochron, as_='phyllochron'),
    #     ],
    #     state_outputs=lambda result, iLC=iLC: [
    #         (result.t_l, f'canopy_component.{iLC}.t_l_estimate'),
    #         (result.t_lem, f'canopy_component.{iLC}.t_lem'),
    #         (result.t_lep, f'canopy_component.{iLC}.t_lep'),
    #         (result.t_lse, f'canopy_component.{iLC}.t_lse'),
    #         (result.t_lma, f'canopy_component.{iLC}.t_lma'),
    #     ]
    # ),
    # TODO: Calculate f_phen

    # TODO: Implement calculate_t_l_from_leaf_f_phen
    # Process(
    #     func=calculate_t_l_from_leaf_f_phen,
    #     gate=config.Land_Cover.parameters[iLC].pn_gsto.life_span_method == LifeSpanMethods.LEAF_F_PHEN,
    #     comment="Use leaf_f_phen input data to estimate the leaf phenology life stages.",
    #     external_state_inputs=lambda e_state, row_index: [
    #         I(e_state.leaf_fphen, as_='leaf_f_phen'),
    #         I(e_state.dd, as_='dd_full'),
    #         I(e_state.td, as_='td_full'),
    #     ],
    #     state_outputs=lambda result, iLC=iLC: [
    #         (result.t_l, f'canopy_component.{iLC}.t_l_estimate'),
    #         (result.t_lem, f'canopy_component.{iLC}.t_lem'),
    #         (result.t_lma, f'canopy_component.{iLC}.t_lma'),
    #         (result.t_lep, f'canopy_component.{iLC}.t_lep'),
    #         (result.t_lse, f'canopy_component.{iLC}.t_lse'),
    #         (result.AStart, f'canopy_component.{iLC}.SGS'),
    #         (result.AStart, f'canopy_component.{iLC}.Astart'),
    #         (-result.t_emerg, f'prev_hour.canopy_component.{iLC}.td_dd'),
    #         (-result.t_emerg, f'canopy_component.{iLC}.td_dd'),
    #     ]
    # ),
    # Process(
    #     func=set_value,
    #     # gate=config.Land_Cover.parameters[iLC].pn_gsto.life_span_method in [
    #     #     LifeSpanMethods.CONSTANT, LifeSpanMethods.JULES] or config.Land_Cover.dvi_method == DVIMethods.JULES,
    #     comment="Set the model life stage values to match config input",
    #     config_inputs=lambda config, iLC=iLC: [
    #         I(config.Land_Cover.parameters[iLC].season.sowing_day,
    #             as_='sowing_day'),
    #     ],
    #     state_outputs=lambda result, iLC=iLC: [
    #         (result['sowing_day'], f'canopy_component.{iLC}.sowing_day'),
    #     ]
    # ),
    # TODO: Implement constant life span method
    # Process(
    #     func=set_value,
    #     gate=config.Land_Cover.parameters[iLC].pn_gsto.life_span_method in [
    #         LifeSpanMethods.CONSTANT, LifeSpanMethods.JULES],
    #     comment="Set the model life stage values to match config input",
    #     config_inputs=lambda config, iLC=iLC: [
    #         I(config.Land_Cover.parameters[iLC].pn_gsto.t_l_estimate, as_='t_l'),
    #         I(config.Land_Cover.parameters[iLC].pn_gsto.t_lem, as_='t_lem'),
    #         I(config.Land_Cover.parameters[iLC].pn_gsto.t_lma, as_='t_lma'),
    #         I(config.Land_Cover.parameters[iLC].pn_gsto.t_lep, as_='t_lep'),
    #         I(config.Land_Cover.parameters[iLC].pn_gsto.t_lse, as_='t_lse'),
    #         I(config.Land_Cover.parameters[iLC].season.SGS,
    #             as_='SGS'),
    #         I(config.Land_Cover.parameters[iLC].season.Astart,
    #           as_='Astart'),
    #     ],
    #     state_outputs=lambda result, iLC=iLC: [
    #         (result['t_l'], f'canopy_component.{iLC}.t_l_estimate'),
    #         (result['t_lem'], f'canopy_component.{iLC}.t_lem'),
    #         (result['t_lma'], f'canopy_component.{iLC}.t_lma'),
    #         (result['t_lep'], f'canopy_component.{iLC}.t_lep'),
    #         (result['t_lse'], f'canopy_component.{iLC}.t_lse'),
    #         (result['SGS'], f'canopy_component.{iLC}.SGS'),
    #         (result['Astart'], f'canopy_component.{iLC}.Astart'),
    #     ]
    # ),

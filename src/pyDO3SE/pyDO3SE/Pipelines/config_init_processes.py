"""Default external state initialization processes.

This takes the config and external state data and returns the
initialized parameters

This includes calculating missing data such as total PAR from direct and diffuse PAR

Note: state_inputs and state_outputs refers to config here
"""

from dataclasses import asdict
from typing import List
from functools import wraps

from data_helpers.list_helpers import flatten_list
from proflow.Objects.Process import Process
from proflow.helpers import NOT_IMPLEMENTED_PROCESS, set_value, skip
from proflow.Objects.Interface import I

from do3se_phenology.switchboard import process_phenology_config
from do3se_phenology.config import ZeroDayOptions
from do3se_met.resistance import calc_displacement_and_roughness_parameters
from do3se_met.soil_moisture.penman_monteith import get_initial_SWC
from do3se_met.soil_moisture.helpers import get_soil_config, calc_ASW
from do3se_met.soil_moisture.enums import SoilMoistureSource
from pyDO3SE.Config import Config_Shape
from pyDO3SE.Config.ConfigLocation import LandCoverType
from pyDO3SE.util.logger import defaultLogger

def _partial(fn, *args, **kwargs):
    @wraps(fn)
    def _inner(*_args, **_kwargs):
        return fn(*[*args, *_args], **{**kwargs, **_kwargs})
    return _inner


def setup_phenology_parameters(config: Config_Shape) -> List[Process]:
    """Setup the phenology parameters from the provided parameters."""
    nLC = config.Land_Cover.nLC

    def get_thermal_time_base_and_offset(zero_day: ZeroDayOptions, sowing_day: int, Astart_day: int):
        if zero_day == ZeroDayOptions.SOWING:
            return sowing_day, 0
        elif zero_day == ZeroDayOptions.ASTART:
            return Astart_day, 0
        elif zero_day == ZeroDayOptions.DATA_START:
            return 0, 0

    return [
        [Process(
            func=_partial(process_phenology_config, logger=defaultLogger),
            comment="Estimate key phenology dates from input config and data.",
            state_inputs=lambda state, iLC=iLC: [
                I(state.Land_Cover.nP, as_="nP"),
                I(state.Land_Cover.phenology_options, as_="model_config"),
                I(state.Land_Cover.parameters[iLC].phenology, as_="species_config"),
                # TODO: Base temperature should be species specifig
                I(state.Met.td_base_temperature, as_='td_base_temperature'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state and asdict(e_state), as_="external_data"),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result[0], 'Land_Cover.phenology_options'),
                (result[1], f'Land_Cover.parameters.{iLC}.phenology'),
            ],
        ) for iLC in range(nLC)],
        # TODO: Below should be per iLC
        Process(
            func=get_thermal_time_base_and_offset,
            comment="Set td start and offset",
            state_inputs=lambda state: [
                I(state.Land_Cover.phenology_options.zero_day, as_="zero_day"),
                I(state.Land_Cover.parameters[0].phenology.key_dates.sowing, as_='sowing_day'),
                I(state.Land_Cover.parameters[0].phenology.key_dates.Astart, as_='Astart_day'),
            ],
            state_outputs=lambda result: [
                (result[0], 'Met.thermal_time_start'),
                (result[1], 'Met.thermal_time_offset'),
            ],
        )
    ]


def init_soil_config(config: Config_Shape):
    """Complete config setup of soil."""
    return [
        Process(
            func=get_soil_config,
            comment="Check that a Soil_t is valid",
            state_inputs=lambda state: [
                I(state.soil_moisture.soil_texture, as_='name'),
                I(state.soil_moisture.soil_config, as_='soil_in'),
            ],
            state_outputs=lambda result: [
                (result, 'soil_moisture.soil_config'),
            ]
        ),
        # TODO: Can we do this without modifiying config
        Process(
            func=calc_ASW,
            comment="Calculate available soil water at field capacity",
            gate=config.soil_moisture.ASW_FC is None, # Only run if not already set
            state_inputs=lambda state: [
                I(state.soil_moisture.soil_config, as_='soil_config'),
                I(state.soil_moisture.PWP, as_='PWP'),
                I(state.soil_moisture.root, as_='root_depth'),
            ],
            state_outputs=lambda result: [
                (result, 'soil_moisture.ASW_FC'),
            ]
        ),
        Process(
            func=get_initial_SWC,
            comment="Set up soil water data source parameters",
            gate=config.soil_moisture.source == SoilMoistureSource.P_M,
            state_inputs=lambda state: [
                I(state.soil_moisture.initial_SWC, as_='SWC_in'),
                I(state.soil_moisture.soil_config.FC, as_='FC'),
            ],
            state_outputs=lambda result: [(result, 'soil_moisture.initial_SWC')]
        )
    ]


def setup_height_parameters(config: Config_Shape) -> List[Process]:
    return [
        Process(
            func=calc_displacement_and_roughness_parameters,
            comment="calculate measured wind canopy displacement parameters",
            state_inputs=lambda state: [
                I(state.Location.h_u, as_="h"),
                I(state.Location.land_cover_type == LandCoverType.FOREST, as_="is_forest"),
            ],
            state_outputs=lambda result: [
                (result[0], 'Location.u_d'),
                (result[1], 'Location.u_z0'),
            ],
        ),
        Process(
            func=calc_displacement_and_roughness_parameters,
            comment="calculate measured ozone canopy displacement parameters",
            state_inputs=lambda state: [
                I(state.Location.h_O3, as_="h"),
                I(state.Location.land_cover_type == LandCoverType.FOREST, as_="is_forest"),
            ],
            state_outputs=lambda result: [
                (result[0], 'Location.O3_d'),
                (result[1], 'Location.O3_z0'),
            ],
        ),
    ]

def config_init_processes(
    config_in: Config_Shape = None,
) -> List[Process]:
    """Get a flattened list of all processes to be ran over an annual cycle."""
    return flatten_list([
        setup_phenology_parameters(config_in),
        init_soil_config(config_in),
        setup_height_parameters(config_in),
    ])

"""Default external state initialization processes.

This takes the config and external state data and returns the
initialized parameters

This includes calculating missing data such as total PAR from direct and diffuse PAR

Note: state_outputs refers to external_state here
"""

import warnings
from typing import List, Union
import numpy as np
from data_helpers.list_helpers import flatten_list
from proflow.helpers import assert_defined, skip
from proflow.Objects.Process import Process
from proflow.Objects.Interface import I
from proflow.Switch import switch

from thermal_time import calcs as td_functions
from do3se_met import helpers as met_helpers
from do3se_met.solar_position import calc_solar_elevation
from do3se_met.irradiance import (
    calc_is_daylight,
    calc_radiation_list,
    get_net_radiation,
)

from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.External_State.External_State_Config import Config_Met, Config_Met_Inputs, InputMethod
from pyDO3SE.constants.physical_constants import T0
from pyDO3SE.plugins.O3 import helpers as O3_helpers


# TODO: Calculate the thermal time gradient instead of the actual thermal time
# That way we are hourly independent
def calc_thermal_time_processes(config_met: Config_Met) -> List[Process]:
    """Calculate the thermal time using mean temps array."""
    return [
        Process(
            func=td_functions.calc_thermal_time_range,
            comment="Get thermal time from accumulated input temperature",
            gate=config_met.inputs.td_method == InputMethod.CALCULATED,
            config_inputs=lambda config: [
                # TODO: Base temperature should be species specifig
                I(config_met.td_base_temperature, as_='t_b'),
                I(config_met.thermal_time_start, as_='thermal_time_start'),
                I(config_met.thermal_time_offset, as_='thermal_time_offset'),
            ],
            state_inputs=lambda state: [
                # TODO: Should we use effective temperature here
                I(state.Ts_C, as_='Ts_C'),
                I(state.dd, as_='dd_list'),
            ],
            state_outputs=lambda result: [
                (result, 'td'),
            ],
        ),
    ]


def calc_met_params_list() -> List[Process]:
    """Run calculations on external state related to met."""
    return [
        Process(
            func=met_helpers.calc_humidity_list,
            comment="Calculate missing humidity data from input data",
            state_inputs=lambda state: [
                I(state.Ts_C, as_='Ts_C_list'),
                I(state.VPD, as_='VPD_list'),
                I(state.RH, as_='RH_list'),
            ],
            state_outputs=lambda result:[
                (result.esat, 'esat'),
                (result.RH, 'RH'),
                (result.eact, 'eact'),
                (result.VPD, 'VPD'),
            ]
        ),
        Process(
            func=lambda lat, lon, dd_list, hr_list: [calc_solar_elevation(
                lat, lon, dd, hr) for dd, hr in zip(dd_list, hr_list)],
            comment="Calculate the solar elevation(sinB) from location and time",
            config_inputs=lambda config: [
                I(config.Location.lat, as_='lat'),
                I(config.Location.lon, as_='lon'),
            ],
            state_inputs=lambda state: [
                I(state.dd, as_='dd_list'),
                I(state.hr, as_='hr_list'),
            ],
            state_outputs=lambda result:[
                (result, 'sinB')
            ],
        ),
        Process(
            func=calc_radiation_list,
            comment="Estimate missing radiation data",
            state_inputs=lambda state: [
                I(state.PAR, as_='PAR_list'),
                I(state.Idrctt, as_='Idrctt_list'),
                I(state.Idfuse, as_='Idfuse_list'),
                I(state.PPFD, as_='PPFD_list'),
                I(state.R, as_='R_list'),
                I(state.Rn, as_='Rn_list'),
                I(state.sinB, as_='sinB_list'),
                I(state.P, as_='P_list'),
                I(state.cloudfrac, as_='cloudfrac_list'),
            ],
            state_outputs=lambda result:[
                (result.PAR, 'PAR'),
                (result.Idrctt, 'Idrctt'),
                (result.Idfuse, 'Idfuse'),
                (result.PPFD, 'PPFD'),
                (result.R, 'R'),
            ],
        ),

        Process(
            func=lambda lat, lon, elev, albedo, Rn_in_list, dd_list, hr_list,
            sinB_list, R_list, Ts_C_list, eact_list:
            [get_net_radiation(Rn_in, lat, lon, elev, albedo, dd, hr, sinB, R, Ts_C, eact)
             for Rn_in, dd, hr, sinB, R, Ts_C, eact in
             zip(Rn_in_list, dd_list, hr_list, sinB_list, R_list, Ts_C_list, eact_list)],
            comment="Get net radiation",
            config_inputs=lambda config: [
                I(config.Location.lat, as_='lat'),
                I(config.Location.lon, as_='lon'),
                I(config.Location.elev, as_='elev'),
                I(config.Location.albedo, as_='albedo'),
            ],
            state_inputs=lambda state: [
                I(state.Rn, as_='Rn_in_list'),
                I(state.dd, as_='dd_list'),
                I(state.hr, as_='hr_list'),
                I(state.sinB, as_='sinB_list'),
                I(state.R, as_='R_list'),
                I(state.Ts_C, as_='Ts_C_list'),
                I(state.eact, as_='eact_list'),
            ],
            state_outputs=lambda result:[
                (result, 'Rn'),
            ],
        ),
        Process(
            func=lambda R_list: [calc_is_daylight(R) for R in R_list],
            comment="Calculate daylight hours",
            state_inputs=lambda state: [
                I(state.R, as_='R_list'),
            ],
            state_outputs=lambda result:[
                (result, 'is_daylight'),
            ]
        ),
        Process(
            func=met_helpers.calc_vpd_daily_accumulation_list,
            comment="Calculate daily accumulated VPD during daylight hours",
            state_inputs=lambda state: [
                I(state.hr, as_='hr_list'),
                I(state.VPD, as_='VPD_list'),
                I(state.is_daylight, as_='is_daylight_list'),
            ],
            state_outputs=lambda result:[
                (result, 'VPD_dd'),
            ],
        ),
    ]


def calc_CO2_params(config_met: Config_Met, row_index: int):
    return switch(
        gate=config_met.inputs.CO2_method,
        comment=' '.join(
            ["Define a CO2 concentration, either using a supplied value",
             "or a configured constant value."]),
        options={
            InputMethod.CONSTANT: Process(
                func=lambda CO2: CO2,
                comment="Set CO2 to constant value",
                config_inputs=lambda config: [I(config.Met.CO2_constant, as_='CO2')],
                state_outputs=lambda result: [(result, f'CO2.{row_index}')]
            ),
            InputMethod.INPUT: Process(
                func=assert_defined,
                state_inputs=lambda state: [I(state.CO2[row_index])]
            )
        }
    )


def calc_CO2_params_list(config_met: Config_Met, row_indexes: int):
    """Calc CO2 state from constant or input."""
    row_count = len(row_indexes)
    return switch(
        gate=config_met.inputs.CO2_method,
        comment=" ".join(["Define a CO2 concentration, either using a supplied",
                          "value or a configured constant value."]),
        options={
            InputMethod.CONSTANT: Process(
                func=lambda CO2: np.full((row_count,), CO2).tolist(),
                comment="Set CO2 to constant value",
                config_inputs=lambda config: [
                    I(config.Met.inputs.CO2_constant, as_='CO2')],
                state_outputs=lambda result: [
                    (result, 'CO2'),
                ]
            ),
            InputMethod.INPUT: Process(
                func=lambda CO2_list: [assert_defined(CO2) for CO2 in CO2_list],
                state_inputs=lambda state: [I('CO2', as_='CO2_list')],
                state_outputs=lambda result: [
                    (result, 'CO2'),
                ]
            )
        }
    )


def convert_units(config_met: Config_Met):
    return [
        Process(
            func=lambda Ts_K: np.array(Ts_K) - T0,
            comment="Convert kelvin to degrees C",
            gate=config_met.inputs.Ts_C_method == InputMethod.CALCULATED and config_met.inputs.Ts_K_method == InputMethod.INPUT,
            state_inputs=lambda state: [
                I(state.Ts_K, as_='Ts_K'),
            ],
            state_outputs=lambda result:[
                (result, 'Ts_C'),
            ]
        ),
    ]


def calc_O3_params_list(config_met: Config_Met, row_indexes: List[int]) -> List[Process]:
    """Run O3 calculations on config and external state."""
    row_count = len(row_indexes)
    return [
        switch(
            gate=config_met.inputs.O3_method,
            comment="Get O3 Concentration",
            options={
                InputMethod.INPUT: Process(
                    func=skip,
                    comment="Skip if input",
                ),
                InputMethod.CONSTANT: Process(
                    func=lambda O3: np.full((row_count,), O3).tolist(),
                    comment="Set O3 as constant",
                    config_inputs=lambda config: [
                        I(config.Met.inputs.O3_constant, as_='O3'),
                    ],
                    state_outputs=lambda result:[
                        (result, 'O3'),
                    ]
                ),
                InputMethod.OFFSET: Process(
                    func=lambda O3_list, O3_const: [max(0.0, O3 + O3_const) for O3 in O3_list],
                    comment="Set as O3 with constant offset",
                    config_inputs=lambda config: [
                        I(config.Met.inputs.O3_constant, as_='O3_const'),
                    ],
                    state_inputs=lambda state: [
                        I(state.O3, as_='O3_list'),
                    ],
                    state_outputs=lambda result:[
                        (result, 'O3'),
                    ]
                ),
            }),
        Process(
            func=lambda Ts_C_list, P_list, O3_ppb_list: [O3_helpers.O3_ppb_to_nmol(
                Ts_C, P, O3_ppb) for Ts_C, P, O3_ppb in zip(Ts_C_list, P_list, O3_ppb_list)],
            comment="Convert O3 ppb to O3 nmol",
            state_inputs=lambda state: [
                I(state.O3, as_='O3_ppb_list'),
                I(state.Ts_C, as_='Ts_C_list'),
                I(state.P, as_='P_list'),
            ],
            state_outputs=lambda result:[
                (result, 'O3_nmol')
            ],
        ),
    ]


def set_constants(
    k: str,
    config_met: Config_Met,
    row_count: int,
    input_vals: List[Union[float, int]],
):
    """If any inputs are set to constant we assign them here."""
    config_ext_keys = Config_Met_Inputs.__annotations__.keys()
    input_config = config_met.inputs
    if f'{k}_method' in config_ext_keys:
        method = getattr(input_config, f'{k}_method', None)
        value = getattr(input_config, f'{k}_constant', None)
        if method == InputMethod.SKIP:
            if input_vals is None:
                return None
                # return np.full((row_count,), None).tolist()
            else:
                return input_vals
        if method == InputMethod.INPUT:
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must supply {k} input')
            return input_vals
        if method == InputMethod.CONSTANT:
            if value is None:
                raise ValueError(f'Must supply {k}_constant')
            return np.full((row_count,), value).tolist()
        if method == InputMethod.OFFSET:
            if value is None:
                raise ValueError(f'Must supply {k}_constant')
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must supply {k} input')
            return (np.array(input_vals) + value).tolist()
        if method == InputMethod.MULTIPLY:
            if value is None:
                raise ValueError(f'Must supply {k}_constant')
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must supply {k} input')
            return (np.array(input_vals) * value).tolist()
        if method == InputMethod.CALCULATED:
            # Currently if it is calculated we set initially to None
            return np.full((row_count,), None).tolist()
        raise ValueError(f'Invalid input method {method} for {k}')
    else:
        if input_vals is None or input_vals[0] is None:
            raise ValueError(f'Must supply {k} input')
        return input_vals


def validate_inputs(
    k: str,
    config_met: Config_Met,
    input_vals: List[Union[float, int]],
):
    """Validate inputs have been provided."""
    config_ext_keys = Config_Met_Inputs.__annotations__.keys()
    input_config = config_met.inputs
    if f'{k}_method' in config_ext_keys:
        method = getattr(input_config, f'{k}_method', None)
        if method == InputMethod.SKIP:
            if input_vals is not None and input_vals[0] is not None:
                warnings.warn(f'Did you mean to provide {k} input as it is marked as skip')
            return
        if method == InputMethod.INPUT:
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must supply {k} input')
            return
        if method == InputMethod.CONSTANT:
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must supply {k} constant')
            return
        if method == InputMethod.OFFSET:
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must supply {k} constant')
            return
        if method == InputMethod.MULTIPLY:
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must supply {k} constant')
            return
        if method == InputMethod.CALCULATED:
            if input_vals is None or input_vals[0] is None:
                raise ValueError(f'Must calculate {k}')
            return
        raise ValueError(f'Invalid input method {method} for {k}')
    else:
        if input_vals is None or input_vals[0] is None:
            raise ValueError(f'Must supply {k} input')
        return input_vals


def validate_PAR_data(PAR_list, PPFD_list, Rn_list):
    if max(PAR_list) > 700:
        warnings.warn(
            f"max PAR is greater than 700. Are you sure this is not PPFD (Max is {max(PAR_list)})?")

    if max(PPFD_list) < 1300:
        warnings.warn(
            f"max PPFD is less than 1300. Are you sure this is not PAR or Rn (Max is {max(PPFD_list)})?")

    if max(Rn_list) > 1300:
        warnings.warn(
            f"max PPFD is less than 1300. Are you sure this is not PAR (Max is {max(Rn_list)})?")


def validate_vpd_data(VPD_list):
    if max(VPD_list) > 10:
        warnings.warn(
            f"max VPD is greater than 10. Check inputs. (Max is {max(VPD_list)})?")
    if min(VPD_list) < 0:
        warnings.warn(
            f"min VPD is less than 0. Check inputs. (Min is {min(VPD_list)})?")

def init_params(row_count: int) -> List[Process]:
    """Set data to input or constant data."""
    keys = [k for k in list(External_State_Shape.__annotations__.keys())]
    return [
        Process(
            func=set_constants,
            comment=f"initialize lists: {k}",
            additional_inputs=lambda k=k: [
                I(k, 'k'),
                I(row_count, 'row_count'),
            ],
            config_inputs=lambda config: [
                I(config.Met, as_='config_met'),
            ],
            state_inputs=lambda state, k=k: [
                I(getattr(state, k), as_='input_vals'),
            ],
            state_outputs=lambda result, k=k:[
                (result, k),
            ],
        ) for k in keys
    ]


def validate_input() -> List[Process]:
    """Validate that each required input has been set."""
    # TODO: Should this be moved to validation processes?
    keys = [k for k in list(External_State_Shape.__annotations__.keys())]
    return [
        [
            Process(
                func=validate_inputs,
                comment=f"validate external inputs: {k}",
                additional_inputs=lambda k=k: [
                    I(k, 'k'),
                ],
                config_inputs=lambda config: [
                    I(config.Met, as_='config_met'),
                ],
                state_inputs=lambda state, k=k: [
                    I(getattr(state, k), as_='input_vals'),
                ],
            ) for k in keys
        ],
        Process(
            func=validate_PAR_data,
            comment="Validate PAR data",
            state_inputs=lambda state:[
                I(state.PAR, as_='PAR_list'),
                I(state.PPFD, as_='PPFD_list'),
                I(state.Rn, as_='Rn_list'),
            ],
        ),
        Process(
            func=validate_vpd_data,
            comment="Validate VPD data",
            state_inputs=lambda state:[
                I(state.VPD, as_='VPD_list'),
            ],
        ),
    ]


def external_state_init_processes(
    start_row: int = 0, end_row: int = (365 - 1) * 24,
    config_met: Config_Met = None,
) -> List[Process]:
    """Get flattened list of all processes to be ran over a annual cycle."""
    row_indexes = list(range(start_row, end_row))
    # TODO: check dd and hr in data are correct
    row_count = len(row_indexes)
    return flatten_list([
        init_params(row_count),
        convert_units(config_met),
        calc_O3_params_list(config_met, row_indexes),
        calc_thermal_time_processes(config_met),
        calc_CO2_params_list(config_met, row_indexes),
        calc_met_params_list(),
        validate_input(),
    ])

"""These processes are run on the current hours external state.

Processes here should only depend on the external data and config.
"""
from proflow.Objects.Process import Process
from proflow.Objects.Interface import I
from proflow.Switch import switch
from proflow.helpers import set_value

from do3se_met import irradiance as met_irrad_helpers
from thermal_time import calcs as td_functions
from do3se_phenology import phyllochron_dvi

from pyDO3SE.External_State.External_State_Config import ThermalTimeMethods


def accumulate_precipitation_process() -> Process:
    return Process(
        func=lambda precip_acc_dd, precip_mm: precip_acc_dd + (precip_mm / 1000),
        comment="Accumulate precipitation",
        group="soil-moisture",
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


def set_thermal_time_process(
    thermal_time_method: ThermalTimeMethods,
    iLC: int,
) -> Process:
    return Process(
        func=set_value,
        comment="Set thermal time from external data",
        group="thermal-time",
        gate=thermal_time_method == ThermalTimeMethods.EXTERNAL,
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.td[row_index], as_='td'),
        ],
        state_outputs=lambda result: [
            (result['td'], f'canopy_component.{iLC}.td'),
        ],
    )


def accumulate_hourly_temperature_process() -> Process:
    return Process(
        func=lambda t, t_acc, hr: t_acc + t if hr > 0 else t,
        comment="Accumulate daily temperature",
        group="thermal-time",
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.Ts_C[row_index], as_='t'),
        ],
        state_inputs=lambda state: [
            I(state.temporal.hr, as_='hr'),
            I(state.external_met.t_acc, as_='t_acc'),
        ],
        state_outputs=lambda result: [
            (result, 'external_met.t_acc'),
        ],
    )


def calculate_daily_thermal_time_process(
    thermal_time_method: ThermalTimeMethods,
    iLC: int,
) -> Process:
    def _inner(**kwargs):
        initial_args = {**kwargs}
        try:
            return (kwargs.pop("dd", "MISSING_DD") >= kwargs.pop("thermal_time_start", "MISSING thermal_time_start")
                        ) and td_functions.calc_daily_thermal_time(**kwargs) \
                            or kwargs.get("previous_thermal_time")
        except Exception as e:
            print("calculate_daily_thermal_time_process failed with initial_args: ", initial_args)
            raise e

    return Process(
        func=_inner,
        comment="Calculate thermal time from hourly accumulated temperature",
        gate=thermal_time_method == ThermalTimeMethods.HOURLY,
        group="thermal-time",
        config_inputs=lambda config: [
            I(config.Met.thermal_time_start, as_="thermal_time_start"),
        ],
        state_inputs=lambda state: [
            I(state.prev_hour.canopy_component[iLC].td, as_="previous_thermal_time"),
            I(state.temporal.dd, as_="dd"),
            I(state.canopy_component[iLC].t_eff, as_='effective_temperature'),
        ],
        state_outputs=lambda result: [
            (result, f'canopy_component.{iLC}.td'),
        ],
    )


def store_accumulate_precipitation_process() -> Process:
    return Process(
        func=set_value,
        comment="Accumulate precipitation",
        group="soil-moisture",
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


def calc_effective_temperature_process(iLC: int) -> Process:
    """Calculate the effective temperature(t_eff).

    Called at start of day.
    """
    return Process(
        # TODO: Manage hourly here
        func=td_functions.calc_effective_temperature,
        comment="Calculate the effective temperature",
        group="thermal-time",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].pn_gsto.t_b, as_='t_b'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.t_o, as_='t_o'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.t_m, as_='t_m'),
        ],
        state_inputs=lambda state: [
            I(state.external_met.t_acc, as_="t_acc"),
        ],
        state_outputs=lambda result: [
            (result, f'canopy_component.{iLC}.t_eff'),
        ],
    )


def calc_photoperiod_process() -> Process:
    """Calculate the day length (photoperiod)."""
    return Process(
        func=met_irrad_helpers.calc_photoperiod,
        comment="Calculate the day length (photoperiod)",
        config_inputs=lambda config: [
                I(config.Location.lat, as_='lat'),
        ],
        state_inputs=lambda state: [
            I(state.temporal.dd, as_='dd'),
        ],
        state_outputs=lambda result: [
            (result, 'external_met.photoperiod'),
        ],
    )


def calc_photoperiod_factor_process(iLC: int) -> Process:
    """Calculate the photoperiod factor."""
    return Process(
        func=met_irrad_helpers.calc_photoperiod_factor,
        comment="Calculate the photoperiod factor",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].PID, as_='PID'),
        ],
        state_inputs=lambda state: [
            I(state.external_met.photoperiod, as_='photoperiod'),
        ],
        state_outputs=lambda result: [
            (result, f'canopy_component.{iLC}.photoperiod_factor'),
        ],
    )


def calculate_relative_photoperiod_process(iLC: int) -> Process:
    return Process(
        func=phyllochron_dvi.calc_rpe,
        comment="Calculate the relative photoperiod",
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.parameters[iLC].pn_gsto.p_crit,
                as_='p_crit'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.p_sens,
                as_='p_sens'),
        ],
        state_inputs=lambda state: [
            I(state.external_met.photoperiod, as_='p'),
        ],
        state_outputs=lambda result: [
            (result, f'canopy_component.{iLC}.rpe'),
        ],
    )

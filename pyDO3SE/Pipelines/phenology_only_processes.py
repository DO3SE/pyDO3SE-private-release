"""Phenology only DO3SE processes.

This runs the phenology only DO3SE model, which calculates td, dvi, fphen etc

"""

from pathlib import Path
from typing import List
from data_helpers.list_helpers import flatten_list

from proflow.ProcessRunnerCls import advance_time_step_process
from proflow.Objects.Process import Process, ProcessType
from proflow.Objects.Interface import I
from proflow.logger import log_values

from do3se_phenology.config import (
    PhenologyMethods,
    TimeTypes,
)
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.External_State.External_State_Config import (
    ThermalTimeMethods,
)

from pyDO3SE.Pipelines.validation_processes import setup_validation_processes
from .es_hour_processes import (
    accumulate_hourly_temperature_process,
    calc_effective_temperature_process,
    calc_photoperiod_process,
    calc_photoperiod_factor_process,
    calculate_daily_thermal_time_process,
    calculate_relative_photoperiod_process,
    set_thermal_time_process,
)


from .default_processes import (
    lget,
    tag_process,
    calc_DVI_process,
    calc_phyllochron_process,
    calc_if_plant_is_sown_process,
    calc_if_plant_has_emerged_process,
    calc_if_flag_leaf_has_emerged_process,
    calc_emergence_rate_process,
    calc_emerged_leaf_count_process,
    get_phenology_stage_process_td,
    get_phenology_stage_process,
    get_growing_populations_process,
    f_phen_method_process,
    leaf_f_phen_process,
    get_vernalisation_factor_process,
    get_td_dd_process,
    calc_td_dd_per_leaf_pop_process,
    set_hour,
    sync_row_index,
    store_prev_state,
    store_prev_state,
    set_day_offset,
    reset_model,
    set_day,
    save_state,
    daily_start_processes,
    log_external_state_values,
)


def daily_start_processes(config: Config_Shape, run_dir: Path, dump_state_n: int) -> List[Process]:
    """Get processes to be ran at start of day."""
    nL = config.Land_Cover.nL
    nLC = config.Land_Cover.nLC
    nP = config.Land_Cover.nP
    thermal_time_method = config.Met.thermal_time_method
    multi_season = config.Location.multi_season
    time_type = config.Land_Cover.phenology_options.time_type
    phenology_method = config.Land_Cover.phenology_options.phenology_method
    f_phen_method = [
        config.Land_Cover.parameters[iLC].phenology.f_phen_method for iLC in range(nLC)]
    leaf_f_phen_method = [
        config.Land_Cover.parameters[iLC].phenology.leaf_f_phen_method for iLC in range(nLC)]
    dvi_method = config.Land_Cover.dvi_method
    plant_emerge_method = config.Land_Cover.phenology_options.plant_emerge_method
    flag_leaf_emerge_method = config.Land_Cover.phenology_options.flag_leaf_emerge_method
    use_vernalisation = config.Land_Cover.phenology_options.use_vernalisation
    use_photoperiod_factor = config.Land_Cover.phenology_options.use_photoperiod_factor
    sparse_data = config.Met.sparse_data

    return [
        tag_process("===== Start of Daily Processes ====="),
        set_day_offset(multi_season),
        reset_model() if multi_season else [],
        set_day(),
        save_state(run_dir, dump_state_n)if dump_state_n is not None else [],

        # MET
        tag_process("== Met Processes =="),
        # perNlc(nLC, calc_effective_temperature_process),
        [calc_effective_temperature_process(iLC) for iLC in range(nLC)] if not sparse_data else [],
        [calculate_daily_thermal_time_process(
            thermal_time_method, iLC) for iLC in range(nLC)] \
        if thermal_time_method != ThermalTimeMethods.EXTERNAL and not sparse_data else [],
        calc_photoperiod_process(),
        [calc_photoperiod_factor_process(iLC)
         for iLC in range(nLC)]if use_photoperiod_factor else [],

        # PHENOLOGY
        tag_process("== Phenology Processes =="),
        # f_phen
        [[
            f_phen_method_process(iL, iLC, f_phen_method[iLC]),
            leaf_f_phen_process(iL, iLC, leaf_f_phen_method[iLC]),
        ] for iL in range(nL) for iLC in range(nLC)] if phenology_method != PhenologyMethods.DISABLED else [],


        [
            [
                calc_DVI_process(iLC, dvi_method),
                calc_phyllochron_process(iLC),
                calculate_relative_photoperiod_process(iLC),
                calc_if_plant_is_sown_process(iLC),
                calc_if_plant_has_emerged_process(iLC, plant_emerge_method),
                get_vernalisation_factor_process(iLC) if use_vernalisation else [],
                get_td_dd_process(iLC) if time_type == TimeTypes.THERMAL_TIME else [],
                calc_if_flag_leaf_has_emerged_process(iLC, flag_leaf_emerge_method),
                calc_emergence_rate_process(nP, iLC),
                calc_emerged_leaf_count_process(nP, iLC, time_type),
                calc_td_dd_per_leaf_pop_process(iLC) if time_type == TimeTypes.THERMAL_TIME else [],
                get_phenology_stage_process_td(
                    nP, iLC) if time_type == TimeTypes.THERMAL_TIME else [],
                get_phenology_stage_process(
                    nP, iLC, leaf_f_phen_method[0]) if time_type == TimeTypes.JULIAN_DAY else [],
                get_growing_populations_process(nP, iLC),
            ]
            for iLC in range(nLC)
        ] if phenology_method != PhenologyMethods.DISABLED else [],

    ]


def log_processes(nL: int, nLC: int, nP: int, fields: List[str]):
    """Log state to output data.

    Where a single layer is used should output top layer
    Where a single population is used should output flag leaf

    """
    flag_index = nP - 1
    top_layer_index = nL - 1
    component_index = 0
    def __SPACER__(name): return I(f'|{name}|', as_='|')
    # TODO: Can we optimize the below logging?
    return [
        # TODO: This should come from config or output_values_map
        log_external_state_values(lambda e_state, row_index: list(filter(lambda f: f, [
            __SPACER__("External Values >"),
            I(row_index, as_='row_index') if 'row_index' in fields else None,
            I(lget(e_state.dd, row_index), as_='dd_e') if 'dd_e' in fields else None,
            I(lget(e_state.hr, row_index), as_='hr') if 'hr' in fields else None,
            I(e_state.Ts_C[row_index], as_='ts_c') if 'ts_c' in fields else None,
            I(lget(e_state.PAR, row_index), as_='par') if 'par' in fields else None,
            I(lget(e_state.P, row_index), as_='p') if 'p' in fields else None,
            I(lget(e_state.u, row_index), as_='uh_zr') if 'uh_zr' in fields else None,
            I(lget(e_state.precip, row_index), as_='precip') if 'precip' in fields else None,
            I(lget(e_state.O3, row_index), as_='o3_ppb_zr') if 'o3_ppb_zr' in fields else None,
            I(lget(e_state.RH, row_index), as_='rh') if 'rh' in fields else None,
            I(lget(e_state.CO2, row_index), as_='co2') if 'co2' in fields else None,
            I(lget(e_state.Hd, row_index), as_='hd') if 'hd' in fields else None,
            I(lget(e_state.R, row_index), as_='r') if 'r' in fields else None,
            I(lget(e_state.Rn, row_index), as_='rn') if 'rn' in fields else None,
            I(lget(e_state.VPD, row_index), as_='vpd') if 'vpd' in fields else None,
            I(lget(e_state.SWC, row_index), as_='swc') if 'swc' in fields else None,
        ]))),
        log_values(lambda state: list(filter(lambda f: f, [
            # Values here should match values in Output_Shape
            __SPACER__("Canopy Level >"),
            I(state.temporal.dd, as_='dd') if 'dd' in fields else None,
            I(state.temporal.dd_offset, as_='dd_offset') if 'dd_offset' in fields else None,
            I(state.canopy_component[component_index].td, as_='td') if 'td' in fields else None,
            I(state.canopy.Vd, as_='canopy_vd') if 'canopy_vd' in fields else None,
            I(state.canopy_component[component_index].dvi, as_='dvi') if 'dvi' in fields else None,
            I(state.canopy_component[component_index].td_dd,
              as_='td_dd') if 'td_dd' in fields else None,
            I(state.canopy_component[component_index].t_eff,
              as_='t_eff') if 't_eff' in fields else None,
            I(state.canopy_component[component_index].td_v,
              as_='td_v') if 'td_v' in fields else None,
            I(state.canopy_component[component_index].phenology.phenology_stage,
              as_='phenology_stage') if 'phenology_stage' in fields else None,
            I(state.canopy_component[0].growing_populations[flag_index],
              as_=f"growing_populations_flag") if 'growing_populations' in fields else [],
            I(state.canopy_component[component_index].phenology.Vf,
              as_='Vf') if 'Vf' in fields else None,
            I(state.canopy_component[component_index].phenology.V_acc,
              as_='V_acc') if 'V_acc' in fields else None,
            I(state.canopy_component[component_index].phenology.V_tot,
              as_='V_tot') if 'V_tot' in fields else None,
            I(state.canopy_component[component_index].phenology.V_pos,
                as_='V_pos') if 'V_pos' in fields else None,
            I(state.canopy_component[component_index].phenology.V_neg,
                as_='V_neg') if 'V_neg' in fields else None,
            I(state.external_met.photoperiod, as_='photoperiod') if 'photoperiod' in fields else None,
            I(state.canopy_component[component_index].photoperiod_factor,
              as_='photoperiod_factor') if 'photoperiod_factor' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.f_phen, as_='f_phen') if 'f_phen' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.leaf_f_phen, as_='leaf_f_phen') if 'leaf_f_phen' in fields else None,
        ]))),
    ]


def hourly_processes(config: Config_Shape, hr: int, run_dir: Path = None, dump_state_n: int = None) -> List[Process]:
    """Take the current hour and returns a list of processes.

    to be ran in this hour
    """
    nL = config.Land_Cover.nL
    nLC = config.Land_Cover.nLC
    nP = config.Land_Cover.nP

    thermal_time_method = config.Met.thermal_time_method

    sparse_data = config.Met.sparse_data
    output_fields = config.output.fields
    log_multilayer = config.output.log_multilayer

    return [
        tag_process(f"===== Start of Hourly Processes (Hour: {hr}) ====="),
        set_hour(hr),
        # dump_state_process(), ## Turn on to debug state at specific hours
        sync_row_index(),  # External data row pointer moved to hr = 0

        [set_thermal_time_process(thermal_time_method, iLC) for iLC in range(nLC)]
        if thermal_time_method == ThermalTimeMethods.EXTERNAL else [],

        # TODO: It would be great to not rely on knowing the hour here.
        daily_start_processes(config, run_dir, dump_state_n) if (
            hr == 0 or sparse_data) else [],


        accumulate_hourly_temperature_process(),
        # LOGGING
        tag_process("== Hourly Logging Processes =="),
        log_processes(nL, nLC, nP, output_fields),
        store_prev_state(),  # Used when we require the data from previous hour.
        # At end of hour advance time step
        advance_time_step_process(),
        tag_process("===== End of Hourly Processes ====="),
    ]


def get_row_processes_hourly(
    config: Config_Shape,
    hours: List[int],
    run_dir: Path = None,
    dump_state_n: int = None,
) -> List[Process]:
    """Get the model processes for each hour in hours.

    Parameters
    ----------
    config : Config_Shape
        _description_
    hours : List[int]
        _description_
    run_dir : Path, optional
        If this and dump_state_n provided will
        periodically dump state to file, by default None
    dump_state_n : int, optional
        _description_, by default None

    Returns
    -------
    List[Process]
        _description_
    """
    return flatten_list([hourly_processes(config, hr, run_dir, dump_state_n) for hr in hours])


def full_model_processes(
    config: Config_Shape,
    hours: List[int],
    run_validation: bool = False,
    run_dir: Path = None,
    dump_state_n: int = None,
) -> List[Process]:
    """Get a flattened list of all processes to be ran between dates inclusive."""
    row_processes = get_row_processes_hourly(config, hours, run_dir, dump_state_n)
    all_processes = [
        setup_validation_processes(config) if run_validation else [],
        row_processes,
    ]
    return flatten_list(all_processes)

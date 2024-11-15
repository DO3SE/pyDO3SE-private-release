"""Model setup function.

Returns the arguments for the run model function.
"""
import os
import re
import math
import pandas as pd
from copy import deepcopy
from collections import namedtuple
from pathlib import Path
from typing import Callable, List, Tuple, Dict
from datetime import datetime
from enum import Enum
from data_helpers.cls_parsing import rsetattr
from data_helpers.list_helpers import flatten_list
from proflow.ProcessRunnerCls import ProcessRunner
from proflow.Objects.Process import Process

from pyDO3SE.External_State.setup_dd import get_day_crop_override
from pyDO3SE.util.logger import Logger
from pyDO3SE.version import config_version
from pyDO3SE.Config import Config_Shape
from pyDO3SE.External_State import External_State_Shape
from pyDO3SE.External_State.external_state_loader import (
    run_init_processes_on_e_state
)
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.Pipelines.default_processes import (
    full_model_processes,
    hourly_processes,
)
from pyDO3SE.Pipelines.validation_processes import setup_validation_processes
from pyDO3SE.Pipelines.state_init_processes import state_init_processes
from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes
from pyDO3SE.Pipelines.config_init_processes import config_init_processes
from pyDO3SE.Output.Output_Shape import default_output_fields
from pyDO3SE.overrides import Main_Overrides

Model = namedtuple('Model', 'config,external_state,initial_state,model_processes')


class LocationMethod(Enum):
    CONFIG_INPUT = "config_input"
    EXTERNAL_STATE_INPUT = "external_state_input"


def setup_config(
    config_in: Config_Shape,
    external_state_data: External_State_Shape = None,
    config_overrides: Dict[str, any] = {},
    logger: Logger = Logger(),
    overrides: Main_Overrides = Main_Overrides(),
) -> Config_Shape:
    """Setup the model config from input config and data.

    Parameters
    ----------
    config_in: Config_Shape,
        loaded raw config
    external_state_data : Path
        External state data.
    config_overrides: Dict[str, any]
        a mapping of config field using dot notation against a value override
        e.g. {"Location.lat": 99}
    overrides : Main_Overrides
        model overrides

    Returns
    -------
    Config
        initialized config

    """
    if overrides.config_override:
        return overrides.config_override
    if config_in.VERSION != config_version:
        raise ValueError(
            f"Config must be updated to latest version {config_version}. Config is version {config_in.VERSION}. Use pyDO3SE config migrate [CONFIG_FILE]. Config description: {config_in.id}")

    process_runner = ProcessRunner(DEBUG_MODE=overrides.debug)
    process_runner.external_state = external_state_data
    process_runner.config = config_in
    init_processes = overrides.init_config_processes or config_init_processes(config_in)

    config_amended = deepcopy(config_in)
    if config_overrides:
        logger(f"Overriding the config with the following values: {config_overrides}")
        for k, v in config_overrides.items():
            if v is not None and not (type(v) == float and math.isnan(v)):
                config_amended = rsetattr(config_amended, k, v)
    config_amended: Config_Shape = process_runner.run_processes(
        init_processes,
        config_amended,
    )
    if config_amended.VERSION != config_version:
        raise ValueError(
            f"Config must be updated to latest version {config_version}. Config is version {config_amended.VERSION}. Use pyDO3SE config migrate [CONFIG_FILE]. Config description: {config_amended.id}")
    if overrides.output_fields is not None:
        config_amended.output.fields = overrides.output_fields
    if not config_amended.output.fields:
        config_amended.output.fields = default_output_fields
    return config_amended


def pad_external_data(external_data, start_day, end_day) -> External_State_Shape:
    raise NotImplementedError()
    # MAX_LENGTH = 600  # Maximum span for data from 0
    # data_start_day = external_data.dd[0] - 1
    # data_end_day = external_data.dd[-1]

    # keys = list(External_State_Shape.__annotations__.keys())
    # data_out = {}
    # for k in keys:
    #     data = list(0 for _ in range(data_start_day)) + \
    #         getattr(external_data, k) + list(0 for _ in range(MAX_LENGTH))
    #     data_out[k] = data[start_day:end_day]
    # return External_State_Shape(**data_out)


def setup_external_state_simple(
    external_state_data: External_State_Shape,
    config: Config_Shape,
    init_processes: List[Process],
    overrides: Main_Overrides = Main_Overrides(),
) -> External_State_Shape:
    return run_init_processes_on_e_state(
        external_state_data,
        config,
        init_processes,
        overrides=overrides,
    )


def setup_external_state(
    config: Config_Shape,
    external_state_data: External_State_Shape,
    pad_data: bool = False,
    logger: Callable[[str], None] = Logger(),
    overrides: Main_Overrides = Main_Overrides(),
) -> Tuple[External_State_Shape, int, int]:
    """Setup the external state from a data location and config.

    NOTE: Mutates external state data input

    Understanding day ranges in the DO3SE model
    ===========================================

    We can override the start and end date either in the config or overrides.
    Input data typically uses dd==1 for row==0. The below setup takes this into account.

    Example (Overrides)
    -------------------
    - Overrides: `start_day = 3, end_day = 10`
    - External state `dd = [1,1,1,1...15,15,15,15]`

    The resulting sliced dd data would be `[3,3,3,3....10,10,10,10,10]`.
    This would be 8 days of data.
    Note that the start and end dates are inclusive so we get `end_day - start_day + 1 = 8 days`.

    Example (No Overrides)
    ----------------------
    - Overrides: `start_day = None, end_day = None`
    - External state `dd = [1,1,1,1...15,15,15,15]`

    The resulting start date and end are `start_day=1, end_day=15` with 15 days of data.


    Parameters
    ----------
    config : Config_Shape
        Config input
    external_state_data : External_State_Shape
        Raw external data
    pad_data: bool, optional
        If true then we pad data to full days if it does not start or end at midnight
    logger: Logger
        logger
    overrides : Main_Overrides, optional
        Overrides that have been passed to main, by default Main_Overrides()


    Returns
    -------
    External_State_Shape
        Processed external data
    start_date: int
        #
    end_date: int
        #

    Raises
    ------
    InputDataError
        Invalid input data supplied
    DayRangeError
        Invalid day range supplied

    """
    if overrides.external_state_override:
        return overrides.external_state_override

    start_day, end_day = get_day_crop_override(
        config.Location.start_day,
        config.Location.end_day,
        overrides,
    )

    if start_day is not None or end_day is not None:
        logger(f"Overriding date range with start_day: {start_day} and end_day: {end_day}")

    # Use the process runner to modify the external state data based on some
    # initialization processes
    process_runner = ProcessRunner(config, DEBUG_MODE=overrides.debug)
    process_runner.external_state = external_state_data

    # TODO: replace external_state_init_processes with get_external_...
    init_processes = overrides.init_external_state_processes or \
        external_state_init_processes(start_day, end_day, pad_data, process_runner.config.Met)
    external_state = process_runner.run_processes(
        init_processes,
        external_state_data,
    )
    assert len(external_state.dd) == len(
        external_state.hr), f"Day and hour lengths do not match. Check data ends at hour 23.\n dd: {len(external_state.dd)} hr: {len(external_state.hr)}"

    start_day = external_state.dd[0]
    end_day = external_state.dd[-1]

    logger(f"External state setup with start_day: {start_day}, end_day: {end_day}")

    return external_state, start_day, end_day


def run_first_hour_on_state(
    state_in: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    debug: bool = False,
) -> Model_State_Shape:
    process_runner = ProcessRunner(config, external_state, DEBUG_MODE=debug)
    init_processes = flatten_list(hourly_processes(config, 0))
    assert external_state is not None, "External state must be supplied to run first hour"
    state = process_runner.run_processes(
        init_processes,
        state_in,
    )
    return state


def setup_initial_state(
    state_in: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    run_init_processes: bool = False,
    override_init_state_processes: List[Process] = None,
    logger: Callable[[str], None] = Logger(),
    overrides: Main_Overrides = Main_Overrides(),
) -> Model_State_Shape:
    try:
        state = deepcopy(state_in)
        if state_in.prev_hour == None:
            # TODO: Check this still works
            state.prev_hour = deepcopy(state_in)
        if run_init_processes:
            logger("Running state init processes")
            process_runner = ProcessRunner(config, external_state, DEBUG_MODE=overrides.debug)
            init_processes = override_init_state_processes or \
                state_init_processes(config)
            state = process_runner.run_processes(
                init_processes,
                state,
            )
        return state
    except Exception as e:
        raise Exception("Failed to run setup initial state") from e


def setup_model_processes(
    config: Config_Shape,
    overrides: Main_Overrides,
    hours: List[int] = None,
    run_validation: bool = False,
    run_dir: Path = None,
    logger: Callable[[str], None] = Logger(),
) -> List[Process]:
    """Get a list of all processes to be run by the model.

    These processes are then run by the process runner.

    Parameters
    ----------
    config : Config_Shape
        model input config
    overrides : Main_Overrides
        Model overrides
    hours: List[int]
        Hours to run
    run_validation: bool
        if true prepends validation processes
    run_dir : Path
        Path to current project run
    logger: Logger
        Log class

    Returns
    -------
    List[Process]
        A list of processes to be ran by the model.

    """
    if overrides.model_processes:
        return overrides.model_processes
    logger("Getting model processes")

    return full_model_processes(config, hours, run_validation, run_dir, overrides.dump_state_n)


def setup_model_processes_hour_map(
    config: Config_Shape,
    overrides: Main_Overrides = Main_Overrides(),
    run_validation: bool = False,
    run_dir: Path = None,
    logger: Callable[[str], None] = Logger(),
) -> List[Process]:
    """Get a list of all processes to be run by the model mapped to hour index.

    I.e. {
        0: List[Process],
        1: List[Process],
        ...
    }

    These processes are then run by the process runner.

    Parameters
    ----------
    config : Config_Shape
        model input config
    overrides : Main_Overrides
        Model overrides
    hours: List[int]
        Hours to run
    run_validation: bool
        if true prepends validation processes
    run_dir : Path
        Path to current project run. Only required if logging state for debugging
    logger: Logger
        Log class

    Returns
    -------
    List[Process]
        A list of processes to be ran by the model.

    """
    if overrides.model_processes:
        return overrides.model_processes
    logger("Getting model processes")
    return {
        hr: full_model_processes(config, [hr], run_validation, run_dir, overrides.dump_state_n) for hr in range(24)
    }


def run_model_validation(
    config: Config_Shape,
    initial_state: Model_State_Shape,
    external_state: External_State_Shape,
):
    process_runner = ProcessRunner(config, external_state)
    model_state = deepcopy(initial_state)
    process_runner.run_processes(
        setup_validation_processes(config),
        model_state,
    )


def setup_model(
    config_in: Config_Shape,
    state_in: Model_State_Shape,
    external_state_in: External_State_Shape,
    # data_location: Path=None,
    run_dir: Path = None,
    config_overrides: Dict[str, any] = None,
    logger: Callable[[str], None] = Logger(),
    overrides: Main_Overrides = Main_Overrides(),
) -> Model:
    """Get the model arguments from the input locations.

    Parameters
    ----------
    config_in: Config_Shape
        Merged and loaded config file
    state_in: Model_State_Shape
        Initial loaded state
    external_state_in: External_State_Shape
        Loaded external data
    run_dir: Path
        Run directory for outputs
    config_overrides: Dict[str, any]
        a mapping of config field using dot notation against a value override
        e.g. {"Location.lat": 99}
    logger: Callable[[str], None], optional
        Log function, by default print
    overrides : Main_Overrides
        model overrides

    Returns
    -------
    (config,
    external_state,
    initial_state,
    model_processes)

    """
    logger("Setting up model")
    start_time = datetime.now()
    # external_state_data = next(load_external_state(
    #     data_location,
    #     logger=logger,
    # ))

    config = setup_config(
        config_in,
        external_state_in,
        config_overrides,
        logger=logger,
        overrides=overrides,
    )
    external_state, start_day, end_day = setup_external_state(
        config, external_state_in, logger=logger, overrides=overrides)
    config.Location.start_day = start_day
    config.Location.end_day = end_day

    assert (
        start_day is None or end_day is None) \
        or start_day <= end_day, f"Start day({start_day}) is after end day({end_day})!"

    overrides = overrides._replace(start_day=start_day, end_day=end_day)
    initial_state = setup_initial_state(
        state_in,
        config,
        external_state,
        run_init_processes=not overrides.skip_state_init,
        override_init_state_processes=overrides.init_state_processes,
        logger=logger,
        overrides=overrides,
    )

    if overrides.run_validation:
        run_model_validation(config, initial_state, external_state)

    # TODO: Implement both methods
    # model_processes = setup_model_processes(config, overrides, hours, logger=logger)
    model_processes = setup_model_processes_hour_map(
        config, overrides=overrides, run_dir=run_dir, logger=logger)
    time_taken = datetime.now() - start_time
    logger(f"Model setup complete in {time_taken}")

    return Model(
        config=config,
        external_state=external_state,
        initial_state=initial_state,
        model_processes=model_processes,
    )


def get_input_files_list(
    dir: Path,
    multi_file_netcdf: bool = False,
    regex_multi_file_filter: str = None,
    start_input_index: int = None,
    end_input_index: int = None,
    logger: Callable[[str], None] = print,
) -> List[str]:
    """Gets a sorted list of files to run.

    When using multi file netcdfs we need to pull out the file groups that are to be ran.
    E.g.

        files_list = [
            'demo_wrf_2014-11_HFX_FORCE',
            'demo_wrf_2014-11_td_2m',
            'demo_wrf_2014-12_rh',
            'demo_wrf_2014-11_SWDOWN',
            'demo_wrf_2014-11_RAINNC',
            'demo_wrf_2014-12_SNOWH',
            'demo_wrf_2014-12_RAINNC',
            'demo_wrf_2014-12_td_2m',
            'demo_wrf_2014-11_pres',
            'demo_wrf_2014-11_SNOWH',
            'demo_wrf_2014-12_wspeed',
            'demo_wrf_2014-12_pres',
            'demo_wrf_2014-11_o3',
            'demo_wrf_2014-12_SWDOWN',
            'demo_wrf_2014-11_wspeed',
            'demo_wrf_2014-12_HFX_FORCE',
            'demo_wrf_2014-11_rh',
            'demo_wrf_2014-12_o3',
        ]

    We need to process the above into:
    files_list_groups = [
        2014-11,
        2014-21,
    ]


    Parameters
    ----------
    dir : Path
        directory to load files from
    multi_file_netcdf: boolean
        If true then we are loading a multi file netcdf
    regex_multi_file_filter : str, optional
        regex pattern to group files for loading, e.g. '[0-9]{4}-[0-9]{2}'
    start_input_index: int
        If set then start at this input index
    end_input_index: int
        If set then end at this input index
    logger : Callable[[str], None], optional
        logger func, by default print

    Returns
    -------
    List[str]
        List of files or file group filters to run.

    """
    files_list = sorted(os.listdir(dir))
    if regex_multi_file_filter is None:
        if multi_file_netcdf:
            files_list = ["*"]
        return files_list
    else:
        files_list = sorted(set([
            getattr(re.search(regex_multi_file_filter, f), 'group', lambda: None)()
            for f in files_list
        ]))

    if start_input_index is None:
        start_input_index = 0
    if end_input_index is None:
        end_input_index = len(files_list) - 1
    files_list = files_list[start_input_index:end_input_index + 1]
    total_files = len(files_list)
    logger(
        f"== Running Model: Total File: {total_files}. Start index: {start_input_index}. End index: {end_input_index} ===")
    return files_list


def get_config_overrides(
    input_id: str,
    per_input_config_overrides: pd.DataFrame,
) -> Dict[str, any]:
    if per_input_config_overrides is None:
        return {}
    try:
        return per_input_config_overrides[per_input_config_overrides.inputid == input_id].iloc[0].to_dict()
    except IndexError as e:
        print(per_input_config_overrides.head())
        raise IndexError(
            f"input id '{input_id}' could not be found in the per_input_config_overrides!")

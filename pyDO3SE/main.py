r"""The pyDO3SE main entrpoint.

This should usually be ran via the :mod:`CLI <pyDO3SE.tools.cli>`, juptyer notebooks or from other libraries.

This file should contain various model entrypoints that link to the UI or CLI.

For information on how to run the model view :mod:`pyDO3SE Home <pyDO3SE>`

-------------------------------------------------------------------------------
Data flow in model
==================

.. code-block:: python

    {config file(json)}                            {external_data(csv)}
                |                                               |
                v                                               v
        (config_loader)                             (load_external_data)
                |                                               |
                v                                               |
            (setup_config)                                      |
                |                                               |
                v                                               |
        <Config_Shape>                                          |
            |   | |                                             |
            |   | |________________________________________     |
            |   |                                          |    |
            |   |                                          v    v
            |   |                                    (setup_external_state)
            |   |                                          |    |
            |   | ________________                         v    v
            |                     |                    <External_State_Shape>
            |                     |    ____________________|    |
            |                     |   |                         |
            |                     v   v                         |
            |               (final_config)                      |
            |                      |                            |
            |                      |___________                 |
            |                      |           |                |
            |                      |           |                |
            |                      |           |                |
            |                      |           |_________       |
            |______________________|____________________ | _____|
                                \|/                     \|/
                                 v                       |
                            (setup_state)                |
                                |                        |
                                v               (get_model_processes)
                <Model_State_Shape(Initial)>             |
                                |                        v
                                |                  <ProcessRunner>
                                |                        |
                                v                        v
                        <Model_State_Shape>----->(Run_Process)
                            |  ^                  |
                            |  |                  |
                            |  |   Loops through  |
                            |  |   processes      |
                            |  |__________________|
                            |
                            V
                        Model Outputs


"""

import sys
import os
import re
import shutil
import warnings
import numpy as np
import pandas as pd
import deprecated
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, List, NamedTuple, Tuple, Union
from datetime import datetime

from pyDO3SE.Config.config_loader import config_loader, grid_config_loader
from pyDO3SE.Model_State.model_state_loader import dump_state_to_file, dump_state_to_file_quick, load_current_cell_state, model_state_loader, model_state_loader_quick
from pyDO3SE.Pipelines.default_processes import get_row_processes_hourly
from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes
from pyDO3SE.util.logger import Logger, generate_run_notes
from pyDO3SE.version import config_version, version as model_version
from pyDO3SE.External_State.external_state_loader import (
    get_date_bounds_from_ext_data,
    load_external_state,
)
from pyDO3SE.External_State.External_State_Config import (
    FileTypes,
    EStateOptions,
)
from pyDO3SE.util.loader import csv_loader, json_loader
from pyDO3SE.setup_model import (
    Main_Overrides,
    get_grid_coords_from_file,
    get_input_files_list,
    initialize_grid_configs,
    setup_config,
    setup_external_state_simple,
    setup_initial_state,
    setup_model,
)
from pyDO3SE.run_model import run_model, run_model_daily
from pyDO3SE.Output.process_outputs import (
    dump_output_to_file_netcdf_grid,
    export_output,
)

from pyDO3SE.Config import Config_Shape
from pyDO3SE.Model_State import Model_State_Shape

# Currently we can only import json configs
CONFIG_FILE_TYPE = "json"
DATA_FILE_TYPE = "csv"


class GridProjectPaths(NamedTuple):
    project_dir: str
    config_dir: str
    input_data_dir: str
    e_state_overrides_file_path: str
    variable_map_path: str
    preprocess_map_path: str
    e_state_overrides_field_map_path: str
    base_config_path: str
    base_state_path: str


class GridRunPaths(NamedTuple):
    run_dir: str
    log_path: str
    config_path: str
    e_state_overrides_field_map_path: str
    initial_state_dir: str
    output_data_dir: str
    processed_configs_dir: str
    prev_state_dir: str
    live_state_dir: str
    final_state_dir: str
    run_mask_path: str


class GridRunFiles(NamedTuple):
    e_state_overrides_field_map: dict
    config: Config_Shape
    state: Model_State_Shape
    preprocess_map: dict
    variable_map: dict


class ProjectPaths(NamedTuple):
    project_dir: str = None
    config_dir: str = None
    input_data_dir: str = None
    output_directory: str = None
    variable_map_path: str = None
    preprocess_map_path: str = None
    base_config_path: str = None
    base_state_path: str = None


class RunPaths(NamedTuple):
    run_id: str = None
    run_dir: str = None
    log_path: str = None
    config_path: str = None
    output_directory: str = None
    input_data_file_path: str = None


class RunFiles(NamedTuple):
    config: Config_Shape
    state: Model_State_Shape


def get_configs(dir: Path):
    return os.listdir(dir)


def get_grid_project_paths(
    project_dir: str,
) -> GridProjectPaths:
    return GridProjectPaths(
        project_dir=project_dir,
        config_dir=f"{project_dir}/configs",
        input_data_dir=f"{project_dir}/inputs",
        e_state_overrides_file_path=f"{project_dir}/e_state_overrides.nc",
        variable_map_path=f"{project_dir}/variable_map.json",
        preprocess_map_path=f"{project_dir}/preprocess_map.json",
        e_state_overrides_field_map_path=f"{project_dir}/e_state_overrides_field_maps",
        base_config_path=f"{project_dir}/base_config.json",
        base_state_path=f"{project_dir}/base_state.json",
    )


def get_grid_run_paths(
    project_paths: GridProjectPaths,
    config_id: str,
    run_id: str,
) -> GridRunPaths:
    run_dir = f"{project_paths.project_dir}/runs/{run_id}/{config_id}"

    return GridRunPaths(
        run_dir=run_dir,
        log_path=f"{run_dir}/run.log",
        config_path=f"{project_paths.project_dir}/configs/{config_id}.json",
        e_state_overrides_field_map_path=f"{project_paths.e_state_overrides_field_map_path}/{config_id}.json",
        initial_state_dir=f"{run_dir}/initial_state",
        live_state_dir=f"{run_dir}/current_state",
        final_state_dir=f"{run_dir}/final_state",
        prev_state_dir=f"{run_dir}/prev_state",
        output_data_dir=f"{run_dir}/outputs_grid",
        processed_configs_dir=f"{run_dir}/processed_configs",
        run_mask_path=f"{project_paths.project_dir}/run_masks/{config_id}.nc",
    )


def get_project_paths(
    project_dir: str,
) -> ProjectPaths:
    return ProjectPaths(
        project_dir=project_dir,
        config_dir=f"{project_dir}/configs",
        input_data_dir=f"{project_dir}/inputs",
        output_directory=f"{project_dir}/outputs",
        variable_map_path=f"{project_dir}/variable_map.json",
        preprocess_map_path=f"{project_dir}/preprocess_map.json",
        base_config_path=f"{project_dir}/base_config.json",
        base_state_path=f"{project_dir}/base_state.json",
    )


def get_run_paths(
    project_paths: ProjectPaths,
    config_id: str,
    input_file_id: str,
    run_id: str,
) -> RunPaths:
    run_dir = f"{project_paths.output_directory}/{run_id}/{config_id}/{input_file_id}"

    return RunPaths(
        run_id=run_id,
        run_dir=run_dir,
        log_path=f"{run_dir}/run.log",
        config_path=f"{project_paths.project_dir}/configs/{config_id}.json",
        input_data_file_path=f"{project_paths.project_dir}/inputs/{input_file_id}.csv",
        output_directory=f"{run_dir}",
    )


def create_grid_run_path_directories(
    run_paths: GridRunPaths,
):
    os.makedirs(run_paths.run_dir, exist_ok=True) if run_paths.run_dir else None
    os.makedirs(run_paths.processed_configs_dir,
                exist_ok=True) if run_paths.processed_configs_dir else None
    os.makedirs(run_paths.prev_state_dir, exist_ok=True) if run_paths.prev_state_dir else None
    os.makedirs(run_paths.live_state_dir, exist_ok=True) if run_paths.live_state_dir else None
    os.makedirs(run_paths.output_data_dir, exist_ok=True) if run_paths.output_data_dir else None
    os.makedirs(run_paths.initial_state_dir, exist_ok=True) if run_paths.initial_state_dir else None
    os.makedirs(run_paths.final_state_dir, exist_ok=True) if run_paths.final_state_dir else None


def create_run_path_directories(
    run_paths: RunPaths,
):
    os.makedirs(run_paths.run_dir, exist_ok=True) if run_paths.run_dir else None
    os.makedirs(run_paths.output_directory, exist_ok=True) if run_paths.output_directory else None


# def get_grid_coords_from_file(
#     run_mask_path: str,
# ) -> Tuple[np.ndarray, int, int]:
#     """Get the grid coordinates and size from a coordinates file

#     Parameters
#     ----------
#     run_mask_path : str
#         Path to a coordinates csv that has headings for x and y

#     Returns
#     -------
#     Tuple[np.ndarray, int, int]
#         grid_coords, grid_x_size, grid_y_size

#     """
#     grid_coords = np.array([[int(i['x']), int(i['y'])] for i in csv_loader(run_mask_path)])
#     grid_x_size, grid_y_size = np.ptp(grid_coords, axis=0) + [1, 1]

#     return grid_coords, grid_x_size, grid_y_size


def init_grid_model_config(
    project_paths: GridProjectPaths,
    run_paths: GridRunPaths,
    loaded_files: GridRunFiles,
    grid_coords,
    logger,
    netcdf_chunks: dict = None,
):
    logger(f'=== Initializing for: {run_paths.run_dir} ===')
    start_time = datetime.now()

    main_partial_initialize(
        config=loaded_files.config,
        state=loaded_files.state,
        processed_config_dir=run_paths.processed_configs_dir,
        state_out_path=run_paths.live_state_dir,
        e_state_overrides_file_path=project_paths.e_state_overrides_file_path,
        e_state_overrides_field_map=loaded_files.e_state_overrides_field_map,
        grid_coords=grid_coords,
        logger=logger,
        # netcdf_chunks=netcdf_chunks,
    )
    logger(f'=== Initialization complete: {run_paths.run_dir} ===')
    end_time = datetime.now()
    setup_duration = end_time - start_time
    return setup_duration


def init_all_grid_model_configs(
    project_paths: GridProjectPaths,
    runid: str = '',
    log_level: int = 0,
    netcdf_chunks: dict = None,
):
    """Initialize grid state for all configs in project dir

    Parameters
    ----------
    project_paths: GridProjectPaths
        file paths specific to project
    runid : str, optional
        Unique id for run, by default ''
    log_level : int, optional
        Log level, by default 0
    netcdf_chunks: dict, optional
        Chunks to use when loading netcdf data

    """
    logger_main = Logger(log_level)
    configs = os.listdir(project_paths.config_dir)
    logger_main(f'== Found {len(configs)} configs to run =====')
    for config_file_path in configs:
        config_name = '.'.join(config_file_path.split('.')[:-1])
        run_paths = get_grid_run_paths(project_paths, config_name, runid)
        create_grid_run_path_directories(run_paths)
        loaded_files: GridRunFiles = load_grid_run_files(project_paths, run_paths)
        logger = Logger(log_level, run_paths.log_path, write_mode='w', set_as_default=True)
        logger(f'== Running config {config_file_path} ==')

        grid_coords, _, _= get_grid_coords_from_file(run_paths.run_mask_path)

        setup_duration = init_grid_model_config(
            project_paths,
            run_paths,
            loaded_files,
            grid_coords,
            logger,
            # netcdf_chunks,
        )
        logger.close()
        return setup_duration


def load_run_files(
    project_paths: GridProjectPaths,
    run_paths: GridRunPaths,
) -> RunFiles:
    return RunFiles(
        config=config_loader(run_paths.config_path, project_paths.base_config_path, 'json'),
        state=project_paths.base_state_path and model_state_loader(
            project_paths.base_state_path, None, 'json', False) or Model_State_Shape(),
    )


def load_grid_run_files(
    project_paths: GridProjectPaths,
    run_paths: GridRunPaths,
) -> GridRunFiles:
    return GridRunFiles(
        e_state_overrides_field_map=json_loader(run_paths.e_state_overrides_field_map_path),
        config=config_loader(run_paths.config_path, project_paths.base_config_path, 'json'),
        state=model_state_loader(project_paths.base_state_path, None, 'json'),
        preprocess_map=json_loader(project_paths.preprocess_map_path),
        variable_map=json_loader(project_paths.variable_map_path),
    )


def main(
    project_paths: ProjectPaths,
    run_paths: RunPaths,
    runnotes: str = '',
    output_results_only: bool = False,
    fields_to_graph: List[str] = [],
    logger: Callable[[str], None] = Logger(),
    use_daily_loop: bool = False,
    *args,
    **kwargs,
) -> Tuple[Model_State_Shape, List[List[Any]], Config_Shape, Model_State_Shape]:
    """Run the full model from a config file location and external data location.

    Takes the config, data and output directory locations.
    Initializes the config then runs the model
    Outputs are saved in the output directory

    Parameters
    ----------
    project_paths: GridProjectPaths
        file paths specific to project
    run_paths : GridRunPaths
        file paths specific to run
    runnotes : str, optional
        run notes to add to log out
    output_results_only: bool
        If true then only outputs the results file
    logger: Callable[[str], None], optional
        Log function, by default print
    use_daily_loop: bool, optional
        If true then loops over daily processes instead of pre generating processes for entire run

    *args, **kwargs
        Passed to Main_Overrides


    Returns
    -------
    Tuple[Model_State_Shape, List[List[Any]], Config_Shape, Model_State_Shape]
       final_state, output_logs, config, initial_state

    """
    overrides = Main_Overrides(*args, **kwargs)

    # MODEL SETUP
    loaded_run_files = load_run_files(
        project_paths=project_paths,
        run_paths=run_paths,
    )

    start_time_setup = datetime.now()
    [
        config,
        external_state,
        initial_state,
        model_processes,
    ] = setup_model(
        config_in=loaded_run_files.config,
        state_in=loaded_run_files.state,
        data_location=run_paths.input_data_file_path,
        use_daily_loop=use_daily_loop,
        overrides=overrides,
    )

    time_taken_setup = datetime.now() - start_time_setup
    # MODEL RUN
    start_time = datetime.now()
    model_runner = run_model_daily if use_daily_loop else run_model
    final_state, output_logs = model_runner(
        initial_state,
        config,
        external_state,
        model_processes,
        DEBUG_MODE=overrides.debug,
    )
    time_taken = datetime.now() - start_time
    log_notes = generate_run_notes(
        runnotes, time_taken, time_taken_setup, config.VERSION, model_version)
    output_filename = f'{run_paths.run_id}.csv' if output_results_only else f'{run_paths.run_id}_out.csv'

    if run_paths.output_directory:
        export_output(
            output_logs,
            final_state,
            external_state,
            run_paths.output_directory,
            final_config=config,
            output_filename=output_filename,
            output_results_only=output_results_only,
            fields_to_graph=fields_to_graph,
            runid=run_paths.run_id,
            log_notes=log_notes,
        )

    return final_state, output_logs, config, initial_state


def main_partial(
    cell_configs: Config_Shape,
    grid_coords: List[Tuple[int, int]],
    output_shape: Tuple[int, int],
    external_data_file_path: Path,
    previous_hour_state_path: Path,
    output_directory: Path,
    output_fields: List[str],
    external_state_options: EStateOptions,
    state_out_path: Path = None,
    logger: Callable[[str], None] = Logger(),
) -> Tuple[Model_State_Shape, List[List[Any]], Config_Shape, Model_State_Shape]:
    """Run the DO3SE model on partial data.

    Use case: we have netcdf data in hourly or monthly batches.


    We pass all the inputs to a setup func here that returns an iterator. Each pass of the
    iterator is a model setup for a single grid cell.
    In the model setup are the config, external_state, initial_state, model_processes and grid coords (x, y).

    Parameters
    ----------
    cell_configs : Config_Shape
        Loaded configs for this run
    grid_coords: List[Tuple[int, int]]
        Grid coordinates to run
    output_shape: Tuple[int,int]
          of the output grid
    external_data_file_path : Path
        The path the the input data file
    previous_hour_state_path: Path
        The path to the directory with previous hour states
    output_directory : Path, optional
        The directory to save outputs, by default None
    output_fields : List[str], optional
        A list of fields to output
    external_state_options: EStateOptions
        Options for loading external data
    state_out_path: Path, optional
        If provided output state is saved to this path else overrides previous hour state
    logger: Callable[[str], None], optional
        Log function, by default print

    Returns
    -------
    Tuple[Model_State_Shape, List[List[Any]], Config_Shape, Model_State_Shape]
        [description]

    """
    [
        start_day,
        end_day,
        start_date,
        end_date,
        row_count,
        time_string,
        hours,
    ] = get_date_bounds_from_ext_data(
        # TODO: Check passing zero year
        external_data_file_path,
        grid_coords,
        external_state_options,
    )
    T = row_count

    # This is an iterator over each grid cell.
    external_state_data_iter = load_external_state(
        external_data_file_path,
        grid_coords,
        external_state_options,
        # TODO: All below being passed correctly in external_state_options
        # file_type=overrides.external_file_type,
        # coords=grid_coords,
        # variable_map=overrides.netcdf_variable_map,
        # preprocess_map=overrides.met_preprocess_map,
        # multi_file_data=overrides.multi_file_netcdf,
        # zero_year=base_config.Location.zero_year,
        # netcdf_chunks=overrides.netcdf_chunks,
    )
    e_state_init_processes = (external_state_init_processes(0, row_count, config.Met)
                              for config in cell_configs)

    prev_hour_states_loaded = (load_current_cell_state(previous_hour_state_path, x, y)
                               for x, y in grid_coords)

    # TODO: Optimize this by just getting a day of processes
    model_processes = (get_row_processes_hourly(config, hours)
                       for config in cell_configs)

    # Process external state ready for model run.
    external_state_loaded = (setup_external_state_simple(
        external_state_data,
        config,
        e_state_processes_cell,
    ) for config, external_state_data, e_state_processes_cell
        in zip(cell_configs, external_state_data_iter, e_state_init_processes))

    # Setup output Data
    grid_coords = np.array(grid_coords)
    grid_x_min, grid_y_min = grid_coords.min(axis=0)
    # grid_x_size, grid_y_size = np.ptp(grid_coords, axis=0) + [1, 1]

    grid_x_size, grid_y_size = output_shape
    full_output_data = {
        k: np.full((grid_x_size, grid_y_size, T), None,
                   dtype=np.float64) for k in output_fields
    } if output_fields else None

    lat_data = np.full((grid_x_size, grid_y_size), None, dtype=float)
    lon_data = np.full((grid_x_size, grid_y_size), None, dtype=float)

    output_state_dir = state_out_path or previous_hour_state_path

    # == Run each cell
    # TODO: can we run in parallel?
    for cell_model in zip(
        grid_coords, cell_configs, external_state_loaded, prev_hour_states_loaded, model_processes
    ):
        # Each iteration runs a single grid cell over all the hours in the input data.
        [(xi, yi), config_cell, external_state_data_cell,
         prev_hour_state_loaded, cell_processes] = cell_model

        logger(f'Running coord ({xi}, {yi} for date: {start_date})', verbose=True)
        xi = xi - grid_x_min
        yi = yi - grid_y_min

        # Make sure start and end day set in config.
        config_cell.Location.start_day = start_day
        config_cell.Location.end_day = end_day

        # Run DO3SE model
        final_state, output_logs = run_model(
            prev_hour_state_loaded,
            config_cell,
            external_state_data_cell,
            cell_processes,
        )

        # Save state at each iteration
        state_file_out_path = f"{output_state_dir}/{xi}_{yi}.state"
        dump_state_to_file_quick(final_state, state_file_out_path)

        # TODO: Can we save to netcdf here rather than storing all output in memory till end of runs?
        for k in output_fields:
            # OPTIMIZE: Can we optimize this?
            # Output logs is List[dict[str, float]]
            full_output_data[k][xi, yi] = [o[k] for o in output_logs]
        lat_data[xi, yi] = 0  # TODO: Get lat and lon data
        lon_data[xi, yi] = 0  # TODO: Get lat and lon data

        yield final_state, output_logs, cell_model

    # After all grid cells ran we save to netcdf
    if output_directory:
        # TODO: Check time of output
        time_data = pd.date_range(start_date, periods=T, freq="1H")
        target_file_name = f'output_data_{time_string}.nc'

        dump_output_to_file_netcdf_grid(
            output_data=full_output_data,
            X=grid_x_size, Y=grid_y_size, T=T,
            lat_data=lat_data,
            lon_data=lon_data,
            time_data=time_data,
            target_path=f"{output_directory}/{target_file_name}",
            output_fields=output_fields,
        )

    return full_output_data


def main_partial_initialize(
    config: Config_Shape,
    state: Model_State_Shape,
    processed_config_dir: Path,
    state_out_path: Path,
    e_state_overrides_file_path: Path = None,
    e_state_overrides_field_map: dict = None,
    grid_coords: List[Tuple[int, int]] = None,
    logger: Callable[[str, str], None] = Logger(),
    *args,
    **kwargs,
):
    """Setup initial state and config for hourly run.


    This should only be ran once before the model iterations.

    The hourly run requires the initial state and config for each grid to be saved to a seperate file.
    It also requires the config to be preprocessed.

    Parameters
    ----------
    config_file_path : Path
        Path of file to run
    processed_config_dir : Path
        Location to save processed configs
    state_out_path : Path
        location to save initial state files
    base_config_path : Path, optional
        path to base config file, by default None
    base_state_path : Path, optional
        path to base initial state file, by default None
    e_state_overrides_file_path : Path, optional
        path to netcdf file with per coord overrides, by default None
    e_state_overrides_field_map : dict, optional
        map of netcdf fields to Config_Shape fields, by default None
    grid_coords: List[Tuple[int, int]]
        grid coordinates
    logger : Callable[[str, str], None], optional
        logger to use, by default print

    """
    overrides = Main_Overrides(*args, **kwargs)

    logger('===== Setting up Grid config =====')

    # NOTE: We initialize all phenology config here even if some values overriden later.
    # This also saves out per cell configs
    # config_loaded = setup_config(run_paths.config_path, None, project_paths.base_config_path, overrides),
    initialized_config = setup_config(config, None, overrides)
    initialize_grid_configs(
        initialized_config,
        processed_config_dir,
        grid_coords,
        e_state_overrides_file_path,
        e_state_overrides_field_map,
        logger,
    )

    # TODO: Need state per config file
    # TODO: Do we need to run this on a per cell/per config basis?
    logger('===== Setting up initial state =====')

    initialized_state = setup_initial_state(state, initialized_config, None, True)

    logger('======= Copying State for each grid ===========')
    for xi, yi in grid_coords:
        state_file_out_path = f"{state_out_path}/{xi}_{yi}.state"
        sys.stdout.write(f"\r {state_file_out_path}")
        sys.stdout.flush()
        dump_state_to_file_quick(initialized_state, state_file_out_path)

    logger('===== Completed setting up initial state =====')

    return initialized_state


# @deprecated.deprecated()
# def main_hour(
#     processed_config_dir: Path,
#     external_data_row_path: Path,
#     previous_hour_state_path: Path,
#     output_data_dir: Path = None,
#     output_fields: List[str] = [],
#     logger: Callable[[str], None] = Logger(),
#     *args,
#     **kwargs,
# ) -> List[List[Model_State_Shape]]:
#     """Run an hour of model.


#     # NOTE: If grid_coords are sparse, output will contain NaN values where cell has not been ran.

#     Parameters
#     ----------
#     processed_config_dir : Path
#         path to pre-processed config
#     external_data_row_path : Path
#         path to data for current hour
#     previous_hour_state_path : Path
#         path to previous hour final state
#     output_data_dir: Path
#         path to file to save output data, by default None
#     output_fields : List[str], optional
#         Fields to graph, by default []
#     logger: Callable[[str], None], optional
#         Log function, by default print

#     Returns
#     -------
#     Tuple[Model_State_Shape, List[List[Any]], Config_Shape, Model_State_Shape]
#        final_state, output_logs, config, initial_state

#     """
#     overrides = Main_Overrides(*args, **kwargs)

#     # == Model setup creates a run setup for each grid square and config option
#     # TODO: This is not returning the model processes required.
#     models = setup_model_single_hour_grid(
#         processed_config_dir,
#         external_data_row_path,
#         previous_hour_state_path,
#         overrides=overrides,
#     )

#     # TODO: Below only works if grid method. Should split grid and none grid hour runs
#     if (overrides.grid_coords is None):
#         raise NotImplementedError("Only implemented for grid setup!")

#     grid_coords = np.array(overrides.grid_coords)
#     grid_x_min, grid_y_min = grid_coords.min(axis=0)
#     grid_x_size, grid_y_size = np.ptp(grid_coords, axis=0) + [1, 1]

#     # TODO: Can make this more efficient by using a predefined schema for the dict to supply to np.full.
#     # TODO: Need to manage T > 1
#     # To implement T then we need to ensure that the models run in order
#     # Models could instead be an iterator for each grid cell. I.e shape would be (nX*nY, nT)
#     T = 1
#     full_output_data = np.full((grid_x_size, grid_y_size, T), dict())
#     lat_data = np.full((grid_x_size, grid_y_size), None)
#     lon_data = np.full((grid_x_size, grid_y_size), None)

#     output_state_dir = overrides.state_out_path or previous_hour_state_path

#     # == Run each cell
#     # TODO: can we run in parallel?
#     for model in models:
#         [
#             config,
#             external_state,
#             initial_state,
#             model_processes,
#             x, y
#             # TODO: Get base netcdf file for output
#         ] = model
#         model_date = external_state.date[0]  # TODO: Get date from inputs
#         logger(f'Running coord ({x}, {y} for date: {model_date})', verbose=True)
#         xi = x - grid_x_min
#         yi = y - grid_y_min
#         final_state, output_logs = run_model(
#             deepcopy(initial_state), config, external_state, model_processes)
#         # Save state at each iteration
#         state_file_out_path = f"{output_state_dir}/{x}_{y}.state"
#         dump_state_to_file_quick(final_state, state_file_out_path)

#         # TODO: Can we save to netcdf here rather than storing all output in memory till end of runs?
#         full_output_data[xi, yi] = output_logs
#         lat_data[xi, yi] = 0  # TODO: Get lat and lon data
#         lon_data[xi, yi] = 0  # TODO: Get lat and lon data

#         yield final_state, output_logs, model

#     # After all grid cells ran we save to netcdf
#     if output_data_dir:
#         time_data = pd.date_range(model_date, periods=1)
#         time_string = f'{time_data[0].year}-{str(time_data[0].month).zfill(2)}-{str(time_data[0].day).zfill(2)}_{str(time_data[0].hour).zfill(2)}'
#         target_file_name = f'output_data_{time_string}.nc'
#         dump_output_to_file_netcdf_grid(
#             full_output_data,
#             lat_data=lat_data,
#             lon_data=lon_data,
#             time_data=time_data,
#             target_path=f"{output_data_dir}/{target_file_name}",
#             output_fields=output_fields,
#         )
#     return full_output_data

# @deprecated.deprecated()
# def multi_run(
#     configs: List[Union[str, Config_Shape]],
#     data_location: str,
#     *args,
#     **kwargs,
# ) -> bool:
#     """Run multiple model passes for each config in list."""
#     overrides = Main_Overrides(*args, **kwargs)
#     for config_in in configs:
#         if not isinstance(config_in, (str, Config_Shape)):
#             raise ValueError('Config must be Config_Shape or string path')
#         config = config_in if isinstance(
#             config_in, Config_Shape) else setup_config(config_in, overrides)
#         # TODO: Cache setup external state
#         external_state, start_day, end_day = setup_external_state(config, data_location, overrides)
#         initial_state = setup_initial_state(
#             config, external_state, run_init_processes=True, overrides=overrides,
#         )
#         model_processes = setup_model_processes(config, overrides)
#         run_model(initial_state, config, external_state, model_processes)

#     return True


def main_grid_seq_per_config(
    project_paths: GridProjectPaths,
    run_paths: GridRunPaths,
    loaded_run_files: GridRunFiles,
    grid_coords: List[Tuple[int, int]],
    output_shape: Tuple[int, int],
    runnotes: str = '',
    output_fields: List[str] = None,
    log_level: int = 0,
    seperate_live_state: bool = False,
    multi_file_netcdf: bool = False,
    logger: Callable[[str], None] = print,
):
    """Run the grid run on a single config.

    Parameters
    ----------
    project_paths: GridProjectPaths
        file paths specific to project
    run_paths : GridRunPaths
        file paths specific to run
    loaded_run_files : GridRunFiles
        loaded run files
    grid_coords : List[Tuple[int, int]]
        grid coordinates to run
    output_shape: Tuple[int,int]
          of the output grid
    runnotes : str, optional
        Additional notes to store with output, by default ''
    output_fields : List[str], optional
        Output state fields to save, by default None
    log_level : int, optional
        log level, by default 0
    seperate_live_state : bool, optional
        If true then don't overwrite prev hour state with new state, by default False
    multi_file_netcdf : bool, optional
        If true use xr.open_mfdataset, by default False
    logger : Callable[[str], None], optional
        logger func, by default print

    """

    log_path = run_paths.log_path  # f"{run_dir}/run.log"
    logger = Logger(log_level, log_path, write_mode='w', set_as_default=True)

    cell_configs: List[Config_Shape] = grid_config_loader(
        run_paths.processed_configs_dir, grid_coords)
    zero_year = cell_configs[0].Location.zero_year  # TODO: Get this from config
    errors = []
    input_files = get_input_files_list(project_paths.input_data_dir, multi_file_netcdf)
    total_files = len(input_files)
    logger(f"== Running Model: Total File: {total_files} ===")

    external_state_options = EStateOptions(
        file_type=FileTypes.NETCDF,
        variable_map=loaded_run_files.variable_map,
        multi_file_data=multi_file_netcdf,
        preprocess_map=loaded_run_files.preprocess_map,
        zero_year=zero_year,
        # file_suffix=input_file_suffix, # May not be needed
        # netcdf_chunks=None,
    )

    duration = None
    try:
        start_time = datetime.now()

        # Iterate over each input file
        # NOTE: Assumes files sorted will be in chronological order.
        # TODO: Fix this for when using multi file loading
        for i, f in enumerate(input_files):
            logger(f"==== Running file: {f} ({i+1}/{total_files}) =======")
            state_out_path = run_paths.live_state_dir
            input_data_file = f"{project_paths.input_data_dir}/{f}" if not multi_file_netcdf else project_paths.input_data_dir
            previous_hour_state_path = run_paths.live_state_dir if not seperate_live_state else run_paths.prev_state_dir

            if multi_file_netcdf:
                external_state_options = external_state_options._replace(
                    # file_suffix=f
                    data_filter=f,
                )

            if seperate_live_state:
                shutil.rmtree(run_paths.prev_state_dir)
                os.rename(run_paths.live_state_dir, run_paths.prev_state_dir)
                os.makedirs(run_paths.live_state_dir, exist_ok=True)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for out in main_partial(
                    cell_configs=cell_configs,
                    grid_coords=grid_coords,
                    output_shape=output_shape,
                    external_data_file_path=input_data_file,
                    previous_hour_state_path=previous_hour_state_path,
                    output_directory=run_paths.output_data_dir,
                    output_fields=output_fields,
                    external_state_options=external_state_options,
                    state_out_path=state_out_path,
                    logger=logger,
                ):
                    final_state, output_logs, model = out

            logger(f"==== Running file complete: {f} ({i+1}/{total_files}) =======")
            logger(f"==== Files saved in {run_paths.output_data_dir} =======")

        end_time = datetime.now()
        duration = end_time - start_time if end_time else 0

    except Exception as e:
        errors.append((f"Run Dir: {run_paths.run_dir} failed", e))

    if log_level >= 0:
        with open(f'{run_paths.output_data_dir}/notes.log', 'w') as f:
            log_notes = generate_run_notes(
                runnotes,
                time_taken=duration,
                time_taken_setup=0,
                config_version=config_version,
                model_version=model_version,
                errors=errors,
            )
            f.write("\n".join(log_notes))

        if len(errors) > 0:
            for m, e in errors:
                logger(m)
            logger(errors)
            raise errors[0][1]
    logger.close


def main_grid_run(
    project_paths: GridProjectPaths,
    output_fields: List[str],
    multi_file_netcdf: bool = False,
    runid: str = 'run',
    runnotes: str = '',
    log_level: int = 0,
    seperate_live_state: bool = False,
    save_final_state: bool = False,
):
    """Run model for grid data split per hour.

    Parameters
    ----------
    project_paths: GridProjectPaths
        file paths specific to project
    output_fields : List[str], optional
        Fields to output
    multi_file_netcdf : bool, optional
        If true use xr.open_mfdataset, by default False
    runid : str, optional
        Unique id for run, by default ''
    runnotes : str, optional
        Additional notes to save with output, by default ''
    log_level : int, optional
        Log level, by default 0
    seperate_live_state : bool, optional
        If true then will save live state in seperate directory(This is more expensive), by default False

    """
    # TODO: Check this is all working still
    logger_main = Logger(log_level)
    # Inputs

    config_file_names = get_configs(project_paths.config_dir)

    logger_main(f'== Found {len(config_file_names)} configs to run =====')

    for config_file_name in config_file_names:
        config_name = '.'.join(config_file_name.split('.')[:-1])
        run_paths = get_grid_run_paths(project_paths, config_name, runid)
        loaded_run_files = load_grid_run_files(project_paths, run_paths)
        create_grid_run_path_directories(run_paths)
        grid_coords, *output_shape = get_grid_coords_from_file(
            run_paths.run_mask_path
        )

        main_grid_seq_per_config(
            project_paths=project_paths,
            run_paths=run_paths,
            loaded_run_files=loaded_run_files,
            grid_coords=grid_coords,
            output_shape=output_shape,
            runnotes=runnotes,
            output_fields=output_fields,
            log_level=log_level,
            seperate_live_state=seperate_live_state,
            multi_file_netcdf=multi_file_netcdf,
            logger=logger_main,
        )
        if save_final_state:
            for x, y in grid_coords:
                file_name = f"{x}_{y}"
                dump_state_to_file(model_state_loader_quick(
                    f"{run_paths.live_state_dir}/{file_name}.state"), f"{run_paths.final_state_dir}/{file_name}.json")

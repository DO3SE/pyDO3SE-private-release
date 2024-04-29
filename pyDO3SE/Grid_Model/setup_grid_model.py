from dataclasses import dataclass
from functools import partial
import os
import math
from copy import deepcopy
import numpy as np
from pathlib import Path
from typing import Generator, Iterator, Tuple, Callable, List, NamedTuple
from data_helpers.cls_parsing import rsetattr
from datetime import datetime
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.util.logger import Logger

from pyDO3SE.util.loader import json_loader
from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.util.error_handling import InputDataError, Do3seSetupError
from pyDO3SE.Config import Config_Shape
from pyDO3SE.Output.process_outputs import dump_config_to_file_binary
from pyDO3SE.External_State.external_state_loader import Coord
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.setup_model import (
    Main_Overrides,
    setup_config,
    setup_initial_state,
)
from pyDO3SE.Model_State.model_state_loader import (
    dump_state_to_file_quick,
    model_state_loader,
)


@dataclass
class ParallelArgs:
    """ Additional args for run in parallel

    Parameters
    ----------
    MAX_PROCESSES: int, default 8
        Max number of processes
    TIMEOUT: float
        NOT IMPLEMENTED
    SLEEP_TIME: float
        Time to wait for process to finish. lower values will increase responsiveness
        to runs finishing but will increase cpu on listener thread

    """
    MAX_PROCESSES: int = 8
    TIMEOUT: float = 1
    SLEEP_TIME: float = 0.01


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
    run_id: str
    run_dir: str
    log_path: str


class GridRunPaths(NamedTuple):
    config_id: str
    config_run_dir: str
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


def pull_config_vars_from_netcdf(
    ds: 'xr.DataSet',
    coord: Coord,
    e_state_overrides_field_map: dict,
) -> dict:
    try:
        x, y = coord
        return {k: ds[v].item((x, y)) for k, v in e_state_overrides_field_map.items()}
    except TypeError as e:
        raise e
    except ValueError as e:
        print(ds)
        raise InputDataError('Check shape of e_state_override netcdf vars') from e
    except KeyError as e:
        raise InputDataError(
            f"Check all keys in e_state_overrides_field_map are in overrides dataset, \nkeys: {e_state_overrides_field_map.keys()}, \n{list(ds.keys())}") from e
    except Exception as e:
        print(coord)
        print(e_state_overrides_field_map)
        print(ds)
        raise Do3seSetupError("Failed to load config override map")


def load_additional_gridded_config_data(
    coords: List[Coord],
    e_state_overrides_dataset: 'xr.DataSet',
    e_state_overrides_field_map: dict,
) -> Iterator[Tuple[Coord, dict]]:
    """Load additional per grid cell config from netcdf file.


    Example e_state_overrides_field_map:

        {
            'Location.lat': "latitude",
            'Location.lon': "longitude",
        }

    Example output dict:

        {
            'Location.lat': 43.1,
            'Location.lon': -2.3,
        }


    Parameters
    ----------
    coords : List[Coord]
        Coordinates in netcdf file to load values from.
    e_state_overrides_dataset : xr.DataSet
        netcdf Dataset
    e_state_overrides_field_map : dict
        Dictionary mapping of field in netcdf to field in Config_Shape

    Yields
    -------
    Iterator[Tuple[Coord, dict]]
        An iterator where each iteration corresponds to a coordinate.
        Each iteration contains a tuple of the coordinate and a map of config var to value.

    """
    for coord in coords:
        field_map = pull_config_vars_from_netcdf(
            e_state_overrides_dataset, coord, e_state_overrides_field_map)
        yield (coord, field_map)


def get_grid_coords_from_dataarray(
    da: "xr.DataArray",
) -> List[Tuple[int, int]]:
    """Get the grid coordinates and size from a coordinates file

    Parameters
    ----------
    da : xr.DataArray
        Xr DataArray with mask of cells to run.

    Returns
    -------
    Tuple[np.ndarray, int, int]
        grid_coords, grid_x_size, grid_y_size

    """
    assert len(
        da.data.shape) == 2, f"Coordinates netcdf file must have shape (x, y) but has shape {da.data.shape}"

    x_indices, y_indices = np.indices(da.data.shape)
    flat_mask = da.data.astype(bool).flatten()

    x_indices_masked = np.ma.masked_array(x_indices.flatten(), ~flat_mask)
    y_indices_masked = np.ma.masked_array(y_indices.flatten(), ~flat_mask)
    x_indices_masked.compressed()[0:10]
    grid_coords = list(zip(x_indices_masked.compressed(), y_indices_masked.compressed()))

    # grid_coords = np.array([[int(i['x']), int(i['y'])] for i in csv_loader(run_mask_path)])
    grid_x_size, grid_y_size = da.shape

    return grid_coords, grid_x_size, grid_y_size


def get_grid_coords_from_file(
    run_mask_path: str,
) -> Tuple[np.ndarray, int, int]:
    """Get the grid coordinates and size from a coordinates file

    Parameters
    ----------
    run_mask_path : str
        Path to a netcdf run mask

    Returns
    -------
    Tuple[np.ndarray, int, int]
        grid_coords, grid_x_size, grid_y_size

    """
    da = xr.open_dataarray(run_mask_path, engine="netcdf4")
    return get_grid_coords_from_dataarray(da)


def save_configs_from_generator(
    configs: Generator[Config_Shape, None, None],
    coords: List[Tuple[int, int]],
    output_dir: Path,
):
    """Save all the configs stored in a generator.

    Also returns the last config processed.

    Parameters
    ----------
    configs : Generator[Config_Shape]
        _description_
    coords : List[Tuple[int, int]]
        _description_
    output_dir : Path
        _description_

    """
    for config, [x, y] in zip(configs, coords):
        config_output_path = f"{output_dir}/{x}_{y}.config"
        dump_config_to_file_binary(config, config_output_path)


def save_state_from_generator(
    states: Generator[Model_State_Shape, None, None],
    coords: List[Tuple[int, int]],
    output_dir: Path,
):
    """Save all the configs stored in a generator.

    Also returns the last config processed.

    Parameters
    ----------
    states : Generator[Model_Shape]
        _description_
    coords : List[Tuple[int, int]]
        _description_
    output_dir : Path
        _description_

    """
    for state, [x, y] in zip(states, coords):
        state_output_path = f"{output_dir}/{x}_{y}.state"
        dump_state_to_file_quick(state, state_output_path)


def initialize_grid_configs(
    base_config: Config_Shape,
    grid_coords: List[Tuple[int, int]],
    e_state_overrides_dataset: 'xr.DataSet' = None,
    e_state_overrides_field_map: dict = None,
    logger: Callable[[str, str], None] = print,
    overrides: Main_Overrides = Main_Overrides(),
) -> Generator[Config_Shape, None, None]:
    """Initialize and save a config file for each grid coord.

    This is for multi grid runs where we are using grid data such as sowing date
    to modify the config parameters.

    Parameters
    ----------
    base_config: Config_Shape
        Config common to all cells
    grid_coords: List[Tuple[int, int]]
        grid coordinates
    e_state_overrides_dataset : xr.DataSet, optional
        netcdf dataset with per coord overrides, by default None
    e_state_overrides_field_map : dict, optional
        map of netcdf fields to Config_Shape fields, by default None
    logger : Callable[[str, str], None], optional
        logger, by default print
    overrides: Main_Overrides
        Any high level model overrides

    Returns
    -------
    Generator[Config_Shape]
        generator of configs

    Raises
    ------
    AssertionError
        No configs found
    NotImplementedError
        Multi config not implemented

    """
    logger("Setting up cell config")

    # = We can supply a netcdf file to overrides per coord configs
    config_overrides = load_additional_gridded_config_data(
        grid_coords,
        e_state_overrides_dataset,
        e_state_overrides_field_map,
    ) if e_state_overrides_dataset is not None \
        else [(coord, {}) for coord in grid_coords]

    initialized_config = None
    for (xi, yi), override_dict in config_overrides:
        logger(f"Setting up config for ({xi},{yi})")
        logger(f"override_dict: {override_dict}")
        config_cell = deepcopy(base_config)
        for k, v in override_dict.items():
            if v is None or math.isnan(v):
                raise ValueError(f'{k} is invalid for ({xi},{yi})')
            config_cell = rsetattr(config_cell, k, v, True)

        # TODO: We only need to set up for each cell if the overrides effect the setup config
        initialized_config = setup_config(
            config_cell, external_state_data=None, overrides=overrides)
        yield initialized_config


def init_all_grid_model_configs(
    project_paths: GridProjectPaths,
    logger: Logger = print,
    sample_size: int = 0,
):
    """Initialize grid state for all configs in project dir

    Parameters
    ----------
    project_paths: GridProjectPaths
        file paths specific to project
    logger: Logger,
        Logger func or class
    sample_size: int, optional
        If greater than 0 then only runs up to sample size number of cells

    """

    logger(f'== Running init_all_grid_model_configs =====')
    configs = os.listdir(project_paths.config_dir)
    logger(f'== Found {len(configs)} configs to run =====')
    setup_duration = datetime.now() - datetime.now()
    e_state_overrides_dataset = xr.open_dataset(project_paths.e_state_overrides_file_path)
    for config_file_path in configs:
        config_name = '.'.join(config_file_path.split('.')[:-1])
        run_paths = get_grid_run_paths(project_paths, config_name)
        try:
            create_grid_run_path_directories(run_paths)
            loaded_files: GridRunFiles = load_grid_run_files(project_paths, run_paths)
            # logger = Logger(log_level, run_paths.log_path, write_mode='w', set_as_default=True)
            logger(f'== Running Init on config {config_file_path} ==')

            grid_coords, _, _ = get_grid_coords_from_file(run_paths.run_mask_path)
            grid_coords_sampled = grid_coords[0:sample_size]if sample_size else grid_coords

            logger(f"Initializing grid setup for: \n{run_paths.config_run_dir}")
            start_time = datetime.now()

            initialized_config_gen, initialized_state_gen = init_grid_model(
                config=loaded_files.config,
                state=loaded_files.state,
                e_state_overrides_dataset=e_state_overrides_dataset,
                e_state_overrides_field_map=loaded_files.e_state_overrides_field_map,
                grid_coords=grid_coords_sampled,
                logger=logger,
                debug=logger.log_level >= 2,
            )
            save_configs_from_generator(
                initialized_config_gen,
                grid_coords_sampled,
                run_paths.processed_configs_dir,
            )
            save_state_from_generator(
                initialized_state_gen,
                grid_coords_sampled,
                run_paths.live_state_dir,
            )
            end_time = datetime.now()
            setup_duration_run = end_time - start_time
            logger(
                f"Initialization complete for: \n{run_paths.config_run_dir}\nSetup took: ({setup_duration})")

        except Exception as e:
            logger.close()
            raise Do3seSetupError(
                f"Failed to run grid init for config: {config_file_path}\nRun Paths: {run_paths}", e)
        setup_duration += setup_duration_run
    logger(
        f'== COMPLETE - Running init_all_grid_model_configs\nSetup took: {setup_duration} ===== \n')
    return setup_duration


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


def get_grid_project_paths(
    project_dir: str,
    run_id: str,
) -> GridProjectPaths:
    run_dir = f"{project_dir}/runs/{run_id}"
    os.makedirs(run_dir, exist_ok=True)

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
        run_dir=run_dir,
        run_id=run_id,
        log_path=f"{run_dir}/run.log",
    )


def get_grid_run_paths(
    project_paths: GridProjectPaths,
    config_id: str,
) -> GridRunPaths:
    config_run_dir = f"{project_paths.run_dir}/{config_id}"

    return GridRunPaths(
        config_id=config_id,
        config_run_dir=config_run_dir,
        config_path=f"{project_paths.project_dir}/configs/{config_id}.json",
        e_state_overrides_field_map_path=f"{project_paths.e_state_overrides_field_map_path}/{config_id}.json",
        initial_state_dir=f"{config_run_dir}/initial_state",
        live_state_dir=f"{config_run_dir}/current_state",
        final_state_dir=f"{config_run_dir}/final_state",
        prev_state_dir=f"{config_run_dir}/prev_state",
        output_data_dir=f"{config_run_dir}/outputs_grid",
        processed_configs_dir=f"{config_run_dir}/processed_configs",
        run_mask_path=f"{project_paths.project_dir}/run_masks/{config_id}.nc",
    )


def create_grid_run_path_directories(
    run_paths: GridRunPaths,
):
    os.makedirs(run_paths.config_run_dir, exist_ok=True) if run_paths.config_run_dir else None
    os.makedirs(run_paths.processed_configs_dir,
                exist_ok=True) if run_paths.processed_configs_dir else None
    os.makedirs(run_paths.prev_state_dir, exist_ok=True) if run_paths.prev_state_dir else None
    os.makedirs(run_paths.live_state_dir, exist_ok=True) if run_paths.live_state_dir else None
    os.makedirs(run_paths.output_data_dir, exist_ok=True) if run_paths.output_data_dir else None
    os.makedirs(run_paths.initial_state_dir, exist_ok=True) if run_paths.initial_state_dir else None
    os.makedirs(run_paths.final_state_dir, exist_ok=True) if run_paths.final_state_dir else None


def init_grid_model(
    config: Config_Shape,
    state: Model_State_Shape,
    e_state_overrides_dataset: 'xr.DataSet' = None,
    e_state_overrides_field_map: dict = None,
    grid_coords: List[Tuple[int, int]] = None,
    logger: Callable[[str, str], None] = Logger(),
    *args,
    **kwargs,
) -> Tuple[
    Generator[Config_Shape, None, None],
    Generator[Model_State_Shape, None, None],
]:
    """Setup initial state and config for hourly run.

    This should only be ran once before the model iterations.

    The hourly run requires the initial state and config for each grid to be saved to a seperate file.
    It also requires the config to be preprocessed.

    Parameters
    ----------
    config_file_path : Path
        Path of file to run
    e_state_overrides_dataset: xr.DataSet, optional
        Netcdf Dataset with config overrides
    e_state_overrides_field_map : dict, optional
        map of netcdf fields to Config_Shape fields, by default None
    grid_coords: List[Tuple[int, int]]
        grid coordinates
    logger : Callable[[str, str], None], optional
        logger to use, by default print

    """
    overrides = Main_Overrides(*args, **kwargs)

    logger('Setting up grid config')

    # NOTE: We are calling this twice so we can use it for state init also
    initialized_config_gen = partial(
        initialize_grid_configs,
        base_config=config,
        grid_coords=grid_coords,
        e_state_overrides_dataset=e_state_overrides_dataset,
        e_state_overrides_field_map=e_state_overrides_field_map,
        logger=logger,
        overrides=overrides,
    )

    logger('Setting up grid state')
    initialized_state_gen = (setup_initial_state(
        state, config, external_state=None, run_init_processes=True, overrides=overrides)
        for config in initialized_config_gen())

    logger('init_grid_model - COMPLETE')

    return initialized_config_gen(), initialized_state_gen

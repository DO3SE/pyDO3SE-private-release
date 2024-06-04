"""Model setup function.

Returns the arguments for the run model function.
"""
import os
import re
from typing import Iterator, Tuple
from copy import deepcopy
import deprecated
import numpy as np
import math
from collections import namedtuple
from dataclasses import asdict
from pathlib import Path
from typing import Callable, List, NamedTuple
from enum import Enum
from math import floor

from data_helpers.cls_parsing import rsetattr
from proflow.ProcessRunnerCls import ProcessRunner
from proflow.Objects.Process import Process

from pyDO3SE.Output.process_outputs import dump_config_to_file_binary
from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.version import config_version
from pyDO3SE.error_handling import DayRangeError, InputDataError
from pyDO3SE.Config import Config_Shape
from pyDO3SE.Config.config_loader import config_loader_pickled
from pyDO3SE.External_State import External_State_Shape
from pyDO3SE.External_State.external_state_loader import Coord, load_external_state
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.Pipelines.default_processes import (
    full_model_processes,
)
from pyDO3SE.Pipelines.state_init_processes import state_init_processes
from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes
from pyDO3SE.Pipelines.config_init_processes import config_init_processes

Model = namedtuple('Output', 'config,external_state,initial_state,model_processes')


GET_INIT_PROCESSES_T = Callable[[Config_Shape], List[Process]]


class LocationMethod(Enum):
    CONFIG_INPUT = "config_input"
    EXTERNAL_STATE_INPUT = "external_state_input"


class Main_Overrides(NamedTuple):
    """Main model overrides and options.

    Parameters
    ----------
    start_day: int
        Override the external data and config start day. default = None
    end_day: int
        Override the external data and config end day. default = None
    init_external_state_processes: List[Process]
        [description] default = None
    init_param_processes: List[Process]
        [description] default = None
    init_state_processes: List[Process]
        [description] default = None
    init_config_processes: List[Process]
        [description] default = None
    processes_in: GET_INIT_PROCESSES_T
        [description] default = None
    config_override: Config_Shape
        [description] default = None
    external_state_override: External_State_Shape
        [description] default = None
    input_state_override: Model_State_Shape
        [description] default = None
    allow_data_padding: bool
        If true then will pad start and end dates in data to match config start and end dates.
        default = False,
    output_to_netcdf: bool, default False
        Output final row data to netcdf if true, csv if false
    state_out_path: str, optional
        Location to save final state
    output_fields: List[str]
        Override the requried output fields
    e_state_overrides_file_path: Path
        Path to netcdf file with gridded external state overrides
    e_state_overrides_field_map: dict
        Dictionary mapping of field in netcdf to field in Config_Shape
    netcdf_chunks: dict, optional
        Chunks to use when loading netcdf data
    debug: bool
        Set debug mode on

    """
    start_day: int = None
    end_day: int = None
    init_external_state_processes: List[Process] = None
    init_param_processes: List[Process] = None
    init_state_processes: List[Process] = None
    init_config_processes: List[Process] = None
    processes_in: GET_INIT_PROCESSES_T = None
    config_override: Config_Shape = None
    external_state_override: External_State_Shape = None
    input_state_override: Model_State_Shape = None
    # If true then will padd start and end dates in data to match config start and end dates
    allow_data_padding: bool = False
    # external_file_type: FileTypes = FileTypes.CSV
    # grid_coords: List[Tuple[int, int]] = None
    # netcdf_variable_map: dict = None
    # met_preprocess_map: dict = None
    # multi_file_netcdf: bool = False
    output_to_netcdf: bool = False
    state_out_path: str = None
    output_fields: List[str] = None
    e_state_overrides_file_path: Path = None
    e_state_overrides_field_map: dict = None
    # netcdf_chunks: dict = None
    debug: bool = False


def setup_config(
    config_in: Config_Shape,
    external_state_data: External_State_Shape,
    overrides: Main_Overrides = Main_Overrides(),
) -> Config_Shape:
    """Setup the model config from input config and data.

    Parameters
    ----------
    config_in: Config_Shape,
        loaded raw config
    external_state_data : Path
        External state data.

    overrides : Main_Overrides
        model overrides

    Returns
    -------
    Config
        initialized config

    """
    process_runner = ProcessRunner()
    process_runner.external_state = external_state_data
    if overrides.config_override:
        return overrides.config_override
    # config = config_loader(config_location, base_config_file, 'json')
    process_runner.config = config_in
    init_processes = overrides.init_config_processes or config_init_processes(config_in)
    config_amended = process_runner.run_processes(
        init_processes,
        config_in)
    if config_amended.VERSION != config_version:
        raise ValueError(
            f"Config must be updated to latest version {config_version}. Use pyDO3SE_cli config-migration [CONFIG_FILE]")
    if overrides.output_fields is not None:
        config_amended.output.fields = overrides.output_fields
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


def setup_dd(
    dd,
    start_day=None,
    end_day=None,
) -> Tuple[List[int], int, int]:
    """Reconstruct dd data to manage wrapping around a year.

    Parameters
    ----------
    dd : List[int]
        external state day data
    start_day : _type_, optional
        start day index inclusive
    end_day : _type_, optional
        end day index inclusive

    Returns
    -------
    Tuple[List[int], int, int]
        _description_

    Raises
    ------
    DayRangeError
        _description_

    """
    dd_out = []
    d_prev = dd[0]
    wrap_d = 0
    start_row = 0
    end_row = None

    if start_day is not None and start_day < dd[0]:
        raise DayRangeError(f"Start day: {start_day} is before start of data: {dd[0]}")

    if end_day is not None and end_day < dd[0]:
        raise DayRangeError(f"End day: {end_day} is before start of data: {dd[0]}")

    for i, d in enumerate(dd):
        if d < d_prev:
            wrap_d = d_prev
        new_dd = d + wrap_d

        if end_day is not None and new_dd > end_day:
            break
        if start_day is None or new_dd >= start_day:
            dd_out.append(new_dd)
        else:
            start_row = i
        end_row = i
        d_prev = d

    if end_day is not None and dd_out[-1] < end_day:
        raise DayRangeError(f"End day: {end_day} is after end of data: {dd_out[-1]}")

    return dd_out, start_row, end_row


def extract_start_and_end_dates(
    config_start_day: int,
    config_end_day: int,
    external_state_dd_data: List[int],
    overrides: Main_Overrides = None,
) -> Tuple[int, int]:
    """Get the start and end dates from ext data taking into account overrides.

    Parameters
    ----------
    config_start_day : int
        config start day overrides
    config_end_day : int
        config start end overrides (inclusive)
    external_state_data : External_State_Shape
        _description_
    overrides : Main_Overrides
        _description_

    Returns
    -------
    Tuple[int, int]
        _description_

    """
    start_day_override = next(int(sd) for sd in [
        overrides and overrides.start_day,
        config_start_day,
        False,
    ] if sd is not None) or None

    end_day_override = next(int(sd) for sd in [
        overrides and overrides.end_day,
        config_end_day,
        False
    ] if sd is not None) or None

    if start_day_override is not None and end_day_override is not None and start_day_override > end_day_override:
        raise DayRangeError(f"Start day: {start_day_override} is after end day {end_day_override}")

    alt_dd, start_row, end_row = setup_dd(
        external_state_dd_data, start_day_override, end_day_override)
    start_day = alt_dd[0]
    end_day = alt_dd[-1]
    assert end_row >= start_row, "No rows!"
    return [start_day, end_day, start_row, end_row, alt_dd]


def setup_external_state_simple(
    external_state_data: External_State_Shape,
    config: Config_Shape,
    init_processes: List[Process],
) -> External_State_Shape:
    process_runner = ProcessRunner(config)
    # process_runner.external_state = external_state_data

    # TODO: replace external_state_init_processes with get_external_...
    external_state = process_runner.run_processes(
        init_processes,
        external_state_data)

    return external_state


def setup_external_state(
    config: Config_Shape,
    external_state_data: External_State_Shape,
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
    -------------------
    - Overrides: `start_day = None, end_day = None`
    - External state `dd = [1,1,1,1...15,15,15,15]`

    The resulting start date and end are `start_day=1, end_day=15` with 15 days of data.


    Parameters
    ----------
    config : Config_Shape
        Config input
    external_state_data : External_State_Shape
        Raw external data
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

    if external_state_data.dd is None:
        raise InputDataError("External state data must include dd")

    [start_day, end_day, start_row, end_row, alt_dd] = extract_start_and_end_dates(
        config.Location.start_day,config.Location.end_day, external_state_data.dd, overrides
    )
    external_state_data.dd = alt_dd
    assert len(alt_dd) > 0, f"Generated alt_dd is empty. Start day is {start_day}, end_day: {end_day}"

    # TODO: Check we are not slicing again here
    external_state_data_sliced = External_State_Shape(
        **{k: (v[start_row:end_row + 1] if v is not None else None)
            for k, v in asdict(external_state_data).items()})

    # Use the process runner to modify the external state data based on some
    # initialization processes
    process_runner = ProcessRunner(config)
    process_runner.external_state = external_state_data_sliced

    # TODO: replace external_state_init_processes with get_external_...
    init_processes = overrides.init_external_state_processes or \
        external_state_init_processes(start_row, end_row + 1, process_runner.config.Met)
    external_state = process_runner.run_processes(
        init_processes,
        external_state_data_sliced)

    return external_state, start_day, end_day

# @deprecated("Use init configs instead")
# def setup_location(
#     method: LocationMethod,
#     config: Config_Shape,
#     external_meta_data: EStateMetaData,
# ) -> Config_Location:
#     if method == LocationMethod.CONFIG_INPUT or method == LocationMethod.CONFIG_INPUT.value:
#         return config.Location
#     if method == LocationMethod.EXTERNAL_STATE_INPUT or method == LocationMethod.EXTERNAL_STATE_INPUT.value:
#         # NOTE: External state type must be NetCDF
#         config_location = deepcopy(config.Location)
#         config_location.elev = external_meta_data.elev
#         config_location.lat = external_meta_data.lat
#         config_location.lon = external_meta_data.lon
#         return config_location
#     raise ConfigError(f"Invalid location method: {method}")


def setup_initial_state(
    state_in: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    run_init_processes: bool = False,
    override_init_state_processes: List[Process] = None,
) -> Model_State_Shape:
    state = deepcopy(state_in)
    if state_in.prev_hour == None:
        # TODO: Check this still works
        state.prev_hour = deepcopy(state_in)
    if run_init_processes:
        process_runner = ProcessRunner(config, external_state)
        init_processes = override_init_state_processes or \
            state_init_processes(config)
        state = process_runner.run_processes(
            init_processes,
            state,
        )
    return state


def setup_model_processes(
    config: Config_Shape,
    overrides: Main_Overrides,
    hours: List[int] = None,
    run_validation: bool = False,
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

    Returns
    -------
    List[Process]
        A list of processes to be ran by the model.

    """
    if overrides.processes_in:
        return overrides.processes_in

    return full_model_processes(config, hours, run_validation)


def setup_model(
    config_in: Config_Shape,
    state_in: Model_State_Shape,
    data_location: Path,
    use_daily_loop: bool = False,
    overrides: Main_Overrides = Main_Overrides(),
) -> Model:
    """Get the model arguments from the input locations.

    Parameters
    ----------
    config_in: Config_Shape
        Merged and loaded config file
    state_in: Model_State_Shape
        Initial loaded state
    data_location: str
        Path to external state file
    use_daily_loop: bool, optional
        If true then loops over daily processes instead of pre generating processes for entire run
    overrides : Main_Overrides
        model overrides

    Returns
    -------
    (config,
    external_state,
    initial_state,
    model_processes)

    """
    # TODO: Implement grid setup here
    external_state_data = next(load_external_state(
        data_location,
    ))
    config = setup_config(
        config_in,
        external_state_data,
        overrides,
    )
    external_state, start_day, end_day = setup_external_state(
        config, external_state_data, overrides)
    config.Location.start_day = start_day
    config.Location.end_day = end_day
    overrides = overrides._replace(start_day=start_day, end_day=end_day)
    initial_state = setup_initial_state(
        state_in,
        config,
        external_state,
        True,
        overrides.init_state_processes,
    )
    hours = external_state.hr if not use_daily_loop else list(range(24))

    model_processes = setup_model_processes(config, overrides, hours)
    return Model(
        config=config,
        external_state=external_state,
        initial_state=initial_state,
        model_processes=model_processes,
    )


# def setup_model_partial_cell(
#     coord: Coord,
#     config: Config_Shape,
#     external_state_data: External_State_Shape,
#     prev_hour_state_loaded: Model_State_Shape,
#     e_state_init_processes: List[Process],
# ) -> Model:
#     """Convert each external state grid iteration into a model run.

#     # TODO: Chech multi hour data works too
#     Parameters
#     ----------
#     coord : Coord
#         The coordinate
#     config: Config_Shape
#         Loaded config for this cell
#     external_state_data : External_State_Shape
#         The external state data for this grid cell
#     previous_hour_state_loaded : Model_State_Shape
#         Path to previous hour state for this grid cell
#     e_state_init_processes: List[Process],
#         List of processes to run on external state

#     Returns
#     -------
#     Model
#         A model setup ready to run

#     """
#     x, y = coord

#     # TODO: Can we optimize this? #OPTIMIZE
#     # TODO: Move some of estate loading higher up.
#     external_state = setup_external_state_simple(
#         external_state_data,
#         config,
#         e_state_init_processes,
#     )
#     # TODO: Can this handle when data wraps around end of year?
#     start_day = external_state.dd[0]
#     end_day = external_state.dd[-1]
#     start_date = external_state.date[0]
#     end_date = external_state.date[-1]

#     row_count = len(external_state.dd)

#     # TODO: Optimize getting model processes.
#     # TODO: Should be same for all hours so can setup higher up.
#     model_processes = get_row_processes(config, start_day, end_day, row_count)
#     config_run = deepcopy(config)
#     config_run.Location.start_day = start_day
#     config_run.Location.end_day = end_day

#     return Model(
#         config=config_run,
#         external_state=external_state,
#         initial_state=prev_hour_state_loaded,
#         model_processes=model_processes,
#         x=x,
#         y=y,
#         ts=start_date,
#         te=end_date,
#     )


# @deprecated.deprecated()
# def setup_model_partial_grid(
#     processed_config_dir: Path,
#     external_data_row_path: Path,
#     previous_hour_state_path: Path = None,
#     overrides: Main_Overrides = Main_Overrides(),
# ) -> Iterator[Model]:
#     """Get the model arguments from the input locations.

#     We pass all the full grid inputs and return an iterator.
#     Each pass of the iterator is a model setup for a single grid cell.
#     In the model setup are the config, external_state, initial_state, model_processes,
#     grid coords (x, y) and time boundaries (ts, te).


#     TODO: This is heavily focused towards grid runs. We should seperate grid and single cell runs.

#     Differences to full data run:
#     - config must be preprocessed so that we do not have to process every hour.
#     - we provide initial state input

#     Parameters
#     ----------
#     processed_config_dir : Path
#         path to processed config file (.json)
#     external_data_row_path : Path
#         path to data location (.csv) TODO: May need to be different format
#     previous_hour_state_path: Path = None,
#         path to state from previous hour (.json)
#     overrides : Main_Overrides
#         model overrides

#     Returns
#     -------
#     Model = namedtuple('Output', 'config,external_state,initial_state,model_processes,x,y,ts,te')

#     """

#     grid_coords = overrides.grid_coords

#     # get first config to use to get e_state_init_processes
#     config = config_loader_pickled(
#         f'{processed_config_dir}/{grid_coords[0][0]}_{grid_coords[0][1]}.config',
#     )

#     # External state iterator
#     external_state_data_iter = load_external_state(
#         external_data_row_path,
#         file_type=overrides.external_file_type,
#         coords=grid_coords,
#         variable_map=overrides.netcdf_variable_map,
#         preprocess_map=overrides.met_preprocess_map,
#         multi_file_data=overrides.multi_file_netcdf,
#         zero_year=config.Location.zero_year,
#         netcdf_chunks=overrides.netcdf_chunks,
#     )

#     e_state_init_processes = external_state_init_processes(0, 0, config.Met)

#     # Iterate over cells in external state data iterator
#     model_iterator = (
#         setup_model_partial_cell(
#             (xi, yi),
#             f'{processed_config_dir}/{xi}_{yi}.config',
#             external_state_data_cell,
#             previous_hour_state_path=previous_hour_state_path,
#             e_state_init_processes=e_state_init_processes,
#         )
#         for (xi, yi), external_state_data_cell in external_state_data_iter
#     )
#     return model_iterator


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
        raise InputDataError('Check shape of e_state_override netcdf vars')
    except Exception as e:
        print(coord)
        print(e_state_overrides_field_map)
        print(ds)


def load_additional_gridded_config_data(
    coords: List[Coord],
    e_state_overrides_file_path: Path,
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
    e_state_overrides_file_path : Path
        Path to external netcdf file
    e_state_overrides_field_map : dict
        Dictionary mapping of field in netcdf to field in Config_Shape

    Yields
    -------
    Iterator[Tuple[Coord, dict]]
        An iterator where each iteration corresponds to a coordinate.
        Each iteration contains a tuple of the coordinate and a map of config var to value.

    """
    ds = xr.open_dataset(e_state_overrides_file_path)

    for coord in coords:
        field_map = pull_config_vars_from_netcdf(ds, coord, e_state_overrides_field_map)
        yield (coord, field_map)


def get_grid_coords_from_dataarray(
    da,
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


def initialize_grid_configs(
    base_config: Config_Shape,
    processed_config_dir: Path,
    grid_coords: List[Tuple[int, int]],
    e_state_overrides_file_path: Path = None,
    e_state_overrides_field_map: dict = None,
    logger: Callable[[str, str], None] = print,
) -> Config_Shape:
    """Initialize and save a config file for each grid coord.

    This is for multi grid runs where we are using grid data such as sowing date
    to modify the config parameters.

    Parameters
    ----------
    base_config: Config_Shape
        Config common to all cells
    processed_config_dir : Path
        location to save processed configs
    grid_coords: List[Tuple[int, int]]
        grid coordinates
    e_state_overrides_file_path : Path, optional
        path to netcdf file with per coord overrides, by default None
    e_state_overrides_field_map : dict, optional
        map of netcdf fields to Config_Shape fields, by default None
    logger : Callable[[str, str], None], optional
        logger, by default print

    Returns
    -------
    Config_Shape
        the last generated config

    Raises
    ------
    AssertionError
        No configs found
    NotImplementedError
        Multi config not implemented

    """
    logger('===== Setting up cell config =====')

    # = We can supply a netcdf file to overrides per coord configs
    config_overrides = load_additional_gridded_config_data(
        grid_coords,
        e_state_overrides_file_path,
        e_state_overrides_field_map,
    ) if e_state_overrides_file_path is not None \
        else [(coord, {}) for coord in grid_coords]

    for (xi, yi), override_dict in config_overrides:
        config_cell = deepcopy(base_config)
        for k, v in override_dict.items():
            if v is None or math.isnan(v):
                raise ValueError(f'{k} is invalid for ({xi},{yi})')
            config_cell = rsetattr(config_cell, k, v, True)

        # TODO: implement overrides
        config_output_path = f"{processed_config_dir}/{xi}_{yi}.config"
        dump_config_to_file_binary(config_cell, config_output_path)


# @deprecated.deprecated()
# def setup_model_single_hour(
#     coord: Coord,
#     config_path: Path,
#     external_state_data: External_State_Shape,
#     previous_hour_state_path: Path,
#     e_state_init_processes: List[Process],
# ) -> Model:
#     """Convert each external state grid iteration into a model run.

#     # TODO: Chech multi hour data works too
#     Parameters
#     ----------
#     coord : Coord
#         The coordinate
#     config_path: Path
#         path of processed config file
#     external_state_data : External_State_Shape
#         The external state data for this grid cell
#     previous_hour_state_path : Path
#         Path to previous hour state for this grid cell
#     e_state_init_processes: List[Process],
#         List of processes to run on external state

#     Returns
#     -------
#     Model
#         A model setup ready to run

#     """
#     x, y = coord
#     config = config_loader_pickled(config_path)

#     # TODO: Can we optimize this? #OPTIMIZE
#     external_state = setup_external_state_simple(
#         external_state_data,
#         config,
#         e_state_init_processes,
#     )
#     # TODO: Can this handle when data wraps around end of year?
#     start_day = external_state.dd[0]
#     end_day = external_state.dd[-1]

#     previous_hour_state_path_tile = f"{previous_hour_state_path}/{x}_{y}.state"

#     # TODO: Can we optimize this. # OPTIMIZE
#     initial_state = model_state_loader_quick(previous_hour_state_path_tile)
#     initial_state.prev_hour = None
#     initial_state.prev_hour = deepcopy(initial_state)

#     hr = external_state_data.hr[0]
#     # TODO: Optimize getting model processes.
#     model_processes = [hourly_run_model_processes(config, hr)]
#     config_run = deepcopy(config)
#     config_run.Location.start_day = start_day
#     config_run.Location.end_day = end_day

#     return Model(
#         config=config_run,
#         external_state=external_state,
#         initial_state=initial_state,
#         model_processes=model_processes,
#         x=x,
#         y=y,
#     )


# @deprecated.deprecated()
# def setup_model_single_hour_grid(
#     processed_config_dir: Path,
#     external_data_row_path: Path,
#     previous_hour_state_path: Path = None,
#     overrides: Main_Overrides = Main_Overrides(),
# ) -> Iterator[Model]:
#     """Get the model arguments from the input locations.

#     We pass all the full grid inputs and return an iterator.
#     Each pass of the iterator is a model setup for a single grid cell.
#     In the model setup are the config, external_state, initial_state, model_processes,
#     grid coords (x, y) and time boundaries (ts, te).


#     TODO: This is heavily focused towards grid runs. We should seperate grid and single cell runs.

#     Differences to full data run:
#     - config must be preprocessed so that we do not have to process every hour.
#     - we provide initial state input

#     Parameters
#     ----------
#     processed_config_dir : Path
#         path to processed config file (.json)
#     external_data_row_path : Path
#         path to data location (.csv) TODO: May need to be different format
#     previous_hour_state_path: Path = None,
#         path to state from previous hour (.json)
#     overrides : Main_Overrides
#         model overrides

#     Returns
#     -------
#     Model = namedtuple('Output', 'config,external_state,initial_state,model_processes,x,y,ts,te')

#     """

#     grid_coords = overrides.grid_coords

#     # get first config to use to get e_state_init_processes
#     config = config_loader_pickled(
#         f'{processed_config_dir}/{grid_coords[0][0]}_{grid_coords[0][1]}.config',
#     )

#     # External state iterator
#     external_state_data_iter = load_external_state(
#         external_data_row_path,
#         file_type=overrides.external_file_type,
#         coords=grid_coords,
#         variable_map=overrides.netcdf_variable_map,
#         preprocess_map=overrides.met_preprocess_map,
#         multi_file_data=overrides.multi_file_netcdf,
#         zero_year=config.Location.zero_year,
#         netcdf_chunks=overrides.netcdf_chunks,
#     )

#     e_state_init_processes = external_state_init_processes(0, 0, config.Met)

#     # Iterate over cells in external state data iterator
#     model_iterator = (
#         setup_model_single_hour(
#             (xi, yi),
#             f'{processed_config_dir}/{xi}_{yi}.config',
#             external_state_data_cell,
#             previous_hour_state_path=previous_hour_state_path,
#             e_state_init_processes=e_state_init_processes,
#         )
#         for (xi, yi), external_state_data_cell in external_state_data_iter
#     )
#     return model_iterator


def get_input_files_list(
    dir: Path,
    multi_file_netcdf: bool = False,
    regex_multi_file_filter: str = '[0-9]{4}-[0-9]{2}',
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
        _description_
    multi_file_netcdf : bool, optional
        _description_, by default False
    regex_multi_file_filter : str, optional
        _description_, by default '[0-9]{4}-[0-9]{2}'

    Returns
    -------
    List[str]
        List of files or file group filters to run.

    """
    files_list = sorted(os.listdir(dir))
    if not multi_file_netcdf:
        return files_list
    else:
        files_list_grouped = sorted(set([
            getattr(re.search(regex_multi_file_filter, f), 'group', lambda: None)()
            for f in files_list
        ]))
        return files_list_grouped

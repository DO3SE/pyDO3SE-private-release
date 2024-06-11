import os
import shutil
import warnings
import numpy as np
import pandas as pd
import time
import queue
from pyDO3SE.optional_dependencies import xarray as xr
from multiprocessing import Process, Queue
from pathlib import Path
from typing import Any, Callable, List, Tuple, Generator, Optional
from datetime import datetime

from pyDO3SE.Config.config_loader import grid_config_loader
from pyDO3SE.setup_model import (
    Main_Overrides,
)
from pyDO3SE.Model_State.model_state_loader import (
    dump_state_to_file,
    dump_state_to_file_quick,
    load_current_cell_state,
    model_state_loader_quick,
)
from pyDO3SE.Pipelines.default_processes import get_row_processes_hourly
from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes
from pyDO3SE.util.logger import Logger, generate_run_notes, wrap_log
from pyDO3SE.version import config_version, version as model_version
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.External_State.external_state_loader import (
    get_date_bounds_from_ext_data,
    load_external_state,
    run_init_processes_on_e_state,
)
from pyDO3SE.External_State.External_State_Config import (
    FileTypes,
    EStateOptions,
)
from pyDO3SE.util.error_handling import Do3seRunError
from pyDO3SE.setup_model import (
    get_input_files_list,
    setup_external_state_simple,
)
from pyDO3SE.Grid_Model.setup_grid_model import (
    GridProjectPaths,
    GridRunFiles,
    GridRunPaths,
    ParallelArgs,
    create_grid_run_path_directories,
    get_grid_run_paths,
    get_grid_coords_from_file,
    load_grid_run_files,
)
from pyDO3SE.run_model import run_model, run_model_daily
from pyDO3SE.Output.process_outputs import (
    dump_output_to_file_netcdf_grid,
)
from pyDO3SE.Output.Output_Shape import OutputData

from pyDO3SE.Config import Config_Shape
from pyDO3SE.Model_State import Model_State_Shape


MISSING_OUTPUT_VALUE = -99999

Coords = Tuple[int, int]
OutputFields = List[str]
CellModel = Tuple[
    int,  # index
    Coords,
    Config_Shape,
    External_State_Shape,
    Model_State_Shape,
    List[Process],
]
# ModelGridCellOutput matrix of output with shape (row_count, len(output_fields))
ModelGridCellOutput = List[List[any]]


def run_cell(
    cell_model: CellModel,
    start_date: str,
    start_day: int,
    end_day: int,
    output_fields: List[str] = [],
    use_daily_loop: bool = False,
    # output_state_dir: str = None,
    row_count: int = None,
    cell_count: int = None,
    logger: Callable[[str], None] = Logger(),
    debug: bool = False,
) -> Tuple[Coords, List[any], Model_State_Shape, OutputFields]:
    """Run a single cell in grid model run.

    Can be ran in parallel.

    Parameters
    ----------
    cell_model : _type_
        Cell data combining all cell specific args
    start_date : str
        _description_
    start_day : int
        _description_
    end_day : int
        _description_
    output_fields : List[str], optional
        _description_, by default []
    use_daily_loop: bool
        If true then only preprocess a day worth of processes and run in a loop from start day to end day
    output_state_dir : str, optional
        _description_, by default None
    row_count: int
        Number of rows to run cell on
    cell_count: int
        Total number of cells for logging
    logger : Callable[[str], None], optional
        _description_, by default Logger()
    debug: bool, optional
        If on then run in debug mode

    Returns
    -------
    _type_
        _description_
    """
    # Each iteration runs a single grid cell over all the hours in the input data.
    # OPTIMIZE: This step taking a long time

    [i, (xi, yi), config_cell, external_state_data_cell,
        prev_hour_state_loaded, cell_processes] = cell_model

    logger(
        f'Preping coord ({xi}, {yi} ({i+1}/{cell_count}) for date: {start_date})\nExtracting Cell Model', verbose=True)

    logger('Model Extracted. Ready to run cell.')
    assert row_count <= len(
        external_state_data_cell.dd), \
        f"External state row count({len(external_state_data_cell.dd)}) is less than requested row count{row_count}"
    # Make sure start and end day set in config.
    config_cell.Location.start_day = start_day
    config_cell.Location.end_day = end_day

    # Run DO3SE model
    # TODO: Check if we should still use daily loop setup here or mapped runner
    model_runner = run_model_daily if use_daily_loop else run_model
    logger(
        f'Running coord ({xi}, {yi} for date: {start_date})\nExtracting Cell Model', verbose=True)
    try:
        final_state, output_logs = model_runner(
            initial_state=prev_hour_state_loaded,
            config=config_cell,
            external_state=external_state_data_cell,
            model_processes=cell_processes,
            start_index=0,
            DEBUG_MODE=debug,
        )
    except Exception as e:
        logger(f"Failed to run file at {xi}_{yi}")
        logger(f"Error: {e}")
        logger(f"Dumping state to tmp/prev_hour_state_loaded.json")
        os.makedirs('tmp', exist_ok=True)
        dump_state_to_file(prev_hour_state_loaded, 'tmp/prev_hour_state_loaded.json')
        raise Do3seRunError(f"Failed to run file at {xi}_{yi}", e)

    # Save state at each iteration
    # state_file_out_path = f"{output_state_dir}/{xi}_{yi}.state"
    # dump_state_to_file_quick(final_state, state_file_out_path)
    output_coords = [xi, yi]
    output_data = [[o.get(k, MISSING_OUTPUT_VALUE)
                    for o in output_logs] for k in output_fields]
    return output_coords, output_data, final_state, output_fields

    # # TODO: Can we save to netcdf here rather than storing all output in memory till end of runs?
    # for k in output_fields:
    #     # OPTIMIZE: Can we optimize this?
    #     # Output logs is List[dict[str, float]]
    #     full_output_data[k][xi, yi] = [o[k] for o in output_logs]
    # lat_data[xi, yi] = 0  # TODO: Get lat and lon data
    # lon_data[xi, yi] = 0  # TODO: Get lat and lon data

    # yield final_state, output_logs, cell_model


def run_in_sequence(
    func,
    input_args: List[List[Any]],
    func_kwargs: dict,
    **kwargs,
):
    return [func(*args, **func_kwargs) for args in input_args]


def run_in_parallel(
    func,
    total_population: int,
    input_args: List[List[Any]],
    func_kwargs: dict,
    parallel_args: ParallelArgs = ParallelArgs(),
) -> List[any]:
    """Run the model in parallel

    Parameters
    ----------
    func : Callable[[Parameters, ModelData], ModelOutput]
        model function.
    population: int
        Number of inputs to run in each iteration
    input_args : List[List[Any]]
        list of dictionaries of arguments for model run
    func_kwargs: dict
        Additional kwargs to pass to model func
    parallel_args: ParallelArgs
        Additional args for parallel processing

    Returns
    -------
    List[ModelOutput]
        The outputs of the Model for each run

    """
    Qu = None
    procs = []

    try:
        Qu = Queue(maxsize=parallel_args.MAX_PROCESSES)

        def _func(i, args, kwargs):
            # We use a middleware func to link the function outputs to the population index
            try:
                result = func(*args, **kwargs)
            except Exception as e:
                return i, e
            return i, result

        def q_wrap(q, args, kwargs):
            # Add the output of _func to the queue to be collected later
            q.put(_func(*args, kwargs))

        outputs = [None for _ in range(total_population)]

        # ==== for each run in population we create a process containing a queue writer.
        # NOTE: Args here is (i, model_args)
        next_index = 0
        running_processes = 0
        completed_processes = 0
        procs = []
        # TODO: Implement TIMEOUT
        while completed_processes < total_population:
            if running_processes > 0:
                time.sleep(parallel_args.SLEEP_TIME)
                try:
                    queue_get = Qu.get(timeout=.1)
                    if queue_get:
                        i, res = queue_get
                        procs[i].join()
                        procs[i].close()
                        if isinstance(res, Exception):
                            raise res
                        outputs[i] = res
                        running_processes -= 1
                        completed_processes += 1
                except queue.Empty:
                    pass

            if running_processes >= parallel_args.MAX_PROCESSES:
                time.sleep(0.1)
                # Already running max number of processes. Wait for some to finish
                continue
            if (next_index < total_population):
                args = next(input_args)
                p = Process(target=q_wrap, args=([Qu, [next_index, args], func_kwargs]))
                procs.append(p)
                p.start()
                next_index += 1
                running_processes += 1
    except Exception as e:
        try:
            for p in procs:
                try:
                    print(p.exitcode)
                    p.terminate()
                    p.join(.1)
                except:
                    p.kill()
        except:
            # TODO: Should warn about fail to teardown multi processes
            pass
        print(e)
        raise Exception("Model run multiprocessing failed: {}".format(str(e)))
    return outputs


Coords = Tuple[int, int]
OutputFields = List[str]


class ExternalStateIterable:
    start_day: int
    end_day: int
    start_date: str
    end_date: str
    row_count: int
    hours: List[int]
    time_string: str
    lat_data: np.ndarray
    lon_data: np.ndarray
    time_data: np.ndarray
    input_shape: List[int]
    external_state_preprocessed: Generator[External_State_Shape, None, None]

    def __init__(
        self,
        start_day: int,
        end_day: int,
        start_date: str,
        end_date: str,
        row_count: int,
        hours: List[int],
        time_string: str,
        lat_data: np.ndarray,
        lon_data: np.ndarray,
        time_data: np.ndarray,
        input_shape: List[int],
        external_state_preprocessed: Generator[External_State_Shape, None, None],
    ):
        self.start_day = start_day
        self.end_day = end_day
        self.start_date = start_date
        self.end_date = end_date
        self.row_count = row_count
        self.hours = hours
        self.time_string = time_string
        self.input_shape = input_shape
        self.lat_data = lat_data
        self.lon_data = lon_data
        self.time_data = time_data
        self.external_state_preprocessed = external_state_preprocessed

    def __iter__(self):
        return self.external_state_preprocessed

    def __next__(self):
        return next(self.external_state_preprocessed)


def setup_external_state_iter(
    cell_configs: List[Config_Shape],
    external_data_file_path: Path,
    external_state_options: EStateOptions,
    e_state_overrides_dataset: xr.Dataset,
    grid_coords: List[Tuple[int, int]],
    logger: Callable[[str], None] = Logger(),
    overrides: Main_Overrides = Main_Overrides(),
    debug: bool = False,
) -> ExternalStateIterable:
    logger(f"File path: {external_data_file_path}")

    # TODO: Should only load 1 cell.
    # OPTIMIZE: This step is taking over 1 min. Seems to be exponential based on number of rows loading
    # e_state_overrides_dataset = xr.open_dataset(project_paths.e_state_overrides_file_path)
    sample_e_state = next(load_external_state(
        external_data_file_path,
        [grid_coords[0]],
        external_state_options,
        logger=logger,
    ))
    sample_e_state_post_processed = run_init_processes_on_e_state(
        sample_e_state,
        config=cell_configs[0],
        init_processes=external_state_init_processes(config_met=cell_configs[0].Met),
        overrides=overrides,
    )
    [
        start_day,
        end_day,
        start_date,
        end_date,
        row_count,
        time_string,
        hours,
    ] = get_date_bounds_from_ext_data(
        sample_e_state_post_processed,
    )

    # This is an iterator over each grid cell.
    external_state_data_iter = load_external_state(
        external_data_file_path,
        grid_coords,
        external_state_options,
    )
    e_state_init_processes = (wrap_log(
        f"Generating e_state processes for config {i}",
        external_state_init_processes,
        coord,
        logger=logger
    )(config_met=config.Met)
        for i, (config, coord) in enumerate(zip(cell_configs, grid_coords)))

    # Process external state ready for model run.
    external_state_preprocessed = (setup_external_state_simple(
        external_state_data,
        config,
        e_state_processes_cell,
    ) for config, external_state_data, e_state_processes_cell
        in zip(cell_configs, external_state_data_iter, e_state_init_processes))

    logger(f"Row count: {row_count}, start_day: {start_day}, end_day: {end_day}")
    assert len(
        hours) == row_count, f"Hours length({len(hours)}) does not match row count ({row_count})"
    try:
        lat_data = e_state_overrides_dataset.lat_data.values
        lon_data = e_state_overrides_dataset.lon_data.values
        time_data = pd.date_range(start_date, periods=row_count, freq="1H")
        input_shape = [len(e_state_overrides_dataset.x), len(
            e_state_overrides_dataset.y)]

    except AttributeError as e:
        raise AttributeError(
            f"e_state_overrides_dataset must be an xarray.Dataset with lat_data, lon_data and x, y dimensions. Dataset has following keys: {list(e_state_overrides_dataset.keys())}") from e


    return ExternalStateIterable(
        start_day,
        end_day,
        start_date,
        end_date,
        row_count,
        hours,
        time_string,
        lat_data,
        lon_data,
        time_data,
        input_shape,
        external_state_preprocessed,
    )


def main_partial(
    cell_configs: List[Config_Shape],
    prev_hour_states_loaded: List[Model_State_Shape],
    grid_coords: List[Tuple[int, int]],
    external_states: ExternalStateIterable,
    output_fields: Optional[List[str]] = None,
    use_daily_loop: bool = False,
    parallel: bool = False,
    logger: Callable[[str], None] = Logger(),
    parallel_args: ParallelArgs = ParallelArgs(),
    overrides: Main_Overrides = Main_Overrides(),
    debug: bool = False,
) -> List[Tuple[Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]]:
    """Run the DO3SE model on partial data for a list of grid coordinates.

    Use case: we have netcdf data in hourly or monthly batches.

    We pass all the inputs to a setup func here that returns an iterator. Each pass of the
    iterator is a model setup for a single grid cell.
    In the model setup are the config, external_state, initial_state, model_processes and grid coords (x, y).

    Note: All inputs are provided as a list with length == len(grid_coords)

    Parameters
    ----------
    cell_configs : Config_Shape
        Loaded configs for this run
    prev_hour_states_loaded : List[Model_State_Shape]
        Loaded previous hour states for each grid cell for this run
    grid_coords: List[Tuple[int, int]]
        Grid coordinates to run
    external_states: ExternalStateIterable
        External state data for this run (already loaded and processed)
    output_fields : Optional[List[str]], optional
        A list of fields to output. Overrides config
    use_daily_loop: bool
        If true then only preprocess a day worth of processes and run in a loop from start day to end day
    parallel: bool
        If true then run cells in parallel
    logger: Callable[[str], None], optional
        log function, by default print
    parallel_args: ParallelArgs
        Args to pass to parallel
    overrides: Main_Overrides
        Overrides for main
    debug: bool, optional
        If on then run in debug mode

    Returns
    -------
    List[Tuple[Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]]
        List of outputs for each grid cell.
        length==len(grid_coords)
        Each item contains a tuple with:
            - Coords: (x, y) grid coordinates
            - ModelGridCellOutput: matrix of output with shape (row_count, len(output_fields))
            - Model_State_Shape: final state of model
            - OutputFields: list of output fields

    """
    logger("== Running Main = SETUP ==")

    model_processes = overrides.model_processes \
        or (get_row_processes_hourly(config, list(range(24)))
            for config in cell_configs) \
        if use_daily_loop \
        else (get_row_processes_hourly(config, external_states.hours)
              for config in cell_configs)

    # setup input args for each cell run
    input_args = (
        [cell] for cell in zip(
            range(len(grid_coords)), grid_coords, cell_configs, external_states, prev_hour_states_loaded, model_processes
        )
    )

    run_count = len(grid_coords)
    logger(
        f"== Running Main = Starting runs\nTotal cells: {len(grid_coords)}\nParallel: {parallel}")
    run_cells = run_in_parallel if parallel else run_in_sequence

    start_time = datetime.now()
    outputs_full_grid: List[Tuple[Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]] = run_cells(
        run_cell,
        total_population=run_count,
        input_args=input_args,
        func_kwargs=dict(
            start_date=external_states.start_date,
            start_day=external_states.start_day,
            end_day=external_states.end_day,
            output_fields=output_fields,
            use_daily_loop=use_daily_loop,
            # output_state_dir=output_data.output_state_dir,
            row_count=external_states.row_count,
            cell_count=len(grid_coords),
            logger=logger,
            debug=debug,
        ),
        parallel_args=parallel_args,
    )
    end_time = datetime.now()
    duration = end_time - start_time if end_time else 0
    logger(f"== Running Main = COMPLETE\nRuns took: ({duration}) ==")

    return outputs_full_grid
    # if output_data.output_state_dir:
    #     for coord, output_data, final_state in outputs_full_grid:
    #         # Save state at each iteration
    #         xi, yi = coord
    #         state_file_out_path = f"{output_data.output_state_dir}/{xi}_{yi}.state"
    #         dump_state_to_file_quick(final_state, state_file_out_path)

    # try:
    #     output_data = output_data.set_output_data_from_outputs_full_grid(
    #         outputs_full_grid
    #     )
    #     # for [xi, yi], outputs_cell in outputs_full_grid:
    #     #     for i, k in enumerate(output_fields):
    #     #         output_data.full_output_data[k][xi, yi] = outputs_cell[i]
    # except ValueError as e:
    #     if "could not convert string to float" in str(e):
    #         raise ValueError(f"Invalid value for key: {k}")
    #     raise e

    # # After all grid cells ran we save to netcdf
    # if output_directory:
    #     # TODO: Check time of output
    #     # time_data = pd.date_range(external_states.start_date,
    #     #                           periods=external_states.row_count, freq="1H")
    #     target_file_name = f'output_data_{external_states.time_string}.nc'
    #     dump_output_to_file_netcdf_grid(
    #         output_data=output_data,
    #         # X=grid_x_size, Y=grid_y_size, T=external_states.row_count,
    #         # lat_data=lat_data,
    #         # lon_data=lon_data,
    #         # time_data=time_data,
    #         target_path=f"{output_directory}/{target_file_name}",
    #     )

    # #   TODO: We have changed the output type here!
    # return output_data


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


def output_data_from_outputs_full_grid(
    outputs_full_grid: List[Tuple[Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]],
    output_shape: Tuple[int, int],
    output_fields: List[str],
    lat_data: np.ndarray,
    lon_data: np.ndarray,
    time_data: np.ndarray,
) -> OutputData:
    grid_x_size, grid_y_size, row_count = output_shape

    full_output_data = {
        k: np.full((grid_x_size, grid_y_size, row_count), None,
                   dtype=np.float64) for k in output_fields
    } if output_fields else None
    try:
        for [xi, yi], outputs_cell, final_state, output_fields in outputs_full_grid:
            for i, k in enumerate(output_fields):
                full_output_data[k][xi, yi] = outputs_cell[i]
        return OutputData(
            output_shape=output_shape,
            full_output_data=full_output_data,
            output_fields=output_fields,
            lat_data=lat_data,
            lon_data=lon_data,
            time_data=time_data,
        )
    except ValueError as e:
        if "could not convert string to float" in str(e):
            raise ValueError(f"Invalid value for key: {k}")
        raise e
    return


def save_model_iteration_output(
    outputs_full_grid: List[Tuple[Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]],
    # external_states_data: ExternalStateIterable,
    output_state_dir: Path,
    output_directory: Path,
    time_string: str,
    output_shape: Tuple[int, int, int],
    output_fields: str,
    lat_data: np.ndarray,
    lon_data: np.ndarray,
    time_data: np.ndarray,
):
    if output_state_dir:
        for coord, output_data, final_state, output_fields in outputs_full_grid:
            # Save state at each iteration
            xi, yi = coord
            state_file_out_path = f"{output_state_dir}/{xi}_{yi}.state"
            dump_state_to_file_quick(final_state, state_file_out_path)

    # After all grid cells ran we save to netcdf
    if output_directory:
        output_data = output_data_from_outputs_full_grid(
            outputs_full_grid,
            output_shape,
            output_fields,
            lat_data,
            lon_data,
            time_data,
        )

        # TODO: Check time of output
        # time_data = pd.date_range(external_states.start_date,
        #                           periods=external_states.row_count, freq="1H")
        target_file_name = f'output_data_{time_string}.nc'
        dump_output_to_file_netcdf_grid(
            full_output_data=output_data.full_output_data,
            output_shape=output_data.output_shape,
            lat_data=output_data.lat_data,
            lon_data=output_data.lon_data,
            time_data=output_data.time_data,
            output_fields=output_data.output_fields,
            target_path=f"{output_directory}/{target_file_name}",
        )


def main_grid_seq_per_config(
    project_paths: GridProjectPaths,
    run_paths: GridRunPaths,
    loaded_run_files: GridRunFiles,
    grid_coords: List[Tuple[int, int]],
    output_shape: Tuple[int, int],
    runnotes: str = '',
    output_fields: List[str] = None,
    seperate_live_state: bool = False,
    multi_file_netcdf: bool = False,
    regex_multi_file_filter: str = None,
    use_daily_loop: bool = False,
    parallel: bool = False,
    logger: Callable[[str], None] = print,
    parallel_args: ParallelArgs = ParallelArgs(),
    netcdf_loader_kwargs: dict = {},
    overrides: Main_Overrides = Main_Overrides(),
    return_outputs: bool = False,
    start_input_index: int = None,
    end_input_index: int = None,
    debug: bool = False,
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
        shape of the output grid
    runnotes : str, optional
        Additional notes to store with output, by default ''
    output_fields : List[str], optional
        Output state fields to save, by default None
    seperate_live_state : bool, optional
        If true then don't overwrite prev hour state with new state, by default False
    multi_file_netcdf : bool, optional
        If true use xr.open_mfdataset, by default False
    regex_multi_file_filter: str = '[0-9]{4}-[0-9]{2}'
        regex filter for grouping input files
    parallel: bool
        If true runs model in parallel
    use_daily_loop: bool
        If true then only preprocess a day worth of processes and run in a loop from start day to end day
    logger : Callable[[str], None], optional
        logger func, by default print
    parallel_args: ParallelArgs
        Args to pass to parallel
    netcdf_loader_kwargs: dict, optional
        kwargs to pass to xr.open_dataset or xr.open_mfdataset
    overrides: Main_Overrides
        Overrides for main
    return_outputs: bool
        If true then return outputs. Leaving as false may conserve memory.
    start_input_index: int
        If set then start at this input index
    end_input_index: int
        If set then end at this input index
    debug: bool, optional
        If on then run in debug mode

    """
    logger(f"== Running main_grid_seq_per_config:\nConfig: {run_paths.config_id} ===")
    cell_configs: List[Config_Shape] = grid_config_loader(
        run_paths.processed_configs_dir, grid_coords)
    for config in cell_configs:
        config.output.fields = (config.output.fields or []) + (output_fields or [])
        assert len(config.output.fields) > 0, "Must supply output fields in config or cli args!"
    zero_year = cell_configs[0].Location.zero_year
    assert zero_year, "Must supply zero year in configuration for grid runs!"
    errors = []
    e_state_overrides_dataset = xr.open_dataset(project_paths.e_state_overrides_file_path)
    input_files = get_input_files_list(
        project_paths.input_data_dir,
        multi_file_netcdf,
        regex_multi_file_filter,
    )
    if start_input_index is None:
        start_input_index = 0
    if end_input_index is None:
        end_input_index = len(input_files) - 1
    input_files = input_files[start_input_index:end_input_index + 1]
    total_files = len(input_files)
    logger(
        f"== Running Model: Total File: {total_files}. Start index: {start_input_index}. End index: {end_input_index} ===")

    external_state_options = EStateOptions(
        file_type=FileTypes.NETCDF,
        variable_map=loaded_run_files.variable_map,
        multi_file_data=multi_file_netcdf,
        preprocess_map=loaded_run_files.preprocess_map,
        zero_year=zero_year,
        netcdf_loader_kwargs=netcdf_loader_kwargs,
    )
    duration = None

    full_outputs: List[List[Tuple[Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]]] = [] if return_outputs else None

    try:
        start_time = datetime.now()

        # Iterate over each input file
        # NOTE: Assumes files sorted will be in chronological order.
        for i, f in enumerate(input_files):
            logger(f"Running file: {f} ({i+1}/{total_files})")
            state_out_path = run_paths.live_state_dir
            input_data_file = f"{project_paths.input_data_dir}/{f}" if not multi_file_netcdf else project_paths.input_data_dir
            previous_hour_state_path = run_paths.live_state_dir if not seperate_live_state else run_paths.prev_state_dir

            if regex_multi_file_filter:
                external_state_options = external_state_options._replace(
                    data_filter=f,
                )

            if seperate_live_state:
                shutil.rmtree(run_paths.prev_state_dir)
                os.rename(run_paths.live_state_dir, run_paths.prev_state_dir)
                os.makedirs(run_paths.live_state_dir, exist_ok=True)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # Iterates over each config for each cell in the grid
                external_states = setup_external_state_iter(
                    cell_configs,
                    external_data_file_path=input_data_file,
                    external_state_options=external_state_options,
                    e_state_overrides_dataset=e_state_overrides_dataset,
                    grid_coords=grid_coords,
                    logger=logger,
                    overrides=overrides,
                    debug=debug,
                )

                prev_hour_states_loaded = (load_current_cell_state(previous_hour_state_path, x, y)
                                           for x, y in grid_coords)
                outputs_full_grid = main_partial(
                    cell_configs=cell_configs,
                    prev_hour_states_loaded=prev_hour_states_loaded,
                    grid_coords=np.array(grid_coords),
                    external_states=external_states,
                    output_fields=output_fields,
                    parallel=parallel,
                    use_daily_loop=use_daily_loop,
                    logger=logger,
                    parallel_args=parallel_args,
                    overrides=overrides,
                    debug=debug,
                )
                save_model_iteration_output(
                    outputs_full_grid=outputs_full_grid,
                    output_state_dir=state_out_path or previous_hour_state_path,
                    output_directory=run_paths.output_data_dir,
                    time_string=external_states.time_string,
                    output_shape=[*output_shape, external_states.row_count],
                    output_fields=output_fields,
                    lat_data=external_states.lat_data,
                    lon_data=external_states.lon_data,
                    time_data=external_states.time_data,
                )
                if return_outputs:
                    full_outputs.append(outputs_full_grid)

            logger(f"==== Running file complete: {f} ({i+1}/{total_files}) =======")
            logger(f"==== Files saved in {run_paths.output_data_dir} =======")

        end_time = datetime.now()
        duration = end_time - start_time if end_time else 0

    except Exception as e:
        errors.append((f"Run Dir: {run_paths.config_run_dir} failed", e))

    if isinstance(logger, Logger) and logger.log_level >= 0:
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
    return full_outputs


def get_configs(dir: Path):
    return os.listdir(dir)


def main_grid_run(
    project_paths: GridProjectPaths,
    output_fields: List[str],
    multi_file_netcdf: bool = False,
    runnotes: str = '',
    logger: Logger = print,
    seperate_live_state: bool = False,
    save_final_state: bool = False,
    use_daily_loop: bool = False,
    parallel: bool = False,
    parallel_args: ParallelArgs = ParallelArgs(),
    netcdf_loader_kwargs: dict = {},
    regex_multi_file_filter: str = '',
    sample_size: int = 0,
    start_input_index: int = None,
    overrides: Main_Overrides = Main_Overrides(),
    debug: bool = False,
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
    logger : Logger, optional
        A Logger instance or print by default
    seperate_live_state : bool, optional
        If true then will save live state in seperate directory(This is more expensive), by default False
    use_daily_loop: bool
        If true then only preprocess a day worth of processes and run in a loop from start day to end day
    parallel: bool
        If true run in parallel
    parallel_args: ParallelArgs
        Args to pass to parallel
    netcdf_loader_kwargs: dict, optional
        kwargs to pass to xr.open_dataset or xr.open_mfdataset
    regex_multi_file_filter: str, optional
        regex filter for grouping input files
    sample_size: int, optional
        If greater than 0 then only runs up to sample size number of cells
    start_input_index: int
        If set then start at this input index
    overrides: Main_Overrides
        Overrides for main
    debug: bool, optional
        If on then run in debug mode

    """
    logger(f'== Running main_grid_run on:\n{project_paths.run_dir}=====')
    config_file_names = get_configs(project_paths.config_dir)

    logger(f'Found {len(config_file_names)} configs to run')

    for config_file_name in config_file_names:
        config_name = '.'.join(config_file_name.split('.')[:-1])
        run_paths = get_grid_run_paths(project_paths, config_name)
        loaded_run_files = load_grid_run_files(project_paths, run_paths)
        create_grid_run_path_directories(run_paths)
        grid_coords, *output_shape = get_grid_coords_from_file(
            run_paths.run_mask_path
        )
        grid_coords_sampled = grid_coords[0:sample_size]if sample_size else grid_coords

        main_grid_seq_per_config(
            project_paths=project_paths,
            run_paths=run_paths,
            loaded_run_files=loaded_run_files,
            grid_coords=grid_coords_sampled,
            output_shape=output_shape,
            runnotes=runnotes,
            output_fields=output_fields,
            seperate_live_state=seperate_live_state,
            multi_file_netcdf=multi_file_netcdf,
            use_daily_loop=use_daily_loop,
            parallel=parallel,
            logger=logger,
            parallel_args=parallel_args,
            netcdf_loader_kwargs=netcdf_loader_kwargs,
            regex_multi_file_filter=regex_multi_file_filter,
            start_input_index=start_input_index,
            overrides=overrides,
            debug=debug,
        )
        if save_final_state:
            for x, y in grid_coords_sampled:
                file_name = f"{x}_{y}"
                dump_state_to_file(model_state_loader_quick(
                    f"{run_paths.live_state_dir}/{file_name}.state"), f"{run_paths.final_state_dir}/{file_name}.json")

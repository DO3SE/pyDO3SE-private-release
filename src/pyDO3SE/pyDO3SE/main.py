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

from typing import get_type_hints
import os
import pandas as pd
from pathlib import Path
from dataclasses import asdict
from shutil import copyfile
from functools import partial
from multiprocessing import Pool
import warnings
import ntpath
from typing import Callable, List, NamedTuple, Tuple, Dict
from datetime import datetime

from pyDO3SE.lib import (
    RunFiles,
    ProjectPaths,
    RunPaths,
)
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Output.Output_Shape import default_output_fields
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.Model_State.model_state_loader import (
    model_state_loader,
)

from pyDO3SE.External_State.external_state_loader import (
    load_external_state,
)
from pyDO3SE.Output.run_comparisons import (
    create_comparison_graphs,
    get_output_files_from_run_list,
)
from pyDO3SE.util.logger import Logger
from pyDO3SE.setup_model import (
    Main_Overrides,
    setup_model,
    get_config_overrides,
)
from pyDO3SE.run_model import run_model_on_mapped_processes
from pyDO3SE.Output.process_outputs import (
    export_output,
    export_failed_run_output,
)
from pyDO3SE.Output.OutputConfig import (
    OutputOptions,
)

from pyDO3SE.Config import Config_Shape
from pyDO3SE.Model_State import Model_State_Shape


class MainOutput(NamedTuple):
    final_state: Model_State_Shape
    output_logs: List[List[any]]
    config_processed: Config_Shape
    initial_state: Model_State_Shape
    external_state: External_State_Shape


def load_run_files(
    project_paths: ProjectPaths,
    run_paths: RunPaths,
    logger: Callable[[str], None] = Logger(),
) -> RunFiles:
    logger("Loading run files")
    return RunFiles(
        config=config_loader(run_paths.config_path,
                             project_paths.base_config_path, 'json', logger=logger),
        state=project_paths.base_state_path and model_state_loader(
            project_paths.base_state_path, None, 'json', False) or Model_State_Shape(),
        per_input_config_overrides=pd.read_csv(project_paths.per_input_config_overrides)
        if project_paths.per_input_config_overrides and os.path.exists(
            project_paths.per_input_config_overrides) else None
    )


def get_project_paths(
    project_dir: str,
    config_dir: str = "configs",
    input_dir: str = "inputs",
    base_config_file_name: str = "base_config.json",
    runs_dir: str = "runs",
) -> ProjectPaths:
    return ProjectPaths(
        project_dir=project_dir,
        config_dir=f"{project_dir}/{config_dir}",
        input_data_dir=f"{project_dir}/{input_dir}",
        preprocess_map_path=f"{project_dir}/preprocess_map.json",
        base_config_path=f"{project_dir}/{base_config_file_name}",
        base_state_path=f"{project_dir}/base_state.json",
        runs_dir=f"{project_dir}/{runs_dir}",
        per_input_config_overrides=f"{project_dir}/per_input_config_overrides.csv"
    )


def get_run_paths(
    run_id: str,
    project_paths: ProjectPaths,
    config_id: str,
    input_file_id: str,
    output_path_format: str = "{runid}/{config}/{input}",
    run_name_format: str = "{config}_{input}_{runid}",
) -> RunPaths:
    run_dir = project_paths.runs_dir + "/" + output_path_format.format(
        config=config_id,
        runid=run_id,
        input=input_file_id,
    )
    run_name = run_name_format.format(
        config=config_id,
        runid=run_id,
        input=input_file_id,
    )

    return RunPaths(
        run_id=run_id,
        run_name=run_name,
        run_dir=run_dir,
        log_path=f"{run_dir}/run.log",
        config_path=f"{project_paths.config_dir}/{config_id}.json",
        input_data_file_path=f"{project_paths.input_data_dir}/{input_file_id}.csv",
        output_directory=f"{run_dir}/outputs",
        comparisons_dir=f"{run_dir}/comparisons",
        output_filename=f"{run_name}_output.csv",
        input_file_id=input_file_id,
        config_id=config_id,
    )


def create_run_path_directories(
    run_paths: RunPaths,
):
    os.makedirs(run_paths.run_dir, exist_ok=True) if run_paths.output_directory else None
    os.makedirs(run_paths.output_directory, exist_ok=True) if run_paths.output_directory else None


def main(
    project_paths: ProjectPaths,
    run_paths: RunPaths,
    fields_to_graph: List[str] = [],
    logger: Callable[[str], None] = Logger(),
    output_options: OutputOptions = OutputOptions(),
    *args,
    **kwargs,
) -> MainOutput:
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
    logger: Callable[[str], None], optional
        Log function, by default print
    output_options: OutputOptions
        Additional options around outputs

    *args, **kwargs
        Passed to Main_Overrides


    Returns
    -------
    Tuple[Model_State_Shape, List[List[Any]], Config_Shape, Model_State_Shape]
       final_state, output_logs, config, initial_state

    """
    overrides = Main_Overrides(*args, **kwargs)

    # MODEL SETUP
    loaded_run_files: RunFiles = load_run_files(
        project_paths=project_paths,
        run_paths=run_paths,
        logger=logger,
    )
    config_overrides: Dict[str, any] = get_config_overrides(
        run_paths.input_file_id, loaded_run_files.per_input_config_overrides)
    start_time_setup = datetime.now()

    external_state_data = next(load_external_state(
        run_paths.input_data_file_path,
        logger=logger,
    ))
    [
        config,
        external_state,
        initial_state,
        model_processes,
    ] = setup_model(
        config_in=loaded_run_files.config,
        state_in=loaded_run_files.state,
        external_state_in=external_state_data,
        run_dir=run_paths.run_dir,
        config_overrides=config_overrides,
        logger=logger,
        overrides=overrides,
    )

    time_taken_setup = datetime.now() - start_time_setup

    # MODEL RUN
    start_time = datetime.now()
    model_runner = run_model_on_mapped_processes
    try:
        final_state, output_logs = model_runner(
            initial_state,
            config,
            external_state,
            model_processes,
            DEBUG_MODE=overrides.debug,
        )
    except Exception as e:
        logger(f"Error running model. Attempting to dump logs to {run_paths.output_directory}.")
        try:
            export_failed_run_output(
                external_state,
                run_paths.output_directory,
                final_config=config,
                initial_state=initial_state,
                model_processes=model_processes,
                logger=logger,
                **asdict(output_options),
            )
        except Exception:
            logger("Error exporting failed run output")
        raise e
    time_taken = datetime.now() - start_time
    logger(f"Model run complete in {time_taken}")

    output_filename = run_paths.output_filename
    if run_paths.output_directory:
        logger("Exporting outputs")
        export_output(
            output_logs,
            final_state,
            external_state,
            run_paths.output_directory,
            final_config=config,
            initial_state=initial_state,
            output_filename=output_filename,
            model_processes=model_processes,
            fields_to_graph=fields_to_graph,
            observed_diurnal_path=project_paths.observed_diurnal_path,
            runid=run_paths.run_id,
            time_taken=time_taken,
            time_taken_setup=time_taken_setup,
            logger=logger,
            **asdict(output_options),
        )

    return MainOutput(final_state, output_logs, config, initial_state, external_state)


Args = get_type_hints(main)


def single(
    config_file: Path,
    data_file: Path,
    output_directory: Path,
    verbose: int = 0,
    runid: str = '',
    base_config_file: Path = None,
    per_input_config_overrides: Path = None,
    initial_state_path: Path = None,
    runnotes: str = '',
    plot_fields: str = None,
    observed_diurnal_data: Path = None,
    output_options: OutputOptions = OutputOptions(),
    overrides: Tuple[str, str] = (),
) -> MainOutput:
    """Run the model with provided config and input data.

    then runs the model for the number of days defined in the config
    outputting the final state
    """
    if not config_file:
        raise ValueError('Missing configfile')

    if not data_file:
        raise ValueError('Missing datafile')

    if not output_directory:
        raise ValueError('Missing outputdirectory')

    logger = Logger(log_level=verbose)
    logger("Running pyDO3SE")
    logger("Logger init")
    os.makedirs(output_directory, exist_ok=True)
    copyfile(config_file, f'{output_directory}/config.json')
    _plot_fields = plot_fields and plot_fields.split(',') or []

    project_paths = ProjectPaths(
        base_config_path=base_config_file,
        base_state_path=initial_state_path,
        observed_diurnal_path=observed_diurnal_data,
        per_input_config_overrides=per_input_config_overrides,
    )
    run_paths = RunPaths(
        run_name=runid,
        run_id=runid,
        config_path=config_file,
        input_data_file_path=data_file,
        output_directory=output_directory,
        output_filename=f"{runid}.csv",
        run_dir=f"{output_directory}/{runid}",
        input_file_id=ntpath.basename(data_file).replace('.csv', ''),
        config_id=ntpath.basename(config_file.replace('.json', ''))
    )
    create_run_path_directories(run_paths)

    _overrides = [o.split('=') for o in overrides]

    final_state, output_logs, config_processed, initial_state, external_state = main(
        project_paths=project_paths,
        run_paths=run_paths,
        fields_to_graph=_plot_fields,
        logger=logger,
        output_options=output_options,
        skip_state_init=initial_state_path is not None,
        debug=verbose > 0,
        **dict(_overrides),
    )
    return MainOutput(final_state, output_logs, config_processed, initial_state, external_state)


class RunOutput(NamedTuple):
    result: any
    error: Exception
    model_output: MainOutput


def run_from_args(args: Args, verbose) -> RunOutput:
    try:
        if not verbose:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out = main(
                    **args,
                )
        else:
            out = main(
                **args,
            )
        return RunOutput(0, None, out)
    except Exception as e:
        if (verbose):
            raise e
        config_file = args['run_paths'].config_path
        input_file = args['run_paths'].input_data_file_path
        warnings.warn(f"Failed to run config {config_file} on data {input_file}")
        return RunOutput(1, e, None)


def batch(
    project_directory: str = None,
    config_dir_override: str = None,
    input_dir_override: str = None,
    verbose: int = 0,
    run_comparisons=False,
    day_range: str = "0,365",
    output_fields: str = None,
    use_base_state: bool = False,
    compare_fields: str = None,
    default_paths: bool = True,
    output_options: OutputOptions = OutputOptions(),
    parallel: bool = True,
    runid: str = None,
    logger: Logger = Logger(),
):
    """Run the model with provided config and input data.

    A batch run allows you to run the model over multiple input(.csv) and config files(.json).

    Please note for the most up to date documentation on the batch cli inputs you
    should run `python -m pyDO3SE run batch --help`.
    This documentation page is intended to provide more indepth understanding of batch runs.

    There are a number of parameters and inputs specific to batch runs:


    Project Directory
    -----------------

    This is the folder that contains all the input and parameter files for your model run.

    Input files
    ^^^^^^^^^^^

    It should be set out as below::

    -- <PROJECT_DIR>
    ---- configs
    ------- config_1.json # can have any name
    ------- config_2.json # can have any name
    ------- ... additional configs ...
    ---- inputs
    ------- input_1.csv # can have any name
    ------- input_2.csv # can have any name
    ------- ... additional inputs...
    -----base_config.json
    -----base_state.json (Optional)
    -----per_input_config_overrides.csv (Optional)



    configs
        The model will run each config on each input file. The values in these
        config files will override the base config

    inputs
        This is where input met data should be located. For more info on input
        data requirements go to :ref:`input_data_setup`

    base_config.json
        This will be used for each run. Parts of the config will be
        overriden by the parameters in the configs directory. For info on setting up configs
        go to :ref:`config_setup`

    base_state.json
        This can be used to set the initial model state. Note that this
        should be taken from the model outputs and not created manually.

    per_input_config_overrides.csv
        This contains a list of config overrides for each
        input file. The headings should correspond to a config parameter e.g. the heading
        "Location.lat" will set that parameter for each input file. There must also be a "runid"
        heading which corresponds to the name of the input file without the .csv. E.g.
        "input1.csv" runid will be "input1".
        Please note that these overrides only apply after a config has been processed.
        Check the processed_configs in the outputs to see if this has been applied correctly


    The output directory will be populated as below::

    - --- runs
    - ----- <config_1>
    - ------- <input_1>
    - --------- config.json
    - --------- final_state.json
    - --------- notes.log
    - --------- processed_config.json
    - --------- external_data.csv
    - --------- <config>_<input>_out.csv
    - --------- <field_1>.png
    - --------- <field_2>.png
    - --------- ... additional field plots



    \b
    config.json
        A copy of the input config file
    final_state.json
        The final state of the model
    notes.log
        Any additional notes on the model run including model version;
        time taken to run etc.
    processed_config.json
        The actual model config that has been ran including any default,
        config processing etc.
    external_data.csv
        The processed external data that has been used for the model run.
        This includes any precalculated values.
    <config>_<input>_out.csv
        The output of the model run
    <field_1>.png
        Plots of fields selected

    """
    logger("Running pyDO3SE multirun")
    if verbose:
        logger("Running with verbose")

    project_paths = get_project_paths(
        project_directory) if default_paths else ProjectPaths()
    project_paths = project_paths._replace(
        config_dir=config_dir_override or project_paths.config_dir,
        input_data_dir=input_dir_override or project_paths.input_data_dir,
    )

    config_file_type = "json"
    config_files = [c for c in os.listdir(project_paths.config_dir) if len(c.split('.')) == 2 and c.split('.')[
        1] == config_file_type]
    config_files_str = "\n-".join(config_files)
    logger(f"Running following configs:\n{config_files_str}")

    input_files = [f for f in os.listdir(project_paths.input_data_dir) if f.split('.')[-1] == 'csv']
    output_fields_to_save = output_fields and output_fields.split(',') or None
    _compare_fields = default_output_fields if compare_fields == "_all" \
        else compare_fields and compare_fields.split(',') or []

    common_args = {}
    common_args['output_fields'] = output_fields_to_save
    common_args['output_options'] = output_options
    common_args['skip_state_init'] = use_base_state
    common_args['fields_to_graph'] = _compare_fields
    common_args['debug'] = verbose
    common_args['logger'] = logger

    run_list = [
        get_run_paths(
            runid,
            project_paths,
            config_file.split('.')[0],
            input_file.replace('.csv', ''),
        )
        for config_file in config_files for input_file in input_files
    ]

    for run_paths_i in run_list:
        create_run_path_directories(run_paths_i)

    # args to pass to main() for each run
    args_to_run = [
        {
            **common_args,
            "project_paths": project_paths,
            "run_paths": run_paths_i,
        } for run_paths_i in run_list
    ]

    start_time = datetime.now()

    # Run each file distributed
    results_info: List[Tuple[Args, RunOutput]] = []
    if parallel:
        with Pool(processes=8) as pool:
            results = pool.map(partial(run_from_args, verbose=verbose), args_to_run)
            results_info = zip(args_to_run, results)
    else:
        for args in args_to_run:
            result = run_from_args(args, verbose)
            results_info.append([args, result])

    runtime = datetime.now() - start_time

    logger(f"Complete! \nRuns took: {runtime}")
    logger(f"Output file located in {project_paths.project_dir}")

    if run_comparisons:
        if compare_fields is None:
            raise ValueError("Must supply compare_fields to run comparisons.")
        data_files = get_output_files_from_run_list(run_list)

        comparisons_dir = f"{project_paths.project_dir}/comparisons/{runid}"
        os.makedirs(comparisons_dir, exist_ok=True)
        output_plot_start_day, output_plot_end_day = [int(d) for d in day_range.split(',')] \
            if day_range is not None else [None, None]
        logger(
            f"Running comparisons on {compare_fields} for files: {[d[2] for d in data_files]}")
        create_comparison_graphs(
            data_files,
            _compare_fields,
            comparisons_dir,
            use_versioned_outdir=False,
            start_day=output_plot_start_day,
            end_day=output_plot_end_day,
            figsize=(20, 10),
        )

    # Check for errors
    failed_runs = [(args, r) for args, r in results_info if r.result == 1]
    if len(failed_runs) > 0:
        logger('Failed Runs (Use --verbose to see more details):', fg='red')
    for args, r in failed_runs:
        config_file = args['run_paths'].config_path
        input_file = args['run_paths'].input_data_file_path
        logger(f"{config_file} | {input_file}", fg='red')
        if verbose:
            logger(r.error)
            if verbose > 1:
                raise r.error

    return results_info

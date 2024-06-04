""" pyDO3SE python cli interface

Uses the click api to create a number of cli endpoints

Run `python pyDO3SE_cli.py --help` to get usage instructions

Available commands
---------------------

 - `run` - Run the model
 - `generate_config`
 - `demo` - runs

"""
from pathlib import Path
from functools import partial
from collections import namedtuple
import os
import warnings
from multiprocessing import Pool
from shutil import copyfile
from typing import Tuple
import click

from datetime import datetime
from pyDO3SE.Output.run_comparisons import create_comparison_graphs, get_output_files
from pyDO3SE import main


@click.command()
@click.option('--loglevel', default=0, help='Log level (0=off, 1=standard, 2=verbose)')
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--base_config_file', default=None, help='The base config file path', type=click.Path(exists=True))
@click.option('--runnotes', default='', help='Add a note to the run')
@click.option('--output_results_only', default=False, help='Only output results file, default False')
@click.option('--plot-fields', default='', help='A comma seperated list of output fields to graph')
@click.argument(
    'config_file',
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.argument(
    'data_file',
    type=click.Path(exists=True),
    # prompt='Enter data location (.csv)',
    # help='The location of the data file. Should be a csv file',
)
@click.argument(
    'output_directory',
    type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    'overrides', nargs=-1,
)
def single(
    config_file: Path,
    data_file: Path,
    output_directory: Path,
    loglevel: int = 0,
    runid: str = '',
    base_config_file: Path = None,
    runnotes: str = '',
    plot_fields: str = 'gsto_l',
    output_results_only: bool = False,
    overrides: Tuple[str, str] = (),
):
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

    click.echo("Running pyDO3SE")
    os.makedirs(output_directory, exist_ok=True)
    copyfile(config_file, f'{output_directory}/config.json')
    _plot_fields = plot_fields and plot_fields.split(',') or ['gsto_l', 'pody', 'lai']

    project_paths = main.ProjectPaths(
        base_config_path=base_config_file,
    )
    run_paths = main.RunPaths(
        config_path=config_file,
        input_data_file_path=data_file,
        output_directory=output_directory,
    )

    final_state, output_logs, config_processed, initial_state = main.main(
        project_paths=project_paths,
        run_paths=run_paths,
        runnotes=runnotes,
        fields_to_graph=_plot_fields,
        use_daily_loop=True,
        output_results_only=output_results_only,
        **dict(overrides),
    )

    click.echo("Complete!")
    click.echo(f"Output files located in {output_directory}")


RunOutput = namedtuple('RunOutput', 'result error')


def run_from_args(args, verbose):
    try:
        if not verbose:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                main.main(
                    **args,
                    # log_level=log_level,
                    output_results_only=True,
                )
        else:
            main.main(
                **args,
                # log_level=log_level,
                output_results_only=True,
            )
        return RunOutput(0, None)
    except Exception as e:
        config_file = args['run_paths'].config_path
        input_file = args['run_paths'].input_data_file_path
        warnings.warn(f"Failed to run config {config_file} on data {input_file}")
        return RunOutput(1, e)


@click.command()
@click.option('-v', '--verbose', count=True, default=False, help='Run full logging output')
@click.option('--runnotes', default='', help='Add a note to the run')
@click.option('--run_comparisons', default=False, help='Generate comparison graphs on outputs')
@click.option('--day_range', default='0,365', help='Set the start and end day')
@click.option('--output_fields', default=None, help='A comma seperated list of output fields to output')
@click.option('--compare_fields', default=None, help='A comma seperated list of output fields to graph')
@click.option('-c', '--default-paths/--custom-paths', default=True, help="If true use default paths")
@click.option(
    '--config-dir',
    required=False,
    type=click.Path(),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.option(
    '--input-dir',
    required=False,
    type=click.Path(),
    # prompt='Enter data location (.csv)',
    # help='The location of the data file. Should be a csv file',
)
@click.option(
    '--output-directory',
    required=False,
    type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    'project_directory',
    required=False,
    type=click.Path(),
)
def batch(
    project_directory: str = None,
    config_dir: str = None,
    input_dir: str = None,
    output_directory: str = None,
    verbose: int = 0,
    runnotes: str = '',
    config_file_type='json',
    run_comparisons=False,
    day_range: str = "0,365",
    output_fields: str = None,
    compare_fields: str = None,
    default_paths: bool = True,
):
    """Run the model with provided config and input data.

    then runs the model for the number of days defined in the config
    outputting the final state

    You can either provide the config, input and output directories seperately or
    define a project-directory. The project directory must have the following structure:

    \b
    -- <PROJECT_DIR>
    ---- configs
    ------- config_1.json # can have any name
    ------- config_2.json # can have any name
    ------- ... additional configs ...
    ---- inputs
    ------- input_1.json # can have any name
    ------- input_2.json # can have any name
    ------- ... additional inputs...
    -----base_config.json
    ---- outputs

    The output directory will be populated as below

    \b
    - --- outputs
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

    For help on creating the config run `python -m pyDO3SE.tools.cli generate-config --help`

    """
    click.echo("Running pyDO3SE multirun")
    if verbose:
        click.echo("Running with verbose")

    # try:
    #     assert config_dir or project_directory
    #     assert input_dir or project_directory
    #     assert output_directory or project_directory
    # except AssertionError:
    #     raise ValueError("Must supply config, input and output directories or project directory!")

    # config_dir = config_dir or f'{project_directory}/configs'
    # input_dir = input_dir or f'{project_directory}/inputs'
    # output_directory = output_directory or f'{project_directory}/outputs'

    project_paths = main.get_project_paths(
        project_directory) if default_paths else main.ProjectPaths()
    project_paths = project_paths._replace(
        config_dir=config_dir or project_paths.config_dir,
        input_data_dir=input_dir or project_paths.input_data_dir,
        output_directory=output_directory or project_paths.output_directory,
    )
    # run_paths = main.get_run_paths(project_paths, config_file, input_file, runid)

    start_day, end_day = [int(d) for d in day_range.split(',')] \
        if day_range is not None else [None, None]

    config_files = [c for c in os.listdir(project_paths.config_dir) if len(c.split('.')) == 2 and c.split('.')[
        1] == config_file_type]

    click.echo(f"Running following configs:\n{config_files}")
    input_files = [f for f in os.listdir(project_paths.input_data_dir) if f.split('.')[-1] == 'csv']

    # os.makedirs(output_directory, exist_ok=True)

    # NOT IMPLEMENTED YET FOR MULTIRUN
    output_fields_to_save = output_fields and output_fields.split(',') or None

    common_args = {}
    common_args['output_fields'] = output_fields_to_save
    run_paths = [
        main.get_run_paths(
            project_paths,
            config_file.split('.')[0],
            input_file.replace('.csv', ''),
            input_file.replace('.csv', ''),
        )
        for config_file in config_files for input_file in input_files
    ]

    for run_paths_i in run_paths:
        main.create_run_path_directories(run_paths_i)

    # args to pass to main() for each run
    args_to_run = [
        {
            **common_args,
            "project_paths": project_paths,
            "run_paths": run_paths_i,
            "runid": run_paths_i.run_id,
            'use_daily_loop': True,
        } for run_paths_i in run_paths
    ]

    start_time = datetime.now()

    # Run each file distributed
    results_info = []
    with Pool(processes=8) as pool:
        results = pool.map(partial(run_from_args, verbose=verbose), args_to_run)
        results_info = zip(args_to_run, results)
    endtime = datetime.now() - start_time

    log_notes = f"{runnotes}\n'Model took: {endtime} \n"
    click.echo(log_notes)
    click.echo("Complete!")
    click.echo(f"Output file located in {project_paths.output_directory}")
    if run_comparisons:
        if compare_fields is None:
            raise ValueError("Must supply compare_fields to run comparisons.")
        _compare_fields = compare_fields.split(',')
        click.echo(f"Running comparisons on {compare_fields}")
        data_files = get_output_files(project_paths.output_directory)
        create_comparison_graphs(
            data_files,
            project_paths.output_directory,
            _compare_fields,
            project_paths.output_directory,
            use_versioned_outdir=False,
            start_day=start_day,
            end_day=end_day,
        )

    # Check for errors
    failed_runs = [(args, r) for args, r in results_info if r.result == 1]
    if len(failed_runs) > 0:
        click.secho('Failed Runs (Use --verbose to see more details):', fg='red')
    for args, r in failed_runs:
        config_file = args['run_paths'].config_path
        input_file = args['run_paths'].input_data_file_path
        click.secho(f"{config_file} | {input_file}", fg='red')
        if verbose:
            click.echo(r.error)
            if verbose > 1:
                raise r.error


@click.command()
def generate_grid_run_setup(
    output_location: Path,
):
    """Generate the directory structure and base files for a grid run.

    """
    raise NotImplementedError("Generate grid run setup not implemented")


@click.command()
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--loglevel', default=0, help='Log level (0=off, 1=standard, 2=verbose)')
@click.argument(
    'project_dir',
    type=click.Path(exists=True),
)
def grid_run_init(
    project_dir: Path,
    runid: str = '',
    loglevel: int = 0,
):
    """Init the DO3SE grid model.

    """
    click.echo("Running Grid Init")
    project_paths = main.get_grid_project_paths(project_dir)
    main.init_all_grid_model_configs(
        project_paths,
        runid,
        log_level=loglevel,
    )
    click.echo("Running Grid Init - Complete")



@click.command()
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--loglevel', default=0, help='Log level (0=off, 1=standard, 2=verbose)')
@click.option('--init-model/--no-init-model', default=False, help='If true will initialize state and configs.')
@click.option('--runnotes', default='', help='Add a note to the run')
@click.option('--output_fields', default='pody', help='A comma seperated list of output fields to graph')
@click.option('--multi-file-netcdf/--single-file-netcdf', help='True if each variable is in its own netcdf file')
@click.argument(
    'project_dir',
    type=click.Path(exists=True),
)
def grid_run(
    project_dir: Path,
    runid: str = '',
    init_model: bool = True,
    output_fields: str = 'pody',
    runnotes: str = '',
    loglevel: int = 0,
    multi_file_netcdf: bool = False,
):
    """Run the DO3SE grid model.

    Inputs must netcdf files.
    """
    project_paths = main.get_grid_project_paths(project_dir)
    output_fields_to_graph = output_fields and output_fields.split(',') or ['gsto_l', 'pody', 'lai']
    if init_model:
        click.echo("Running grid init")
        main.init_all_grid_model_configs(
            project_paths,
            runid,
            log_level=loglevel,
        )
    click.echo("Running grid model")
    main.main_grid_run(
        project_paths,
        output_fields_to_graph,
        multi_file_netcdf,
        runid,
        runnotes,
        log_level=loglevel,
        save_final_state=True,
    )
    click.echo("==== Complete ===")

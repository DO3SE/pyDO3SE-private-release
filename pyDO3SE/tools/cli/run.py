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
from typing import Tuple
import click


from pyDO3SE import main
from pyDO3SE.Grid_Model import setup_grid_model, run_grid_model
from pyDO3SE.util.logger import Logger
from pyDO3SE.Output.Output_Shape import default_grid_output_fields
from pyDO3SE.Output.OutputConfig import OutputOptions, output_results_only_options


@click.command()
@click.option('-v', '--verbose', count=True, default=False, help='Run full logging output')
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--base_config_file', default=None, help='The base config file path', type=click.Path(exists=True))
@click.option('--per-input-config-overrides', default=None, help='Per input overrides csv file', type=click.Path(exists=True))
@click.option('--initial_state_path', default=None, help='The initial state file path', type=click.Path(exists=True))
@click.option('--runnotes', default='', help='Add a note to the run')
@click.option('--output_results_only', default=False, help='Only output results file, default False')
@click.option('--plot-fields', default='', help='A comma seperated list of output fields to graph')
@click.option('--observed-diurnal-data', default=None, help='Monthly average diurnal data', type=click.Path(exists=True))
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
    verbose: int = 0,
    runid: str = '',
    base_config_file: Path = None,
    per_input_config_overrides: Path = None,
    initial_state_path: Path = None,
    runnotes: str = '',
    plot_fields: str = 'gsto_l',
    observed_diurnal_data: Path = None,
    output_results_only: bool = False,
    overrides: Tuple[str, str] = (),
):
    """Run the model with provided config and input data.

    then runs the model for the number of days defined in the config
    outputting the final state
    """
    main.single(
        config_file=config_file,
        data_file=data_file,
        output_directory=output_directory,
        verbose=verbose,
        runid=runid,
        base_config_file=base_config_file,
        per_input_config_overrides=per_input_config_overrides,
        initial_state_path=initial_state_path,
        runnotes=runnotes,
        plot_fields=plot_fields,
        observed_diurnal_data=observed_diurnal_data,
        output_options=output_results_only and output_results_only_options() or OutputOptions(
            save_model_processes=True,
            save_model_processes_detailed=True,
        ),
        overrides=overrides,
    )
    click.echo("Complete!")
    click.echo(f"Output files located in {output_directory}")


@click.command()
@click.option('-v', '--verbose', count=True, default=False, help='Run full logging output')
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--run_comparisons/--skip-comparisons', default=False, help='Generate comparison graphs on outputs')
@click.option('--parallel/--single-core', default=True, help='If parallel then use multiprocessing to distribute run tasks.')
@click.option('--day_range', default='0,365', help='Set the start and end day')
@click.option('--output_results_only', default=False, help='Only output results file, default False')
@click.option('--output_fields', default=None, help='A comma seperated list of output fields to output')
@click.option('--compare_fields', default=None, help='A comma seperated list of output fields to graph. Use _all to plot all outputs.')
@click.option('--use-base-state/--no-base-state', default=False, help='If true will use base state input and skip state init')
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
@click.argument(
    'project_directory',
    required=False,
    type=click.Path(),
)
def batch(
    project_directory: str = None,
    config_dir: str = None,
    input_dir: str = None,
    verbose: int = 0,
    run_comparisons=False,
    day_range: str = "0,365",
    output_fields: str = None,
    use_base_state: bool = False,
    output_results_only: bool = False,
    compare_fields: str = None,
    default_paths: bool = True,
    parallel: bool = True,
    runid: str = None,
):
    """Run the model with provided config and input data.

    Full docs can be found at

    https://pydo3se-docs.onrender.com/source/pyDO3SE.main.html#pyDO3SE.main.batch

    """
    logger = Logger(verbose)
    main.batch(
        project_directory=project_directory,
        config_dir_override=config_dir,
        input_dir_override=input_dir,
        verbose=verbose,
        run_comparisons=run_comparisons,
        day_range=day_range,
        output_fields=output_fields,
        use_base_state=use_base_state,
        compare_fields=compare_fields,
        output_options=output_results_only and output_results_only_options() or OutputOptions(),
        default_paths=default_paths,
        parallel=parallel,
        runid=runid,
        logger=logger,
    )


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
    project_paths = setup_grid_model.get_grid_project_paths(project_dir, runid)
    init_log_path = f"{project_paths.run_dir}/init.log"
    logger_main = Logger(loglevel, init_log_path, set_as_default=True, write_mode='w')
    setup_grid_model.init_all_grid_model_configs(
        project_paths,
        logger_main,
    )
    logger_main.close()
    click.echo("Running Grid Init - Complete")


@click.command()
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--loglevel', default=0, help='Log level (0=off, 1=standard, 2=verbose)')
@click.option('--sample-size', default=0, help='If greater than 0 then only runs up to sample size')
@click.option('--start-input-index', default=None, help='If set will start at this input index')
@click.option('--init-model/--no-init-model', default=False, help='If true will initialize state and configs.')
@click.option('--runnotes', default='', help='Add a note to the run')
@click.option('--netcdf-regex-multi-file-filter', default=None, help='regex filter for grouping input files')
@click.option('--netcdf-concat-dim', default="Time", help='Passed to xarray.open_mfdataset')
@click.option('--output-fields', default='pody', help='A comma seperated list of output fields to graph')
@click.option('--multi-file-netcdf/--single-file-netcdf', help='True if each variable is in its own netcdf file')
@click.option('--parallel/--no-parallel', default=False, help='Use parallelized paramga')
@click.option('--e_state_overrides_field_map_path', default=None, help='use_per_config_e_state_overrides_field_map_path')
@click.option('--run_mask_path', default=None, help='use_per_config_run_mask_path')
@click.option('--max-processes', default=8, help='Maximum processes when running in parallel')
@click.option('--cpu-sleep-time', default=0.01, help='Cpu sleep time when running in parallel')
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
    sample_size: int = 0,
    loglevel: int = 0,
    start_input_index: int = None,
    multi_file_netcdf: bool = False,
    netcdf_concat_dim: str = "Time",
    parallel: bool = False,
    max_processes: int = 8,
    cpu_sleep_time: float = 0.01,
    netcdf_regex_multi_file_filter=None,
    e_state_overrides_field_map_path=None,
    run_mask_path=None,
):
    """Run the DO3SE grid model.

    Inputs must netcdf files.
    """
    project_paths = setup_grid_model.get_grid_project_paths(project_dir, runid)
    if not output_fields:
        raise ValueError("--output-fields must be set")
    _output_fields = default_grid_output_fields if output_fields == "_all" else output_fields.split(
        ',')

    logger_main = Logger(loglevel, project_paths.log_path, set_as_default=True,
                         write_mode='w', flush_per_log=loglevel >= 2)

    if init_model:
        click.echo("Running grid init")
        setup_grid_model.init_all_grid_model_configs(
            project_paths,
            logger_main,
            sample_size=sample_size,
            e_state_overrides_field_map_path=e_state_overrides_field_map_path,
            run_mask_path=run_mask_path,
        )
    click.echo("Running grid model")
    parallel_args = setup_grid_model.ParallelArgs(
        MAX_PROCESSES=max_processes,
        SLEEP_TIME=cpu_sleep_time,
    )

    netcdf_loader_kwargs = dict(
        parallel=parallel,
        concat_dim=netcdf_concat_dim,
        combine="nested",
    ) if multi_file_netcdf else {}
    try:
        run_grid_model.main_grid_run(
            project_paths=project_paths,
            output_fields=_output_fields,
            multi_file_netcdf=multi_file_netcdf,
            runnotes=runnotes,
            logger=logger_main,
            save_final_state=True,
            parallel=parallel,
            parallel_args=parallel_args,
            netcdf_loader_kwargs=netcdf_loader_kwargs,
            regex_multi_file_filter=netcdf_regex_multi_file_filter,
            start_input_index=start_input_index and int(start_input_index) or None,
            sample_size=sample_size,
            e_state_overrides_field_map_path=e_state_overrides_field_map_path,
            run_mask_path=run_mask_path,
        )
    except Exception as e:
        logger_main.log(f"Error running grid model: {e}")
        raise e
    finally:
        logger_main.close()
    click.echo("=== Running Grid Model Complete ===")

# """ pyDO3SE python cli interface

# Uses the click api to create a number of cli endpoints

# Run `python pyDO3SE_cli.py --help` to get usage instructions

# Available commands
# ---------------------

#  - `run` - Run the model
#  - `generate_config`
#  - `demo` - runs

# """
# from pathlib import Path
# from pprint import pprint
# from functools import partial
# import json
# from collections import namedtuple
# import os
# import warnings
# from multiprocessing import Pool
# from shutil import copyfile
# from typing import List, Tuple
# import click

# from datetime import datetime
# from data_helpers.list_helpers import flatten_list

# from do3se_phenology.plots import plot_phenology_from_config
# from pyDO3SE.setup_model import setup_config

# from pyDO3SE.version import config_version
# from pyDO3SE.Config.config_loader import config_loader
# from pyDO3SE.Config.generate_config import generate_config as generate_config_fn
# from pyDO3SE.Config.config_migration import Migrations
# from pyDO3SE.Config.Config_Shape import Config_Shape
# from pyDO3SE.Output.Output_Shape import output_fields_map
# from pyDO3SE.Output.run_comparisons import create_comparison_graphs, get_output_files
# from pyDO3SE.Output.Output_Shape import output_fields as model_output_fields
# from pyDO3SE.Output.process_outputs import dump_config_to_file_json, extract_final_results
# from pyDO3SE import version
# from pyDO3SE import main
# # from pyDO3SE.Pipelines.default_processes import (
# #     hourly_run_model_processes,
# # )


# @click.group()
# def cli():
#     """Main cli entrypoint."""
#     click.echo("Welcome to pyDO3SE")


# @cli.command()
# @click.option('--loglevel', default=0, help='Log level (0=off, 1=standard, 2=verbose)')
# @click.option('--runid', default='0', help='Set a run id for output naming')
# @click.option('--base_config_file', default=None, help='The base config file path', type=click.Path(exists=True))
# @click.option('--runnotes', default='', help='Add a note to the run')
# @click.option('--output_results_only', default=False, help='Only output results file, default False')
# @click.option('--plot-fields', default='', help='A comma seperated list of output fields to graph')
# @click.argument(
#     'config_file',
#     type=click.Path(exists=True),
#     # prompt='Enter config location(.json)',
#     # help='The location of the config file. Should be a json file',
# )
# @click.argument(
#     'data_file',
#     type=click.Path(exists=True),
#     # prompt='Enter data location (.csv)',
#     # help='The location of the data file. Should be a csv file',
# )
# @click.argument(
#     'output_directory',
#     type=click.Path(),
#     # prompt='Enter output directory',
#     # help='The location to save outputs',
# )
# @click.argument(
#     'overrides', nargs=-1,
# )
# def run(
#     config_file: Path,
#     data_file: Path,
#     output_directory: Path,
#     loglevel: int = 0,
#     runid: str = '',
#     base_config_file: Path = None,
#     runnotes: str = '',
#     plot_fields: str = 'gsto_l',
#     output_results_only: bool = False,
#     overrides: Tuple[str, str] = (),
# ):
#     """Run the model with provided config and input data.

#     then runs the model for the number of days defined in the config
#     outputting the final state
#     """
#     if not config_file:
#         raise ValueError('Missing configfile')

#     if not data_file:
#         raise ValueError('Missing datafile')

#     if not output_directory:
#         raise ValueError('Missing outputdirectory')

#     click.echo("Running pyDO3SE")
#     os.makedirs(output_directory, exist_ok=True)
#     copyfile(config_file, f'{output_directory}/config.json')
#     _plot_fields = plot_fields and plot_fields.split(',') or ['gsto_l', 'pody', 'lai']

#     final_state, output_logs, config_processed, initial_state = main.main(
#         config_file,
#         data_file,
#         output_directory,
#         fields_to_graph=_plot_fields,
#         runid=runid,
#         runnotes=runnotes,
#         log_level=loglevel,
#         base_config_file=base_config_file,
#         output_results_only=output_results_only,
#         **dict(overrides),
#     )

#     click.echo("Complete!")
#     click.echo(f"Output files located in {output_directory}")


# RunOutput = namedtuple('RunOutput', 'result error')


# def run_from_args(args, verbose):
#     log_level = 1 if verbose else 0
#     try:
#         if not verbose:
#             with warnings.catch_warnings():
#                 warnings.simplefilter("ignore")
#                 main.main(
#                     **args,
#                     log_level=log_level,
#                     output_results_only=True,
#                 )
#         else:
#             main.main(
#                 **args,
#                 log_level=log_level,
#                 output_results_only=True,
#             )
#         return RunOutput(0, None)
#     except Exception as e:
#         config_file = args['config_location']
#         input_file = args['data_location']
#         warnings.warn(f"Failed to run config {config_file} on data {input_file}")
#         return RunOutput(1, e)


# @cli.command()
# @click.option('--verbose', default=False, help='Run full logging output')
# @click.option('--runnotes', default='', help='Add a note to the run')
# @click.option('--base_config_file', default=None, help='The base config file path', type=click.Path(exists=True))
# @click.option('--run_comparisons', default=False, help='Generate comparison graphs on outputs')
# @click.option('--day_range', default='0,365', help='Set the start and end day')
# @click.option('--output_fields', default=None, help='A comma seperated list of output fields to output')
# @click.option('--compare_fields', default=None, help='A comma seperated list of output fields to graph')
# @click.option('--project_directory', default=None, help='Project directory', type=click.Path())
# @click.argument(
#     'config_dir',
#     required=False,
#     type=click.Path(),
#     # prompt='Enter config location(.json)',
#     # help='The location of the config file. Should be a json file',
# )
# @click.argument(
#     'input_dir',
#     required=False,
#     type=click.Path(),
#     # prompt='Enter data location (.csv)',
#     # help='The location of the data file. Should be a csv file',
# )
# @click.argument(
#     'output_directory',
#     required=False,
#     type=click.Path(),
#     # prompt='Enter output directory',
#     # help='The location to save outputs',
# )
# def multi_run(
#     config_dir: str = None,
#     input_dir: str = None,
#     output_directory: str = None,
#     project_directory: str = None,
#     verbose: bool = False,
#     base_config_file: Path = None,
#     runnotes: str = '',
#     config_file_type='json',
#     run_comparisons=False,
#     day_range: str = "0,365",
#     output_fields: str = None,
#     compare_fields: str = None,
# ):
#     """Run the model with provided config and input data.

#     then runs the model for the number of days defined in the config
#     outputting the final state

#     You can either provide the config, input and output directories seperately or
#     define a project-directory. The project directory must have the following structure:

#     \b
#     -- <PROJECT_DIR>
#     ---- configs
#     ------- config_1.json # can have any name
#     ------- config_2.json # can have any name
#     ------- ... additional configs ...
#     ---- inputs
#     ------- input_1.json # can have any name
#     ------- input_2.json # can have any name
#     ------- ... additional inputs...
#     -----base_config.json
#     ---- outputs

#     The output directory will be populated as below

#     \b
#     - --- outputs
#     - ----- <config_1>
#     - ------- <input_1>
#     - --------- config.json
#     - --------- final_state.json
#     - --------- notes.log
#     - --------- processed_config.json
#     - --------- external_data.csv
#     - --------- <config>_<input>_out.csv
#     - --------- <field_1>.png
#     - --------- <field_2>.png
#     - --------- ... additional field plots

#     \b
#     config.json
#         A copy of the input config file
#     final_state.json
#         The final state of the model
#     notes.log
#         Any additional notes on the model run including model version;
#         time taken to run etc.
#     processed_config.json
#         The actual model config that has been ran including any default,
#         config processing etc.
#     external_data.csv
#         The processed external data that has been used for the model run.
#         This includes any precalculated values.
#     <config>_<input>_out.csv
#         The output of the model run
#     <field_1>.png
#         Plots of fields selected

#     For help on creating the config run `python -m pyDO3SE.tools.cli generate-config --help`

#     """
#     click.echo("Running pyDO3SE multirun")
#     if verbose:
#         click.echo("Running with verbose")

#     try:
#         assert config_dir or project_directory
#         assert input_dir or project_directory
#         assert output_directory or project_directory
#     except AssertionError:
#         raise ValueError("Must supply config, input and output directories or project directory!")

#     config_dir = config_dir or f'{project_directory}/configs'
#     input_dir = input_dir or f'{project_directory}/inputs'
#     output_directory = output_directory or f'{project_directory}/outputs'

#     log_level = 1 if verbose else 0
#     start_day, end_day = [int(d) for d in day_range.split(',')] \
#         if day_range is not None else [None, None]

#     config_files = [c for c in os.listdir(config_dir) if len(c.split('.')) == 2 and c.split('.')[
#         1] == config_file_type]

#     click.echo(f"Running following configs:\n{config_files}")
#     input_files = [f for f in os.listdir(input_dir) if f.split('.')[-1] == 'csv']

#     os.makedirs(output_directory, exist_ok=True)

#     # NOT IMPLEMENTED YET FOR MULTIRUN
#     output_fields_to_save = output_fields and output_fields.split(',') or None

#     common_args = {}
#     common_args['output_fields'] = output_fields_to_save

#     # args to pass to main() for each run
#     args_to_run = [
#         {
#             **common_args,
#             "config_location": config_dir + '/' + config_file,
#             "data_location": input_dir + '/' + input_file,
#             "runid": input_file.replace('.csv', ''),
#             "output_directory": output_directory + '/' + config_file.split('.')[0],
#             "base_config_file": base_config_file,
#         } for config_file in config_files for input_file in input_files
#     ]

#     start_time = datetime.now()

#     # Run each file distributed
#     results_info = []
#     with Pool(processes=8) as pool:
#         results = pool.map(partial(run_from_args, verbose=verbose), args_to_run)
#         results_info = zip(args_to_run, results)
#     endtime = datetime.now() - start_time

#     log_notes = f"{runnotes}\n'Model took: {endtime} \n"
#     click.echo(log_notes)
#     click.echo("Complete!")
#     click.echo(f"Output file located in {output_directory}")
#     if run_comparisons:
#         if compare_fields is None:
#             raise ValueError("Must supply compare_fields to run comparisons.")
#         _compare_fields = compare_fields.split(',')
#         click.echo(f"Running comparisons on {compare_fields}")
#         data_files = get_output_files(output_directory)
#         create_comparison_graphs(
#             data_files,
#             output_directory,
#             _compare_fields,
#             output_directory,
#             use_versioned_outdir=False,
#             start_day=start_day,
#             end_day=end_day,
#         )

#     # Check for errors
#     failed_runs = [(args, r) for args, r in results_info if r.result == 1]
#     if len(failed_runs) > 0:
#         click.secho('Failed Runs (Use --verbose to see more details):', fg='red')
#     for args, r in failed_runs:

#         click.secho(f"{args['config_location']} | {args['data_location']}", fg='red')
#         if verbose:
#             click.echo(r.error)


# @cli.command()
# def generate_grid_run_setup(
#     output_location: Path,
# ):
#     """Generate the directory structure and base files for a grid run.

#     """
#     raise NotImplementedError("Generate grid run setup not implemented")


# @cli.command()
# @click.option('--runid', default='0', help='Set a run id for output naming')
# @click.option('--loglevel', default=0, help='Log level (0=off, 1=standard, 2=verbose)')
# @click.argument(
#     'project_dir',
#     type=click.Path(exists=True),
# )
# def grid_run_init(
#     project_dir: Path,
#     runid: str = '',
#     loglevel: int = 0,
# ):
#     """Init the DO3SE grid model.

#     """
#     main.init_all_grid_model_configs(
#         project_dir,
#         runid,
#         log_level=loglevel,
#     )


# @cli.command()
# @click.option('--runid', default='0', help='Set a run id for output naming')
# @click.option('--loglevel', default=0, help='Log level (0=off, 1=standard, 2=verbose)')
# @click.option('--init-model', default=True, help='If true will initialize state and configs.')
# @click.option('--runnotes', default='', help='Add a note to the run')
# @click.option('--output_fields', default='pody', help='A comma seperated list of output fields to graph')
# @click.argument(
#     'project_dir',
#     type=click.Path(exists=True),
# )
# def grid_run(
#     project_dir: Path,
#     runid: str = '',
#     init_model: bool = True,
#     output_fields: str = 'pody',
#     runnotes: str = '',
#     loglevel: int = 0,
# ):
#     """Run the DO3SE grid model.

#     Inputs must netcdf files.
#     """
#     output_fields_to_graph = output_fields and output_fields.split(',') or ['gsto_l', 'pody', 'lai']
#     if init_model:
#         main.init_all_grid_model_configs(
#             project_dir,
#             runid,
#             log_level=loglevel,
#         )
#     main.main_grid_seq(
#         project_dir,
#         runid,
#         output_fields_to_graph,
#         runnotes,
#         log_level=loglevel,
#     )


# @cli.command()
# @click.option("--gsto-method", default=None, help="Pick gsto method for config")
# @click.argument("out_location", type=click.Path())
# def generate_config(
#     out_location: Path,
#     gsto_method: str = None,
# ):
#     """Generate a template config file based on the response to some simple questions.



#     """
#     click.echo("Welcome to pyDO3SE config generator")
#     click.echo(gsto_method)
#     # TODO: Fix prompts
#     if gsto_method is None:
#         gsto_method = click.prompt('Please enter photosynthesis method', type=str)
#     click.echo(f"Config set to {gsto_method} defaults")
#     new_config = generate_config_fn(**{
#         "Land_Cover.parameters.0.gsto.method": gsto_method,
#     })
#     dump_config_to_file_json(new_config, out_location)


# @cli.command()
# @click.argument('config_location')
# def migrate_config(config_location: str):
#     """Migrate an old config file to the latest version."""
#     click.echo("Welcome to pyDO3SE config migrator!")
#     # Load old config
#     if os.path.isdir(config_location):
#         raise NotImplementedError("config_location must point to a config file not directory")
#     else:
#         with open(config_location) as config_file_data:
#             config = json.load(config_file_data)
#             input_version = config.get('VERSION', 0)
#             migrated_config = Migrations.run_migrations(config, input_version)
#             out_dir = os.path.dirname(config_location)
#             out_file_name = os.path.splitext(os.path.basename(config_location))[
#                 0] + f"_{config_version}.json"
#             with open(f"{out_dir}/{out_file_name}", 'w') as out_location:
#                 json.dump(migrated_config, out_location, indent=4)


# @cli.command()
# @click.option('--verbose', default=False, help='Run full logging output')
# @click.option('--runid', default='0', help='Set a run id for output naming')
# @click.option('--runnotes', default='', help='Add a note to the run')
# def demo(verbose=False, runid=0, runnotes='Demo run'):
#     """Run pyDO3SE in demo mode with example inputs and config."""
#     click.echo("Running pyDO3SE in demo mode")

#     demo_config = "./examples/demo/configs/demo_config.json"
#     demo_data = "./examples/demo/inputs/demo_data.csv"
#     demo_output_dir = f"./examples/demo/outputs/{runid or version}"
#     log_level = 1 if verbose else 0

#     # Create output dir
#     os.makedirs(demo_output_dir, exist_ok=True)

#     # Copy config to this directory
#     copyfile(demo_config, f'{demo_output_dir}/config.json')

#     output_fields_to_graph = ['gsto_l', 'A_n']

#     final_state, output_logs, config_processed, initial_state = main.main(
#         demo_config,
#         demo_data,
#         demo_output_dir,
#         output_fields=output_fields_to_graph,
#         runid=runid,
#         runnotes=runnotes,
#         log_level=log_level,
#         start_day=1,
#         end_day=365,
#         debug=True,
#     )

#     click.echo("Complete!")
#     click.echo(f"Output file located in {demo_output_dir}")


# @cli.command()
# @click.option('--verbose', default=False, help='Run full logging output')
# def demo_multirun(verbose=False):
#     """Run pyDO3SE in demo mode with example inputs and config."""
#     # TODO: Make sure this works
#     click.echo("Running pyDO3SE in multirun demo mode")
#     demo_config_directory = "./examples/spanish_wheat/configs"
#     demo_data_directory = "./examples/spanish_wheat/data/"
#     demo_output_dir = "./examples/spanish_wheat/outputs"
#     log_level = 1 if verbose else 0
#     click.echo(f'Started at: {datetime.now()}')
#     config_files = (f for f in os.listdir(demo_config_directory) if f[-5:] == '.json')
#     data_files = (f for f in os.listdir(demo_data_directory) if f[-4:] == '.csv')
#     for config in config_files:
#         for data in data_files:
#             click.echo(f'Running config: {config} on datafile: {data}')
#             main.main(config, data, demo_output_dir, log_level=log_level)
#     click.echo(f'Finished at: {datetime.now()}')

#     click.echo("Complete!")
#     click.echo(f"Output file located in {demo_output_dir}")


# @click.option('--day_range', default='0,365', help='Set the start and end day')
# @click.option('--filter', default='', help='Filter directory')
# @click.argument(
#     'fields_to_graph',
#     # prompt='Enter fields to compare',
#     # help='The fields to run comparisons',
#     nargs=-1,
# )
# @click.argument(
#     'output_directory',
#     # prompt='Enter output directory',
#     # help='The location to save outputs',
#     nargs=1,
# )
# @click.argument(
#     'input_directory',
#     # prompt='Enter model outputs directory',
#     # help='The location of model outputs',
#     nargs=1,
# )
# @cli.command()
# def compare_outputs(
#     input_directory: Path,
#     output_directory: Path,
#     fields_to_graph: List[str] = [''],
#     day_range: str = None,  # "0,365",
#     filter: str = '',
# ):
#     """Create comparison graphs of all files in a directory.

#     Each output to compare should be in a subdirectory. The name of the directory will
#     be the name of the series.
#     """
#     click.echo(f"Comparing outputs from {input_directory}")
#     start_day, end_day = [int(d) for d in day_range.split(
#         ',')] if day_range is not None else [None, None]
#     click.echo(f"Day range: {start_day} to {end_day}")
#     _fields_to_graph = fields_to_graph if len(fields_to_graph) > 0 else [
#         o.id for o in model_output_fields]
#     click.echo(f"Comparing fields {_fields_to_graph}")

#     data_files = get_output_files(input_directory, filter)
#     click.echo(f"Running Comparisons on outputs")
#     create_comparison_graphs(
#         data_files,
#         input_directory,
#         _fields_to_graph,
#         output_directory,
#         use_versioned_outdir=False,
#         start_day=start_day,
#         end_day=end_day,
#         use_versioned_outfile=False,
#     )

#     click.echo("Complete")


# @click.argument(
#     'field_name',
#     # prompt='field name to extract',
#     nargs=1,
# )
# @click.argument(
#     'output_directory',
#     # prompt='Enter output directory',
#     # help='The location to save outputs',
#     nargs=1,
# )
# @click.argument(
#     'input_directory',
#     # prompt='Enter model outputs directory',
#     # help='The location of model outputs',
#     nargs=1,
# )
# @cli.command()
# def multi_run_final_value(
#     input_directory: Path,
#     output_directory: Path,
#     field_name: str,
# ):
#     """Get the final value from each run (i.e. pody)"""
#     extract_final_results(input_directory, output_directory, field_name)


# @click.option('--input-data-file', type=click.Path(exists=True), default=None, help='Input data csv file')
# @click.option('--plot-dd', default=False, help='Plot day data')
# @click.argument(
#     'output_directory',
#     # type=click.Path(),
#     # prompt='Enter output directory',
#     # help='The location to save outputs',
# )
# @click.argument(
#     'config_file',
#     type=click.Path(exists=True),
#     # prompt='Enter config location(.json)',
#     # help='The location of the config file. Should be a json file',
# )
# @cli.command()
# def plot_phenology(
#     config_file: Path,
#     output_directory: Path,
#     input_data_file: Path = None,
#     plot_dd: bool = False,
# ):
#     """Plot the phenology from config file."""
#     print(config_file)
#     config: Config_Shape = config_loader(config_file, 'json')
#     os.makedirs(output_directory, exist_ok=True)
#     day_count = config.Location.end_day - config.Location.start_day

#     if plot_dd:
#         raise NotImplementedError('Plotting day data is not implemented')
#         assert input_data_file is not None

#     plot_phenology_from_config(
#         config.Land_Cover.parameters[0].phenology,
#         config.Land_Cover.phenology_options,
#         nP=config.Land_Cover.nP,
#         output_location=f"{output_directory}/phenology.png",
#         # TODO: Add external data input
#         day_count=day_count,
#         plot_dd=plot_dd,
#     )


# @cli.command()
# def available_outputs():
#     """Print available outputs for output data and graphs.

#     NOTE: This is not fully implemented and only a subset are
#     outputed to the output file.
#     """
#     print("WARNING: Only a subset of these fields are actually available currently")
#     pprint(output_fields_map)


# if __name__ == "__main__":

#     cli()


# @cli.command()
# @click.option('--out-location', default=None, type=click.Path(exists=True))
# @click.option('--base-config-path', default=None, type=click.Path(exists=False))
# @click.argument('config-location', type=click.Path(exists=True))
# def output_process_list(
#     config_location: Path,
#     out_location: Path = None,
#     base_config_path: Path = None,
# ):
#     config = config_loader(config_location, base_config_path, 'json')

#     # hourly_processes_out = hourly_processes(config, 0)
#     # daily_processes_out = daily_process_list(config)
#     # full_model_processes_out = full_model_processes(config, 1, 1)

#     full_model_processes_out = hourly_run_model_processes(config, 0)

#     flattened_process_comments = [
#         p.comment or p.func.__name__ for p in flatten_list(full_model_processes_out)]

#     if out_location:
#         with open(out_location, 'w') as f:
#             f.write('\n'.join(flattened_process_comments))
#     else:
#         print('\n'.join(flattened_process_comments))

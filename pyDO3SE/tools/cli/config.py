from pathlib import Path
from pprint import pprint
import json
import os
import click

from data_helpers.list_helpers import flatten_list


from do3se_phenology.plots import plot_phenology_from_config
from pyDO3SE.version import config_version
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Config.generate_config import generate_config
from pyDO3SE.Config.config_migration import Migrations
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.Output.process_outputs import dump_config_to_file_json
from pyDO3SE.Pipelines.default_processes import (
    full_model_processes,
)


@click.command()
@click.option("--gsto-method", default=None, help="Pick gsto method for config")
@click.argument("out_location", type=click.Path())
def generate(
    out_location: Path,
    gsto_method: str = None,
):
    """Generate a template config file based on the response to some simple questions.



    """
    click.echo("Welcome to pyDO3SE config generator")
    click.echo(gsto_method)
    # TODO: Fix prompts
    if gsto_method is None:
        gsto_method = click.prompt('Please enter photosynthesis method', type=str)
    click.echo(f"Config set to {gsto_method} defaults")
    new_config = generate_config(**{
        "Land_Cover.parameters.0.gsto.method": gsto_method,
    })
    dump_config_to_file_json(new_config, out_location)


@click.command()
@click.argument('config_location')
def migrate(config_location: str):
    """Migrate an old config file to the latest version."""
    click.echo("Welcome to pyDO3SE config migrator!")
    # Load old config
    if os.path.isdir(config_location):
        raise NotImplementedError("config_location must point to a config file not directory")
    else:
        with open(config_location) as config_file_data:
            config = json.load(config_file_data)
            input_version = config.get('VERSION', 0)
            migrated_config = Migrations.run_migrations(config, input_version)
            out_dir = os.path.dirname(config_location)
            out_file_name = os.path.splitext(os.path.basename(config_location))[
                0] + f"_{config_version}.json"
            with open(f"{out_dir}/{out_file_name}", 'w') as out_location:
                json.dump(migrated_config, out_location, indent=4)


@click.command()
def available_outputs():
    """Print available outputs for output data and graphs.

    NOTE: This is not fully implemented and only a subset are
    outputed to the output file.
    """
    print("WARNING: Only a subset of these fields are actually available currently")
    pprint(output_fields_map)


@click.command()
@click.option('--detailed', default=False)
@click.option('--out-location', default=None, type=click.Path(exists=True))
@click.option('--base-config-path', default=None, type=click.Path(exists=False))
@click.argument('config-location', type=click.Path(exists=True))
def output_process_list(
    config_location: Path,
    out_location: Path = None,
    base_config_path: Path = None,
    detailed: bool = False,
):
    """Output all the processes ran with current config.

    Parameters
    ----------
    config_location : Path
        location of config file
    out_location : Path, optional
        Location to save output, by default None
    base_config_path : Path, optional
        Location of base config file, by default None
    detailed: bool
        If true will print detailed info on each process

    """
    config = config_loader(config_location, base_config_path, 'json')

    hours= list(range(24))
    full_model_processes_out = full_model_processes(config, hours)
    flattened_processes = [p.human() for p in flatten_list(full_model_processes_out)]
    if detailed:
        flattened_process_comments = [
            "\n".join(filter(lambda a: a, [
                # "________________________________________",
                " ",
                (f"# " if "=====" in p.get("comment") else f"## " if "==" in p.get("comment") else "### ") +
                f"{p.get('comment') or p.get('func')}",
                "\tconfig_inputs: " + str(p.get('config_inputs')) if p.get('config_inputs') else "",
                "\texternal_state_inputs: " + str(p.get('external_state_inputs')) if p.get('external_state_inputs') else "",
                "\tstate_inputs: " + str(p.get('state_inputs')) if p.get('state_inputs') else "",
                "\tstate_outputs: " + str(p.get('state_outputs')) if p.get('state_outputs') else "",
            ]))
            for p in flattened_processes]
    else:
        flattened_process_comments = [
            p.get('comment') or p.get('func') for p in flattened_processes]

    if out_location:
        with open(out_location, 'w') as f:
            f.write('\n'.join(flattened_process_comments))
    else:
        print('\n'.join(flattened_process_comments))


@click.option('--input-data-file', type=click.Path(exists=True), default=None, help='Input data csv file')
@click.option('--plot-dd', default=False, help='Plot day data')
@click.argument(
    'output_directory',
    # type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    'config_file',
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.command()
def plot_phenology(
    config_file: Path,
    output_directory: Path,
    input_data_file: Path = None,
    plot_dd: bool = False,
):
    """Plot the phenology from config file."""
    print(config_file)
    config: Config_Shape = config_loader(config_file, 'json')
    os.makedirs(output_directory, exist_ok=True)
    day_count = config.Location.end_day - config.Location.start_day

    if plot_dd:
        raise NotImplementedError('Plotting day data is not implemented')
        assert input_data_file is not None

    plot_phenology_from_config(
        config.Land_Cover.parameters[0].phenology,
        config.Land_Cover.phenology_options,
        nP=config.Land_Cover.nP,
        output_location=f"{output_directory}/phenology.png",
        # TODO: Add external data input
        day_count=day_count,
        plot_dd=plot_dd,
    )

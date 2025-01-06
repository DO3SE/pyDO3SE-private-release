from pathlib import Path
from pprint import pprint
import json
import os
import click

from do3se_phenology.plots import plot_phenology_from_config
from pyDO3SE.version import config_version
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Config.generate_config import generate_config
from pyDO3SE.Config.config_migration import Migrations
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.Output.process_outputs import dump_config_to_file_json, dump_model_processes_info_to_string
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
    dump_config_to_file_json()(new_config, out_location)


def migrate_config(config_dir: Path, filename: Path, override):
    config_location = f"{config_dir}/{filename}"
    with open(config_location) as config_file_data:
        print(f"Running migrations on \"{config_location}/{filename}\"")
        config = json.load(config_file_data)
        input_version = config.get('VERSION', 0)
        try:
            migrated_config = Migrations.run_migrations(config, input_version)
        except Exception as e:
            raise Exception(f"Failed to migrate {filename}. {e}")
        out_dir = config_dir  # os.path.dirname(config_location)
        out_file_name = os.path.splitext(filename)[0]
        out_file_name = f"{out_file_name}.json" if override else f"{out_file_name}_{config_version}.json"
        with open(f"{out_dir}/{out_file_name}", 'w') as out_location:
            json.dump(migrated_config, out_location, indent=4)


@click.command()
@click.option("--override/--no-override", default=False, help="If true then will override the input config. Otherwise a copy is created.")
@click.argument('config_location')
def migrate(config_location: str, override: bool):
    """Migrate an old config file to the latest version.

    CONFIG_LOCATION can be a file or directory.

    """
    click.echo("Welcome to pyDO3SE config migrator!")
    # Load old config
    if os.path.isdir(config_location):
        click.echo(f"Migrating all files in '{config_location}'")
        config_files = [f for f in os.listdir(config_location) if ".json" in f]
        for f in config_files:
            migrate_config(config_location, f, override)
    else:
        config_dir = os.path.dirname(config_location)
        click.echo(f"Migrating the following file: '{config_location}' in dir '{config_dir}'")
        migrate_config(config_dir, os.path.basename(config_location), override)


@click.command()
def available_outputs():
    """Print available outputs for output data and graphs.

    NOTE: This is not fully implemented and only a subset are
    outputed to the output file.
    """
    print("WARNING: Only a subset of these fields are actually available currently")
    pprint(output_fields_map)


@click.command()
@click.option('--detailed/--short', default=False)
@click.option('--allow-errors/--strict', default=False)
@click.option('--out-location', default=None, type=click.Path(exists=False))
@click.option('--base-config-path', default=None, type=click.Path(exists=False))
@click.argument('config-location', type=click.Path(exists=True))
def output_process_list(
    config_location: Path,
    out_location: Path = None,
    base_config_path: Path = None,
    detailed: bool = False,
    allow_errors: bool = False,
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

    hours = list(range(24))
    full_model_processes_out = full_model_processes(config, hours)
    flattened_process_comments = dump_model_processes_info_to_string(
        full_model_processes_out, detailed, allow_errors)

    if out_location:
        with open(out_location, 'w') as f:
            f.write('\n'.join(flattened_process_comments))
    else:
        print('\n'.join(flattened_process_comments))


@click.option('--input-data-file', type=click.Path(exists=True), default=None, help='Input data csv file')
@click.option('--base-config-file', default=None, help='The base config file path', type=click.Path(exists=True))
@click.option('--plot-dd', default=False, help='Plot day data')
@click.argument(
    'output-directory',
    # type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    'config-file',
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.command()
def plot_phenology(
    config_file: Path,
    output_directory: Path,
    input_data_file: Path = None,
    base_config_file: Path = None,
    plot_dd: bool = False,
    day_count: int = 365,
):
    """Plot the phenology from config file."""
    config: Config_Shape = config_loader(config_file, base_config_file, 'json')
    os.makedirs(output_directory, exist_ok=True)

    day_count = (config.Location.end_day or 0) - (config.Location.start_day or 0) if config.Location.start_day \
        is not None and config.Location.end_day is not None else day_count

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




@click.option('--input-data-file', type=click.Path(exists=True), default=None, help='Input data csv file')
@click.option('--base-config-file', default=None, help='The base config file path', type=click.Path(exists=True))
@click.argument(
    'output-directory',
    # type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    'config-file',
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.command()
def process_config(
    config_file: Path,
    output_directory: Path,
    input_data_file: Path = None,
    base_config_file: Path = None,
):
    """Process the config by combining the config and base config.

    There are also some setup processes applied to the config such as setting up
    the phenology parameters.

    Note this does not include parameters defined by the external state.

    """

    from pyDO3SE.setup_model import setup_config, get_config_overrides
    from pyDO3SE.Output.process_outputs import dump_config_to_file_json
    config_in=config_loader(config_file, base_config_file, 'json')

    external_state_in = None # TODO: Implement external state if provided
    input_file_id = None # TODO: Implement input file id if provided
    per_input_config_overrides = {} # TODO: Implement per input config overrides if provided

    # config_overrides: Dict[str, any] = get_config_overrides(
    #     input_file_id, per_input_config_overrides)
    config_overrides = {}
    config = setup_config(
        config_in,
        external_state_in,
        config_overrides,
    )
    dump_config_to_file_json(config, f'{output_directory}/processed_config.json')

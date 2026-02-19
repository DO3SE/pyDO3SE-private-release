from pathlib import Path
from pprint import pprint
import json
import os
import click
from typing import Optional

from do3se_phenology.plots import plot_phenology_from_config
from do3se_phenology.units import TimeTypes
from pyDO3SE.version import config_version
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Config.generate_config import generate_config
from pyDO3SE.Config.config_migration import Migrations
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.Output.process_outputs import (
    dump_config_to_file_json,
    dump_model_processes_info_to_string,
)
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
    """Generate a template config file based on the response to some simple questions."""
    click.echo("Welcome to pyDO3SE config generator")
    click.echo(gsto_method)
    # TODO: Fix prompts
    if gsto_method is None:
        gsto_method = click.prompt("Please enter photosynthesis method", type=str)
    click.echo(f"Config set to {gsto_method} defaults")
    new_config = generate_config(
        **{
            "Land_Cover.parameters.0.gsto.method": gsto_method,
        }
    )
    dump_config_to_file_json()(new_config, out_location)


def migrate_config(config_dir: Path, filename: Path, override):
    config_location = f"{config_dir}/{filename}"
    with open(config_location) as config_file_data:
        print(f'Running migrations on "{config_location}/{filename}"')
        config = json.load(config_file_data)
        input_version = config.get("VERSION", 0)
        try:
            migrated_config = Migrations.run_migrations(config, input_version)
        except Exception as e:
            raise Exception(f"Failed to migrate {filename}. {e}")
        out_dir = config_dir  # os.path.dirname(config_location)
        out_file_name = os.path.splitext(filename)[0]
        out_file_name = (
            f"{out_file_name}.json" if override else f"{out_file_name}_{config_version}.json"
        )
        with open(f"{out_dir}/{out_file_name}", "w") as out_location:
            json.dump(migrated_config, out_location, indent=4)


@click.command()
@click.option(
    "--override/--no-override",
    default=False,
    help="If true then will override the input config. Otherwise a copy is created.",
)
@click.argument("config_location")
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
@click.option("--detailed/--short", default=False)
def available_outputs(
    detailed: bool = False,
):
    """Print available outputs for output data and graphs."""
    if detailed:
        pprint(output_fields_map)
    else:
        print("ID\tShort\tType")
        print("------\t------")
        print("\n".join([f"{o.id}\t{o.short}\t{o.type}" for o in output_fields_map.values()]))


@click.command()
@click.option("--detailed/--short", default=False)
@click.option("--allow-errors/--strict", default=False)
@click.option("--out-location", default=None, type=click.Path(exists=False))
@click.option("--base-config-path", default=None, type=click.Path(exists=False))
@click.argument("config-location", type=click.Path(exists=True))
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
    config = config_loader(config_location, base_config_path, "json")

    hours = list(range(24))
    full_model_processes_out = full_model_processes(config, hours)
    flattened_process_comments = dump_model_processes_info_to_string(
        full_model_processes_out, detailed, allow_errors
    )

    if out_location:
        with open(out_location, "w") as f:
            f.write("\n".join(flattened_process_comments))
    else:
        print("\n".join(flattened_process_comments))


@click.option(
    "--input-data-file", type=click.Path(exists=True), default=None, help="Input data csv file"
)
@click.option(
    "--base-config-file",
    default=None,
    help="The base config file path",
    type=click.Path(exists=True),
)
@click.argument(
    "output-directory",
    # type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    "config-file",
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.command()
def plot_phenology(
    config_file: Path,
    output_directory: Path,
    input_data_file: Optional[Path] = None,
    base_config_file: Optional[Path] = None,
    day_count: int = 365,
):
    """Plot the phenology from config file."""
    config: Config_Shape = config_loader(config_file, base_config_file, "json")
    os.makedirs(output_directory, exist_ok=True)

    day_count = (
        int(config.Location.end_day or 0) - int(config.Location.start_day or 0)
        if config.Location.start_day is not None and config.Location.end_day is not None
        else day_count
    )

    plot_phenology_from_config(
        config.Land_Cover.parameters[0].phenology,
        config.Land_Cover.phenology_options,
        nP=config.Land_Cover.nP,
        output_location=Path(f"{output_directory}/phenology.png"),
        # TODO: Add external data input
        day_count=day_count,
        plot_dd=config.Land_Cover.phenology_options.time_type == TimeTypes.JULIAN_DAY,
        plot_td=config.Land_Cover.phenology_options.time_type == TimeTypes.THERMAL_TIME,
        plot_f_phen=True,
        plot_lengths=True,
        plot_carbon=config.Land_Cover.phenology_options.time_type == TimeTypes.THERMAL_TIME,
        plot_growing=config.Land_Cover.phenology_options.time_type == TimeTypes.THERMAL_TIME,
    )


@click.option(
    "--grid-coord", type=click.STRING, default="0,0", help="Grid coordinate to use"
)
@click.option(
    "--grid-overrides-file", type=click.Path(exists=True), default=None, help="Input data csv file"
)
@click.option(
    "--grid-overrides-field-map-path",
    type=click.Path(exists=True),
    default=None,
    help="Input data csv file",
)
@click.option(
    "--base-config-file",
    default=None,
    help="The base config file path",
    type=click.Path(exists=True),
)
@click.argument(
    "output-directory",
    # type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    "config-file",
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.command()
def plot_grid_phenology(
    config_file: Path,
    output_directory: Path,
    grid_overrides_file: Path,
    grid_overrides_field_map_path: Path,
    base_config_file: Optional[Path] = None,
    day_count: int = 365,
    grid_coord: str = "0,0",
):
    """Plot the phenology from config file."""
    from pyDO3SE.Grid_Model.setup_grid_model import process_grid_config as process_grid_config_fn

    processed_grid_config = process_grid_config_fn(
        config_file,
        grid_overrides_file=grid_overrides_file,
        grid_overrides_field_map_path=grid_overrides_field_map_path,
        base_config_file=base_config_file,
        grid_coord=tuple(map(int, grid_coord.split(","))),
    )
    os.makedirs(output_directory, exist_ok=True)

    day_count = (
        int(processed_grid_config.Location.end_day or 0)
        - int(processed_grid_config.Location.start_day or 0)
        if processed_grid_config.Location.start_day is not None
        and processed_grid_config.Location.end_day is not None
        else day_count
    )

    plot_phenology_from_config(
        processed_grid_config.Land_Cover.parameters[0].phenology,
        processed_grid_config.Land_Cover.phenology_options,
        nP=processed_grid_config.Land_Cover.nP,
        output_location=Path(f"{output_directory}/phenology.png"),
        day_count=day_count,
        plot_dd=processed_grid_config.Land_Cover.phenology_options.time_type
        == TimeTypes.JULIAN_DAY,
        plot_td=processed_grid_config.Land_Cover.phenology_options.time_type
        == TimeTypes.THERMAL_TIME,
        plot_f_phen=True,
        plot_lengths=True,
        plot_carbon=processed_grid_config.Land_Cover.phenology_options.time_type
        == TimeTypes.THERMAL_TIME,
        plot_growing=processed_grid_config.Land_Cover.phenology_options.time_type
        == TimeTypes.THERMAL_TIME,
    )



@click.option(
    "--grid-coord", type=click.STRING, default="0,0", help="Grid coordinate to use"
)
@click.option(
    "--grid-overrides-file", type=click.Path(exists=True), default=None, help="Input data csv file"
)
@click.option(
    "--grid-overrides-field-map-path",
    type=click.Path(exists=True),
    default=None,
    help="Input data csv file",
)
@click.option(
    "--base-config-file",
    default=None,
    help="The base config file path",
    type=click.Path(exists=True),
)
@click.argument(
    "output-directory",
)
@click.argument(
    "config-file",
    type=click.Path(exists=True),
)
@click.command()
def plot_grid_f_sw_curve(
    config_file: Path,
    output_directory: Path,
    grid_overrides_file: Path,
    grid_overrides_field_map_path: Path,
    base_config_file: Optional[Path] = None,
    grid_coord: str = "0,0",
):
    """Plot the f_sw curve from config file for single grid coord."""
    from pyDO3SE.Grid_Model.setup_grid_model import process_grid_config as process_grid_config_fn
    from do3se_met.plots import plot_f_SW_curve

    processed_grid_config = process_grid_config_fn(
        config_file,
        grid_overrides_file=grid_overrides_file,
        grid_overrides_field_map_path=grid_overrides_field_map_path,
        base_config_file=base_config_file,
        grid_coord=tuple(map(int, grid_coord.split(","))),
    )
    os.makedirs(output_directory, exist_ok=True)

    fig, ax = plot_f_SW_curve(
        f_SW_method=processed_grid_config.Land_Cover.parameters[0].gsto.f_SW_method,
        ASW_FC=processed_grid_config.soil_moisture.ASW_FC,
        fmin=processed_grid_config.Land_Cover.parameters[0].gsto.fmin,
        ASW_min=processed_grid_config.Land_Cover.parameters[0].gsto.ASW_min,
        ASW_max=processed_grid_config.Land_Cover.parameters[0].gsto.ASW_max,
        SWP_min=processed_grid_config.Land_Cover.parameters[0].gsto.SWP_min,
        SWP_max=processed_grid_config.Land_Cover.parameters[0].gsto.SWP_max,
    )
    fig.savefig(f"{output_directory}/f_SW_curve.png")


@click.option(
    "--input-data-file", type=click.Path(exists=True), default=None, help="Input data csv file"
)
@click.option(
    "--base-config-file",
    default=None,
    help="The base config file path",
    type=click.Path(exists=True),
)
@click.argument(
    "output-directory",
    # type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    "config-file",
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.command()
def process_config(
    config_file: Path,
    output_directory: Path,
    base_config_file: Optional[Path] = None,
):
    """Process the config by combining the config and base config.

    There are also some setup processes applied to the config such as setting up
    the phenology parameters.

    Note this does not include parameters defined by the external state.

    """

    from pyDO3SE.setup_model import setup_config, get_config_overrides
    from pyDO3SE.Output.process_outputs import dump_config_to_file_json

    config_in = config_loader(config_file, base_config_file, "json")

    external_state_in = None  # TODO: Implement external state if provided
    input_file_id = None  # TODO: Implement input file id if provided
    per_input_config_overrides = {}  # TODO: Implement per input config overrides if provided

    # config_overrides: Dict[str, any] = get_config_overrides(
    #     input_file_id, per_input_config_overrides)
    config_overrides = {}
    config = setup_config(
        config_in,
        external_state_in,
        config_overrides,
    )
    dump_config_to_file_json(config, f"{output_directory}/processed_config.json")


@click.option(
    "--grid-coord", type=click.STRING, default="0,0", help="Grid coordinate to use"
)
@click.option(
    "--grid-overrides-file", type=click.Path(exists=True), default=None, help="Input data csv file"
)
@click.option(
    "--grid-overrides-field-map-path",
    type=click.Path(exists=True),
    default=None,
    help="Input data csv file",
)
@click.option(
    "--base-config-file",
    default=None,
    help="The base config file path",
    type=click.Path(exists=True),
)
@click.argument(
    "output-directory",
    # type=click.Path(),
    # prompt='Enter output directory',
    # help='The location to save outputs',
)
@click.argument(
    "config-file",
    type=click.Path(exists=True),
    # prompt='Enter config location(.json)',
    # help='The location of the config file. Should be a json file',
)
@click.command()
def process_grid_config(
    config_file: Path,
    output_directory: Path,
    grid_overrides_file: Path,
    grid_overrides_field_map_path: Path,
    base_config_file: Optional[Path] = None,
    grid_coord: str = "0,0",
):
    """Process a config with grid overrides by combining the config and base config.

    There are also some setup processes applied to the config such as setting up
    the phenology parameters.

    This copies parts of initialize_grid_configs

    """
    from pyDO3SE.Grid_Model.setup_grid_model import process_grid_config as process_grid_config_fn

    processed_config = process_grid_config_fn(
        config_file,
        grid_overrides_file,
        grid_overrides_field_map_path,
        base_config_file,
        tuple(map(int, grid_coord.split(","))),
    )
    os.makedirs(output_directory, exist_ok=True)
    dump_config_to_file_json(processed_config, Path(f"{output_directory}/processed_config.json"))

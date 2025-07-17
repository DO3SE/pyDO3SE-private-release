from dataclasses import replace
from typing import Optional
import json
import click
from pathlib import Path
from collections import namedtuple

from data_helpers.encoders import AdvancedJsonEncoder
from data_helpers.data_loaders import load_json_to_cls
from thermal_time.calcs import calc_thermal_time_range

from do3se_phenology.switchboard import process_phenology_config
from do3se_phenology.phyllochron_dvi import calculate_t_l_from_leaf_f_phen
from do3se_phenology.data_access import get_var_from_file, get_td_from_file, get_vars_from_file
from do3se_phenology.config import PhenologyConfig


@click.group()
def cli():
    """Main cli entrypoint."""
    click.echo("Welcome to do3se_phenology")


ExternalData = namedtuple('ExternalData', ['td', 'dd', 'leaf_fphen'])


@cli.command()
@click.option('--t-b', default=0, help='Base temperature for thermal time')
@click.option('-o', '--output-file', default=None, help='Location to save processed config', type=click.Path())
@click.argument(
    'config-path',
    type=click.Path(exists=True),
)
@click.argument(
    'data-path',
    type=click.Path(exists=True),
)
def convert_phenology(
    data_path: Path,
    config_path: Path,
    t_b: float = 0,
    output_file: Optional[Path] = None
):
    click.echo(f"Generating phenology for {config_path} using data in {data_path}")
    # 1. load config into dataclass
    config = load_json_to_cls(config_path, PhenologyConfig)
    Ts_C, dd, leaf_fphen = get_vars_from_file(data_path, ['Ts_C', 'dd', 'leaf_fphen'])
    external_data = ExternalData(
        td=calc_thermal_time_range(Ts_C, t_b=t_b),
        dd=dd,
        leaf_fphen=leaf_fphen,
    )
    processed_configs = [process_phenology_config(
        config.model_config, species_config, external_data) for species_config in config.species_config]

    processed_config = replace(
        config,
        model_config=processed_configs[0][0],
        species_config=[p[1] for p in processed_configs]
    )

    json_data = json.dumps(
        processed_config,
        cls=AdvancedJsonEncoder, indent=4, sort_keys=True)

    if output_file:
        with open(output_file, 'w') as out_file_raw:
            out_file_raw.write(json_data)
    else:
        click.echo(json_data)


@cli.command()
def td_phenology():
    """Calculate the phenology variables from the thermal time data."""
    raise NotImplementedError()


@cli.command()
@click.option('--t-b', default=0, help='Base temperature for thermal time')
@click.option('--tsc-header', default="Ts_C", help='temperature header in file')
@click.option('--dd-header', default="dd", help='julian day header in file')
@click.option('--leaf-fphen-header', default="leaf_fphen", help='leaf fphen header in file')
@click.argument(
    'data-path',
    type=click.Path(exists=True),
)
def leaf_fphen_to_tl(
    data_path: Path,
    t_b: float = 0,
    tsc_header: str = "Ts_C",
    dd_header: str = "dd",
    leaf_fphen_header: str = "leaf_f_phen",
):
    """Convert fphen data to t_l values."""

    td_data = get_td_from_file(data_path, t_b, tsc_header)
    dd_data = get_var_from_file(data_path, dd_header)
    leaf_fphen_data = get_var_from_file(data_path, leaf_fphen_header)
    out = calculate_t_l_from_leaf_f_phen(td_data, leaf_fphen_data, dd_data)
    click.echo(out._asdict())


if __name__ == "__main__":
    cli()


def test_cli():
    cli()

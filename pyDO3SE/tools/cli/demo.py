import os
from shutil import copyfile
import click

from datetime import datetime
from pyDO3SE import version
from pyDO3SE import main


@click.command()
@click.option('--verbose', default=False, help='Run full logging output')
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--runnotes', default='', help='Add a note to the run')
def run(verbose=False, runid=0, runnotes='Demo run'):
    """Run pyDO3SE in demo mode with example inputs and config."""
    click.echo("Running pyDO3SE in demo mode")

    demo_config = "./examples/demo/configs/demo_config.json"
    demo_data = "./examples/demo/inputs/demo_data.csv"
    demo_output_dir = f"./examples/demo/outputs/{runid or version}"
    log_level = 1 if verbose else 0

    # Create output dir
    os.makedirs(demo_output_dir, exist_ok=True)

    # Copy config to this directory
    copyfile(demo_config, f'{demo_output_dir}/config.json')

    output_fields_to_graph = ['gsto_l', 'A_n']

    final_state, output_logs, config_processed, initial_state, external_state = main.main(
        demo_config,
        demo_data,
        demo_output_dir,
        output_fields=output_fields_to_graph,
        runid=runid,
        runnotes=runnotes,
        log_level=log_level,
        start_day=1,
        end_day=365,
        debug=True,
    )

    click.echo("Complete!")
    click.echo(f"Output file located in {demo_output_dir}")


@click.command()
@click.option('--verbose', default=False, help='Run full logging output')
def batch(verbose=False):
    """Run pyDO3SE in demo mode with example inputs and config."""
    # TODO: Make sure this works
    click.echo("Running pyDO3SE in multirun demo mode")
    demo_config_directory = "./examples/spanish_wheat/configs"
    demo_data_directory = "./examples/spanish_wheat/data/"
    demo_output_dir = "./examples/spanish_wheat/outputs"
    log_level = 1 if verbose else 0
    click.echo(f'Started at: {datetime.now()}')
    config_files = (f for f in os.listdir(demo_config_directory) if f[-5:] == '.json')
    data_files = (f for f in os.listdir(demo_data_directory) if f[-4:] == '.csv')
    for config in config_files:
        for data in data_files:
            click.echo(f'Running config: {config} on datafile: {data}')
            main.main(config, data, demo_output_dir, log_level=log_level)
    click.echo(f'Finished at: {datetime.now()}')

    click.echo("Complete!")
    click.echo(f"Output file located in {demo_output_dir}")

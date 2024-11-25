import os
from shutil import copyfile
import click

from datetime import datetime
from pyDO3SE import version
from pyDO3SE import main
from pyDO3SE.Output.OutputConfig import OutputOptions, output_results_only_options


@click.command()
@click.option('--verbose', default=False, help='Run full logging output')
@click.option('--runid', default='0', help='Set a run id for output naming')
@click.option('--runnotes', default='', help='Add a note to the run')
@click.option('--example-dir', default="./examples/demo", help='The path to the example directory', type=click.Path(exists=True))
def run(verbose=False, runid=0, runnotes='Demo run', example_dir=None):
    """Run pyDO3SE in demo mode with example inputs and config."""
    click.echo("Running pyDO3SE in demo mode")

    demo_config = f"{example_dir}/configs/demo_config.json"
    demo_data = f"{example_dir}/inputs/demo_data.csv"
    demo_output_dir = f"{example_dir}/outputs/{runid or version}"

    # Create output dir
    os.makedirs(demo_output_dir, exist_ok=True)

    # Copy config to this directory
    copyfile(demo_config, f'{demo_output_dir}/config.json')

    output_fields_to_graph = ",".join(['gsto_l', 'A_n'])

    main.single(
        config_file=demo_config,
        data_file=demo_data,
        output_directory=demo_output_dir,
        verbose=verbose,
        runid=runid,
        plot_fields=output_fields_to_graph,
        runnotes=runnotes,
        output_options=OutputOptions(
            save_model_processes=True,
            save_model_processes_detailed=True,
        ),
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

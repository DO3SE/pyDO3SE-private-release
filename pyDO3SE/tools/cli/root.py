
import click

from pyDO3SE.version import version
from . import run as run_commands
from . import config as config_commands
from . import analysis as analysis_commands
from . import demo as demo_commands
from . import setup as setup_commands

@click.group()
def cli():
    """Main cli entrypoint."""
    click.echo(f"Welcome to pyDO3SE version: {version}")
    click.echo(f"==== Update Info ====")
    click.echo(f"\tCLI major update. Commands split into run, config and analysis.")
    click.echo(f"\trun is now `run single`.")
    click.echo(f"\tmulti-run is now `run batch`.")
    click.echo(f"=====================\n")

@cli.group()
def run():
    pass

@cli.group()
def config():
    pass


@cli.group()
def analysis():
    pass

@cli.group()
def demo():
    pass


@cli.group()
def setup():
    pass

run.add_command(run_commands.single)
run.add_command(run_commands.batch)
run.add_command(run_commands.grid_run)
run.add_command(run_commands.grid_run_init)

config.add_command(config_commands.available_outputs)
config.add_command(config_commands.generate)
config.add_command(config_commands.migrate)
config.add_command(config_commands.output_process_list)
config.add_command(config_commands.plot_phenology)

analysis.add_command(analysis_commands.compare_outputs)
analysis.add_command(analysis_commands.multi_run_final_value)

demo.add_command(demo_commands.run)
demo.add_command(demo_commands.batch)

setup.add_command(setup_commands.generate_project_dir)

if __name__ == '__main__':
    cli()

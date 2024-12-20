
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
    from do3se_met.version import version as do3se_met_version
    from do3se_phenology.version import version as do3se_phenology_version
    from proflow.version import version as proflow_version
    click.echo(f"Welcome to pyDO3SE version: {version}")
    click.echo(f"==== Update Info ====")

    click.echo(f"=====================\n")
    click.echo(f"===== dependency versions ======\n")

    click.echo(f"do3se_met: {do3se_met_version}\n")
    click.echo(f"do3se_phenology: {do3se_phenology_version}\n")
    click.echo(f"proflow: {proflow_version}\n")

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
config.add_command(config_commands.process_config)

analysis.add_command(analysis_commands.compare_outputs)
analysis.add_command(analysis_commands.multi_run_final_value)

demo.add_command(demo_commands.run)
demo.add_command(demo_commands.batch)

setup.add_command(setup_commands.generate_project_dir)
setup.add_command(setup_commands.print_outputs)

if __name__ == '__main__':
    cli()

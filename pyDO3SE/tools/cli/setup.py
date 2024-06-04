from pathlib import Path
import os
import click

from pyDO3SE.Config.generate_config import generate_config
from pyDO3SE.Output.process_outputs import dump_config_to_file_json


@click.command()
@click.argument("path", type=click.Path())
def generate_project_dir(
    path: Path,
):
    """Generate the directory structure for DO3SE runs.

    Parameters
    ----------
    path : Path
        location to setup directory

    """
    click.echo(f"Generating project directory in {path}")
    os.makedirs(path, exist_ok=True)
    os.makedirs(f"{path}/data", exist_ok=True)
    os.makedirs(f"{path}/scripts", exist_ok=True)
    os.makedirs(f"{path}/runs", exist_ok=True)
    os.makedirs(f"{path}/runs/01", exist_ok=True)
    os.makedirs(f"{path}/runs/01/configs", exist_ok=True)
    os.makedirs(f"{path}/runs/01/inputs", exist_ok=True)
    base_config_path = f"{path}/runs/01/base_config.json"
    with open(f"{path}/runs/01/run.sh", 'w') as f:
        f.write("""python -m pyDO3SE run batch \
--project_directory={path}/runs/01 \
--base_config_file={base_config_path}
""")

    with open(f"{path}/README.md", 'w') as f:
        f.write("""# PROJECT TITLE

Project directory auto generated.

## Setup
...

## Run
...
""")
    new_config = generate_config()
    dump_config_to_file_json(new_config, base_config_path)

    click.echo(f"Generating project directory in {path} - COMPLETE")

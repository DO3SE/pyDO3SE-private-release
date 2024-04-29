from pathlib import Path
import os
import click
import pandas as pd

from pyDO3SE.Config.generate_config import generate_config
from pyDO3SE.Output.process_outputs import dump_config_to_file_json
from pyDO3SE.Output.Output_Shape import output_fields


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


@click.command()
def print_outputs():
    df = pd.DataFrame([field._asdict() for field in output_fields])
    df.to_csv('do3se_output_fields.csv')
    print("outputs saved to do3se_output_fields.csv")

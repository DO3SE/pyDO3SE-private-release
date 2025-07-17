from pathlib import Path
import click
from typing import List

from pyDO3SE.Output.run_comparisons import create_comparison_graphs, get_output_files
from pyDO3SE.Output.Output_Shape import output_fields as model_output_fields
from pyDO3SE.Output.process_outputs import extract_final_results


@click.option('--day_range', default='0,365', help='Set the start and end day')
@click.option('--filter', default='', help='Filter directory')
@click.option('--additional-input-dirs', default='', help='Additional directories to compare as comma seperated list')
@click.argument(
    'fields_to_graph',
    # prompt='Enter fields to compare',
    # help='The fields to run comparisons',
    nargs=-1,
)
@click.argument(
    'output_directory',
    # prompt='Enter output directory',
    # help='The location to save outputs',
    nargs=1,
)
@click.argument(
    'input_directory',
    # prompt='Enter model outputs directory',
    # help='The location of model outputs',
    nargs=1,
)
@click.command()
def compare_outputs(
    input_directory: Path,
    output_directory: Path,
    fields_to_graph: List[str] = [''],
    day_range: str = None,  # "0,365",
    filter: str = '',
    additional_input_dirs: str = '',
):
    """Create comparison graphs of all files in a directory.

    Each output to compare should be in a subdirectory. The name of the directory will
    be the name of the series.
    """
    start_day, end_day = [int(d) for d in day_range.split(
        ',')] if day_range is not None else [None, None]
    click.echo(f"Day range: {start_day} to {end_day}")
    _fields_to_graph = fields_to_graph if len(fields_to_graph) > 0 else [
        o.id for o in model_output_fields]
    click.echo(f"Comparing fields {_fields_to_graph}")

    all_input_dirs = (additional_input_dirs and additional_input_dirs.split(',') or []) + [input_directory]
    click.echo(f"Comparing outputs from {all_input_dirs}")
    data_files = get_output_files(all_input_dirs, filter)
    if len(data_files) == 0:
        raise ValueError("No output files found!")
    click.echo(f"Running Comparisons on outputs")
    create_comparison_graphs(
        data_files,
        _fields_to_graph,
        output_directory,
        use_versioned_outdir=False,
        start_day=start_day,
        end_day=end_day,
        use_versioned_outfile=False,
    )

    click.echo("Complete")


@click.argument(
    'field_name',
    # prompt='field name to extract',
    nargs=1,
)
@click.argument(
    'output_directory',
    # prompt='Enter output directory',
    # help='The location to save outputs',
    nargs=1,
)
@click.argument(
    'input_directory',
    # prompt='Enter model outputs directory',
    # help='The location of model outputs',
    nargs=1,
)
@click.command()
def multi_run_final_value(
    input_directory: Path,
    output_directory: Path,
    field_name: str,
):
    """Get the final value from each run (i.e. pody)"""
    extract_final_results(input_directory, output_directory, field_name)

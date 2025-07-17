"""Test running a few lines of data to make sure ozone deposition is correct."""

from pathlib import Path
from typing import NamedTuple
import warnings
import pandas as pd

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main


def run_with_config(runid: str, project_dir: Path, config_file: str, input_file: str, project_path_overrides: dict,) -> main.MainOutput:
    project_paths = main.get_project_paths(project_dir, **project_path_overrides)
    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)

    # Create output dir
    main.create_run_path_directories(run_paths)

    output_options = output_results_only_options()
    output_options.save_hourly_output_data = True
    output_options.save_processed_config = True
    output_options.save_final_state = True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = main.single(
            config_file=run_paths.config_path,
            data_file=run_paths.input_data_file_path,
            output_directory=run_paths.output_directory,
            base_config_file=project_paths.base_config_path,
            plot_fields=None,
            runid=runid,
            verbose=0,
            output_options=output_options,
        )
    return out


class Setup(NamedTuple):
    runid: str
    config_file: str
    input_file: str
    project_path_overrides: dict|None=None
    overrides: dict|None=None


class TestRunOutput(NamedTuple):
    out: main.MainOutput
    hourly_output: pd.DataFrame


def get_test_run(setups: list[Setup], run_outputs: dict[str, TestRunOutput], project_dir: str):

    def get_setup(runid: str):
        return next(s for s in setups if s[0] == runid)


    def test_run(runid: str) -> TestRunOutput:
        # global run_outputs
        setup = get_setup(runid)
        runid, config_file, input_file, project_path_overrides, overrides = setup

        if runid in run_outputs:
            return run_outputs[runid]
        try:
            out = run_with_config(
                runid=runid,
                project_dir=Path(project_dir),
                config_file=config_file,
                input_file=input_file,
                project_path_overrides=project_path_overrides or {},
                **(overrides or {}),
            )
        except Exception as e:
            print(f"Failed to run {runid}")
            raise e

        final_state, output_logs, final_config, initial_state, external_state = out
        run_outputs[runid] = TestRunOutput(
            out=out,
            hourly_output=pd.DataFrame(output_logs),
        )
        return run_outputs[runid]

    return test_run

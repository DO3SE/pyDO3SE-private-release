"""Contains common types and classes"""

from typing import NamedTuple
import pandas as pd
from pathlib import Path

from pyDO3SE.Config import Config_Shape
from pyDO3SE.Model_State import Model_State_Shape


class RunFiles(NamedTuple):
    config: Config_Shape
    state: Model_State_Shape
    per_input_config_overrides: pd.DataFrame


class ProjectPaths(NamedTuple):
    """Paths specific to the project run.

    This includes paths to multiple configs and inputs that will be iterated over.

    Parameters
    ----------

    project_dir: str
        _description_
    config_dir: str
        _description_
    input_data_dir: str
        _description_
    preprocess_map_path: str
        _description_
    base_config_path: str
        _description_
    base_state_path: str
        _description_
    observed_diurnal_path: str
        _description_
    runs_dir: str
        _description_
    per_input_config_overrides: str
        _description_

    """

    project_dir: Path | None = None
    config_dir: Path | None = None
    input_data_dir: Path | None = None
    preprocess_map_path: Path | None = None
    base_config_path: Path | None = None
    base_state_path: Path | None = None
    observed_diurnal_path: Path | None = None
    runs_dir: Path | None = None
    per_input_config_overrides: Path | None = None


class RunPaths(NamedTuple):
    """Paths specific to this config and input file.

    Parameters
    ----------
    run_id: str
        _description_
    run_name: str
        _description_
    run_dir: str
        ?
    log_path: str
        path to save log outputs from this run
    config_path: str
        path to the config file to run
    output_directory: str
        path to save outputs for this run
    input_data_file_path: str
        _description_
    config_run_dir: str
        _description_
    comparisons_dir: str
        _description_
    output_filename: str
        _description_
    input_file_id: str
        input file name without extension
    config_id: str
        config file name without extension

    """

    run_id: str
    run_name: str
    run_dir: str
    config_path: Path
    output_directory: Path
    input_data_file_path: Path
    output_filename: str
    input_file_id: str
    config_id: str
    comparisons_dir: Path | None = None
    log_path: Path | None = None
    config_run_dir: Path | None = None

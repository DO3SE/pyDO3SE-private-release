"""Contains common types and classes"""

from typing import NamedTuple
import pandas as pd

from pyDO3SE.Config import Config_Shape
from pyDO3SE.Model_State import Model_State_Shape


class RunFiles(NamedTuple):
    config: Config_Shape
    state: Model_State_Shape
    per_input_config_overrides: pd.DataFrame


class ProjectPaths(NamedTuple):
    project_dir: str = None
    config_dir: str = None
    input_data_dir: str = None
    preprocess_map_path: str = None
    base_config_path: str = None
    base_state_path: str = None
    observed_diurnal_path: str = None
    runs_dir: str = None
    per_input_config_overrides: str = None


class RunPaths(NamedTuple):
    run_id: str = None
    run_name: str = None
    run_dir: str = None
    log_path: str = None
    config_path: str = None
    output_directory: str = None
    input_data_file_path: str = None
    config_run_dir: str = None
    comparisons_dir: str = None
    output_filename: str = None
    input_file_id: str = None
    config_id: str = None

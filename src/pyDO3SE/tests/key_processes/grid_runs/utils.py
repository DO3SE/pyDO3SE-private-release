"""A set of tests that run the hourly model then compare the output against the previous version."""

import shutil
import os
from pathlib import Path
from typing import Callable
from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.Grid_Model import setup_grid_model
from pyDO3SE.util.logger import Logger
from pyDO3SE.Grid_Model.setup_grid_model import (
    GridProjectPaths,
    GridRunFiles,
    GridRunPaths,
)


class TestSetup:
    run_paths: GridRunPaths
    project_paths: GridProjectPaths
    multi_file_netcdf: bool
    runid: str
    project_dir: str | Path
    seperate_state_path: bool
    runnotes: list[str] | str
    loaded_run_files: GridRunFiles
    log_level: int
    outputs: list
    logs: list
    output_shape: tuple[int, int]
    output_fields: list[str]
    grid_coords: list[tuple[int, int]]
    grid_x_size: int
    grid_y_size: int
    logger_main: Callable
    regex_multi_file_filter: str
    netcdf_loader_kwargs: dict


def _assertTestSetup(self: TestSetup):
    assert self.multi_file_netcdf != "SETME", "Must set multi_file_netcdf in test class"
    assert self.runid != "SETME", "Must set runid in test class"
    assert self.project_dir != "SETME", "Must set project_dir in test class"


def _setup(self: TestSetup):
    print("Setting up")
    runid = self.runid
    project_dir = self.project_dir
    config_id = "bangor_wheat"
    # multi_file_netcdf = self.multi_file_netcdf
    self.seperate_state_path = True

    self.runnotes = ""
    self.log_level = 2

    project_paths = setup_grid_model.get_grid_project_paths(project_dir, runid)
    run_paths = setup_grid_model.get_grid_run_paths(project_paths, config_id)
    loaded_run_files = setup_grid_model.load_grid_run_files(project_paths, run_paths)

    # Clean up previous test output
    try:
        print(f"Removing {project_paths.run_dir}")
        shutil.rmtree(project_paths.run_dir)
    except FileNotFoundError:
        pass
    setup_grid_model.create_grid_run_path_directories(run_paths)

    self.loaded_run_files = loaded_run_files
    self.project_paths = project_paths
    self.run_paths = run_paths
    self.outputs = []
    self.logs = []

    grid_coords, grid_x_size, grid_y_size = setup_grid_model.get_grid_coords_from_file(run_paths.run_mask_path)

    self.output_shape = (grid_x_size, grid_y_size)
    self.grid_coords = grid_coords
    self.grid_x_size = grid_x_size
    self.grid_y_size = grid_y_size
    # self.logger_main = Logger(self.log_level, project_paths.log_path,
    #                           write_mode='w', set_as_default=True)
    self.logger_main = Logger(self.log_level, None, write_mode="w", set_as_default=True)


def _run_initialization(self):
    print("Running Init")
    try:
        e_state_overrides_dataset = xr.open_dataset(self.project_paths.e_state_overrides_file_path)
        initialized_config_gen, initialized_state_gen = setup_grid_model.init_grid_model(
            config=self.loaded_run_files.config,
            state=self.loaded_run_files.state,
            e_state_overrides_dataset=e_state_overrides_dataset,
            e_state_overrides_field_map=self.loaded_run_files.e_state_overrides_field_map,
            grid_coords=self.grid_coords,
            logger=self.logger_main,
            debug=True,
        )
        setup_grid_model.save_configs_from_generator(
            initialized_config_gen,
            self.grid_coords,
            self.run_paths.processed_configs_dir,
        )
        setup_grid_model.save_state_from_generator(
            initialized_state_gen,
            self.grid_coords,
            self.run_paths.live_state_dir,
        )
    except Exception as e:
        raise Exception("Failed to run init_grid_model")

    assert os.path.exists(f"{self.run_paths.live_state_dir}/0_0.state")
    assert os.path.exists(f"{self.run_paths.processed_configs_dir}/0_0.config")

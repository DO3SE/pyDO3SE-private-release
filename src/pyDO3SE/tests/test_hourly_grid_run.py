"""A set of tests that run the hourly model then compare the output against the previous version."""
import math
import numpy as np
from distutils.dir_util import copy_tree
import shutil
import warnings
import os
from datetime import datetime
import pytest
from pyDO3SE.Config.config_loader import config_loader_pickled
from pyDO3SE.Model_State.model_state_loader import dump_state_to_file, model_state_loader_quick
from pyDO3SE.setup_model import get_grid_coords_from_file

from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE import main
from pyDO3SE.util.logger import Logger, generate_run_notes

example_folders_netcdf = [
    'test_hourly_grid_single_file_hour+examples/net_cdf/single_file_hour',
    'test_hourly_grid_multi_file_range+examples/net_cdf/multi_file_range',
    # 'wrfchem+examples/net_cdf/wrfchem_post',
]
output_fields = ['dd', 'hr', 'gsto_canopy', 'td_dd',
                 'lai', 'pody', 'fst', 'canopy_height', 'micro_u']


# =============== Full grid sequential ======== #
# TODO: Only currently set up for WRFChem
@pytest.fixture(scope="class", params=example_folders_netcdf)
def setup_test_run_from_netcdf_grid_seq(request):
    [runid, project_dir] = request.param.split('+')
    # Inputs
    config_id = "bangor_wheat"
    multi_file_netcdf = "multi" in runid
    seperate_state_path = True

    runnotes = []
    log_level = 0

    project_paths = main.get_grid_project_paths(project_dir)
    run_paths = main.get_grid_run_paths(project_paths, config_id, runid)
    loaded_run_files = main.load_grid_run_files(project_paths, run_paths)

    # Clean up previous test output
    try:
        shutil.rmtree(run_paths.run_dir)
    except FileNotFoundError:
        pass
    main.create_grid_run_path_directories(run_paths)

    request.cls.project_paths = project_paths
    request.cls.run_paths = run_paths
    request.cls.outputs = []
    request.cls.logs = []

    grid_coords, grid_x_size, grid_y_size = get_grid_coords_from_file(
        run_paths.run_mask_path)

    output_shape = (grid_x_size, grid_y_size)
    request.cls.grid_coords = grid_coords
    request.cls.grid_x_size = grid_x_size
    request.cls.grid_y_size = grid_y_size

    logger_main = Logger(log_level)

    # Initialize the state and configs first.
    # In live runs this should only be ran for the first hour.
    # TODO: Get zero year from here
    try:
        main.main_partial_initialize(
            config=loaded_run_files.config,
            state=loaded_run_files.state,
            processed_config_dir=run_paths.processed_configs_dir,
            state_out_path=run_paths.live_state_dir,
            e_state_overrides_file_path=project_paths.e_state_overrides_file_path,
            e_state_overrides_field_map=loaded_run_files.e_state_overrides_field_map,
            grid_coords=grid_coords,
            logger=logger_main,
        )
    except Exception as e:
        raise Exception("Failed to run main_partial_initialize")
#
    assert os.path.exists(f"{run_paths.live_state_dir}/0_0.state")
    assert os.path.exists(f"{run_paths.processed_configs_dir}/0_0.config")
    copy_tree(run_paths.live_state_dir, run_paths.initial_state_dir)

    errors = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            start_time = datetime.now()
            main.main_grid_seq_per_config(
                project_paths=project_paths,
                run_paths=run_paths,
                loaded_run_files=loaded_run_files,
                grid_coords=grid_coords,
                output_shape=output_shape,
                runnotes=runnotes,
                output_fields=output_fields,
                log_level=log_level,
                seperate_live_state=seperate_state_path,
                multi_file_netcdf=multi_file_netcdf,
                logger=logger_main,
            )
            end_time = datetime.now()
            duration = end_time - start_time

            with open(f'{run_paths.output_data_dir}/notes.log', 'w') as f:
                log_notes = generate_run_notes(
                    runnotes='',
                    time_taken=duration,
                    # time_taken_setup=setup_duration,
                    # config_version=config_version,
                    # model_version=model_version,
                    errors=errors,
                )
                f.write("\n".join(log_notes))

        # Save a human readable copy of the final state
        for x, y in grid_coords:
            file_name = f"{x}_{y}"
            dump_state_to_file(model_state_loader_quick(
                f"{run_paths.live_state_dir}/{file_name}.state"), f"{run_paths.final_state_dir}/{file_name}.json")

    except Exception as e:
        errors.append((f"Project dir: {project_dir} failed", e))

    if len(errors) > 0:
        for m, e in errors:
            print(m)
        print(errors)
        raise errors[0][1]
    # run_type = model.processed_config.Land_Cover.parameters[0].gsto.method
    # request.cls.run_type = run_type


@pytest.mark.usefixtures('setup_test_run_from_netcdf_grid_seq')
class TestHourlyRunNetCDFGridSeq:

    def test_should_run_without_errors(self):
        pass

    def test_should_have_combined_all_grid_outputs_into_a_single_file(self):
        output_files = sorted(os.listdir(self.run_paths.output_data_dir))
        output_data_file = f"{self.run_paths.output_data_dir}/{output_files[-1]}"
        try:
            with open(output_data_file) as f:
                assert f is not None
        except FileNotFoundError:
            raise Exception("final output not saved")

    def test_should_contain_correct_coordinates(self):
        output_files = sorted(os.listdir(self.run_paths.output_data_dir))
        output_data_file = f"{self.run_paths.output_data_dir}/{output_files[-1]}"
        ds = xr.open_dataset(output_data_file)
        assert set(ds.coords) == {'time', 'lat', 'lon'}

    def test_should_have_only_outputed_required_fields(self):
        output_files = sorted(os.listdir(self.run_paths.output_data_dir))
        output_data_file = f"{self.run_paths.output_data_dir}/{output_files[-1]}"
        ds = xr.open_dataset(output_data_file)
        assert list(ds.keys()) == output_fields

        for f in output_fields:
            assert ds[f] is not None

    def test_output_should_be_correct_shape(self):
        output_files = sorted(os.listdir(self.run_paths.output_data_dir))
        output_data_file = f"{self.run_paths.output_data_dir}/{output_files[-1]}"
        ds = xr.open_dataset(output_data_file)
        [x, y, t] = ds[output_fields[0]].shape
        assert x == self.grid_x_size
        assert y == self.grid_y_size
        assert t > 0  # This should match row count

    def test_combined_outputs(self):
        output_data_dir = self.run_paths.output_data_dir
        output_files = [f for f in sorted(os.listdir(output_data_dir)) if ".log" not in f]
        ds_files = [xr.open_dataset(f"{output_data_dir}/{o}") for o in output_files]
        ds = xr.merge(ds_files)
        [x, y, t] = ds.dd.shape
        assert x == self.grid_x_size
        assert y == self.grid_y_size
        assert t == 24 * 4 * 2

    def test_should_output_state_for_each_grid_tile(self):
        output_state_files_dir = f"{self.run_paths.live_state_dir}"
        state_files = os.listdir(output_state_files_dir)
        assert len(state_files) == len(self.grid_coords)

    def test_should_store_prev_state_for_each_grid_tile(self):
        output_state_files_dir = f"{self.run_paths.prev_state_dir}"
        state_files = os.listdir(output_state_files_dir)
        print(self.run_paths.prev_state_dir)
        assert len(state_files) == len(self.grid_coords)

    def test_should_create_processed_config_for_each_grid_tile(self):
        output_config_files_dir = f"{self.run_paths.processed_configs_dir}"
        config_files = os.listdir(output_config_files_dir)
        assert len(config_files) == len(self.grid_coords)

    def test_should_have_set_external_config_override(self):
        output_config_files_dir = f"{self.run_paths.processed_configs_dir}"
        config_files = os.listdir(output_config_files_dir)
        config_file = config_loader_pickled(f'{output_config_files_dir}/{config_files[0]}')
        assert config_file.Land_Cover.parameters[0].phenology.key_dates.sowing is not None
        assert not math.isnan(config_file.Land_Cover.parameters[0].phenology.key_dates.sowing)

    def test_should_have_advanced_hour_in_state(self):
        # output_state_files_hr_1_dir = f"{self.run_paths.initial_state_dir}"
        output_state_files_hr_1_dir = f"{self.run_paths.prev_state_dir}"
        output_state_files_hr_2_dir = f"{self.run_paths.live_state_dir}"
        state_files = os.listdir(output_state_files_hr_1_dir)
        # input_state = None

        # input_state = model_state_loader_quick(f"{input_state_files_dir}/{state_files[0]}")

        state_data_hr_01 = model_state_loader_quick(
            f"{output_state_files_hr_1_dir}/{state_files[0]}")

        state_data_hr_02 = model_state_loader_quick(
            f"{output_state_files_hr_2_dir}/{state_files[0]}")

        # TODO: Get better way of testing this
        # assert state_data_hr_01.temporal.dd == 274
        assert state_data_hr_02.temporal.dd == 368
        # assert state_data_hr_01.temporal.hr == 22
        assert state_data_hr_02.temporal.hr == 23

    def test_should_handle_multiple_years(self):
        """When using multiple years the day of year should increase beyond 365."""
        # TODO: Implement this test
        pass

    # def test_compare_initial_to_final(self, snapshot):
    #     final_state, output_logs, processed_config, initial_state = self.output
    #     initial_state_copy = {k: v for k, v in asdict(initial_state).items() if k != "prev_hour"}
    #     final_state_copy = {k: v for k, v in asdict(final_state).items() if k != "prev_hour"}
    #     compared = get_dict_differences("model_state", initial_state_copy, final_state_copy)
    #     snapshot.assert_match({
    #         "total_changes": len(compared),
    #         "changes": sorted(compared)
    #     }, "Changes")

# ====== FUNCTIONAL TESTS ========== #

project_directories=[
    'examples/net_cdf/multi_file_range',
]

class TestIntegrated:
    """Full integrated test"""

    # @pytest.mark.skip(reason="Out of date")
    @pytest.mark.parametrize('project_directory', project_directories)
    def test_grid_run(self, project_directory):
        # cli args
        runid = 'TestIntegrated'
        log_level = 2
        run_notes=[]
        multi_file_netcdf=True
        config_name = 'bangor_wheat'

        project_paths = main.get_grid_project_paths(project_directory)
        output_fields_to_graph = ['gsto_l', 'pody', 'lai']
        main.init_all_grid_model_configs(
            project_paths,
            runid,
            log_level=log_level,
        )
        main.main_grid_run(
            project_paths,
            output_fields_to_graph,
            multi_file_netcdf,
            runid,
            run_notes,
            log_level=log_level,
            save_final_state=True,
        )

        # Check outputs
        ds = xr.load_dataset(
            f'{project_directory}/runs/{runid}/bangor_wheat/outputs_grid/output_data_2017-10-01_01.nc')
        assert ds['pody'] is not None

        run_paths = main.get_grid_run_paths(project_paths, config_name, runid)
        grid_coords, grid_x_size, grid_y_size = get_grid_coords_from_file(run_paths.run_mask_path)
        T = 1
        X = grid_x_size
        Y = grid_y_size
        assert ds['pody'].shape == (X, Y, T)
        assert ds['lai'].shape == (X, Y, T)

        assert ds['time'][0].values.astype(str) == '2017-10-01T01:00:00.000000000'

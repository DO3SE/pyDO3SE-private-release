"""A set of tests that run the hourly model then compare the output against the previous version."""
import math
import pandas as pd
import numpy as np
import warnings
import os
from pathlib import Path
from datetime import datetime
import pytest
from data_helpers.list_helpers import flatten_list
from pyDO3SE.Output.utils import get_multi_dimension_output_fields
from pyDO3SE.Config.config_loader import config_loader_pickled
from pyDO3SE.Model_State.model_state_loader import dump_state_to_file, model_state_loader_quick
from pyDO3SE.Grid_Model.setup_grid_model import get_grid_coords_from_file

from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.Grid_Model import setup_grid_model
from pyDO3SE.Grid_Model import run_grid_model
from pyDO3SE.util.logger import Logger, generate_run_notes

from .utils import _assertTestSetup, _setup, _run_initialization, TestSetup


def _default_run(self: TestSetup, **kwargs):
    self.logger_main("Running Model")
    errors = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            start_time = datetime.now()
            run_grid_model.main_grid_seq_per_config(
                project_paths=self.project_paths,
                run_paths=self.run_paths,
                loaded_run_files=self.loaded_run_files,
                grid_coords=self.grid_coords,
                output_shape=self.output_shape,
                runnotes=self.runnotes,
                output_fields=self.output_fields,
                seperate_live_state=self.seperate_state_path,
                multi_file_netcdf=self.multi_file_netcdf,
                regex_multi_file_filter=self.regex_multi_file_filter,
                logger=self.logger_main,
                netcdf_loader_kwargs=self.netcdf_loader_kwargs,
                debug=True
            )
            end_time = datetime.now()
            duration = end_time - start_time

            with open(f'{self.run_paths.config_run_dir}/notes.log', 'w') as f:
                log_notes = generate_run_notes(
                    runnotes='',
                    time_taken=str(duration),
                    # time_taken_setup=setup_duration,
                    # config_version=config_version,
                    # model_version=model_version,
                    errors=errors,
                )
                f.write("\n".join(log_notes))

        # Save a human readable copy of the final state
        for x, y in self.grid_coords:
            file_name = f"{x}_{y}"
            dump_state_to_file(model_state_loader_quick(
                Path(f"{self.run_paths.live_state_dir}/{file_name}.state")), Path(f"{self.run_paths.final_state_dir}/{file_name}.json"))

    except Exception as e:
        errors.append((f"Project dir: {self.project_paths.project_dir} failed", e))

    if len(errors) > 0:
        for m, e in errors:
            print(m)
        print(errors)
        raise errors[0][1]
    print("Complete")





@pytest.fixture(scope="class")
def before_all(request):
    _self = request.cls
    assert not _self.complete
    print("Before all")
    _assertTestSetup(_self)
    _setup(_self)
    _run_initialization(_self)
    _self.logger_main("Test logger")
    _default_run(_self)
    _self.complete = True


class TestHourlyGridRun(TestSetup):
    class Base(TestSetup):
        config_path = None
        config_id = "bangor_wheat"
        output_fields = ['dd', 'hr', 'gsto_canopy', 'td_dd',
                         'canopy_lai', 'pody', 'fst', 'canopy_height', 'micro_u', 'par', 'ts_c']
        input_fields = ['SWDOWN', 'HFX_FORCE', 'td_2m',
                        'rh', 'o3', 'wspeed', 'pres', 'RAINNC', 'SNOWH']
        multi_file_netcdf = "SETME" # type: ignore
        runid = "SETME"
        project_dir = "SETME"
        regex_multi_file_filter = None # type: ignore
        netcdf_loader_kwargs = {}
        complete = False
        expected_total_days = 8

        def test_should_run_without_errors(self):
            assert self.complete is True

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

        # No longer valid if we include multi layer fields
        # def test_should_have_only_outputed_required_fields(self):
        #     output_files = sorted(os.listdir(self.run_paths.output_data_dir))
        #     output_data_file = f"{self.run_paths.output_data_dir}/{output_files[-1]}"
        #     ds = xr.open_dataset(output_data_file)
        #     assert list(ds.keys()) == self.output_fields

        #     for f in self.output_fields:
        #         assert ds[f] is not None

        def test_should_have_only_outputed_required_multidimensional_fields(self):
            output_files = sorted(os.listdir(self.run_paths.output_data_dir))
            output_data_file = f"{self.run_paths.output_data_dir}/{output_files[-1]}"
            ds = xr.open_dataset(output_data_file)

            multi_dimensional_fields =flatten_list([[f, *get_multi_dimension_output_fields(f)] for f in self.output_fields])
            assert len(multi_dimensional_fields) > len(self.output_fields)
            for f in multi_dimensional_fields:
                assert ds[f] is not None
                for x, y in self.grid_coords:
                    assert all((not np.isnan(v) and v is not None) for v in ds[f].sel(y=y, x=x, drop=True).values)


        def test_output_should_be_correct_shape(self):
            output_files = sorted(os.listdir(self.run_paths.output_data_dir))
            output_data_file = f"{self.run_paths.output_data_dir}/{output_files[-1]}"
            ds = xr.open_dataset(output_data_file)
            [x, y, t] = ds[self.output_fields[0]].shape
            assert x == self.grid_x_size
            assert y == self.grid_y_size
            assert t > 0  # This should match row count

        def test_output_should_match_input_shape(self):
            output_files = sorted(os.listdir(self.run_paths.output_data_dir))
            input_files = sorted(os.listdir(f"{self.project_dir}/inputs"))

            ds_out_files = [xr.open_dataset(
                f"{self.run_paths.output_data_dir}/{f}") for f in output_files]
            ds_in_files = [xr.open_dataset(f"{self.project_dir}/inputs/{f}") for f in input_files]
            ds_in = xr.concat(ds_in_files, dim="Time")
            ds_out = xr.concat(ds_out_files, dim="time")  # Note: Output is time not Time
            output_field = self.output_fields[0]
            # NOTE: time location in output is different
            [x_out, y_out, t_out] = ds_out[output_field].shape
            [t_in, x_in, y_in] = ds_in[self.input_fields[0]].shape
            assert x_out == x_in, f"Output x {x_out} != input x {x_in}"
            assert y_out == y_in, f"Output y {y_out} != input y {y_in}"
            assert t_out == t_in, f"Output t {t_out} != input t {t_in}"

        def test_combined_outputs(self):
            output_data_dir = self.run_paths.output_data_dir
            output_files = [f for f in sorted(os.listdir(output_data_dir)) if ".log" not in f]
            ds_files = [xr.open_dataset(f"{output_data_dir}/{o}") for o in output_files]
            ds = xr.merge(ds_files)
            [x, y, t] = ds.dd.shape
            assert x == self.grid_x_size
            assert y == self.grid_y_size
            assert t == (24 * self.expected_total_days)

        def test_should_output_state_for_each_grid_tile(self):
            output_state_files_dir = f"{self.run_paths.live_state_dir}"
            state_files = os.listdir(output_state_files_dir)
            assert len(state_files) == len(self.grid_coords)

        def test_should_store_prev_state_for_each_grid_tile(self):
            output_state_files_dir = f"{self.run_paths.prev_state_dir}"
            state_files = os.listdir(output_state_files_dir)
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

        def test_time_in_output_should_match_time_in_input(self):
            output_files = sorted(os.listdir(self.run_paths.output_data_dir))
            input_files = sorted(os.listdir(f"{self.project_dir}/inputs"))

            ds_out_files = [xr.open_dataset(
                f"{self.run_paths.output_data_dir}/{f}") for f in output_files]
            ds_in_files = [xr.open_dataset(f"{self.project_dir}/inputs/{f}") for f in input_files]
            ds_in = xr.concat(ds_in_files, dim="Time")
            ds_out = xr.concat(ds_out_files, dim="time")
            assert len(ds_out.time.values) == len(
                ds_in.XTIME.values), f"Output time {len(ds_out.time.values)} != input time {len(ds_in.XTIME.values)}"
            assert ds_out.time.values.min() == ds_in.XTIME.values.min(
            ), f"Output time min {ds_out.time.values.min()} != input time min {ds_in.XTIME.values.min()}"
            assert ds_out.time.values.max() == ds_in.XTIME.values.max(
            ), f"Output time max {ds_out.time.values.max()} != input time max {ds_in.XTIME.values.max()}"
            assert (ds_out.time.values == ds_in.XTIME.values).all()

        def test_should_have_advanced_hour_in_state(self):
            """After running a model input netcdf file the hour should be advanced by 1.

            Assuming the input data is hourly and the model is run hourly.
            """

            output_files = sorted(os.listdir(self.run_paths.output_data_dir))

            ds_out_file_last = xr.open_dataset(
                f"{self.run_paths.output_data_dir}/{output_files[-1]}")
            assert ds_out_file_last.hr.sel(y=0, x=0, drop=True).values[-1] == 2.0
            assert ds_out_file_last.dd.sel(y=0, x=0, drop=True).values[-1] == 369.0
            output_state_files_dir = f"{self.run_paths.live_state_dir}"

            state_files = os.listdir(output_state_files_dir)
            final_state = model_state_loader_quick(
                Path(f"{output_state_files_dir}/{state_files[0]}"))
            assert final_state.temporal.hr == 2.0
            assert final_state.temporal.dd == 369.0

        def test_should_handle_multiple_years(self):
            """When using multiple years the day of year should increase beyond 365."""
            # TODO: Implement this test
            pass

        def test_should_handle_skipping_rows_in_data(self):
            """"""
            pass

        def test_should_handle_offsetting_time_in_data(self):
            """Sometimes the input data time is in GMT rather than local timezone"""
            pass

        def test_input_temperature_data_should_match_output_temperature_data(self):
            """The input temperature data should not be modified by the model before output.

            This is a good way of ensuring that all data was merged and loaded correctly.

            """
            output_files = sorted(os.listdir(self.run_paths.output_data_dir))
            input_files = sorted(os.listdir(f"{self.project_dir}/inputs"))

            ds_out_files = [xr.open_dataset(
                f"{self.run_paths.output_data_dir}/{f}") for f in output_files]
            ds_in_files = [xr.open_dataset(f"{self.project_dir}/inputs/{f}") for f in input_files]
            ds_in = xr.concat(ds_in_files, dim="Time")
            ds_out = xr.concat(ds_out_files, dim="time")

            for x, y in self.grid_coords:
                assert (ds_out.ts_c.sel(y=y, x=x, drop=True).values ==
                        ds_in.td_2m.sel(y=y, x=x, drop=True).values).all()

    # @pytest.mark.skip(reason="Takes too long")
    # @pytest.mark.usefixtures("before_all")
    # class TestSingleFileHour(Base):
    #     runid = "test_hourly_grid_single_file_hour"
    #     project_dir = "examples/net_cdf/single_file_hour"
    #     multi_file_netcdf = False

    @pytest.mark.usefixtures("before_all")
    class TestSingleFileRangeCombined(Base):
        runid = "test_hourly_grid_single_file_range_combined"
        project_dir = "examples/net_cdf/single_file_range"
        multi_file_netcdf = True
        # regex_multi_file_filter = '[0-9]{4}'
        netcdf_loader_kwargs = dict(
            parallel=False,
            concat_dim="Time",
            combine="nested",
        )

        # days cropped to 7
        expected_total_days = 8

        def test_single_file_range_worked(self):
            pass

        def test_should_have_merged_all_files(self):
            """Should have merged all inputs into 1 before running the model."""
            output_files = os.listdir(
                f"{self.project_dir}/runs/{self.runid}/{self.config_id}/outputs_grid")
            assert len(
                output_files) == 1, f"Model should have only produced 1 nc file but {len(output_files)} files"

    @pytest.mark.usefixtures("before_all")
    class TestSingleFileRange(Base):
        """Test grid run when the input data is split across multiple files per day.

        In this case we do not merge the files before running the model.
        """

        runid = "test_hourly_grid_single_file_range"
        project_dir = "examples/net_cdf/single_file_range"
        multi_file_netcdf = False
        # regex_multi_file_filter = '[0-9]{4}'
        netcdf_loader_kwargs = dict(
            # parallel=False,
            # concat_dim="Time",
            # combine="nested",
        )

        # days cropped to 7
        expected_total_days = 8

        def test_single_file_range_worked(self):
            pass

    @pytest.mark.usefixtures("before_all")
    class TestSingleFileRangeGroupByYear(Base):
        """Test grid run when the input data is split across multiple files per day.

        In this we merge the files by year and then run the model.
        """

        runid = "test_hourly_grid_single_file_range_group_year"
        project_dir = "examples/net_cdf/single_file_range"
        multi_file_netcdf = True
        regex_multi_file_filter = '[0-9]{4}'
        expected_total_days = 8
        netcdf_loader_kwargs = dict(
            parallel=False,
            concat_dim="Time",
            combine="nested",
        )

        def test_single_file_range_worked(self):
            pass

        def test_should_have_merged_all_files(self):
            """Should have merged all inputs into 1 before running the model."""
            output_files = os.listdir(
                f"{self.project_dir}/runs/{self.runid}/{self.config_id}/outputs_grid")
            assert len(
                output_files) == 2, f"Model should have only produced 2 nc files but {len(output_files)} files"

    # TODO: Fix this test
    # @pytest.mark.usefixtures("before_all")
    # class TestMultiFileRange(Base):
    #     """Test grid run when the input data is split across multiple files per
    #     variable and each files represents a range of data."""
    #     runid = "test_hourly_grid_multi_file_range"
    #     project_dir = "examples/net_cdf/multi_file_range"
    #     regex_multi_file_filter = '[0-9]{4}-[0-9]{2}'
    #     multi_file_netcdf = True


project_directories = [
    'examples/net_cdf/multi_file_range',
]


class TestIntegrated:
    """Full integrated test"""

    def test_grid_run(self):
        import warnings
        warnings.warn("TEST")
        project_directory = 'examples/net_cdf/multi_file_range'
        # cli args
        runid = 'TestIntegrated'
        log_level = 2
        run_notes = ""
        multi_file_netcdf = True
        config_name = 'bangor_wheat'
        regex_multi_file_filter = "[0-9]{4}-[0-9]{2}"
        netcdf_loader_kwargs = {

        }

        project_paths = setup_grid_model.get_grid_project_paths(project_directory, runid)
        output_fields = ['gsto_l', 'pody', 'canopy_lai',
                         'td', 'ts_c', 'o3_ppb_zr', 'vpd', 'uh_zr', 'par']
        logger_main = Logger(log_level, log_to_file=None,
                             set_as_default=True, write_mode='w', flush_per_log=True)
        logger_main(f"=== Running Integration test on {project_directory}")
        print(f"=== Running Integration test on {project_directory}")
        input_files = os.listdir(f'{project_directory}/inputs')
        ds_in = xr.open_dataset(f'{project_directory}/inputs/{input_files[0]}')
        t = ds_in.XTIME.values.astype(str)[0:19]
        time_data = pd.date_range(str(t[0]), periods=1, freq="1H")
        time_string = f'{time_data[0].year}-{str(time_data[0].month).zfill(2)}-{str(time_data[0].day).zfill(2)}_{str(time_data[0].hour).zfill(2)}'
        target_file_name = f'output_data_{time_string}.nc'
        input_time_length = ds_in.XTIME.size

        setup_grid_model.init_all_grid_model_configs(
            project_paths,
            logger_main,
        )
        run_grid_model.main_grid_run(
            project_paths,
            output_fields,
            multi_file_netcdf,
            run_notes,
            logger=logger_main,
            save_final_state=True,
            parallel=False,
            parallel_args=setup_grid_model.ParallelArgs(8, 0, .01),
            regex_multi_file_filter=regex_multi_file_filter,
            netcdf_loader_kwargs=netcdf_loader_kwargs,
            debug=True,
        )

        # Check outputs
        ds_out = xr.load_dataset(
            f'{project_directory}/runs/{runid}/bangor_wheat/outputs_grid/{target_file_name}')
        assert ds_out['pody'] is not None

        run_paths = setup_grid_model.get_grid_run_paths(project_paths, config_name)
        grid_coords, grid_x_size, grid_y_size = get_grid_coords_from_file(run_paths.run_mask_path)

        K = list(ds_in.keys())[0]
        T, X, Y = ds_in[K].shape
        assert grid_x_size == X
        assert grid_y_size == Y
        assert input_time_length == T
        assert ds_out['pody'].shape == (X, Y, T)
        assert ds_out['canopy_lai'].shape == (X, Y, T)

        # Check no missing vars in outputs
        for k in ds_out.keys():
            assert ds_out[k][0, 0, 0] != run_grid_model.MISSING_OUTPUT_VALUE, f"{k} is Missing"
            assert ds_out[k][0, 0, 0] is not None, f"{k} is Missing"
            assert not math.isnan(ds_out[k][0, 0, 0]), f"{k} is Missing"

        assert ds_out['time'][0].values.astype(str) == ds_in['XTIME'][0].values.astype(str)


class TestGridRunTimezoneShift:

    def test_can_offset_hours(self):
        pass

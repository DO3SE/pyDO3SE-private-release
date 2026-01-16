"""A set of tests that run the grid model on global data.

This tests that we can run with gridded hourly offset.
"""
import math
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

from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.Grid_Model import run_grid_model
from pyDO3SE.util.logger import generate_run_notes

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


class TestGlobalGridRun(TestSetup):
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
    class TestGlobalGridRun(Base):
        runid = "TestGlobalGridRun"
        project_dir = "examples/net_cdf/single_file_global"
        multi_file_netcdf = False
        # regex_multi_file_filter = '[0-9]{4}'
        netcdf_loader_kwargs = dict(
            parallel=False,
            # concat_dim="Time",
            # combine="nested",
        )

        # days cropped to 7
        expected_total_days = 1

        def test_single_file_global_worked(self):
            pass

        def test_should_have_merged_all_files(self):
            """Should have merged all inputs into 1 before running the model."""
            output_files = os.listdir(
                f"{self.project_dir}/runs/{self.runid}/{self.config_id}/outputs_grid")
            assert len(
                output_files) == 1, f"Model should have only produced 1 nc file but got {len(output_files)} files"


        def test_should_have_offset_hr_in_some_grid_cells(self):
            output_files = sorted(os.listdir(self.run_paths.output_data_dir))
            ds = xr.open_dataset(self.run_paths.output_data_dir + "/" + output_files[0])
            assert ds.hr.shape == (3, 2, 24)
            # Check at time = 0
            np.testing.assert_array_equal(
                ds.hr[:,:,0].values,
                [[0,0],[np.nan,np.nan],[12.0,np.nan]]
            )
            # Check at time = 23
            np.testing.assert_array_equal(
                ds.hr[:,:,23].values,
                [[23,23],[np.nan,np.nan],[(12.0 + 23)%24,np.nan]]
            )







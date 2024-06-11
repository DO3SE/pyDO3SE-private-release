# """A set of tests that run the hourly model then compare the output against the previous version."""
# import json
# import numpy as np
# import shutil
# import importlib
# import itertools
# import warnings
# import pandas as pd
# import os
# from datetime import datetime
# import pytest
# from dataclasses import asdict
# from data_helpers.diff import diff
# from pyDO3SE.Model_State.model_state_loader import model_state_loader_quick

# from pyDO3SE.util.loader import csv_loader, json_loader
# from pyDO3SE.optional_dependencies import xarray as xr
# from pyDO3SE.External_State.external_state_loader import FileTypes
# from pyDO3SE import main
# from pyDO3SE.setup_model import LocationMethod
# from pyDO3SE.main import generate_run_notes

# example_folders = [f.path for f in os.scandir('examples/hourly_runs') if f.is_dir()]
# output_fields = ['dd', 'hr', 'gsto_canopy', 'td_dd',
#                  'lai', 'pody', 'fst', 'canopy_height', 'micro_u']


# @pytest.fixture(scope="class", params=example_folders)
# def setup_test_run(request):
#     project_dir = request.param
#     config_file = f"{project_dir}/processed_config.json"
#     initial_state_file = f"{project_dir}/initial_state.json"
#     input_data_file = f"{project_dir}/hourly_data.csv"
#     output_data_file = f"{project_dir}/hourly_output.csv"
#     output_state_file = f"{project_dir}/final_state.json"
#     request.cls.project_dir = project_dir
#     request.cls.outputs = []
#     request.cls.logs = []
#     errors = []
#     try:
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             start_time = datetime.now()
#             for out in main.main_hour(
#                 processed_config_path=config_file,
#                 external_data_row_path=input_data_file,
#                 previous_hour_state_path=initial_state_file,
#                 output_fields=output_fields,
#                 runid=f"{project_dir}",
#                 runnotes="Running from test",
#                 log_level=0,
#                 state_out_path=output_state_file,
#             ):
#                 final_state, output_logs, processed_config, initial_state = out
#                 request.cls.outputs.append(out)
#                 request.cls.logs.append(pd.DataFrame(output_logs))

#     except Exception as e:
#         errors.append((f"Project dir: {project_dir} failed", e))

#     if len(errors) > 0:
#         for m, e in errors:
#             print(m)
#         print(errors)
#         raise errors[0][1]
#     run_type = processed_config.Land_Cover.parameters[0].gsto.method
#     request.cls.run_type = run_type


# @pytest.mark.skip(reason="Broken!")
# @pytest.mark.usefixtures('setup_test_run')
# class TestHourlyRun:

#     def test_should_run_without_errors(self):
#         pass

#     def test_only_single_output(self):
#         """Test that when we only pass single x and y values that we get a single output"""
#         assert len(self.outputs) == 1

#     def test_compare_initial_to_final(self, snapshot):
#         final_state, output_logs, processed_config, initial_state = self.outputs[0]
#         initial_state_copy = {k: v for k, v in asdict(initial_state).items() if k != "prev_hour"}
#         final_state_copy = {k: v for k, v in asdict(final_state).items() if k != "prev_hour"}
#         compared = diff("model_state", initial_state_copy, final_state_copy)
#         snapshot.assert_match({
#             "total_changes": len(compared),
#             "changes": sorted(compared)
#         }, "Changes")

#         with open(f"{self.project_dir}/final_state.json") as outfile:
#             outfile_json = json.load(outfile)
#             with open(f"{self.project_dir}/expected_state.json") as expected_outfile:
#                 expected_outfile_json = json.load(expected_outfile)
#                 compared = diff("model_state", expected_outfile_json, outfile_json)
#                 if len(compared) > 0:
#                     compared_message = "\n".join(compared)
#                     raise AssertionError(
#                         f"Output and expected output do not match: \n expected -> actual \n\n {compared_message}")

#     def test_should_be_able_to_run_following_hour(self):
#         project_dir = self.project_dir
#         config_file = f"{project_dir}/processed_config.json"
#         initial_state_file = f"{project_dir}/final_state.json"
#         input_data_file = f"{project_dir}/hourly_data.csv"  # TODO: Replace with next hour of data
#         output_data_file = f"{project_dir}/hourly_output_next.csv"
#         output_state_file = f"{project_dir}/final_state_next.json"

#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             out = main.main_hour(
#                 processed_config_path=config_file,
#                 external_data_row_path=input_data_file,
#                 previous_hour_state_path=initial_state_file,
#                 output_data_file=output_data_file,  # Fix outputs
#                 output_state_file=output_state_file,
#                 output_fields=output_fields,
#                 runid=f"{project_dir}",
#                 runnotes="Running from test",
#                 log_level=0,
#             )


# example_folders_netcdf = ['wrfchem+examples/net_cdf/wrfchem']
# # example_folders_netcdf = [f.path for f in os.scandir('examples/net_cdf') if f.is_dir()]


# @pytest.fixture(scope="class", params=example_folders_netcdf)
# def setup_test_run_from_netcdf(request):
#     [runid, project_dir] = request.param.split('+')
#     config_file = f"{project_dir}/processed_config.json"
#     initial_state_file = f"{project_dir}/initial_state.json"
#     input_data_folder = f"{project_dir}/inputs"
#     multi_file_netcdf = len(os.listdir(input_data_folder)) > 1
#     input_data_file = input_data_folder if multi_file_netcdf else f"{input_data_folder}/{os.listdir(input_data_folder)[0]}"
#     output_data_dir = f"{project_dir}/outputs"
#     request.cls.project_dir = project_dir
#     request.cls.outputs = []
#     request.cls.logs = []

#     # Clean up previous test output
#     try:
#         os.remove(output_data_dir)
#     except FileNotFoundError:
#         pass

#     maps = importlib.import_module(f"{project_dir}/maps".replace('/', '.'))
#     netcdf_variable_map = maps.variable_map
#     met_preprocess_map = maps.preprocess_map

#     errors = []
#     try:
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             for out in main.main_hour(
#                 processed_config_path=config_file,
#                 external_data_row_path=input_data_file,
#                 previous_hour_state_path=initial_state_file,
#                 output_data_dir=output_data_dir,  # Fix outputs
#                 output_fields=output_fields,
#                 runid=runid,
#                 runnotes="Running from test",
#                 log_level=0,

#                 external_file_type=FileTypes.NETCDF,
#                 grid_coords=[(0,0)],
#                 netcdf_variable_map=netcdf_variable_map,
#                 met_preprocess_map=met_preprocess_map,
#                 multi_file_netcdf=multi_file_netcdf,
#                 output_to_netcdf=True,
#             ):
#                 final_state, output_logs, processed_config, initial_state = out
#                 request.cls.outputs.append(out)
#                 request.cls.logs.append(pd.DataFrame(output_logs))
#         # Make a copy for the version
#     except Exception as e:
#         errors.append((f"Project dir: {project_dir} failed", e))

#     if len(errors) > 0:
#         for m, e in errors:
#             print(m)
#         print(errors)
#         raise errors[0][1]
#     run_type = processed_config.Land_Cover.parameters[0].gsto.method
#     request.cls.run_type = run_type


# @pytest.mark.skip(reason="Broken!")
# @pytest.mark.usefixtures('setup_test_run_from_netcdf')
# class TestHourlyRunNetCDF:

#     def test_should_run_without_errors(self):
#         pass

#     def test_only_single_output(self):
#         """Test that when we only pass single x and y values that we get a single output"""
#         assert len(self.outputs) == 1

#     def test_compare_initial_to_final(self, snapshot):
#         final_state, output_logs, processed_config, initial_state = self.outputs[0]
#         initial_state_copy = {k: v for k, v in asdict(initial_state).items() if k != "prev_hour"}
#         final_state_copy = {k: v for k, v in asdict(final_state).items() if k != "prev_hour"}
#         compared = diff("model_state", initial_state_copy, final_state_copy)
#         snapshot.assert_match({
#             "total_changes": len(compared),
#             "changes": sorted(compared)
#         }, "Changes")

#     def test_should_output_netcdf_file(self):
#         output_data_file = f"{self.project_dir}/outputs/hourly_output.nc"
#         try:
#             with open(output_data_file) as f:
#                 assert f is not None
#         except FileNotFoundError:
#             raise Exception("final output not saved")
#         ds = xr.open_dataset(output_data_file)
#         assert ds.f_phen is not None
#         assert ds.pody is not None

#     def test_should_output_final_state(self):
#         output_state_file = f"{self.project_dir}/outputs/final_state.json"
#         try:
#             with open(output_state_file) as f:
#                 assert f is not None
#         except FileNotFoundError:
#             raise Exception("final state not saved")

# # =============== Full grid runs =====================

# @pytest.mark.skip(reason="Replaced by sequential run below")
# @pytest.fixture(scope="class", params=example_folders_netcdf)
# def setup_test_run_from_netcdf_grid(request):
#     [runid, project_dir] = request.param.split('+')
#     config_file = f"{project_dir}/processed_config.json"
#     initial_state_path = f"{project_dir}/state"
#     input_data_folder = f"{project_dir}/inputs"
#     multi_file_netcdf = len(os.listdir(input_data_folder)) > 1
#     input_data_file = input_data_folder if multi_file_netcdf else f"{input_data_folder}/{os.listdir(input_data_folder)[0]}"
#     output_data_dir = f"{project_dir}/outputs_grid"
#     state_out_path = f"{project_dir}/state_out"
#     processed_configs_path = f"{project_dir}/processed_configs"
#     request.cls.project_dir = project_dir
#     request.cls.outputs = []
#     request.cls.logs = []

#     # Clean up previous test output
#     # TODO: Clean up out dir
#     # try:
#     #     os.remove(output_data_dir)
#     # except FileNotFoundError:
#     #     pass

#     # == initialize directory == #

#     maps = importlib.import_module(f"{project_dir}/maps".replace('/', '.'))
#     netcdf_variable_map = maps.variable_map
#     met_preprocess_map = maps.preprocess_map

#     grid_x = list(range(49))
#     grid_y = list(range(49))
#     request.cls.grid_x_size = len(grid_x)
#     request.cls.grid_y_size = len(grid_y)

#     grid_coords = list(itertools.product(grid_x, grid_y))

#     # Initialize the state and configs first.
#     # In live runs this should only be ran for the first hour.
#     main_hour_initialize(
#         f"{project_dir}/configs",
#         processed_configs_path,
#         initial_state_path,
#         base_config_path=f"{project_dir}/base_config.json",
#         base_state_path=f"{project_dir}/base_state.json",
#         runid=runid,
#         grid_coords=grid_coords,
#     )

#     config_file = os.listdir(processed_configs_path)[0]

#     errors = []
#     try:
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             start_time = datetime.now()
#             for out in main.main_hour(
#                 processed_config_path=f"{processed_configs_path}/{config_file}",
#                 external_data_row_path=input_data_file,
#                 previous_hour_state_path=initial_state_path,
#                 output_data_dir=output_data_dir,  # Fix outputs
#                 output_fields=output_fields,
#                 runid=runid,
#                 runnotes="Running from test",
#                 log_level=0,

#                 external_file_type=FileTypes.NETCDF,
#                 grid_coords=grid_coords,
#                 netcdf_variable_map=netcdf_variable_map,
#                 met_preprocess_map=met_preprocess_map,
#                 multi_file_netcdf=multi_file_netcdf,
#                 output_to_netcdf=True,
#                 state_out_path=state_out_path,
#                 location_method=LocationMethod.EXTERNAL_STATE_INPUT,
#             ):

#                 final_state, output_logs, model = out
#                 # run_type = model.processed_config.Land_Cover.parameters[0].gsto.method
#                 # request.cls.run_type = run_type

#                 request.cls.outputs.append(out)
#                 request.cls.logs.append(pd.DataFrame(output_logs))
#             end_time = datetime.now()
#             duration = end_time - start_time
#             with open('notes.log', 'w') as f:
#                 f.write(f'\nRun took: {duration} for project_dir: {project_dir}')
#         # Make a copy for the version
#     except Exception as e:
#         errors.append((f"Project dir: {project_dir} failed", e))

#     if len(errors) > 0:
#         for m, e in errors:
#             print(m)
#         print(errors)
#         raise errors[0][1]
#     # run_type = model.processed_config.Land_Cover.parameters[0].gsto.method
#     # request.cls.run_type = run_type


# @pytest.mark.skip(reason="Broken!")
# @pytest.mark.usefixtures('setup_test_run_from_netcdf_grid')
# class TestHourlyRunNetCDFGrid:

#     def test_should_run_without_errors(self):
#         pass

#     def test_should_have_combined_all_grid_outputs_into_a_single_file(self):
#         output_data_file = f"{self.project_dir}/outputs_grid/output_data.nc"
#         try:
#             with open(output_data_file) as f:
#                 assert f is not None
#         except FileNotFoundError:
#             raise Exception("final output not saved")

#     def test_should_contain_correct_coordinates(self):
#         output_data_file = f"{self.project_dir}/outputs_grid/output_data.nc"
#         ds = xr.open_dataset(output_data_file)
#         assert set(ds.coords) == {'time', 'lat', 'lon'}

#     def test_should_have_only_outputed_required_fields(self):
#         output_data_file = f"{self.project_dir}/outputs_grid/output_data.nc"
#         ds = xr.open_dataset(output_data_file)
#         assert list(ds.keys()) == output_fields

#         for f in output_fields:
#             assert ds[f] is not None

#     def test_output_should_be_correct_shape(self):
#         output_data_file = f"{self.project_dir}/outputs_grid/output_data.nc"
#         ds = xr.open_dataset(output_data_file)
#         assert ds[output_fields[0]].shape == (self.grid_x_size, self.grid_y_size, 1)

#     def test_should_output_state_for_each_grid_tile(self):
#         output_state_files_dir = f"{self.project_dir}/state_out"
#         state_files = os.listdir(output_state_files_dir)
#         assert len(state_files) == self.grid_x_size * self.grid_y_size

#     def test_should_have_advanced_hour_in_state(self):
#         output_state_files_dir = f"{self.project_dir}/state_out"
#         state_files = os.listdir(output_state_files_dir)
#         input_state = None

#         with open(f"{self.project_dir}/initial_state.json") as f:
#             input_state = json.load(f)

#         with open(f"{self.project_dir}/state_out/{state_files[0]}") as f:
#             state_data = json.load(f)
#             assert state_data.get('temporal').get('hr') == input_state.get('temporal').get('hr') + 1

#     # def test_compare_initial_to_final(self, snapshot):
#     #     final_state, output_logs, processed_config, initial_state = self.output
#     #     initial_state_copy = {k: v for k, v in asdict(initial_state).items() if k != "prev_hour"}
#     #     final_state_copy = {k: v for k, v in asdict(final_state).items() if k != "prev_hour"}
#     #     compared = diff("model_state", initial_state_copy, final_state_copy)
#     #     snapshot.assert_match({
#     #         "total_changes": len(compared),
#     #         "changes": sorted(compared)
#     #     }, "Changes")


# # =============== Full grid sequential ======== #
# # TODO: Only currently set up for WRFChem
# @pytest.fixture(scope="class", params=example_folders_netcdf)
# def setup_test_run_from_netcdf_grid_seq(request):
#     [runid, project_dir] = request.param.split('+')
#     # Inputs
#     config_names = ['bangor_wheat']
#     input_data_folder = f"{project_dir}/inputs"
#     e_state_overrides_file_path = f"{project_dir}/e_state_overrides.nc"
#     multi_file_netcdf = False
#     variable_map_path = f"{project_dir}/variable_map.json"
#     preprocess_map_path = f"{project_dir}/preprocess_map.json"
#     e_state_overrides_field_map_path = f"{project_dir}/e_state_overrides_field_map.json"

#     base_config_path=f"{project_dir}/base_config.json"
#     base_state_path=f"{project_dir}/base_state.json"

#     variable_map = json_loader(variable_map_path)
#     preprocess_map = json_loader(preprocess_map_path)
#     e_state_overrides_field_map = json_loader(e_state_overrides_field_map_path)

#     run_dir = f"{project_dir}/runs/{runid}_test/{config_names[0]}"

#     # Clean up previous test output
#     try:
#         shutil.rmtree(run_dir)
#     except FileNotFoundError:
#         pass
#     finally:
#         os.makedirs(run_dir, exist_ok=True)
#     # == initialize directory == #

#     # Processed file locations
#     initial_state_path = f"{run_dir}/current_state/0"
#     output_data_dir = f"{run_dir}/outputs_grid"
#     processed_configs_dir = f"{run_dir}/processed_configs"
#     prev_state_dir = f"{run_dir}/current_state"
#     run_mask_path = f"{project_dir}/coords/{config_names[0]}.csv"

#     request.cls.project_dir = project_dir
#     request.cls.output_data_dir = output_data_dir
#     request.cls.state_path = prev_state_dir
#     request.cls.processed_configs_dir = processed_configs_dir
#     request.cls.outputs = []
#     request.cls.logs = []

#     os.makedirs(processed_configs_dir, exist_ok=True)
#     os.makedirs(prev_state_dir, exist_ok=True)
#     os.makedirs(output_data_dir, exist_ok=True)


#     grid_coords = np.array([[int(i['x']), int(i['y'])] for i in csv_loader(run_mask_path)])
#     grid_x_size, grid_y_size = np.ptp(grid_coords, axis=0) + [1,1]
#     request.cls.grid_coords = grid_coords
#     request.cls.grid_x_size = grid_x_size
#     request.cls.grid_y_size = grid_y_size

#     os.makedirs(initial_state_path, exist_ok=True)

#     # Initialize the state and configs first.
#     # In live runs this should only be ran for the first hour.
#     main_hour_initialize(
#         f"{project_dir}/configs/{config_names[0]}.json",
#         processed_configs_dir,
#         initial_state_path,
#         base_config_path=base_config_path,
#         base_state_path=base_state_path,
#         e_state_overrides_file_path=e_state_overrides_file_path,
#         e_state_overrides_field_map=e_state_overrides_field_map,
#         grid_coords=grid_coords,
#     )

#     errors = []
#     try:
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             start_time = datetime.now()
#             for i, f in enumerate(sorted(os.listdir(input_data_folder))):
#                 print(f"==== Running file: {f} =======")
#                 input_data_file = f"{input_data_folder}/{f}"

#                 # for tests we save state in seperate folders.
#                 # In practice we would save out override state
#                 previous_hour_state_path = f"{prev_state_dir}/{i}"
#                 state_out_path = f"{prev_state_dir}/{i+1}"
#                 os.makedirs(previous_hour_state_path, exist_ok=True)
#                 os.makedirs(state_out_path, exist_ok=True)

#                 for out in main.main_hour(
#                     processed_config_dir=f"{processed_configs_dir}",
#                     external_data_row_path=input_data_file,
#                     previous_hour_state_path=previous_hour_state_path,
#                     output_data_dir=output_data_dir,  # Fix outputs
#                     output_fields=output_fields,
#                     external_file_type=FileTypes.NETCDF,
#                     grid_coords=grid_coords,
#                     netcdf_variable_map=variable_map,
#                     met_preprocess_map=preprocess_map,
#                     multi_file_netcdf=multi_file_netcdf,
#                     output_to_netcdf=True,
#                     state_out_path=state_out_path,
#                 ):

#                     final_state, output_logs, model = out
#                     # run_type = model.processed_config.Land_Cover.parameters[0].gsto.method
#                     # request.cls.run_type = run_type

#                     request.cls.outputs.append(out)
#                     request.cls.logs.append(pd.DataFrame(output_logs))
#             end_time = datetime.now()
#             duration = end_time - start_time

#             with open(f'{output_data_dir}/notes.log', 'w') as f:
#                 log_notes = generate_run_notes(
#                     runnotes='',
#                     time_taken=duration,
#                     # time_taken_setup=setup_duration,
#                     # config_version=config_version,
#                     # model_version=model_version,
#                     errors=errors,
#                 )
#                 f.write("\n".join(log_notes))

#         # Make a copy for the version
#     except Exception as e:
#         errors.append((f"Project dir: {project_dir} failed", e))

#     if len(errors) > 0:
#         for m, e in errors:
#             print(m)
#         print(errors)
#         raise errors[0][1]
#     # run_type = model.processed_config.Land_Cover.parameters[0].gsto.method
#     # request.cls.run_type = run_type


# @pytest.mark.usefixtures('setup_test_run_from_netcdf_grid_seq')
# class TestHourlyRunNetCDFGridSeq:

#     def test_should_run_without_errors(self):
#         pass

#     def test_should_have_combined_all_grid_outputs_into_a_single_file(self):
#         output_data_file = f"{self.output_data_dir}/output_data.nc"
#         try:
#             with open(output_data_file) as f:
#                 assert f is not None
#         except FileNotFoundError:
#             raise Exception("final output not saved")

#     def test_should_contain_correct_coordinates(self):
#         output_data_file = f"{self.output_data_dir}/output_data.nc"
#         ds = xr.open_dataset(output_data_file)
#         assert set(ds.coords) == {'time', 'lat', 'lon'}

#     def test_should_have_only_outputed_required_fields(self):
#         output_data_file = f"{self.output_data_dir}/output_data.nc"
#         ds = xr.open_dataset(output_data_file)
#         assert list(ds.keys()) == output_fields

#         for f in output_fields:
#             assert ds[f] is not None

#     def test_output_should_be_correct_shape(self):
#         output_data_file = f"{self.output_data_dir}/output_data.nc"
#         ds = xr.open_dataset(output_data_file)
#         assert ds[output_fields[0]].shape == (self.grid_x_size, self.grid_y_size, 1)

#     def test_should_output_state_for_each_grid_tile(self):
#         output_state_files_dir = f"{self.state_path}/1"
#         state_files = os.listdir(output_state_files_dir)
#         assert len(state_files) == len(self.grid_coords)


#     def test_should_create_processed_config_for_each_grid_tile(self):
#         output_config_files_dir = f"{self.processed_configs_dir}"
#         config_files = os.listdir(output_config_files_dir)
#         assert len(config_files) == len(self.grid_coords)

#     def test_should_have_advanced_hour_in_state(self):
#         input_state_files_dir = f"{self.state_path}/0"
#         output_state_files_hr_1_dir = f"{self.state_path}/1"
#         output_state_files_hr_2_dir = f"{self.state_path}/2"
#         state_files = os.listdir(output_state_files_hr_1_dir)
#         input_state = None

#         input_state = model_state_loader_quick(f"{input_state_files_dir}/{state_files[0]}")

#         state_data_hr_01 = model_state_loader_quick(f"{output_state_files_hr_1_dir}/{state_files[0]}")

#         state_data_hr_02 = model_state_loader_quick(f"{output_state_files_hr_2_dir}/{state_files[0]}")

#         assert state_data_hr_01.temporal.hr == 0
#         assert state_data_hr_02.temporal.hr == 1

#     def test_should_handle_multiple_years(self):
#         """When using multiple years the day of year should increase beyond 365."""


#     # def test_compare_initial_to_final(self, snapshot):
#     #     final_state, output_logs, processed_config, initial_state = self.output
#     #     initial_state_copy = {k: v for k, v in asdict(initial_state).items() if k != "prev_hour"}
#     #     final_state_copy = {k: v for k, v in asdict(final_state).items() if k != "prev_hour"}
#     #     compared = diff("model_state", initial_state_copy, final_state_copy)
#     #     snapshot.assert_match({
#     #         "total_changes": len(compared),
#     #         "changes": sorted(compared)
#     #     }, "Changes")

# # ====== FUNCTIONAL TESTS ========== #
# class TestCli:


#     def test_grid_run(self):

#         main.main_grid_seq(
#             project_dir='examples/net_cdf/wrfchem',
#             runid='0',
#             init_model=True,
#             output_fields=['pody'],
#             runnotes='',
#             log_level=2,
#         )

#         ds = xr.load_dataset('examples/net_cdf/wrfchem/runs/0/bangor_wheat/outputs_grid/output_data.nc')
#         assert ds['pody'] is not None

#         run_mask_path = 'examples/net_cdf/wrfchem/coords/bangor_wheat.csv'
#         grid_coords = np.array([[int(i['x']), int(i['y'])] for i in csv_loader(run_mask_path)])
#         grid_x_size, grid_y_size = np.ptp(grid_coords, axis=0) + [1,1]
#         T = 1
#         X = grid_x_size
#         Y = grid_y_size
#         assert ds['pody'].shape == (X,Y,T)

#         assert ds['time'][0].values.astype(str) == '2017-10-01T01:00:00.000000000'



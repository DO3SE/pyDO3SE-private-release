"""Unit tests for main.py ."""
from copy import deepcopy
import numpy as np
import pandas as pd
from typing import List, Tuple
from unittest.mock import MagicMock, Mock
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.External_State.external_state_loader import (
    get_date_bounds_from_ext_data,
)

from pyDO3SE.Model_State.demo_data import example_state
from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes
from pyDO3SE.setup_model import (
    setup_config,
)
from .run_grid_model import (
    main_partial,
    ExternalStateIterable,
    setup_external_state_simple
)


__module_loc__ = __package__ + ".run_grid_model"


# def Obj(name: str, *args, **kwargs) -> object:
#     argsAsKwargs = {k._extract_mock_name(): k for k in args}
#     return type(name, (object,), {**kwargs, **argsAsKwargs})()


# == Mocked return values == #
# mock_config = MagicMock()
# mock_initial_state = MagicMock()
# mock_model_processes = MagicMock()
# mock_external_state = MagicMock()
# mock_parameters = MagicMock()
# mock_state = MagicMock()
# mock_state_out = MagicMock()
# mock_run_processes = MagicMock(return_value=mock_state_out)
# mock_loaded_run_files = MagicMock()
# mock_state_logs = [
#     {"a": 1, "b": 2},
#     {"a": 3, "b": 4},
#     {"a": 5, "b": 6},
# ]
# mock_run_notes = MagicMock()
# mock_processes = [1, 2, 3]
# mock_overrides = MagicMock()
# mock_overrides.debug = False
# mock_process_runner = Mock(
#     run_processes=mock_run_processes,
#     state_logs=mock_state_logs,
# )
# config_location = "demo_config_location"
# mock_data_location = "demo_data_location"
# mock_base_config_file = "demo_base_config.json"
# mock_initial_state_file = None


# class TestMainPartial:

#     class TestMainPartialPointLocation:
#         """Test that we can run the DO3SE model as a partial run on single location csv data.

#         Should handle the following scenarios:
#         - csv file provided as a single file covering all hours
#         - multiple csv files that can be concatenated

#         """
#         # TODO: Implement tests
#         pass

#     class TestMainPartialGrid:
#         """Test that we can run the DO3SE model as a partial run on netcdf gridded data.

#         Should handle the following scenarios:
#         - netcdf data provided as a file per hour where a file contains all required vars
#         - netcdf data provided as a file covering a time range which contains all required vars
#         - netcdf data provided as multiple files each containing a single variable over a time range

#         """

#         def set_defaults(self):
#             """Set the default parameters for a fully mocked run. """
#             # Default Input Params
#             defaults.config_name = "demo_config"
#             defaults.output_fields = ['a']
#             defaults.config_names = [defaults.config_name]
#             defaults.run_id = 0
#             defaults.grid_coords = [(0, 0), (1, 0), (2, 1)]
#             defaults.grid_x_size = 5  # this is range of x coords
#             defaults.grid_y_size = 6  # this is range of y coords
#             defaults.output_shape = (defaults.grid_x_size, defaults.grid_y_size)
#             defaults.cell_configs = [MagicMock() for _ in defaults.grid_coords]
#             defaults.cell_initial_states = [MagicMock() for _ in defaults.grid_coords]
#             defaults.cell_output_states = [MagicMock() for _ in defaults.grid_coords]
#             defaults.cell_external_states = [MagicMock() for _ in defaults.grid_coords]
#             # defaults.external_state_options = MagicMock()
#             # defaults.external_state_init_processes = (MagicMock())
#             defaults.model_processes = [MagicMock]
#             defaults.mock_overrides = mock_overrides
#             defaults.mock_overrides.state_out_path = None
#             defaults.mock_overrides.grid_coords = defaults.grid_coords
#             defaults.start_day = 274
#             defaults.end_day = 274
#             defaults.row_count = 1
#             defaults.e_state_dates = dict(
#                 # TODO: replace date strings with correct values below
#                 start_day=defaults.start_day,
#                 end_day=defaults.end_day,
#                 start_date="01-02-2017",
#                 end_date="28-02-2017",
#                 row_count=defaults.row_count,
#                 time_string="01-02-2017",
#                 hours=[hr for _ in range(defaults.start_day - defaults.end_day) for hr in range(24)]
#             )
#             # Expected Results
#             defaults.output_data = {
#                 "a": [row['a'] for row in mock_state_logs],
#                 "b": [row['b'] for row in mock_state_logs],
#             }

#             defaults.output_logs = [
#                 {"a": i + 1} for i in range(defaults.row_count)
#             ]

#             # Iterators
#             defaults.cell_output_states_i = cycle(defaults.cell_output_states)
#             defaults.cell_initial_states_i = cycle(defaults.cell_initial_states)
#             defaults.e_states_i = cycle(defaults.cell_external_states)

#         def set_integration_test_defaults(self):
#             """Setup the defaults when only mocking model run and outputs. """
#             # Full run without any mocks
#             defaults.variable_map = {
#                 "time": "XTIME",
#                 "_SHAPE": "td_2m",
#                 "Ts_C": "td_2m",
#                 "P": "pres",
#                 "PAR": "SWDOWN",
#                 "precip": "RAINNC",
#                 "RH": "rh",
#                 "u": "wspeed",
#                 "O3": "o3",
#                 "Hd": "HFX_FORCE",
#                 "snow_depth": "SNOWH"
#             }
#             defaults.data_location = 'examples/net_cdf/single_file_hour/inputs/demo_wrf_2017-12-27-00-00-00.nc'
#             defaults.multi_file_data = False
#             defaults.preprocess_map = {}
#             defaults.zero_year = 2017
#             defaults.netcdf_chunks = None
#             defaults.prev_hour_state_location = "demo_state_location"
#             defaults.output_fields = ['pody']

#             demo_config = config_loader(
#                 'examples/net_cdf/single_file_hour/configs/bangor_wheat.json',
#                 'examples/net_cdf/single_file_hour/base_config.json',
#             )

#             defaults.cell_configs = [demo_config for _ in defaults.grid_coords]

#             defaults.cell_states = [deepcopy(example_state) for _ in defaults.grid_coords]
#             defaults.cell_states_i = cycle(defaults.cell_states)

#             defaults.mock_state_out = MagicMock()
#             defaults.output_logs = [
#                 {"pody": 1},
#             ]

#         def mock_for_unit_test(defaults, mocker):
#             """Setup mocks for a unit test."""
#             # setup_funcs_to_mock = [
#             #     ('get_date_bounds_from_ext_data', defaults.e_state_dates.values()),
#             #     ('run_model', None, lambda *_, **__: (next(defaults.cell_output_states_i), defaults.output_logs)),
#             #     ('dump_state_to_file_quick', None),
#             #     ('Main_Overrides', mock_overrides),
#             #     ('dump_output_to_file_netcdf_grid', None),
#             #     ('external_state_init_processes', defaults.external_state_init_processes),
#             #     ('load_current_cell_state', None, lambda *_, **__: next(defaults.cell_initial_states_i)),
#             #     ('get_row_processes_hourly', defaults.model_processes),
#             #     ('load_external_state', defaults.e_states_i),
#             # ]

#             # mock_funcs(defaults, mocker, __module_loc__, setup_funcs_to_mock)
#             run_funcs_to_mock = [
#                 ('get_date_bounds_from_ext_data', defaults.e_state_dates.values()),
#                 ('run_model', None, lambda *_, **__: (next(defaults.cell_output_states_i), defaults.output_logs)),
#                 ('dump_state_to_file_quick', None),
#                 # ('Main_Overrides', mock_overrides),
#                 ('dump_output_to_file_netcdf_grid', None),
#                 ('external_state_init_processes', defaults.external_state_init_processes),
#                 ('load_current_cell_state', None, lambda *_, **__: next(defaults.cell_initial_states_i)),
#                 ('get_row_processes_hourly', defaults.model_processes),
#                 ('load_external_state', defaults.e_states_i),
#             ]

#             mock_funcs(defaults, mocker, __module_loc__, run_funcs_to_mock)

#         def mock_for_integration_test(defaults, mocker):
#             """Setup mocks for integration test.
#             # TODO: Ideally we should only mock IO here.

#             """

#             defaults.external_state_options = EStateOptions(
#                 file_type=FileTypes.NETCDF,
#                 variable_map=defaults.variable_map,
#                 multi_file_data=defaults.multi_file_data,
#                 preprocess_map=defaults.preprocess_map,
#                 zero_year=defaults.zero_year,
#                 data_filter=getattr(defaults, 'data_filter', None),
#                 netcdf_loader_kwargs=dict(
#                     chunks=defaults.netcdf_chunks,
#                 ),
#             )
#             funcs_to_mock = [
#                 ('dump_state_to_file_quick', None),
#                 ('dump_output_to_file_netcdf_grid', None),
#                 ('load_current_cell_state', None, lambda *_, **__: next(defaults.cell_initial_states_i)),
#                 # We mock run model so we don't need to set up the entire config/init state setup.
#                 ('run_model', None, lambda *_, **__: (next(defaults.cell_output_states_i), defaults.output_logs)),
#             ]
#             mock_funcs(defaults, mocker, __module_loc__, funcs_to_mock)

#         def default_run(defaults, kwargs={}):
#             """Default main_partial run with some presets.

#             kwargs will override the presets.
#             """
#             args = {
#                 **dict(
#                     # cell_configs=defaults.cell_configs,
#                     # grid_coords=defaults.grid_coords,
#                     # output_shape=defaults.output_shape,
#                     # external_data_file_path=defaults.data_location,
#                     # previous_hour_state_path=defaults.prev_hour_state_location,
#                     # output_directory=defaults.output_directory,
#                     # output_fields=defaults.output_fields,
#                     # external_state_options=defaults.external_state_options,
#                     cell_configs=defaults.cell_configs,
#                     prev_hour_states_loaded=defaults.prev_hour_states_loaded,
#                     grid_coords=defaults.grid_coords,
#                     external_states=defaults.external_states,
#                     output_fields=defaults.output_fields,
#                     # use_daily_loop=defaults.use_daily_loop,
#                     # parallel=defaults.parallel,
#                     # logger=defaults.logger,
#                     # parallel_args=defaults.parallel_args,
#                     # overrides=defaults.overrides,
#                     # debug=defaults.debug,
#                 ),
#                 **kwargs,
#             }
#             return main_partial(
#                 **args,
#             )

#         def _setup(self):
#             defaults.row_count = getattr(defaults, 'row_count', (defaults.end_day - defaults.start_day + 1) * 24)
#             defaults.output_logs = [
#                 {"a": i + 1} for i in range(defaults.row_count)
#             ]

#             defaults.e_state_dates = dict(
#                 start_day=defaults.start_day,
#                 end_day=defaults.end_day,
#                 start_date="01-02-2017",
#                 end_date="28-02-2017",
#                 row_count=defaults.row_count,
#                 time_string="01-02-2017",
#                 hours=[hr for _ in range(defaults.start_day - defaults.end_day + 1)
#                        for hr in range(24)][0:defaults.row_count + 1]
#             )
#             defaults.output_logs = [
#                 {"a": i + 1} for i in range(defaults.row_count)
#             ]

#         def _override_defaults(self):
#             # Set in sub test groups to override defaults
#             pass

#         def setup_for_unit_test(self):
#             TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#             defaults._override_defaults()
#             defaults._setup()

#         def setup_for_integration_test(self):
#             TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#             TestMainPartial.TestMainPartialGrid.set_integration_test_defaults(defaults)
#             defaults._setup()
#             # Expected Results
#             defaults.output_logs = [
#                 {"pody": i + 1} for i in range(defaults.row_count)
#             ]

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_runs_unit_test(self):
#             """Check that the default run works with unit test mocks."""
#             defaults.default_run()

#         @setup_test(ETestTypes.INTEGRATION_TEST)
#         def test_runs_integration_test(self):
#             """Check that the default run works with integration test mocks."""
#             defaults.default_run()

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_call_dump_state_to_file_for_each_coord(self):
#             list(defaults.default_run())

#             assert defaults.mock_dump_state_to_file_quick.call_count == len(defaults.grid_coords)
#             expected_calls = [
#                 call.dump_state_to_file(
#                     defaults.cell_output_states[i],
#                     f"{defaults.prev_hour_state_location}/{x}_{y}.state",
#                 ) for i, (x, y) in enumerate(defaults.grid_coords)

#             ]
#             assert defaults.mock_dump_state_to_file_quick.mock_calls == expected_calls

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_call_run_model(self):
#             defaults.default_run()
#             assert defaults.mock_run_model.call_count == 3
#             defaults.mock_run_model.assert_any_call(
#                 defaults.cell_initial_states[0],
#                 defaults.cell_configs[0],
#                 defaults.cell_external_states[0],
#                 defaults.model_processes,
#             )

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_call_dump_output_file_netcdf_grid_with_target_path(self):
#             time_string = defaults.e_state_dates['time_string']
#             target_file_name = f'{defaults.output_directory}/output_data_{time_string}.nc'

#             defaults.default_run()
#             calls = defaults.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
#             assert calls['target_path'] == target_file_name

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_call_dump_output_file_netcdf_grid_with_output_fields(self):
#             defaults.default_run()
#             calls = defaults.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
#             assert calls['output_fields'] == defaults.output_fields

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_call_dump_output_file_netcdf_grid_with_output_data(self):
#             defaults.output_data = {
#                 k: np.full(
#                     (defaults.grid_x_size, defaults.grid_y_size, defaults.row_count),
#                     None,
#                     dtype=np.float64) for k in defaults.output_fields
#             }

#             for x, y in defaults.grid_coords:
#                 defaults.output_data['a'][x][y] = [i + 1 for i in range(defaults.row_count)]

#             defaults.default_run()
#             calls = defaults.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
#             assert str(calls['output_data']['a'].shape) == str(
#                 (defaults.grid_x_size, defaults.grid_y_size, defaults.row_count))

#             np.testing.assert_array_equal(
#                 calls['output_data']['a'][0][0], defaults.output_data['a'][0][0])
#             np.testing.assert_array_equal(calls['output_data']['a'][0], defaults.output_data['a'][0])
#             np.testing.assert_array_equal(calls['output_data']['a'], defaults.output_data['a'])

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_call_dump_output_file_netcdf_grid_with_time_data(self):
#             # defaults.time_data = np.full(
#             #     (defaults.row_count),
#             #     None,
#             #     dtype=np.float64)

#             defaults.time_data = pd.date_range(
#                 defaults.e_state_dates['start_date'], periods=defaults.row_count, freq="1H")

#             defaults.default_run()
#             calls = defaults.mock_dump_output_to_file_netcdf_grid.call_args.kwargs

#             assert str(calls['time_data'].shape) == str((defaults.row_count,))

#             np.testing.assert_array_equal(calls['time_data'], defaults.time_data)

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_call_dump_output_file_netcdf_grid_with_location_data(self):
#             defaults.lat_data = np.full(
#                 (defaults.grid_x_size, defaults.grid_y_size),
#                 None,
#                 dtype=np.float64)
#             defaults.lon_data = np.full(
#                 (defaults.grid_x_size, defaults.grid_y_size),
#                 None,
#                 dtype=np.float64)

#             # TODO: Need to set lat and lon data
#             # for x, y in defaults.grid_coords:
#             #     defaults.lat_data[x][y] = 0
#             #     defaults.lon_data[x][y] = 0

#             defaults.default_run()
#             calls = defaults.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
#             np.testing.assert_array_equal(calls['lat_data'], defaults.lat_data)
#             np.testing.assert_array_equal(calls['lon_data'], defaults.lon_data)

#     class TestSingleHourOfData(TestMainPartialGrid):
#         """Test that the partial model can be ran when data is provided 1 hour at a time."""

#         def _setup(self):
#             defaults.row_count = 1

#         def setup_for_unit_test(self):
#             TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#             defaults._setup()

#             # Expected Results
#             defaults.output_logs = [
#                 {"a": 1},
#             ]

#         def setup_for_integration_test(self):
#             TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#             TestMainPartial.TestMainPartialGrid.set_integration_test_defaults(defaults)
#             defaults.data_location = 'examples/net_cdf/single_file_hour/inputs/demo_wrf_2017-12-27-00-00-00.nc'
#             defaults._setup()
#             # Expected Results
#             defaults.output_logs = [
#                 {"pody": 1},
#             ]

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_should_run_full_model(self):
#             output = defaults.default_run()
#             assert list(output.keys()) == defaults.output_fields
#             assert output[defaults.output_fields[0]].shape == (
#                 defaults.grid_x_size, defaults.grid_y_size, defaults.row_count)

#     class TestMultiFileTimeRange(TestMainPartialGrid):
#         """Netcdf data provided as multiple files each containing a single variable over a time range.

#         The below setup functions setup the tests for where they vary from the default setup in TestMainPartialGrid.

#         """

#         def _override_defaults(self):
#             defaults.start_day = 361
#             defaults.end_day = 364
#             defaults.multi_file_data = True
#             defaults.row_count = (defaults.end_day - defaults.start_day + 1) * 24

#         def setup_for_unit_test(self):
#             TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#             defaults._override_defaults()
#             defaults._setup()

#         def setup_for_integration_test(self):
#             TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#             TestMainPartial.TestMainPartialGrid.set_integration_test_defaults(defaults)
#             defaults.multi_file_data = True
#             defaults.data_filter = 'demo_wrf_2017-10'
#             defaults.data_location = 'examples/net_cdf/multi_file_range/inputs'

#             defaults._override_defaults()
#             # Expected Results
#             defaults.output_logs = [
#                 {"pody": i + 1} for i in range(defaults.row_count)
#             ]

#         @setup_test(ETestTypes.UNIT_TEST)
#         def test_runs_unit_test(self):
#             defaults.default_run()

#         @setup_test(ETestTypes.INTEGRATION_TEST)
#         def test_runs_integration_test(self):
#             defaults.default_run()

#     # TODO: Implement test
#     # class TestSingleFileTimeRange(TestMainPartialGrid):

#     #     def _setup(self):
#     #         defaults.row_count = 24 * 2

#     #     def setup_for_unit_test(self):
#     #         TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#     #         defaults._setup()

#     #         # Expected Results
#     #         defaults.output_logs = [
#     #             {"a": i+1} for i in range(defaults.row_count)
#     #         ]

#     #     def setup_for_integration_test(self):
#     #         TestMainPartial.TestMainPartialGrid.set_defaults(defaults)
#     #         TestMainPartial.TestMainPartialGrid.set_integration_test_defaults(defaults)
#     #         defaults.data_location = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'
#     #         defaults._setup()
#     #         # Expected Results
#     #         defaults.output_logs = [
#     #             {"pody": i+1} for i in range(defaults.row_count)
#     #         ]

#     class TestMultiHourData:
#         def test_should_run_without_errors(self):
#             # TODO: Implement test
#             pass
#             # raise NotImplementedError("Implement Test")

# TODO: Move these tests to correct func
# class TestMainGridSeq:

#     @pytest.fixture(autouse=True)
#     def _setup(defaults, mocker):
#         # Default params
#         defaults.project_dir = 'project_dir'
#         defaults.runid = 'test_run'
#         defaults.output_fields = 'pody'
#         defaults.runnotes = 'runnotes'
#         # Expected params
#         defaults.config_names = ['config']
#         defaults.configs = [f'{c}.json' for c in defaults.config_names]
#         defaults.inputs = ['input.json']
#         defaults.runs = [('final_state', 'output_logs', 'model')]
#         defaults.variable_map = MagicMock()
#         defaults.preprocess_map = MagicMock()
#         defaults.first_run_dir = f"{defaults.project_dir}/runs/{defaults.runid}/{defaults.config_names[0]}"

#         defaults.processed_configs_dir = f'{defaults.first_run_dir}/processed_configs'
#         defaults.base_config_file = f"{defaults.project_dir}/base_config.json"
#         defaults.grid_coords = [
#             [0, 0], [0, 1], [0, 2],
#             [1, 0], [1, 1], [1, 2],
#             [2, 0], [2, 1], [2, 2],
#         ]
#         defaults.input_data_files = [f"{defaults.project_dir}/inputs/{f}" for f in defaults.inputs]
#         defaults.previous_hour_state_path = f"{defaults.first_run_dir}/prev_state"
#         defaults.output_data_dir = f"{defaults.first_run_dir}/outputs_grid"
#         defaults.live_state_dir = f"{defaults.first_run_dir}/current_state"
#         defaults.logger = MagicMock()

#         # Mocks
#         defaults.spy_os_makedirs = mocker.patch('os.makedirs')
#         mocker.patch('os.rename')
#         mocker.patch('shutil.rmtree')
#         mocker.patch('builtins.open', return_value=MagicMock())

#         mocker.patch(__module_loc__ + '.Logger', return_value=defaults.logger)
#         mocker.patch(__module_loc__ + '.get_configs', return_value=defaults.configs)
#         mocker.patch(__module_loc__ + '.get_input_files_list', return_value=defaults.inputs)

#         mocker.patch(
#             __module_loc__ + '.json_loader',
#             side_effect=lambda p: defaults.variable_map if 'variable_map' in p
#             else defaults.preprocess_map if 'preprocess_map' in p
#             else None)
#         mocker.patch(
#             __module_loc__ + '.csv_loader', side_effect=lambda _: [{'x': x, 'y': y} for x, y in defaults.grid_coords])
#         mocker.patch(
#             __module_loc__ + '.generate_run_notes', return_value=[])
#         defaults.spy_main_hour = mocker.patch(__module_loc__ + '.main_hour', return_value=defaults.runs)

#     def default_run(self):
#         return main_grid_seq(
#             defaults.project_dir,
#             defaults.runid,
#             defaults.output_fields,
#             defaults.runnotes,
#             seperate_live_state=True,
#         )

#     def test_runs_without_errors(self):
#         defaults.default_run()

#     def test_calls_main_hour(self):
#         defaults.default_run()
#         defaults.spy_main_hour.assert_called_once_with(
#             processed_config_dir=f"{defaults.processed_configs_dir}",
#             external_data_row_path=defaults.input_data_files[0],
#             previous_hour_state_path=defaults.previous_hour_state_path,
#             output_data_dir=defaults.output_data_dir,
#             output_fields=defaults.output_fields,
#             logger=defaults.logger,
#             external_file_type=FileTypes.NETCDF,
#             grid_coords=defaults.grid_coords,
#             netcdf_variable_map=defaults.variable_map,
#             met_preprocess_map=defaults.preprocess_map,
#             multi_file_netcdf=False,
#             output_to_netcdf=True,
#             state_out_path=defaults.live_state_dir,
#         )

#     def test_calls_makedirs(self):
#         defaults.default_run()
#         assert defaults.spy_os_makedirs.mock_calls == [
#             call(f'{defaults.first_run_dir}/current_state', exist_ok=True),
#         ]


def generate_external_state_data(
    start_day: int,
    end_day: int,
    grid_coords: List[Tuple[int, int]],
    configs: List[Config_Shape]
) -> External_State_Shape:
    """Generates external data for hour range runs"""
    # TODO: Generate data

    DX = 3
    DY = 2
    DT = 24 * 8
    daily_sqdown = [3, 4, 5, 6, 7, 8, 9, 9.5, 9.8, 10, 9.8, 9.5, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2]
    demo_data_SWDOWN = np.array([[[daily_sqdown[i] for _ in range(DT)
                                for i in range(24)] for _ in range(DX)] for _ in range(DY)]).transpose()
    demo_data_HFX_FORCE = np.zeros((DT, DX, DY))
    demo_data_td_2m = np.ones((DT, DX, DY)) * 10
    demo_data_rh = np.zeros((DT, DX, DY)) + 0.3
    demo_data_o3 = np.ones((DT, DX, DY)) * 10
    demo_data_wspeed = np.ones((DT, DX, DY)) * 1.4
    demo_data_pres = np.ones((DT, DX, DY)) * 101
    demo_data_RAINNC = np.ones((DT, DX, DY)) * 4
    demo_data_SNOWH = np.zeros((DT, DX, DY))
    time_data = pd.date_range("2017-12-27", periods=DT, freq="1H")
    # Julian days
    dd = [dd for dd in range(start_day, end_day + 1) for _ in range(24)]
    hr = [hr for _ in range(start_day, end_day + 1) for hr in range(24)]
    external_states_output = []
    for (x, y), config in zip(grid_coords, configs):

        e_state_init_processes = external_state_init_processes(config_met=config.Met)
        external_states_output.append(setup_external_state_simple(External_State_Shape(
            time=time_data,
            dd=dd,
            hr=hr,
            Ts_C=demo_data_td_2m[:, x, y],
            P=demo_data_pres[:, x, y],
            PAR=demo_data_SWDOWN[:, x, y],
            precip=demo_data_RAINNC[:, x, y],
            RH=demo_data_rh[:, x, y],
            u=demo_data_wspeed[:, x, y],
            O3=demo_data_o3[:, x, y],
            Hd=demo_data_HFX_FORCE[:, x, y],
            snow_depth=demo_data_SNOWH[:, x, y],
        ), config, e_state_init_processes))
    return external_states_output


demo_config = config_loader(
    'examples/net_cdf/single_file_hour/configs/bangor_wheat.json',
    'examples/net_cdf/single_file_hour/base_config.json',
)


class TestMainPartialGrid:

    def set_defaults(self, config=demo_config):
        defaults = MagicMock()
        defaults.grid_coords = [(0, 0), (1, 0), (2, 1)]
        defaults.grid_x_size = 5  # this is range of x coords
        defaults.grid_y_size = 6  # this is range of y coords
        defaults.start_day = 1  # "01-01-2017"
        defaults.end_day = 4  # "04-01-2017"
        defaults.start_date = "01-01-2017"
        defaults.end_date = "04-01-2017"
        defaults.cell_configs = [demo_config for _ in defaults.grid_coords]
        defaults.prev_hour_states_loaded = [
            deepcopy(example_state) for _ in defaults.grid_coords]

        defaults.row_count = 24 * (defaults.end_day - defaults.start_day)
        defaults.hours = [hr for _ in range(
            defaults.start_day - defaults.end_day) for hr in range(24)]
        defaults.time_string = defaults.start_date
        defaults.lat_data = np.arange(defaults.grid_x_size)
        defaults.lon_data = np.arange(defaults.grid_y_size)
        for i, (c, (x, y)) in enumerate(zip(defaults.cell_configs, defaults.grid_coords)):
            c.Location.lat = defaults.lat_data[x]
            c.Land_Cover.phenology_options.latitude = defaults.lat_data[x]
            c = setup_config(c)
            defaults.cell_configs[i] = c
        defaults.time_data = pd.date_range(
            defaults.start_date, periods=defaults.row_count, freq="1H")
        defaults.input_shape = (defaults.grid_x_size, defaults.grid_y_size, defaults.row_count)
        defaults.external_state_preprocessed = generate_external_state_data(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            grid_coords=defaults.grid_coords,
            configs=defaults.cell_configs,
        )

        [
            start_day,
            end_day,
            start_date,
            end_date,
            row_count,
            time_string,
            hours,
        ] = get_date_bounds_from_ext_data(
            defaults.external_state_preprocessed[0],
        )
        defaults.start_day = start_day
        defaults.end_day = end_day
        defaults.start_date = start_date
        defaults.end_date = end_date
        defaults.row_count = row_count
        defaults.hours = hours

        external_state_iter = ExternalStateIterable(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            start_date=defaults.start_date,
            end_date=defaults.end_date,
            row_count=defaults.row_count,
            hours=defaults.hours,
            time_string=defaults.time_string,
            lat_data=defaults.lat_data,
            lon_data=defaults.lon_data,
            time_data=defaults.time_data,
            input_shape=defaults.input_shape,
            external_state_preprocessed=iter(defaults.external_state_preprocessed),
        )

        defaults.external_states = external_state_iter
        defaults.output_fields = ['pody', 'par']
        demo_config.output.fields = defaults.output_fields
        return defaults

    def default_run(self, defaults, kwargs={}):
        """Default main_partial run with some presets.

        kwargs will override the presets.
        """
        args = {
            **dict(
                cell_configs=defaults.cell_configs,
                prev_hour_states_loaded=defaults.prev_hour_states_loaded,
                grid_coords=defaults.grid_coords,
                external_states=defaults.external_states,
                output_fields=defaults.output_fields,
                use_daily_loop=True,
            ),
            **kwargs,
        }
        return main_partial(
            **args,
        )

    def test_runs_without_errors(self):
        defaults = self.set_defaults()
        self.default_run(defaults)

    def test_returns_model_outputs(self):
        defaults = self.set_defaults()
        output = self.default_run(defaults)
        coords, data, state, output_fields = output[0]
        assert coords == list(defaults.grid_coords[0])
        print(data)
        assert len(data) == len(defaults.output_fields)
        assert output_fields == defaults.output_fields

    def test_returns_model_outputs_for_all_coords(self):
        defaults = self.set_defaults()
        output = self.default_run(defaults)
        assert len(output) == len(defaults.grid_coords)

    def test_returns_model_outputs_for_all_coords_and_fields(self):
        defaults = self.set_defaults()
        output = self.default_run(defaults)
        assert len(output) == len(defaults.grid_coords)
        assert len(output[0][1]) == len(defaults.output_fields)

    def test_returns_model_outputs_for_all_coords_and_fields_and_time(self):
        defaults = self.set_defaults()
        output = self.default_run(defaults)
        coords, data, state, output_fields = output[0]
        assert len(output) == len(defaults.grid_coords)
        assert len(data) == len(defaults.output_fields)
        assert len(data[0]) == defaults.row_count

    def test_should_run_ok_with_offset_date(self):
        config = demo_config
        config.Location.hr_offset = 3
        config.Location.crop_to_day_start = True
        config.Location.crop_to_day_end = True
        defaults = self.set_defaults(config)
        # We crop remaining first day
        assert len(defaults.external_state_preprocessed[0].hr) == (
            (4 - 1) * 24) == defaults.row_count
        output = self.default_run(defaults)
        coords, data, state, output_fields = output[0]
        assert len(data) == len(defaults.output_fields)
        assert len(data[0]) == defaults.row_count
        assert output_fields[1] == 'par'
        assert data[1][0] == 0


class TestGridExternalState:

    def set_defaults(self, config=demo_config):
        class Foo:
            pass
        defaults = Foo()
        defaults.config = config
        defaults.grid_coords = [(0, 0), (1, 0), (2, 1)]
        defaults.grid_x_size = 5  # this is range of x coords
        defaults.grid_y_size = 6  # this is range of y coords
        defaults.start_day = 1  # "01-01-2017"
        defaults.end_day = 4  # "04-01-2017"
        defaults.start_date = "01-01-2017"
        defaults.end_date = "04-01-2017"
        defaults.cell_configs = [demo_config for _ in defaults.grid_coords]
        defaults.prev_hour_states_loaded = [deepcopy(example_state) for _ in defaults.grid_coords]

        defaults.row_count = 24 * (defaults.end_day - defaults.start_day)
        defaults.hours = [hr for _ in range(
            defaults.start_day - defaults.end_day) for hr in range(24)]
        defaults.time_string = defaults.start_date
        defaults.lat_data = np.arange(defaults.grid_x_size)
        defaults.lon_data = np.arange(defaults.grid_y_size)
        for i, (c, (x, y)) in enumerate(zip(defaults.cell_configs, defaults.grid_coords)):
            c.Location.lat = defaults.lat_data[x]
            c.Land_Cover.phenology_options.latitude = defaults.lat_data[x]
            c = setup_config(c)
            defaults.cell_configs[i] = c
        defaults.time_data = pd.date_range(
            defaults.start_date, periods=defaults.row_count, freq="1H")
        defaults.input_shape = (defaults.grid_x_size, defaults.grid_y_size, defaults.row_count)
        return defaults

    def test_can_setup_external_grid_data(self):
        config = demo_config
        config.Location.hr_offset = 0
        config.Location.crop_to_day_start = False
        config.Location.crop_to_day_end = False
        defaults = self.set_defaults()
        print(defaults.config.Location.hr_offset)
        external_state_preprocessed = generate_external_state_data(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            grid_coords=defaults.grid_coords,
            configs=defaults.cell_configs,
        )
        external_state_iter = ExternalStateIterable(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            start_date=defaults.start_date,
            end_date=defaults.end_date,
            row_count=defaults.row_count,
            hours=defaults.hours,
            time_string=defaults.time_string,
            lat_data=defaults.lat_data,
            lon_data=defaults.lon_data,
            time_data=defaults.time_data,
            input_shape=defaults.input_shape,
            external_state_preprocessed=iter(external_state_preprocessed),
        )
        assert external_state_iter is not None
        external_state_processed = next(external_state_iter)
        assert external_state_processed.hr == [
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23
        ]

        assert external_state_processed.dd == [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        ]
        assert external_state_processed.PAR[0] == 3.0

    def test_can_offset_hours(self):
        config = demo_config
        config.Location.hr_offset = 3
        config.Location.crop_to_day_start = False
        config.Location.crop_to_day_end = False
        defaults = self.set_defaults(config)
        external_state_preprocessed = generate_external_state_data(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            grid_coords=defaults.grid_coords,
            configs=defaults.cell_configs,
        )
        external_state_iter = ExternalStateIterable(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            start_date=defaults.start_date,
            end_date=defaults.end_date,
            row_count=defaults.row_count,
            hours=defaults.hours,
            time_string=defaults.time_string,
            lat_data=defaults.lat_data,
            lon_data=defaults.lon_data,
            time_data=defaults.time_data,
            input_shape=defaults.input_shape,
            external_state_preprocessed=iter(external_state_preprocessed),
        )
        external_state_processed = next(external_state_iter)
        assert external_state_processed.hr == [
            3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2
        ]
        assert external_state_processed.dd == [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5,
        ]

    def test_can_offset_hours_with_crop(self):
        config = demo_config
        config.Location.hr_offset = 3
        config.Location.crop_to_day_start = True
        config.Location.crop_to_day_end = True
        defaults = self.set_defaults(config)
        external_state_preprocessed = generate_external_state_data(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            grid_coords=defaults.grid_coords,
            configs=defaults.cell_configs,
        )
        external_state_iter = ExternalStateIterable(
            start_day=defaults.start_day,
            end_day=defaults.end_day,
            start_date=defaults.start_date,
            end_date=defaults.end_date,
            row_count=defaults.row_count,
            hours=defaults.hours,
            time_string=defaults.time_string,
            lat_data=defaults.lat_data,
            lon_data=defaults.lon_data,
            time_data=defaults.time_data,
            input_shape=defaults.input_shape,
            external_state_preprocessed=iter(external_state_preprocessed),
        )
        external_state_processed = next(external_state_iter)
        assert external_state_processed.hr == [
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
        ]
        assert external_state_processed.dd == [
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        ]
        # TODO: Test that temperature data etc is also cropped correctly
        assert external_state_processed.PAR[0] == 0

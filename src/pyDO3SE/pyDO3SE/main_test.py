"""Unit tests for main.py ."""
from copy import deepcopy
import numpy as np
import pandas as pd
from itertools import cycle
import pytest

from unittest.mock import MagicMock, Mock, call
from pyDO3SE.Config.config_loader import config_loader

from pyDO3SE.External_State.external_state_loader import EStateOptions, FileTypes
from pyDO3SE.Model_State.demo_data import example_state
from pyDO3SE.util.test_utils import (
    mock_funcs, setup_test, TestTypes,
)

from .main import (
    ProjectPaths,
    RunPaths,
    main,
    main_partial_initialize,
    main_partial,
)


__module_loc__ = __package__ + '.main'


def Obj(name: str, *args, **kwargs) -> object:
    argsAsKwargs = {k._extract_mock_name(): k for k in args}
    return type(name, (object,), {**kwargs, **argsAsKwargs})()


# == Mocked return values == #
mock_config = MagicMock()
mock_initial_state = MagicMock()
mock_model_processes = MagicMock()
mock_external_state = MagicMock()
mock_parameters = MagicMock()
mock_state = MagicMock()
mock_state_out = MagicMock()
mock_run_processes = MagicMock(return_value=mock_state_out)
mock_loaded_run_files = MagicMock()
mock_state_logs = [
    {"a": 1, "b": 2},
    {"a": 3, "b": 4},
    {"a": 5, "b": 6},
]
mock_run_notes = MagicMock()
mock_processes = [1, 2, 3]
mock_overrides = MagicMock()
mock_overrides.debug = False
mock_process_runner = Mock(
    run_processes=mock_run_processes,
    state_logs=mock_state_logs,
)
config_location = "demo_config_location"
mock_data_location = "demo_data_location"
mock_base_config_file = "demo_base_config.json"
mock_initial_state_file = None
mock_project_paths = ProjectPaths(config_location, mock_data_location,
                                  base_config_path=mock_base_config_file)
mock_run_paths = RunPaths()

mock_processes = [1, 2, 3]


# @pytest.fixture(autouse=True)
# def mock_imports(mocker: MockerFixture) -> object:
#     """Mock imported functions."""
#     mock_setup_model = mocker.patch(__module_loc__ + '.setup_model', return_value=[
#         mock_config,
#         mock_external_state,
#         mock_initial_state,
#         mock_model_processes,
#     ])
#     return Obj('MockedImports', mock_setup_model)


# class TestRunModel:
#     @pytest.fixture(autouse=True)
#     def _setup(self, mocker) -> None:
#         self.mock_setup_state = mocker.patch(
#             __module_loc__ + '.setup_initial_state', return_value=(mock_state))
#         self.mock_get_processes = mocker.patch(
#             __module_loc__ + '.full_model_processes', return_value=(mock_processes))

#     def test_creates_a_process_runner_with_config_input(self, mock_imports):
#         out = run_model(mock_config, mock_external_state)
#         mock_imports.ProcessRunner.assert_called_once_with(mock_config, mock_external_state)
#         self.mock_setup_state.assert_called_once_with(mock_process_runner, Main_Overrides())
#         self.mock_get_processes.assert_called_once_with(
#             mock_config, 0, 365)
#         assert mock_process_runner.run_processes.call_count == 1
#         assert out[0] == mock_state_out
#         assert out[1] == mock_state_logs


class TestMain:

    @pytest.fixture(autouse=True)
    def _setup(self, mocker):
        self.mock_run_model = mocker.patch(
            __module_loc__ + '.run_model_daily', return_value=(mock_state_out, mock_state_logs))
        # self.mock_setup_config = mocker.patch(
        #     __module_loc__ + '.setup_config', return_value=mock_config)
        self.mock_setup_model = mocker.patch(
            __module_loc__ + '.setup_model', return_value=[
                mock_config,
                mock_external_state,
                mock_initial_state,
                mock_model_processes,
            ])
        mocker.patch(__module_loc__ + '.load_run_files', return_value=mock_loaded_run_files)
        mocker.patch(__module_loc__ + '.Main_Overrides', return_value=mock_overrides)
        mocker.patch(__module_loc__ + '.generate_run_notes', return_value=mock_run_notes)
        mocker.patch(__module_loc__ + '.export_output')

        # self.mock_setup_external_state = mocker.patch(
        #     __module_loc__ + '.setup_external_state', return_value=mock_external_state)

    def test_should_receive_config_and_data_location_and_run_model(self, mocker):
        final_state, output_logs, config, initial_state = main(
            project_paths=mock_project_paths,
            run_paths=mock_run_paths,
        )
        self.mock_setup_model.assert_called_once_with(
            config_in=mock_loaded_run_files.config,
            state_in=mock_loaded_run_files.state,
            data_location=mock_project_paths.input_data_dir,
            overrides=mock_overrides,
        )

        assert final_state == mock_state_out
        assert output_logs == mock_state_logs

        self.mock_run_model.assert_called_once_with(
            mock_initial_state,
            mock_config,
            mock_external_state,
            mock_model_processes,
            DEBUG_MODE=mock_overrides.debug,
        )

# @Deprecated
# class TestMultiPass:
#     @pytest.fixture(autouse=True)
#     def _setup(self, mocker) -> None:
#         self.mock_run_model = mocker.patch(
#             __module_loc__ + '.run_model', return_value=(mock_state_logs))
#         self.mock_setup_config = mocker.patch(
#             __module_loc__ + '.setup_config', return_value=mock_config)
#         self.mock_setup_external_state = mocker.patch(
#             __module_loc__ + '.setup_external_state', return_value=(mock_external_state, 0, 365))
#         self.mock_setup_initial_state = mocker.patch(
#             __module_loc__ + '.setup_initial_state', return_value=mock_initial_state)
#         self.mock_setup_model_processses = mocker.patch(
#             __module_loc__ + '.setup_model_processes', return_value=mock_model_processes)
#         mocker.patch(__module_loc__ + '.Main_Overrides', return_value=mock_overrides)


#     @pytest.fixture(autouse=True)
#     def configs_input(self):
#         return ['a', 'b', 'c']

#     def test_should_raise_value_error_if_invalid_config_supplied(self):
#         with pytest.raises(ValueError) as ve:
#             data_location = "demo_data_location"
#             configs_input = [1, 2, 3]
#             multi_run(configs_input, data_location)
#         assert 'Config must be' in str(ve.value)

#     def test_should_run_config_setup_for_each_config(self, configs_input):
#         multi_run(configs_input, mock_data_location)
#         assert self.mock_setup_config.call_count == len(configs_input)
#         assert self.mock_setup_config.call_args_list == [
#             call(configs_input[0], mock_overrides),
#             call(configs_input[1], mock_overrides),
#             call(configs_input[2], mock_overrides),
#         ]

#     def test_should_setup_external_state_for_each_config(self, configs_input):
#         multi_run(configs_input, mock_data_location)

#         assert self.mock_setup_external_state.call_count == len(configs_input)
#         # TODO: Should really check that each mock_config is different
#         assert self.mock_setup_external_state.call_args_list == [
#             call(mock_config, mock_data_location, mock_overrides),
#             call(mock_config, mock_data_location, mock_overrides),
#             call(mock_config, mock_data_location, mock_overrides),
#         ]

#     def test_should_run_a_model_run_for_each_config_input(self, configs_input):
#         multi_run(configs_input, mock_data_location)

#         assert self.mock_run_model.call_count == len(configs_input)
#         assert self.mock_run_model.call_args_list == [
#             call(mock_initial_state, mock_config, mock_external_state, mock_model_processes),
#             call(mock_initial_state, mock_config, mock_external_state, mock_model_processes),
#             call(mock_initial_state, mock_config, mock_external_state, mock_model_processes),
#         ]

#     def test_should_return_true_when_complete(self, configs_input):
#         output = multi_run(configs_input, mock_data_location)
#         assert output is True

#     @pytest.mark.skip(reason="Not implemented!")
#     def test_should_create_output_folder_for_each_config(self):
#         assert False

#     @pytest.mark.skip(reason="Not implemented!")
#     def test_should_export_results_for_each_config(self):
#         assert False

#     @pytest.mark.skip(reason="Not implemented!")
#     def test_should_run_additional_analysis_processes_on_each_output(self):
#         assert False


class TestMainhourInitialize:

    @pytest.fixture(autouse=True)
    def _setup(self, mocker):
        # Default params
        self.processed_config_dir = 'processed_config_dir'
        self.state_out_path = 'state_out_path'

        self.e_state_overrides_file_path = 'e_state_overrides_file_path'
        self.e_state_overrides_field_map = {}
        self.logger = MagicMock()
        self.config = MagicMock()
        self.state = MagicMock()
        self.init_state = MagicMock()
        self.grid_coords = MagicMock()
        self.init_config = MagicMock()

        self.mock_initialize_grid_configs = mocker.patch(
            __module_loc__ + '.initialize_grid_configs',
        )
        mocker.patch(__module_loc__ + '.Logger', return_value=self.logger)
        mocker.patch(__module_loc__ + '.setup_config', return_value=self.init_config)
        mocker.patch(__module_loc__ + '.setup_initial_state', return_value=self.init_state)
        mocker.patch(__module_loc__ + '.dump_state_to_file_quick')

    def default_run(self):
        return main_partial_initialize(
            config=self.config,
            state=self.state,
            processed_config_dir=self.processed_config_dir,
            state_out_path=self.state_out_path,
            e_state_overrides_file_path=self.e_state_overrides_file_path,
            e_state_overrides_field_map=self.e_state_overrides_field_map,
            grid_coords=self.grid_coords,
            logger=self.logger,
        )

    def test_should_call_initialize_grid_configs(self):
        self.default_run()
        self.mock_initialize_grid_configs.assert_called_once_with(
            self.init_config,
            self.processed_config_dir,
            self.grid_coords,
            self.e_state_overrides_file_path,
            self.e_state_overrides_field_map,
            self.logger,
        )

# TODO: Check if these tests should move to MainPartial
# class TestMainHour:

#     def setup(self, mocker):
#         # Default Params
#         self.config_location = "demo_config_location"
#         self.data_location = "demo_data_location"
#         self.prev_hour_state_location = "demo_state_location"
#         self.config_name = "demo_config"
#         self.processed_config_dir = 'processed_config_dir'
#         self.config_names = [self.config_name]
#         self.run_id = 0
#         self.x = 0
#         self.y = 0
#         self.grid_coords = [(0, 0), (1, 0)]

#         # Mocks
#         self.mock_run_model = mocker.patch(
#             __module_loc__ + '.run_model', return_value=(mock_state_out, [mock_state_logs]))

#         self.mock_setup_model_hour = mocker.patch(
#             __module_loc__ + '.setup_model_single_hour_grid', return_value=[[
#                 mock_config,
#                 mock_external_state,
#                 mock_initial_state,
#                 mock_model_processes,
#                 0,  # x
#                 0,  # y
#             ]])
#         self.mock_dump_state_to_file_quick = mocker.patch(
#             __module_loc__ + '.dump_state_to_file_quick'
#         )
#         mock_overrides.state_out_path = None
#         mocker.patch(__module_loc__ + '.Main_Overrides', return_value=mock_overrides)

#     @pytest.fixture(autouse=True)
#     def _setup(self, mocker):
#         self.setup(mocker)

#     def default_run(self, kwargs={}):
#         args = {
#             **{
#                 'config_names': self.config_names,
#                 'processed_config_dir': self.config_location,
#                 'external_data_row_path': self.data_location,
#                 'previous_hour_state_path': self.prev_hour_state_location,
#                 'state_out_path': None,
#             },
#             **kwargs,
#         }
#         return main_hour(
#             **args,
#         )

#     @pytest.mark.skip(reason="Broken as only works if supplied with grid info")
#     def test_should_run_a_single_hour_of_data(self):
#         output = next(self.default_run())
#         x = 0
#         y = 0
#         self.mock_setup_model_hour.assert_called_once_with(
#             self.config_location,
#             self.data_location,
#             self.config_name,
#             self.prev_hour_state_location,
#             self.run_id,
#             mock_overrides,
#         )
#         self.mock_dump_state_to_file_quick.assert_called_once_with(
#             mock_state_out,
#             f"{self.prev_hour_state_location}/{self.run_id}_{x}_{y}.state",
#         )

#         assert output[0] == mock_state_out
#         assert output[1] == mock_state_logs

#         self.mock_run_model.assert_called_once_with(
#             mock_initial_state,
#             mock_config,
#             mock_external_state,
#             mock_model_processes,
#         )

#     class TestGrid:

#         def default_run(self, kwargs={}):
#             return TestMainHour.default_run(self, {
#                 **dict(grid_coords=self.grid_coords),
#                 **kwargs,
#             })

#         @pytest.fixture(autouse=True)
#         def _setup(self, mocker):
#             TestMainHour.setup(self, mocker)
#             mock_overrides.grid_coords = self.grid_coords
#             self.mock_setup_model_hour = mocker.patch(
#                 __module_loc__ + '.setup_model_single_hour_grid', return_value=[[
#                     mock_config,
#                     mock_external_state,
#                     mock_initial_state,
#                     mock_model_processes,
#                     x,
#                     y,
#                 ] for x, y in self.grid_coords])

#         def test_runs_ok(self):
#             next(TestMainHour.default_run(self))

#         def test_should_run_a_single_hour(self):
#             output = next(self.default_run())
#             assert output[0] == mock_state_out
#             assert output[1] == [mock_state_logs]

#         def test_should_call_setup_model_single_hour_grid(self):
#             next(self.default_run())

#             self.mock_setup_model_hour.assert_called_once_with(
#                 self.config_location,
#                 self.data_location,
#                 self.prev_hour_state_location,
#                 overrides=mock_overrides,
#             )

#         def test_should_call_dump_state_to_file_for_each_coord(self):
#             list(self.default_run())
#             assert self.mock_dump_state_to_file_quick.call_count == len(self.grid_coords)
#             expected_calls = [
#                 call.dump_state_to_file(
#                     mock_state_out,
#                     f"{self.prev_hour_state_location}/{x}_{y}.state",
#                 ) for x, y in self.grid_coords

#             ]
#             assert self.mock_dump_state_to_file_quick.mock_calls == expected_calls

#         def test_should_call_run_model(self):
#             next(self.default_run())
#             self.mock_run_model.assert_called_once_with(
#                 mock_initial_state,
#                 mock_config,
#                 mock_external_state,
#                 mock_model_processes,
#             )

#         # def test_should_return_the_state_after_run(self):
#         #     pass
#         # def test_should_return_the_logged_fields_to_save(self):
#         #     pass
#         # def test_should_save_state_to_output_location(self):
#         #     pass


class TestMainPartial:

    class TestMainPartialPointLocation:
        """Test that we can run the DO3SE model as a partial run on single location csv data.

        Should handle the following scenarios:
        - csv file provided as a single file covering all hours
        - multiple csv files that can be concatenated

        """
        # TODO: Implement tests
        pass

    class TestMainPartialGrid:
        """Test that we can run the DO3SE model as a partial run on netcdf gridded data.

        Should handle the following scenarios:
        - netcdf data provided as a file per hour where a file contains all required vars
        - netcdf data provided as a file covering a time range which contains all required vars
        - netcdf data provided as multiple files each containing a single variable over a time range

        """

        def set_defaults(self):
            """Set the default parameters for a fully mocked run. """
            # Default Input Params
            self.data_location = "demo_data_location"
            self.prev_hour_state_location = "demo_state_location"
            self.config_name = "demo_config"
            self.output_directory = "output_directory"
            self.output_fields = ['a']
            self.config_names = [self.config_name]
            self.run_id = 0
            self.grid_coords = [(0, 0), (1, 0), (4, 5)]
            self.grid_x_size = 5  # this is range of x coords
            self.grid_y_size = 6  # this is range of y coords
            self.output_shape = (self.grid_x_size,self.grid_y_size)
            self.cell_configs = [MagicMock() for _ in self.grid_coords]
            self.cell_initial_states = [MagicMock() for _ in self.grid_coords]
            self.cell_output_states = [MagicMock() for _ in self.grid_coords]
            self.cell_external_states = [MagicMock() for _ in self.grid_coords]
            self.external_state_options = MagicMock()
            self.external_state_init_processes = (MagicMock())
            self.model_processes = [MagicMock]
            self.mock_overrides = mock_overrides
            self.mock_overrides.state_out_path = None
            self.mock_overrides.grid_coords = self.grid_coords
            self.start_day = 274
            self.end_day = 274
            self.row_count = 1
            self.e_state_dates = dict(
                # TODO: replace date strings with correct values below
                start_day=self.start_day,
                end_day=self.end_day,
                start_date="01-02-2017",
                end_date="28-02-2017",
                row_count=self.row_count,
                time_string="01-02-2017",
                hours=[hr for _ in range(self.start_day - self.end_day) for hr in range(24)]
            )
            # Expected Results
            self.output_data = {
                "a": [row['a'] for row in mock_state_logs],
                "b": [row['b'] for row in mock_state_logs],
            }

            self.output_logs = [
                {"a": i + 1} for i in range(self.row_count)
            ]

            # Iterators
            self.cell_output_states_i = cycle(self.cell_output_states)
            self.cell_initial_states_i = cycle(self.cell_initial_states)
            self.e_states_i = cycle(self.cell_external_states)

        def set_integration_defaults(self):
            """Setup the defaults when only mocking model run and outputs. """
            # Full run without any mocks
            self.variable_map = {
                "time": "XTIME",
                "Ts_C": "td_2m",
                "P": "pres",
                "PAR": "SWDOWN",
                "precip": "RAINNC",
                "RH": "rh",
                "u": "wspeed",
                "O3": "o3",
                "Hd": "HFX_FORCE",
                "snow_depth": "SNOWH"
            }
            self.data_location = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'
            self.multi_file_data = False
            self.preprocess_map = {}
            self.zero_year = 2017
            self.netcdf_chunks = None
            self.prev_hour_state_location = "demo_state_location"
            self.output_fields = ['pody']

            demo_config = config_loader(
                'examples/net_cdf/wrfchem/configs/bangor_wheat.json',
                'examples/net_cdf/wrfchem/base_config.json',
            )

            self.cell_configs = [demo_config for _ in self.grid_coords]

            self.cell_states = [deepcopy(example_state) for _ in self.grid_coords]
            self.cell_states_i = cycle(self.cell_states)

            self.mock_state_out = MagicMock()
            self.output_logs = [
                {"pody": 1},
            ]


        def mock_for_unit_test(self, mocker):
            """Setup mocks for a unit test."""
            funcs_to_mock = [
                ('get_date_bounds_from_ext_data', self.e_state_dates.values()),
                ('run_model', None, lambda *_, **__: (next(self.cell_output_states_i), self.output_logs)),
                ('dump_state_to_file_quick', None),
                ('Main_Overrides', mock_overrides),
                ('dump_output_to_file_netcdf_grid', None),
                ('external_state_init_processes', self.external_state_init_processes),
                ('load_current_cell_state', None, lambda *_, **__: next(self.cell_initial_states_i)),
                ('get_row_processes_hourly', self.model_processes),
                ('load_external_state', self.e_states_i),
            ]

            mock_funcs(self, mocker, __module_loc__, funcs_to_mock)

        def mock_for_integration_test(self, mocker):
            """Setup mocks for integration test.
            # TODO: Ideally we should only mock IO here.

            """

            self.external_state_options = EStateOptions(
                file_type=FileTypes.NETCDF,
                variable_map=self.variable_map,
                multi_file_data=self.multi_file_data,
                preprocess_map=self.preprocess_map,
                zero_year=self.zero_year,
                netcdf_chunks=self.netcdf_chunks,
                data_filter=getattr(self, 'data_filter', None),
            )
            funcs_to_mock = [
                ('dump_state_to_file_quick', None),
                ('dump_output_to_file_netcdf_grid', None),
                ('load_current_cell_state', None, lambda *_, **__: next(self.cell_initial_states_i)),
                # We mock run model so we don't need to set up the entire config/init state setup.
                ('run_model', None, lambda *_, **__: (next(self.cell_output_states_i), self.output_logs)),
            ]
            mock_funcs(self, mocker, __module_loc__, funcs_to_mock)

        def default_run(self, kwargs={}):
            """Default main_partial run with some presets.

            kwargs will override the presets.
            """
            args = {
                **dict(
                    cell_configs=self.cell_configs,
                    grid_coords=self.grid_coords,
                    output_shape=self.output_shape,
                    external_data_file_path=self.data_location,
                    previous_hour_state_path=self.prev_hour_state_location,
                    output_directory=self.output_directory,
                    output_fields=self.output_fields,
                    external_state_options=self.external_state_options,
                ),
                **kwargs,
            }
            return main_partial(
                **args,
            )

        def _setup(self):
            self.row_count = getattr(self, 'row_count', (self.end_day - self.start_day + 1) * 24)
            self.output_logs = [
                {"a": i + 1} for i in range(self.row_count)
            ]

            self.e_state_dates = dict(
                start_day=self.start_day,
                end_day=self.end_day,
                start_date="01-02-2017",
                end_date="28-02-2017",
                row_count=self.row_count,
                time_string="01-02-2017",
                hours=[hr for _ in range(self.start_day - self.end_day) for hr in range(24)]
            )
            self.output_logs = [
                {"a": i + 1} for i in range(self.row_count)
            ]


        def _override_defaults(self):
            # Set in sub test groups to override defaults
            pass

        def setup_for_unit_test(self):
            TestMainPartial.TestMainPartialGrid.set_defaults(self)
            self._override_defaults()
            self._setup()

        def setup_for_integration_test(self):
            TestMainPartial.TestMainPartialGrid.set_defaults(self)
            TestMainPartial.TestMainPartialGrid.set_integration_defaults(self)
            self._setup()
            # Expected Results
            self.output_logs = [
                {"pody": i + 1} for i in range(self.row_count)
            ]

        @setup_test(TestTypes.UNIT_TEST)
        def test_runs_unit_test(self):
            """Check that the default run works with unit test mocks."""
            next(self.default_run())

        @setup_test(TestTypes.INTEGRATION_TEST)
        def test_runs_integration_test(self):
            """Check that the default run works with integration test mocks."""
            next(self.default_run())

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_call_dump_state_to_file_for_each_coord(self):
            list(self.default_run())

            assert self.mock_dump_state_to_file_quick.call_count == len(self.grid_coords)
            expected_calls = [
                call.dump_state_to_file(
                    self.cell_output_states[i],
                    f"{self.prev_hour_state_location}/{x}_{y}.state",
                ) for i, (x, y) in enumerate(self.grid_coords)

            ]
            assert self.mock_dump_state_to_file_quick.mock_calls == expected_calls

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_call_run_model(self):
            next(self.default_run())

            self.mock_run_model.assert_called_once_with(
                self.cell_initial_states[0],
                self.cell_configs[0],
                self.cell_external_states[0],
                self.model_processes,
            )

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_call_dump_output_file_netcdf_grid_with_target_path(self):
            time_string = self.e_state_dates['time_string']
            target_file_name = f'{self.output_directory}/output_data_{time_string}.nc'

            list(self.default_run())
            calls = self.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
            assert calls['target_path'] == target_file_name

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_call_dump_output_file_netcdf_grid_with_output_fields(self):
            list(self.default_run())
            calls = self.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
            assert calls['output_fields'] == self.output_fields

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_call_dump_output_file_netcdf_grid_with_output_data(self):
            self.output_data = {
                k: np.full(
                    (self.grid_x_size, self.grid_y_size, self.row_count),
                    None,
                    dtype=np.float64) for k in self.output_fields
            }

            for x, y in self.grid_coords:
                self.output_data['a'][x][y] = [i + 1 for i in range(self.row_count)]

            list(self.default_run())
            calls = self.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
            assert str(calls['output_data']['a'].shape) == str(
                (self.grid_x_size, self.grid_y_size, self.row_count))

            np.testing.assert_array_equal(calls['output_data']['a'][0][0], self.output_data['a'][0][0])
            np.testing.assert_array_equal(calls['output_data']['a'][0], self.output_data['a'][0])
            np.testing.assert_array_equal(calls['output_data']['a'], self.output_data['a'])

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_call_dump_output_file_netcdf_grid_with_time_data(self):
            # self.time_data = np.full(
            #     (self.row_count),
            #     None,
            #     dtype=np.float64)

            self.time_data = pd.date_range(
                self.e_state_dates['start_date'], periods=self.row_count, freq="1H")

            list(self.default_run())
            calls = self.mock_dump_output_to_file_netcdf_grid.call_args.kwargs

            assert str(calls['time_data'].shape) == str((self.row_count,))

            np.testing.assert_array_equal(calls['time_data'], self.time_data)

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_call_dump_output_file_netcdf_grid_with_location_data(self):
            self.lat_data = np.full(
                (self.grid_x_size, self.grid_y_size),
                None,
                dtype=np.float64)
            self.lon_data = np.full(
                (self.grid_x_size, self.grid_y_size),
                None,
                dtype=np.float64)

            for x, y in self.grid_coords:
                self.lat_data[x][y] = 0
                self.lon_data[x][y] = 0

            list(self.default_run())
            calls = self.mock_dump_output_to_file_netcdf_grid.call_args.kwargs
            np.testing.assert_array_equal(calls['lat_data'], self.lat_data)
            np.testing.assert_array_equal(calls['lon_data'], self.lon_data)

    class TestSingleHourOfData(TestMainPartialGrid):
        """Test that the partial model can be ran when data is provided 1 hour at a time."""

        def _setup(self):
            self.row_count = 1

        def setup_for_unit_test(self):
            TestMainPartial.TestMainPartialGrid.set_defaults(self)
            self._setup()

            # Expected Results
            self.output_logs = [
                {"a": 1},
            ]

        def setup_for_integration_test(self):
            TestMainPartial.TestMainPartialGrid.set_defaults(self)
            TestMainPartial.TestMainPartialGrid.set_integration_defaults(self)
            self.data_location = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'
            self._setup()
            # Expected Results
            self.output_logs = [
                {"pody": 1},
            ]

        @setup_test(TestTypes.UNIT_TEST)
        def test_should_run_a_single_hour(self):
            output = next(self.default_run())

            final_state, output_logs, model = output
            assert final_state == self.cell_output_states[0]
            assert len(output_logs) == 1
            assert list(output_logs[0].keys()) == self.output_fields

    class TestMultiFileTimeRange(TestMainPartialGrid):
        """Netcdf data provided as multiple files each containing a single variable over a time range.

        The below setup functions setup the tests for where they vary from the default setup in TestMainPartialGrid.

        """


        def _override_defaults(self):
            self.start_day = 361
            self.end_day = 364
            self.multi_file_data = True
            self.row_count = (self.end_day - self.start_day + 1) * 24

        def setup_for_unit_test(self):
            TestMainPartial.TestMainPartialGrid.set_defaults(self)
            self._override_defaults()
            self._setup()


        def setup_for_integration_test(self):
            TestMainPartial.TestMainPartialGrid.set_defaults(self)
            TestMainPartial.TestMainPartialGrid.set_integration_defaults(self)
            self.multi_file_data = True
            self.data_filter = 'demo_wrf_2017-10'
            self.data_location = 'examples/net_cdf/multi_file_range/inputs'

            self._override_defaults()
            # Expected Results
            self.output_logs = [
                {"pody": i + 1} for i in range(self.row_count)
            ]

        @setup_test(TestTypes.UNIT_TEST)
        def test_runs_unit_test(self):
            next(self.default_run())

        @setup_test(TestTypes.INTEGRATION_TEST)
        def test_runs_integration_test(self):
            next(self.default_run())

    # TODO: Implement test
    # class TestSingleFileTimeRange(TestMainPartialGrid):

    #     def _setup(self):
    #         self.row_count = 24 * 2

    #     def setup_for_unit_test(self):
    #         TestMainPartial.TestMainPartialGrid.set_defaults(self)
    #         self._setup()

    #         # Expected Results
    #         self.output_logs = [
    #             {"a": i+1} for i in range(self.row_count)
    #         ]

    #     def setup_for_integration_test(self):
    #         TestMainPartial.TestMainPartialGrid.set_defaults(self)
    #         TestMainPartial.TestMainPartialGrid.set_integration_defaults(self)
    #         self.data_location = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'
    #         self._setup()
    #         # Expected Results
    #         self.output_logs = [
    #             {"pody": i+1} for i in range(self.row_count)
    #         ]

    class TestMultiHourData:
        def test_should_run_without_errors(self):
            # TODO: Implement test
            pass
            # raise NotImplementedError("Implement Test")

# TODO: Move these tests to correct func
# class TestMainGridSeq:

#     @pytest.fixture(autouse=True)
#     def _setup(self, mocker):
#         # Default params
#         self.project_dir = 'project_dir'
#         self.runid = 'test_run'
#         self.output_fields = 'pody'
#         self.runnotes = 'runnotes'

#         # Expected params
#         self.config_names = ['config']
#         self.configs = [f'{c}.json' for c in self.config_names]
#         self.inputs = ['input.json']
#         self.runs = [('final_state', 'output_logs', 'model')]
#         self.variable_map = MagicMock()
#         self.preprocess_map = MagicMock()
#         self.first_run_dir = f"{self.project_dir}/runs/{self.runid}/{self.config_names[0]}"

#         self.processed_configs_dir = f'{self.first_run_dir}/processed_configs'
#         self.base_config_file = f"{self.project_dir}/base_config.json"
#         self.grid_coords = [
#             [0, 0], [0, 1], [0, 2],
#             [1, 0], [1, 1], [1, 2],
#             [2, 0], [2, 1], [2, 2],
#         ]
#         self.input_data_files = [f"{self.project_dir}/inputs/{f}" for f in self.inputs]
#         self.previous_hour_state_path = f"{self.first_run_dir}/prev_state"
#         self.output_data_dir = f"{self.first_run_dir}/outputs_grid"
#         self.live_state_dir = f"{self.first_run_dir}/current_state"
#         self.logger = MagicMock()

#         # Mocks
#         self.spy_os_makedirs = mocker.patch('os.makedirs')
#         mocker.patch('os.rename')
#         mocker.patch('shutil.rmtree')
#         mocker.patch('builtins.open', return_value=MagicMock())

#         mocker.patch(__module_loc__ + '.Logger', return_value=self.logger)
#         mocker.patch(__module_loc__ + '.get_configs', return_value=self.configs)
#         mocker.patch(__module_loc__ + '.get_input_files_list', return_value=self.inputs)

#         mocker.patch(
#             __module_loc__ + '.json_loader',
#             side_effect=lambda p: self.variable_map if 'variable_map' in p
#             else self.preprocess_map if 'preprocess_map' in p
#             else None)
#         mocker.patch(
#             __module_loc__ + '.csv_loader', side_effect=lambda _: [{'x': x, 'y': y} for x, y in self.grid_coords])
#         mocker.patch(
#             __module_loc__ + '.generate_run_notes', return_value=[])
#         self.spy_main_hour = mocker.patch(__module_loc__ + '.main_hour', return_value=self.runs)

#     def default_run(self):
#         return main_grid_seq(
#             self.project_dir,
#             self.runid,
#             self.output_fields,
#             self.runnotes,
#             seperate_live_state=True,
#         )

#     def test_runs_without_errors(self):
#         self.default_run()

#     def test_calls_main_hour(self):
#         self.default_run()
#         self.spy_main_hour.assert_called_once_with(
#             processed_config_dir=f"{self.processed_configs_dir}",
#             external_data_row_path=self.input_data_files[0],
#             previous_hour_state_path=self.previous_hour_state_path,
#             output_data_dir=self.output_data_dir,
#             output_fields=self.output_fields,
#             logger=self.logger,
#             external_file_type=FileTypes.NETCDF,
#             grid_coords=self.grid_coords,
#             netcdf_variable_map=self.variable_map,
#             met_preprocess_map=self.preprocess_map,
#             multi_file_netcdf=False,
#             output_to_netcdf=True,
#             state_out_path=self.live_state_dir,
#         )

#     def test_calls_makedirs(self):
#         self.default_run()
#         assert self.spy_os_makedirs.mock_calls == [
#             call(f'{self.first_run_dir}/current_state', exist_ok=True),
#         ]

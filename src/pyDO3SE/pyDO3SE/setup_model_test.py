"""Unit tests for setup_model.py ."""
# External Libraries
from math import floor
import itertools
import numpy as np
from dataclasses import replace, asdict
import pytest
from unittest.mock import MagicMock, Mock, call
from copy import deepcopy
import xarray as xr
import pandas as pd
from data_helpers.diff import diff

# Internal Libraries
from pyDO3SE.error_handling import DayRangeError
from pyDO3SE.Config.ConfigLocation import Config_Location
from pyDO3SE.External_State.External_State_Config import InputMethod, ThermalTimeMethods
from pyDO3SE.Model_State.model_state_loader import dump_state_to_string
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.Model_State.Model_State import Model_State_Shape

# Libraries we are testing
from pyDO3SE.setup_model import (
    Main_Overrides,
    extract_start_and_end_dates,
    get_grid_coords_from_dataarray,
    initialize_grid_configs,
    load_additional_gridded_config_data,
    pull_config_vars_from_netcdf,
    setup_dd,
    setup_external_state,
    setup_initial_state,
    setup_model,
)

__module_loc__ = __package__ + '.setup_model'

# == Mocked return values == #
mock_config = MagicMock()
mock_initial_state = MagicMock()
mock_model_processes = MagicMock()
mock_external_state = MagicMock()
mock_parameters = MagicMock()
mock_state = MagicMock()
mock_state_out = MagicMock()
mock_run_processes = MagicMock(return_value=mock_state_out)
mock_state_logs = [1, 2, 3]
mock_processes = [1, 2, 3]
mock_overrides = MagicMock()
mock_process_runner = Mock(
    run_processes=mock_run_processes,
    state_logs=mock_state_logs,
)

mock_processes = [1, 2, 3]


class TestSetupModel:

    @pytest.fixture(autouse=True)
    def _setup(self, mocker):
        mocker.patch(__module_loc__ + '.setup_config', return_value=mock_config)
        mocker.patch(__module_loc__ + '.setup_external_state',
                     return_value=(mock_external_state, 0, 100))
        mocker.patch(__module_loc__ + '.setup_initial_state', return_value=mock_initial_state)
        mocker.patch(__module_loc__ + '.setup_model_processes', return_value=mock_model_processes)

        mocker.patch(__module_loc__ + '.load_external_state',
                     return_value=(iter(([None, None], mock_external_state) for _ in range(20))))

    def test_should_call_each_setup_function(self, mocker):

        config_location = "demo_config_location"
        data_location = "demo_data_location"
        (config,
         external_state,
         initial_state,
         model_processes) = setup_model(config_location, data_location, mock_overrides)

        assert initial_state == mock_initial_state
        assert config == mock_config
        assert external_state == mock_external_state
        assert model_processes == mock_model_processes


external_state_data_demo = External_State_Shape(
    dd=np.array([floor(i / 24) for i in range(0, 120 * 24)]),
    hr=[i % 24 for i in range(0, 120 * 24)],
    Ts_C=[abs((i % 24) - 12) for i in range(0, 120 * 24)],
    P=[102 for _ in range(0, 120 * 24)],
    precip=[abs((i % 24) - 12) for i in range(0, 120 * 24)],
    u=[4 for _ in range(0, 120 * 24)],
    O3=[34 for _ in range(0, 120 * 24)],
    PAR=[30 * abs((i % 24) - 12) for i in range(0, 120 * 24)],
    VPD=[0.1 for _ in range(0, 120 * 24)],
)
external_state_data_hour_demo = External_State_Shape(
    dd=np.array([32]),
    hr=[0],
    Ts_C=[12],
    P=[102],
    precip=[0],
    u=[4],
    O3=[34],
    PAR=[90],
    VPD=[0.1],
)
demo_config = Config_Shape()
demo_config.Location.lon = 1
demo_config.Location.lat = 1
demo_config.Location.elev = 1
demo_config.Location.albedo = 1
demo_config_hourly = replace(demo_config)
demo_config_hourly.Met.thermal_time_method = ThermalTimeMethods.HOURLY
demo_config_hourly.Met.inputs.td_method = InputMethod.SKIP


setup_external_state_test_params = [
    # start_day, end_day, external_state_data, config
    (101, 103, external_state_data_demo, demo_config),
    (98, 103, external_state_data_demo, demo_config),
    (101, 103, external_state_data_demo, demo_config),
    (None, None, external_state_data_hour_demo, demo_config_hourly),
]


@pytest.fixture(scope="class", params=setup_external_state_test_params)
def setup_external_state_fixture(request):
    try:
        start_day, end_day, external_state_data, config = request.param
        request.cls.runner = None
        request.cls.start_day = start_day
        request.cls.end_day = end_day
        request.cls.config = config
        external_state_data=deepcopy(external_state_data)
        external_state_data.dd += start_day or 0
        print(start_day)
        config.Location.start_day = start_day
        config.Location.end_day = end_day
        external_state, start_day_out, end_day_out = setup_external_state(
            config=config,
            external_state_data=deepcopy(external_state_data),
            overrides=Main_Overrides(start_day=start_day, end_day=end_day)
        )
        request.cls.external_state_out = external_state
    except Exception as e:
        print("setup_external_state_fixture FAILED")
        raise e


@pytest.mark.usefixtures('setup_external_state_fixture')
class TestSetupExternalState:

    def test_should_have_correct_dd(self):
        assert self.external_state_out.dd is not None
        assert self.external_state_out.dd[0] is not None
        assert self.external_state_out.dd[-1] is not None
        assert self.external_state_out.dd[0] >= 0
        assert self.external_state_out.dd[-1] > 1

    def test_should_have_correct_hr(self):
        assert self.external_state_out.hr is not None
        assert self.external_state_out.hr[0] is not None
        assert self.external_state_out.hr[-1] is not None
        assert self.external_state_out.hr[0] >= 0
        assert self.external_state_out.hr[-1] >= 0
        if len(self.external_state_out.hr) > 1:
            assert self.external_state_out.hr[0] < self.external_state_out.hr[1]


    def test_should_have_validated_inputs(self):
        pass

    def test_setup_constant_values(self):
        pass

    def test_should_setup_O3_input(self):
        # TODO: Set up this test
        pass

    def test_should_setup_CO2_input(self):
        # TODO: Set up this test
        pass

    def test_should_calculate_humidity(self):
        pass

    def test_should_calculate_solar_elevation(self):
        pass

    def test_should_calculate_radiation(self):
        pass

    def test_should_calculate_is_daylight(self):
        pass

    def test_should_calculate_vpd(self):
        pass

    def test_should_calculate_thermal_time_diff(self):
        if self.config.Met.inputs.td_method != InputMethod.SKIP:
            # TODO: At the moment we are calculating thermal time acc
            assert self.external_state_out.td is not None
            assert self.external_state_out.td[0] is not None
            assert self.external_state_out.td[-1] is not None
            assert self.external_state_out.td[0] >= 0
            assert self.external_state_out.td[-1] > 1


# TODO: Merge with above tests
# class TestSetupExternalStateB:

#     @pytest.fixture(autouse=True)
#     def _setup(self, mocker):
#         mocker.patch(__module_loc__ + '.setup_config', return_value=mock_config)
#         mocker.patch(__module_loc__ + '.setup_initial_state', return_value=mock_initial_state)
#         mocker.patch(__module_loc__ + '.setup_model_processes', return_value=mock_model_processes)
#         mocker.patch(__module_loc__ + '.external_state_init_processes', return_value=[])
#         mocker.patch(__module_loc__ + '.asdict',
#                      side_effect=lambda x: {"dd": x.dd} if x == mock_external_state else {})

#     def test_should_slice_data_to_fit_start_and_end_days(self):
#         data_start_day = 99
#         start_day = 101
#         end_day = 103
#         data_end_day = 108

#         mock_external_state.dd = [data_start_day +
#                                   floor(i / 24) for i in range(24 * (data_end_day - data_start_day))]
#         external_state, start_day_out, end_day_out = setup_external_state(
#             config=mock_config,
#             external_state_data=mock_external_state,
#             overrides=Main_Overrides(start_day=start_day, end_day=end_day)
#         )
#         assert len(external_state.dd) == len(external_state.hr)
#         assert len(external_state.dd) == (end_day - start_day) * 24
#         assert external_state.dd[0] == start_day
#         assert external_state.dd[-1] == end_day
#         assert all([len(v) == len(external_state.dd)
#                     if v is not None else True for v in asdict(external_state).values()])

#     def test_should_handle_dd_spanning_multiple_years(self):
#         start_day = 361
#         end_day = 365 + 3

#         mock_external_state.dd = [360 + floor(i / 24) for i in range(24 * 10)]
#         mock_external_state.dd = [dd if dd < 365 else dd - 365 for dd in mock_external_state.dd]
#         external_state, start_day_out, end_day_out = setup_external_state(
#             config=mock_config,
#             external_state_data=mock_external_state,
#             overrides=Main_Overrides(start_day=start_day, end_day=end_day)
#         )
#         assert len(external_state.dd) == len(external_state.hr)
#         assert len(external_state.dd) == (end_day - start_day) * 24
#         assert external_state.dd[0] == start_day
#         assert external_state.dd[-1] == end_day
#         assert all([len(v) == len(external_state.dd)
#                     if v is not None else True for v in asdict(external_state).values()])


class TestSetupInitialState:

    @pytest.fixture(autouse=True)
    def DEMO_CONFIG(self):
        return Config_Shape()

    @pytest.fixture(autouse=True)
    def DEMO_INITIAL_STATE(self):
        # TODO: Modify this
        return Model_State_Shape()

    @pytest.fixture(autouse=True)
    def DEMO_EXTERNAL_STATE(self):
        return External_State_Shape(
            dd=[4, 5, 6],
            Ts_C=[3, 2, 1],
            td=[3, 3, 3],
        )

    @pytest.mark.skip(reason="not implemented")
    def test_should_work_with_empty_initial_state(
        self,
        DEMO_CONFIG,
        DEMO_EXTERNAL_STATE,
        DEMO_INITIAL_STATE,
    ):
        overrides = Main_Overrides()
        initial_state_file = None
        initial_state: Model_State_Shape = setup_initial_state(
            DEMO_CONFIG, DEMO_EXTERNAL_STATE, Model_State_Shape(), initial_state_file, overrides=overrides)
        assert initial_state.temporal.dd == DEMO_CONFIG.Location.start_day
        # assert initial_state == DEMO_INITIAL_STATE
        compared = diff("model_state", dump_state_to_string(
            initial_state), dump_state_to_string(DEMO_INITIAL_STATE))
        if len(compared) > 0:
            compared_message = "\n".join(compared)
            raise AssertionError(
                f"Output and expected output do not match: \n expected -> actual \n\n {compared_message}")

    # def test_should_work_with_none_empty_initial_state(
    #     self,
    #     DEMO_CONFIG,
    #     DEMO_EXTERNAL_STATE,
    #     DEMO_INITIAL_STATE,
    # ):
    #     overrides = Main_Overrides()
    #     initial_state: Model_State_Shape = setup_initial_state(
    #         DEMO_CONFIG, DEMO_EXTERNAL_STATE, DEMO_INITIAL_STATE, overrides)
    #     assert initial_state.temporal.dd == 4


class TestPullConfigVarsFromNetcdf:

    def test_should_return_values_at_coord(self):
        e_state_overrides_field_map = {
            'Location.lat': "latitude",
            'Location.lon': "longitude",
            'hi': 'hello',
        }
        coords = [
            (0, 0),
            (0, 1),
            (1, 1),
            (1, 0),
        ]
        index = pd.MultiIndex.from_tuples(coords, names=["x", "y"])
        ds = xr.Dataset.from_dataframe(
            pd.DataFrame(
                [{
                    # 'x': coord[0],
                    # 'y': coord[1],
                    'longitude': coord[0],
                    'latitude': coord[1],
                    'hello': "world",
                } for coord in coords],
                index=index,
            )
        )

        for coord in coords:
            field_map = pull_config_vars_from_netcdf(
                ds,
                coord,
                e_state_overrides_field_map
            )
            assert field_map == {
                'Location.lat': coord[1],
                'Location.lon': coord[0],
                'hi': 'world',
            }


class TestInitializeGridConfigs:

    @pytest.fixture(autouse=True)
    def _setup(self, mocker):
        # = default params
        self.config_dir = "config_dir"
        self.processed_config_dir = "processed_config_dir"
        self.base_config_path = "base_config_path"
        self.grid_coords = [(0,0), (0,1)]
        self.overrides = Main_Overrides(
        )
        self.mock_config = deepcopy(mock_config)

        self.mock_config.Location = Config_Location()

        self.logger = lambda *args, **kwargs: None
        self.config_name = 'config'
        self.config_path = f'{self.config_name}.json'

        self.e_state_overrides_file_path = 'e_state_overrides_file_path'
        self.e_state_overrides_field_map = 'e_state_overrides_field_map'

        # = Mocks

        self.spy_os_listdir = mocker.patch(
            'os.listdir', return_value=[self.config_path])
        self.spy_dump_config_to_file_binary = mocker.patch(
            __module_loc__ + '.dump_config_to_file_binary')
        self.spy_setup_config = mocker.patch(
            __module_loc__ + '.setup_config', return_value=self.mock_config)
        self.spy_load_additional_gridded_config_data = mocker.patch(
            __module_loc__ + '.load_additional_gridded_config_data',
            return_value=[
                ((xi, yi), {
                    "Location.lat": yi,
                    "Location.lon": xi,
                }) for (xi, yi) in self.grid_coords
            ],
        )

    def run_default(self):
        initialize_grid_configs(
            self.mock_config,
            self.processed_config_dir,
            self.grid_coords,
            self.e_state_overrides_file_path,
            self.e_state_overrides_field_map,
            self.logger,
        )

    # def test_should_call_setup_config(self):
    #     self.run_default()
    #     self.spy_setup_config.assert_called_once_with(
    #         self.config_path,
    #         None,
    #         self.base_config_path,
    #         self.overrides,
    #     )

    def test_should_call_dump_config_file_to_pickle_with_each_coord(self, mocker):
        # Patch deep copy so mock config has same id
        self.spy_deepcopy = mocker.patch(__module_loc__ + '.deepcopy', side_effect=lambda x: x)
        self.run_default()
        assert self.spy_deepcopy.call_count == len(self.grid_coords)
        expected_calls = [
            call.self.spy_dump_config_to_file_pickle(
                self.mock_config,
                f"{self.processed_config_dir}/{xi}_{yi}.config",
            )
            for xi, yi in self.grid_coords
        ]
        assert self.spy_dump_config_to_file_binary.mock_calls == expected_calls

    def test_should_called_load_additional_grid_data(self):
        self.run_default()
        self.spy_load_additional_gridded_config_data.assert_called_once_with(
            self.grid_coords,
            self.e_state_overrides_file_path,
            self.e_state_overrides_field_map,
        )

    def test_should_have_set_coord_specific_params(self):
        self.run_default()
        config_1 = self.spy_dump_config_to_file_binary.mock_calls[0].args[0]
        assert config_1.Location.lon == self.grid_coords[0][0]
        assert config_1.Location.lat == self.grid_coords[0][1]

class TestLoadAdditionalGriddedConfigData:

    @pytest.fixture(autouse=True)
    def _setup(self, mocker):
        self.grid_x = [0, 1, 2]
        self.grid_y = [0, 1, 2]
        self.grid_coords = list(itertools.product(self.grid_x, self.grid_y))

        self.e_state_overrides_file_path = 'e_state_overrides_file_path'
        self.e_state_overrides_field_map = {
            'Location.lat': "latitude",
            'Location.lon': "longitude",
        }
        self.mock_xr = MagicMock
        self.mock_dataset = dict()

        mock_lat_data = np.array([[i for _ in range(len(self.grid_y))]
                                  for i in range(len(self.grid_x))])
        mock_lon_data = np.array([[i for i in range(len(self.grid_y))]
                                  for _ in range(len(self.grid_x))])

        self.mock_dataset['latitude'] = mock_lat_data
        self.mock_dataset['longitude'] = mock_lon_data
        self.mock_xr.open_dataset = MagicMock(return_value=self.mock_dataset)
        mocker.patch(__module_loc__ + '.xr', return_value=self.mock_xr)

    def run_default(self):
        return load_additional_gridded_config_data(
            self.grid_coords,
            self.e_state_overrides_file_path,
            self.e_state_overrides_field_map,
        )

    def test_should_return_an_iterator(self):
        out = self.run_default()
        # = returns an iterator of each model variation
        assert len(list(out)) == len(self.grid_coords)

    def test_should_output_the_coordinate(self):
        out = self.run_default()
        coord, field_map = next(out)
        assert coord == self.grid_coords[0]

    def test_should_output_a_config_object(self):
        out = self.run_default()
        coord, field_map = next(out)
        assert isinstance(field_map, dict)

    def test_should_have_overriden_values_in_config(self):
        out = self.run_default()
        coord, field_map = next(out)
        # NOTE:We set the lat and lon to match grid coord
        assert field_map['Location.lat'] == self.grid_coords[0][0] == coord[0]
        assert field_map['Location.lon'] == self.grid_coords[0][1] == coord[1]

        coord, field_map = next(out)
        assert field_map['Location.lat'] == self.grid_coords[1][0] == coord[0]
        assert field_map['Location.lon'] == self.grid_coords[1][1] == coord[1]

    def test_should_load_external_state_file(self):
        next(self.run_default())
        self.mock_xr.open_dataset.assert_called_once_with(self.e_state_overrides_file_path)


# class TestSetupModelSingleHour:

#     @pytest.fixture(autouse=True)
#     def _setup(self, mocker):
#         # == default parameters
#         self.data_start_day = 1
#         self.data_end_day = 1
#         self.date = '2017-10-01T01:00:00'
#         self.config_path = 'config_path'
#         self.base_config_path = 'base_config_path'
#         self.base_state_path = 'base_state_path'
#         self.extenal_data_path = 'external_data_path'
#         self.e_state_init_processes = MagicMock()
#         self.previous_hour_state_path = 'prev_state_path'
#         self.coord = (5, 7)
#         self.overrides = Main_Overrides(
#             # external_file_type=FileTypes.NETCDF,
#         )
#         self.run_id = "runid"
#         self.mock_processed_external_state = deepcopy(mock_external_state)
#         self.mock_processed_external_state.hr = [0]
#         self.mock_processed_external_state.dd = [self.data_start_day]
#         self.mock_processed_external_state.date = [self.date]

#         # == Mocks
#         self.mock_config_loader_pickled = mocker.patch(
#             __module_loc__ + '.config_loader_pickled',
#             return_value=mock_config,
#         )
#         self.mock_setup_external_state_simple = mocker.patch(
#             __module_loc__ + '.setup_external_state_simple',
#             return_value=self.mock_processed_external_state,
#         )
#         self.mock_hourly_run_model_processes = mocker.patch(
#             __module_loc__ + '.hourly_run_model_processes',
#             return_value=mock_model_processes,
#         )
#         # = Setup mock manager
#         self.mock_manager = Mock()
#         self.mock_manager.attach_mock(self.mock_config_loader_pickled, 'config_loader_pickled')
#         self.mock_manager.attach_mock(self.mock_setup_external_state_simple, 'setup_external_state_simple')
#         self.mock_manager.attach_mock(
#             self.mock_hourly_run_model_processes, 'hourly_run_model_processes')


#     def run_default(self):
#         model = setup_model_single_hour(
#             self.coord,
#             self.config_path,
#             external_state_data_hour_demo,
#             self.previous_hour_state_path,
#             self.e_state_init_processes,
#         )
#         return model

#     def test_should_return_model(self):
#         model = self.run_default()
#         # = returns an iterator of each model variation
#         config, external_state, initial_state, model_processes, x, y = model
#         assert config.id != mock_config.id  # Has been deep copied
#         assert external_state == self.mock_processed_external_state
#         assert initial_state == mock_initial_state
#         assert model_processes == mock_model_processes
#         assert x == self.coord[0]
#         assert y == self.coord[1]

#     def test_setup_external_state_simple_called(self):
#         self.run_default()
#         self.mock_setup_external_state_simple.assert_called_once_with(
#             external_state_data_hour_demo,
#             mock_config,
#             self.e_state_init_processes,
#         )


#     def test_should_call_all_mocked_functions(self):
#         self.run_default()
#         expected_calls = [
#             call.config_loader_pickled(
#                 self.config_path,
#             ),
#             call.setup_external_state_simple(
#                 external_state_data_hour_demo,
#                 mock_config,
#                 self.e_state_init_processes,
#             ),
#             call.hourly_run_model_processes(
#                 mock_config,
#                 self.mock_processed_external_state.hr[0],
#             ),
#         ]
#         assert self.mock_manager.mock_calls == expected_calls

#     def test_should_have_set_config_location(self):
#         model = self.run_default()
#         config, external_state, initial_state, model_processes, x, y= model
#         # assert config.Location.id == self.mock_config_location.id
#         assert config.Location.start_day == self.data_start_day
#         assert config.Location.end_day == self.data_end_day


# class TestSetupModelSingleHourGrid:

#     @pytest.fixture(autouse=True)
#     def _setup(self, mocker):
#         self.grid_x = [0, 1, 2]
#         self.grid_y = [0, 1, 2]
#         self.grid_coords = list(itertools.product(self.grid_x, self.grid_y))
#         self.mock_model = MagicMock()

#         self.config_name = 'config'
#         self.config_path = f'{self.config_name}.config'
#         self.extenal_data_path = 'external_data_path'
#         self.previous_hour_state_path = 'prev_state_path'
#         self.base_config_path = 'base_config_path'
#         self.base_state_path = 'base_state_path'
#         self.e_state_init_processes = MagicMock()
#         self.e_state_overrides_file_path = 'e_state_overrides_file_path'
#         self.e_state_overrides_field_map = {
#             'Location.lat': "latitude",
#             'Location.lon': "longitude",
#         }
#         self.run_id = "runid"

#         self.overrides = Main_Overrides(
#             # external_file_type=FileTypes.NETCDF,
#             # grid_coords=self.grid_coords,
#             # netcdf_variable_map={},
#             # met_preprocess_map={},
#             e_state_overrides_file_path=self.e_state_overrides_file_path,
#             e_state_overrides_field_map=self.e_state_overrides_field_map,
#         )

#         # == Mocks
#         self.mock_config_loader_pickled = mocker.patch(
#             __module_loc__ + '.config_loader_pickled',
#             return_value=mock_config,
#         )
#         self.mock_external_state_init_processes = mocker.patch(
#             __module_loc__ + '.external_state_init_processes',
#             return_value=self.e_state_init_processes,
#         )
#         mocker.patch(__module_loc__ + '.load_external_state',
#                      return_value=([coord, mock_external_state] for coord in self.grid_coords))

#         self.spy_setup_model_single_hour = mocker.patch(
#             __module_loc__ + '.setup_model_single_hour', return_value=self.mock_model)


#     def run_default(self):
#         models = setup_model_single_hour_grid(
#             self.config_path,
#             self.extenal_data_path,
#             self.previous_hour_state_path,
#             overrides=self.overrides,
#         )
#         return models

#     def test_should_init_grid_model(self):
#         models = self.run_default()
#         # = returns an iterator of each model variation
#         assert len(list(models)) == len(self.grid_coords)

#     def test_should_call_setup_model_hour_for_each_coord(self):
#         for x, y in self.grid_coords:
#             next(self.run_default())
#             self.spy_setup_model_single_hour.assert_called_with(
#                 (x, y),
#                 f'{self.config_path}/{x}_{y}.config',
#                 mock_external_state,
#                 previous_hour_state_path=self.previous_hour_state_path,
#                 e_state_init_processes=self.e_state_init_processes,
#             )


# @pytest.mark.parametrize(['start_day', 'end_day', 'row_count'], [
#     [0,1,12],
#     [340,455,2808], # Should be max 2760
#     [0,365,8760],
#     [0,366,8784],
#     [0,367,8808],
# ])
# def test_setup_dd(start_day, end_day, row_count):
#     out = setup_dd(start_day, end_day, row_count)
#     assert len(out) == row_count

# @pytest.mark.parametrize(['dd', 'is_true'], [
#     [range(0,1), False],
#     [range(0,365), False],
#     [range(0,366), False],
#     [range(0,367), True],
#     [range(100,367), True],
#     [range(100,400), True],
#     [range(100,800), True],
#     [range(0,730), True],
#     [range(0,733), True],
#     # [range(365, 400), False],
# ])
# def test_get_years_spanned_count(dd, is_true):
#     out = get_years_spanned_count(list(dd))
#     assert out == is_true


class TestSetupDD:

    def test_some_days(self):
        day_count = 10
        dd_in = [_dd for _dd in range(day_count) for _ in range(24)]
        start_day = 0
        end_day = 9
        dd_out, start_row, end_row = setup_dd(
            dd=dd_in,
            start_day=start_day,
            end_day=end_day,
        )
        assert len(dd_out) > 0
        assert len(dd_out) == end_row - start_row + 1
        assert dd_out[0] == start_day
        assert dd_out[-1] == end_day

    @pytest.mark.parametrize(['dd', 'last_val'], [
        [range(1,2), 1],
        [range(1,365), 364],
        [range(1,366), 365],
        [range(1,367), 366],
        [range(100,367), 366],
        [range(100,400), 399],
        [range(100,800), 799],
        [range(0,730), 729],
        [range(0,733), 732],
        [list(range(0,366)) + list(range(0,1)), 365],
        [list(range(0,365)) + list(range(0,10)), 373],
        [list(range(0,366)) + list(range(0,10)), 374],
        [list(range(365,366)) + list(range(0,1)), 365],
        [list(range(365,366)) + list(range(0,2)), 366],
        [list(range(365,367)) + list(range(0,1)), 366],
        [list(range(365,367)) + list(range(0,2)), 367],
    ])
    @pytest.mark.parametrize(['start_day', 'end_day'], [
        [None, None],
        [1, None],
        [2, None],
        [None, 1],
        [None, 10],
    ])
    def test_setup_dd(self, dd, last_val, start_day, end_day):
        dd_in = [dd for dd in dd for _ in range(24)]
        if start_day is None and end_day is None:
            out, start_row, end_row = setup_dd(dd_in)
            assert len(out) == len(dd_in)
            assert out[-1] == last_val
            assert start_row == 0
            assert end_row == len(out) - 1
        elif start_day is not None and end_day is None:
            if start_day > dd_in[0]:
                # Test that input is clipped by start day
                out, start_row, end_row = setup_dd(dd_in, start_day, end_day)
                assert len(out) == len(dd_in) - (start_day - dd_in[0]) * 24
                assert start_row == (start_day - dd_in[0]) * 24 - 1
                assert end_row == len(dd_in) - 1

            elif start_day < dd_in[0]:
                with pytest.raises(DayRangeError):
                    out, start_row, end_row = setup_dd(dd_in, start_day, end_day)
            elif start_day == dd_in[0]:
                out, start_row, end_row = setup_dd(dd_in, start_day, end_day)
                assert len(out) == len(dd_in)
                assert out[-1] == last_val
                assert end_row == len(out) -1
            else:
                raise NotImplementedError("This condition is not tested!")

        elif start_day is None and end_day is not None:
            if end_day < dd_in[0]:
                with pytest.raises(DayRangeError):
                    out, start_row, end_row = setup_dd(dd_in, start_day, end_day)
                return

            try:
                out, start_row, end_row = setup_dd(dd_in, start_day, end_day)
            except DayRangeError as e:
                # TODO: Assert should only be if end day outside data
                return
            if len(out) > 0:
                assert out[-1] == end_day
                assert start_row == 0
                assert end_row == len(out) - 1

            assert len(out) <= len(dd_in)

        else:
            raise NotImplementedError("This condition is not tested!")


class TestGetGridCoordsFromFile:

    def test_works(self, mocker):
        data = np.zeros((9,8))
        data[0,0] = 1
        data[0,4] = 1
        ds = xr.DataArray(data, coords={"lat": np.arange(9),"lon": np.arange(8)})
        grid_coords, grid_x_size, grid_y_size = get_grid_coords_from_dataarray(ds)
        assert len(grid_coords) == 2
        assert grid_coords[0] == (0,0)
        assert grid_coords[1] == (0,4)

        assert grid_x_size == 9
        assert grid_y_size == 8

class TestExtractStartAndEndDates:

    def test_works(self):
        dd_data = np.array([dd for dd in range(100) for _ in range(24)])
        [start_day, end_day, start_row, end_row, alt_dd] = extract_start_and_end_dates(
            config_start_day=0,
            config_end_day=9,
            external_state_dd_data=dd_data,
            overrides=None,
        )
        assert start_day == 0
        assert end_day == 9
        assert start_row == 0
        assert end_row == 10*24 - 1
        assert len(alt_dd) == end_row - start_row + 1
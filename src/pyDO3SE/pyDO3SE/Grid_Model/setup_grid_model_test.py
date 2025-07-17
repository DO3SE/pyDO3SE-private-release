"""Unit tests for setup_model.py ."""
# External Libraries
from dataclasses import replace
import itertools
import numpy as np
import pytest
from unittest.mock import MagicMock, Mock, call
from copy import deepcopy
import pandas as pd

# Internal Libraries

from pyDO3SE.Config.ConfigLocation import Config_Location
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Config.config_loader import process_json_config

# Libraries we are testing
from pyDO3SE.Grid_Model.setup_grid_model import (
    get_grid_coords_from_dataarray,
    init_grid_model,
    initialize_grid_configs,
    load_additional_gridded_config_data,
    pull_config_vars_from_netcdf,
)
from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.setup_model import (
    Main_Overrides,
)
from pyDO3SE.util.loader import json_loader

__module_loc__ = __package__ + '.setup_grid_model'

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
        self.grid_coords = [(0, 0), (0, 1)]
        self.overrides = Main_Overrides(
        )
        self.mock_config = deepcopy(mock_config)

        self.mock_config.Location = Config_Location()

        self.logger = lambda *args, **kwargs: None
        self.config_name = 'config'
        self.config_path = f'{self.config_name}.json'

        self.e_state_overrides_dataset = MagicMock()
        self.e_state_overrides_field_map = 'e_state_overrides_field_map'

        # = Mocks

        self.spy_os_listdir = mocker.patch(
            'os.listdir', return_value=[self.config_path])
        self.spy_dump_config_to_file_binary = mocker.patch(
            __module_loc__ + '.dump_config_to_file_binary')
        self.spy_setup_config = mocker.patch(
            __module_loc__ + '.setup_config', side_effect=lambda c, *args: c)
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
        return initialize_grid_configs(
            self.mock_config,
            self.grid_coords,
            self.e_state_overrides_dataset,
            self.e_state_overrides_field_map,
            self.logger,
        )

    # TODO: Move to init_all_grid_model_configs
    # def test_should_call_dump_config_file_to_pickle_with_each_coord(self, mocker):
    #     # Patch deep copy so mock config has same id
    #     self.spy_deepcopy = mocker.patch(__module_loc__ + '.deepcopy', side_effect=lambda x: x)
    #     self.run_default()
    #     assert self.spy_deepcopy.call_count == len(self.grid_coords)
    #     expected_calls = [
    #         call.self.spy_dump_config_to_file_pickle(
    #             self.mock_config,
    #             f"{self.processed_config_dir}/{xi}_{yi}.config",
    #         )
    #         for xi, yi in self.grid_coords
    #     ]
    #     assert self.spy_dump_config_to_file_binary.mock_calls == expected_calls

    def test_should_called_load_additional_grid_data(self):
        next(self.run_default())
        self.spy_load_additional_gridded_config_data.assert_called_once_with(
            self.grid_coords,
            self.e_state_overrides_dataset,
            self.e_state_overrides_field_map,
        )

    def test_should_have_set_coord_specific_params(self):
        config_1 = next(self.run_default())
        assert config_1.Location.lon == self.grid_coords[0][0]
        assert config_1.Location.lat == self.grid_coords[0][1]


class TestLoadAdditionalGriddedConfigData:

    @pytest.fixture(autouse=True)
    def _setup(self, mocker):
        self.grid_x = [0, 1, 2]
        self.grid_y = [0, 1, 2]
        self.grid_coords = list(itertools.product(self.grid_x, self.grid_y))

        self.e_state_overrides_field_map = {
            'Location.lat': "latitude",
            'Location.lon': "longitude",
        }
        self.mock_xr = MagicMock
        self.e_state_overrides_dataset = dict()

        mock_lat_data = np.array([[i for _ in range(len(self.grid_y))]
                                  for i in range(len(self.grid_x))])
        mock_lon_data = np.array([[i for i in range(len(self.grid_y))]
                                  for _ in range(len(self.grid_x))])

        self.e_state_overrides_dataset['latitude'] = mock_lat_data
        self.e_state_overrides_dataset['longitude'] = mock_lon_data
        mocker.patch(__module_loc__ + '.xr', return_value=self.mock_xr)

    def run_default(self):
        return load_additional_gridded_config_data(
            self.grid_coords,
            self.e_state_overrides_dataset,
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
        print(field_map['Location.lat'])
        assert field_map['Location.lat'] == self.grid_coords[0][0] == coord[0]
        assert field_map['Location.lon'] == self.grid_coords[0][1] == coord[1]

        coord, field_map = next(out)
        assert field_map['Location.lat'] == self.grid_coords[1][0] == coord[0]
        assert field_map['Location.lon'] == self.grid_coords[1][1] == coord[1]


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

class TestGetGridCoordsFromFile:

    def test_works(self, mocker):
        data = np.zeros((9, 8))
        data[0, 0] = 1
        data[0, 4] = 1
        ds = xr.DataArray(data, coords={"lat": np.arange(9), "lon": np.arange(8)})
        grid_coords, grid_x_size, grid_y_size = get_grid_coords_from_dataarray(ds)
        assert len(grid_coords) == 2
        assert grid_coords[0] == (0, 0)
        assert grid_coords[1] == (0, 4)

        assert grid_x_size == 9
        assert grid_y_size == 8


class TestLatitudeSowing:

    def test_should_get_sowing_date_from_latitude(self):
        base_config = process_json_config(json_loader("examples/net_cdf/full_season/base_config.json"))
        assert base_config.Land_Cover.phenology_options.sowing_day_method == ""
        base_state = Model_State_Shape()
        grid_coords = [[0, 0]]
        configs_out, states_out = init_grid_model(
            base_config,
            base_state,
            e_state_overrides_dataset=None,
            e_state_overrides_field_map={},
            grid_coords=grid_coords,
        )
        config_out = next(configs_out)
        # state_out = next(states_out)
        assert config_out.Land_Cover.parameters[0].phenology.key_dates.sowing == 99

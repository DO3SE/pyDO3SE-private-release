"""Unit tests for setup_model.py ."""

# External Libraries
from math import floor
import numpy as np
from dataclasses import replace
import pytest
from unittest.mock import MagicMock, Mock
from copy import deepcopy
from data_helpers.diff import diff

# Internal Libraries
from pyDO3SE.External_State.External_State_Config import InputMethod, ThermalTimeMethods
from pyDO3SE.Model_State.model_state_loader import dump_state_to_string
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.Model_State.Model_State import Model_State_Shape

# Libraries we are testing
from pyDO3SE.setup_model import (
    Main_Overrides,
    setup_external_state,
    setup_initial_state,
    setup_model,
)

__module_loc__ = __package__ + ".setup_model"

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
        external_state_data = deepcopy(external_state_data)
        external_state_data.dd += start_day or 0
        config.Location.start_day = start_day
        config.Location.end_day = end_day
        external_state, start_day_out, end_day_out = setup_external_state(
            config=config,
            external_state_data=deepcopy(external_state_data),
            overrides=Main_Overrides(start_day=start_day, end_day=end_day),
        )
        request.cls.external_state_out = external_state
    except Exception as e:
        print("setup_external_state_fixture FAILED")
        raise e


@pytest.mark.usefixtures("setup_external_state_fixture")
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
            DEMO_CONFIG,
            DEMO_EXTERNAL_STATE,
            Model_State_Shape(),
            initial_state_file,
            overrides=overrides,
        )
        assert initial_state.temporal.dd == DEMO_CONFIG.Location.start_day
        # assert initial_state == DEMO_INITIAL_STATE
        compared = diff(
            "model_state",
            dump_state_to_string(initial_state),
            dump_state_to_string(DEMO_INITIAL_STATE),
        )
        if len(compared) > 0:
            compared_message = "\n".join(compared)
            raise AssertionError(
                f"Output and expected output do not match: \n expected -> actual \n\n {compared_message}"
            )

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

import pytest
import json
from data_helpers.encoders import AdvancedJsonEncoder
from data_helpers.time_data import get_row_index
from pyDO3SE.External_State.External_State_Config import Config_Met, Config_Met_Inputs, InputMethod
from pyDO3SE.External_State.external_state_loader import FileTypes, load_external_state
from pyDO3SE.Config.Config_Shape import Config_Location, Config_Shape
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from proflow.Objects.Process import Process
from proflow.ProcessRunnerCls import ProcessRunner
from pyDO3SE.util.test_utils import process_snapshot
from .es_init_processes import (
    calc_thermal_time_processes,
    external_state_init_processes,
    init_params,
)


def process_snapshot(data):
    return json.dumps(data, cls=AdvancedJsonEncoder, indent=4, sort_keys=True)


def test_es_init_processes(snapshot):
    config = Config_Shape()
    config.Met = Config_Met()

    config.Location.elev = 20
    processes = external_state_init_processes(start_day=0, end_day=10 * 24, config_met=config.Met)
    assert isinstance(processes[0], Process)
    # assert len(processes) == 81
    snapshot.assert_match(process_snapshot(processes), "es init processes")
    # assert len(flatten_list(processes)) == 81


def test_es_init_params(snapshot):
    config = Config_Shape()
    config.Location.elev = 20
    config.Met.inputs.row_index_method = InputMethod.INPUT
    process_runner = ProcessRunner(config_in=config)
    row_count = 10
    init_processes = init_params()
    initial_state = External_State_Shape(
        row_index=[i for i in range(row_count)],
        dd=[i for i in range(row_count)],
        hr=[i for i in range(row_count)],
        Ts_C=[i for i in range(row_count)],
        O3=[i for i in range(row_count)],
        P=[i for i in range(row_count)],
        precip=[i for i in range(row_count)],
        u=[i for i in range(row_count)],
        PAR=[i for i in range(row_count)],
        VPD=[i for i in range(row_count)],
        Hd=[i for i in range(row_count)],
    )
    final_state = process_runner.run_processes(init_processes, initial_state)
    # snapshot.assert_match(asdict(final_state))
    assert len(final_state.sinB) == row_count
    assert len(final_state.RH) == row_count


def test_throws_error_if_missing_input():
    process_runner = ProcessRunner()
    config_in = Config_Shape()
    config_in.Location.elev = 20
    config_in.Met.h_method = "input"
    process_runner.config = config_in
    init_processes = init_params()
    with pytest.raises(ValueError) as e:
        process_runner.run_processes(
            init_processes, External_State_Shape(row_index=[i for i in range(10)])
        )
    assert "Must supply" in str(e.value)


def test_es_init_processes_run(snapshot):
    EXT_DATA_COLS = [
        # TODO: This should be based on the config
        "PAR",
        "VPD",
        "Ts_C",
        "u",
        "P",
        "O3",
        "dd",
        "hr",
        "precip",
    ]
    start_day = 0
    end_day = 10

    row_indexes = [get_row_index(dd, hr) for dd in range(start_day, end_day) for hr in range(24)]
    # TODO: Replace below loader with static initial state
    data_location = "examples/spanish_wheat/inputs/spanish_wheat_data.csv"
    external_state_data = next(
        load_external_state(
            data_location,
            file_type=FileTypes.CSV,
            row_indexes=row_indexes,
        )
    )
    assert len(external_state_data.O3) == len(row_indexes)
    assert len(external_state_data.O3) == len(row_indexes)

    config = Config_Shape(Location=Config_Location(lat=50, lon=1.3, elev=20, albedo=1))
    config.Met = Config_Met(
        inputs=Config_Met_Inputs(
            row_index_method=InputMethod.CALCULATED, Rn_method=InputMethod.CALCULATED
        )
    )

    process_runner = ProcessRunner()
    process_runner.config = config
    process_runner.external_state = external_state_data
    external_state = process_runner.run_processes(
        external_state_init_processes(
            start_day=start_day + 1,
            end_day=end_day,
            config_met=config.Met,
        ),
        external_state_data,
    )
    assert len(external_state.eact) == (end_day - start_day) * 24
    assert len(external_state.sinB) == (end_day - start_day) * 24

    snapshot.assert_match(process_snapshot(external_state), "external_state ran")


def test_calc_thermal_time_processes():
    config = Config_Shape()
    config.Met.inputs.td_method = InputMethod.CALCULATED
    process_runner = ProcessRunner(config_in=config)
    processes = calc_thermal_time_processes(config.Met)
    row_count = 24 * 5
    initial_state = External_State_Shape(
        dd=[i for i in range(row_count)],
        hr=[i for i in range(row_count)],
        Ts_C=[i for i in range(row_count)],
        O3=[i for i in range(row_count)],
        P=[i for i in range(row_count)],
        precip=[i for i in range(row_count)],
        u=[i for i in range(row_count)],
        PAR=[i for i in range(row_count)],
        VPD=[i for i in range(row_count)],
        Hd=[i for i in range(row_count)],
    )

    final_state = process_runner.run_processes(processes, initial_state)
    assert final_state.td is not None
    assert len(final_state.td) == row_count
    assert final_state.td[0] == 0.0
    assert final_state.td[-1] == 286.0
    assert final_state.td[2 * 24] == 95.0

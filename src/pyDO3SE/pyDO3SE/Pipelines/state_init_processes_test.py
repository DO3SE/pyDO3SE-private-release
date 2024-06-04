from data_helpers.list_helpers import flatten_list

from data_helpers.cls_parsing import _replace_recursive

from proflow.ProcessRunnerCls import ProcessRunner

from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.Config.demo_data.demo_config import DEMO_CONFIG
from pyDO3SE.External_State.External_State_Shape import External_State_Shape

from .state_init_processes import state_init_processes

# import proflow as _Process_Runner

series_array = _d = [i for i in range(365 * 24)]
demo_external_shape = External_State_Shape(
    _d, _d, _d, _d, _d, _d, _d, td=_d,
)

demo_config_b = _replace_recursive(
    DEMO_CONFIG, 'Land_Cover.parameters.0.phenology.key_dates.sowing', 2)


def test_state_init_processes(snapshot):
    processes = state_init_processes(DEMO_CONFIG)
    flattened_process_comments = [p.comment or p.func.__name__ for p in flatten_list(processes)]
    snapshot.assert_match(flattened_process_comments, 'Process_comments')
    # assert len(flattened_process_comments) == 7


def test_state_init_processes_run(snapshot):
    state: Model_State_Shape = Model_State_Shape()
    process_runner = ProcessRunner(config_in=demo_config_b, external_state_in=demo_external_shape)
    run_processes = process_runner.initialize_processes(state_init_processes(config=DEMO_CONFIG))
    state_2: Model_State_Shape = run_processes(initial_state=state)
    snapshot.assert_match(state_2)

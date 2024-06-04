from pyDO3SE.Config.demo_data.demo_config import DEMO_CONFIG
from data_helpers.list_helpers import flatten_list
from proflow.Objects.Process import Process
from pyDO3SE.Config import Config_Shape
from pyDO3SE.Config.Config_Shape import Config_Land_Cover, Config_Land_Cover_Parameters
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.Output.Output_Shape import output_fields

from .default_processes import (
    full_model_processes,
    hourly_processes,
    log_processes,
    set_hour,
)


def test_hourly_processes(snapshot):
    number_of_layers = 3
    number_of_components = 2
    config = Config_Shape(
        Land_Cover=Config_Land_Cover(
            nL=number_of_layers,
            nLC=number_of_components,
            parameters=[
                Config_Land_Cover_Parameters(

                ),
                Config_Land_Cover_Parameters(

                ),
                Config_Land_Cover_Parameters(

                )
            ]
        ),
    )
    config.carbon_allocation.use_carbon_allocation = True
    processes = hourly_processes(config, 0)
    assert isinstance(processes[0], Process)
    assert processes[1].func == set_hour(0).func
    snapshot.assert_match(processes, 'proccesses')
    flattened_process_comments = [p.comment or p.func.__name__ for p in flatten_list(processes)]
    snapshot.assert_match(flattened_process_comments, 'Process_comments')
    # assert len(processes) == 31
    # assert len(flatten_list(processes)) == 140


def test_all_processes_first_24_hours(snapshot):
    hours = list(range(24))
    processes = full_model_processes(DEMO_CONFIG, hours)
    flattened_process_comments = [p.comment or p.func.__name__ for p in flatten_list(processes)]
    snapshot.assert_match(flattened_process_comments, 'Process_comments')
    # assert len(flattened_process_comments) == 876


class TestLogProcesses:

    def test_should_log_all_outputs(self):
        # output_fields
        fields = [f.id for f in output_fields]
        external, state, _ = log_processes(
            nL=1,
            nLC=1,
            nP=1,
            fields=fields,
            log_multilayer=False,
        )
        example_state = Model_State_Shape()
        example_e_state = External_State_Shape(
            *[[1] for _ in External_State_Shape.__annotations__]
        )
        state_inputs = list(filter(lambda f: f.as_ != '|', state.state_inputs(example_state)))
        e_state_inputs = list(external.external_state_inputs(example_e_state, 0))

        log_fields = set([f.as_ for f in state_inputs] + [f.as_ for f in e_state_inputs])
        assert log_fields == set(fields)
        assert len(state_inputs) + len(e_state_inputs) == len(fields)

    def test_should_handle_multi_layer_output(self):
        # output_fields
        fields = [f.id for f in output_fields]
        external, state, _ = log_processes(
            nL=2,
            nLC=1,
            nP=1,
            fields=fields,
            log_multilayer=False,
        )
        example_state = Model_State_Shape()
        example_e_state = External_State_Shape(
            *[[1] for _ in External_State_Shape.__annotations__]
        )
        state_inputs = list(filter(lambda f: f.as_ != '|', state.state_inputs(example_state)))
        e_state_inputs = list(external.external_state_inputs(example_e_state, 0))

        log_fields = set([f.as_ for f in state_inputs] + [f.as_ for f in e_state_inputs])
        assert log_fields == set(fields)
        assert len(state_inputs) + len(e_state_inputs) == len(fields)

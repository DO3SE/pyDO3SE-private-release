import pytest
from data_helpers.list_helpers import flatten_list
from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.Config.demo_data.demo_config import DEMO_CONFIG
from pyDO3SE.External_State.External_State_Shape import External_State_Shape

from proflow.ProcessRunnerCls import ProcessRunner
from .validation_processes import *

series_array = _d = [i for i in range(365 * 24)]
demo_external_shape = External_State_Shape(
    _d, _d, _d, _d, _d, _d, _d
)


class TestValidateExternalState:

    def test_should_throw_error_if_required_values_are_not_truthy(self):
        state: Model_State_Shape = Model_State_Shape()
        e_state: External_State_Shape = External_State_Shape()
        process_runner = ProcessRunner(config_in=DEMO_CONFIG, external_state_in=e_state)
        run_processes = process_runner.initialize_processes(
            flatten_list(validate_external_state(config=DEMO_CONFIG)))
        with pytest.raises(ValueError) as e:
            run_processes(initial_state=state)
        assert "sinB is not defined" in repr(e)

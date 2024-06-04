"""The model run function.

This is where the actual model starts!
The inputs should come from setup_model.
"""

from typing import Any, List, Tuple
from copy import deepcopy
from proflow.ProcessRunnerCls import ProcessRunner
from proflow.Objects.Process import Process

from pyDO3SE.Config import Config_Shape
from pyDO3SE.External_State import External_State_Shape
from pyDO3SE.Model_State import Model_State_Shape


def run_model_daily(
    initial_state: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    model_processes: List[Process],
    # TODO: Implement log level
    # log_level: int = 0,
    DEBUG_MODE: bool = False,
) -> Tuple[Model_State_Shape, List[List[Any]]]:
    """Run the model between day range from the point we have received all external data."""
    process_runner = ProcessRunner(config, external_state, DEBUG_MODE=DEBUG_MODE)

    model_state = deepcopy(initial_state)
    assert config.Location.start_day is not None, "Missing config.Location.start_day. Possibly not setup correctly."
    assert config.Location.end_day is not None, "Missing config.Location.end_day. Possibly not setup correctly."
    assert (config.Location.end_day - config.Location.start_day + 1) * 24 == len(external_state.dd), "Day range and external data length do not match!"
    for _ in range(config.Location.start_day, config.Location.end_day + 1):
        model_state = process_runner.run_processes(
            model_processes,
            model_state,
        )

    output_logs = process_runner.state_logs
    return model_state, output_logs


def run_model(
    initial_state: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    model_processes: List[Process],
    # TODO: Implement log level
    # log_level: int = 0,
    DEBUG_MODE: bool = False,
) -> Tuple[Model_State_Shape, List[List[Any]]]:
    """Run the model from the point we have received all external data."""
    process_runner = ProcessRunner(config, external_state, DEBUG_MODE=DEBUG_MODE)

    model_state = deepcopy(initial_state)
    final_state = process_runner.run_processes(
        model_processes,
        model_state,
    )

    output_logs = process_runner.state_logs
    return final_state, output_logs

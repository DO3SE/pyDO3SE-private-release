"""The model run function.

This is where the actual model starts!
The inputs should come from setup_model.
"""

from typing import Any, Callable, List, Tuple, Dict
from copy import deepcopy
from proflow.ProcessRunnerCls import ProcessRunner
from proflow.Objects.Process import Process

from pyDO3SE.Config import Config_Shape
from pyDO3SE.External_State import External_State_Shape
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.util.logger import Logger


def run_model_daily(
    initial_state: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    model_processes: List[Process],
    start_index: int = None,
    logger: Callable[[str], None] = Logger(),
    DEBUG_MODE: bool = False,
) -> Tuple[Model_State_Shape, List[List[Any]]]:
    """Run the model between day range from the point we have received all external data."""
    process_runner = ProcessRunner(config, external_state, DEBUG_MODE=DEBUG_MODE)
    model_state = deepcopy(initial_state)

    start_index = start_index if start_index is not None else model_state.temporal.row_index or 0
    process_runner.tm.row_index = model_state.temporal.row_index or 0

    start_day = config.Location.start_day
    end_day = config.Location.end_day

    logger(
        f"Running daily model from start day:{start_day} to end day{end_day}")

    # Validate Input
    assert config.output.fields and len(
        config.output.fields) > 0, "Must supply output fields in config or cli args!"
    assert start_day is not None, "Missing config.Location.start_day. Possibly not setup correctly."
    assert end_day is not None, "Missing config.Location.end_day. Possibly not setup correctly."
    assert (end_day - start_day + 1) * \
        24 == len(
            external_state.dd), f"Day range and external data length do not match!\nstart_day:{start_day}, end_day: {end_day}, ext length= {len(external_state.dd)}"

    # Run for day range
    for _ in range(start_day, end_day + 1):
        model_state = process_runner.run_processes(
            model_processes,
            model_state,
        )

    output_logs = process_runner.state_logs
    return model_state, output_logs


def run_model_on_mapped_processes(
    initial_state: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    model_processes: Dict[int, List[Process]],
    start_index: int = None,
    logger: Callable[[str], None] = Logger(),
    DEBUG_MODE: bool = False,
) -> Tuple[Model_State_Shape, List[List[Any]]]:
    """Run the model between day range from the point we have received all external data."""
    process_runner = ProcessRunner(config, external_state, DEBUG_MODE=DEBUG_MODE)
    model_state = deepcopy(initial_state)

    start_index = start_index if start_index is not None else model_state.temporal.row_index or 0
    process_runner.tm.row_index = model_state.temporal.row_index or 0

    # Run for day range
    for hr in external_state.hr:
        model_state = process_runner.run_processes(
            model_processes[hr],
            model_state,
        )

    output_logs = process_runner.state_logs
    return model_state, output_logs


def run_model(
    initial_state: Model_State_Shape,
    config: Config_Shape,
    external_state: External_State_Shape,
    model_processes: List[Process],
    start_index: int = None,
    logger: Callable[[str], None] = Logger(),
    DEBUG_MODE: bool = False,
) -> Tuple[Model_State_Shape, List[List[Any]]]:
    """Run the model from the point we have received all external data."""
    logger("Running model")

    # Validate Input
    assert len(config.output.fields) > 0, "Must supply output fields in config or cli args!"

    process_runner = ProcessRunner(config, external_state, DEBUG_MODE=DEBUG_MODE)
    model_state = deepcopy(initial_state)
    start_index = start_index if start_index is not None else model_state.temporal.row_index or 0
    process_runner.tm.row_index = model_state.temporal.row_index or 0

    final_state = process_runner.run_processes(
        model_processes,
        model_state,
    )

    output_logs = process_runner.state_logs
    return final_state, output_logs

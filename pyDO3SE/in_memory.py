"""A version of the main file that has no file I/O."""

from typing import Callable, List, NamedTuple, Dict
from datetime import datetime

from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE.util.logger import Logger
from pyDO3SE.run_model import run_model_on_mapped_processes
from pyDO3SE.Config import Config_Shape
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.setup_model import (
    Main_Overrides,
    setup_model,
)


class MainOutput(NamedTuple):
    final_state: Model_State_Shape
    output_logs: List[List[any]]
    config_processed: Config_Shape
    initial_state: Model_State_Shape
    external_state: External_State_Shape


def main(
    config: Config_Shape,
    state: Model_State_Shape,
    external_state: External_State_Shape,
    config_overrides: Dict[str, any] = {},
    logger: Callable[[str], None] = Logger(),
    *args,
    **kwargs,
) -> MainOutput:
    """Run the full model from a config file location and external data location.

    Takes the config, data and output directory locations.
    Initializes the config then runs the model
    Outputs are saved in the output directory

    Parameters
    ----------
    project_paths: GridProjectPaths
        file paths specific to project
    run_paths : GridRunPaths
        file paths specific to run
    logger: Callable[[str], None], optional
        Log function, by default print
    output_options: OutputOptions
        Additional options around outputs

    *args, **kwargs
        Passed to Main_Overrides


    Returns
    -------
    Tuple[Model_State_Shape, List[List[Any]], Config_Shape, Model_State_Shape]
       final_state, output_logs, config, initial_state

    """
    overrides = Main_Overrides(*args, **kwargs)

    # MODEL SETUP
    [
        config_processed,
        external_state_processed,
        initial_state,
        model_processes,
    ] = setup_model(
        config_in=config,
        state_in=state,
        external_state_in=external_state,
        # TODO: Allow inputing already loaded data here
        # data_location=run_paths.input_data_file_path,
        # run_dir=run_paths.run_dir,
        config_overrides=config_overrides,
        logger=logger,
        overrides=overrides,
    )

    # MODEL RUN
    start_time = datetime.now()
    model_runner = run_model_on_mapped_processes
    try:
        final_state, output_logs = model_runner(
            initial_state,
            config_processed,
            external_state_processed,
            model_processes,
            DEBUG_MODE=overrides.debug,
        )
    except Exception as e:
        logger(f"Error running model")
        raise e
    time_taken = datetime.now() - start_time
    logger(f"Model run complete in {time_taken}")
    return MainOutput(final_state, output_logs, config, initial_state, external_state)

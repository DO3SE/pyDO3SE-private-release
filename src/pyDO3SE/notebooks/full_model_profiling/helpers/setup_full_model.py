"""Test the full model."""
from helpers.named_tuple_helpers import unpack
from proflow.ProcessRunnerCls import ProcessRunner


from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.External_State.external_state_loader import load_external_state
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Defaults.state_init_processes import state_init_processes
from pyDO3SE.Defaults.es_init_processes import external_state_init_processes
from pyDO3SE.Defaults.config_init_processes import config_init_processes


def setup_initial_model(
    start_day=0,
    end_day=40,
    config_location='examples/spanish_wheat/spanish_wheat_config_for_short_test.json',
    data_location='examples/spanish_wheat/spanish_wheat_data.csv',
):

    EXT_DATA_COLS = [
        # TODO: This should be based on the config
        'PAR',
        'VPD',
        'Ts_C',
        'u',
        'P',
        'O3',
        'dd',
        'hr',
        'precip',
    ]

    # %%
    """SETUP CONFIG"""

    config = config_loader(config_location, 'json')
    process_runner = ProcessRunner(config, DEBUG_MODE=True)

    config_amended = process_runner.run_processes(
        config_init_processes(config),
        config)
    process_runner.config = config_amended
    unpack(config_amended)

    # %%
    """Setup External state"""

    external_state_data = load_external_state(data_location, 'csv', EXT_DATA_COLS)
    process_runner.external_state = external_state_data
    external_state = process_runner.run_processes(
        external_state_init_processes(start_day, end_day, config),
        external_state_data)
    process_runner.external_state = external_state

    assert external_state.sinB[0] is not None
    assert process_runner.external_state.sinB[0] is not None

    """Setup initial state"""

    empty_state = Model_State_Shape()
    state_init = process_runner.run_processes(
        state_init_processes(config, start_day, end_day),
        empty_state,
    )
    return state_init, config, process_runner


if __name__ == "__main__":
    state, conf = setup_initial_model()
    print(state.temporal)

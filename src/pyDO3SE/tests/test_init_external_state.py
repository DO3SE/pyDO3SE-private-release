"""Test the full model."""

from dataclasses import asdict
import json

from data_helpers.snapshot_helpers import prep_data_for_snapshot
from data_helpers.encoders import AdvancedJsonEncoder
from data_helpers.cls_parsing import unpack
from proflow.ProcessRunnerCls import ProcessRunner


from pyDO3SE.External_State.external_state_loader import FileTypes, load_external_state
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes
from pyDO3SE.Pipelines.config_init_processes import config_init_processes

DEMO_START_DAY = 0
DEMO_END_DAY = 40
HOURS_IN_DAY = 24
DAYS_IN_YEAR = 365
DAYS_TO_RUN = DEMO_END_DAY - DEMO_START_DAY
HOURS_IN_YEAR = HOURS_IN_DAY * DAYS_IN_YEAR
HOURS_TO_RUN = HOURS_IN_DAY * DAYS_TO_RUN
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


# def test_init_exteranl_state_optimized(snapshot):
#     start_day = DEMO_START_DAY
#     end_day = DEMO_END_DAY
#     config_location = 'examples/spanish_wheat/configs/spanish_wheat_config.json'
#     data_location = 'examples/spanish_wheat/inputs/spanish_wheat_data.csv'
#     output_directory = "demo_output"

#     # == 1. SETUP CONFIG

#     config = config_loader(config_location=config_location, config_type='json')
#     snapshot.assert_match(unpack(config), 'Config')
#     process_runner = ProcessRunner(config, DEBUG_MODE=True)

#     config_amended = process_runner.run_processes(
#         config_init_processes(config),
#         config)
#     process_runner.config = config_amended
#     snapshot.assert_match(unpack(config_amended), 'Config-amended')

#     external_state_data = load_external_state(data_location, FileTypes.CSV, EXT_DATA_COLS)


def test_init_external_state(snapshot):
    start_day = DEMO_START_DAY
    end_day = DEMO_END_DAY
    config_location = 'examples/spanish_wheat/configs/spanish_wheat_config.json'
    data_location = 'examples/spanish_wheat/inputs/spanish_wheat_data.csv'

    # == 1. SETUP CONFIG

    config = config_loader(config_location=config_location, config_type='json')
    process_runner = ProcessRunner(config, DEBUG_MODE=True)

    config_amended = process_runner.run_processes(
        config_init_processes(config),
        config)
    process_runner.config = config_amended

    # == 2. SETUP EXTERNAL STATE

    external_state_data = next(load_external_state(data_location, file_type=FileTypes.CSV))
    process_runner.external_state = external_state_data
    external_state = process_runner.run_processes(
        external_state_init_processes(start_day * 24, end_day * 24, config.Met),
        external_state_data)

    assert external_state.sinB and sum(external_state.sinB) > 0
    assert external_state.VPD and sum(external_state.VPD) > 0
    assert external_state.P and sum(external_state.P) > 0
    assert external_state.PAR and sum(external_state.PAR) > 0
    assert external_state.Idrctt and sum(external_state.Idrctt) > 0
    assert external_state.Idfuse and sum(external_state.Idfuse) > 0
    assert external_state.R and sum(external_state.R) > 0
    assert external_state.PPFD and sum(external_state.PPFD) > 0
    assert external_state.is_daylight and external_state.is_daylight[0] is not None
    assert external_state.VPD_dd and sum(external_state.VPD_dd) > 0
    assert external_state.CO2 and sum(external_state.CO2) > 0
    assert external_state.RH and sum(external_state.RH) > 0
    assert external_state.Ts_C and sum(external_state.Ts_C) > 0
    assert external_state.eact and sum(external_state.eact) > 0
    # assert external_state.uh and sum(external_state.uh) > 0 # Not currently set. In state instead?
    # assert external_state.Rn and sum(external_state.Rn) > 0 # currently always 0.0
    snapshot.assert_match(json.dumps(prep_data_for_snapshot(asdict(process_runner.external_state)),
                                     cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'external_state')
    process_runner.external_state = external_state
    assert process_runner.external_state.sinB[0] is not None

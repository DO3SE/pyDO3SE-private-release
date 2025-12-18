# """Test the full model."""

# from dataclasses import asdict
# import pytest
# import os
# import json
# import pandas as pd

# from data_helpers.encoders import AdvancedJsonEncoder
# from data_helpers.cls_parsing import unpack
# from proflow.ProcessRunnerCls import ProcessRunner


# from pyDO3SE.Model_State.Model_State import Model_State_Shape
# from pyDO3SE.External_State.external_state_loader import FileTypes, load_external_state
# from pyDO3SE.Config.config_loader import config_loader
# from pyDO3SE.Pipelines.state_init_processes import state_init_processes
# from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes
# from pyDO3SE.Pipelines.default_processes import full_model_processes
# from pyDO3SE.Pipelines.config_init_processes import config_init_processes

# SAVE_OUTPUTS_LOCATION = os.path.dirname(os.path.realpath(__file__)) + '/output/'

# DEMO_START_DAY = 1
# DEMO_END_DAY = 365
# HOURS_IN_DAY = 24
# DAYS_IN_YEAR = 365
# DAYS_TO_RUN = DEMO_END_DAY - DEMO_START_DAY + 1
# HOURS_IN_YEAR = HOURS_IN_DAY * DAYS_IN_YEAR
# HOURS_TO_RUN = HOURS_IN_DAY * DAYS_TO_RUN
# EXT_DATA_COLS = [
#     # TODO: This should be based on the config
#     'PAR',
#     'VPD',
#     'Ts_C',
#     'u',
#     'P',
#     'O3',
#     'dd',
#     'hr',
#     'precip',
# ]


# def test_processes(snapshot, benchmark):
#     if os.environ.get('TQUICK', 'False') == 'True':
#         pytest.skip('Skip test in TQUICK mode as it takes a while to run')

#     start_day = DEMO_START_DAY
#     end_day = DEMO_END_DAY
#     config_location = 'examples/spanish_wheat/configs/spanish_wheat_config.json'
#     data_location = 'examples/spanish_wheat/inputs/spanish_wheat_data.csv'

#     # == 1. SETUP CONFIG

#     config = config_loader(config_path=config_location, config_type='json')
#     snapshot.assert_match(unpack(config), 'Config')
#     process_runner = ProcessRunner(config, DEBUG_MODE=True)

#     config_amended = process_runner.run_processes(
#         config_init_processes(config),
#         config)
#     process_runner.config = config_amended
#     snapshot.assert_match(unpack(config_amended), 'Config-amended')

#     # == 2. SETUP EXTERNAL STATE
#     _, external_state_data = next((data_location, FileTypes.CSV))
#     process_runner.external_state = external_state_data
#     external_state = process_runner.run_processes(
#         external_state_init_processes((start_day - 1) * 24, (end_day) * 24, config.Met),
#         external_state_data)

#     assert external_state.sinB and sum(external_state.sinB) > 0
#     assert external_state.VPD and sum(external_state.VPD) > 0
#     assert external_state.P and sum(external_state.P) > 0
#     assert external_state.PAR and sum(external_state.PAR) > 0
#     assert external_state.Idrctt and sum(external_state.Idrctt) > 0
#     assert external_state.Idfuse and sum(external_state.Idfuse) > 0
#     assert external_state.R and sum(external_state.R) > 0
#     assert external_state.PPFD and sum(external_state.PPFD) > 0
#     assert external_state.is_daylight and external_state.is_daylight[0] is not None
#     assert external_state.VPD_dd and sum(external_state.VPD_dd) > 0
#     assert external_state.CO2 and sum(external_state.CO2) > 0
#     assert external_state.RH and sum(external_state.RH) > 0
#     assert external_state.Ts_C and sum(external_state.Ts_C) > 0
#     assert external_state.eact and sum(external_state.eact) > 0
#     assert len(external_state.PAR) == 8760
#     assert len(external_state.Idrctt) == 8760
#     # assert external_state.uh and sum(external_state.uh) > 0 # Not currently set. In state instead?
#     # assert external_state.Rn and sum(external_state.Rn) > 0 # currently always 0.0
#     # snapshot.assert_match(json.dumps(prep_data_for_snapshot(asdict(process_runner.external_state)),
#     #                                  cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'external_state')
#     process_runner.external_state = external_state
#     assert process_runner.external_state.sinB[0] is not None

#     # == SETUP INITIAL STATE
#     initial_state = Model_State_Shape()
#     state = process_runner.run_processes(
#         state_init_processes(config),
#         initial_state,
#     )
#     assert state.canopy_component[0].season_Astart_td is not None

#     # == RUN PROCESSES
#     final_state = process_runner.run_processes(
#         full_model_processes(config),
#         state,
#     )
#     # assert final_state.canopy.LAI_total == 3.5

#     snapshot.assert_match(
#         json.dumps(
#             asdict(final_state.temporal),
#             cls=AdvancedJsonEncoder,
#             indent=4,
#             sort_keys=True,
#         ), 'final_state - temporal',
#     )

#     snapshot.assert_match(
#         json.dumps(
#             asdict(final_state.canopy),

#             cls=AdvancedJsonEncoder,
#             indent=4,
#             sort_keys=True,
#         ), 'final_state - canopy',
#     )
#     snapshot.assert_match(
#         json.dumps(
#             asdict(final_state.canopy_component[0]),
#             cls=AdvancedJsonEncoder,
#             indent=4,
#             sort_keys=True,
#         ), 'final_state - canopy_component',
#     )
#     snapshot.assert_match(
#         json.dumps(
#             asdict(final_state.canopy_layer_component[0][0]),
#             cls=AdvancedJsonEncoder,
#             indent=4,
#             sort_keys=True,
#         ), 'final_state - canopy_layer_component',
#     )
#     snapshot.assert_match(
#         json.dumps(
#             asdict(final_state.canopy_layers[0]),
#             cls=AdvancedJsonEncoder,
#             indent=4,
#             sort_keys=True,
#         ), 'final_state - canopy_layers',
#     )

#     # snapshot.assert_match(process_runner.state_logs, 'logs')

#     # print output logs to csv
#     df_out = pd.DataFrame(process_runner.state_logs)
#     df_out.to_csv(SAVE_OUTPUTS_LOCATION + 'logs_full_year.csv')

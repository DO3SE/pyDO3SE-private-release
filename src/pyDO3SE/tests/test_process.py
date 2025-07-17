# """Test the full model."""

# from dataclasses import asdict
# import json
# import pytest
# import os
# import pandas as pd

# from data_helpers.snapshot_helpers import prep_data_for_snapshot
# from data_helpers.encoders import AdvancedJsonEncoder
# from data_helpers.cls_parsing import unpack
# from proflow.ProcessRunnerCls import ProcessRunner


# from pyDO3SE.Model_State.Model_State import Model_State_Shape
# from pyDO3SE.Config.config_loader import config_loader
# from pyDO3SE.Pipelines.state_init_processes import state_init_processes
# from pyDO3SE.Pipelines.default_processes import full_model_processes
# from pyDO3SE.Pipelines.config_init_processes import config_init_processes
# from pyDO3SE.setup_model import Main_Overrides, setup_external_state

# SAVE_OUTPUTS_LOCATION = os.path.dirname(os.path.realpath(__file__)) + '/output/'

# DEMO_START_DAY = 1
# DEMO_END_DAY = 40

# @pytest.mark.skip('Out of date.')
# def test_processes_photo(snapshot):
#     start_day = DEMO_START_DAY
#     end_day = DEMO_END_DAY
#     config_location = 'examples/spanish_wheat/old_configs/spanish_wheat_config_for_short_test.json'
#     data_location = 'examples/spanish_wheat/inputs/spanish_wheat_data.csv'

#     # == 1. SETUP CONFIG

#     config = config_loader(config_location=config_location, config_type='json')
#     snapshot.assert_match(unpack(config), 'Config')
#     process_runner = ProcessRunner(config, DEBUG_MODE=True)

#     config_amended = process_runner.run_processes(
#         config_init_processes(config),
#         config)
#     process_runner.config = config_amended
#     snapshot.assert_match(unpack(config_amended), 'Config-amended')

#     # == 2. SETUP EXTERNAL STATE
#     external_state, start_day_out, end_day_out = setup_external_state(
#         config, data_location, Main_Overrides(start_day=start_day, end_day=end_day))
#     ROW_COUNT = (end_day - start_day + 1) * 24
#     assert external_state.sinB and len(
#         external_state.sinB) == ROW_COUNT and sum(external_state.sinB) > 0
#     assert external_state.VPD and len(
#         external_state.VPD) == ROW_COUNT and sum(external_state.VPD) > 0
#     assert external_state.P and len(external_state.P) == ROW_COUNT and sum(external_state.P) > 0
#     assert external_state.PAR and len(
#         external_state.PAR) == ROW_COUNT and sum(external_state.PAR) > 0
#     assert external_state.Idrctt and len(
#         external_state.Idrctt) == ROW_COUNT and sum(external_state.Idrctt) > 0
#     assert external_state.Idfuse and len(
#         external_state.Idfuse) == ROW_COUNT and sum(external_state.Idfuse) > 0
#     assert external_state.R and len(external_state.R) == ROW_COUNT and sum(external_state.R) > 0
#     assert external_state.PPFD and len(
#         external_state.PPFD) == ROW_COUNT and sum(external_state.PPFD) > 0
#     assert external_state.is_daylight and len(
#         external_state.is_daylight) == ROW_COUNT and external_state.is_daylight[0] is not None
#     assert external_state.VPD_dd and len(
#         external_state.VPD_dd) == ROW_COUNT and sum(external_state.VPD_dd) > 0
#     assert external_state.CO2 and len(
#         external_state.CO2) == ROW_COUNT and sum(external_state.CO2) > 0
#     assert external_state.RH and len(external_state.RH) == ROW_COUNT and sum(external_state.RH) > 0
#     assert external_state.Ts_C and len(
#         external_state.Ts_C) == ROW_COUNT and sum(external_state.Ts_C) > 0
#     assert external_state.eact and len(
#         external_state.eact) == ROW_COUNT and sum(external_state.eact) > 0
#     # assert external_state.uh and sum(external_state.uh) > 0 # Not currently set. In state instead?
#     # assert external_state.Rn and sum(external_state.Rn) > 0 # currently always 0.0
#     snapshot.assert_match(json.dumps(prep_data_for_snapshot(asdict(external_state)),
#                                      cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'external_state')
#     process_runner.external_state = external_state
#     assert process_runner.external_state.sinB[0] is not None
#     assert type(external_state.dd[0]) == int

#     # TODO: init params

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

#     snapshot.assert_match(json.dumps(asdict(
#         final_state.temporal), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - temporal')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_component[0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy_component')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_layer_component[0][0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy_layer_component')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_layers[0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy_layers')

#     snapshot.assert_match(process_runner.state_logs, 'logs')

#     # print output logs to csv
#     df_out = pd.DataFrame(process_runner.state_logs)
#     df_out.to_csv(SAVE_OUTPUTS_LOCATION + 'logs.csv')


# def test_processes_multiplicative(snapshot):
#     start_day = DEMO_START_DAY
#     end_day = DEMO_END_DAY
#     config_location = 'examples/spanish_wheat_multiplicative/old_configs/spanish_wheat_multiplicative_config_short_tests.json'
#     data_location = 'examples/spanish_wheat/inputs/spanish_wheat_data.csv'

#     # == 1. SETUP CONFIG

#     config = config_loader(config_location=config_location, config_type='json')
#     snapshot.assert_match(unpack(config), 'Config')
#     process_runner = ProcessRunner(config, DEBUG_MODE=True)

#     config_amended = process_runner.run_processes(
#         config_init_processes(config),
#         config)
#     process_runner.config = config_amended
#     snapshot.assert_match(unpack(config_amended), 'Config-amended')

#     # == 2. SETUP EXTERNAL STATE
#     external_state, start_day_out, end_day_out = setup_external_state(
#         config, data_location, Main_Overrides(start_day=start_day, end_day=end_day))

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
#     # assert external_state.uh and sum(external_state.uh) > 0 # Not currently set. In state instead?
#     # assert external_state.Rn and sum(external_state.Rn) > 0 # currently always 0.0
#     snapshot.assert_match(json.dumps(prep_data_for_snapshot(asdict(external_state)),
#                                      cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'external_state')
#     process_runner.external_state = external_state
#     assert process_runner.external_state.sinB[0] is not None

#     # TODO: init params

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

#     snapshot.assert_match(json.dumps(asdict(
#         final_state.temporal), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - temporal')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_component[0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy_component')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_layer_component[0][0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy_layer_component')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_layers[0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True),
#         'final_state - canopy_layers')

#     snapshot.assert_match(process_runner.state_logs, 'logs')

#     # print output logs to csv
#     df_out = pd.DataFrame(process_runner.state_logs)
#     df_out.to_csv(SAVE_OUTPUTS_LOCATION + 'logs_multiplicative.csv')

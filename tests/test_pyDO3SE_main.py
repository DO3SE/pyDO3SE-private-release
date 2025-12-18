# ''' testing pyDO3SE main'''

# from dataclasses import asdict
# import pytest
# import json
# from unittest.mock import patch
# import numpy as np
# from proflow import ProcessRunnerCls

# from data_helpers.cls_parsing import unpack
# from data_helpers.snapshot_helpers import prep_data_for_snapshot
# from data_helpers.encoders import AdvancedJsonEncoder


# from pyDO3SE.main import main
# from pyDO3SE.Config import Config_Shape
# from pyDO3SE.Config.Config_Shape import Config_Location
# from pyDO3SE.External_State import External_State_Shape


# from pyDO3SE.Pipelines import es_init_processes
# from pyDO3SE.Pipelines.es_init_processes import external_state_init_processes

# # from pyDO3SE.Pipelines.default_processes import hourly_processes, daily_start_processes, \
# #     daily_end_processes

# DEMO_START_DAY = 0
# DEMO_END_DAY = 50
# HOURS_IN_DAY = 24
# DAYS_IN_YEAR = 365
# DAYS_TO_RUN = DEMO_END_DAY - DEMO_START_DAY
# HOURS_IN_YEAR = HOURS_IN_DAY * DAYS_IN_YEAR
# HOURS_TO_RUN = HOURS_IN_DAY * DAYS_TO_RUN

# # @patch('pyDO3SE.plugins.gsto.photosynthesis.run_gsto_pn')


# # @patch('pyDO3SE.Pipelines.default_processes.ewert', side_effect=ewert.ewert)
# # @patch('pyDO3SE.Pipelines.default_processes.multiplicative',
# #        side_effect=multiplicative.multiplicative)
# # @patch('pyDO3SE.main.initialize_processes', side_effect=ProcessRunner.initialize_processes)
# # @patch('proflow.ProccessRunner.ProcessRunner.run_process', side_effect=ProcessRunner.run_process)
# @pytest.mark.skip(reason='Replaced with test_process_full_year.py')
# def test_main(
#         # mock_proccessRunner_run_process,
#         # mock_proccessRunner_initialize_processes,
#         # mock_process_calc_gsto_multiplicative,
#         # mock_process_calc_gsto_ewert,
#         snapshot,
# ):
#     ''' Functional test to ensure we can run the model from start to finish'''
#     output = main(
#         start_day=DEMO_START_DAY,
#         end_day=DEMO_END_DAY,
#         config_path='examples/spanish_wheat/configs/spanish_wheat_config.json',
#         data_location='examples/spanish_wheat/inputs/spanish_wheat_data.csv',
#         output_directory="demo_output",
#     )
#     final_state, process_runner = output

#     # Check that config was set

#     snapshot.assert_match(json.dumps(unpack(process_runner.config),
#                                      cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'Config')
#     snapshot.assert_match(json.dumps(prep_data_for_snapshot(asdict(process_runner.external_state)),
#                                      cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'external_state')

#     external_state = process_runner.external_state
#     snapshot.assert_match(
#         f'max - {max(external_state.O3_nmol)}, min - {min(external_state.O3_nmol)}',
#         'external state O3_nmol')
#     # snapshot.assert_match(
#     #     f'max - {max(external_state.td)}, min - {min(external_state.td)}', 'external state td')

#     # Check that runner was called
#     # assert mock_proccessRunner_initialize_processes.called

#     # Check the args that the runner was called with
#     # args, kwargs = mock_proccessRunner_initialize_processes.call_args
#     # assert len(args) == 1
#     # assert list(kwargs.keys()) == []

#     # all_processes_to_run = args[0]
#     # TODO: Failing because we need config to get hourly processes
#     # expected_total_process_count = \
#     #     len(hourly_processes(1)) * HOURS_IN_YEAR + \
#     #     len(daily_start_processes(1)) * DAYS_IN_YEAR + \
#     #     len(daily_end_processes(1)) * DAYS_IN_YEAR

#     # # assert expected_total_process_count == 18615
#     # assert len(all_processes_to_run) == expected_total_process_count
#     # assert len(all_processes_to_run) >= HOURS_TO_RUN + DAYS_TO_RUN + DAYS_TO_RUN
#     # assert all(isinstance(p, ProcessRunner.Process) for p in all_processes_to_run)

#     # Check that the function of the process is callable
#     # assert all(isinstance(p.func, Callable) for p in all_processes_to_run)

#     # Check that the processes were called
#     # assert mock_process_calc_gsto_ewert.called
#     # mock_process_calc_gsto_multiplicative.assert_not_called()

#     # check that the output is saved to the output directory
#     # TODO: Check output
#     assert final_state.temporal.hr == 23
#     assert final_state.temporal.dd == DAYS_TO_RUN + -1
#     snapshot.assert_match(json.dumps(asdict(final_state.temporal),
#                                      cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'temporal')
#     snapshot.assert_match(json.dumps(asdict(final_state.canopy),
#                                      cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'canopy')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_component[0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'canopy_component')
#     snapshot.assert_match(json.dumps(asdict(
#         final_state.canopy_layer_component[0][0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'canopy_layer_component')
#     snapshot.assert_match(json.dumps(
#         asdict(final_state.canopy_layers[0]), cls=AdvancedJsonEncoder, indent=4, sort_keys=True), 'canopy_layers')


# # @patch('pyDO3SE.Pipelines.es_init_processes.calc_O3_params',
# #        side_effect=es_init_processes.calc_O3_params)
# # # @patch('pyDO3SE.Pipelines.param_init_processes.O3_helpers',
# # #        side_effect=O3_helpers.O3_ppb_to_nmol)
# # def test_setup_external_state(
# #     # mock_calc_O3_params,
# #     mock_O3_pbnmol,
# # ):
# #     config = Config_Shape(Location=Config_Location(10, 10))
# #     external_state: External_State_Shape = External_State_Shape(
# #         O3=np.arange(130),
# #         dd=np.arange(130),
# #         hr=np.arange(130),
# #         Ts_C=np.arange(130),
# #         P=np.full_like(np.zeros(130), 94.1),
# #         precip=np.arange(130),
# #         u=np.arange(130),
# #         VPD=np.arange(130),
# #         PAR=np.arange(130),
# #     )
# #     process_runner = ProcessRunnerCls.ProcessRunner(
# #         external_state_in=external_state, config_in=config)

# #     start_day = 0
# #     end_day = 5
# #     processes_to_init = external_state_init_processes(start_day, end_day, config)
# #     run_init = process_runner.initialize_processes(processes_to_init)
# #     external_state_initialized: External_State_Shape = run_init(external_state)

# #     # assert mock_calc_O3_params.called
# #     assert mock_O3_pbnmol.called

# #     assert external_state_initialized.O3_nmol[0:3] == [0, 41.131684654895416, 81.66973142915009]
# #     assert external_state_initialized.eact[0:3] == [
# #         0.611, -0.34307580354071987, -1.2941275110314323]
# #     assert external_state_initialized.RH[0:3] == [1.0, -0.5222456493303876, -1.8333729267766652]
# #     assert external_state_initialized.esat[0:3] == [0.611, 0.6569241964592801, 0.7058724889685677]

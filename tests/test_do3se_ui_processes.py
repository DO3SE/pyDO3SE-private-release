# import os
# import pytest
# from tests.test_helpers import long_test
# from pyDO3SE.main import export_output
# from pyDO3SE.setup_model import setup_model
# from pyDO3SE.run_model import run_model
# from pyDO3SE.Pipelines.do3se_ui_processes import full_model_processes

# SAVE_OUTPUTS_LOCATION = os.path.dirname(os.path.realpath(__file__)) + '/output/'


# @long_test
# @pytest.mark.skip(reason="No longer relevant")
# def test_do3se_ui_processes():
#     demo_config = "./examples/spanish_wheat_multiplicative/configs/spanish_wheat_multiplicative_config.json"
#     demo_data = "./examples/spanish_wheat/inputs/spanish_wheat_data.csv"

#     # MODEL SETUP
#     [
#         config,
#         external_state,
#         initial_state,
#         model_processes,
#     ] = setup_model(demo_config, demo_data)
#     model_processes = full_model_processes(config, 1, 365)
#     final_state, output_logs = run_model(initial_state, config, external_state, model_processes)

#     output_fields = ['gsto_l', 'A_n', 'PARshade', 'PARsun']

#     export_output(
#         output_logs,
#         final_state,
#         external_state,
#         SAVE_OUTPUTS_LOCATION,
#         config,
#         output_fields=output_fields,
#         runid=0,
#         log_notes="test_do3se_ui_processes",
#     )

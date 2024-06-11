# """A set of tests that run the hourly model then compare the output against the previous version."""
# import warnings
# import pandas as pd
# import os
# import pytest
# from dataclasses import asdict
# from data_helpers.diff import diff

# from pyDO3SE import main

# example_folders = [f"examples/parial_runs/{i}" for i in os.listdir('examples/hourly_runs')]


# @pytest.fixture(scope="class", params=example_folders)
# def setup_test_run(request):
#     project_dir = request.param
#     config_file = f"{project_dir}/processed_config.json"
#     config_file = f"{project_dir}/processed_config.json"
#     initial_state_file = f"{project_dir}/initial_state.json"
#     input_data_file = f"{project_dir}/hourly_data.csv"
#     # output_data_file = f"{project_dir}/hourly_output.csv"
#     # output_state_file = f"{project_dir}/final_state.json"
#     output_directory = f"{project_dir}/outputs"
#     request.cls.project_dir = project_dir

#     errors = []
#     try:
#         # == Run First pass
#         output_fields = ['gsto_canopy', 'td_dd', 'lai', 'pody', 'fst', 'canopy_height']
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             out = main.main_partial(
#                 processed_config_path=config_file,
#                 data_path=input_data_file,
#                 output_directory=output_directory,
#                 initial_state_file=initial_state_file,
#                 output_fields=output_fields,
#                 runid=f"{project_dir}",
#                 runnotes="Running from test",
#                 log_level=0,
#                 # base_config_file=base_demo_config,
#             )

#         final_state, output_logs, processed_config, initial_state = out
#         request.cls.output = out
#         request.cls.logs = pd.DataFrame(output_logs)
#         # Make a copy for the version
#     except Exception as e:
#         errors.append((f"Project dir: {project_dir} failed", e))

#     if len(errors) > 0:
#         for m, e in errors:
#             print(m)
#         print(errors)
#         raise errors[0][1]
#     run_type = processed_config.Land_Cover.parameters[0].gsto.method
#     request.cls.run_type = run_type


# @pytest.mark.usefixtures('setup_test_run')
# class TestHourlyRun:

#     def test_should_run_without_errors(self):
#         pass

#     def test_compare_initial_to_final(self, snapshot):
#         final_state, output_logs, processed_config, initial_state = self.output
#         initial_state_copy = {k: v for k, v in asdict(initial_state).items() if k != "prev_hour"}
#         final_state_copy = {k: v for k, v in asdict(final_state).items() if k != "prev_hour"}
#         compared = diff("model_state", initial_state_copy, final_state_copy)
#         snapshot.assert_match({
#             "total_changes": len(compared),
#             "changes": compared
#             }, "Changes")



# """A set of tests that run a short version of the model to check timings."""
# from timeit import timeit
# import pytest
# import os
# from shutil import copyfile
# from tests.test_helpers import long_test

# from pyDO3SE.Output.run_comparisons import create_comparison_graphs, get_output_files

# from pyDO3SE import main

# FIELDS_TO_SAVE = [
#     'gsto_l',
#     'gsto_bulk',
#     'gsto',
#     'rsto',
#     'rsto_c',
#     'rsto_l',
#     "td_dd",
#     'dvi',
#     'PARsun',
#     'PARshade',
#     'fO3_d',
#     'fst',
#     'A_n',
#     'lai',
#     'sai',
#     'o3_ppb',
#     'o3_ppb_i',
#     'o3_nmol_m3',
#     'canopy_height',
#     'swp',
#     'ustar',
#     'ra',
#     'rb',
#     'rsur',
#     'rinc',
#     'f_LS',
#     'f_VPD',
#     'LAIsunfrac',
# ]

# @pytest.mark.skip(reason="Needs updating")
# def test_demo_run_short_multi():
#     config_dir = "./examples/short/configs"
#     input_dir = "./examples/short/inputs"
#     output_dir = f"./examples/short/outputs"

#     config_file_type = 'json'

#     config_files = [c for c in os.listdir(config_dir) if len(c.split('.')) == 2 and c.split('.')[
#         1] == config_file_type]
#     input_files = [f for f in os.listdir(input_dir) if f.split('.')[-1] == 'csv']

#     for config_file in config_files:
#         output_dir_run = output_dir + '/' + config_file.split('.')[0]
#         os.makedirs(output_dir_run, exist_ok=True)
#         for input_file in input_files:
#             final_state, output_logs, config_final, initial_state = main.main(
#                 config_location=config_dir + '/' + config_file,
#                 data_location=input_dir + '/' + input_file,
#                 runid=input_file.replace('.csv', ''),
#                 output_directory=output_dir_run,
#                 output_fields=[],
#                 runnotes="Running from test",
#                 # output_results_only=True,
#                 start_day=100,
#                 end_day=104,
#             )

#     data_files = get_output_files(output_dir)
#     create_comparison_graphs(
#         data_files,
#         output_dir,
#         FIELDS_TO_SAVE,
#         './examples/short/comparisons',
#         use_daily_average=False,
#     )


# config_dir = "./examples/short/configs"
# input_dir = "./examples/short/inputs"

# config_file_type = 'json'

# config_files = [c for c in os.listdir(config_dir) if len(c.split('.')) == 2 and c.split('.')[
#     1] == config_file_type]
# input_files = [f for f in os.listdir(input_dir) if f.split('.')[-1] == 'csv']


# @pytest.mark.mark(reason="Needs updating")
# @pytest.mark.parametrize('config_file', config_files)
# @pytest.mark.parametrize('input_file', input_files)
# def test_benchmark_run_short_multi(benchmark, config_file, input_file):
#     result_time = benchmark(lambda: main.main(
#         config_location=config_dir + '/' + config_file,
#         data_location=input_dir + '/' + input_file,
#         runid=input_file.replace('.csv', ''),
#         output_fields=[],
#         runnotes="Running from test",
#         start_day=100,
#         end_day=104,
#     ))

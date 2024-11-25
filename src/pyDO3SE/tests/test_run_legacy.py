# """A set of tests that run the full model then compare the output against the previous version."""
# from tests.test_helpers import long_test
# import warnings
# import pandas as pd
# import os
# import pytest
# from shutil import copyfile
# from proflow.errors import Run_Process_Error
# from do3se_phenology.state import LeafPhenologyStage, PhenologyStage
# from pyDO3SE.Config.config_loader import config_loader
# from pyDO3SE.Output.process_outputs import (
#     dump_model_processes_info_to_string,
#     dump_config_to_file_json,
#     dump_state_to_file,
# )
# from pyDO3SE.Pipelines.default_processes import (
#     full_model_processes,
# )
# from pyDO3SE import main


# example_folders = [
#     './examples/legacy',
# ]
# setups = [
#     ','.join([project_dir, cf.split('.json')[0], f.split('.csv')[0]])
#     for project_dir in example_folders
#     for cf in sorted(os.listdir(f"{project_dir}/configs"))
#     for f in os.listdir(f"{project_dir}/inputs")
# ]


# @pytest.fixture(scope="class", params=setups)
# def setup_test_run_and_compare(request):
#     project_dir, config_file, input_file = request.param.split(',')
#     skip_outputs = not request.config.getoption("--skip-outputs") == "False"
#     request.cls.project_dir = project_dir
#     errors = []
#     try:
#         runid = "demo_test_output_comparison"
#         project_paths = main.get_project_paths(project_dir)
#         run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)

#         # Create output dir
#         main.create_run_path_directories(run_paths)

#         demo_config = run_paths.config_path
#         base_demo_config = project_paths.base_config_path
#         demo_output_dir = run_paths.output_directory

#         # Copy config to this directory
#         copyfile(demo_config, f'{demo_output_dir}/config.json')

#         # Print out model processes
#         hours = list(range(24))
#         full_model_processes_out = full_model_processes(
#             config_loader(demo_config, base_demo_config), hours)
#         flattened_process_comments = dump_model_processes_info_to_string(
#             full_model_processes_out, True)

#         with open(f"{demo_output_dir}/processes.txt", 'w') as f:
#             f.write('\n'.join(flattened_process_comments))
#         dump_config_to_file_json(config_loader(demo_config, base_demo_config),
#                                  f"{demo_output_dir}/processed_config_init.json")
#         # ======

#         fields_to_graph = ','.join(['gsto_canopy', 'lai', 'pody', 'fst', 'f_phen', 'leaf_f_phen'])

#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             out = main.single(
#                 config_file=run_paths.config_path,
#                 data_file=run_paths.input_data_file_path,
#                 output_directory=demo_output_dir,
#                 base_config_file=base_demo_config,
#                 plot_fields=fields_to_graph,
#                 runid="test_run_legacy",
#                 verbose=2,
#             )

#         final_state, output_logs, final_config, initial_state = out
#         request.cls.output = out
#         request.cls.logs = pd.DataFrame(output_logs)

#     except Exception as e:
#         errors.append((f"Config: {config_file} with input {input_file} failed", e))

#     if len(errors) > 0:

#         for m, e in errors:
#             print(m)
#             if(type(e) == Run_Process_Error):
#                 dump_state_to_file(e.state, 'failed_state.json')
#         print(errors)

#         raise errors[0][1]
#     run_type = final_config.Land_Cover.parameters[0].gsto.method
#     start_day = final_config.Location.start_day or 0
#     end_day = final_config.Location.end_day or 365
#     request.cls.run_type = run_type

#     # if not skip_outputs:
#     #     # Compare outputs
#     #     if run_type == "multiplicative":
#     #         fields_to_graph = multip_output_fields
#     #     elif run_type == "photosynthesis":
#     #         fields_to_graph = pn_output_fields
#     #     else:
#     #         raise Exception("Invalid run type")
#     #     output_dir = f"{project_dir}/outputs"
#     #     data_files = list(sorted(get_output_files(output_dir, "latest"), key=lambda f: f[1]))
#     #     with warnings.catch_warnings():
#     #         warnings.simplefilter("ignore")
#     #         create_comparison_graphs(
#     #             data_files,
#     #             fields_to_graph,
#     #             f"{project_dir}/comparisons",
#     #             use_versioned_outfile=False,
#     #             use_versioned_outdir=False,
#     #             log_level=0,
#     #             start_day=start_day,
#     #             end_day=end_day,
#     #         )


# @pytest.mark.usefixtures('setup_test_run_and_compare')
# class TestRunAndCompare:

#     @long_test
#     def test_main(self):
#         pass

#     @long_test
#     def test_should_have_increased_pody(self):
#         final_state, output_logs, final_config, initial_state = self.output
#         assert final_state.canopy_component_population[0][-1].POD_Y > 0

#     @long_test
#     def test_plant_should_have_emerged(self):
#         final_state, output_logs, final_config, initial_state = self.output
#         # assert output_logs[0]['phenology_stage'] < PhenologyStage.EMERGED
#         assert final_state.canopy_component[0].phenology.phenology_stage > PhenologyStage.NOT_SOWN
#         assert final_state.canopy_component[0].phenology.phenology_stage > PhenologyStage.SOWN
#         assert final_state.canopy_component[0].phenology.phenology_stage > PhenologyStage.EMERGED
#         assert final_state.canopy_component[0].phenology.phenology_stage > PhenologyStage.ASTART
#         assert final_state.canopy_component[0].phenology.phenology_stage == PhenologyStage.HARVEST

#     @long_test
#     def test_flag_leaf_should_have_finished_growing(self):
#         final_state, output_logs, final_config, initial_state = self.output
#         assert final_state.canopy_component_population[0][-1].phenology.phenology_stage > LeafPhenologyStage.GROWING

#     @long_test
#     def test_should_have_increased_gsto(self):
#         output_logs = self.logs
#         assert max(output_logs['gsto']) > 0

#     @long_test
#     def test_flag_leaf_should_have_emerged(self):
#         final_state, output_logs, final_config, initial_state = self.output
#         assert final_state.canopy_component[0].total_emerged_leaf_populations > 0

#     @long_test
#     def test_should_have_increased_f_phen(self):
#         output_logs = self.logs
#         assert max(output_logs['f_phen']) > 0

#     @long_test
#     def test_should_have_increased_leaf_f_phen(self):
#         output_logs = self.logs
#         assert max(output_logs['leaf_f_phen']) > 0

#     @long_test
#     def test_should_have_increased_lai(self):
#         output_logs = self.logs
#         assert max(output_logs['canopy_lai']) > 0

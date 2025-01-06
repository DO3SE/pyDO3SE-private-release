"""Generate initial state for partial runs."""

# %%
import warnings
from pyDO3SE.Model_State.Model_State import Model_State_Shape

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main

source_dir = "examples/bangor_2015_multilayer"
runid = "get_states"

project_paths = main.get_project_paths(source_dir)
run_paths = main.get_run_paths(runid, project_paths, "bangor_wheat_single_pop", "bangor_2015_hb_ww")

# Create output dir
main.create_run_path_directories(run_paths)

output_options = output_results_only_options()
output_options.save_hourly_output_data = False

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    out = main.single(
        config_file=run_paths.config_path,
        data_file=run_paths.input_data_file_path,
        output_directory=run_paths.output_directory,
        base_config_file=project_paths.base_config_path,
        plot_fields=None,
        runid=runid,
        verbose=2,
        output_options=output_options,
        overrides = ["dump_state_n=100"],
    )


# %%
import os
from pyDO3SE.Model_State.model_state_loader import model_state_loader
from pyDO3SE.Output.process_outputs import dump_state_to_file
state_in_dir = "examples/bangor_2015_multilayer/runs/get_states/bangor_wheat/bangor_2015_hb_ww/outputs/get_states/running_state"
state_out_dir = "tests/key_processes/flag_only_pngsto/initial_state"
for sf in os.listdir(state_in_dir):
    model_state: Model_State_Shape = model_state_loader(f"{state_in_dir}/{sf}")
    model_state.temporal.row_index = 0
    model_state.temporal.dd = 0
    dvi = model_state.canopy_component[0].dvi
    dvi_string = "{:.2f}".format(dvi).replace('.', '_')
    print(dvi_string)
    out_filename = f"{state_out_dir}/dump_state_dvi_{dvi_string}.json"
    dump_state_to_file(model_state, out_filename)

# %%

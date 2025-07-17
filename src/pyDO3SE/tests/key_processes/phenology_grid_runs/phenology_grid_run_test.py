"""We can run just the phenology module on gridded data"""
import warnings
import os
import json
import pytest
from typing import Optional, Dict
from data_helpers.encoders import AdvancedJsonEncoder
from dataclasses import asdict
from datetime import datetime
from pyDO3SE.Model_State.model_state_loader import dump_state_to_file, model_state_loader_quick
from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.setup_model import (
    Main_Overrides,
)
from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.Grid_Model import setup_grid_model
from pyDO3SE.Grid_Model import run_grid_model
from pyDO3SE.util.logger import Logger, generate_run_notes
from pyDO3SE.Pipelines.phenology_only_processes import get_row_processes_hourly

output_fields = ['dd', 'hr', 'gsto_canopy', 'td_dd',
                 'canopy_lai', 'pody', 'fst', 'canopy_height', 'micro_u', 'dvi', 'td_v']
multi_file_netcdf = False
regex_multi_file_filter: Optional[str] = None
netcdf_loader_kwargs = {}
seperate_state_path = True


def run_phenology_on_gridded_data(
    project_dir: str,
    runid: str,
    config_id: str,
    start_index: Optional[int] = None,
    end_index: Optional[int] = None,
    init: bool = True,
):
    project_paths = setup_grid_model.get_grid_project_paths(project_dir, runid)
    run_paths = setup_grid_model.get_grid_run_paths(project_paths, config_id)
    loaded_run_files = setup_grid_model.load_grid_run_files(project_paths, run_paths)
    setup_grid_model.create_grid_run_path_directories(run_paths)
    grid_coords, grid_x_size, grid_y_size = setup_grid_model.get_grid_coords_from_file(
        run_paths.run_mask_path)
    e_state_overrides_dataset = xr.open_dataset(project_paths.e_state_overrides_file_path)
    logger_main = Logger(0, None, set_as_default=True)

    output_shape = (grid_x_size, grid_y_size)

    if init:
        initialized_config_gen, initialized_state_gen = setup_grid_model.init_grid_model(
            config=loaded_run_files.config,
            state=loaded_run_files.state,
            e_state_overrides_dataset=e_state_overrides_dataset,
            e_state_overrides_field_map=loaded_run_files.e_state_overrides_field_map,
            grid_coords=grid_coords,
            logger=logger_main,
            debug=True,
        )
        setup_grid_model.save_configs_from_generator(
            initialized_config_gen,
            grid_coords,
            run_paths.processed_configs_dir,
        )
        setup_grid_model.save_state_from_generator(
            initialized_state_gen,
            grid_coords,
            run_paths.live_state_dir,
        )
    model_processes = get_row_processes_hourly(loaded_run_files.config, list(range(24)))

    errors = []
    model_output = None
    final_states: Dict[str, Model_State_Shape] = {}

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            start_time = datetime.now()
            model_output = run_grid_model.main_grid_seq_per_config(
                project_paths=project_paths,
                run_paths=run_paths,
                loaded_run_files=loaded_run_files,
                grid_coords=grid_coords,
                output_shape=output_shape,
                runnotes="",
                output_fields=output_fields,
                seperate_live_state=seperate_state_path,
                multi_file_netcdf=multi_file_netcdf,
                regex_multi_file_filter=regex_multi_file_filter,
                logger=logger_main,
                netcdf_loader_kwargs=netcdf_loader_kwargs,
                overrides=Main_Overrides(
                    model_processes=model_processes
                ),
                start_input_index=start_index,
                end_input_index=end_index,
                return_outputs=True,
                debug=True
            )
            end_time = datetime.now()
            duration = end_time - start_time

            with open(f'{run_paths.output_data_dir}/notes.log', 'w') as f:
                log_notes = generate_run_notes(
                    runnotes='',
                    time_taken=duration,
                    # time_taken_setup=setup_duration,
                    # config_version=config_version,
                    # model_version=model_version,
                    errors=errors,
                )
                f.write("\n".join(log_notes))

        # Save a human readable copy of the final state
        for x, y in grid_coords:
            file_name = f"{x}_{y}"
            final_state = model_state_loader_quick(
                f"{run_paths.live_state_dir}/{file_name}.state")
            final_states[file_name] = final_state
            dump_state_to_file(model_state_loader_quick(
                f"{run_paths.live_state_dir}/{file_name}.state"), f"{run_paths.final_state_dir}/{file_name}.json")

    except Exception as e:
        errors.append((f"Project dir: {project_paths.project_dir} failed", e))

    if len(errors) > 0:
        for m, e in errors:
            print(m)
        print(errors)
        raise errors[0][1]
    print("Complete")
    return model_output, final_states


def test_run_phenology_on_gridded_data():
    project_dir = "examples/net_cdf/full_season_monthly"
    runid = "phenology_grid_run"
    model_output, final_states = run_phenology_on_gridded_data(project_dir, runid, "bangor_wheat")
    input_file_count = len(os.listdir(f"{project_dir}/inputs"))
    assert model_output is not None
    time_len = 743 # hours per file?
    input_file_count = 7
    grid_coord_length = 2
    assert len(model_output) == input_file_count
    assert len(model_output[0]) == grid_coord_length
    # 4 == [Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]
    Coords, ModelGridCellOutput, Model_State_Shape, OutputFields = model_output[0][0]
    assert len(ModelGridCellOutput) == len(OutputFields)
    assert len(ModelGridCellOutput[0]) == time_len
    # DVI should reach 2.0
    # final_state: Model_State_Shape = final_states.values()[0]
    # assert final_state.canopy_component[0].dvi >= 2.0
    output_fields = model_output[0][0][3]
    dvi_index = output_fields.index('dvi')
    assert model_output[-1][0][1][-1][dvi_index] >= 2.0


@pytest.mark.skip(reason="test currently fails")
def test_continue_run_phenology_on_gridded_data():
    """Test that we can continue a run that has already completed some of the input files."""
    project_dir = "examples/net_cdf/full_season_monthly"
    runid = "phenology_grid_run_continue_full"
    model_output, final_states = run_phenology_on_gridded_data(project_dir, runid, "bangor_wheat")
    time_len = 743 # hours per file?
    output_field_count = 11
    input_file_count = 7
    grid_coord_length = 2
    assert len(model_output) == input_file_count
    assert len(model_output[0]) == grid_coord_length
    # 4 == [Coords, ModelGridCellOutput, Model_State_Shape, OutputFields]
    assert len(model_output[0][0]) == 4
    assert len(model_output[0][0][1]) == output_field_count
    assert len(model_output[0][0][1][0]) == time_len

    # assert out[0][0].get('dd').shape == (grid_len_x, grid_len_y, time_len)
    start_index = 2
    end_index = 3
    runid = "phenology_grid_run_continue_partial"
    out_partial_a, final_states_a = run_phenology_on_gridded_data(
        project_dir, runid, "bangor_wheat", start_index=0, end_index=start_index - 1, init=True)
    out_partial_b, final_states_b = run_phenology_on_gridded_data(
        project_dir, runid, "bangor_wheat", start_index=start_index, init=False)

    final_state_left = json.dumps(
        asdict(final_states['0_0']),
        cls=AdvancedJsonEncoder, indent=4, sort_keys=True)
    final_state_right = json.dumps(
        asdict(final_states_b['0_0']),
        cls=AdvancedJsonEncoder, indent=4, sort_keys=True)
    assert final_state_left == final_state_right

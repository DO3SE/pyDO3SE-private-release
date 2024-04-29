"""Test running a few lines of data to make sure soil moisture is correct."""
from collections import namedtuple
from math import isclose
from pathlib import Path
from copy import deepcopy
import pytest
import warnings
import pandas as pd
from pyDO3SE.Config.Config_Shape import Config_Shape

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main

project_dir = "tests/key_processes/ozone_damage"


def run_with_config(runid: str, project_dir: Path, config_file: str, input_file: str, input_state_file: str, **kwargs):
    project_paths = main.get_project_paths(project_dir)
    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)
    project_paths = project_paths._replace(
        base_state_path=f"{project_dir}/initial_state/{input_state_file}.json")

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
            initial_state_path=project_paths.base_state_path,
            plot_fields=None,
            runid=runid,
            verbose=2,
            output_options=output_options,
        )
    return out


RunSetup = namedtuple('RunSetup', ['runid', 'config_file',
                                   'input_file', 'input_state', 'overrides'])

setups = [
    RunSetup("start_second_emergence", "second_emergence", "three_days", 'start', {}),
    RunSetup("dvi_0_91_second_emergence", "second_emergence",
             "three_days", 'dump_state_dvi_0_91', {}),
    RunSetup("default_dvi_0_91", "default", "three_days", 'dump_state_dvi_0_91', {}),
    RunSetup("default_dvi_1_41", "default", "three_days", 'dump_state_dvi_1_41', {}),
    RunSetup("default_dvi_1_86", "default", "three_days", 'dump_state_dvi_1_86', {}),
    RunSetup("long_ozone_damage_dvi_1_86", "long_ozone_damage",
             "three_days", 'dump_state_dvi_1_86', {}),
    RunSetup("short_ozone_damage_dvi_1_86", "short_ozone_damage",
             "three_days", 'dump_state_dvi_1_86', {}),
    RunSetup("default_dvi_2_00", "default", "three_days", 'dump_state_dvi_2_00', {}),
    RunSetup("long_ozone_damage_dvi_2_00", "long_ozone_damage",
             "three_days", 'dump_state_dvi_2_00', {}),
    RunSetup("short_ozone_damage_dvi_2_00", "short_ozone_damage",
             "three_days", 'dump_state_dvi_2_00', {}),
]

setupsMap = {s.runid: s for s in setups}
all_setups = [s.runid for s in setups]


@pytest.fixture(scope="class")
def ozone_damage_test_run(request):
    request.cls.output = {}

    for runid, config_file, input_file, input_state, overrides in setups:
        try:
            out = run_with_config(
                runid=runid,
                project_dir=project_dir,
                config_file=config_file,
                input_file=input_file,
                input_state_file=input_state,
                **overrides,
            )
        except Exception as e:
            print(f"Failed to run {runid}")
            raise e

        final_state, output_logs, final_config, initial_state, external_state = out
        request.cls.output[runid] = {}
        request.cls.output[runid]['out'] = out
        request.cls.output[runid]['hourly_output'] = pd.DataFrame(output_logs)
        request.cls.output[runid]['final_config'] = final_config
        request.cls.output[runid]['final_state'] = final_state
        request.cls.output[runid]['initial_state'] = initial_state


@pytest.mark.usefixtures('ozone_damage_test_run')
class TestRunAndCompare:

    def test_preruns_run_without_error(self):
        pass

    @pytest.mark.parametrize('runid', ['default_dvi_0_91'])
    def test_should_be_1_when_dvi_less_than_1(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        f_LS = hourly_output['f_LS'].values
        assert f_LS[0] == 1

    @pytest.mark.parametrize('runid', ['default_dvi_0_91'])
    def test_should_increase_td(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        td = hourly_output['td'].values
        assert td[0] < td[-1]

    @pytest.mark.parametrize('runid', ['default_dvi_1_86', 'long_ozone_damage_dvi_1_86', 'short_ozone_damage_dvi_1_86'])
    def test_should_reduce_f_ls(self, runid):
        """f_LS should reduce as DVI increases after DVI = approx 1.7."""
        hourly_output = self.output[runid]['hourly_output']
        f_LS = hourly_output['f_LS'].values
        td_dd = hourly_output['td_dd'].values
        print(td_dd)
        assert f_LS[0] > f_LS[-1]

    @pytest.mark.parametrize('runid', ['default_dvi_2_00'])
    def test_thermal_time_should_reach_end_of_plant_dvi(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_config: Config_Shape = self.output[runid]['final_config']
        td_end = hourly_output['td_dd'].values[-1]
        td_harv = final_config.Land_Cover.parameters[0].phenology.key_lengths_td.emerg_to_end
        assert td_end > td_harv

    @pytest.mark.parametrize('runid', ['default_dvi_2_00', 'long_ozone_damage_dvi_2_00', 'short_ozone_damage_dvi_2_00'])
    def test_f_ls_should_reach_0_at_dvi_2(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        f_LS = hourly_output['f_LS'].values
        assert f_LS[-1] == 0

    @pytest.mark.parametrize('runid', ['start_second_emergence'])
    def test_should_not_accumulate_fst_until_second_emergence(self, runid):
        final_state = self.output[runid]['final_state']
        final_config = self.output[runid]['final_config']
        leaf_emerg_to_leaf_fst_acc = final_config.Land_Cover.parameters[
            0].phenology.key_lengths_flag_leaf_td.leaf_emerg_to_leaf_fst_acc
        assert leaf_emerg_to_leaf_fst_acc > 0

        assert final_state.canopy_component[0].td_dd_leaf_pops[0] < leaf_emerg_to_leaf_fst_acc
        assert final_state.canopy_component_population[0][0].O3up_acc == 0

    @pytest.mark.parametrize('runid', ['dvi_0_91_second_emergence',
                                       "default_dvi_0_91",
                                       "default_dvi_1_41",
                                       "default_dvi_1_86",
                                       "long_ozone_damage_dvi_1_86",
                                       "short_ozone_damage_dvi_1_86",
                                       "default_dvi_2_00",
                                       "long_ozone_damage_dvi_2_00",
                                       "short_ozone_damage_dvi_2_00",
                                       ])
    def test_should_accumulate_fst_after_second_emergence(self, runid):
        final_state = self.output[runid]['final_state']
        initial_state = self.output[runid]['initial_state']
        final_config = self.output[runid]['final_config']
        leaf_emerg_to_leaf_fst_acc = final_config.Land_Cover.parameters[
            0].phenology.key_lengths_flag_leaf_td.leaf_emerg_to_leaf_fst_acc

        assert final_state.canopy_component[0].td_dd_leaf_pops[0] > leaf_emerg_to_leaf_fst_acc
        assert final_state.canopy_component_population[0][0].O3up_acc > initial_state.canopy_component_population[0][0].O3up_acc

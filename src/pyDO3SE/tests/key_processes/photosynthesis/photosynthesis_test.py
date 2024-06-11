"""Test running a few lines of data to make sure soil moisture is correct."""
from collections import namedtuple
from pathlib import Path
import pytest
import warnings
import pandas as pd

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main

project_dir = "tests/key_processes/photosynthesis"


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
    RunSetup("default_dvi_1_03", "default", "three_days", 'dump_state_dvi_1_03', {}),
    RunSetup("default_dvi_1_94", "default", "three_days", 'dump_state_dvi_1_94', {}),
]

setupsMap = {s.runid: s for s in setups}
all_setups = [s.runid for s in setups]


@pytest.fixture(scope="class")
def photosynthesis_test_run(request):
    request.cls.output = {}

    for runid, config_file, input_file, input_state, overrides in setups:
        out = run_with_config(
            runid=runid,
            project_dir=project_dir,
            config_file=config_file,
            input_file=input_file,
            input_state_file=input_state,
            **overrides,
        )

        final_state, output_logs, final_config, initial_state, external_state = out
        request.cls.output[runid] = {}
        request.cls.output[runid]['out'] = out
        request.cls.output[runid]['hourly_output'] = pd.DataFrame(output_logs)


@pytest.mark.usefixtures('photosynthesis_test_run')
class TestRunAndCompare:

    def test_preruns_run_without_error(self):
        pass


    @pytest.mark.parametrize('runid', ['default_dvi_1_03'])
    def test_should_have_sun_frac_below_0(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        LAIsunfrac_top = hourly_output['LAIsunfrac_top'].values
        assert any(i < 1.0 for i in LAIsunfrac_top)


    @pytest.mark.parametrize('runid', ['default_dvi_1_03'])
    def test_should_output_canopy_sunlit_gsto(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        gsto_l_sunlit = hourly_output['gsto_l_sunlit'].values
        gsto_l = hourly_output['gsto_l'].values
        assert any(a < b for a, b in zip(gsto_l, gsto_l_sunlit))

    @pytest.mark.parametrize('runid', ['default_dvi_1_03'])
    def test_should_output_canopy_sunlit_An(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        A_n_sunlit = hourly_output['A_n_sunlit'].values
        A_n = hourly_output['A_n'].values
        in_sun_index = next(i for i, v in enumerate(A_n) if v > 0)
        assert A_n_sunlit[in_sun_index] > A_n[in_sun_index]

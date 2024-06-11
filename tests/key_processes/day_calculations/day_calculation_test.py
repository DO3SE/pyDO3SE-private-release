"""Test running a few lines of data to make sure ozone deposition is correct."""
from collections import namedtuple
from math import isclose
from pathlib import Path
import pytest
import warnings
import pandas as pd

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE import main


def run_with_config(runid: str, project_dir: Path, config_file: str, input_file: str, **kwargs):
    project_paths = main.get_project_paths(project_dir)
    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)

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
        )
    return out


project_dir = "tests/key_processes/day_calculations"

RunSetup = namedtuple('RunSetup', ['runid', 'config_file', 'input_file', 'overrides'])

setups = [
    RunSetup("use_input_dd", "use_input_dd", "three_days", {}),
    RunSetup("bangor_wheat", "bangor_wheat", "three_days", {}),
    RunSetup("sparse_simple", "sparse", "sparse", {}),

]


@pytest.fixture(scope="class")
def day_calculation_test_run(request):
    request.cls.output = {}

    for runid, config_file, input_file, overrides in setups:
        out = run_with_config(
            runid=runid,
            project_dir=project_dir,
            config_file=config_file,
            input_file=input_file,
            **overrides,
        )

        final_state, output_logs, final_config, initial_state, external_state = out
        request.cls.output[runid] = {}
        request.cls.output[runid]['out'] = out
        request.cls.output[runid]['hourly_output'] = pd.DataFrame(output_logs)


@pytest.mark.usefixtures('day_calculation_test_run')
class TestRunAndCompare:

    def test_preruns_run_without_error(self):
        pass

    @pytest.mark.parametrize('runid', ['sparse_simple'])
    def test_should_calculate_dd_correctly(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        dd = hourly_output['dd'].values
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        assert dd[0] is not None
        assert dd[-1] is not None
        assert max(dd) > 0
        assert min(dd) < 365
        assert final_state.temporal.dd < 365


    @pytest.mark.parametrize('runid', ['use_input_dd', 'bangor_wheat'])
    def test_should_calculate_dd_correctly(self, runid):
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        external_state: External_State_Shape = external_state
        hourly_output = self.output[runid]['hourly_output']
        dd = hourly_output['dd'].values
        assert external_state.dd[0] < external_state.dd[-1]
        assert dd[0] < dd[-1]
        assert all([isclose(a, b, abs_tol=1e1) for a, b in zip(dd, external_state.dd)])
        assert final_state.temporal.dd < 365

    @pytest.mark.parametrize('runid', ['use_input_dd', 'bangor_wheat'])
    def test_should_return_correct_number_of_rows_of_data(self, runid):
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        hourly_output = self.output[runid]['hourly_output']
        dd = hourly_output['dd'].values
        assert len(dd) == len(external_state.dd) == 3 * 24

    @pytest.mark.parametrize('runid', ['use_input_dd', 'bangor_wheat'])
    def test_should_calculate_hr_correctly(self, runid):
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        external_state: External_State_Shape = external_state
        hourly_output = self.output[runid]['hourly_output']
        hr_values = hourly_output['hr'].values
        dd = external_state.dd
        assert external_state.hr[0] < external_state.hr[-1]
        assert len([hr for hr in external_state.hr if hr == 0]) == dd[-1] - dd[0] + 1
        assert len([hr for hr in hr_values if hr == 0]) == dd[-1] - dd[0] + 1

"""Test running a few lines of data to make sure ozone deposition is correct."""
from pathlib import Path
import pytest
import warnings
import pandas as pd

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main


def run_with_config(runid: str, project_dir: Path, config_file: str, input_file: str, **kwargs):
    project_paths = main.get_project_paths(project_dir)
    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)

    # Create output dir
    main.create_run_path_directories(run_paths)

    output_options = output_results_only_options()
    output_options.save_hourly_output_data = False

    kwargs_all = {
        **dict(
            config_file=run_paths.config_path,
            data_file=run_paths.input_data_file_path,
            output_directory=run_paths.output_directory,
            base_config_file=project_paths.base_config_path,
            plot_fields=None,
            runid=runid,
            verbose=2,
            output_options=output_options,
        ),
        **kwargs,
    }
    with warnings.catch_warnings():
        out = main.single(**kwargs_all)
    return out


project_dir = "tests/key_processes/hourly_sparse_runs"

setups = [
    # ["hourly_sparse_simple", "hourly", "hourly", {}],
    ["bihourly_sparse_simple", "bihourly", "bihourly", dict()],
    ["bihourly_sparse_simple_cropped", "bihourly", "bihourly",
        dict(overrides=['end_day=212', 'dump_state_n=3'])],
]


@pytest.fixture(scope="class")
def hourly_sparse_runs_test_run(request):
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


@pytest.mark.usefixtures('hourly_sparse_runs_test_run')
class TestRunAndCompare:

    def test_preruns_run_without_error(self):
        pass

    @pytest.mark.parametrize('runid', ['bihourly_sparse_simple'])
    def test_should_calculate_dd_correctly(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        dd = hourly_output['dd'].values
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        assert dd[0] is not None
        assert dd[-1] is not None
        assert max(dd) > 0
        assert min(dd) < 365
        assert final_state.temporal.dd < 365

    @pytest.mark.parametrize('runid', ['bihourly_sparse_simple'])
    def test_should_use_multi_hour_hr_correctly(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        hr = hourly_output['hr'].values
        assert hr[0] is not None
        assert hr[-1] is not None

        assert hr[0] == hr[1]

    @pytest.mark.parametrize('runid', ['bihourly_sparse_simple'])
    @pytest.mark.skip(reason="not working yet")
    def test_should_calculate_f_phen(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        f_phen = hourly_output['f_phen'].values
        assert max(f_phen) == 1.0
        assert min(f_phen) == 0.0

    @pytest.mark.parametrize('runid', ['bihourly_sparse_simple'])
    def test_should_reset_SAI(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        sai = hourly_output['canopy_sai'].values
        lai = hourly_output['canopy_lai'].values
        assert all([a == b for a, b in zip(sai, lai)])


"""Test running a few lines of data to make sure multiplicative runs work."""
from pathlib import Path
import pytest
import warnings
import pandas as pd

from pyDO3SE.Output.OutputConfig import (
    OutputOptions
)
from pyDO3SE import main


def run_with_config(runid: str, project_dir: Path, config_file: str, input_file: str, **kwargs):
    project_paths = main.get_project_paths(project_dir)
    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)

    # Create output dir
    main.create_run_path_directories(run_paths)

    output_options = OutputOptions()

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


project_dir = Path("tests/key_processes/multiplicative")

setups = [
    # runid, config_file, input_file, overrides
    ["default_sparse_simple", "default", "three_days", dict()],
    ["bihourly_sparse_simple", "bihourly", "bihourly", dict()],
    ["alt", "alt", "A_03", dict()],
]


@pytest.fixture(scope="class")
def legacy_fphen_test_run(request):
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


def test_run_alt():
    out = run_with_config(
        runid="alt",
        project_dir=project_dir,
        config_file="alt",
        input_file="A_03",
    )
    final_state, output_logs, final_config, initial_state, external_state = out

@pytest.mark.usefixtures('legacy_fphen_test_run')
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
    def test_should_calculate_fphen_correctly(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        f_phen = hourly_output['f_phen'].values
        assert f_phen[0] is not None
        assert f_phen[-1] is not None
        assert all(f is not None for f in f_phen)
        assert max(f_phen) == 1
        assert min(f_phen) == 0
        assert f_phen[-1] == 1

    @pytest.mark.parametrize('runid', ['alt'])
    def test_should_calculate_leaf_f_phen_correctly(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        leaf_f_phen = hourly_output['leaf_f_phen'].values
        assert leaf_f_phen[0] is not None
        assert leaf_f_phen[-1] is not None
        assert all(f is not None for f in leaf_f_phen)
        assert max(leaf_f_phen) == 1
        assert min(leaf_f_phen) == 0
        assert leaf_f_phen[-1] == 0
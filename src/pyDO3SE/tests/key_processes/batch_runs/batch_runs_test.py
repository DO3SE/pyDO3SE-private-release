"""Test running a few lines of data to make sure soil moisture is correct."""
from collections import namedtuple
import pytest
import warnings
from typing import List, Tuple
import math

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main

project_dir = "tests/key_processes/batch_runs"


def run_with_config(runid: str, **kwargs):
    output_options = output_results_only_options()
    output_options.save_hourly_output_data = False
    output_options.save_processed_config = True
    output_options.save_model_processes_detailed = True
    logger = main.Logger(2, log_to_file=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = main.batch(
            project_directory=project_dir,
            runid=runid,
            verbose=2,
            output_options=output_options,
            logger=logger,
            **kwargs,
        )
    return out


RunSetup = namedtuple('RunSetup', ['runid', 'overrides'])

setups = [
    RunSetup("parallel", dict(parallel=True)),
    RunSetup("series", dict(parallel=False)),
    RunSetup("per_run_config", dict(parallel=False)),
]

setupsMap = {s.runid: s for s in setups}
all_setups = [s.runid for s in setups]


@pytest.fixture(scope="class")
def batch_test_run(request):
    request.cls.output = {}

    for runid, overrides in setups:
        results_info: List[Tuple[main.Args, main.RunOutput]] = run_with_config(
            runid=runid,
            **overrides,
        )

        request.cls.output[runid] = {}
        request.cls.output[runid]['out'] = results_info


@pytest.mark.parametrize(['runid', 'overrides'], setups)
def test_runs_without_error(runid, overrides):
    run_with_config(
        runid=runid,
        **overrides,
    )


@pytest.mark.usefixtures('batch_test_run')
class TestRunAndCompare:

    @pytest.mark.parametrize('runid', all_setups)
    def test_preruns_run_without_error(self, runid):
        results_info: List[Tuple[main.Args, main.RunOutput]] = self.output[runid]['out']
        assert all([result == 0 for args, [result, error, model_output] in results_info])

    @pytest.mark.parametrize('runid', ['per_run_config'])
    def test_can_include_per_input_file_config(self, runid):
        results_info: List[Tuple[main.Args, main.RunOutput]] = self.output[runid]['out']
        elevations: List = [model_output.config_processed.Location.elev for args, [
            result, error, model_output] in results_info]
        base_config_elevation = 5.0
        print(elevations)
        assert all(not math.isclose(e, base_config_elevation, abs_tol=1e-2) for e in elevations)

    @pytest.mark.parametrize('runid', ['per_run_config'])
    def test_can_set_per_input_sowing_day(self, runid):
        results_info: List[Tuple[main.Args, main.RunOutput]] = self.output[runid]['out']
        sowing_days: List = [model_output.config_processed.Land_Cover.parameters[0].phenology.key_dates.sowing for args, [
            result, error, model_output] in results_info]
        base_config_sowing = 1
        print(sowing_days)
        assert all(not math.isclose(e, base_config_sowing, abs_tol=1e-2) for e in sowing_days)
        assert sowing_days[0] != sowing_days[1]

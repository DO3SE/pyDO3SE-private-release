"""Test running a few lines of data to make sure soil moisture is correct."""

from collections import namedtuple
from pathlib import Path
import pytest
import warnings
import os
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.in_memory import main, MainOutput
from pyDO3SE.main import get_project_paths, get_run_paths
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.External_State.external_state_loader import (
    load_external_state,
)

project_dir = os.path.dirname(os.path.realpath(__file__))


def run_with_config(
    runid: str,
    project_dir: Path,
    config_file: str,
    input_file: str,
):
    project_paths = get_project_paths(project_dir)
    run_paths = get_run_paths(runid, project_paths, config_file, input_file)

    # Create output dir
    config = config_loader(
        run_paths.config_path, project_paths.base_config_path, "json"
    )

    external_state_data = next(
        load_external_state(
            run_paths.input_data_file_path,
        )
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = main(
            config=config,
            state=Model_State_Shape(),
            external_state=external_state_data,
        )
    return out


RunSetup = namedtuple("RunSetup", ["runid", "config_file", "input_file"])

setups = [
    # RunSetup("default_three_days", "default", "three_days"),
    RunSetup("default_full_season", "default", "full_season"),
]

setupsMap = {s.runid: s for s in setups}
all_setups = [s.runid for s in setups]


@pytest.fixture(scope="class")
def in_memory_test_run(request):
    request.cls.output = {}

    for runid, config_file, input_file in setups:
        out = run_with_config(
            runid=runid,
            project_dir=project_dir,
            config_file=config_file,
            input_file=input_file,
        )

        request.cls.output[runid] = {}
        request.cls.output[runid]["out"] = out


@pytest.mark.usefixtures("in_memory_test_run")
class TestRunAndCompare:
    def test_preruns_run_without_error(self):
        pass


    @pytest.mark.parametrize("runid", all_setups)
    def test_should_output_canopy_sunlit_An(self, runid):
        model_output = self.output[runid]["out"]
        assert isinstance(model_output, MainOutput)

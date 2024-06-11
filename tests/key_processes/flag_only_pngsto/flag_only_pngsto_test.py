"""Test running a few lines of data to make sure ozone deposition is correct."""
from pathlib import Path
from dataclasses import replace
from collections import namedtuple
import pytest
import warnings
import pandas as pd
from math import isclose

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE import main
from pyDO3SE import setup_model
from pyDO3SE.Config.ConfigLandCover import LAIMethods


def run_with_config(
    runid: str,
    project_dir: Path,
    config_file: str,
    input_file: str,
    input_state_file: str = None,
    **kwargs,
):
    print(kwargs)
    project_paths = main.get_project_paths(project_dir)
    if input_state_file is not None:
        project_paths = project_paths._replace(
            base_state_path=f"{project_dir}/initial_state/{input_state_file}.json")

    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)

    # Create output dir
    main.create_run_path_directories(run_paths)

    output_options = output_results_only_options()
    # output_options = OutputOptions()
    output_options.save_hourly_output_data = False
    # output_options.save_processed_config = True

    kwargs_all = {
        **dict(
            config_file=run_paths.config_path,
            data_file=run_paths.input_data_file_path,
            output_directory=run_paths.output_directory,
            base_config_file=project_paths.base_config_path,
            initial_state_path=project_paths.base_state_path,
            plot_fields=None,
            runid=runid,
            verbose=2,
            output_options=output_options,
            overrides=[
                'run_validation=True'
            ]
        ),
        **kwargs,
    }
    print(kwargs_all)
    with warnings.catch_warnings():
        out = main.single(**kwargs_all)
    return out


project_dir = "tests/key_processes/flag_only_pngsto"
RunSetup = namedtuple('RunSetup', ['runid', 'config_file',
                      'input_file', 'initial_state', 'overrides'])

setups = [
    RunSetup("three_days_default", "three_days", "three_days", None, dict()),
    RunSetup("three_days_default_dvi_1_94", "three_days",
             "three_days", "dump_state_dvi_1_94", dict()),
    # RunSetup("three_days_no_ozone", "three_days_no_ozone", "three_days", None, dict()),
    # RunSetup("test1", "lisa_base_config_pngsto", "sparse", None, dict(base_config_file='')),
]

all_setups = [s.runid for s in setups]


@pytest.fixture(scope="class")
def flag_only_pngsto_test_run(request):
    request.cls.output = {}

    for runid, config_file, input_file, initial_state_path, overrides in setups:
        print(overrides)
        out = run_with_config(
            runid=runid,
            project_dir=project_dir,
            config_file=config_file,
            input_file=input_file,
            input_state_file=initial_state_path,
            **overrides,
        )

        final_state, output_logs, final_config, initial_state, external_state = out
        request.cls.output[runid] = {}
        request.cls.output[runid]['out'] = out
        request.cls.output[runid]['hourly_output'] = pd.DataFrame(output_logs)


@pytest.mark.usefixtures('flag_only_pngsto_test_run')
class TestRunAndCompare:

    def test_preruns_run_without_error(self):
        pass

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_gsto(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        gsto = hourly_output['gsto'].values
        config: Config_Shape = self.output[runid]['out'][2]
        min_gsto = config.Land_Cover.parameters[0].pn_gsto.g_sto_0
        assert gsto[0] is not None
        assert gsto[-1] is not None
        assert max(gsto) > 0
        assert isclose(min(gsto), min_gsto / 1000, abs_tol=1)

    @pytest.mark.parametrize('runid', ['three_days_default_dvi_1_94'])
    def test_should_calculate_pody(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        pody = hourly_output['pody'].values
        assert pody[0] is not None
        assert pody[-1] is not None
        assert max(pody) > 0

    @pytest.mark.parametrize('runid', ['three_days_default_dvi_1_94'])
    def test_should_calculate_sunlit_pody(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        # leaf_rmodel_O3
        # leaf_rmodel_O3
        # leaf_rmodel_O3
        # mean_gsto_per_layer
        # fst_sunlit
        fst_sun = hourly_output['fst_sun'].values
        print(fst_sun)
        pody_sun = hourly_output['pody_sun'].values
        print(pody_sun)
        assert pody_sun[0] is not None
        assert pody_sun[-1] is not None
        assert max(pody_sun) > 0


def test_should_not_change_pody_if_lai_changes():
    config = main.config_loader("tests/key_processes/flag_only_pngsto/base_config.json")
    external_state_data = next(setup_model.load_external_state(
        "tests/key_processes/flag_only_pngsto/inputs/three_days.csv",
    ))

    processed_config = setup_model.setup_config(
        config,
        external_state_data,
    )

    # TODO: make LAI constant
    processed_config.Land_Cover.LAI_method = LAIMethods.CONSTANT
    processed_config.Land_Cover.LAI = 3
    initial_state = main.model_state_loader(
        "tests/key_processes/flag_only_pngsto/initial_state/dump_state_dvi_1_94.json")
    initial_state = replace(initial_state, prev_hour=initial_state)
    external_state, start_day, end_day = setup_model.setup_external_state(
        processed_config, external_state_data)
    model_processes = setup_model.setup_model_processes_hour_map(processed_config)

    final_state, output_logs = main.run_model_on_mapped_processes(
        initial_state,
        processed_config,
        external_state,
        model_processes,
    )

    assert isclose(final_state.canopy_component_population[0][0].POD_Y, 0.227126866, abs_tol=1e-4)
    assert final_state.canopy.LAI_total == 3
    assert final_state.canopy_component[0].LAI == 3
    assert output_logs[0]['pody'] < output_logs[-1]['pody']

    processed_config.Land_Cover.LAI_method = LAIMethods.CONSTANT
    processed_config.Land_Cover.LAI = 6
    final_state_b, output_logs = main.run_model_on_mapped_processes(
        initial_state,
        processed_config,
        external_state,
        model_processes,
    )
    assert isclose(final_state_b.canopy_component_population[0][0].POD_Y, 0.227126866, abs_tol=1e-4)
    assert final_state_b.canopy_component[0].LAI == 6

    assert output_logs[0]['pody'] < output_logs[-1]['pody']

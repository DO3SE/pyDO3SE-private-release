"""Test running a few lines of data to make sure soil moisture is correct."""
from collections import namedtuple
from math import isclose
from pathlib import Path
from copy import deepcopy
import pytest
import warnings
import pandas as pd
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Model_State.Model_State import Model_State_Shape

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main

project_dir = "tests/key_processes/carbon_allocation"


def run_with_config(
    runid: str,
    project_dir: Path,
    config_file: str,
    input_file: str,
    input_state_file: str,
    **kwargs,
):
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
    RunSetup("default_dvi_1_49", "default", "three_days", 'dump_state_dvi_1_49', {}),
    RunSetup("default_dvi_1_94", "default", "three_days", 'dump_state_dvi_1_94', {}),
]

setupsMap = {s.runid: s for s in setups}
all_setups = [s.runid for s in setups]


@pytest.fixture(scope="class")
def carbon_allocation_test_run(request):
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


@pytest.mark.usefixtures('carbon_allocation_test_run')
class TestRunAndCompare:

    def test_preruns_run_without_error(self):
        pass

    @pytest.mark.parametrize('runid', ['default_dvi_1_03'])
    def test_should_preserve_mass(self, runid):
        """Total NPP should equal the sum of all carbon pools."""
        hourly_output = self.output[runid]['hourly_output']
        c_leaf = hourly_output['c_leaf'].values
        c_harv = hourly_output['c_harv'].values
        c_stem = hourly_output['c_stem'].values
        c_lbrn = hourly_output['c_lbrn'].values
        c_root = hourly_output['c_root'].values
        c_resv = hourly_output['c_resv'].values
        npp = hourly_output['npp'].values
        npp_acc = hourly_output['npp_acc'].values
        carbon_pool_start = sum([
            c_leaf[0],
            c_harv[0],
            c_stem[0],
            c_lbrn[0],
            c_root[0],
            c_resv[0],
        ])
        carbon_pool_end = sum([
            c_leaf[-1],
            c_harv[-1],
            c_stem[-1],
            c_lbrn[-1],
            c_root[-1],
            c_resv[-1],
        ])
        carbon_pool_change = carbon_pool_end - carbon_pool_start
        total_carbon_added = sum(npp) - npp_acc[-1]
        assert isclose(carbon_pool_change, total_carbon_added, abs_tol=1e-6)

    @pytest.mark.parametrize('runid', ['default_dvi_1_03'])
    def test_should_increase_carbon_pools(self, runid):
        """Total NPP should equal the sum of all carbon pools."""
        hourly_output = self.output[runid]['hourly_output']
        c_leaf = hourly_output['c_leaf'].values
        c_harv = hourly_output['c_harv'].values
        c_stem = hourly_output['c_stem'].values
        c_lbrn = hourly_output['c_lbrn'].values
        c_root = hourly_output['c_root'].values
        c_resv = hourly_output['c_resv'].values

        assert max(c_leaf) > 0
        assert max(c_harv) > 0
        assert max(c_stem) > 0
        assert max(c_lbrn) > 0
        assert max(c_root) > 0
        assert max(c_resv) > 0

    @pytest.mark.parametrize('runid', ['default_dvi_1_94'])
    def test_should_distribute_green_leaf_carbon_to_brown_leaf_carbon_after_dvi_1_5(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        hourly_output = self.output[runid]['hourly_output']
        c_leaf = hourly_output['c_leaf'].values
        c_lbrn = hourly_output['c_lbrn'].values
        assert c_leaf[-1] < c_leaf[0]
        assert c_lbrn[-1] > c_lbrn[0]

    @pytest.mark.parametrize('runid', ['default_dvi_1_94'])
    def test_should_distribute_carbon_to_harvest_after_dvi_1_5(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        hourly_output = self.output[runid]['hourly_output']
        c_harv = hourly_output['c_harv'].values
        c_resv = hourly_output['c_resv'].values
        assert c_resv[-1] < c_resv[0]
        assert c_harv[-1] > c_harv[0]

    @pytest.mark.parametrize('runid', ['default_dvi_1_94'])
    def test_should_have_increased_brown_lai_after_start_of_senescence(self, runid):
        # - split leaf carbon between green and brown leaf. 85% of green goes to brown
        # - split decrease 85/15 grain brown leaf

        hourly_output = self.output[runid]['hourly_output']
        canopy_lai_brown = hourly_output['canopy_lai_brown'].values
        assert canopy_lai_brown[0] < canopy_lai_brown[-1]

    @pytest.mark.parametrize('runid', ['default_dvi_1_94'])
    def test_should_set_sai_to_green_plus_brown_lai(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        canopy_sai = hourly_output['canopy_sai'].values
        canopy_lai = hourly_output['canopy_lai'].values
        canopy_lai_brown = hourly_output['canopy_lai_brown'].values
        assert canopy_sai[-1] == canopy_lai[-1] + canopy_lai_brown[-1]

    @pytest.mark.parametrize('runid', ['default_dvi_1_94'])
    def test_should_distribute_sai_between_layers(self, runid):
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        hourly_output = self.output[runid]['hourly_output']
        canopy_sai = hourly_output['canopy_sai'].values
        final_state: Model_State_Shape = final_state
        config: Config_Shape = final_config
        nL = config.Land_Cover.nL
        total_sai = 0
        for iL in range(nL):
            total_sai += final_state.canopy_layer_component[iL][0].SAI
            assert final_state.canopy_layer_component[iL][0].SAI > 0
        assert total_sai > 0
        assert canopy_sai[-1] > 0
        assert final_state.canopy.SAI_total > 0

        assert total_sai == final_state.canopy.SAI_total
        assert final_state.canopy.SAI_total == canopy_sai[-1]

        # TODO: Check why we this is 0
        # assert final_state.canopy_component[0].SAI > 0

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_yield_outputs(self, runid):
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        final_state: Model_State_Shape = final_state
        stem_dm = final_state.canopy_component[0].stem_dm
        assert stem_dm > 0, "stem_dm has not increased above 0!"
        leaf_dm = final_state.canopy_component[0].leaf_dm
        assert leaf_dm > 0, "leaf_dm has not increased above 0!"
        lbrn_dm = final_state.canopy_component[0].lbrn_dm
        assert lbrn_dm > 0, "lbrn_dm has not increased above 0!"
        total_leaf_dm = final_state.canopy_component[0].total_leaf_dm
        assert total_leaf_dm > 0, "total_leaf_dm has not increased above 0!"
        straw_dm = final_state.canopy_component[0].straw_dm
        assert straw_dm > 0, "straw_dm has not increased above 0!"
        ear_dm = final_state.canopy_component[0].ear_dm
        assert ear_dm > 0, "ear_dm has not increased above 0!"
        aboveground_dm = final_state.canopy_component[0].aboveground_dm
        assert aboveground_dm > 0, "aboveground_dm has not increased above 0!"
        belowground_dm = final_state.canopy_component[0].belowground_dm
        assert belowground_dm > 0, "belowground_dm has not increased above 0!"
        grain_dm = final_state.canopy_component[0].grain_dm
        assert grain_dm > 0, "grain_dm has not increased above 0!"
        harvest_index = final_state.canopy_component[0].harvest_index
        assert harvest_index > 0, "harvest_index has not increased above 0!"
        yield_ha = final_state.canopy_component[0].yield_ha
        assert yield_ha > 0, "yield_ha has not increased above 0!"

    @pytest.mark.parametrize('runid', ['default_dvi_1_49'])
    def test_should_output_respiration_and_npp(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        A_n_canopy = hourly_output['A_n_canopy'].values
        gpp = hourly_output['gpp'].values
        npp = hourly_output['npp'].values
        R_pg = hourly_output['R_pg'].values
        R_pm = hourly_output['R_pm'].values
        assert any(i > 0 for i in A_n_canopy)
        assert any(i > 0 for i in gpp)
        assert any(i > 0 for i in npp)
        assert any(i > 0 for i in R_pg)
        assert any(i > 0 for i in R_pm)


def run_single_override_config(runid, config_override=None):
    runid, config_file, input_file, input_state, overrides = setupsMap[runid]
    project_paths = main.get_project_paths(project_dir)

    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)
    project_paths = project_paths._replace(
        base_state_path=f"{project_dir}/initial_state/{input_state}.json")

    loaded_run_files = main.load_run_files(
        project_paths=project_paths,
        run_paths=run_paths,
    )

    overrides_main = main.Main_Overrides(**overrides)
    overrides_main = overrides_main._replace(skip_state_init=True)

    [
        config,
        external_state,
        initial_state,
        model_processes,
    ] = main.setup_model(
        config_in=loaded_run_files.config,
        state_in=loaded_run_files.state,
        data_location=run_paths.input_data_file_path,
        run_dir=run_paths.run_dir,
        # logger=logger,
        overrides=overrides_main,
    )

    config_in = config if not config_override else config_override(deepcopy(config))

    final_state, output_logs = main.run_model_on_mapped_processes(
        initial_state,
        config_in,
        external_state,
        model_processes,
        DEBUG_MODE=True,
    )
    return final_state, config_in


class TestIndividual:
    @pytest.mark.parametrize('runid', ['default_dvi_1_94'])
    def test_setting_carbon_leaf_fractions(self, runid):
        final_state, config = run_single_override_config(runid)

        assert final_state.canopy_component[0].c_leaf > 0, "c_leaf is 0"
        assert final_state.canopy_component[0].c_lbrn > 0, "c_lbrn is 0"
        assert final_state.canopy_component[0].c_harv > 0, "c_harv is 0"
        assert final_state.canopy_component[0].c_resv > 0, "c_resv is 0"
        assert final_state.canopy_component[0].c_stem > 0, "c_stem is 0"
        assert final_state.canopy_component[0].c_root > 0, "c_root is 0"

        # # Decreasing f_brown_leaf increases leaf carbon and reduces harvest and brown leaf carbon
        def config_override(c):
            c.carbon_allocation.f_green_leaf = 0.5
            return c

        final_state_b, config_b = run_single_override_config(runid, config_override)

        assert final_state.canopy_component[0].c_leaf > final_state_b.canopy_component[0].c_leaf
        assert final_state.canopy_component[0].c_lbrn < final_state_b.canopy_component[0].c_lbrn
        assert final_state.canopy_component[0].c_harv < final_state_b.canopy_component[0].c_harv

        # Decreasing f_green_leaf increases brown leaf carbon and reduces harvest
        def config_override(c):
            c.carbon_allocation.f_brown_leaf = 0.5
            return c

        final_state_c, config_c = run_single_override_config(runid, config_override)

        assert final_state.canopy_component[0].c_leaf == final_state_c.canopy_component[0].c_leaf
        assert final_state.canopy_component[0].c_lbrn > final_state_c.canopy_component[0].c_lbrn
        assert final_state.canopy_component[0].c_harv < final_state_c.canopy_component[0].c_harv

    @pytest.mark.parametrize('runid', ['default_dvi_1_94'])
    def test_should_increase_grain_dm_if_we_increase_grain_to_ear(self, runid):
        final_state, config = run_single_override_config(runid)

        assert final_state.canopy_component[0].c_leaf > 0, "c_leaf is 0"
        assert final_state.canopy_component[0].c_lbrn > 0, "c_lbrn is 0"
        assert final_state.canopy_component[0].c_harv > 0, "c_harv is 0"
        assert final_state.canopy_component[0].c_resv > 0, "c_resv is 0"
        assert final_state.canopy_component[0].c_stem > 0, "c_stem is 0"
        assert final_state.canopy_component[0].c_root > 0, "c_root is 0"

        # increasing grain_to_ear should increase grain_dm
        def config_override(c):
            c.carbon_allocation.grain_to_ear = c.carbon_allocation.grain_to_ear * 1.1
            return c

        final_state_b, config_b = run_single_override_config(runid, config_override)

        assert final_state.canopy_component[0].grain_dm < final_state_b.canopy_component[0].grain_dm

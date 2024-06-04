"""A set of tests that run the full model then compare the output against the previous version."""
from tests.test_helpers import long_test
import warnings
import pandas as pd
import os
import pytest
from shutil import copyfile

from do3se_phenology.state import LeafPhenologyStage

from pyDO3SE.Output.run_comparisons import (
    create_comparison_graphs,
    get_output_files,
    multip_output_fields,
    pn_output_fields,
)
from pyDO3SE.Output.Output_Shape import default_output_fields
from pyDO3SE.version import version
from pyDO3SE import main


example_folders = [
    # './examples/spanish_wheat_multiplicative', TODO: Fix this
    # TODO: Change this back
    # './examples/spanish_wheat',
    './examples/bangor_2015',
    # './examples/bangor_2015_multiplicative',
    './examples/bangor_2015_multilayer',
    # './examples/bangor_2015_multilayer_simple',
    './examples/icp_veg',
    './examples/xiaoji_2008',
    # './examples/icp_veg_bangor',
]
setups = [
    ','.join([project_dir, cf.split('.json')[0], f.split('.csv')[0]])
    for project_dir in example_folders
    for cf in sorted(os.listdir(f"{project_dir}/configs"))
    for f in os.listdir(f"{project_dir}/inputs")
]


@pytest.fixture(scope="class", params=setups)
def setup_test_run_and_compare(request):
    project_dir, config_file, input_file = request.param.split(',')
    skip_outputs = not request.config.getoption("--skip-outputs") == "False"
    request.cls.project_dir = project_dir
    errors = []
    try:
        runid = "demo_test_output_comparison"
        project_paths = main.get_project_paths(project_dir)
        run_paths = main.get_run_paths(project_paths, config_file, input_file, runid)

        # Create output dir
        main.create_run_path_directories(run_paths)

        demo_config = run_paths.config_path
        base_demo_config = project_paths.base_config_path
        demo_output_dir = run_paths.output_directory

        # Copy config to this directory
        copyfile(demo_config, f'{demo_output_dir}/config.json')
        output_fields = default_output_fields  # pn_output_fields # TODO: Fix for multip
        fields_to_graph = ['gsto_canopy', 'td_dd', 'lai', 'pody', 'fst', 'canopy_height']
        if skip_outputs:
            print("skipping outputs")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = main.main(
                project_paths=project_paths,
                run_paths=run_paths,
                runnotes="Running from test",
                fields_to_graph=fields_to_graph,
                use_daily_loop=True,
                output_fields=output_fields,
            )

        final_state, output_logs, final_config, initial_state = out
        request.cls.output = out
        request.cls.logs = pd.DataFrame(output_logs)
        # Make a copy for the version
        # NOTE: Disabled as creating too much data!
        # if not skip_outputs:
        #     versioned_dir = f"{project_dir}/outputs/{input_file}/{config_file}/{version}"
        #     os.makedirs(versioned_dir, exist_ok=True)
        #     copytree(demo_output_dir, versioned_dir, dirs_exist_ok=True)
    except Exception as e:
        errors.append((f"Config: {config_file} with input {input_file} failed", e))

    if len(errors) > 0:

        for m, e in errors:
            print(m)
        print(errors)
        raise errors[0][1]
    run_type = final_config.Land_Cover.parameters[0].gsto.method
    start_day = final_config.Location.start_day or 0
    end_day = final_config.Location.end_day or 365
    request.cls.run_type = run_type

    if not skip_outputs:
        # Compare outputs
        if run_type == "multiplicative":
            fields_to_graph = multip_output_fields
        elif run_type == "photosynthesis":
            fields_to_graph = pn_output_fields
        else:
            raise Exception("Invalid run type")
        output_dir = f"{project_dir}/outputs"
        data_files = list(sorted(get_output_files(output_dir, "latest"), key=lambda f: f[1]))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            create_comparison_graphs(
                data_files,
                output_dir,
                fields_to_graph,
                f"{project_dir}/comparisons",
                use_versioned_outfile=False,
                use_versioned_outdir=False,
                log_level=0,
                start_day=start_day,
                end_day=end_day,
            )


@pytest.mark.usefixtures('setup_test_run_and_compare')
class TestRunAndCompare:

    @long_test
    def test_main(self):
        pass

    @long_test
    def test_should_have_increased_pody(self):
        final_state, output_logs, final_config, initial_state = self.output
        assert final_state.canopy_component_population[0][-1].POD_Y > 0

    @long_test
    def test_should_have_increased_td_dd(self):
        final_state, output_logs, final_config, initial_state = self.output
        assert final_state.canopy_component[0].td_dd > 0
        assert final_state.canopy_component[0].td_dd_leaf_pops[-1] > 0

    @long_test
    def test_flag_leaf_should_have_finished_growing(self):
        final_state, output_logs, final_config, initial_state = self.output
        assert final_state.canopy_component_population[0][-1].phenology.phenology_stage > LeafPhenologyStage.GROWING

    @long_test
    def test_should_have_set_f_LS_to_zero(self):
        if self.run_type == "photosynthesis":
            final_state, output_logs, final_config, initial_state = self.output
            # Note some runs don't end exactly at 0
            assert final_state.canopy_component_population[0][-1].f_LS < 0.025

    @long_test
    def test_should_have_increased_gsto(self):
        output_logs = self.logs
        assert max(output_logs['gsto']) > 0

    @long_test
    def test_flag_leaf_should_have_emerged(self):
        final_state, output_logs, final_config, initial_state = self.output
        assert final_state.canopy_component[0].total_emerged_leaf_populations > 0

    @long_test
    def test_lai_distributed_correctly(self):
        final_state, output_logs, final_config, initial_state = self.output
        np = final_config.Land_Cover.nP == 0
        if np == 0:
            return


        # TODO: Assert that the sum of layer LAI == sum of population LAI == total canopy lai
        # NOTE This is difficult with current setup as when canopy LAI goes down at end of season the
        # layer and pop LAI does not.
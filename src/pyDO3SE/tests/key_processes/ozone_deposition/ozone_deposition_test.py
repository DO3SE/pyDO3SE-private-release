"""Test running a few lines of data to make sure ozone deposition is correct."""
from math import isclose
from pathlib import Path
import numpy as np
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
    output_options.save_processed_config = True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = main.single(
            config_file=run_paths.config_path,
            data_file=run_paths.input_data_file_path,
            output_directory=run_paths.output_directory,
            base_config_file=project_paths.base_config_path,
            plot_fields=None,
            runid=runid,
            verbose=0,
            output_options=output_options,
        )
    return out


project_dir = "tests/key_processes/ozone_deposition"

setups = [
    ["legacy_three_days", "legacy", "three_days", {}],
    ["legacy_two_hours", "legacy", "two_hours", dict()],
    ["canopy_measure_height_two_hours", "canopy_measure_height",
        "two_hours", dict()],
    ["target_h_equals_model_h_two_hours", "target_h_equals_model_h",
        "two_hours", dict()],
    ["measured_h_eq_izr_two_hours", "measured_h_eq_izr",
        "two_hours", dict()],
    ["input_ustar_ref_two_hours", "input_ustar_ref", "two_hours", dict()],
    ["input_ustar_two_hours", "input_ustar", "two_hours", dict()],
    ["otc_three_days", "using_otc", "three_days", {}],
]

all_setups = [s[0] for s in setups]


def otherThan(ignore_these):
    return [s for s in all_setups if s not in ignore_these]


@pytest.fixture(scope="class")
def ozone_deposition_test_run(request):
    request.cls.output = {}

    for runid, config_file, input_file, overrides in setups:
        try:
            out = run_with_config(
                runid=runid,
                project_dir=project_dir,
                config_file=config_file,
                input_file=input_file,
                **overrides,
            )
        except Exception as e:
            print(f"Failed to run {runid}")
            raise e

        final_state, output_logs, final_config, initial_state, external_state = out
        request.cls.output[runid] = {}
        request.cls.output[runid]['out'] = out
        request.cls.output[runid]['hourly_output'] = pd.DataFrame(output_logs)


@pytest.mark.usefixtures('ozone_deposition_test_run')
class TestRunAndCompare:

    def test_preruns_run_without_error(self):
        pass

    @pytest.mark.parametrize('runid', ['legacy_three_days', 'otc_three_days'])
    def test_should_calculate_ustar_correctly(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        ustar = hourly_output['ustar'].values
        ustar_ref = hourly_output['ustar_ref'].values
        u_i = hourly_output['u_i'].values
        assert all(ustar)
        assert all(ustar_ref)
        assert all(u_i)

    @pytest.mark.parametrize('runid', ['canopy_measure_height_two_hours', 'otc_three_days'])
    def test_ozone_in_should_be_the_same_as_ozone_out(self, runid):
        """Test that ozone remains the same if not transfered up or down."""
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        canopy_top_o3 = hourly_output['canopy_top_o3'].values
        micro_o3 = hourly_output['micro_O3'].values
        # assert final_config.Location.h_O3 == initial_state.canopy.canopy_height
        if not final_config.Location.OTC:
            assert final_config.Location.h_O3 == final_config.Location.z_O3
            assert final_config.Location.h_O3 == final_state.canopy.canopy_height
        external_state: External_State_Shape = external_state
        assert all(canopy_top_o3)
        assert all([isclose(a, b, abs_tol=1e-1) for a, b in zip(canopy_top_o3, external_state.O3)])
        assert all([isclose(a, b, abs_tol=1e-1) for a, b in zip(canopy_top_o3, micro_o3)])

    @pytest.mark.parametrize('runid', ['measured_h_eq_izr_two_hours'])
    def test_o3_i_should_equal_o3_zr_if_h_o3_equals_izr(self, runid):
        # TODO: Is this right
        """Test that ozone changes if transfered up or down with different canopy heights."""
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        assert final_config.Location.z_O3 == final_config.Location.izr

        o3_ppb_i = hourly_output['o3_ppb_i'].values
        o3_ppb_zr = hourly_output['o3_ppb_zr'].values
        print(o3_ppb_i)
        print(o3_ppb_zr)
        external_state: External_State_Shape = external_state
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(o3_ppb_i, external_state.O3)])
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(o3_ppb_zr, external_state.O3)])

    @pytest.mark.parametrize('runid', otherThan(['measured_h_eq_izr_two_hours', 'otc_three_days']))
    def test_ozone_transfered_up_should_increase_ozone_at_izr(self, runid):
        # TODO: Is this right
        """Test that ozone changes if transfered up or down with different canopy heights."""
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        assert final_config.Location.z_O3 < final_config.Location.izr
        print(final_config.Location.z_O3, final_config.Location.izr)
        o3_ppb_i = hourly_output['o3_ppb_i'].values

        external_state: External_State_Shape = external_state
        assert any([a > b for a, b in zip(o3_ppb_i, external_state.O3)])

    @pytest.mark.parametrize('runid', otherThan(['canopy_measure_height_two_hours', 'target_h_equals_model_h_two_hours', 'otc_three_days', 'input_ustar_ref_two_hours']))
    def test_ozone_in_should_not_be_the_same_as_ozone_out_when_transfered(self, runid):
        # TODO: Is this right
        """Test that ozone changes if transfered up or down with different canopy heights."""
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        micro_O3 = hourly_output['micro_O3'].values
        external_state: External_State_Shape = external_state
        assert all(micro_O3)
        assert any([not isclose(a, b, abs_tol=1e-1) for a, b in zip(micro_O3, external_state.O3)])

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_ra(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        ra = hourly_output['ra'].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in ra))

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_rb(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        rb = hourly_output['rb'].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rb))

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_rsur(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        rsur = hourly_output['rsur'].values
        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rsur))

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_rinc(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        rinc = hourly_output['rinc'].values
        canopy_sai = hourly_output['canopy_sai'].values

        print(rinc)
        print(canopy_sai)
        external_state: External_State_Shape = external_state
        assert all((sai is not None and sai > 0 for sai in canopy_sai))
        assert all((r is not None and r > 0 for r in rinc))

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_rext(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        rext = hourly_output['rext'].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rext))

    @pytest.mark.parametrize('runid', all_setups)
    def test_should_calculate_rsto_c(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        rsto_c = hourly_output['rsto_c'].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rsto_c))

    @pytest.mark.parametrize('runid', ['canopy_measure_height_two_hours', 'input_ustar_ref_two_hours', 'otc_three_days'])
    def test_windspeed_in_should_be_the_same_as_windspeed_out(self, runid):
        """Test that windspeed remains the same if not transfered up or down."""
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        micro_u = hourly_output['micro_u'].values

        external_state: External_State_Shape = external_state
        assert all(micro_u)
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(micro_u, external_state.u)])

    @pytest.mark.parametrize('runid', ['canopy_measure_height_two_hours', 'input_ustar_ref_two_hours'])
    def test_ustar_ref_in_should_be_the_same_as_ustar_out_if_using_observed_canopy_height(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        ustar = hourly_output['ustar'].values

        external_state: External_State_Shape = external_state
        ustar_observed = external_state.ustar_ref
        assert all(ustar_observed)
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(ustar, ustar_observed)])

    @pytest.mark.parametrize('runid', ['input_ustar_two_hours'])
    def test_ustar_in_should_be_the_same_as_ustar_out_if_using_observed_canopy_height(self, runid):
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        ustar = hourly_output['ustar'].values

        external_state: External_State_Shape = external_state
        ustar_observed = external_state.ustar
        assert all(ustar_observed)
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(ustar, ustar_observed)])

    @pytest.mark.parametrize('runid', ['measured_h_eq_izr_two_hours'])
    def test_that_if_measured_height_equal_to_izr_that_O3_i_equals_O3_zr(self, runid):
        """Test that when the measured height is equal to the izr(Mixed height) that O3_i(At izr) equals O3_zr(measured O3)."""
        hourly_output = self.output[runid]['hourly_output']
        final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
        o3_ppb_i = hourly_output['o3_ppb_i'].values
        o3_ppb_zr = hourly_output['o3_ppb_zr'].values
        external_state: External_State_Shape = external_state
        O3_observed = external_state.O3
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(o3_ppb_i, o3_ppb_zr)])
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(O3_observed, o3_ppb_zr)])

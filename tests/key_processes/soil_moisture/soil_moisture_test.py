"""Test running a few lines of data to make sure soil moisture is correct."""

from pathlib import Path
import pytest
import warnings
import numpy as np
from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from do3se_met.soil_moisture.helpers import SWP_to_SWC
from do3se_met.soil_moisture.enums import SoilMoistureSource
from pyDO3SE import main
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from tests.utils import get_test_run, Setup

project_dir = "tests/key_processes/soil_moisture"


def run_with_config(runid: str, project_dir: Path, config_file: str, input_file: str, **kwargs):
    project_paths = main.get_project_paths(project_dir)
    run_paths = main.get_run_paths(runid, project_paths, config_file, input_file)

    # Create output dir
    main.create_run_path_directories(run_paths)

    output_options = output_results_only_options()
    output_options.save_hourly_output_data = False
    output_options.save_processed_config = True
    output_options.save_model_processes_detailed = True

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


# TODO: Fix precip tests
setups: list[Setup] = [
    Setup("precip_three_days", "precip", "three_days"),
    Setup("swc_input_three_days", "swc_input", "three_days"),
    Setup("asw_input_three_days", "asw_input", "three_days_asw"),
]

all_setups = [s[0] for s in setups]


def otherThan(ignore_these):
    return [s for s in all_setups if s not in ignore_these]


def get_setup(runid: str):
    return next(s for s in setups if s[0] == runid)


run_outputs = {}


model_run = get_test_run(setups, run_outputs, project_dir)

# @pytest.fixture(scope="class")
# def soil_moisture_test_run(request):
#     request.cls.output = {}

#     for runid, config_file, input_file, overrides in setups:
#         out = run_with_config(
#             runid=runid,
#             project_dir=project_dir,
#             config_file=config_file,
#             input_file=input_file,
#             **overrides,
#         )

#         final_state, output_logs, final_config, initial_state, external_state = out
#         request.cls.output[runid] = {}
#         request.cls.output[runid]['out'] = out
#         request.cls.output[runid]['hourly_output'] = pd.DataFrame(output_logs)


class TestRunAndCompare:
    @pytest.mark.parametrize("runid", all_setups)
    def test_preruns_run_without_error(self, runid):
        model_run(runid)

    @pytest.mark.parametrize("runid", ["precip_three_days"])
    def test_should_accumulate_precip(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        precip_acc = hourly_output["precip_acc"].values
        assert any([p > 0 for p in precip_acc])

    @pytest.mark.parametrize("runid", ["precip_three_days"])
    def test_should_accumulate_penman_monteith_daily_values(self, runid):
        output = model_run(runid)
        final_state, output_logs, final_config, initial_state, external_state = output.out
        assert initial_state.canopy.PM.Sn_diff == 0
        assert final_state.canopy.PM.Sn_diff < 0
        assert final_state.canopy.SMD.Sn > 0

    @pytest.mark.parametrize("runid", ["precip_three_days"])
    def test_should_have_correct_penman_monteith_setup(self, runid):
        """Check that Sn is decreased from initial value.
        For P-M setup this relies on sn_diff which is calculated in penman_monteith_daily.

        """
        output = model_run(runid)
        final_state, output_logs, final_config, initial_state, external_state = output.out
        # config_processed: Config_Shape = self.output[runid]['out'].config_processed
        assert final_config.soil_moisture.source == SoilMoistureSource.P_M
        # NOTE: We set a high SWP_max to force Es_blocked = False
        assert final_config.Land_Cover.parameters[0].gsto.SWP_max == -0.9

    @pytest.mark.parametrize("runid", ["swc_input_three_days", "precip_three_days"])
    def test_should_decrease_Sn(self, runid):
        """Check that Sn is decreased from initial value.

        For P-M setup this relies on sn_diff which is calculated in penman_monteith_daily.

        """

        # TODO: Fix for precip
        output = model_run(runid)
        hourly_output = output.hourly_output
        sn = hourly_output["sn"].values
        final_state, output_logs, final_config, initial_state, external_state = output.out
        # config_processed: Config_Shape = self.output[runid]['out'].config_processed
        sn_max = final_config.soil_moisture.soil_config.FC
        soil_config = final_config.soil_moisture.soil_config
        sn_min = SWP_to_SWC(soil_config.SWP_AE, soil_config.b, final_config.soil_moisture.PWP)
        assert sn_max > 0.2, f"sn_max is too low"

        assert min(sn) < sn_max, f"Sn never goes below max cutoff value: {max(sn)}"

        assert min(sn) < 0.25, f"Sn does not go low enough to impact f_SW min: {min(sn)}"
        assert all([0 < s <= sn_max for s in sn]), "Sn goes outside threshold"

    @pytest.mark.parametrize("runid", ["swc_input_three_days"])
    # @pytest.mark.parametrize('runid', ['precip_three_days', 'swc_input_three_days'])
    def test_swp_output_should_go_down(self, runid):
        """Check that SWP calculated correctly.

        NOTE: swp is directly tied to Sn.

        """
        output = model_run(runid)
        hourly_output = output.hourly_output
        swp = hourly_output["swp"].values
        assert all([s <= 0 for s in swp])
        assert min(swp) < -0.15, f"swp does not drop low enough to impact f_SW, min: {min(swp)}"

    @pytest.mark.parametrize("runid", ["swc_input_three_days"])
    # @pytest.mark.parametrize('runid', ['precip_three_days', 'swc_input_three_days'])
    def test_f_sw_output_should_go_down(self, runid):
        """Here we check it calcs correct f_SW.

        NOTE: f_SW calcs take SWP as input.

        """
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out

        assert final_config.Land_Cover.parameters[0].gsto.f_SW_method != "disabled"
        f_SW = hourly_output["f_SW"].values
        assert min(f_SW) < 1, "f_SW never drops below 1. Check SWP values."
        assert all([0 < f <= 1 for f in f_SW])

    @pytest.mark.parametrize("runid", ["swc_input_three_days", "asw_input_three_days"])
    def test_f_sw_uses_correct_method(self, runid):
        """Check that f_SW uses correct method."""
        output = model_run(runid)
        final_state, output_logs, final_config, initial_state, external_state = output.out
        # Check that the f_SW method is not the default which is disabled
        assert final_config.Land_Cover.parameters[0].gsto.f_SW_method != "disabled"


    @pytest.mark.parametrize("runid", ["asw_input_three_days"])
    def test_asw_input_calculates_f_sw_correctly(self, runid):
        """Check that f_SW is calculated correctly when using asw_input."""
        output = model_run(runid)
        final_state, output_logs, final_config, initial_state, external_state = output.out

        hourly_output = output.hourly_output

        fmin = final_config.Land_Cover.parameters[0].gsto.fmin
        ASW_FC = final_config.soil_moisture.ASW_FC
        ASW_min = final_config.Land_Cover.parameters[0].gsto.ASW_min
        ASW_max = final_config.Land_Cover.parameters[0].gsto.ASW_max
        ASW = np.array(external_state.ASW or [])
        assert ASW_FC == 1.0
        f_PAW_val = fmin + (1 - fmin) * ((100 * ((ASW) / ASW_FC)) - ASW_min) / (ASW_max - ASW_min)
        f_PAW_out_calced = [np.min([1.0, np.max([fmin, v])]) for v in f_PAW_val]

        f_sw_model_out = hourly_output["f_SW"].values
        assert np.allclose(f_PAW_out_calced, f_sw_model_out, atol=0.01)

        assert min(f_sw_model_out) < 1, "f_SW never drops below 1. Check SWP values."

    # @pytest.mark.parametrize('runid', ['precip_three_days'])
    # def test_smd_output_should_go_up(self, runid):
    #     output = model_run(runid)
    #     hourly_output = output.hourly_output
    #     smd = hourly_output['smd'].values
    #     assert max(smd) > 0, f"smd never increases beyond 0"
    #     assert all([s >= 0 for s in smd])

    # TODO: Finish implementing this test
    # @pytest.mark.parametrize('runid', ['precip_three_days'])
    # def test_should_reduce_swp_when_precip_increases(self, runid):
    #     output = model_run(runid)
    #     hourly_output = output.hourly_output
    #     final_state, output_logs, final_config, initial_state, external_state = self.output[runid]['out']
    #     swp = hourly_output['swp'].values
    #     precip_acc = hourly_output['precip_acc'].values
    #     swp_daily = [s for i, s in enumerate(swp) if i % 24 == 0]

    #     external_state: External_State_Shape = external_state
    #     precip = external_state.precip
    #     daily_precip = [sum(precip[i:i + 24]) for i, _ in enumerate(precip[0:-24]) if i % 24 == 0]
    #     daily_precip_grad = [p - pp for p, pp in zip(daily_precip, daily_precip[1:])]
    #     daily_precip_is_increasing = [d > 0 for d in daily_precip_grad]
    #     assert all([s > 0 if d else s <= 0 for s, d in zip(swp_daily, daily_precip_is_increasing)])

    @pytest.mark.parametrize("runid", ["precip_three_days"])
    def test_should_calculate_rsto_h2o(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        rsto_h2o = hourly_output["rsto_h2o"].values

        external_state: External_State_Shape = external_state
        print(rsto_h2o)
        assert all((r is not None and r > 0 for r in rsto_h2o))

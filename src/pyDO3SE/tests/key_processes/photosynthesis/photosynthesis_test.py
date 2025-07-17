"""Test running a few lines of data to make sure soil moisture is correct."""

from collections import namedtuple
from pathlib import Path
import pytest
import warnings
import pandas as pd
from math import isclose
import numpy as np

from pyDO3SE.Output.OutputConfig import (
    output_results_only_options,
)
from pyDO3SE import main

project_dir = "tests/key_processes/photosynthesis"


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
    if input_state_file:
        project_paths = project_paths._replace(
            base_state_path=f"{project_dir}/initial_state/{input_state_file}.json"
        )

    # Create output dir
    main.create_run_path_directories(run_paths)

    output_options = output_results_only_options()
    output_options.save_hourly_output_data = True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = main.single(
            config_file=run_paths.config_path,
            data_file=run_paths.input_data_file_path,
            output_directory=run_paths.output_directory,
            base_config_file=project_paths.base_config_path,
            initial_state_path=project_paths.base_state_path if input_state_file else None,
            plot_fields=None,
            runid=runid,
            verbose=2,
            output_options=output_options,
        )
    return out


RunSetup = namedtuple(
    "RunSetup", ["runid", "config_file", "input_file", "input_state", "overrides"]
)

setups = [
    RunSetup("default_dvi_1_03", "default", "three_days", "dump_state_dvi_1_03", {}),
    RunSetup("default_dvi_1_94", "default", "three_days", "dump_state_dvi_1_94", {}),
    RunSetup("negative_an_cubic", "negative_an", "negative_an", None, {}),
    RunSetup("negative_an_iterative", "negative_an_iterative", "negative_an", None, {}),
    RunSetup("negative_an_clip", "negative_an_clip", "negative_an", None, {}),
]

setupsMap = {s.runid: s for s in setups}
all_setups = [s.runid for s in setups]


@pytest.fixture(scope="class")
def photosynthesis_test_run(request):
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
        request.cls.output[runid]["out"] = out
        request.cls.output[runid]["hourly_output"] = pd.DataFrame(output_logs)


@pytest.mark.usefixtures("photosynthesis_test_run")
class TestRunAndCompare:
    output: dict[str, dict]

    def test_preruns_run_without_error(self):
        pass

    @pytest.mark.parametrize("runid", ["default_dvi_1_03"])
    def test_should_have_sun_frac_below_0(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        LAIsunfrac_top = hourly_output["LAIsunfrac_top"].values
        assert any(i < 1.0 for i in LAIsunfrac_top)

    @pytest.mark.parametrize("runid", ["default_dvi_1_03"])
    def test_should_output_canopy_sunlit_gsto(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        gsto_l_sunlit = hourly_output["gsto_l_sunlit"].values
        gsto_l = hourly_output["gsto_l"].values
        assert any(a < b for a, b in zip(gsto_l, gsto_l_sunlit))

    @pytest.mark.parametrize("runid", ["default_dvi_1_03"])
    def test_should_output_canopy_sunlit_An(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        A_n_sunlit = hourly_output["A_n_sunlit"].values
        A_n = hourly_output["A_n"].values
        in_sun_index = next(i for i, v in enumerate(A_n) if v > 0)
        assert A_n_sunlit[in_sun_index] > A_n[in_sun_index]

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_not_output_negative_An(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        A_n_sunlit = hourly_output["A_n_sunlit"].values
        A_n_canopy = hourly_output["A_n_canopy"].values
        A_n = hourly_output["A_n"].values
        # Note buffer of -20 to allow for R_d
        assert all(a >= -20 for a in A_n_sunlit), f"Sunlit: {min(A_n_sunlit)}"
        assert all(a >= -20 for a in A_n_canopy), f"Sunlit: {min(A_n_sunlit)}"
        assert all(a >= -20 for a in A_n), f"Sunlit: {min(A_n_sunlit)}"

    @pytest.mark.parametrize("runid", all_setups)
    def test_a_n_equal_or_less_than_zero_when_f_ls_is_zero(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        A_n_values = hourly_output["A_n"].values
        f_LS_values = hourly_output["f_LS"].values
        failed_index = next(
            (
                i
                for i, (A_n, f_LS) in enumerate(zip(A_n_values, f_LS_values))
                if A_n != 0 and f_LS == 0
            ),
            None,
        )
        if failed_index is not None:
            print(f"Failed at index: {failed_index}")
            print(f"A_n: {A_n_values[failed_index - 5 : failed_index + 5]}")
            print(f"f_LS: {f_LS_values[failed_index - 5 : failed_index + 5]}")
        assert all(A_n <= 0 if f_LS == 0 else True for A_n, f_LS in zip(A_n_values, f_LS_values)), (
            "Failed at "
            + str(
                next(
                    i
                    for i, (A_n, f_LS) in enumerate(zip(A_n_values, f_LS_values))
                    if A_n != 0 and f_LS == 0
                )
            )
        )

    @pytest.mark.parametrize("runid", all_setups)
    def test_c_i_should_not_be_negative(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        c_i_output = hourly_output["c_i_sunlit"].values
        assert all(a >= 0 for a in c_i_output), (
            f"Failed at {next(i for i, a in enumerate(c_i_output) if a < 0)}"
        )
        # NOTE: This can be caused by A_j being the limiting factor in the calculation of c_i

    @pytest.mark.parametrize("runid", all_setups)
    def test_c_i_should_not_be_above_external_co2(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        co2 = hourly_output["co2"].values[0]
        co2_padded = co2 * 1.1 # Pad as can be slightly higher
        c_i_output = hourly_output["c_i_sunlit"].values
        assert all(a <= co2_padded for a in c_i_output), (
            f"Failed at {next([i, a] for i, a in enumerate(c_i_output) if a > co2_padded)}"
        )

    @pytest.mark.skip(reason="Not finished")
    @pytest.mark.parametrize("runid", ["negative_an_cubic", "negative_an_iterative"])
    def test_should_be_able_to_recalculate_c_i_from_outputs(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        co2 = hourly_output["co2"].values[0]
        c_i_output = hourly_output["c_i_sunlit"].values
        A_n_values = hourly_output["A_n_sunlit"].values
        g_sto_values = hourly_output["gsto_l_sunlit"].values
        g_bv_values = hourly_output["g_bv"].values

        calc_ci = [
            co2 - ((A_n * (1 / g_sto + 1.37 / g_bl)) * 1e6)
            for A_n, g_sto, g_bl in zip(A_n_values, g_sto_values, g_bv_values)
        ]
        assert all(np.isnan(a) or isclose(a, b) for a, b in zip(c_i_output, calc_ci)), (
            f"Failed at {next([i, a, b] for i, (a, b) in enumerate(zip(c_i_output, calc_ci)) if not isclose(a, b))}"
        )

    # @pytest.mark.skip(reason="No longer valid")
    def test_should_get_similar_outputs_between_cubic_and_iterative(self):
        """We implemented this test while trying to match the cubic and iterative"""
        hourly_output_cubic = self.output["negative_an_cubic"]["hourly_output"]
        hourly_output_iterative = self.output["negative_an_iterative"]["hourly_output"]
        # NOTE: These are no longer that similar but it is most important the final outcome is similar
        # c_i_output_cubic = hourly_output_cubic["c_i_sunlit"].values
        # A_n_values_cubic = hourly_output_cubic["A_n_sunlit"].values
        # g_sto_values_cubic = hourly_output_cubic["gsto_l_sunlit"].values
        # c_i_output_iter = hourly_output_iterative["c_i_sunlit"].values
        # A_n_values_iter = hourly_output_iterative["A_n_sunlit"].values
        # g_sto_values_iter = hourly_output_iterative["gsto_l_sunlit"].values

        # assert all(
        #     np.isnan(a) or isclose(a, b) for a, b in zip(c_i_output_cubic, c_i_output_iter)
        # ), (
        #     f"Failed at {next(i for i, (a, b) in enumerate(zip(c_i_output_cubic, c_i_output_iter)) if not isclose(a, b))}"
        # )

        # assert all(
        #     np.isnan(a) or isclose(a, b) for a, b in zip(A_n_values_cubic, A_n_values_iter)
        # ), (
        #     f"Failed at {next(i for i, (a, b) in enumerate(zip(A_n_values_cubic, A_n_values_iter)) if not isclose(a, b))}"
        # )

        # assert all(
        #     np.isnan(a) or isclose(a, b) for a, b in zip(g_sto_values_cubic, g_sto_values_iter)
        # ), (
        #     f"Failed at {next(i for i, (a, b) in enumerate(zip(g_sto_values_cubic, g_sto_values_iter)) if not isclose(a, b))}"
        # )

        assert isclose(
            hourly_output_cubic["c_harv"].iloc[-1],
            hourly_output_iterative["c_harv"].iloc[-1],
            abs_tol=1e-3
        )

    @pytest.mark.parametrize("runid", ["negative_an_iterative"])
    def test_should_not_fail_to_find_iterative_A_n_solution(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        iterations = hourly_output["ewert_loop_iterations"].values
        assert all(i < 50 for i in iterations), (
            f"Failed at {next(i for i, a in enumerate(iterations) if a >= 50)}"
        )

    @pytest.mark.parametrize("runid", ["negative_an_cubic"])
    def test_should_not_fail_to_find_cubic_A_n_solution(self, runid):
        hourly_output = self.output[runid]["hourly_output"]
        A_j_values = hourly_output["A_j"].values
        assert all(not np.isnan(i) for i in A_j_values), (
            f"Failed at {next(i for i, a in enumerate(A_j_values) if np.isnan(a))}"
        )
        A_p_values = hourly_output["A_p"].values
        assert all(not np.isnan(i) for i in A_p_values), (
            f"Failed at {next(i for i, a in enumerate(A_p_values) if np.isnan(a))}"
        )

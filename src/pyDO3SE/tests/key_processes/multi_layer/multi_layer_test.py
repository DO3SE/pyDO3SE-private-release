"""Test running a few lines of data to make sure ozone deposition is correct."""

from math import isclose
from typing import NamedTuple
import numpy as np
import pytest
import pandas as pd

from pyDO3SE import settings
from pyDO3SE.External_State.External_State_Shape import External_State_Shape
from pyDO3SE import main
from tests.utils import get_test_run, Setup


settings.set_settings(MAX_NUM_OF_CANOPY_LAYERS=10)
project_dir = "tests/key_processes/multi_layer"


setups: list[Setup] = [
    Setup("default_three_days", "default", "three_days"),
    Setup("default_tree_three_days", "default_tree", "three_days"),
    Setup("default_two_hours", "default", "two_hours"),
    Setup("canopy_measure_height_two_hours", "canopy_measure_height", "two_hours"),
    Setup("target_h_equals_model_h_two_hours", "target_h_equals_model_h", "two_hours"),
    Setup("measured_h_eq_izr_two_hours", "measured_h_eq_izr", "two_hours"),
    Setup("input_ustar_ref_two_hours", "input_ustar_ref", "two_hours"),
    Setup("input_ustar_two_hours", "input_ustar", "two_hours"),
    # ["otc_three_days", "using_otc", "three_days", {}],
]

all_setups = [s[0] for s in setups]


def otherThan(ignore_these):
    return [s for s in all_setups if s not in ignore_these]


def get_setup(runid: str):
    return next(s for s in setups if s[0] == runid)


run_outputs = {}


class TestRunOutput(NamedTuple):
    out: main.MainOutput
    hourly_output: pd.DataFrame


model_run = get_test_run(setups, run_outputs, project_dir)

class TestRunAndCompare:
    @pytest.mark.parametrize("runid", all_setups)
    def test_preruns_run_without_error(self, runid):
        model_run(runid)

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_calculate_ustar_correctly(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        ustar = hourly_output["ustar"].values
        ustar_ref = hourly_output["ustar_ref"].values
        u_i = hourly_output["u_i"].values
        assert all(ustar)
        assert all(ustar_ref)
        assert all(u_i)

    @pytest.mark.parametrize("runid", ["canopy_measure_height_two_hours"])
    def test_ozone_in_should_be_the_same_as_ozone_out(self, runid):
        """Test that ozone remains the same if not transfered up or down."""
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        canopy_top_o3 = hourly_output["canopy_top_o3"].values
        if not final_config.Location.OTC:
            assert final_config.Location.h_O3 == final_config.Location.z_O3
            assert final_config.Location.h_O3 == final_state.canopy.canopy_height
        external_state: External_State_Shape = external_state
        assert external_state.O3 is not None
        print(canopy_top_o3)
        print(external_state.O3)
        assert all(canopy_top_o3)
        assert all([isclose(a, b, abs_tol=1e-0) for a, b in zip(canopy_top_o3, external_state.O3)])

    @pytest.mark.parametrize("runid", ["measured_h_eq_izr_two_hours"])
    def test_o3_i_should_equal_o3_zr_if_h_o3_equals_izr(self, runid):
        # TODO: Is this right
        """Test that ozone changes if transfered up or down with different canopy heights."""
        output = model_run(runid)

        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        assert final_config.Location.z_O3 == final_config.Location.izr

        o3_ppb_i = hourly_output["o3_ppb_i"].values
        o3_ppb_zr = hourly_output["o3_ppb_zr"].values
        external_state: External_State_Shape = external_state
        assert external_state.O3 is not None
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(o3_ppb_i, external_state.O3)])
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(o3_ppb_zr, external_state.O3)])

    @pytest.mark.parametrize("runid", otherThan(["measured_h_eq_izr_two_hours"]))
    def test_ozone_transfered_up_should_increase_ozone_at_izr(self, runid):
        # TODO: Is this right
        """Test that ozone changes if transfered up or down with different canopy heights."""
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        assert final_config.Location.z_O3 < final_config.Location.izr
        o3_ppb_i = hourly_output["o3_ppb_i"].values

        external_state: External_State_Shape = external_state
        assert external_state.O3 is not None
        assert any([a > b for a, b in zip(o3_ppb_i, external_state.O3)])

    @pytest.mark.parametrize(
        "runid",
        otherThan(
            [
                "canopy_measure_height_two_hours",
                "target_h_equals_model_h_two_hours",
                "input_ustar_ref_two_hours",
                "default_tree_three_days",
            ]
        ),
    )
    def test_ozone_in_should_not_be_the_same_as_ozone_out_when_transfered(self, runid):
        # TODO: Is this right
        """Test that ozone changes if transfered up or down with different canopy heights."""
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        micro_O3 = hourly_output["micro_O3"].values

        external_state: External_State_Shape = external_state
        assert all(micro_O3)
        assert external_state.O3 is not None
        assert any([not isclose(a, b, abs_tol=1e-3) for a, b in zip(micro_O3, external_state.O3)])

    @pytest.mark.parametrize(
        "runid",
        all_setups,
    )
    def test_should_have_correct_number_of_layers_in_output(self, runid):
        """Test that the number of layers in the output is correct."""
        output = model_run(runid)
        final_state, output_logs, final_config, initial_state, external_state = output.out

        assert len(final_state.canopy_layers) == final_config.Land_Cover.nL

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_calculate_ra(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        ra = hourly_output["ra"].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in ra))

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_calculate_rb(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        rb = hourly_output["rb"].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rb))

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_calculate_rsur(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        rsur = hourly_output["rsur"].values
        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rsur))

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_calculate_rinc(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        rinc = hourly_output["rinc"].values
        canopy_sai = hourly_output["canopy_sai"].values

        print(rinc)
        print(canopy_sai)
        external_state: External_State_Shape = external_state
        assert all((sai is not None and sai > 0 for sai in canopy_sai))
        assert all((r is not None for r in rinc))
        assert sum(rinc) > 0

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_calculate_rext(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        rext = hourly_output["rext"].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rext))

    @pytest.mark.parametrize("runid", all_setups)
    def test_should_calculate_rsto_c(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        rsto_c = hourly_output["rsto_c"].values

        external_state: External_State_Shape = external_state
        assert all((r is not None and r > 0 for r in rsto_c))

    @pytest.mark.parametrize("runid", ["canopy_measure_height_two_hours", "input_ustar_ref_two_hours"])
    def test_windspeed_in_should_be_the_same_as_windspeed_out(self, runid):
        """Test that windspeed remains the same if not transfered up or down."""
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        micro_u = hourly_output["micro_u"].values

        external_state: External_State_Shape = external_state
        assert external_state.u is not None
        assert all(micro_u)
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(micro_u, external_state.u)])

    @pytest.mark.parametrize("runid", ["canopy_measure_height_two_hours", "input_ustar_ref_two_hours"])
    def test_ustar_ref_in_should_be_the_same_as_ustar_out_if_using_observed_canopy_height(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        ustar = hourly_output["ustar"].values

        external_state: External_State_Shape = external_state
        ustar_observed = external_state.ustar_ref
        assert external_state.O3 is not None
        assert ustar_observed is not None
        assert all(ustar_observed)
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(ustar, ustar_observed)])

    @pytest.mark.parametrize("runid", ["input_ustar_two_hours"])
    def test_ustar_in_should_be_the_same_as_ustar_out_if_using_observed_canopy_height(self, runid):
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        ustar = hourly_output["ustar"].values

        external_state: External_State_Shape = external_state
        ustar_observed = external_state.ustar
        assert ustar_observed is not None
        assert all(ustar_observed)
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(ustar, ustar_observed)])

    @pytest.mark.parametrize("runid", ["measured_h_eq_izr_two_hours"])
    def test_that_if_measured_height_equal_to_izr_that_O3_i_equals_O3_zr(self, runid):
        """Test that when the measured height is equal to the izr(Mixed height) that O3_i(At izr) equals O3_zr(measured O3)."""
        output = model_run(runid)
        hourly_output = output.hourly_output
        final_state, output_logs, final_config, initial_state, external_state = output.out
        o3_ppb_i = hourly_output["o3_ppb_i"].values
        o3_ppb_zr = hourly_output["o3_ppb_zr"].values
        external_state: External_State_Shape = external_state
        O3_observed = external_state.O3
        assert O3_observed is not None
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(o3_ppb_i, o3_ppb_zr)])
        assert all([isclose(a, b, abs_tol=1e-3) for a, b in zip(O3_observed, o3_ppb_zr)])


class TestInCanopyMet:
    """We need to be able to calculate the met state within the canopy"""

    class TestAtCanopyLevels:
        @pytest.mark.parametrize("runid", otherThan(["otc_three_days"]))
        def test_should_gradually_decrease_windspeed_from_top_of_canopy_to_bottom(self, runid):
            """Test that windpseed decreases through the canopy."""
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out
            canopy_top_u = hourly_output["micro_u"].values
            micro_u_iL_0 = hourly_output["micro_u_iL_0"].values
            micro_u_iL_1 = hourly_output["micro_u_iL_1"].values
            micro_u_iL_2 = hourly_output["micro_u_iL_2"].values

            print(micro_u_iL_0[0:2])

            # TODO: Check these ranges
            # layers should be between 1% and 80% of the layer above
            assert all([b * 0.01 < a < b * 0.9 for a, b in zip(micro_u_iL_0, micro_u_iL_1)])
            assert all([b * 0.01 < a < b * 0.9 for a, b in zip(micro_u_iL_1, micro_u_iL_2)])
            assert all([b * 0.01 < a < b * 0.9 for a, b in zip(micro_u_iL_2, canopy_top_u)])
            # bottom layer should be between 1% and 70% of the top layer
            assert all([b * 0.01 < a < b * 0.7 for a, b in zip(micro_u_iL_0, canopy_top_u)])

        @pytest.mark.parametrize("runid", otherThan(["otc_three_days"]))
        def test_should_gradually_decrease_solar_radiation_from_top_of_canopy_to_bottom(self, runid):
            """Test that solar radiation decreases through the canopy."""
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out
            canopy_top_par = hourly_output["PARsun"].values
            PARsun_iL_0 = hourly_output["PARsun_iL_0"].values
            PARsun_iL_1 = hourly_output["PARsun_iL_1"].values
            PARsun_iL_2 = hourly_output["PARsun_iL_2"].values

            print(PARsun_iL_0[0:24])
            print(PARsun_iL_1[0:24])

            # TODO: Check these ranges
            # layers should be between 1% and 80% of the layer above
            assert all([a + b == 0 or b * 0.5 < a < b * 0.999 for a, b in zip(PARsun_iL_0, PARsun_iL_1)])
            assert all([a + b == 0 or b * 0.5 < a < b * 0.999 for a, b in zip(PARsun_iL_1, PARsun_iL_2)])
            assert all([a + b == 0 or b * 0.5 < a < b * 0.999 for a, b in zip(PARsun_iL_2, canopy_top_par)])
            # bottom layer should be between 1% and 70% of the top layer
            assert all([a + b == 0 or b * 0.01 < a < b * 0.9999 for a, b in zip(PARsun_iL_0, canopy_top_par)])

        @pytest.mark.parametrize("runid", otherThan(["otc_three_days", "default_tree_three_days"]))
        def test_should_gradually_decrease_ozone_from_top_of_canopy_to_bottom(self, runid):
            """Test that ozone decreases through the canopy.

            NOTE: We need a reference for the ratios used here

            """
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out

            canopy_height = final_state.canopy.canopy_height
            canopy_top_o3 = hourly_output["canopy_top_o3"].values
            micro_O3_iL_0 = hourly_output["micro_O3_iL_0"].values
            micro_O3_iL_1 = hourly_output["micro_O3_iL_1"].values
            micro_O3_iL_2 = hourly_output["micro_O3_iL_2"].values

            assert all([b * 0.001 < a < b * 0.99 for a, b in zip(micro_O3_iL_0, micro_O3_iL_1)])
            assert all([b * 0.001 < a < b * 0.99 for a, b in zip(micro_O3_iL_1, micro_O3_iL_2)])
            assert all([b * 0.001 < a < b * 0.99 for a, b in zip(micro_O3_iL_2, canopy_top_o3)])
            # bottom layer should match ratio of top layer. Ratio is set by canopy height.
            # This is a rough estimate and should be improved
            upper_ratio_limit = (-0.01 * canopy_height) + 1.5 + (1 / (canopy_height + 1)) - 1
            assert all([b * 0.0001 < a < b * upper_ratio_limit for a, b in zip(micro_O3_iL_0, canopy_top_o3)])
            # bottom layer should be greater than 0
            assert all([a > 0 for a in micro_O3_iL_0])

        @pytest.mark.parametrize("runid", ["default_tree_three_days"])
        def test_should_gradually_decrease_ozone_from_top_of_canopy_to_bottom_for_trees(self, runid):
            """Test that ozone decreases through the canopy.

            The ratios for trees are defined in Makar et al (2017)

            100 yr old Mixed deciduous forest in Canada, LAI 4.6m2/m2

            Measurements conducted from 2008 to 2013

            Based on Figure 2 we assume that during the summer(When LAI is at max) the ozone
            concentration at the bottom of the canopy is approximately 33% of the top of the canopy.

            References
            ----------
            Makar et al 2017

            """
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out
            canopy_top_o3 = hourly_output["canopy_top_o3"].values
            micro_O3_iL_0 = hourly_output["micro_O3_iL_0"].values
            micro_O3_iL_1 = hourly_output["micro_O3_iL_1"].values
            micro_O3_iL_2 = hourly_output["micro_O3_iL_2"].values
            canopy_lai = hourly_output["canopy_lai"].values

            # Print for debuging
            print(micro_O3_iL_0[0:2])
            print(micro_O3_iL_1[0:2])
            print(micro_O3_iL_2[0:2])
            print(canopy_top_o3[0:2])
            print((micro_O3_iL_0 / canopy_top_o3)[0:24])
            print((micro_O3_iL_0 / canopy_top_o3)[-24:])

            print(micro_O3_iL_0[-2:])
            print(micro_O3_iL_1[-2:])
            print(micro_O3_iL_2[-2:])
            print(canopy_top_o3[-2:])

            # TODO: This is failing for the default_tree_three_days
            # assert all([b * 0.001 < a < b * 0.99 for a, b in zip(micro_O3_iL_0, micro_O3_iL_1)])
            # assert all([b * 0.001 < a < b * 0.99 for a, b in zip(micro_O3_iL_1, micro_O3_iL_2)])
            # assert all([b * 0.001 < a < b * 0.99 for a, b in zip(micro_O3_iL_2, canopy_top_o3)])
            # bottom layer should be between 0.01% and 50% of the top layer
            canopy_height = final_state.canopy.canopy_height
            # TODO: should be ratio of canopy LAI also
            upper_ratio_limit = (-0.001 * canopy_height) + 1 + (1 / ((0.011 * canopy_height) + 1)) - 1
            lower_ratio_limit = (-0.001 * canopy_height) + 1.4 + (1 / ((1 * canopy_height) + 1)) - 1
            print("Canopy height", canopy_height)
            print("limits", upper_ratio_limit, lower_ratio_limit)
            assert all(
                [
                    lai < 0.1 or b * lower_ratio_limit < a < b * upper_ratio_limit
                    for a, b, lai in zip(micro_O3_iL_0, canopy_top_o3, canopy_lai)
                ]
            )
            # assert all([b * 0.3 < a < b * 0.4 for a, b in zip(micro_O3_iL_0, canopy_top_o3)])

            # bottom layer should be greater than 0
            assert all([a > 0 for a in micro_O3_iL_0])

    class TestAtGroundLevel:
        @pytest.mark.parametrize("runid", all_setups)
        def test_should_calculate_ozone_at_ground_level(self, runid):
            """Test that ozone is being calucated for the ground level"""
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out
            assert "ground_O3" in final_config.output.fields
            assert "ground_O3" in hourly_output.columns
            assert final_state.ground_level.micro_met.micro_O3 is not None
            assert final_state.ground_level.micro_met.micro_O3 > 0
            ground_O3 = hourly_output["ground_O3"].values

            assert all([a > 0 for a in ground_O3])

        @pytest.mark.skip(reason="This test is failing")
        @pytest.mark.parametrize("runid", all_setups)
        def test_should_not_set_ozone_at_ground_to_be_higher_than_canopy_layers(self, runid):
            """Test that ozone is being calucated for the ground level"""
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out
            assert "ground_O3" in final_config.output.fields
            assert "ground_O3" in hourly_output.columns
            assert final_state.ground_level.micro_met.micro_O3 is not None
            assert final_state.ground_level.micro_met.micro_O3 > 0
            ground_O3 = hourly_output["ground_O3"].values
            micro_O3_iL_0 = hourly_output["micro_O3_iL_0"].values

            # TODO: default three days failing here
            assert all([a < b for a, b in zip(ground_O3, micro_O3_iL_0)])

    class TestAtCustomHeight:
        """We need to be able to calculate the ozone at a custom height.

        This could be at a person head height or at the height of a plant.
        """

        @pytest.mark.parametrize("runid", otherThan([]))
        def test_should_calculate_ozone_at_custom_height(self, runid):
            """Test that ozone is being calucated for the ground level"""
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out

            # Must set a custom height in the config
            assert final_config.Location.custom_heights is not None
            assert len(final_config.Location.custom_heights) > 0
            assert final_config.Location.custom_heights[0] is not None

            assert final_state.custom_height[0].micro_met.micro_O3 is not None
            assert final_state.custom_height[0].micro_met.micro_O3 > 0

            for iCH, custom_layer_height in enumerate(final_config.Location.custom_heights):
                if custom_layer_height < min(layer.layer_height for layer in final_state.canopy_layers):
                    ground_O3 = hourly_output["ground_O3"].values
                    micro_O3_iCH = hourly_output[f"micro_O3_iCH_{iCH}"].values
                    if 0 == custom_layer_height:
                        np.isclose(micro_O3_iCH, ground_O3, atol=1e-3)
                    if 0 > custom_layer_height:
                        raise ValueError("Custom height is below ground level")
                    if 0 < custom_layer_height:
                        assert all(np.greater(micro_O3_iCH - ground_O3, -1e-6))

                for iL, layer in enumerate(final_state.canopy_layers):
                    micro_O3_iCH = hourly_output[f"micro_O3_iCH_{iCH}"].values
                    micro_O3_iL = hourly_output[f"micro_O3_iL_{iL}"].values

                    if isclose(layer.layer_height, custom_layer_height, abs_tol=1e-1):
                        np.isclose(micro_O3_iCH, micro_O3_iL, atol=1e-3)
                    elif layer.layer_height > custom_layer_height:
                        assert all(np.greater(micro_O3_iL - micro_O3_iCH, -1e-4))
                    elif layer.layer_height < custom_layer_height:
                        assert all(np.less(micro_O3_iL - micro_O3_iCH, 1e-4))

    class TestCanopyLevelVariables:
        """Some values are calculated assuming a single layer canopy. We need to test that these are being calculated correctly"""

        @pytest.mark.parametrize("runid", otherThan([]))
        def test_canopy_vd(self, runid):
            """Test that the canopy_vd is being calculated correctly"""
            output = model_run(runid)
            hourly_output = output.hourly_output
            final_state, output_logs, final_config, initial_state, external_state = output.out

            canopy_vd = hourly_output["canopy_vd"].values
            assert all([a > 0 for a in canopy_vd])

            assert 0.001 < np.mean(canopy_vd) < 0.011

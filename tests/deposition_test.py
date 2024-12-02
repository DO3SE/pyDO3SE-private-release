import pytest
from math import isclose
import numpy as np
from do3se_met.deposition import (
    calc_canopy_ozone_concentration,
    calc_multi_layer_O3_ozone_concentration,
    calc_ozone_at_custom_height_linear,
)

# TODO: Fix and compare tests for multilayer setup
# @pytest.mark.parametrize(
#     ['O3', 'canopy_height', 'u_i', 'h_O3_in', 'z_O3', 'Rtotal_top_layer', 'expected_out'],
#     [
#         (33.904, 1.0, 2.9349, None, 40.0, 99.9001, [34.871, 19.589]),
#     ],
# )
# def test_calc_canopy_ozone_concentration_multilayer(
#     O3,
#     canopy_height,
#     u_i,
#     h_O3_in,
#     z_O3,
#     Rtotal_top_layer,
#     expected_out,
# ):
#     """Test the output of calc_canopy_ozone_concentration."""
#     out = calc_canopy_ozone_concentration_multilayer(O3, canopy_height, u_i, h_O3_in, z_O3, Rtotal_top_layer)

#     assert out.O3_i > O3
#     assert out.micro_O3 < O3
#     assert isclose(out.O3_i, expected_out[0], abs_tol=1e-4)
#     assert isclose(out.micro_O3, expected_out[1], abs_tol=1e-4)


# @pytest.mark.skip(reason="Currently does not align with ui model")
# def test_calc_canopy_ozone_concentration_compare_ui_alt():
#     """Test the output of calc_canopy_ozone_concentration matches ui outp."""
#     out = calc_canopy_ozone_concentration_multilayer(
#         O3=33.904,
#         canopy_height=1.0,
#         u_i=2.9348522711,
#         h_O3_in=None,  # Use canopy height
#         z_O3=50.0,
#         Rtotal_top_layer=99.900131118
#     )
#     # TODO: Get these values from ui
#     assert isclose(out.O3_i, 33.904, abs_tol=1e-4)
#     assert isclose(out.micro_O3, 28.27369, abs_tol=1e-4)


class TestCalcCanopyOzoneConcentration:
    @pytest.mark.parametrize(
        [
            "O3",
            "Rsur_top_layer",
            "Rb_top_layer",
            "Rb_ref",
            "Ra_ref_canopy",
            "Ra_tar_canopy",
            "Ra_ref_measured",
            "expected_o3_i",
            "expected_micro_o3",
        ],
        [(26.961, 454.6, 7.449, 7.449, 9.304, 9.304, 9.304, 27.503, 26.961)],
    )
    def test_calc_canopy_ozone_concentration_compare_ui(
        self,
        O3,
        Rsur_top_layer,
        Rb_top_layer,
        Ra_ref_canopy,
        Ra_tar_canopy,
        Ra_ref_measured,
        Rb_ref,
        expected_o3_i,
        expected_micro_o3,
    ):
        """Test the output of calc_canopy_ozone_concentration matches ui outp."""

        out = calc_canopy_ozone_concentration(
            O3_ppb_zR=O3,
            Rsur_top_layer=Rsur_top_layer,
            Ra_ref_canopy=Ra_ref_canopy,
            Ra_tar_canopy=Ra_tar_canopy,
            Ra_ref_measured=Ra_ref_measured,
            Rb_ref=Rb_ref,
            Rb_top_layer=Rb_top_layer,
            Ra_tar_canopy_top=Ra_tar_canopy,
            Rsur_ref=Rsur_top_layer,
        )
        # TODO: Get these values from ui
        assert isclose(out.O3_i, expected_o3_i, abs_tol=1e-3)
        assert isclose(out.micro_O3, expected_micro_o3, abs_tol=1e-3)


# @pytest.mark.parametrize(['nL', 'O3_in', 'rm_Ra', 'rm_Rinc', 'rm_Rsur', 'rm_Rgs', 'expected_output'], [
#     (3, 33.3, 1.1, [99, 99, 99], [100, 100, 100], 10, [32.9358, 16.4754, 1.5115]),
#     (4, 23.401, 249.54, [0.168, 0.168, 0.168, 0.168], [
#      2551.67, 2551.67, 2551.67, 2551.67], 200, [9.2233, 9.2143, 9.2060, 9.1983]),
#     (4, 23.561, 283.95, [0.1452420, 0.1452420, 0.1452420, 0.1452420], [
#      2558.797, 2558.797, 2558.797, 2558.797], 200, [8.571, 8.564, 8.557, 8.551]),
#     (4, 30.0, 0, [0.0, 0.0, 0.0, 0.0], [10000, 10000, 10000, 10000], 200, [30.0, 30.0, 30.0, 30.0]),
# ])
# def test_calc_multi_layer_ozone_concentration(
#     nL, O3_in, rm_Ra, rm_Rinc, rm_Rsur, rm_Rgs, expected_output
# ):
#     """Test the output of calc_multi_layer_ozone_concentration."""
#     O3_out = calc_multi_layer_O3_ozone_concentration(
#         nL=nL,
#         O3_in=O3_in,
#         rm_Ra=rm_Ra,
#         rm_Rinc=rm_Rinc,
#         rm_Rsur=rm_Rsur,
#         rm_Rgs=rm_Rgs,
#     )
#     print(O3_out)
#     for iL in range(nL):
#         assert isclose(O3_out[iL], expected_output[iL], abs_tol=1e-3)
#     assert all([a <= b for a, b in zip(O3_out[1:], O3_out[0:-1])])


class TestCalcMultiLayerOzoneConcentration:
    def test_concentration_should_be_highest_at_top_layer(self):
        nL = 3
        O3_in = 10
        rm_Ra = 10  # Assume O3 already at canopy level
        rm_Rinc = [1000.0 for _ in range(nL)]
        rm_Rsur = [6000.0 for _ in range(nL)]
        rm_Rgs = 200  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        assert all([a <= b for a, b in zip(O3_out[1:], O3_out[0:-1])])

    def test_concentration_at_top_of_canopy_should_be_close_to_O3_in(self):
        nL = 3
        O3_in = 10
        rm_Ra = 10  # Assume O3 already at canopy level
        rm_Rinc = [1000.0 for _ in range(nL)]
        rm_Rsur = [6000.0 for _ in range(nL)]
        rm_Rgs = 200  # Soil resistance is constant


        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        assert isclose(O3_out[0], O3_in, abs_tol=10)

    def test_concentration_at_each_layer(self):
        nL = 3
        O3_in = 10
        rm_Ra = 0  # Assume O3 already at canopy level
        rm_Rinc = [1000.0 for _ in range(nL)]
        rm_Rsur = [6000.0 for _ in range(nL)]
        rm_Rgs = 200  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        assert isclose(O3_out[0], O3_in, abs_tol=10)
        assert isclose(O3_out[1], 6.00, abs_tol=1e-3)
        assert isclose(O3_out[2], 2.999, abs_tol=1e-3)

    def test_concentration_at_each_layer_low_resistance(self):
        nL = 3
        O3_in = 10
        rm_Ra = 10  # Assume O3 already at canopy level
        rm_Rinc = [0.0001 for _ in range(nL)]
        rm_Rsur = [0.0001 for _ in range(nL)]
        rm_Rgs = 0.0001  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        assert isclose(O3_out[0], O3_in, abs_tol=10)
        assert isclose(O3_out[1], 0, abs_tol=1e-3)
        assert isclose(O3_out[2], 0, abs_tol=1e-3)

    def test_concentration_at_each_layer_high_surface_resistance(self):
        nL = 3
        O3_in = 10
        rm_Ra = 0  # Assume O3 already at canopy level
        rm_Rinc = [0.00000001 for _ in range(nL)]
        rm_Rsur = [100900000000.0 for _ in range(nL)]
        rm_Rgs = 200  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        # if Rinc is low and Rsur is high, O3 should be close to input value
        print(O3_out)
        assert isclose(O3_out[0], O3_in, abs_tol=1e-3)
        assert isclose(O3_out[1], O3_in, abs_tol=1e-3)
        assert isclose(O3_out[2], O3_in, abs_tol=1e-3)
        assert isclose(O3_out[3], O3_in, abs_tol=1e-3)

    def test_concentration_at_each_layer_high_internal_resistance(self):
        nL = 3
        O3_in = 10
        rm_Ra = 0  # Assume O3 already at canopy level
        rm_Rinc = [1000000000 for _ in range(nL)]
        rm_Rsur = [0.00000001 for _ in range(nL)]
        rm_Rgs = 200  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        # if Rinc is low and Rsur is high, O3 should be close to input value
        print(O3_out)
        assert isclose(O3_out[0], O3_in, abs_tol=1e-3)
        assert isclose(O3_out[1], 0, abs_tol=1e-3)
        assert isclose(O3_out[2], 0, abs_tol=1e-3)
        assert isclose(O3_out[3], 0, abs_tol=1e-3)

    def test_should_contain_value_at_ground_level(self):
        nL = 5
        O3_out = calc_multi_layer_O3_ozone_concentration(
            **{
                "nL": nL,
                "O3_in": 30.167430449675667,
                "rm_Ra": 34.260757032281624,
                "rm_Rinc": [
                    898.2900869477606,
                    1275.278962517186,
                    1810.4802177712204,
                    2570.2914540917577,
                    3648.9756110728913,
                ],
                "rm_Rsur": [
                    6207.1675125913,
                    6220.793784823008,
                    6240.1386474315605,
                    6267.602043361151,
                    6306.591108894859,
                ],
                "rm_Rgs": 200,
            }
        )
        assert len(O3_out) == nL + 1
        bottom_layer_index = nL
        assert O3_out[bottom_layer_index] < O3_out[bottom_layer_index - 1]


    def test_should_have_correct_ratio_between_top_and_bottom_of_canopy(self):
        nL = 5
        O3_out = calc_multi_layer_O3_ozone_concentration(
            **{
                "nL": nL,
                "O3_in": 30.167430449675667,
                "rm_Ra": 34.260757032281624,
                "rm_Rinc": [
                    898.2900869477606,
                    1275.278962517186,
                    1810.4802177712204,
                    2570.2914540917577,
                    3648.9756110728913,
                ],
                "rm_Rsur": [
                    6207.1675125913,
                    6220.793784823008,
                    6240.1386474315605,
                    6267.602043361151,
                    6306.591108894859,
                ],
                "rm_Rgs": 200,
            }
        )
        assert len(O3_out) == nL + 1
        gound_layer_index = nL
        top_layer_index = 0
        print(O3_out)
        assert O3_out[top_layer_index] > O3_out[gound_layer_index]
        ratio = O3_out[gound_layer_index] /O3_out[top_layer_index]
        # TODO: This ratio is not correct
        assert 0.2 < ratio < 0.3


class TestCalcOzoneAtCustomHeight:
    def test_should_return_O3_at_custom_height(self):
        x = np.array([0, 1, 2])
        y = np.array([10, 5, 2.5])

        val = calc_ozone_at_custom_height_linear(y, x, 0.5)
        assert val < 10
        assert val > 5

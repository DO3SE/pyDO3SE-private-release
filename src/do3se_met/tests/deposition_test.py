import pytest
import numpy as np
from math import isclose

from do3se_met.deposition import (
    calc_canopy_ozone_concentration,
    calc_canopy_ozone_concentration_alt,
    calc_multi_layer_O3_ozone_concentration,
)

# TODO: Fix and compare tests for alt setup
# @pytest.mark.parametrize(
#     ['O3', 'canopy_height', 'u_i', 'h_O3_in', 'z_O3', 'Rtotal_top_layer', 'expected_out'],
#     [
#         (33.904, 1.0, 2.9349, None, 40.0, 99.9001, [34.871, 19.589]),
#     ],
# )
# def test_calc_canopy_ozone_concentration_alt(
#     O3,
#     canopy_height,
#     u_i,
#     h_O3_in,
#     z_O3,
#     Rtotal_top_layer,
#     expected_out,
# ):
#     """Test the output of calc_canopy_ozone_concentration."""
#     out = calc_canopy_ozone_concentration_alt(O3, canopy_height, u_i, h_O3_in, z_O3, Rtotal_top_layer)

#     assert out.O3_i > O3
#     assert out.micro_O3 < O3
#     assert isclose(out.O3_i, expected_out[0], abs_tol=1e-4)
#     assert isclose(out.micro_O3, expected_out[1], abs_tol=1e-4)


# @pytest.mark.skip(reason="Currently does not align with ui model")
# def test_calc_canopy_ozone_concentration_compare_ui_alt():
#     """Test the output of calc_canopy_ozone_concentration matches ui outp."""
#     out = calc_canopy_ozone_concentration_alt(
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

    @pytest.mark.parametrize([
        'O3',
        'Rsur_top_layer',
        'Rb_top_layer',
        'Rb_ref',
        'Ra_ref_canopy',
        'Ra_tar_canopy',
        'Ra_ref_measured',
        'expected_o3_i',
        'expected_micro_o3',
    ], [
        (26.961, 454.6, 7.449, 7.449, 9.304, 9.304, 9.304, 27.503, 26.961)
    ],
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
        rm_Ra = 0  # Assume O3 already at canopy level
        rm_Rinc = [10 for _ in range(nL)]
        rm_Rsur = [10 for _ in range(nL)]
        rm_Rgs = 100  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        assert all([a <= b for a, b in zip(O3_out[1:], O3_out[0:-1])])

    def test_concentration_at_top_of_canopy_should_equal_O3_in(self):
        nL = 3
        O3_in = 10
        rm_Ra = 0  # Assume O3 already at canopy level
        rm_Rinc = [10 for _ in range(nL)]
        rm_Rsur = [10 for _ in range(nL)]
        rm_Rgs = 100  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        assert O3_out[0] == O3_in

    def test_concentration_at_each_layer(self):
        nL = 3
        O3_in = 10
        rm_Ra = 0  # Assume O3 already at canopy level
        rm_Rinc = [10 for _ in range(nL)]
        rm_Rsur = [10 for _ in range(nL)]
        rm_Rgs = 100  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        assert O3_out[0] == O3_in
        assert isclose(O3_out[1], 5.0, abs_tol=1e-3)
        assert isclose(O3_out[2], 4.545, abs_tol=1e-3)

    def test_concentration_at_each_layer_low_resistance(self):
        nL = 3
        O3_in = 10
        rm_Ra = 0  # Assume O3 already at canopy level
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
        assert O3_out[0] == O3_in
        assert isclose(O3_out[1], 5.0, abs_tol=1e-3)
        assert isclose(O3_out[2], 2.5, abs_tol=1e-3)

    def test_concentration_at_each_layer_high_surface_resistance(self):
        nL = 3
        O3_in = 10
        rm_Ra = 0  # Assume O3 already at canopy level
        rm_Rinc = [0.0 for _ in range(nL)]
        rm_Rsur = [100000 for _ in range(nL)]
        rm_Rgs = 0.0001  # Soil resistance is constant

        O3_out = calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=O3_in,
            rm_Ra=rm_Ra,
            rm_Rinc=rm_Rinc,
            rm_Rsur=rm_Rsur,
            rm_Rgs=rm_Rgs,
        )
        # if Rinc is 0 then all layers have input O3
        assert isclose(O3_out[0], O3_in, abs_tol=1e-3)
        assert isclose(O3_out[1], O3_in, abs_tol=1e-3)
        assert isclose(O3_out[2], O3_in, abs_tol=1e-3)


class TestCalcCanopyOzoneConcentrationAlt:

    @pytest.mark.parametrize([
        'O3',
        'z_O3',
        'Rsur_top_layer',
        'Rb_top_layer',
        'Ra_top_layer',
        'ustar_ref_O3',
        'L',
        'u_i',
        'O3_h',
        'expected_o3_i',
        'expected_micro_o3',
    ], [
        # (33.904, 2.9348522711, 50.0, 276.901275635,
        # 31.2906856537, 64.121711731, 0.19406299293, 816.99, [33.904, 28.254]),
        (26.961, 20, 454.6, 7.449, 9.304, 0.5486, 1 / -0.002034,
         10.0, 20, 27.109, 26.566),  # micro O3 should be 35.5
        # This should match emep. Check Ra calcs
        (26.961, 45, 454.6, 7.449, 9.304, 0.5486, 1 / -0.002034,
         10.0, 20, 26.961, 26.421),  # micro O3 should be 35.5
    ],
    )
    def test_calc_canopy_ozone_concentration_compare_ui(
        self,
        O3,
        z_O3,
        Rsur_top_layer,
        Rb_top_layer,
        Ra_top_layer,
        ustar_ref_O3,
        L,
        u_i,
        O3_h,
        expected_o3_i,
        expected_micro_o3,
    ):
        """Test the output of calc_canopy_ozone_concentration matches ui outp."""

        O3_d = O3_h * 0.78
        O3_z0 = min(O3_h * 0.07, 1.0)
        canopy_height = O3_h
        Rtotal_top_layer = 400  # This should be worked out from resistance of all layers!
        out = calc_canopy_ozone_concentration_alt(
            O3=O3,
            ustar_ref=ustar_ref_O3,
            canopy_height=canopy_height,
            u_i=u_i,
            z_O3=z_O3,
            Rtotal_top_layer=Rtotal_top_layer,
            O3_d=O3_d,
            O3_z0=O3_z0,
            L=L,
            izr=45,
            ra_method="heat_flux",
        )
        # NOTE: Low precision because matching not alt outputs
        assert isclose(out.O3_i, expected_o3_i, abs_tol=1e-0)
        assert isclose(out.micro_O3, expected_micro_o3, abs_tol=1e-0)


# # @pytest.mark.parametrize(['nL', 'O3_in', 'rm_Ra', 'rm_Rinc', 'rm_Rsur', 'rm_Rgs', 'expected_output'], [
# #     (3, 33.3, 1.1, [99,99,99],[100,100,100],10,[32.714, 12.3684, 4.2675]),
# #     (4, 23.401, 249.54, [0.168,0.168,0.168,0.168],[2551.67,2551.67,2551.67,2551.67],200,[8.885,8.876,8.867,8.859]),
# # ])
# # def test_calc_multi_layer_ozone_concentration(
# #     nL, O3_in, rm_Ra, rm_Rinc, rm_Rsur, rm_Rgs, expected_output
# # ):
# #     """Test the output of calc_multi_layer_ozone_concentration."""
# #     O3_out = calc_multi_layer_O3_ozone_concentration(
# #         nL=nL,
# #         O3_in=O3_in,
# #         rm_Ra=rm_Ra,
# #         rm_Rinc=rm_Rinc,
# #         rm_Rsur=rm_Rsur,
# #         rm_Rgs=rm_Rgs,
# #     )
# #     for iL in range(nL):
# #         assert isclose(O3_out[iL], expected_output[iL], abs_tol=1e-3)

# #     assert all([a > b for a, b in zip(O3_out[1:], O3_out[0:-1])])

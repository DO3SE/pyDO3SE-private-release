"""Tests for met modules wind helpers."""

import pytest
from math import isclose


from do3se_met.wind import (
    calc_monin_obukhov_length,
    ustar_from_velocity,
    calc_windspeed,
    calc_layer_windspeed,
    calc_ustar_and_L,
    ustar_from_velocity_simple,
)


def test_ustar_from_velocity():
    """Test the output of ustar_from_velocity."""
    ustar = ustar_from_velocity(20, 45, 3, -0.00204)
    assert isclose(ustar, 51.619, abs_tol=1e-3)


# def test_velocity_from_ustar():
#     """Test the output of velocity_from_ustar."""
#     u = velocity_from_ustar(3, 25, 5)
#     assert u == 11.776374969030002


class TestCalcMoninObukhovLength:

    def test_outputs_correct_value(self):
        out = calc_monin_obukhov_length(
            ustar=0.5386,
            Hd=-28.25,
            P=102.1,
            Tk=0.428 + 273.15,
        )
        assert isclose(1 / out, -0.002034, abs_tol=1e-3)

    def test_outputs_correct_value_last_row(self):
        out = calc_monin_obukhov_length(
            ustar=0.573,
            Hd=-13.5,
            P=99.2,
            Tk=4.76 + 273.15,
        )
        assert isclose(1 / out, -0.00083, abs_tol=1e-5)


class TestCalcUstarRefAndL:

    @pytest.mark.parametrize(['u', 'Hd', 'P', 'Tk', 'ustar_ref_target'], [
        # Below are EMEP grid%ustar values
        # (6.076,-28.25, 102.1, 0.428+273.15, 0.5386),
        # (6.09,-11.99, 102.1, 0.351+273.15, 0.4796),
        # (3.696,-16.04, 102.6, -0.11+273.15, 0.1663),
        (6.076, -28.25, 102.1, 0.428 + 273.15, 0.75378),
        (3.696, -16.04, 102.6, -0.11 + 273.15, 0.4699),

    ])
    def test_should_output_correct_values(self, u, Hd, P, Tk, ustar_ref_target):
        h_u = 20
        canopy_d = 0.78
        canopy_z0 = 0.07
        u_d = h_u * canopy_d
        u_z0 = min(1.0, h_u * canopy_z0)

        ustar_ref, L = calc_ustar_and_L(
            u=u,
            Hd=Hd,
            P=P,
            Tk=Tk,
            z_u=45,
            u_d=u_d,
            u_z0=u_z0,
            MAX_ITERATIONS=20,
            initial_L=1 / 0.000000000001
        )
        assert isclose(ustar_ref, ustar_ref_target, abs_tol=1e-3)
        # assert isclose(L, 1/-0.002034, abs_tol=1e-3) # Tested with monin obukhov length tests

    @pytest.mark.parametrize('u', [1, 6, 20])
    @pytest.mark.parametrize('h', [1, 5, 20])
    @pytest.mark.parametrize('Hd', [-28.1, 30])
    def test_should_converge_on_values_multiple(self, h, u, Hd):
        h_u = h
        canopy_d = 0.78
        canopy_z0 = 0.07
        u_d = h_u * canopy_d
        u_z0 = min(1.0, h_u * canopy_z0)

        ustar_ref_19, L = calc_ustar_and_L(
            u=u,
            Hd=Hd,
            P=102.1,
            Tk=0.428 + 273.15,
            z_u=45,
            u_d=u_d,
            u_z0=u_z0,
            MAX_ITERATIONS=19,
            initial_L=1 / 0.000000000001
        )
        ustar_ref_20, L = calc_ustar_and_L(
            u=u,
            Hd=Hd,
            P=102.1,
            Tk=0.428 + 273.15,
            z_u=45,
            u_d=u_d,
            u_z0=u_z0,
            MAX_ITERATIONS=20,
            initial_L=1 / 0.000000000001
        )
        assert isclose(ustar_ref_19, ustar_ref_20, abs_tol=1e-2)


class TestCalcWindspeed:

    # TODO: Update these outputs to match UI
    def test_calc_windspeed_parameters(self):
        """Test the output of calc_windspeed_parameters."""
        out = calc_windspeed(
            h=1, u=20, L=0.1, z_u=10, u_z0=1.0, d=0.78, z0=0.07, min_windspeed=0.1)
        assert isclose(out.ustar, 0.0194, abs_tol=1e-3)
        assert isclose(out.u_i, 116.881, abs_tol=1e-3)
        assert isclose(out.micro_u, 0.410, abs_tol=1e-3)

    def test_matches_emep_53_29_row_10(self):
        h_u = 20
        h = 20
        canopy_d = 0.78
        canopy_z0 = 0.07
        z_u = 45
        # Forest specific
        u_z0 = min(1.0, h_u * canopy_z0)
        d = h * canopy_d
        z0 = min(1.0, h * canopy_z0)
        L = 1 / -0.0020217653363943
        u = 6.116
        assert isclose(d, 15.599, abs_tol=1e-2)
        assert z0 == 1.0
        assert u_z0 == 1.0
        out = calc_windspeed(
            h=h,
            u=u,
            L=L,
            z_u=z_u,
            u_z0=u_z0,
            d=d,
            z0=z0,
            izr=45,
        )
        assert isclose(out.u_i, 6.116, abs_tol=1e-3)
        assert isclose(out.ustar, 0.7873, abs_tol=1e-2)  # EMEP OUTPUT
        assert isclose(out.micro_u, 2.795, abs_tol=1e-1)

    def test_calc_windspeed_parameters_compare_ui_row_253(self):
        """Test the output of calc_windspeed_parameters matches ui."""
        out = calc_windspeed(
            h=1, u=0.73, L=0.1, z_u=10, u_z0=1.0, d=0.78, z0=0.07, min_windspeed=0.1)
        assert isclose(out.ustar, 0.0007, abs_tol=1e-4)
        assert isclose(out.u_i, 4.266, abs_tol=1e-3)
        assert isclose(out.micro_u, 0.1, abs_tol=1e-3)

    def test_calc_with_simple_method(self):
        h_u = 20
        h = 20
        canopy_d = 0.78
        canopy_z0 = 0.07
        z_u = 45
        # Forest specific
        u_z0 = min(1.0, h_u * canopy_z0)
        d = h * canopy_d
        u_d = h_u * canopy_d
        z0 = min(1.0, h * canopy_z0)
        u = 6.116
        assert isclose(d, 15.599, abs_tol=1e-2)
        assert z0 == 1.0
        assert u_z0 == 1.0
        out = calc_windspeed(
            h=h,
            u=u,
            L=None,  # Invokes simple method
            ustar_ref=0.7537817236092231,
            z_u=z_u,
            u_z0=u_z0,
            u_d=u_d,
            d=d,
            z0=z0,
            izr=45,
        )
        # Lower precision for simple method
        assert isclose(out.u_i, 6.116, abs_tol=1e-1)
        assert isclose(out.ustar, 0.7873, abs_tol=1e-1)  # EMEP OUTPUT
        assert isclose(out.micro_u, 2.795, abs_tol=1e-1)

    def test_calc_with_simple_method_measured_height_equal_canopy_height(self):
        h_u = 20
        h = 20
        canopy_d = 0.78
        canopy_z0 = 0.07
        z_u = 20
        # Forest specific
        u_z0 = min(1.0, h_u * canopy_z0)
        d = h * canopy_d
        u_d = h_u * canopy_d
        z0 = min(1.0, h * canopy_z0)
        u = 6.116
        izr = 45
        assert isclose(d, 15.599, abs_tol=1e-2)
        assert z0 == 1.0
        assert u_z0 == 1.0
        ustar_ref = ustar_from_velocity_simple(u, z_u - u_d, u_z0)
        out = calc_windspeed(
            h=h,
            u=u,
            L=None,  # Invokes simple method
            ustar_ref=ustar_ref,
            z_u=z_u,
            u_z0=u_z0,
            u_d=u_d,
            d=d,
            z0=z0,
            izr=izr,
        )
        # Lower precision for simple method
        # assert isclose(out.u_i, 6.116, abs_tol=1e-1)
        assert isclose(out.ustar, ustar_ref, abs_tol=1e-3)
        assert isclose(out.micro_u, u, abs_tol=1e-1)

    # Should not cause errors
    # def test_calc_windspeed_parameters_error(self):
    #     """Should return a math domain error if invalid h and z values.

    #     i.e if h is too much higher than z_u
    #     """
    #     with pytest.raises(ValueError) as e:
    #         calc_windspeed_parameters(
    #             h=10, u=0.52, L=0.1)
    #     assert 'math domain error' in str(e)

    #     with pytest.raises(ValueError) as e:
    #         calc_windspeed_parameters(
    #             h=2, u=0.52, L=0.1)
    #     assert 'math domain error' in str(e)

    #     with pytest.raises(ValueError) as e:
    #         calc_windspeed_parameters(
    #             h=2, u=0.52, L=0.1, z_u=1, u_z0=)
    #     assert 'math domain error' in str(e)

    #     with pytest.raises(ValueError) as e:
    #         z = 1
    #         h = z / CANOPY_D
    #         calc_windspeed_parameters(
    #             h=h, u=0.52, L=0.1, z_u=1, u_z0=)
    #     assert 'math domain error' in str(e)

    #     # with pytest.raises(ValueError) as e:
    #     #     calc_windspeed_parameters(
    #     #         h=2, u=0.52, o_top_chamber=False, loc_z_u=10, loc_h_u=None)
    #     # assert 'math domain error' in str(e)


def test_calc_multi_layer_windspeed():
    """Test the output of calc_multi_layer_windspeed."""
    u_at_canopy_top = 21.1
    z_u = calc_layer_windspeed(h=10, w=0.02, SAI=3.0, u_at_canopy_top=u_at_canopy_top, z=5)
    assert isclose(z_u, 2.181, abs_tol=1e-3)

    # When using OTC windspeeds are equal throughout canopy
    z_u_open_top = calc_layer_windspeed(
        h=10, w=0.02, SAI=3.0, u_at_canopy_top=u_at_canopy_top, z=5, o_top_chamber=True)
    z_u_open_top_b = calc_layer_windspeed(
        h=15, w=0.04, SAI=5.0, u_at_canopy_top=u_at_canopy_top, z=5, o_top_chamber=True)
    assert z_u_open_top == z_u_open_top_b == u_at_canopy_top


class TestIntegrated:

    def test_should_output_same_wind_value_if_measured_at_canopy(self):
        for u in range(1, 5):
            for h in range(1, 20):
                try:
                    h_u = h  # Canopy height for windspeed measurements
                    # u_h = h_u # Canopy height for windspeed measurements
                    izr = 20
                    z_u = h  # Windspeed measured height
                    canopy_d = 0.78  # 0.78
                    canopy_z0 = 0.07  # 0.07
                    # u_z0 is the height at which u=0[m] I.e. 7% of canopy height
                    u_z0 = min(1.0, h_u * canopy_z0)
                    # d is height 78% of height of canopy ( Canopy displacement height)
                    h_d = h * canopy_d
                    # u_d is 78% of wind measured height
                    u_d = h_u * canopy_d
                    z0 = min(1.0, h * canopy_z0)
                    ustar_ref = ustar_from_velocity_simple(u, (z_u - u_d), u_z0)

                    out = calc_windspeed(
                        h=h,
                        u=u,
                        L=None,  # Invokes simple method
                        ustar_ref=ustar_ref,
                        z_u=z_u,
                        u_z0=u_z0,
                        u_d=u_d,
                        d=h_d,
                        z0=z0,
                        izr=izr,
                    )
                    # Lower precision for simple method
                    assert isclose(out.micro_u, u, abs_tol=1e-1)

                except Exception as e:
                    print(e)
                    print(f"Failed with {h}")
                    raise e

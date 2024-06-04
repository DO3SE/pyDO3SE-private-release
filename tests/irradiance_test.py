"""Tests for met modules helpers."""

from math import isclose
from typing import NamedTuple
import pytest


from do3se_met.irradiance import (
    MLMC_sunlit_LAI,
    calc_Idrctt_Idfuse,
    calc_Ir_beam_sun,
    calc_Ir_diffuse,
    calc_Ir_scattered_b,
    calc_PAR_shade,
    calc_PAR_sun,
    calc_PAR_sun_shade,
    calc_PAR_sun_shade_UI,
    calc_PAR_sun_shade_farq,
    calc_PAR_sun_shade_farq_b,
    calc_beam_irradiance_horiz,
    calc_beam_irradiance_uad,
    calc_diffuse_irradiance_refl,
    calc_is_daylight,
    calc_photoperiod,
    calc_photoperiod_factor,
    calc_radiation,
    get_net_radiation,
    get_parSunShade_from_par,
    sunlit_LAI,
    calc_net_radiation,
)


def test_calc_net_radiation():
    out = calc_net_radiation(
        lat=52.2,
        lon=-1.12,
        elev=20.0,
        albedo=0.2,
        dd=135,
        hr=1,
        sinB=0.3,
        R=900.0,
        Ts_C=20.1,
        eact=0.8,
    )
    assert isclose(out, 2.2675, abs_tol=1e-3)


def test_calc_daylight():
    """Test the output of calc_daylight."""
    is_daylight = calc_is_daylight(30)
    assert is_daylight is False
    is_daylight = calc_is_daylight(80)
    assert is_daylight is True


@pytest.mark.parametrize(['dd', 'lat', 'out'], [
    (96, 52, 14.392),
    (1, 52, 8.988),
    (0, 52, 8.969),
    (0,0, 12.736),
    (364,52, 8.951),
    (365,52, 8.969),
    (366,52, 8.988),
    (147,63.625, 24), # photoperiod is > 100% I.e 24 hours
    (147,-63.625, 24), # photoperiod is > 100% I.e 24 hours
])
def test_calc_photoperiod(dd, lat, out):
    """Test the output of calc_photoperiod."""
    pr = calc_photoperiod(dd, lat)
    assert isclose(pr, out, abs_tol=1e-3)


def test_calc_photoperiod_factor():
    """Test the output of calc_photoperiod_factor."""
    prf = calc_photoperiod_factor(
        photoperiod=12.684,
        PID=40,
    )
    assert isclose(prf, 0.7859, abs_tol=1e-3)


def test_calc_PAR_sun_shade():
    """Test the output of calc_PAR_sun_shade."""
    PARsun, PARshade = calc_PAR_sun_shade(50, 44, 1.3, 0.5, 1.5)
    assert isclose(PARsun, 38.3252, rel_tol=1e-4)
    assert isclose(PARshade, 22.9406, rel_tol=1e-4)
    PARsun, PARshade = calc_PAR_sun_shade(99, 20, 1.8, 0.4, 3.2)
    assert isclose(PARsun, 24.1218, rel_tol=1e-4)
    assert isclose(PARshade, 6.5218, rel_tol=1e-4)


def test_calc_Idrctt_Idfuse():
    """Test the output of calc_Idrctt_Idfuse."""
    Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(1.8, 1.9, PAR=10.0)
    assert isclose(Idrctt, 0.1402, abs_tol=1e-4)
    assert isclose(Idfuse, 9.8597, abs_tol=1e-4)

    Idrctt, Idfuse, _ = calc_Idrctt_Idfuse(0.384, 91.2, PAR=300)
    assert isclose(Idrctt, 204.9826, abs_tol=1e-4)
    assert isclose(Idfuse, 95.0173, abs_tol=1e-4)


def test_calc_Idrctt_Idfuse_from_cloud_frac():
    """Test the output of calc_Idrctt_Idfuse when using cloudFrac input."""
    Idrctt, Idfuse, PAR = calc_Idrctt_Idfuse(
        cloudFrac=0.60194445, sinB=0.656637609004974, P=99.5498733520507)
    assert isclose(Idrctt, 224.700506994, abs_tol=1e-4)
    assert isclose(Idfuse, 102.700037123, abs_tol=1e-4)
    assert isclose(PAR, 327.40054411, abs_tol=1e-4)


class TestCalcRadiation():
    ExpectedOut = NamedTuple(
        'ExpectedOut',
        [
            ('PAR', float),
            ('PPFD', float),
            ('Idrctt', float),
            ('Idfuse', float),
            ('R', float),
            ('Rn', float),
        ])
    expected_outputs = ExpectedOut(
        350.0,
        1599.5,
        128.574,
        221.426,
        777.777,
        2.7997,
    )

    def matches_expected_out(self, out):
        assert isclose(out.PAR, self.expected_outputs.PAR, abs_tol=1e-1)
        assert isclose(out.PPFD, self.expected_outputs.PPFD, abs_tol=1e-1)
        assert isclose(out.Idrctt, self.expected_outputs.Idrctt, abs_tol=1e-1)
        assert isclose(out.Idfuse, self.expected_outputs.Idfuse, abs_tol=1e-1)
        assert isclose(out.R, self.expected_outputs.R, abs_tol=1e-1)
        assert isclose(out.Rn, self.expected_outputs.Rn, abs_tol=1e-1)

    def test_calc_radiation_from_ppfd(self):
        """Test the output of calc_radiation."""

        P = 1.9  # Always need P!
        sinB = 1.3  # Always need P!
        out = calc_radiation(
            PAR_in=None,
            Idrctt_in=None,
            Idfuse_in=None,
            PPFD_in=self.expected_outputs.PPFD,
            R_in=None,
            sinB=sinB,
            P=P,
            cloudfrac=None,
        )
        self.matches_expected_out(out)

    def test_calculate_from_Idrc_in(self):
        out = calc_radiation(
            PAR_in=None,
            Idrctt_in=self.expected_outputs.Idrctt,
            Idfuse_in=self.expected_outputs.Idfuse,
            PPFD_in=None,
            R_in=None,
            sinB=None,
            P=None,
            cloudfrac=None,
        )
        self.matches_expected_out(out)

    def test_calculate_from_PAR_inputs(self):
        P = 1.9  # Always need P!
        sinB = 1.3  # Always need sinB!
        out = calc_radiation(
            PAR_in=self.expected_outputs.PAR,
            Idrctt_in=None,
            Idfuse_in=None,
            PPFD_in=None,
            R_in=None,
            sinB=sinB,
            P=P,
            cloudfrac=None,
        )
        self.matches_expected_out(out)

    def test_calc_from_cloud_frac(self):
        P = 1.9  # Always need P
        sinB = 1.3  # Always need sinB
        cloudfrac = 0.881985
        out = calc_radiation(
            PAR_in=None,
            Idrctt_in=None,
            Idfuse_in=None,
            PPFD_in=None,
            R_in=None,
            sinB=sinB,
            P=P,
            cloudfrac=cloudfrac,
        )
        self.matches_expected_out(out)

    def test_calc_radiation_from_Rn(self):
        """Test the output of calc_radiation."""

        P = 1.9  # Always need P!
        sinB = 1.3  # Always need P!
        out = calc_radiation(
            PAR_in=None,
            Idrctt_in=None,
            Idfuse_in=None,
            PPFD_in=None,
            R_in=None,
            Rn_in=self.expected_outputs.Rn,
            sinB=sinB,
            P=P,
            cloudfrac=None,
        )
        self.matches_expected_out(out)

    def test_throws_error_if_not_enough_inputs(self):
        with pytest.raises(ValueError):
            calc_radiation(None, None, None, None)

    def test_match_ui(self):
        P = 102  # Always need P!
        sinB = 0.91134  # Always need P!
        cloudfrac = 0.98600
        out = calc_radiation(
            PAR_in=None,
            Idrctt_in=None,
            Idfuse_in=None,
            PPFD_in=None,
            R_in=None,
            sinB=sinB,
            P=P,
            cloudfrac=cloudfrac,
        )
        assert isclose(out.PAR, 143.1154, abs_tol=1e-3)
        assert isclose(out.Idrctt, 10.521, abs_tol=1e-3)
        assert isclose(out.Idfuse, 132.594, abs_tol=1e-3)


def test_get_net_radiation():
    """Test the output of get_net_radiation."""
    Rn = get_net_radiation(1.1)
    assert Rn == 1.1
    Rn = get_net_radiation(
        lat=52.2,
        lon=-1.12,
        elev=20.0,
        albedo=0.2,
        dd=135,
        hr=1,
        sinB=0.3,
        R=900.0,
        Ts_C=20.1,
        eact=0.8,
    )
    # TODO: why is this always 0.0?
    assert isclose(Rn, 2.2675, abs_tol=1e-3)


def test_sunlit_LAI():
    """Test the output of sunlit_LAI."""
    out = sunlit_LAI(
        LAI=2.1,
        sinB=0.7,
    )
    assert isclose(out, 1.087617, abs_tol=1e-3)


def test_MLMC_sunlit_LAI_single_layer():
    """Test the output of MLMC_sunlit_LAI."""
    out = MLMC_sunlit_LAI(
        nL=1,
        nLC=1,
        LAI=[[1.0]],
        sinB=0.3,
    )
    assert isclose(out[0][0], 0.48667, abs_tol=1e-5)


def test_MLMC_sunlit_LAI_single_component():
    """Test the output of MLMC_sunlit_LAI."""
    out = MLMC_sunlit_LAI(
        nL=3,
        nLC=1,
        LAI=[[1.0], [1.0], [1.0]],
        sinB=0.3,
    )
    assert isclose(out[0][0], 0.48667, abs_tol=1e-5)
    assert isclose(out[1][0], 0.091921, abs_tol=1e-5)
    assert isclose(out[2][0], 0.017361, abs_tol=1e-5)


def test_MLMC_sunlit_LAI():
    """Test the output of MLMC_sunlit_LAI."""
    out = MLMC_sunlit_LAI(
        nL=2,
        nLC=3,
        LAI=[[2.3, 2.4, 2.5], [1.4, 1.5, 1.6]],
        sinB=0.7,
    )
    assert isclose(out[0][0], 0.19330, abs_tol=1e-5)
    assert isclose(out[1][0], 0.00174, abs_tol=1e-5)


# Farquhar 1997 leaf Irradiance Equations

def test_calc_beam_irradiance_horiz():
    """Test the output of calc_beam_irradiance_horiz."""
    out = calc_beam_irradiance_horiz(
        sigma=0.15,
    )

    assert isclose(out, 0.04060, abs_tol=1e-3)


def test_calc_beam_irradiance_uad():
    """Test the output of calc_beam_irradiance_uad."""
    out = calc_beam_irradiance_uad(
        P_h=0.04,  # value for horizontal leaves
        k_b=0.5 / 0.87,  # 0.5/sinb
    )
    assert isclose(out, 0.0287749, abs_tol=1e-4)
    # assert isclose(out, -0.0296, abs_tol=1e-4)


def test_calc_diffuse_irradiance_refl():
    """Test the output of calc_diffuse_irradiance_refl."""
    out = calc_diffuse_irradiance_refl(
        sinB=0.87,
        P_h=0.04,
        Ir_dfuse_0=300,
    )

    assert isclose(out, 0.0074, abs_tol=1e-3)


def test_calc_Ir_scattered_b():
    """Test the output of calc_Ir_scattered_b."""
    out = calc_Ir_scattered_b(
        P_cb=0.1,
        Ir_beam_0=300,
        LAI_c=0.3,
        k_b=0.25,
        k_b_alt=0.46 / 0.87,
        sigma=0.15,
    )
    assert isclose(out, 73.299, abs_tol=1e-3)


def test_calc_Ir_beam_sun():
    """Test the output of calc_Ir_beam_sun."""
    out = calc_Ir_beam_sun(
        sinB=0.87,
        cosA=0.5,
        Ir_beam_0=300,
        sigma=0.15,
    )

    assert isclose(out, 146.551, abs_tol=1e-3)


def test_calc_Ir_diffuse():
    """Test the output of calc_Ir_diffuse."""
    out = calc_Ir_diffuse(
        P_cd=0.1,
        Ir_dfuse_0=300,
        LAI_c=0.3,
        k_d_alt=0.719,
    )
    assert isclose(out, 156.464, abs_tol=1e-3)


def test_calc_PAR_shade():
    """Test the output of calc_PAR_shade."""
    out = calc_PAR_shade(
        Ir_diffuse=400,
        Ir_scattered_b=300,
    )

    assert isclose(out, 700, abs_tol=1e-3)


def test_calc_PAR_sun():
    """Test the output of calc_PAR_sun."""
    out = calc_PAR_sun(
        PAR_shade=400,
        Ir_beam_sun=300,
    )
    assert isclose(out, 700, abs_tol=1e-3)


def test_calc_PAR_sun_shade_farq():
    """Test the output of calc_PAR_sun_shade_farq."""
    out_sun, out_shade = calc_PAR_sun_shade_farq(
        Ir_beam_0=300.0,
        Ir_dfuse_0=300.0,
        sinB=0.87,
        cosA=0.5,
        LAI_c=0.3,
        sigma=0.15,
    )

    assert isclose(out_sun, 325.390321, abs_tol=1e-3)
    assert isclose(out_shade, 178.83859, abs_tol=1e-3)


def test_calc_PAR_sun_shade_farq_b():
    """Test the output of calc_PAR_sun_shade_farq_b.

    Should be exactly the same output as test_calc_PAR_sun_shade_farq
    """
    out_sun, out_shade = calc_PAR_sun_shade_farq_b(
        Ir_beam_0=300.0,
        Ir_dfuse_0=300.0,
        sinB=0.87,
        cosA=0.5,
        LAI_c=0.3,
        sigma=0.15,
    )

    assert isclose(out_sun, 325.390321, abs_tol=1e-3)
    assert isclose(out_shade, 178.83859, abs_tol=1e-3)


def test_calc_PAR_sun_shade_farq_b_irrad_0():
    """Test the output of calc_PAR_sun_shade_farq_b when irrad is 0."""
    out_sun, out_shade = calc_PAR_sun_shade_farq_b(
        Ir_beam_0=0.0,
        Ir_dfuse_0=0.0,
        sinB=0.87,
        cosA=0.5,
        LAI_c=0.3,
        sigma=0.15,
    )

    assert isclose(out_sun, 0.0, abs_tol=1e-3)
    assert isclose(out_shade, 0.0, abs_tol=1e-3)


def test_get_parSunShade_from_par():
    PARsun, PARshade = get_parSunShade_from_par(
        PAR=200,
        sinB=0.87,
        cosA=0.5,
        P=99.54,
        LAI=1,
    )
    assert isclose(PARsun, 116.91231679798183, abs_tol=1e-3)
    assert isclose(PARshade, 99.70992235136127, abs_tol=1e-3)


def test_calc_PAR_sun_shade_UI():
    PARsun, PARshade = calc_PAR_sun_shade_UI(
        PAR=200,
        sinB=0.87,
        cosA=0.5,
        P=99.54,
        LAI=1,
    )

    assert isclose(PARsun, 116.91231679798183, abs_tol=1e-3)
    assert isclose(PARshade, 99.70992235136127, abs_tol=1e-3)

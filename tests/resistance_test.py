import pytest
from math import isclose, inf
from decimal import Decimal as D

from do3se_met.resistance import (
    calc_PsiH,
    calc_displacement_and_roughness_parameters,
    calc_PsiM,
    calc_Ra_with_heat_flux,
    calc_Ra_with_heat_flux_old,
    calc_Ra_simple,
    calc_Rb,
    calc_leaf_rb,
    calc_leaf_gb,
    calc_Rtotal,
    calc_Rsto,
    calc_Rext,
    calc_Rinc,
    calc_Rgs,
    calc_Rsur,
    calc_Rsur_multilayer,
    calc_deposition_velocity,
    calc_resistance_model,
    Resistance_Model,
    R_INF,
)


class TestCalcDisplacementAndRoughnessParameters:

    @pytest.mark.parametrize(['h', 'is_forest', 'd', 'z0'], [
        (20, True, 15.600, 1.0),
        (20, False, 14.0, 2.0),
    ])
    def test_outputs_correct_values(self, h, is_forest, d, z0):
        d_out, z0_out = calc_displacement_and_roughness_parameters(h, is_forest)
        assert isclose(d_out, d, abs_tol=1e-3)
        assert isclose(z0_out, z0, abs_tol=1e-3)


class TestCalcPsiH:

    @pytest.mark.parametrize(['length', 'e_out'], [
        (0, 0),
        (0.1, -0.5),
        (1, -5.0),
        (10, -50.0),
    ])
    def test_outputs(self, length, e_out):
        out = calc_PsiH(length)
        assert isclose(out, e_out, abs_tol=1e-3)


class TestCalcPsiM:

    @pytest.mark.parametrize(['length', 'e_out'], [
        (0, 0),
        (D(0.1), -0.5),
        (D(1), -5.0),
        (D(10), -50.0),
    ])
    def test_correct_outputs(self, length, e_out):
        out = calc_PsiM(length)
        assert isclose(out, e_out, abs_tol=1e-3)


class TestCalcRaWithHeatFlux:

    @pytest.mark.parametrize(['ustar', 'z1', 'z2', 'L', 'e_out'], [
        (0.3, 10, 50, -0.001, 0.0224),
        (0.1, 40, 50, -0.001, 0.0064),
    ])
    def test_correct_outputs(self, ustar, z1, z2, L, e_out):
        out = calc_Ra_with_heat_flux(ustar, z1, z2, L)
        assert isclose(out, e_out, abs_tol=1e-3)


# class TestCalcRa_with_heat_flux_old:

#     @pytest.mark.parametrize([], [
#         (),
#     ])
#     def test_correct_outputs(self):
#         out = calc_Ra_with_heat_flux_old()
#         assert isclose(out, 99.9, abs_tol=1e-3)


class TestCalcRaSimple:

    @pytest.mark.parametrize(['ustar', 'z1', 'z2', 'd', 'e_out'], [
        (0.3, 10, 50, 7.8, 24.015),
        # Ra simple should be similar to Ra_heat flux below:
        # (0.3, 10,50,7.8, 0.0224),
        # (0.1, 40,50,0.2, 0.0064),
    ])
    def test_correct_outputs(self, ustar, z1, z2, d, e_out):
        out = calc_Ra_simple(ustar, z1 - d, z2, d)
        assert isclose(out, e_out, abs_tol=1e-3)


class TestCalcRb:

    @pytest.mark.parametrize(['ustar', 'diff', 'e_out'], [
        (0.3, 0.1, 0.05714),
    ])
    def test_correct_outputs(self, ustar, diff, e_out):
        out = calc_Rb(ustar, diff)
        assert isclose(out, e_out, abs_tol=1e-3)


class TestCalcDepositionVelocity:

    @pytest.mark.parametrize(['Ra_c', 'Rtotal', 'e_out'], [
        # TODO: Get correct inputs here
        (99.9, 99.9, 0.00500),
    ])
    def test_correct_outputs(self, Ra_c, Rtotal, e_out):
        out = calc_deposition_velocity(Ra_c, Rtotal)
        assert isclose(out, e_out, abs_tol=1e-3)


def test_init_resistance_model():
    rmodel = Resistance_Model(3)
    assert isinstance(rmodel, Resistance_Model)
    assert rmodel.nL == 3


# def test_calc_Ra_with_heat_flux():
#     Ra = calc_Ra_with_heat_flux(
#         Ts_C=24.4,
#         Hd=0.1,
#         P=101.0,
#         ustar=24.4,
#         d=1.23,
#         zo=20,
#     )
#     assert isclose(Ra, 0.0891, abs_tol=1e-4)


# def test_calc_ra_simple():
#     Ra = calc_Ra_simple(1.1, 2, 5, 1)
#     assert isclose(Ra, 3.0738, abs_tol=1e-4)


# def test_calc_Rb():
#     Rb = calc_Rb(1.1, 1.2)
#     assert isclose(Rb, 0.00297, abs_tol=1e-4)


class TestLeafGb:

    def test_calc_leaf_gb(self):
        leaf_gb = calc_leaf_gb(
            G=1.23,
            Lm=3.45,
            u=3.33,
        )
        assert isclose(leaf_gb, 2.4168, abs_tol=1e-4)

    def test_calc_leaf_gb_zero_wind(self):
        leaf_gb = calc_leaf_gb(
            G=1.23,
            Lm=3.45,
            u=0,
        )
        assert isclose(leaf_gb, 0.0, abs_tol=1e-4)


class TestCalcLeafRb:
    def test_calc_leaf_rb(self):
        leaf_rb = calc_leaf_rb(
            gb=1.23,
        )
        assert isclose(leaf_rb, 33.3333, abs_tol=1e-4)

    def test_calc_leaf_rb_low_gb(self):
        leaf_rb = calc_leaf_rb(
            gb=0.000001,
        )
        assert isclose(leaf_rb, 41000000.0, abs_tol=1e-4)

    def test_calc_leaf_rb_zero_gb(self):
        leaf_rb = calc_leaf_rb(
            gb=0.0,
        )
        assert isclose(leaf_rb, inf, abs_tol=1e-4)


def test_calc_Rinc():
    Rinc = calc_Rinc(1.5, 10, 1.1)
    assert isclose(Rinc, 190.9090, abs_tol=1e-4)


def test_calc_Rext():
    Rext = calc_Rext(2500, 1.2)
    assert isclose(Rext, 2083.3333, abs_tol=1e-4)


def test_calc_Rsto():
    Rsto = calc_Rsto(1.2)
    assert isclose(Rsto, 34166.6666, abs_tol=1e-4)


def test_calc_Rsto_compare_with_ui_row_2816():
    Rsto = calc_Rsto(169.590332031)
    assert isclose(Rsto, 241.759063721, abs_tol=1e-4)


def test_calc_Rsur_multilayer():
    Rsur = calc_Rsur_multilayer(
        nL=3,
        Rb=0.3,
        Rsto=[100, 200, 300],
        Rext=[100, 200, 300],
        LAI=[0, 0.2, 0.3],
        SAI=[0.0, 0.0, 0.5])
    assert isclose(Rsur[0], R_INF, abs_tol=1e9)
    assert isclose(Rsur[1], 100.3, abs_tol=1e-9)
    assert isclose(Rsur[2], 150.3, abs_tol=1e-9)


def test_calc_Rsur_match_ui():
    """Check Rsur calc matches ui."""
    Rsur = calc_Rsur(
        Rb=31.2906856537,
        Rsto_c=100000,
        Rext=2500 / 3,
        Rinc=216.424575806,
        Rsoil=200,
        LAI=3,
        SAI=3,
    )

    assert isclose(Rsur, 276.901275635, abs_tol=1e-4)


class TestCalcRtotal:

    def test_const_rsur(self):
        nL = 4
        Rtotal = calc_Rtotal(
            nL=nL,
            Rsur=[100 for _ in range(nL)],
            Rinc=[3 for _ in range(nL)],
            Rgs=200,
        )
        # Resistance should increase as you descend levels if other resistances are the same
        assert all([a > b for a, b in zip(Rtotal, Rtotal[1:])])

        assert isclose(Rtotal[0], 66.997, abs_tol=1e-3)
        assert isclose(Rtotal[1], 41.175, abs_tol=1e-3)
        assert isclose(Rtotal[2], 30.640, abs_tol=1e-3)
        assert isclose(Rtotal[3], 25.172, abs_tol=1e-3)

    def test_const_rsur_low_rinc(self):
        Rtotal = calc_Rtotal(
            nL=3,
            Rsur=[100, 100, 100],
            Rinc=[0.1, 0.1, 0.1],
            Rgs=333,
        )
        # Resistance should increase as you descend levels if other resistances are the same
        assert all([a > b for a, b in zip(Rtotal, Rtotal[1:])])

        assert isclose(Rtotal[2], 30.365, abs_tol=1e-3)
        assert isclose(Rtotal[1], 43.506, abs_tol=1e-3)
        assert isclose(Rtotal[0], 76.910, abs_tol=1e-3)

    def test(self):
        Rtotal = calc_Rtotal(
            nL=3,
            Rsur=[300, 200, 100],
            Rinc=[103, 102, 101],
            Rgs=333,
        )

        assert isclose(Rtotal[2], 68.515, abs_tol=1e-3)
        assert isclose(Rtotal[1], 116.617, abs_tol=1e-3)
        assert isclose(Rtotal[0], 177.717, abs_tol=1e-3)

    def test_vary_rgs(self):
        Rtotal = calc_Rtotal(
            nL=3,
            Rsur=[10, 10, 10],
            Rinc=[10, 10, 10],
            Rgs=100,
        )

        assert isclose(Rtotal[2], 6.236, abs_tol=1e-3)
        assert isclose(Rtotal[1], 6.571, abs_tol=1e-3)
        assert isclose(Rtotal[0], 9.166, abs_tol=1e-3)
        # Rtotal for layer below layer[0] is the soil therefore will == Rgs


def test_calc_deposition_velocity():
    d_vel = calc_deposition_velocity(1.0, 1.1)
    assert isclose(d_vel, 0.4761, abs_tol=1e-4)


def test_calc_Rgs():
    Rgs = calc_Rgs(Rsoil=200)
    assert isclose(Rgs, 200, abs_tol=1e-4)


def test_calc_Rgs_with_snow_cover():
    Rgs = calc_Rgs(Rsoil=200, snow_depth=0.5, Ts_C=0.5)
    assert isclose(Rgs, 160.8435465029701, abs_tol=1e-4)


class TestCalcResistanceModel:
    def test_works_without_error(self, snapshot):
        rmodel = calc_resistance_model(
            nL=2,
            nLC=2,
            ustar=1.1,
            canopy_height=3.3,
            SAI_values=[[1, 2], [2, 3]],
            LAI_values=[[1, 2], [2, 3]],
            mean_gsto_values=[[1, 2], [2, 3]],
            Rsoil=123,
            rsur_calc_method="multi_layer",
            ra_calc_method="simple",
        )
        assert isinstance(rmodel, Resistance_Model)
        snapshot.assert_match(rmodel)

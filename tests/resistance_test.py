import pytest
from math import isclose
from decimal import Decimal as D

from do3se_met.resistance import (
    calc_PsiH,
    calc_displacement_and_roughness_parameters,
    calc_PsiM,
    calc_Ra_with_heat_flux,
    calc_Ra_with_heat_flux_old,
    calc_Ra_simple,
    calc_Rb,
    calc_deposition_velocity,
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
        out = calc_Ra_simple(ustar, z1, z2, d)
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

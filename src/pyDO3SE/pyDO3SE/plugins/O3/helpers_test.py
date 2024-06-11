import pytest
from math import inf, isclose
from .helpers import (
    calc_FO3_eff,
    calc_O3up_accumulation,
    calc_OT_acc,
    calc_POD,
    calc_fst,
    O3_ppb_to_nmol,
    calc_fst_leaf,
    stomatal_flux_rate,
)


class TestCalcFst:

    def test_calc_fst(self):
        fst = calc_fst(
            Gsto_l=495.529266357,
            Rb_l=10,
            Rsto_l=82.7398147583,
            Rext=2500,
            O3_nmol_m3=11007.807805097667,
        )
        assert isclose(fst, 118.273517, abs_tol=1e-3)

    def test_calc_fst_match_do3se_ui_row_3655(self):
        fst = calc_fst(
            Gsto_l=8.73888683319,
            Rb_l=35,
            Rsto_l=4691.67285156,
            Rext=2500,
            O3_nmol_m3=1087.22546387,
        )
        assert isclose(fst, 0.226157203317, abs_tol=1e-3)


class TestCalcFstLeaf:

    @pytest.mark.parametrize(['Gsto_l', 'expected_output'], [
        ([99.9, 99.9, 99.9], 0.003113),
    ])
    def test_calculates_fst(self, Gsto_l, expected_output):
        nL = 3
        Rext = [300 for _ in range(nL)]
        O3_nmol_m3 = [3, 2, 1]
        fLAI = [0.1, 0.5, 0.4]
        Rsto_l = [41000 / Gsto_l[iL] for iL in range(nL)]
        Rb_l = [3, 3, 3]
        out = calc_fst_leaf(
            Gsto_l,
            Rb_l,
            Rsto_l,
            Rext,
            O3_nmol_m3,
            fLAI,
            nL,
        )

        assert isclose(out, expected_output, abs_tol=1e-3)

    def test_handles_infinite_Rb(self):
        nL = 1
        Rext = [300 for _ in range(nL)]
        O3_nmol_m3 = [3]
        fLAI = [1]
        Rsto_l = [41000 / 1]
        Rb_l = [inf]
        Gsto_l = [1]
        out = calc_fst_leaf(
            Gsto_l,
            Rb_l,
            Rsto_l,
            Rext,
            O3_nmol_m3,
            fLAI,
            nL,
        )

        expected_output = 0.0
        assert isclose(out, expected_output, abs_tol=1e-3)

    @pytest.mark.parametrize(['Gsto_l', 'expected_output'], [
        ([99.9], 0.0065704880),
    ])
    def test_calculates_fst_single_layer(self, Gsto_l, expected_output):
        nL = 1
        Rsto_l = [41000 / Gsto_l[iL] for iL in range(nL)]
        Rext = [300 for _ in range(nL)]
        O3_nmol_m3 = [3]
        fLAI = [1.0]
        Rb_l = [10]
        out = calc_fst_leaf(
            Gsto_l,
            Rb_l,
            Rsto_l,
            Rext,
            O3_nmol_m3,
            fLAI,
            nL,
        )

        assert isclose(out, expected_output, abs_tol=1e-3)

    @pytest.mark.parametrize(['Gsto_l', 'expected_output'], [
        ([99.9, 99.9, 99.9], 0.0065704880),
    ])
    def test_increasing_O3_increases_fst(self, Gsto_l, expected_output):
        nL = 3
        Rsto_l = [41000 / Gsto_l[iL] for iL in range(nL)]
        Rext = [300 for _ in range(nL)]
        O3_nmol_m3 = [3, 2, 1]
        fLAI = [0.1, 0.5, 0.4]
        Rb_l = [10, 20, 30]
        out_high = calc_fst_leaf(
            Gsto_l,
            Rb_l,
            Rsto_l,
            Rext,
            [o * 10 for o in O3_nmol_m3],
            fLAI,
            nL,
        )
        out_low = calc_fst_leaf(
            Gsto_l,
            Rb_l,
            Rsto_l,
            Rext,
            O3_nmol_m3,
            fLAI,
            nL,
        )
        assert out_high > out_low


def test_O3up_accumulation():
    """Test O3up_accumulation output."""
    O3up_acc = calc_O3up_accumulation(
        O3up=0.533,
        O3up_prev=0.577,
        O3up_acc=2.1,
        td_dd=23.22,
        td_dd_prev=11.42,
    )
    assert isclose(O3up_acc, 8.649, abs_tol=1e-3)


def test_O3_ppb_to_nmol():
    O3_nmol = O3_ppb_to_nmol(20, 91.0, 40)
    assert isclose(O3_nmol, 1493.3715602, abs_tol=1e-3)


def test_O3_ppb_to_nmol_match_do3se_ui_row_0():
    O3_nmol = O3_ppb_to_nmol(3.27999997139, 91.8600006104, 28.27368927)
    assert isclose(O3_nmol, 1130.00439453, abs_tol=1e-3)


def test_O3_ppb_to_nmol_match_do3se_ui_row_3679():
    O3_nmol = O3_ppb_to_nmol(11.5200004578, 92.9100036621, 15.3285636902)
    assert isclose(O3_nmol, 601.698059082, abs_tol=1e-3)


def test_stomatal_flux_rate():
    flux_rate = stomatal_flux_rate(
        leaf_r_Rsto=100,
        leaf_r_Rext=100,
        leaf_r_Rb=100,
    )
    assert isclose(flux_rate, 0.0033, abs_tol=1e-3)


def test_calc_POD_flag_leaf_has_not_emerged():
    out = calc_POD(1, 2, False, False, 3, 4)
    assert isclose(out.POD_0, 0.0, abs_tol=1e-3)
    assert isclose(out.POD_Y, 0.0, abs_tol=1e-3)


def test_calc_POD():
    out = calc_POD(1, 2, True, False, 3, 4)
    assert isclose(out.POD_0, 1.0144, abs_tol=1e-3)
    assert isclose(out.POD_Y, 2.0036, abs_tol=1e-3)


def test_calc_POD_match_DO3SE_UI_row_3655():
    out = calc_POD(0, 0, True, False, 6, 0.226157203317)
    assert isclose(out.POD_0, 0.000814165920019, abs_tol=1e-8)
    assert isclose(out.POD_Y, 0, abs_tol=1e-3)


def test_calc_FO3_eff():
    FO3_eff = calc_FO3_eff(
        fst=99.9,
        FO3_eff_prev=99.9,
        F_0=99.9,
    )
    assert isclose(FO3_eff, 99.9, abs_tol=1e-3)


def test_calc_OT():
    out = calc_OT_acc(
        is_daylight=True,
        f_phen=99.9,
        leaf_f_phen=99.9,
        micro_O3=99.9,
        AOT_0_prev=99.9,
        AOT_40_prev=99.9,
    )

    assert isclose(out.OT_0, 0.0999, abs_tol=1e-3)
    assert isclose(out.OT_40, 0.0599, abs_tol=1e-3)
    assert isclose(out.AOT_0, 99.999, abs_tol=1e-3)
    assert isclose(out.AOT_40, 99.959, abs_tol=1e-3)

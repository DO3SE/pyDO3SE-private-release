import pytest
from math import isclose
from pyDO3SE.plugins.resistance.model import Leaf_Resistance_Model, Resistance_Model
from .helpers import (
    calc_FO3_eff,
    calc_O3up_accumulation,
    calc_OT_acc,
    calc_POD,
    calc_fst,
    O3_ppb_to_nmol,
    calc_fst_leaf,
    calc_leaf_resistance_model,
    calc_resistance_model,
    stomatal_flux_rate,
)


def test_calc_fst():
    fst = calc_fst(
        Gsto_l=495.529266357,
        Rb_l=10,
        Rsto_l=82.7398147583,
        Rext=2500,
        O3_nmol_m3=11007.807805097667,
    )
    assert isclose(fst, 118.273517, abs_tol=1e-3)


def test_calc_fst_match_do3se_ui_row_3655():
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
        ([99.9,99.9,99.9], 0.0065704880),
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


def test_calc_resistance_model(snapshot):
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
        rext_calc_method="calculate_SAI",
        ra_calc_method="simple",
    )
    assert isinstance(rmodel, Resistance_Model)
    snapshot.assert_match(rmodel)
    # assert isclose(rmodel.Ra_c, 11.0274, abs_tol=1e-4)
    # assert isclose(rmodel.Ra, 8.5915, abs_tol=1e-4)
    # assert isclose(rmodel.Rb, 5.5203, abs_tol=1e-4)
    # assert isclose(rmodel.Rinc[0], 125.9999, abs_tol=1e-4)
    # assert isclose(rmodel.Rinc[1], 209.9999, abs_tol=1e-4)
    # assert isclose(rmodel.Rext[0], 833.3333, abs_tol=1e-4)
    # assert isclose(rmodel.Rext[1], 500.0, abs_tol=1e-4)
    # assert isclose(rmodel.Rsto[0], 13666.6666, abs_tol=1e-4)
    # assert isclose(rmodel.Rsto[1], 8200.0, abs_tol=1e-4)
    # assert isclose(rmodel.Rgs, 123.0, abs_tol=1e-4)
    # assert isclose(rmodel.Rsur[0], 790.9609, abs_tol=1e-4)
    # assert isclose(rmodel.Rsur[1], 476.7846, abs_tol=1e-4)
    # assert isclose(rmodel.Rtotal[0], 228.8716, abs_tol=1e-4)
    # assert isclose(rmodel.Rtotal[1], 196.0636, abs_tol=1e-4)


def test_calc_resistance_model_match_ui_row_0():
    """Check calc_resistance_model outputs match DO3SE UI."""
    rmodel = calc_resistance_model(
        nL=1,
        nLC=1,
        ustar=0.19406299293,
        canopy_height=1,
        SAI_values=[[3]],
        LAI_values=[[3]],
        mean_gsto_values=[[0]],  # CHCEK THIS
        Rsoil=200,
        ra_calc_method="simple",
        rsur_calc_method="single_layer",
    )
    assert isinstance(rmodel, Resistance_Model)
    assert isclose(rmodel.Ra, 64.121711731, abs_tol=1e-4)
    assert isclose(rmodel.Rb, 31.2906856537, abs_tol=1e-4)
    assert isclose(rmodel.Rinc[0], 216.424575806, abs_tol=1e-4)
    assert isclose(rmodel.Rsto[0], 100000, abs_tol=1e-4)
    assert isclose(rmodel.Rgs, 200, abs_tol=1e-4)
    assert isclose(rmodel.Rsur[0], 276.901275635, abs_tol=1e-4)
    assert isclose(rmodel.Rext[0], 2500, abs_tol=1e-4)


def test_calc_resistance_model_match_ui_row_20():
    """Check calc_resistance_model outputs match DO3SE UI."""
    rmodel = calc_resistance_model(
        nL=1,
        nLC=1,
        ustar=0.0373198054731,
        canopy_height=1,
        SAI_values=[[3]],
        LAI_values=[[3]],
        mean_gsto_values=[[0]],  # CHCEK THIS
        Rsoil=200,
        ra_calc_method="simple",
        rsur_calc_method="single_layer",
    )
    assert isinstance(rmodel, Resistance_Model)
    assert isclose(rmodel.Ra, 333.432891846, abs_tol=1e-4)
    assert isclose(rmodel.Rb, 162.71156311, abs_tol=1e-4)
    assert isclose(rmodel.Rinc[0], 1125.40783691, abs_tol=1e-4)
    assert isclose(rmodel.Rsto[0], 100000, abs_tol=1e-4)
    assert isclose(rmodel.Rgs, 200, abs_tol=1e-4)
    assert isclose(rmodel.Rsur[0], 509.039337158, abs_tol=1e-4)
    assert isclose(rmodel.Rext[0], 2500, abs_tol=1e-4)


def test_calc_leaf_resistance_model():
    rlmodel = calc_leaf_resistance_model(
        nL=3,
        Lm=0.01,
        u_per_layer=[10.1, 10.2, 10.3],
        leaf_gsto_per_layer=[100, 200, 300],
    )
    assert isinstance(rlmodel, Leaf_Resistance_Model)
    assert isclose(rlmodel.Rb[0], 6.1433, abs_tol=1e-4)
    assert isclose(rlmodel.Rext[0], 2500.0, abs_tol=1e-4)
    assert isclose(rlmodel.Rsto[0], 410.0, abs_tol=1e-4)

    assert isclose(rlmodel.Rb[2], 6.0833, abs_tol=1e-4)
    assert isclose(rlmodel.Rext[2], 2500.0, abs_tol=1e-4)
    assert isclose(rlmodel.Rsto[2], 136.6666, abs_tol=1e-4)


# Deprecated as function is for single component
def test_calc_leaf_resistance_model_single_layer():
    rlmodel = calc_leaf_resistance_model(
        nL=1,
        Lm=0.01,
        u_per_layer=[10.1],
        leaf_gsto_per_layer=[100]
    )
    assert isinstance(rlmodel, Leaf_Resistance_Model)
    assert isclose(rlmodel.Rb[0], 6.1433, abs_tol=1e-4)
    assert isclose(rlmodel.Rext[0], 2500.0, abs_tol=1e-4)
    assert isclose(rlmodel.Rsto[0], 410.0, abs_tol=1e-4)


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

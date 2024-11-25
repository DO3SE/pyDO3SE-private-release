from math import isclose
from pyDO3SE.plugins.gsto.multiplicative import apply_VPD_crit, calc_gsto_leaf, calc_gsto_mean, multiplicative


def test_multiplicative(snapshot):
    """Test multiplicative at day 150 hour 12."""
    out = multiplicative(
        VPD_crit=8.0,
        VPD_dd=0.72929,
        initial_leaf_gsto=100.0,
        initial_mean_gsto=80.0,

        fmin=1.0,
        gmax=450.0,
        gmorph=1.0,
        f_phen=1.0,
        f_light=1.0,
        f_temp=1.0,
        f_VPD=1.0,
        f_SW=1.0,
        leaf_f_phen=1.0,
        f_O3=1.0,
        leaf_f_light=1.0,
    )

    assert out.new_leaf_gsto == 450.0
    assert out.new_mean_gsto == 450.0


def test_apply_VPD_crit():
    out = apply_VPD_crit(
        VPD_crit=8.0,
        VPD_dd=0.72929,
        old_gsto=100.0,
        new_gsto=110.0,
    )

    assert isclose(out, 110.0, abs_tol=1e-3)


def test_apply_VPD_crit_above_limit():
    out = apply_VPD_crit(
        VPD_crit=8.0,
        VPD_dd=10.0,
        old_gsto=100.0,
        new_gsto=110.0,
    )

    assert isclose(out, 100.0, abs_tol=1e-3)


def test_compare_mean_and_leaf_gsto():
    # When all fraction values are 1.0 then both gsto_l and mean_gsto == gmax
    gmax = 2000
    leaf_f_phen = 1.0
    f_O3 = 1.0
    leaf_f_light = 1.0
    fmin = 0.2
    f_temp = 1.0
    f_VPD = 1.0
    f_SW = 1.0

    gmorph = 1.0
    f_phen = 1.0
    f_light = 1.0

    gsto_leaf = calc_gsto_leaf(
        gmax,
        leaf_f_phen,
        f_O3,
        leaf_f_light,
        fmin,
        f_temp,
        f_VPD,
        f_SW,
    )
    gsto_mean = calc_gsto_mean(
        gmax,
        gmorph,
        f_phen,
        f_light,
        fmin,
        f_temp,
        f_VPD,
        f_SW,
    )

    assert gsto_leaf == gmax
    assert gsto_mean == gmax

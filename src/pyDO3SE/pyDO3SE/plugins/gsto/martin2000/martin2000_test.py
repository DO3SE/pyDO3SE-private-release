from cmath import isclose
from pyDO3SE.plugins.gsto.martin2000.martin2000 import calc_O3_effect_on_V_xmax_25


def test_calc_O3_effect_on_V_xmax_25():
    out = calc_O3_effect_on_V_xmax_25(
        K_z=99.9,
        FO3_eff=9999.9,
        V_cmax_25_in=180,
        J_max_25_in=400,
    )
    assert isclose(out.V_cmax_25, 178.2018, abs_tol=1e-3)
    assert isclose(out.J_max_25, 396.0040, abs_tol=1e-3)

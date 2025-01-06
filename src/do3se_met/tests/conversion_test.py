from math import isclose
from do3se_met.conversion import convert_ppb_to_nmol


def test_convert_ppb_to_umol_O3():
    o3 = convert_ppb_to_nmol(
        value_ppb=40,
        Ts_C=20.0,
        P=91.0,
        gas_molecular_weight=48.0,
    )

    assert isclose(o3, 1493.3715, abs_tol=1e-3)

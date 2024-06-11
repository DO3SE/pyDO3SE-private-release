from .v_cmax import multilayer_vcmax25

def test_multilayer_vcmax():
    nL = 5
    out = multilayer_vcmax25(
       layer_LAI=[0.2 for _ in range(nL)],
    )
    assert len(out) == nL
    assert out[0] > out[-1]
    assert all([a > b for a, b in zip(out[:-1], out[1:])])

def test_multilayer_vcmax_zero_lai():
    nL = 5
    out = multilayer_vcmax25(
       layer_LAI=[0.0 for _ in range(nL)],
    )
    assert len(out) == nL


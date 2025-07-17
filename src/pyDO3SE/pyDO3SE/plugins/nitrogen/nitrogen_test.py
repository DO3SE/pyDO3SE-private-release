from pyDO3SE.plugins.nitrogen.nitrogen import calc_multilayer_nitrogen

class TestCalcNitrogenMultilayer:
    nL = 5
    out = calc_multilayer_nitrogen(
       layer_LAI=[0.2 for _ in range(nL)],
    )
    assert len(out) == nL
    assert out[0] > out[-1]
    assert all([a > b for a, b in zip(out[:-1], out[1:])])

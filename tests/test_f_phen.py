import pytest
from functools import partial
from itertools import accumulate
from math import isclose

from do3se_phenology.f_phen import (
    tt_f_phen_simple_PLF_value,
    tt_leaf_f_phen_PLF_value,
    f_phen_simple_PLF_value,
    f_phen_complex_PLF_value,
    leaf_f_phen_PLF_value,
    calc_leaf_f_phen_effect_on_V_cmax_25,
    tt_leaf_f_phen_PLF_range,
    f_phen_simple_PLF_range,
    get_fphen_PLF_fn,
    get_leaf_fphen_PLF_fn,
)


@pytest.mark.parametrize(
    "td,expected",
    [
        (0, 0.2),  # before SGS set to f_min
        (10, 0.2),  # at SGS
        (10 + 70, 0.2),
        (10 + 100, 0.282),
        (10 + 360, 1.0),
        (10 + 1145, 1.0),
        (10 + 1146, 0.998),
        (10 + 1845, 0.0),
        (10 + 1846, 0.0),
    ],
)
def test_tt_f_phen_simple_PLF(td, expected):
    """Tests the tt_f_phen_simple_PLF function returns the correct value."""
    f_phen = tt_f_phen_simple_PLF_value(
        td=td,
        t_f_phen_a=70,
        t_f_phen_b=360,
        t_f_phen_c=1145,
        t_f_phen_d=1845,
        f_phen_min=0.2,
        td_at_sgs=10,
    )

    assert isclose(f_phen, expected, abs_tol=1e-3)


def test_tt_f_phen_simple_PLF_snapshot():
    """Tests the tt_f_phen_simple_PLF function returns the correct value."""
    temperature_data = [2 for _ in range(1253)]
    accumulated_temperatures = list(accumulate(temperature_data))
    _get_f_phen = partial(
        tt_f_phen_simple_PLF_value,
        t_f_phen_a=70,
        t_f_phen_b=360,
        t_f_phen_c=1145,
        t_f_phen_d=1845,
        f_phen_min=0.2,
        td_at_sgs=100,
    )
    f_phen = list(map(_get_f_phen, accumulated_temperatures))

    assert isclose(f_phen[0], 0.2, abs_tol=1e-3)
    assert isclose(f_phen[10], 0.2, abs_tol=1e-3)
    assert isclose(f_phen[100], 0.2883, abs_tol=1e-3)
    assert isclose(f_phen[500], 1, abs_tol=1e-3)
    assert isclose(f_phen[1000], 0, abs_tol=1e-3)
    assert isclose(f_phen[-1], 0, abs_tol=1e-3)


@pytest.mark.parametrize(
    "td,expected",
    [
        (0, 0.0),  # before SGS set to f_min
        (100, 0.0),  # at SGS
        (100 + 950 - 1, 0.0),  # between SGS and Astart
        (100 + 950 + 1, 1.0),  # between Astart and mid_anthesis
        (100 + 950 + 200, 1.0),  # at mid_anthesis
        (100 + 950 + 200 + 250, 1.0),  # at start of senescence
        (100 + 950 + 200 + 525, 0.7),  # at mid of senescence
        (100 + 950 + 200 + 700, 0.0),  # at end of senescence
    ],
)
def test_tt_leaf_f_phen_PLF(td, expected):
    """Tests the tt_leaf_f_phen_PLF function returns the correct value."""
    leaf_f_phen = tt_leaf_f_phen_PLF_value(
        td=td,
        t_leaf_f_phen_a=0.3,
        t_leaf_f_phen_b=0.7,
        t_leaf_f_phen_e=200,
        t_leaf_f_phen_g=250,
        t_leaf_f_phen_h=525,
        t_leaf_f_phen_i=700,
        t_astart=950,
        td_at_sgs=100,
    )

    assert isclose(leaf_f_phen, expected, abs_tol=1e-3)


def test_tt_leaf_f_phen_PLF_snapshot():
    """Tests the tt_leaf_f_phen_PLF function returns the correct value."""
    temperature_data = [2 for _ in range(1253)]
    accumulated_temperatures = list(accumulate(temperature_data))
    _get_f_phen = partial(
        tt_leaf_f_phen_PLF_value,
        t_leaf_f_phen_a=0.3,
        t_leaf_f_phen_b=0.7,
        t_leaf_f_phen_e=200,
        t_leaf_f_phen_g=100,
        t_leaf_f_phen_h=525,
        t_leaf_f_phen_i=700,
        t_astart=950,
        td_at_sgs=100,
    )
    leaf_f_phen = list(map(_get_f_phen, accumulated_temperatures))

    assert isclose(leaf_f_phen[0], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[10], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[600], 1, abs_tol=1e-3)
    assert isclose(leaf_f_phen[800], 0.8221, abs_tol=1e-3)
    assert isclose(leaf_f_phen[900], 0.592, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1000], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[-1], 0, abs_tol=1e-3)


def test_tt_leaf_f_phen_PLF_range():
    """Tests the tt_leaf_f_phen_PLF function returns the correct value."""
    temperature_data = [2 for _ in range(1253)]
    accumulated_temperatures = list(accumulate(temperature_data))

    leaf_f_phen = tt_leaf_f_phen_PLF_range(
        td_list=accumulated_temperatures,
        t_leaf_f_phen_a=0.3,
        t_leaf_f_phen_b=0.7,
        t_leaf_f_phen_e=200,
        t_leaf_f_phen_g=100,
        t_leaf_f_phen_h=525,
        t_leaf_f_phen_i=700,
        t_astart=950,
        td_at_sgs=100,
    )

    assert isclose(leaf_f_phen[0], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[10], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[600], 1, abs_tol=1e-3)
    assert isclose(leaf_f_phen[800], 0.8221, abs_tol=1e-3)
    assert isclose(leaf_f_phen[900], 0.592, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1000], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[1100], 0, abs_tol=1e-3)
    assert isclose(leaf_f_phen[-1], 0, abs_tol=1e-3)


@pytest.mark.parametrize(
    "dd,expected",
    [
        (0, 0.0),  # at dd = 0
        (19, 0.0),  # Just before SGS
        (20, 0.2),  # at dd = SGS
        (20 + 1, 0.226),  # just after SGS
        (20 + 30, 1.0),  # at dd = SGS + f_phen_1
        (20 + 30 + 1, 1.0),  # just after SGS + f_phen_1
        (200 - 80, 1.0),  # at dd = SGS + f_phen_1 + f_phen_4
        (200 - 80 + 1, 0.9925),  # at dd = SGS + f_phen_1 + f_phen_4
        (200 + 1, 0.0),  # just after EGS
        (364, 0.0),  # Should wrap around to 0 - 365
        (400, 0.6),  # Should wrap around to 0 - 365
        (365, 0.0),  # Should wrap around to 0 - 365
        (365 + 20 + 30, 1.0),  # Should wrap around to 0 - 365
    ],
)
def test_f_phen_simple_PLF(dd, expected):
    """Tests the f_phen_simple_PLF function returns the correct value."""

    f_phen_1 = 30
    f_phen_4 = 80
    f_phen_a = 0.2
    f_phen_c = 1.0
    f_phen_e = 0.4
    SGS = 20
    EGS = 200
    out = f_phen_simple_PLF_value(
        f_phen_1=f_phen_1,
        f_phen_4=f_phen_4,
        f_phen_a=f_phen_a,
        f_phen_c=f_phen_c,
        f_phen_e=f_phen_e,
        SGS=SGS,
        EGS=EGS,
        dd=dd,
        f_phen_min=0.0,
    )
    assert isclose(out, expected, abs_tol=1e-3)


def test_f_phen_complex_PLF():
    """Tests the f_phen_complex_PLF function returns the correct value.

    The values to check have been matched against the DO3SE UI model
    """

    def calc_fphen_demo_p(dd):
        return f_phen_complex_PLF_value(
            f_phen_1=0,
            f_phen_2=1,
            f_phen_3=1,
            f_phen_4=45,
            f_phen_limA=0,
            f_phen_limB=0,
            f_phen_a=0.1,
            f_phen_b=1.0,
            f_phen_c=1.0,
            f_phen_d=1.0,
            f_phen_e=0.1,
            SGS=118,
            EGS=210,
            dd=dd,
        )

    assert isclose(calc_fphen_demo_p(dd=117), 0.0)  # ddadj = 364
    assert isclose(calc_fphen_demo_p(dd=118), 1.0)  # ddadj = 0 (365)
    assert isclose(calc_fphen_demo_p(dd=165), 1.0)
    assert isclose(calc_fphen_demo_p(dd=166), 0.98)
    assert isclose(calc_fphen_demo_p(dd=210), 0.1)
    assert isclose(calc_fphen_demo_p(dd=211), 0.0)


def test_f_phen_range():
    """Test that f_phen_simple_PLF_range produces outputs that match f_phen_simple"""
    dd_values = [i for i in range(500)]
    f_phen_1 = 30
    f_phen_4 = 80
    f_phen_a = 0.2
    f_phen_c = 1.0
    f_phen_e = 0.4
    SGS = 20
    EGS = 200
    f_phen_min = 0.0
    out = f_phen_simple_PLF_range(
        dd_values,
        f_phen_1,
        f_phen_4,
        f_phen_a,
        f_phen_c,
        f_phen_e,
        SGS,
        EGS,
        f_phen_min,
    )
    values_to_check = [
        (0, 0.0),  # at dd = 0
        (19, 0.0),  # Just before SGS
        (20, 0.2),  # at dd = SGS
        (20 + 1, 0.226),  # just after SGS
        (20 + 30, 1.0),  # at dd = SGS + f_phen_1
        (20 + 30 + 1, 1.0),  # just after SGS + f_phen_1
        (200 - 80, 1.0),  # at dd = SGS + f_phen_1 + f_phen_4
        (200 - 80 + 1, 0.9925),  # at dd = SGS + f_phen_1 + f_phen_4
        (200 + 1, 0.0),  # just after EGS
        (364, 0.0),  # Should wrap around to 0 - 365
        (400, 0.6),  # Should wrap around to 0 - 365
        (365, 0.0),  # Should wrap around to 0 - 365
        (365 + 20 + 30, 1.0),  # Should wrap around to 0 - 365
    ]
    for dd_target, expected in values_to_check:
        f_phen_out = next(f_phen for dd, f_phen in zip(dd_values, out) if dd == dd_target)
        assert isclose(f_phen_out,  expected, abs_tol=1e-3), f"f_phen_out({f_phen_out}) != expected({expected}) for dd({dd_target})"



@pytest.mark.parametrize(
    "dd,expected",
    [
        (0, 0.0),  # before SGS set to f_min
        (50 - 1, 0.0),  # before Astart
        (50, 0.5),  # at Astart
        (50 + 1, 0.533),  # after Astart
        (50 + 15, 1.0),  # Astart + leaf_f_phen_1
        (100 - 30, 1.0),  # AEnd - leaf_f_phen_2
        (100 - 1, 0.13),  # Before AEnd
        (100, 0.1),  # at AEnd
        (100 + 1, 0.0),  # after AEnd
        # Check wrap around
        (365 + 0, 0.0),  # before SGS set to f_min
        (365 + 50 - 1, 0.0),  # before Astart
        (365 + 50, 0.5),  # at Astart
        (365 + 50 + 1, 0.533),  # after Astart
        (365 + 50 + 15, 1.0),  # Astart + leaf_f_phen_1
        (365 + 100 - 30, 1.0),  # AEnd - leaf_f_phen_2
        (365 + 100 - 1, 0.13),  # Before AEnd
        (365 + 100, 0.1),  # at AEnd
        (365 + 100 + 1, 0.0),  # after AEnd
    ],
)
def test_leaf_f_phen_PLF(
    dd,
    expected,
):
    """Tests the leaf_f_phen_PLF function returns the correct value."""
    out = leaf_f_phen_PLF_value(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.5,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.1,
        Astart=50,
        Aend=100,
        dd=dd,
    )
    assert isclose(out, expected, abs_tol=1e-3)


@pytest.mark.parametrize(
    "dd,expected",
    [
        (340 - 10, 0.0),  # before Astart set to f_min
        (340 - 1, 0.0),  # before Astart
        (340, 0.5),  # at Astart
        (340 + 1, 0.533),  # after Astart
        (340 + 15, 1.0),  # Astart + leaf_f_phen_1
        (35 - 30, 1.0),  # AEnd - leaf_f_phen_2
        (35 - 1, 0.13),  # Before AEnd
        (35, 0.1),  # at AEnd
        (35 + 1, 0.0),  # after AEnd
        # Check wrap around
        (365 + 340-10, 0.0),  # before SGS set to f_min
        (365 + 340 - 1, 0.0),  # before Astart
        (365 + 340, 0.5),  # at Astart
        (365 + 340 + 1, 0.533),  # after Astart
        (365 + 340 + 15, 1.0),  # Astart + leaf_f_phen_1
        (365 + 35 - 30, 1.0),  # AEnd - leaf_f_phen_2
        (365 + 35 - 1, 0.13),  # Before AEnd
        (365 + 35, 0.1),  # at AEnd
        (365 + 35 + 1, 0.0),  # after AEnd
    ],
)
def test_leaf_f_phen_PLF_inversed_aend(
    dd,
    expected,
):
    """Tests the leaf_f_phen_PLF function returns the correct value."""
    out = leaf_f_phen_PLF_value(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.5,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.1,
        Astart=340,
        Aend=35,
        dd=dd,
    )
    assert isclose(out, expected, abs_tol=1e-3)


def test_leaf_f_phen_PLF_snapshot():
    """Tests the leaf_f_phen_PLF function returns the correct value."""
    out = leaf_f_phen_PLF_value(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=12.0,
        leaf_f_phen_b=30.0,
        leaf_f_phen_c=8.0,
        Astart=1,
        Aend=100,
        dd=74,
    )
    assert isclose(out, 27.066, abs_tol=1e-3)


def test_leaf_f_phen_PLF_Astart():
    """Tests the leaf_f_phen_PLF function returns the correct value.

    Make sure that the fphen starts to increase at dd==Astart.
    """
    out = leaf_f_phen_PLF_value(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=152,
    )
    assert isclose(out, 0.0, abs_tol=1e-3)
    out = leaf_f_phen_PLF_value(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=153,
    )
    assert isclose(out, 0.8, abs_tol=1e-3)


def test_leaf_f_phen_PLF_Aend():
    """Tests the leaf_f_phen_PLF function returns the correct value.

    Make sure that the fphen hits 0 at dd==Aend + 1.
    """
    out = leaf_f_phen_PLF_value(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=208,
    )
    assert isclose(out, 0.2, abs_tol=1e-3)

    out = leaf_f_phen_PLF_value(
        leaf_f_phen_1=15,
        leaf_f_phen_2=30,
        leaf_f_phen_a=0.8,
        leaf_f_phen_b=1.0,
        leaf_f_phen_c=0.2,
        Astart=153,
        Aend=208,
        dd=209,
    )
    assert isclose(out, 0.0, abs_tol=1e-3)


def test_calc_leaf_f_phen_effect_on_V_cmax_25():
    out = calc_leaf_f_phen_effect_on_V_cmax_25(
        V_cmax_25_in=180,
        J_max_25_in=400,
        leaf_f_phen=1.23,
    )
    assert isclose(out.V_cmax_25, 221.4, abs_tol=1e-3)
    assert isclose(out.J_max_25, 492.0, abs_tol=1e-3)



@pytest.mark.parametrize(
    "SGS, EGS, expected",
    [
        (1,365, [0,1,16,320,365, 366]),
        # Note that if the season length is the same we should always get the same PLF function
        (50,200, [0,1,16,106,151, 152]),
        (100,250, [0,1,16,106,151, 152]),
        (300,85, [0,1,16,106,151, 152]),
    ],
)
def test_get_fphen_PLF_fn(SGS: int, EGS: int, expected: list[int]):
    f_phen_1 = 15
    f_phen_4 = 45
    f_phen_a = 0.1
    f_phen_c = 1.0
    f_phen_e = 0.1
    gs_offset, fphen_values = get_fphen_PLF_fn(
        f_phen_1=f_phen_1,
        f_phen_4=f_phen_4,
        f_phen_a=f_phen_a,
        f_phen_c=f_phen_c,
        f_phen_e=f_phen_e,
        SGS=SGS,
        EGS=EGS,
        f_phen_min=0.0,
        offset_value=SGS-1,
    )

    assert gs_offset == expected
    assert fphen_values == [0.0, 0.1, 1.0, 1.0, 0.1, 0.0]


@pytest.mark.parametrize(
    "Astart, Aend, expected",
    [
        (50,200,  [0, 1, 16, 121, 151, 152]),
        # Note that if the season length is the same we should always get the same PLF function
        (100,160,  [0, 1, 16, 31.0, 61.0, 62.0]),
        (130,190,  [0, 1, 16, 31.0, 61.0, 62.0]),
        (340,35,  [0, 1, 16, 31.0, 61.0, 62.0]),
    ],
)
def test_get_leaf_fphen_PLF_fn(Astart: int,Aend: int, expected: list[int]):
    leaf_f_phen_1=15
    leaf_f_phen_2=30
    leaf_f_phen_a=0.5
    leaf_f_phen_b=1.0
    leaf_f_phen_c=0.1
    gs_offset, leaf_f_phen_values= get_leaf_fphen_PLF_fn(
        leaf_f_phen_1=leaf_f_phen_1,
        leaf_f_phen_2=leaf_f_phen_2,
        leaf_f_phen_a=leaf_f_phen_a,
        leaf_f_phen_b=leaf_f_phen_b,
        leaf_f_phen_c=leaf_f_phen_c,
        Astart=Astart,
        Aend=Aend,
        offset_value=Astart - 1,
    )

    assert gs_offset == expected
    assert leaf_f_phen_values == [0.0, 0.5, 1.0, 1.0, 0.1, 0.0]
from math import isclose
from do3se_phenology.latitude_function import (
    lat_function,
    lat_function_winter_wheat_china,
    lat_function_spring_wheat_europe,
)


def test_latitiude_function():
    out = lat_function(53, -2.6, 0, 376.9)
    assert isclose(out, 239, abs_tol=1)

    out = lat_function(118, -2.6, 0, 376.9)
    assert isclose(out, 71, abs_tol=1)


def test_latitiude_function_winter_wheat_china():
    out = lat_function_winter_wheat_china(53)
    assert isclose(out, 239, abs_tol=1)

    out = lat_function_winter_wheat_china(118)
    assert isclose(out, 71, abs_tol=1)


def test_latitiude_function_spring_wheat_europe():
    # These values from table 5.1 in Simpson et al., (2003)
    out = lat_function_spring_wheat_europe(50)
    assert isclose(out, 105, abs_tol=1)

    out = lat_function_spring_wheat_europe(50 + 1)
    assert isclose(out, 105 + 3, abs_tol=1)

    out = lat_function_spring_wheat_europe(50 - 1)
    assert isclose(out, 105 - 3, abs_tol=1)

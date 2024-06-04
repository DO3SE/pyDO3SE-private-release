from math import isclose
from .latitude_function import lat_function_winter_wheat
def test_latitiude_function():
    out = lat_function_winter_wheat(53)
    assert isclose(out, 239, abs_tol=1)

    out = lat_function_winter_wheat(118)
    assert isclose(out, 71, abs_tol=1)
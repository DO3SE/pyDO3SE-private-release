"""A module for calculating geospatial phenology."""
#
#  Latitude function for calculating start and end of growing season
#
from math import floor
from typing import Tuple


def estimate_latitude_SGS_EGS(
    lat: float,
    elev: float,
    a: float = 105.0,
    b: float = 297.0,
    c: float = 1.5,
    d: float = 2.0,
    e: float = 1000.0,
    f: float = 50,

) -> Tuple[int, int]:

    SGS = floor(a + ((lat - f) * c) + ((elev / e) * 10.0))
    EGS = floor(b - ((lat - f) * d) - ((elev / e) * 10.0))
    return SGS, EGS

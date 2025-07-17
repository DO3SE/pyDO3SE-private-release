"""Met helpers for calculating sun position."""

from math import cos, radians, sin
from typing import List, Tuple


def solar_noon(
    lon: float,
    dd: int,
) -> float:
    """Calculate solar noon (hour of day).

    Parameters
    ----------
    lon: float
        Longitude (degrees East)
    dd: int
        Day of year

    Returns
    -------
    solar_noon: float
        [description][unit]
    """
    # Solar noon correction for day of year
    f = radians(279.575 + (0.9856 * dd))
    e = (-104.7 * sin(f) + 596.2 * sin(2 * f) + 4.3 * sin(3 * f) - 12.7 * sin(4 * f)
         - 429.3 * cos(f) - 2.0 * cos(2 * f) + 19.3 * cos(3 * f)) / 3600  # noqa: W503
    # Calculate the longitudinal meridian
    lonm = int(lon / 15.0) * 15.0
    # Solar noon, with day of year and longitudinal correction
    LC = (lon - lonm) / 15
    solar_noon_value = 12 - LC - e
    return solar_noon_value


def solar_declination(dd: int) -> float:
    """Calculate solar declination (radians).

    Parameters
    ----------
    dd: int
        Day of year

    Returns
    -------
    solar_declination: float
        [description][unit]
    """
    return radians(-23.4 * cos((360 * radians((dd + 10) / 365.0))))


def calc_solar_elevation_list(
    lat: float,
    lon: float,
    day_range: Tuple[int, int]
) -> List[float]:
    """Calculate solar elevation for a range of days.

    Parameters
    ----------
    lat: float
        Latitude [degrees]
    lon: float
        Longitude [degrees]
    day_range: Tuple[int, int]
        Start and end day

    Returns
    -------
    sinB: List[float]
        SinB for each day in range
    """
    sinB_list = [calc_solar_elevation(lat, lon, dd, hr)
                 for dd in range(*day_range) for hr in range(24)]
    return sinB_list


def calc_solar_elevation(
    lat: float,
    lon: float,
    dd: int,
    hr: int,
) -> float:
    """Calculate solar elevation from lat, lon and time.

    Parameters
    ----------
    lat: float
        Latitude (degrees North)
    lon: float
        Longitude (degrees East)
    dd: int
        Day of year (1--365)
    hr: int
        Hour of day (0--23)

    Returns
    -------
    sinB: float
        sin() of solar elevation angle
    """
    t0_ = solar_noon(lon, dd)
    dec = solar_declination(dd)

    # Hour-angle of the sun
    h = radians(15 * (hr - t0_))

    # sin() of solar elevation angle
    sinB = sin(radians(lat)) * sin(dec) + cos(radians(lat)) * cos(dec) * cos(h)
    # TODO: should this line be removed? what effect does it have?  Does any
    #       use of sinB happen when sinB < 0?
    sinB = max(0.0, sinB)
    return sinB

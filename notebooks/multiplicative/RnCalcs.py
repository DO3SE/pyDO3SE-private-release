# %%
from math import cos, exp, pi, radians, sin, sqrt
from pyDO3SE.constants.physical_constants import GSC, SBC, T0
# %%
# PYTHON


def solar_noon(
    lon: float,
    dd: int,
) -> float:
    # Solar noon correction for day of year
    lonm = int(lon / 15.0) * 15.0

    f = radians(279.575 + (0.9856 * dd))  # OK
    e = (-104.7 * sin(f) + 596.2 * sin(2 * f) + 4.3 * sin(3 * f) - 12.7 * sin(4 * f)
         - 429.3 * cos(f) - 2.0 * cos(2 * f) + 19.3 * cos(3 * f)) / 3600  # noqa: W503
    # Calculate the longitudinal meridian
    # Solar noon, with day of year and longitudinal correction
    LC = (lon - lonm) / 15
    solar_noon_value = 12 - LC - e  # OK
    return solar_noon_value


solar_noon(56, 10)
# %%


def solar_declination(dd: int) -> float:
    return radians(-23.4 * cos((360 * radians((dd + 10) / 365.0))))

# %%


def calc_net_radiation(
    lat: float,
    lon: float,
    elev: float,
    albedo: float,
    dd: int,
    hr: int,
    sinB: float,
    R: float,
    Ts_C: float,
    eact: float,
) -> float:
    if sinB <= 0:
        net_radiation = 0.0
    else:
        # Latitude in radians
        lat_rad = radians(lat)

        # Convert global radiation W m-2 to MJ m-2 s-1
        R_MJ = R * 0.0036

        # Hour-angle of the sun
        t0_ = solar_noon(lon, dd)
        h = radians(15 * (hr - t0_))
        h1 = h - (pi / 24)
        h2 = h + (pi / 24)

        dr = 1 + (0.033 * cos(((2 * pi) / 365) * dd))
        dec = solar_declination(dd)
        # External radiation (with fix to stop div by zero)
        # TODO: fix this to be less hackish
        Re = max(0.00000000001,
                 ((12 * 60) / pi) * GSC * dr
                  * ((h2 - h1) * sin(lat_rad) * sin(dec)
                  + cos(lat_rad) * cos(dec) * (sin(h2) - sin(h1))))  # noqa:E501
        # TODO: what was this for?
        # Re = max(0.0, ((12*60)/pi)*Gsc*dr*sinB)

        # Calculate net longwave radiation
        assert elev is not None
        assert Re is not None

        pR = (0.75 + (2e-5 * elev)) * Re  # OK

        Rnl = max(0.0, (SBC * ((Ts_C + T0)**4)) * (0.34 - (0.14 * sqrt(eact)))
                  * ((1.35 * (min(1.0, R_MJ / pR))) - 0.35))  # noqa:W503 # OK
        Rns = (1 - albedo) * R_MJ  # OK

        net_radiation = max(0.0, Rns - Rnl)
    return net_radiation


# TODO: Why is this 0.0???
calc_net_radiation(
    lat=52.2,
    lon=-1.12,
    elev=1.1,
    albedo=0.2,
    dd=135,
    hr=1,
    sinB=0.3,
    R=900.0,
    Ts_C=20.1,
    eact=0.8,
)

# %%
# FORTRAN


def solar_declination(dd):
    return radians(-23.4 * cos(radians(360 * ((dd + 10) / 365.0))))


def calc_solar_noon_f(lon, dd):
    # Calculate the longitudinal meridian
    lonm = int(lon / 15.0) * 15.0

    # Solar noon correction for day of year
    f = radians(279.575 + (0.9856 * dd))
    e = (-104.7 * sin(f) + 596.2 * sin(2 * f) + 4.3 * sin(3 * f) - 12.7 * sin(4 * f)
         - 429.3 * cos(f) - 2.0 * cos(2 * f) + 19.3 * cos(3 * f)) / 3600

    # Solar noon, with day of year and longitudinal correction
    LC = (lon - lonm) / 15
    t0 = 12 - LC - e
    return t0


calc_solar_noon_f(56, 10)
# %%


def Calc_Rn(
    lat: float,
    lon: float,
    elev: float,
    albedo: float,
    dd: int,
    hr: int,
    sinB: float,
    R: float,
    Ts_C: float,
    eact: float,
):
    Gsc = 0.082            # Solar constant (MJ/m^2/min)
    SBC = 4.903e-9 / 24    # Stephan Boltzman constant

    lat_rad = radians(lat)

    if (sinB > 0):
        # Unit conversions
        R_MJ = R * 0.0036
        Ts_K = Ts_C + 273.15

        t0_ = solar_noon(lon, dd)
        h = radians(15 * (hr - t0_))  # OK

        # TODO: CHECK THIS IS CORRECT
        dec = solar_declination(dd)

        # Hour angle stuff
        h1 = h - (pi / 24)
        h2 = h + (pi / 24)

        dr = 1 + (0.033 * cos(((2 * pi) / 365) * dd))
        # External radiation (with fix to stop div by zero)
        # TODO: fix this to be less hackish
        Re = max(0.00000000001,
                 ((12 * 60) / pi) * Gsc * dr
                 * ((h2 - h1) * sin(lat_rad) * sin(dec)
                    + cos(lat_rad) * cos(dec) * (sin(h2) - sin(h1))))
        # Re = max(0.0, ((12*60)/pi)*Gsc*dr*sinB)

        # Calculate net longwave radiation
        pR = (0.75 + (2e-5 * elev)) * Re

        Rnl = max(0.0, (SBC * (Ts_K**4)) * (0.34 - (0.14 * sqrt(eact)))
                  * ((1.35 * (min(1.0, R_MJ / pR))) - 0.35))
        Rns = (1 - albedo) * R_MJ

        Rn = max(0.0, Rns - Rnl)
    else:
        Rn = 0

    # Calculate Rn in W/m2
    Rn_W = Rn * 277.8
    return Rn


Calc_Rn(
    lat=52.2,
    lon=-1.12,
    elev=20.0,
    albedo=0.2,
    dd=1,
    hr=12,
    sinB=0.8,
    R=900.0,
    Ts_C=20.1,
    eact=0.8,
)

def lat_function(
    lat: float,
    k: float,
    b: float,
    c: float
):
    """Latitude based phenology sowing date model.

    # TODO: May not work when lat is negative

    k * (lat - b) + c
    """
    return k * (lat - b) + c


def lat_function_winter_wheat_china(
    lat: float,
    k: float = -2.6,
    b: float = 376.9,
):
    """Latitude based phenology sowing date model for Chinese Winter Wheat.

    Based on figure 2a in Xiao et all(2015)

    Uses :function:`lat_function_wh
    lat_function(lat, k, b=0, c=b)

    References
    ----------
    Xiao et all 2015 - Spatiotemporal variability of winter wheat phenology
                       in response to weather and climate variability in China

    """
    return lat_function(lat, k, b=0, c=b)


def lat_function_spring_wheat_europe(
    lat: float,
    k: float = 3,
    b: float = 50,
    c: float = 105
):
    """Latitude based phenology sowing date model for Chinese Winter Wheat.

    Based on table 5.1 in Simpson et al., (2003)

    Uses :function:`lat_function`

    References
    ----------
    Simpson et al., (2003). Transboundary Acidification, Eutrophication and Ground Level Ozone in Europe PART I. Unified EMEP Model Description. EMEP Report 1/2003

    """
    return lat_function(lat, k, b, c)



def lat_function_forest_europe(
    lat: float,
    k: float = 1.5,
    b: float = 50,
    c: float = 105
):
    """Latitude based phenology sowing date model for Chinese Winter Wheat.

    Based on table 5.1 in Simpson et al., (2003)

    Uses :function:`lat_function`

    References
    ----------
    Simpson et al., (2003). Transboundary Acidification, Eutrophication and Ground Level Ozone in Europe PART I. Unified EMEP Model Description. EMEP Report 1/2003

    """
    return lat_function(lat, k, b, c)

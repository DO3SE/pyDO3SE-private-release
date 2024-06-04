# %%
def lat_function_winter_wheat(
    lat: float,
    k: float = -2.6,
    b: float = 376.9,
):
    """Latitude based phenology sowing date model for Chinese Winter Wheat.

    Based on figure 2a in Xiao et all(2015)

    References
    ----------
    Xiao et all 2015 - Spatiotemporal variability of winter wheat phenology
                       in response to weather and climate variability in China

    """
    return k * lat + b

# No longer used
# # > Canopy displacement (fraction of canopy height)
# CANOPY_D = 0.7
# # > Canopy roughness length (fraction of canopy height)
# CANOPY_Z0 = 0.1

izR = 50    # Intermediate "decoupled" height for transfer of O3 and windspeed

MIN_WINDSPEED = 0.01

MIN_DAYLIGHT_R: float = 50 # R values lower than this are assumed to by night.


# Gas constants [mol m-2 s-1]
LEAF_G_HEAT = 0.135
LEAF_G_H2O = 0.147
LEAF_G_CO2 = 0.110
LEAF_G_O3 = 0.105

######################################
# Fortran Conversion
# DO3SE_ModelConstants_ml.F90
######################################

# UNDEF = -999.0
# IUNDEF = -999

# TODO: see how much we can remove these
MAX_LC = 3  # < Maximum number of land covers (used in some static allocations)
MAX_LAYERS = 5  # < Maximum number of layers (used in some static allocations)
DEFAULT_LC = 1
DEFAULT_LAYERS = 1


DT = 60 * 60  # < Number of seconds in a timestep

# > Canopy displacement (fraction of canopy height)
CANOPY_D = 0.7
# > Canopy roughness length (fraction of canopy height)
CANOPY_Z0 = 0.1


izR = 50    # Intermediate "decoupled" height for transfer of O3 and windspeed

MIN_CANOPY_HEIGHT = 0.01

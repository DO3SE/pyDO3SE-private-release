from dataclasses import dataclass
from .ConfigEnums import LandCoverType


@dataclass(frozen=False)
class Config_Location:
    """Location properties.


    """

    lat: float = None     #: Latitude (degrees North)
    lon: float = None     #: Longitude (degrees East)
    elev: float = None    #: Elevation (m above sea level)
    albedo: float = None  #: Surface albedo [fraction]

    #: Soil resistances are set using values in Table S19 of Simpson, David, et al. 2012
    #: The EMEP MSC-W chemical transport model–technical description.
    Rsoil: float = 100     #: Base Soil resistance(O3) [s m-1]

    OTC: bool = False       #: Is open top chamber experiment

    #: Heights for transfering wind and ozone
    izr: float = 50         #: Decoupled height [m]
    z_u: float = None       #: Measurement height for windspeed[m]
    z_O3: float = None      #: Measurement height for O3[m]
    h_u: float = None       #: Canopy height for windspeed measurement, default = target canopy[m]
    h_O3: float = None      #: Canopy height for O3 measurement, default = target canopy[m]

    #: Calculated heights
    O3_d: float = None  # :  Ozone measured canopy displacement height [m]
    O3_z0: float = None  # :  Ozone measured canopy roughness length [m]

    u_d: float = None  # :  Wind measured canopy displacement height [m]
    u_z0: float = None  # :  Wind measured canopy roughness length [m]

    #: Override the start and end dates of the model. Otherwise the first and last dd values
    #: from the external data will be used. These days are inclusive.
    start_day: float = None
    end_day: float = None

    #: Year for DOY offset
    #: DOY (where year == zero_year + 1) == DOY + 365
    zero_year: int = None

    #: Hour offset (for timezone correction)
    hr_offset: float = 0
    #: Crop the data to whole days only
    crop_to_day_start: bool = False #: Crop to day start
    crop_to_day_end: bool = False #: Crop to day end

    land_cover_type: LandCoverType = LandCoverType.CROP

    #: If true then we reset all model values at a set day
    multi_season: bool = False
    multi_season_day: int = None

    #: CANOPY_D
    #: 0.78 for forest
    #: 0.70 for crops
    canopy_d: float = 0.7

    #: CANOPY_Z0
    #: 0.07 for forest
    #: 0.1 for crops
    canopy_zo: float = 0.1

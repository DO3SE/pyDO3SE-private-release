"""Contains the config data.

Config data is inputed at the start of the model and is constant throughout the model

It consists of dataclasses that can be accessed via dot notation
i.e.
Land_Cover.parameters[0].gsto.VPD_max

Developer Comments
------------------
Any changes to this file need to have a migration included in config_migration.py


Config Shape Module
-------------------

"""
from dataclasses import dataclass, field
from pyDO3SE.plugins.carbon_allocation.config import CarbonAllocationConfig
from pyDO3SE.External_State.External_State_Config import Config_Met
from pyDO3SE.plugins.soil_moisture.config import Soil_Moisture_Config
from pyDO3SE.Output.OutputConfig import OutputConfig

from .ConfigEnums import *
from .ConfigLocation import Config_Location
from .ConfigLandCover import Config_Land_Cover, Config_Land_Cover_Parameters
from .ConfigResistance import ResistanceConfig


@dataclass(frozen=False)
class Config_Shape:
    """The shape of the config data.

    Input config files should be a json file that matches the structure of this class and
    all subclasses.

    In the example below "Location" should match :class:`pyDO3SE.Config.ConfigLocation.Config_Location`.

    Example
    -------
    .. code-block:: json

        {
            "Location": {
                "lat": 53.2,
                "lon": -2.3
            },
            "Land_Cover": {
                "nL": 1,
                "nLC": 2,
                "parameters": [
                    {
                        "name": "demo land cover",
                        "season": ...,
                        "gsto":...
                        ...
                    }
                ]
                ...
            }
            ...
        }

    """
    VERSION: int = None  # This must be set by user to ensure using correct version

    id: str = "MISSING_ID"

    notes: str = ""

    # Location Config
    Location: Config_Location = field(default_factory=lambda: Config_Location())
    # Land Cover Config
    Land_Cover: Config_Land_Cover = field(default_factory=lambda: Config_Land_Cover())

    # Meteorological Config
    Met: Config_Met = field(default_factory=lambda: Config_Met())

    # SMD module configuration
    soil_moisture: Soil_Moisture_Config = field(default_factory=lambda: Soil_Moisture_Config())

    # Resistance Config
    resistance: ResistanceConfig = field(default_factory=lambda: ResistanceConfig())

    # Carbon Allocation Config
    carbon_allocation: CarbonAllocationConfig = field(
        default_factory=lambda: CarbonAllocationConfig())

    output: OutputConfig = field(default_factory=lambda: OutputConfig())

import json
from do3se_met.soil_moisture.config import SOIL_SANDY_LOAM, Soil_Moisture_Config
from do3se_met.soil_moisture.enums import SoilMoistureSource
from pyDO3SE.Config.Config_Shape import Config_Land_Cover, Config_Land_Cover_Parameters, Config_Shape, Config_Location


DEMO_CONFIG_INPUT = {
    "Location": {
        "lat": 52.2,
        "lon": -1.12,
    },
    "Land_Cover": {
        "parameters": [
            {
                "name": "p01",
                "gsto": {
                    "method": "photosynthesis"
                }
            }
        ]
    }
}


# DEMO_CONFIG_INPUT_ALT = {
#     "Location": {
#         "lat": 52.2,
#         "long": -1.12,
#     },
# }


DEMO_CONFIG_JSON = json.dumps(DEMO_CONFIG_INPUT)
# DEMO_CONFIG_JSON_ALT = json.dumps(DEMO_CONFIG_INPUT_ALT)

number_of_layers = 1
number_of_components = 1

DEMO_CONFIG = Config_Shape(
    Location=Config_Location(
        lat=52.2,
        lon=-1.12,
    ),
    Land_Cover=Config_Land_Cover(
        nL=number_of_layers,
        nLC=number_of_components,
        parameters=[
            Config_Land_Cover_Parameters(

            ),
        ]
    ),
    soil_moisture=Soil_Moisture_Config(
        source=SoilMoistureSource.P_M,
        soil_config=SOIL_SANDY_LOAM,
        initial_SWC=SOIL_SANDY_LOAM.FC,
    )
)

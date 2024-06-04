from data_helpers.cls_parsing import unpack
import json
import pytest

from pyDO3SE.Config import Config_Shape
from pyDO3SE.Config.config_loader import config_loader, process_json_config

# from pyDO3SE.Config.Config_Shape import Config_Location, Config_Land_Cover, \
#     Config_Land_Cover_Parameters, Config_Gsto


def test_config_loader(snapshot):
    config = config_loader('examples/spanish_wheat/configs/spanish_wheat_config.json', config_type='json')
    assert isinstance(config, Config_Shape)
    assert config.Location.lat == 40.43
    assert config.Location.lon == -3.7
    assert config.Land_Cover.parameters[0].gsto.method == "photosynthesis"
    snapshot.assert_match(unpack(config), 'Config')


def test_config_loader_invalid_type():
    with pytest.raises(Exception):
        config_loader('examples/spanish_wheat/configs/spanish_wheat_config.json', 'pdf')


def test_json_to_config(snapshot):
    DEMO_CONFIG_INPUT = {
        "Location": {
            "lat": 40.43,
            "lon": -3.7,
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
    DEMO_CONFIG_JSON = json.dumps(DEMO_CONFIG_INPUT)
    config: Config_Shape = process_json_config(DEMO_CONFIG_INPUT)
    assert isinstance(config, Config_Shape)
    snapshot.assert_match(unpack(config), 'Config from json')

    # assert config.Location == Config_Location(
    #     lat=52.2,
    #     lon=-1.12
    # )
    # assert config.Land_Cover.parameters[0] == Config_Land_Cover_Parameters(
    #     name="p01",
    #     gsto=Config_Gsto(
    #         method="photosynthesis"
    #     )
    # )
    # assert config.Land_Cover == Config_Land_Cover(
    #     parameters=[
    #         Config_Land_Cover_Parameters(
    #             name="p01",
    #             gsto=Config_Gsto(
    #                 method="photosynthesis"
    #             )
    #         )
    #     ]
    # )
    # assert config == Config_Shape(
    #     Location=Config_Location(
    #         lat=52.2,
    #         lon=-1.12
    #     ),
    #     Land_Cover=Config_Land_Cover(
    #         parameters=[
    #             Config_Land_Cover_Parameters(
    #                 name="p01",
    #                 gsto=Config_Gsto(
    #                     method="photosynthesis"
    #                 )
    #             )
    #         ]
    #     )
    # )

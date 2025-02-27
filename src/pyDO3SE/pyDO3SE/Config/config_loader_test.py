import json
import pytest
from data_helpers.dictionary_helpers import ListMergeMethods, merge_dictionaries
from pyDO3SE.util.test_utils import process_snapshot
from pyDO3SE.Config import Config_Shape
from pyDO3SE.Config.config_loader import config_loader, process_json_config
from pyDO3SE.Config.ConfigEnums import GstoMethods

# from pyDO3SE.Config.Config_Shape import Config_Location, Config_Land_Cover, \
#     Config_Land_Cover_Parameters, Config_Gsto


def test_config_loader(snapshot):
    config = config_loader(
        'examples/spanish_wheat/configs/spanish_wheat_config.json', config_type='json')
    assert isinstance(config, Config_Shape)
    assert config.Location.lat == 40.43
    assert config.Location.lon == -3.7
    assert config.Land_Cover.parameters[0].gsto.method == GstoMethods.PHOTOSYNTHESIS.value
    snapshot.assert_match(process_snapshot(config), 'Config')


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
                    },
                    "phenology": {
                        "dvi_interval": [[0, -1.0], [73.55, 0.0], [632.53, 1.0], [1471.0, 2.0]]
                    }
                }
            ]
        }
    }
    DEMO_CONFIG_JSON = json.dumps(DEMO_CONFIG_INPUT)
    config: Config_Shape = process_json_config(DEMO_CONFIG_INPUT)
    assert isinstance(config, Config_Shape)
    snapshot.assert_match(process_snapshot(config), 'Config from json')

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


def test_json_to_config_merge(snapshot):
    DEMO_CONFIG_INPUT_BASE = {
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
                    },
                    "phenology": {
                        "PRESET": None,
                        "dvi_interval": [[0, -1.0], [73.55, 0.0], [632.53, 1.0], [1471.0, 2.0]]
                    }
                },

            ]
        }
    }
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
                    },
                    "phenology": {
                        "PRESET": None,
                        "dvi_interval": [[0, -1.0], [73.55, 0.0], [632.53, 1.0], [1471.0, 2.0]]
                    }
                },

            ]
        }
    }
    merged_config_data = merge_dictionaries(
        DEMO_CONFIG_INPUT_BASE, DEMO_CONFIG_INPUT, ListMergeMethods.ZIP)
    assert merged_config_data['Land_Cover']['parameters'][0]['phenology']['PRESET'] == None
    # config: Config_Shape = process_json_config(DEMO_CONFIG_INPUT)
    # assert isinstance(config, Config_Shape)
    # snapshot.assert_match(unpack(config), 'Config from json')

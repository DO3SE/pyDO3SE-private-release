"""Use this file to import all needed env vars and other settings"""

import os
from typing import NamedTuple


class Settings(NamedTuple):
    DEBUG_MODE: bool
    MAX_NUM_OF_CANOPY_LAYERS: int
    MAX_NUM_OF_CUSTOM_LAYERS: int
    MAX_NUM_OF_CANOPY_COMPONENTS: int
    MAX_NUM_OF_LEAF_POPULATIONS: int


global_settings = Settings(
    os.environ.get("DEBUG_MODE") and True or False,
    int(os.environ.get("MAX_NUM_OF_CANOPY_LAYERS", 4)),
    int(os.environ.get("MAX_NUM_OF_CUSTOM_LAYERS", 4)),
    int(os.environ.get("MAX_NUM_OF_CANOPY_COMPONENTS", 1)),
    int(os.environ.get("MAX_NUM_OF_LEAF_POPULATIONS", 3)),
)


def settings():
    global global_settings
    return global_settings


def set_settings(**kwargs):
    global global_settings
    global_settings = global_settings._replace(
        **kwargs,
    )
    return global_settings

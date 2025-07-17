"""Generate config from options."""
from data_helpers.cls_parsing import rsetattr
from .Config_Shape import Config_Shape


required_config = [
    "Location.lat",
    "Location.elev",
]


def generate_config(**kwargs):
    base_config = Config_Shape()
    for k, v in kwargs.items():
        base_config = rsetattr(base_config, k, v)
    for k in required_config:
        if k not in kwargs.keys():
            base_config = rsetattr(base_config, k, "SET_ME!")
    return base_config

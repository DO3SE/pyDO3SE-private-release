from typing import List, Tuple, Callable, Union
from data_helpers.dictionary_helpers import ListMergeMethods, merge_dictionaries
from pathlib import Path
import json
from math import isclose
import pickle
import numpy as np

from data_helpers.cls_parsing import check_types, dict_to_cls

from pyDO3SE.Config import Config_Shape
from pyDO3SE.util.error_handling import ConfigError


def process_json_config(json_config_data: dict) -> Config_Shape:
    config: Union[Config_Shape, None] = dict_to_cls(json_config_data, Config_Shape)
    if not config:
        raise ConfigError("Failed to convert config to Config_Shape")
    # TODO: Should move these to config validation
    assert isclose(np.sum(config.Land_Cover.fLAI), 1, abs_tol=1e-5), (
        "fLAI must sum to 1 but got " + str(np.sum(config.Land_Cover.fLAI))
    )
    assert isclose(sum(config.Land_Cover.layer_height_frac), 1, abs_tol=1e-5), (
        "layer_height_frac must sum to 1 but got " + str(sum(config.Land_Cover.layer_height_frac))
    )
    assert config.Land_Cover.parameters is not None
    check_types(config)
    return config


def config_loader(
    config_path: Path | str,
    base_config_file_path: Path | str | None = None,
    config_type: str = "json",
    logger: Callable[[str], None] = print,
) -> Config_Shape:
    """loads a datafile into a config object"""
    logger("Loading config:", str(config_path))
    logger("Loading base config:", str(base_config_file_path))
    with open(config_path) as config_file:
        read_data = config_file.read()
        if config_type == "json":
            try:
                config_data = json.loads(read_data)
            except:
                raise ConfigError(f"Failed to load config from {config_path}")
        else:
            raise Exception("Invalid config file type")
    if base_config_file_path:
        with open(base_config_file_path) as base_config_file:
            read_data = base_config_file.read()
            if config_type == "json":
                try:
                    base_config_data = json.loads(read_data)
                except Exception:
                    raise ConfigError(f"Failed to load base config from{base_config_file}")
            else:
                raise Exception("Invalid config file type")
    merged_config_data = (
        base_config_file_path
        and merge_dictionaries(base_config_data, config_data, ListMergeMethods.ZIP)
        or config_data
    )
    config_object = process_json_config(merged_config_data)
    return config_object


def config_loader_pickled(
    config_location: Path | str,
) -> Config_Shape:
    """Load binary config file.

    NOTE: Laoding data with Pickle is not secure!


    Parameters
    ----------
    config_location : Path
        Location of config binary file

    Returns
    -------
    Config_Shape
        Config Shape Object

    """
    with open(config_location, "rb") as f:
        config_loaded = pickle.load(f)
    return config_loaded


def grid_config_loader(
    processed_config_dir: Path,
    grid_coords: List[Tuple[int, int]],
) -> List[Config_Shape]:
    """Get all the configs that will be ran for each grid cell.

    Parameters
    ----------
    processed_config_dir : Path
        path to configs that have already been created for each grid cell
    grid_coords: List[Tuple[int, int]]
        coordinates to run for this config

    Returns
    -------
    List[Config_Shape]
        A list of configs matching each coordinate

    """
    try:
        configs_loaded = [
            config_loader_pickled(
                f"{processed_config_dir}/{x}_{y}.config",
            )
            for x, y in grid_coords
        ]
    except FileNotFoundError:
        raise ConfigError(
            f"Could not find config for grid cells: {grid_coords}. Please run grid init first."
        )
    return configs_loaded

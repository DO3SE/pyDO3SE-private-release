from dataclasses import asdict
from copy import deepcopy
import json
from pathlib import Path
import pickle
from typing import IO
import warnings
from data_helpers.encoders import AdvancedJsonEncoder
from data_helpers.dictionary_helpers import ListMergeMethods, merge_dictionaries
from data_helpers.cls_parsing import dict_to_cls
from .Model_State import Model_State_Shape


def process_json_state(json_state_data: dict) -> Model_State_Shape:
    prev_hour = json_state_data.pop('prev_hour', None) # TODO: Manage previous hour data
    state = dict_to_cls(json_state_data, Model_State_Shape)

    # check_types(state)
    return state


def model_state_loader(
    state_location: Path,
    base_state_file: Path = None,
    file_type: str = 'json',
    strict: bool = True,
) -> Model_State_Shape:
    '''loads a datafile into a config object'''
    try:
        with open(state_location) as state_file:
            read_data = state_file.read()
            if file_type == 'json':
                state_data = json.loads(read_data)
            else:
                raise Exception('Invalid state file type')
        if base_state_file:
            with open(base_state_file) as state_file:
                read_data = state_file.read()
                if file_type == 'json':
                    base_state_data = json.loads(read_data)
                else:
                    raise Exception('Invalid state file type')
        merged_state_data = base_state_file and merge_dictionaries(
            base_state_data, state_data, ListMergeMethods.ZIP) or state_data
        state_object = process_json_state(merged_state_data)
        return state_object
    except IOError as e:
        if "No such file or directory" in str(e):
            if strict:
                raise e
            else:
                warnings.warn(f"Failed to load state file from: {state_location}\n {e}")
                return None
        else:
            raise e



def model_state_loader_quick(
    state_location: Path,
) -> Model_State_Shape:
    """Dump state to binary file with Pickle. DO NOT USE THIS ON API.

    NOTE: Loading data with Pickle is not secure!

    Parameters
    ----------
    state_location : Path
        Location of binary state file

    Returns
    -------
    Model_State_Shape
        Model State Object

    """
    with open(state_location, 'rb') as f:
        state_loaded = pickle.load(f)
    return state_loaded


def dump_state_to_string(
    state: Model_State_Shape,
) -> str:
    return json.dumps(
        asdict(state),
        cls=AdvancedJsonEncoder, indent=4, sort_keys=True)


def dump_state_to_file(
    state: Model_State_Shape,
    target_path: Path,
) -> IO:
    with open(target_path, 'w') as statefile:
        statefile.write(dump_state_to_string(state))


def dump_state_to_file_quick(
    state: Model_State_Shape,
    target_path: Path,
):
    """Pickle the model state.

    Parameters
    ----------
    state : Model_State_Shape
        [description]
    target_path : Path
        [description]
    """
    with open(target_path, 'wb') as statefile:
        pickle.dump(state, statefile)


def load_current_cell_state(
    prev_hour_state_path: Path,
    x: int,
    y: int,
) -> Model_State_Shape:
    previous_hour_state_path_tile = f"{prev_hour_state_path}/{x}_{y}.state"

    # TODO: Can we optimize this. # OPTIMIZE
    initial_state = model_state_loader_quick(previous_hour_state_path_tile)
    initial_state.prev_hour = None
    initial_state.prev_hour = deepcopy(initial_state)
    initial_state.temporal.row_index = 0
    return initial_state
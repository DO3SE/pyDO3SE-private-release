"""Define Model overrides.
"""
from pathlib import Path
from typing import List, NamedTuple, Callable
from proflow.Objects.Process import Process
from pyDO3SE.Config import Config_Shape
from pyDO3SE.External_State import External_State_Shape
from pyDO3SE.Model_State import Model_State_Shape

GET_INIT_PROCESSES_T = Callable[[Config_Shape], List[Process]]


class Main_Overrides(NamedTuple):
    """Main model overrides and options.

    Parameters
    ----------
    start_day: int
        Override the external data and config start day. default = None
    end_day: int
        Override the external data and config end day. default = None
    init_external_state_processes: List[Process]
        [description] default = None
    init_param_processes: List[Process]
        [description] default = None
    init_state_processes: List[Process]
        [description] default = None
    init_config_processes: List[Process]
        [description] default = None
    model_processes: GET_INIT_PROCESSES_T
        [description] default = None
    config_override: Config_Shape
        [description] default = None
    external_state_override: External_State_Shape
        [description] default = None
    input_state_override: Model_State_Shape
        [description] default = None
    allow_data_padding: bool
        If true then will pad start and end dates in data to match config start and end dates.
        default = False,
    output_to_netcdf: bool, default False
        Output final row data to netcdf if true, csv if false
    state_out_path: str, optional
        Location to save final state
    output_fields: List[str]
        Override the requried output fields
    e_state_overrides_file_path: Path
        Path to netcdf file with gridded external state overrides
    e_state_overrides_field_map: dict
        Dictionary mapping of field in netcdf to field in Config_Shape
    netcdf_chunks: dict, optional
        Chunks to use when loading netcdf data
    skip_state_init: bool, default False
        If true will not run state init processes
    dump_state_n: int, optional
        If set to a number then model will dump the state to file when row_index % n == 0
    run_validation: boolean,
        If true runs validation step on config
    debug: bool
        Set debug mode on

    """
    start_day: int = None
    end_day: int = None
    init_external_state_processes: List[Process] = None
    init_param_processes: List[Process] = None
    init_state_processes: List[Process] = None
    init_config_processes: List[Process] = None
    model_processes: GET_INIT_PROCESSES_T = None
    config_override: Config_Shape = None
    external_state_override: External_State_Shape = None
    input_state_override: Model_State_Shape = None
    allow_data_padding: bool = False
    output_to_netcdf: bool = False
    state_out_path: str = None
    output_fields: List[str] = None
    e_state_overrides_file_path: Path = None
    e_state_overrides_field_map: dict = None
    skip_state_init: bool = False
    dump_state_n: int = None
    run_validation: bool = True
    debug: bool = False

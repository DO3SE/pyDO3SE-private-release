

import numpy as np
import csv
import json
import os
import pickle
from pathlib import Path
from typing import IO, Callable, Dict, List
import pandas as pd
from dataclasses import asdict
from data_helpers.cls_parsing import unpack
from data_helpers.encoders import AdvancedJsonEncoder

from pyDO3SE.Model_State.model_state_loader import dump_state_to_file

from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.External_State.External_State_Config import FileTypes
from pyDO3SE.Analysis.util import output_log_to_field_data
from pyDO3SE.Analysis.charts import annual_graph
from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.External_State import External_State_Shape
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.error_handling import OutputError


def dump_config_to_string(
    config: Config_Shape,
) -> str:
    return json.dumps(unpack(config), indent=4, cls=AdvancedJsonEncoder)


def dump_config_to_file_json(
    config: Config_Shape,
    target_path: Path,
) -> IO:
    with open(target_path, 'w') as cf:
        json.dump(unpack(config), cf, indent=4, cls=AdvancedJsonEncoder)


def dump_config_to_file_binary(
    config: Config_Shape,
    target_path: Path,
) -> IO:
    """Pickle the model config.

    Parameters
    ----------
    state : Config_Shape
        Config Shape Object
    target_path : Path
        Location to save config

    """
    with open(target_path, 'wb') as configfile:
        pickle.dump(config, configfile)


def dump_output_to_file_csv(
    output_data: List[dict],
    target_path: Path,
) -> IO:
    # TODO: Could filter fields here
    full_logs = pd.DataFrame(output_data)
    full_logs.to_csv(target_path)


def dump_output_to_file_netcdf(
    output_data: List[dict],
    target_path: Path,
) -> IO:
    # TODO: Could filter fields here
    full_logs = pd.DataFrame(output_data)
    ds = full_logs.to_xarray()
    # ds = xr.DataSet.from_dataframe(full_logs)
    xr.save_mfdataset([ds], paths=[target_path])


def merge_netcdf_grid_data(
    output_data: Dict[str, np.ndarray],  # data with shape(x, y ,t) where each element is a dict
    lat_data: List[List[float]],
    lon_data: List[List[float]],
    time_data: List[np.datetime64],
) -> IO:
    variables = output_data[0][0][0].keys()
    X = len(output_data)
    Y = len(output_data[0])
    T = len(output_data[0][0])
    data_vars = {k: (['x', 'y', 'time'], [[[output_data[y][x][t][k] for t in range(T)]
                                           for x in range(X)] for y in range(Y)]) for k in variables}

    # TODO: Could filter fields here
    return xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            lon=(["x", "y"], lon_data),
            lat=(["x", "y"], lat_data),
            time=time_data,
            # reference_time=reference_time,
        ),
        attrs=dict(description="DO3SE outputs"),
    )


def dump_output_to_netcdf_grid(
    output_data: Dict[str, np.ndarray],  # data with shape(x, y ,t) where each element is a dict
    X: int,
    Y: int,
    T: int,
    lat_data: List[List[float]],  # shape(x, y)
    lon_data: List[List[float]],  # shape(x, y)
    time_data: List[np.datetime64],  # shape (t)
    output_fields: List[str] = None,
) -> IO:
    """Dump output to netcdf file.


    Output data should be a dictionary where each value has shape (x, y, t)
    The lat and lon data should have shape (x, y)
    Time length should be t


    Parameters
    ----------
    output_data : Dict[str, np.ndarray]
        Output data should be a dictionary where each value has shape (x, y, t)
    X: int
        gris X size
    Y: int
        gris Y size
    T: int
        gris T size
    lat_data : List[List[float]]
        matrix of lon data for each grid point, shape(x, y)
    lon_data : List[List[float]]
        matrix of lon data for each grid point, shape(x, y)
    time_data : List[np.datetime64]
        matrix of time data with shape(t)
    output_fields : List[str]
        Fields to save

    Returns
    -------
    IO
        [description]
    """
    variables = output_fields or list(output_data.keys())

    # TODO: Pass coords here
    try:
        data_vars = {
            k: (
                ['x', 'y', 'time'],
                output_data[k]
            )
            for k in variables}
        return xr.Dataset(
            data_vars=data_vars,
            coords=dict(
                lon=(["x", "y"], lon_data),
                lat=(["x", "y"], lat_data),
                time=time_data,
                # reference_time=reference_time,
            ),
            attrs=dict(description="DO3SE outputs"),
        )
    except KeyError as e:
        print(e)
        raise OutputError(f"Error saving to netcdf. Check output shape: X:{X} Y: {Y} T: {T}")


def dump_output_to_file_netcdf_grid(
    output_data: Dict[str, np.ndarray],  # data with shape(x, y ,t) where each element is a dict
    X: int,
    Y: int,
    T: int,
    lat_data: List[List[float]],
    lon_data: List[List[float]],
    time_data: List[np.datetime64],
    target_path: Path,
    output_fields: List[str],
) -> IO:
    """Dump the output data to a NETCDF file.

    Output data should be a dictionary where each value has shape (x, y, t)
    The lat and lon data should have shape (x, y)
    Time length should be t


    Parameters
    ----------
    output_data : Dict[str, np.ndarray]
        Output data should be a dictionary where each value has shape (x, y, t)
    X: int
        gris X size
    Y: int
        gris Y size
    T: int
        gris T size
    lon_data : List[List[float]]
        matrix of lon data for each grid point
    time_data : List[np.datetime64]
        matrix of lat data for each grid point
    target_path : Path
        Location to save nc file
    output_fields : List[str]
        Fields to save

    Returns
    -------
    IO
        [description]
    """
    ds = dump_output_to_netcdf_grid(
        output_data,
        X, Y, T,
        lat_data,
        lon_data,
        time_data,
        output_fields,
    )
    xr.save_mfdataset([ds], paths=[target_path])


def export_partial_output(
    output_data: List[dict],
    final_state: Model_State_Shape,
    output_directory: str,
    output_fields: List[str] = [],
    runid: str = None,
    log_notes: List[str] = None,
    dump_output_to_file: Callable[[List[dict], Path], IO] = dump_output_to_file_csv,
):
    """Export the results of a partial run.

    These runs are intened to be followed by another partial run so we should
    export final state ready to use as initial state to next run.

    """
    output_data_file = f"{output_directory}/{runid}_results.csv"
    output_state_file = f"{output_directory}/final_state.json"
    dump_output_to_file(output_data, output_data_file)
    dump_state_to_file(final_state, output_state_file)
    # Save notes
    if log_notes:
        with open(f'{output_directory}/log.txt', 'w') as notesfile:
            notesfile.write("\n".join(log_notes))


def export_output(
    output_data: List[dict],
    final_state: Model_State_Shape,
    external_state: External_State_Shape,
    output_directory: str,
    final_config: Config_Shape,
    output_filename: str = "pyDO3SE_output.csv",
    output_results_only: bool = False,
    fields_to_graph: List[str] = [],
    runid: str = None,
    log_notes: List[str] = None,
    dump_output_to_file: Callable[[List[dict], Path], IO] = dump_output_to_file_csv,
):
    """Output log data to a csv file."""
    # full_logs = pd.DataFrame(output_data)
    row_count = len(external_state.dd)
    external_state_out = pd.DataFrame({k: v[0:row_count] if v is not None else None
                                       for k, v in asdict(external_state).items()})

    os.makedirs(output_directory, exist_ok=True)
    dump_output_to_file(output_data, f"{output_directory}/{output_filename}")

    if output_results_only:
        return

    external_state_out.to_csv(f"{output_directory}/external_data.csv")

    # Save out processed config
    dump_config_to_file_json(final_config, f"{output_directory}/processed_config.json")

    start_day = final_config.Location.start_day
    start_day = start_day if start_day is not None else external_state.dd[0]
    end_day = final_config.Location.end_day
    end_day = end_day if end_day is not None else external_state.dd[-1]

    # Export graphs
    for f in fields_to_graph:
        try:
            field = output_fields_map[f]
        except KeyError:
            raise KeyError(
                f"{f} not in output fields map. Run cli command `pyDO3SE_cli available-outputs")

        data = output_log_to_field_data(output_data, field)

        annual_graph(
            data,
            field,
            output_dir=output_directory,
            chart_id=f'{field.id}-{runid}',
            start_day=start_day,
            end_day=end_day,
            label_x_days=20,
        )

    # Save notes
    if log_notes:
        with open(f'{output_directory}/notes.log', 'w') as notesfile:
            notesfile.write("\n".join(log_notes))

    dump_state_to_file(final_state, f'{output_directory}/final_state.json')


def export_single_hour_output(
    output_data_dir: Path,
    output_logs: List[List[any]],
    final_state: Model_State_Shape,
    output_type: FileTypes = FileTypes.CSV,
):
    os.makedirs(output_data_dir, exist_ok=True)
    if output_type == FileTypes.NETCDF:
        output_data_file = f"{output_data_dir}/output_data.nc"
        dump_output_to_file_netcdf(output_logs, output_data_file)
    else:
        output_data_file = f"{output_data_dir}/output_data.csv"
        dump_output_to_file_csv(output_logs, output_data_file)
    output_state_file = f"{output_data_dir}/output_state.nc"
    dump_state_to_file(final_state, output_state_file)


def get_headers(file):
    with open(file) as f:
        return f.readline().replace('\n', '').replace('\r\n', '').split(',')


def extract_state_val(line_data, field):
    return line_data[field]


def get_final_state_values(input_file_directory, field):
    files = [f for f in os.listdir(input_file_directory) if f.split('.')[-1] == 'csv']
    headers = get_headers(input_file_directory + "/" + files[0])
    # print(headers)
    last_lines = get_last_line_of_files(input_file_directory, files, headers)
    data = [line_data.get(field, None) for line_data in last_lines]
    file_names = [f.split('.')[0] for f in files]
    return dict(zip(file_names, data))


def get_last_line_in_file(filename, headers):
    try:
        with open(filename, 'rb') as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
            last_line = list(csv.reader([f.readline().decode()]))[0]
            last_line_data = {h: d for h, d in zip(headers, last_line)}
            return last_line_data
    except OSError:
        Warning(f'Os Error while getting last line in file {filename}')
        return {h: None for h in headers}


def get_last_line_of_files(dir, files, headers):
    last_lines = [get_last_line_in_file(dir + '/' + f, headers) for f in files]
    return last_lines


def get_and_save_final_state_data(input_file_directory, output_directory, field_name):
    data = get_final_state_values(input_file_directory, field_name)
    name = os.path.basename(input_file_directory)
    print(f"Getting pod for{name}")
    os.makedirs(output_directory, exist_ok=True)
    with open(f"{output_directory}/{name}.json", 'w') as out_file:
        json.dump(data, out_file)
    with open(f"{output_directory}/{name}.csv", 'w', newline='') as out_file:
        spamwriter = csv.writer(out_file, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for k, v in data.items():
            spamwriter.writerow([k, v])


def get_and_save_final_state_data_multi_dir(input_file_directory, output_directory, field_name):
    dirs = [d for d in os.listdir(input_file_directory)
            if os.path.isdir(f"{input_file_directory}/{d}")]
    data = [get_final_state_values(f'{input_file_directory}/{d}', field_name) for d in dirs]
    print(dirs)
    for dir, data_group in zip(dirs, data):
        name = os.path.basename(dir)
        print(f"Getting pod for{name}")
        os.makedirs(output_directory, exist_ok=True)
        with open(f"{output_directory}/{name}.json", 'w') as out_file:
            json.dump(data_group, out_file)
        with open(f"{output_directory}/{name}.csv", 'w', newline='') as out_file:
            spamwriter = csv.writer(out_file, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for k, v in data_group.items():
                spamwriter.writerow([k, v])


def extract_final_results(
    input_directory: Path,
    output_directory: Path,
    field_name: str,
):
    get_and_save_final_state_data(input_directory, output_directory, field_name)

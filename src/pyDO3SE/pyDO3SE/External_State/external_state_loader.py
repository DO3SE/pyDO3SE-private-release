import csv
import os
import numpy as np
from copy import deepcopy
import warnings
import re
from io import StringIO
from typing import Callable, Iterator, List, NamedTuple, Tuple
from functools import partial
from collections import OrderedDict
from pathlib import Path
import pandas as pd
from proflow.Objects.Process import Process
from proflow.ProcessRunnerCls import ProcessRunner

from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.util.error_handling import InputDataError
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.overrides import Main_Overrides

from .External_State_Shape import External_State_Shape, InputField, INPUT_FIELDS
from .External_State_Config import EStateOptions, FileTypes

CSV_DELIMITER = ','
HEADER_REGEX = re.compile(r'^(?P<name>.+?)(\s*,\s*(?P<units>.+))?$')
Coord = Tuple[int, int]


def match_heading_to_field(heading: str, fields: List[InputField]):
    """Matches a heading from the data to a defined field type
    Uses regex HEADER_REGEX to extract name and units from the header string

    e.g. 'P, kPa' matches InputField('P'...),
    """
    if not heading or not heading.replace(' ', ''):
        warnings.warn(f"Found empty heading in input")
        return heading, heading

    # Match information in header string
    match = HEADER_REGEX.match(heading)

    if not match or match is None:
        warnings.warn(f"Failed to find header: {heading}")
        return heading, heading

    match = match.groupdict()
    name, units = match['name'], match['units']

    matching_fields = [f for f in fields if (
        f.name.lower() == name.lower() or name.lower() in f.alt_names)]
    matching_fields_units = [f for f in matching_fields if f.units == units]
    ignore_field = None

    if len(matching_fields) > 1 and len(matching_fields_units) > 0:
        # attempt to match by unit if too many matching fields
        field_out = matching_fields_units[0].name
    elif len(matching_fields) >= 1:
        field_out = matching_fields[0].name
    else:
        ignore_field = heading
        field_out = heading

    return field_out, ignore_field


def match_heading_to_heading(heading: str, headings: List[str]):
    """ matches a heading from the data to a defined field type
    Uses regex HEADER_REGEX to extract name and units from the header string

    e.g. 'P, kPa' matches InputField('P'...),
    """
    # Match information in header string
    match = HEADER_REGEX.match(heading).groupdict()
    name = match['name']

    matching_fields = [h for h in headings if h == name]
    ignore_field = None

    if len(matching_fields) >= 1:
        field_out = matching_fields[0]
    else:
        ignore_field = heading
        field_out = heading

    return field_out, ignore_field


read_header_csv = partial(csv.reader, skipinitialspace=True,
                          delimiter=CSV_DELIMITER, quotechar='"')


def process_csv_data(
    csv_data_raw: str,
    fields: List[InputField] = None,
    has_header_row: bool = True,
    input_headers: List[str] = None,
    row_indexes: List[int] = None,
) -> OrderedDict:
    """processes a csv file into a External State Object

    Inputs:
    csv_data: direct file output from open()
    fields: list of fields to use
    has_header_row: skip first row or use first row as headers depending on headers input
    input_headers: list of headers to use if not using first row or fields

    """
    if not has_header_row and input_headers is None:
        raise Exception('Missing headers')

    # 1. get headers from data or input
    data_headings = next(read_header_csv(csv_data_raw)) \
        if has_header_row else input_headers

    # 3. remove whitespace in headers
    header_row = [h.strip() for h in data_headings] if data_headings else None
    # 4. match headings from data to requested fields or input_headings list
    match_headings_fn = \
        partial(match_heading_to_field, fields=fields) \
        if fields is not None else \
        partial(match_heading_to_heading, headings=input_headers)
    header_names, ignore_fields = list(zip(*[match_headings_fn(h) for h in header_row]))
    ignore_fields = [f for f in ignore_fields if f is not None]
    if len(ignore_fields) > 0:
        warnings.warn(UserWarning(
            f'\nThe following ext data columns have been ignored \n {ignore_fields}'))

    # 1. read csv data
    data = pd.read_csv(csv_data_raw, skipinitialspace=True,
                       delimiter=CSV_DELIMITER, quotechar='"',
                       skiprows=lambda x: x not in row_indexes if row_indexes else False,
                       names=header_names,
                       #    header=0
                       )
    data_filtered = data.drop(ignore_fields, axis=1)
    data_remove_nans = data_filtered.replace({np.nan: None})
    data_out = data_remove_nans.to_dict('list')
    return data_out


def convert_csv_to_external_state(
    csv_data: pd.DataFrame,
    fields_data: List[InputField] = INPUT_FIELDS,
) -> External_State_Shape:
    # TODO: WIP
    data_headings = csv_data.columns
    header_row = [h.strip() for h in data_headings] if data_headings else None
    match_headings_fn = \
        partial(match_heading_to_field, fields=fields_data)
    header_names, ignore_fields = list(zip(*[match_headings_fn(h) for h in header_row]))
    ignore_fields = [f for f in ignore_fields if f is not None]
    if len(ignore_fields) > 0:
        warnings.warn(UserWarning(
            f'\nThe following ext data columns have been ignored \n {ignore_fields}'))
    d = {kk: csv_data[k].values for k, kk in zip(data_headings, header_names)}
    external_state_object = External_State_Shape(**d)
    return external_state_object


def load_external_state_csv(
    external_state_file_location: str,
    has_header_row: bool = True,
    row_indexes: List[int] = None,
    logger: Callable[[str, str], None] = print,
) -> External_State_Shape:
    try:
        with open(external_state_file_location, encoding='utf-8-sig') as external_state_file:
            logger(f"Loading external data from {external_state_file_location}")
            read_data = external_state_file.read()
            if os.path.basename(external_state_file_location).split('.')[-1] != 'csv':
                raise ValueError(f"Input data is not a csv file")
            # get the data associated with each requested field
            fields_data = [fieldData for fieldData in INPUT_FIELDS]
            # fields_data = [fieldData for fieldData in INPUT_FIELDS if fieldData.name in fields]
            external_state_data = process_csv_data(
                StringIO(read_data),
                fields=fields_data,
                has_header_row=has_header_row,
                row_indexes=row_indexes,
            )
            external_state_object = External_State_Shape(**external_state_data)

            return external_state_object
    except UnicodeDecodeError as e:
        print(f"Failed to import data due to unicode error")
        # return External_State_Shape()
        raise e
    except FileNotFoundError as e:
        full_dir= os.path.abspath(external_state_file_location)
        print(f"Failed to import data from '{full_dir}' due to file not found")
        # list all files in the directory
        files = os.listdir(os.path.dirname(full_dir))
        print(f"Files in directory: {files}")
        raise e
    except Exception as e:
        raise e


def extract_cell_data_from_netcdf(
    data,  # Xarray dataset
    xi: int,
    yi: int,
    T: int,
    time_key: str,
    variable_map: dict,
    index_counts: List[int],
    use_dask: bool,
) -> dict:
    """Extract data for a single cell from a netcdf file.

    We handle the issue that accessing a none dask array with .item is quicker than indexing
    but dask arrays do not allow accessing with .item.

    # To manage multi level variables we assume the second index is for the level.

    Parameters
    ----------
    data : [type]
        The input Xarray dataset
    xi : int
        grid x index
    yi : int
        grid y index
    T : int
        Size of Time axis
    time_key: str
        Key for time variable in data
    variable_map : dict(Should be ordered)
        map of netcdf variables to DO3SE variabls
    index_counts : List[int]
        The size of the shape for each variable.

    use_dask : bool
        if true use Xarray Dask

    Returns
    -------
    dict
        a dictionary where each value is a list of the variable values per hour.

    """
    vert_index = 0  # We assume we want the data at layer 0

    if use_dask:
        # We can't currently use .item when using dask
        # TODO: This assumes that the indexing is (T, X, Y)
        data_out = dict([
            (ki, data[kj][:, vert_index, xi, yi].values) if icount == 4 else
            (ki, data[kj][:, xi, yi].values) if icount == 3 else
            (ki, data[kj][:].values)
            for icount, (ki, kj) in zip(index_counts, variable_map.items())])

        data_base = {
            'dd': data['dd'].values,
            'hr': data['hr'].values,
            'time': [data[time_key][t].values.astype(str)[0:19] for t in range(T)]
        }
    else:
        data_base = {
            'dd': [data['dd'].item(t) for t in range(T)],
            'hr': [data['hr'].item(t) for t in range(T)],
            'time': [data[time_key][t].values.astype(str)[0:19] for t in range(T)]
        }
        if T > 1:
            # TODO: Make getting multiple time indexes more efficient
            data_out = dict([
                (ki, [data[kj].item((t, vert_index, xi, yi)) for t in range(T)]) if icount == 4 else
                (ki, [data[kj].item((t, xi, yi)) for t in range(T)]) if icount == 3 else
                (ki, [data[kj].item(t) for t in range(T)])
                for icount, (ki, kj) in zip(index_counts, variable_map.items())])
        else:
            t = 0
            data_out = dict([
                (ki, [data[kj].item((t, vert_index, xi, yi))]) if icount == 4 else
                (ki, [data[kj].item((t, xi, yi))]) if icount == 3 else
                (ki, [data[kj].item(t)])
                for icount, (ki, kj) in zip(index_counts, variable_map.items())])

    return {**data_base, **data_out}


def load_external_state_netcdf(
    data_location: Path,
    coords: List[Coord],
    variable_map: dict = None,
    multi_file_data: bool = False,
    preprocess_map: dict = {},
    zero_year: int = None,
    data_filter: str = '',
    chunks: dict = None,
    parallel: bool = False,
    **kwargs,
) -> Iterator[External_State_Shape]:
    """Load NetCDF data to external state object.

    TODO: Replace variable map with match functions above.

    Parameters
    ----------
    data_location : Path
        Path to NetCDF file(s). Should be directory if multi_file_data is true
    coords: List[Coord]
        coordinate to rnu
    variable_map : dict, optional
        mapping of DO3SE variables to input variables, by default None
    multi_file_data: bool, defaul False
        If true then path is to a directory of NetCDF files
    preprocess_map: dict, default {}
        dictionary mapping output variables to a function that modifies the variable
    zero_year: int, optional
        DOY (where year == zero_year + 1) == DOY + 365
    data_filter: str, optional
        When using multi file mode can use filename filter.
    file_suffix: str, optional
        Allow overrideing input file suffix. This is after chosen filetype
    chunks: dict, optional
        Chunk options to pass to Xarray
    parallel: bool
        if true use parallel loading

    Returns
    -------
    Iterator[External_State_Shape]
        External state shape for model runs

    """
    variable_map = deepcopy(variable_map)  # Copy so we can pop elements

    # Get time and meta data
    assert "time" in variable_map, f"Must include \"time\" in variable_map.json for data: \"{data_location}\""
    assert "_SHAPE" in variable_map, f"Must include \"_SHAPE\" in variable_map.json for data: \"{data_location}\""
    time_key = variable_map.pop('time')
    shape_key = variable_map.pop('_SHAPE')
    use_dask = parallel or chunks is not None or multi_file_data
    # TODO: Handle dropping of layers we don't need
    if multi_file_data:
        # file_filter = f"*{data_filter}*{file_suffix}" if data_filter or file_suffix else '*'
        file_filter = f"*{data_filter}*" if data_filter else "*"
        try:
            data = xr.open_mfdataset(f'{data_location}/{file_filter}',
                                     engine="netcdf4",
                                     parallel=parallel,
                                     **kwargs,
                                     )
        except:
            raise ValueError(
                f'Failed to load external data files in multi file mode: {data_location} using file_filter: {file_filter}')

    else:
        try:
            data = xr.open_dataset(data_location, chunks=chunks, engine="netcdf4", **kwargs)
        except Exception as e:
            raise ValueError(
                f'Failed to load external data files in single file mode: {data_location}') from e

    assert time_key in list(data.keys()) or time_key in list(
        data.dims) or time_key in list(data.coords), \
        f'Time key: "{time_key}" not in data. Got keys: {data.keys()} '

    data_processed = data

    try:
        index_counts = [len(data_processed[k].shape) for k in variable_map.values()]
    except KeyError as e:
        raise InputDataError(
            f"Input data missing required key: {e}\n{','.join(list(data_processed.keys()))}")

    data_processed = data_processed.assign(hr=lambda d: getattr(d, time_key).dt.hour)
    data_processed = data_processed.assign(year=lambda d: getattr(d, time_key).dt.year)
    data_processed = data_processed.assign(
        dd=lambda d: getattr(d, time_key).dt.strftime('%j').astype(int) + (d.year.astype(int) - zero_year) * 365)

    T = data_processed[time_key].shape[0]

    # Check that the first dim is the time key
    # This assumes that we have Ts_C as an input
    data_processed[time_key].dims[0] == data_processed[shape_key].dims[0]

    for xi, yi in coords:
        data_out = extract_cell_data_from_netcdf(
            data_processed,
            xi, yi, T,
            time_key,
            variable_map,
            index_counts,
            use_dask=use_dask,
        )
        for k, func in preprocess_map.items():
            data_out[k] = [eval(func)(d) for d in data_out[k]]

        external_state_object = External_State_Shape(**data_out)
        yield external_state_object


class EStateDates(NamedTuple):
    start_day: int = None
    end_day: int = None
    start_date: str = None
    end_date: str = None
    row_count: int = None
    time_string: str = None
    hours: List[int] = None


def run_init_processes_on_e_state(
    external_state_data: External_State_Shape,
    config: Config_Shape,
    init_processes: List[Process],
    overrides: Main_Overrides = Main_Overrides(),
) -> External_State_Shape:
    """Run process runner on external state.

    Parameters
    ----------
    external_state_data : External_State_Shape
        _description_
    config : Config_Shape
        _description_
    init_processes : List[Process]
        _description_
    overrides : Main_Overrides, optional
        _description_, by default Main_Overrides()

    Returns
    -------
    External_State_Shape
        _description_
    """
    process_runner = ProcessRunner(config, DEBUG_MODE=overrides.debug)

    external_state = process_runner.run_processes(
        init_processes,
        external_state_data,
    )

    return external_state


def get_date_bounds_from_ext_data(
    external_state: External_State_Shape,
) -> EStateDates:
    """Get the time boundaries from external data.

    Parameters
    ----------
    external_state : External_State_Shape
        External state data for single cell

    Returns
    -------
    EStateDates
        Bounding dates from external data

    """
    # TODO: Can this handle when data wraps around end of year?
    start_day = external_state.dd[0]
    end_day = external_state.dd[-1]
    row_count = len(external_state.dd)

    if external_state.time is not None and external_state.time[0]:
        start_date = external_state.time[0]
        end_date = external_state.time[-1]
        time_data = pd.date_range(start_date, periods=len(
            external_state.time), freq="1H")
        # time_data = external_state.time
        time_string = f'{time_data[0].year}-{str(time_data[0].month).zfill(2)}-{str(time_data[0].day).zfill(2)}_{str(time_data[0].hour).zfill(2)}'

    else:
        start_date = None
        end_date = None
        time_string = None
    return EStateDates(
        start_day=start_day,
        end_day=end_day,
        start_date=start_date,
        end_date=end_date,
        row_count=row_count,
        time_string=time_string,
        hours=external_state.hr,
    )


def get_lat_and_lon_from_external_data(
    external_data_file_path: Path,
    external_state_options: EStateOptions,
    **kwargs,
) -> Tuple[float, float]:
    """Get the lat and lon from external data.

    Parameters
    ----------
    external_state : External_State_Shape
        External state data for single cell

    Returns
    -------
    Tuple[float, float]
        lat and lon from external data

    """
    print("External state path", external_data_file_path)
    _external_state_options = external_state_options._replace(
        **kwargs,
    )
    print("External state options", _external_state_options)
    variable_map = deepcopy(_external_state_options.variable_map)  # Copy so we can pop elements

    # Get time and meta data
    assert "lat" in variable_map, f"Must include \"lat\" in variable_map.json for data: \"{external_data_file_path}\""
    assert "lon" in variable_map, f"Must include \"lon\" in variable_map.json for data: \"{external_data_file_path}\""
    lat_key = variable_map.pop('lat')
    lon_key = variable_map.pop('lon')

    if _external_state_options.multi_file_data:
        # file_filter = f"*{data_filter}*{file_suffix}" if data_filter or file_suffix else '*'
        try:
            first_file = os.listdir(external_data_file_path)[0]
            data = xr.open_dataset(f'{external_data_file_path}/{first_file}',
                                   engine="netcdf4",
                                   **kwargs,
                                   )
        except:
            raise ValueError(
                f'Failed to load external data files in multi file mode: {external_data_file_path} using file_filter: {first_file}')

    else:
        try:
            data = xr.open_dataset(external_data_file_path, engine="netcdf4", **kwargs)
        except Exception as e:
            print(e)
            raise ValueError(
                f'Failed to load external data files in single file mode: {external_data_file_path}')
    lat = data[lat_key].values
    lon = data[lon_key].values
    return [lat, lon]


def load_external_state(
    external_state_file_location: Path,
    grid_coords: List[Tuple[int, int]] = None,
    external_state_options: EStateOptions = EStateOptions(),
    logger: Callable[[str, str], None] = print,
    **kwargs,
) -> Iterator[External_State_Shape]:
    """loads external state from a file.

    Checks against the fields in external_state_shape.

    Parameters
    ----------
    external_state_file_location : Path
        Path to external data file
    grid_coords: List[Tuple[int, int]]
        Grid coordinates to load
    external_state_options: EStateOptions
        Options for loading external data

    **kwargs
        Passed to EStateOptions constructor

    Returns
    -------
    Iterator[External_State_Shape]
        External state shape for model runs
        External state meta data

    """
    # _external_state_options = external_state_options or EStateOptions(**kwargs)
    _external_state_options = external_state_options._replace(
        **kwargs,
    )
    if _external_state_options.file_type == FileTypes.CSV:
        if grid_coords is not None:
            raise NotImplementedError("Grid workflow not setup for csv")
        if (_external_state_options.preprocess_map and len(_external_state_options.preprocess_map.keys())):
            raise NotImplementedError("preprocess not implemented for csv")

        try:
            ext_data = load_external_state_csv(
                external_state_file_location,
                _external_state_options.has_header_row,
                _external_state_options.row_indexes,
                logger=logger,
            )
        except Exception:
            raise InputDataError(
                f"Failed to load external data from \"{external_state_file_location}\"")

        yield ext_data

    elif _external_state_options.file_type == FileTypes.NETCDF:
        if grid_coords is None or _external_state_options.variable_map is None:
            raise ValueError("Must supply coords and variable map for NetCDF loading.")

        try:
            loaded_data = load_external_state_netcdf(
                external_state_file_location,
                grid_coords,
                _external_state_options.variable_map,
                _external_state_options.multi_file_data,
                _external_state_options.preprocess_map,
                _external_state_options.zero_year,
                data_filter=_external_state_options.data_filter,
                **_external_state_options.netcdf_loader_kwargs,
            )
        except Exception:
            raise InputDataError(
                f"Failed to load external data from \"{external_state_file_location}\"")

        # loaded data is flattened (x, y) array
        # for each coordinate we yield external data
        for _ in grid_coords:
            ext_data = next(loaded_data)
            yield ext_data
    else:
        raise Exception(
            f'Invalid data file type. Must be csv or netcdf, got {_external_state_options.file_type}')

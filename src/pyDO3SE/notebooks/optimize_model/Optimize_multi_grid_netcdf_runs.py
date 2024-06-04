# %%
import os
from math import isclose
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from timeit import repeat
import warnings
import cProfile
import pstats
from pstats import SortKey
import xarray as xr
# %%

%load_ext autoreload
%autoreload 2

from pyDO3SE.util.loader import csv_loader, json_loader
from pyDO3SE.main import main_hour
from pyDO3SE.External_State.external_state_loader import FileTypes, load_external_state_netcdf
from pyDO3SE.util.logger import Logger
# %%
# %%
log_level = 0
logger_main = Logger(log_level)
# Inputs
project_dir = 'examples/net_cdf/wrfchem'
multi_file_netcdf = False
config_dir = f"{project_dir}/configs"
inputs_dir = f"{project_dir}/inputs"

base_config_file = f"{project_dir}/base_config.json"
base_state_file = f"{project_dir}/base_state.json"
e_state_overrides_file_path = f"{project_dir}/e_state_overrides.nc"

variable_map_path = f"{project_dir}/variable_map.json"
preprocess_map_path = f"{project_dir}/preprocess_map.json"
e_state_overrides_field_map_path = f"{project_dir}/e_state_overrides_field_map.json"

coordinates_dir = f"{project_dir}/coords"

configs = os.listdir(config_dir)
config_file_path = configs[0]
runid = 99
logger_main(f'== Found {len(configs)} configs to run =====')
config_name = '.'.join(config_file_path.split('.')[:-1])
full_config_path = f"{project_dir}/configs/{config_file_path}"
run_dir = f"{project_dir}/runs/{runid}/{config_name}"
os.makedirs(run_dir, exist_ok=True)
log_path = f"{run_dir}/run.log"
logger = Logger(log_level, log_path)
logger(f'== Running config {config_file_path} ==')

coordinates_file_path = f"{coordinates_dir}/{config_name}.csv"

# Processed file locations
live_state_dir = f"{run_dir}/current_state"
prev_state_dir = f"{run_dir}/prev_state"

output_data_dir = f"{run_dir}/outputs_grid"
processed_configs_dir = f"{run_dir}/processed_configs"

# loaded data
variable_map = json_loader(variable_map_path)
preprocess_map = json_loader(preprocess_map_path)
e_state_overrides_field_map = json_loader(e_state_overrides_field_map_path)

# TODO: Get grid from args or input file
grid_coords = [[int(i['x']), int(i['y'])] for i in csv_loader(coordinates_file_path)]
f = sorted(os.listdir(inputs_dir))[0]
input_data_file = f"{inputs_dir}/{f}"
state_out_path = live_state_dir
previous_hour_state_path = live_state_dir
output_fields = ['pody']


# %%

coords = [
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
    (2, 2),
]
# logger= lambda*args, **kwargs: lambda *args, **kwargs: None


def run():

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        next(main_hour(
            processed_config_dir=processed_configs_dir,
            external_data_row_path=input_data_file,
            previous_hour_state_path=previous_hour_state_path,
            output_data_dir=output_data_dir,
            output_fields=output_fields,
            external_file_type=FileTypes.NETCDF,
            grid_coords=coords,
            netcdf_variable_map=variable_map,
            met_preprocess_map=preprocess_map,
            multi_file_netcdf=multi_file_netcdf,
            output_to_netcdf=True,
            state_out_path=state_out_path,
            logger=logger,
        ))


run()
# %%


t = min(repeat(run, number=10, repeat=10))
t
# %%


def run():

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        list(main_hour(
            processed_config_dir=processed_configs_dir,
            external_data_row_path=input_data_file,
            previous_hour_state_path=previous_hour_state_path,
            output_data_dir=None,
            output_fields=output_fields,
            external_file_type=FileTypes.NETCDF,
            grid_coords=grid_coords,
            netcdf_variable_map=variable_map,
            met_preprocess_map=preprocess_map,
            multi_file_netcdf=multi_file_netcdf,
            output_to_netcdf=True,
            state_out_path=state_out_path,
            logger=logger,
        ))


run()
# %%


t = min(repeat(run, number=10, repeat=10))
t

# %%
profiler_output = 'profiler_output'
cProfile.run('run()', profiler_output)

# %%
# %%
p = pstats.Stats(profiler_output)
p.sort_stats(SortKey.CUMULATIVE).print_stats(15, 'numpy')
# %%
p.print_callees(30, 'numpy')
# %%

variable_map = json_loader(variable_map_path)

print(variable_map)
print(variable_map['time'])


def runb():
    data_location = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'
    coords = [[0, 0]]
    data = xr.open_dataset(data_location)
    data_processed = data
    variable_map_b = {**variable_map}
    # Get time and meta data
    time_key = variable_map_b.pop('time')
    data_processed = data_processed.assign(hr=lambda d: getattr(d, time_key).dt.hour)
    data_processed = data_processed.assign(
        dd=lambda d: getattr(d, time_key).dt.strftime('%j').astype(int))

    # We assume a variable is multi level if it has a shape size of 4
    index_counts = [len(data_processed[k].shape) for k in variable_map_b.values()]

    for xi, yi in coords:
        # To manage multi level variables we assume the second index is for the level.
        data_out = dict([
            (ki, data_processed[kj][:, 0, xi, yi].values) if icount == 4 else
            (ki, data_processed[kj][:, xi, yi].values) if icount == 3 else
            (ki, data_processed[kj][:].values)
            for icount, (ki, kj) in zip(index_counts, variable_map_b.items())])

        for k, func in preprocess_map.items():
            data_out[k] = [func(d) for d in data_out[k]]

        data_out['dd'] = data_processed.dd.values
        data_out['hr'] = data_processed.hr.values
        # external_state_object = External_State_Shape(**data_out)


t = min(repeat(runb, number=10, repeat=10))
t
# %%
profiler_output = 'profiler_output'
cProfile.run('runb()', profiler_output)
p = pstats.Stats(profiler_output)
p.sort_stats(SortKey.CUMULATIVE).print_stats(20)
# %%

inputs_dir = f"{project_dir}/inputs"
data_file = f'{inputs_dir}/wrfchem_demo_out_01_10_01.nc'
ds = xr.open_dataset(data_file)
ds
# %%
ds.isel(south_north=0, west_east=0, bottom_top=0)['td_2m']

# %%
ds['td_2m'].item((0,1,1))
# %%

coords = [
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
    (2, 2),
]
print(variable_map)
variable_map = {
    'time': 'XTIME',
    'Ts_C': 'td_2m',
    # 'P': 'pres',
    # 'PAR': 'SWDOWN',
    # 'precip': 'RAINNC',
    # 'RH': 'rh',
    # 'u': 'wspeed',
    # 'O3': 'o3',
    # 'Hd': 'HFX_FORCE',
    # 'snow_depth': 'SNOWH',
}
from copy import deepcopy

def extract_cell_data_from_netcdf(
    data, # Xarray dataset
    xi: int,
    yi: int,
    variable_map: dict,
) -> dict:
    T = 0
    index_counts = [len(data[k].shape) for k in variable_map.values()]
    # data_out = dict([
    #         (ki, data[kj][:, 0, xi, yi].values) if icount == 4 else
    #         (ki, data[kj][:, xi, yi].values) if icount == 3 else
    #         (ki, data[kj][:].values)
    #             for icount, (ki, kj) in zip(index_counts, variable_map.items())])
    data_out = dict([
            (ki, [data[kj].item((T, 0, xi, yi))]) if icount == 4 else
            (ki, [data[kj].item((T, xi, yi))]) if icount == 3 else
            (ki, [data[kj].item(T)])
                for icount, (ki, kj) in zip(index_counts, variable_map.items())])
    return data_out

data = xr.open_dataset(data_file)
def test_load_netcdf():
    variable_map_copy = deepcopy(variable_map) # Copy so we can pop elements

    # Get time and meta data
    time_key = variable_map_copy.pop('time')
    data_processed = data
    # data_processed = data[list(variable_map_copy.values())]

    # data_processed = data_processed.assign(hr=lambda d: getattr(d, time_key).dt.hour)
    # data_processed = data_processed.assign(
    #     dd=lambda d: getattr(d, time_key).dt.strftime('%j').astype(int))

    # We assume a variable is multi level if it has a shape size of 4
    # index_counts = [len(data_processed[k].shape) for k in variable_map_copy.values()]

    xi, yi = coords[0]
    # To manage multi level variables we assume the second index is for the level.
    data_out = extract_cell_data_from_netcdf(
            data_processed,
            xi,yi,
            variable_map,
        )
    # data_out = dict([
    #     (ki, data_processed[kj][:, 0, xi, yi].values.tolist()) if icount == 4 else
    #     (ki, data_processed[kj][:, xi, yi].values.tolist()) if icount == 3 else
    #     (ki, data_processed[kj][:].values.tolist())
    #         for icount, (ki, kj) in zip(index_counts, variable_map_copy.items())])

    # for k, func in preprocess_map.items():
    #     data_out[k] = [func(d) for d in data_out[k]]

    # data_out['dd'] = data_processed.dd.values
    # data_out['hr'] = data_processed.hr.values
    # print(data_out)
    # assert isclose(data_out['Ts_C'], expected['Ts_C'])
    # expected = {
    #     'Ts_C': [294.94],
    #     'dd': [274],
    #     'hr': [1],
    # }
    # print(expected)
    return data_out
print(test_load_netcdf())
# Best is 2.75
t = min(repeat(test_load_netcdf, number=2, repeat=3))
t
# %%
data = xr.open_dataset(data_file)
variable_map_copy = deepcopy(variable_map) # Copy so we can pop elements
data_processed = data[list(variable_map_copy.values())]
data_processed

# %%
data_processed.to_netcdf('temp.nc')
# %%
data_file_small = 'temp.nc'
def extract_cell_data_from_netcdf(
    data, # Xarray dataset
    xi: int,
    yi: int,
    variable_map: dict,
) -> dict:
    T = 0
    data_sel = data# data.isel(south_north=xi, west_east=yi)
    # data_sel = data.isel(south_north=xi, west_east=yi, bottom_top=0)
    index_counts = [len(data[k].shape) for k in variable_map.values()]
    # data_out = dict([
    #         (ki, data[kj][:, 0, xi, yi].values) if icount == 4 else
    #         (ki, data[kj][:, xi, yi].values) if icount == 3 else
    #         (ki, data[kj][:].values)
    #             for icount, (ki, kj) in zip(index_counts, variable_map.items())])
    data_out = dict([
            (ki, [data[kj].item((T, 0, xi, yi))]) if icount == 4 else
            (ki, [data[kj].item((T, xi, yi))]) if icount == 3 else
            (ki, [data[kj].item(T)])
                for icount, (ki, kj) in zip(index_counts, variable_map.items())])
    # data_out = dict([(ki, data_sel[kj].item(0)) for ki, kj in variable_map.items()])
    # data_out = data['td_2m'].item(0) # data_sel[list(variable_map.values())].to_dataframe().to_dict('records')
    return data_out

def run():
    data = xr.open_dataset(data_file)
    variable_map_copy = deepcopy(variable_map) # Copy so we can pop elements
    data_processed = data[list(variable_map_copy.values())]
    # data_file_b = 'tempb.nc'
    # data_processed.to_netcdf(data_file_b)
    # data_small = xr.open_dataset(data_file_b)
    for i in range(48):
        extract_cell_data_from_netcdf(data_processed, 0,i, variable_map)
print(run())
# Best is 3.8
t = min(repeat(run, number=3, repeat=5))
t
# %%

# %%
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
# %%
ds = xr.open_mfdataset(
  "examples/net_cdf/single_file_hour/inputs/*",
  # "examples/net_cdf/multi_file_range/inputs/demo_wrf_2017-10_*",
  engine="netcdf4",
  concat_dim="Time",
#   concat_dim=["X", "Y", "Time"],
  combine="nested",
  parallel=True,
)
ds


# %%
# Method A. Split by date filter then combine by time dim
import re
def get_combined_ds_A(loc, re_filter, time_key="Time"):
    files_list = sorted(os.listdir(loc))
    files_list_groups = [
        getattr(re.search(re_filter, f), 'group', lambda: None)()
        for f in files_list
    ]
    groups = sorted(set(files_list_groups))
    files_list_grouped= [[f"{loc}/{f}" for f in files_list if k in f] for k in groups]
    dss = [xr.open_mfdataset(
        files,
        engine="netcdf4",
        parallel=True,
    ) for files in files_list_grouped]
    ds = xr.combine_nested(dss, concat_dim=time_key)
    return ds

get_combined_ds_A(loc="examples/net_cdf/multi_file_range/inputs", re_filter='[0-9]{4}-[0-9]{2}')
get_combined_ds_A(loc="examples/net_cdf/single_file_range/inputs", re_filter = '[0-9]{4}-[0-9]{2}-[0-9]{2}')
get_combined_ds_A(loc="examples/net_cdf/single_file_hour/inputs", re_filter = '[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{2}')


# %%
# Method B - Combine by variable then merge by time - Only works if split and named by variable
def get_combined_ds_b(loc, time_key="Time"):
    vars = list(xr.open_mfdataset(f"{loc}/*", combine="nested" ,compat="override").keys())
    files_list = sorted(os.listdir(loc))
    files_list_grouped= [[f"{loc}/{f}" for f in files_list if k in f] for k in vars]
    dss = [xr.open_mfdataset(
        files,
        engine="netcdf4",
        concat_dim=time_key,
        combine="nested",
        parallel=True,
    ) for files in files_list_grouped]
    ds = xr.combine_by_coords(dss)
    return ds
# Only works when nc split by variable
get_combined_ds_b(loc="examples/net_cdf/multi_file_range/inputs")
# get_combined_ds_b(loc="examples/net_cdf/single_file_range/inputs")
# get_combined_ds_b(loc="examples/net_cdf/single_file_hour/inputs")

# %%
# Method C - Load all in MF - Only works if not split by variable
def get_combined_ds_c(loc, time_key="Time"):
    ds = xr.open_mfdataset(
        f"{loc}/*",
        # "examples/net_cdf/multi_file_range/inputs/demo_wrf_2017-10_*",
        engine="netcdf4",
        concat_dim=time_key,
        #   concat_dim=["X", "Y", "Time"],
        combine="nested",
        parallel=True,
    )
    return ds

# get_combined_ds_c(loc="examples/net_cdf/multi_file_range/inputs")
get_combined_ds_c(loc="examples/net_cdf/single_file_range/inputs")
get_combined_ds_c(loc="examples/net_cdf/single_file_hour/inputs")

# %%
# Catch all method A
def get_combined_ds(loc, vars=None, re_filter=None):
    if vars:
        files_list = sorted(os.listdir(loc))
        files_list_grouped= [[f"{loc}/{f}" for f in files_list if k in f] for k in vars]
        dss = [xr.open_mfdataset(
            files,
            engine="netcdf4",
            concat_dim="Time",
            combine="nested",
            parallel=True,
        ) for files in files_list_grouped]
        ds = xr.combine_by_coords(dss)
    elif re_filter:
        files_list = sorted(os.listdir(loc))
        files_list_groups = [
            getattr(re.search(re_filter, f), 'group', lambda: None)()
            for f in files_list
        ]
        groups = sorted(set(files_list_groups))
        files_list_grouped= [[f"{loc}/{f}" for f in files_list if k in f] for k in groups]
        dss = [xr.open_mfdataset(
            files,
            engine="netcdf4",
            parallel=True,
        ) for files in files_list_grouped]
        ds = xr.combine_nested(dss, concat_dim="Time")
    else:
        raise ValueError("Must supply vars or a re_filter")

    return ds
vars = list(xr.open_mfdataset("examples/net_cdf/multi_file_range/inputs/*", combine="nested" ,compat="override").keys())
get_combined_ds(loc="examples/net_cdf/multi_file_range/inputs", vars=vars)
get_combined_ds(loc="examples/net_cdf/single_file_range/inputs", re_filter = '[0-9]{4}-[0-9]{2}-[0-9]{2}')
get_combined_ds(loc="examples/net_cdf/single_file_hour/inputs", re_filter = '[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{2}')

# %%
#  === Compare times for A, B and C
from timeit import repeat
grid_coords= [
    [0,1],
    [2,1],
    [2,0],
]
def check_method(fn):
    def _inner():
        ds = fn()
        variable_map= {
            "Ts_C": "td_2m",
            "O3": "o3",
            "u": "wspeed",
        }
        for xi, yi in grid_coords:
            data_out = dict([(ki, ds[kj][:, xi, yi].values)
                for ki, kj in variable_map.items()])
    return _inner

# %%
# multi file range

opt_a = lambda: get_combined_ds_A(loc="examples/net_cdf/multi_file_range/inputs", re_filter='[0-9]{4}-[0-9]{2}')
opt_b = lambda: get_combined_ds_b(loc="examples/net_cdf/multi_file_range/inputs")

t_a = min(repeat(check_method(opt_a), number=5, repeat=10))
t_b = min(repeat(check_method(opt_b), number=5, repeat=10))
print(t_a, t_b)

# %%
# single file range
opt_a = lambda: get_combined_ds_A(loc="examples/net_cdf/single_file_range/inputs", re_filter = '[0-9]{4}-[0-9]{2}-[0-9]{2}')
opt_b = lambda: get_combined_ds_c(loc="examples/net_cdf/single_file_range/inputs")

t_a = min(repeat(check_method(opt_a), number=5, repeat=10))
t_b = min(repeat(check_method(opt_b), number=5, repeat=10))
print(t_a, t_b)

# # %%
# # single file hour
# # THIS IS VERY SLOW!
# # opt_a = lambda: get_combined_ds_A(loc="examples/net_cdf/single_file_hour/inputs", re_filter = '[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{2}')
# # opt_b = lambda: get_combined_ds_c(loc="examples/net_cdf/single_file_hour/inputs")

# # t_a = min(repeat(check_method(opt_a), number=5, repeat=10))
# # t_b = min(repeat(check_method(opt_b), number=5, repeat=10))
# # print(t_a, t_b)

# # %%
# # Work out the fastest way of opening up hourly data

# ds = xr.open_mfdataset(
#   "examples/net_cdf/single_file_hour/inputs/*",
#   # "examples/net_cdf/multi_file_range/inputs/demo_wrf_2017-10_*",
#   engine="netcdf4",
#   concat_dim="Time",
# #   concat_dim=["X", "Y", "Time"],
#   combine="nested",
#   parallel=True,
#   chunks={'x': 1, 'y': 1, "Time": 24},
#     coords="minimal",
# )
# ds


# # %%
# from netCDF4 import Dataset, MFDataset
# # %%
# rootgrp = Dataset("examples/net_cdf/single_file_hour/inputs/demo_wrf_2017-12-27-00-00-00", "r", format="NETCDF4")
# print(rootgrp.data_model)
# # rootgrp.close()
# # %%
# rootgrp.dimensions
# # %%
# f = MFDataset("examples/net_cdf/single_file_hour/inputs/*nc",aggdim="Time")
# f
# # %%
# import dask
# from glob import glob
# paths = glob("examples/net_cdf/single_file_hour/inputs/*nc")
# paths
# # %%
# open_dataset= xr.open_dataset
# datasets = [dask.delayed(open_dataset)(p, "r", format="NETCDF4") for p in paths]
# closers = [dask.delayed(getattr)(ds, "_close") for ds in datasets]
# datasets, closers = dask.compute(datasets, closers)
# # %%
# # combined = _nested_combine(
# #     datasets,
# #     concat_dims=concat_dim,
# #     compat=compat,
# #     data_vars=data_vars,
# #     coords=coords,
# #     ids=ids,
# #     join=join,
# #     combine_attrs=combine_attrs,
# # )
# combined_ds = xr.concat(datasets, 'Time')
# combined_ds

# # %%
# from xarray.core.combine import _infer_concat_order_from_positions
# combined_ids_paths = _infer_concat_order_from_positions(paths)
# ids, paths = (
#     list(combined_ids_paths.keys()),
#     list(combined_ids_paths.values()),
# )
# open_kwargs = dict(engine="netcdf4", chunks={"Time": 24})


# open_dataset= xr.open_dataset
# # wrap the open_dataset, getattr, and preprocess with delayed
# open_ = dask.delayed(open_dataset)
# getattr_ = dask.delayed(getattr)

# datasets = [open_(p, **open_kwargs) for p in paths]
# closers = [getattr_(ds, "_close") for ds in datasets]

# datasets, closers = dask.compute(datasets, closers)
# # %%
# concat_dim= ""
# engine="netcdf4"
# concat_dim=["Time"]
# #   concat_dim=["X", "Y", "Time"]
# combine="nested"
# parallel=True
# chunks={'x': 1, 'y': 1, "Time": 24}
# compat = "no_conflicts"
# data_vars = "all"
# coords="minimal"
# combined = xr.core.combine._nested_combine(
#                 datasets,
#                 concat_dims=concat_dim,
#                 compat=compat,
#                 data_vars=data_vars,
#                 coords=coords,
#                 ids=ids,
#                 # join=join,
#                 # combine_attrs=combine_attrs,
#             )
# combined
# # %%
# from xarray.core.combine import _infer_concat_order_from_positions
# def open_manually():
#     paths = glob("examples/net_cdf/single_file_range/inputs/*nc")
#     combined_ids_paths = _infer_concat_order_from_positions(paths)
#     ids, paths = (
#         list(combined_ids_paths.keys()),
#         list(combined_ids_paths.values()),
#     )
#     open_kwargs = dict(engine="netcdf4", chunks={"Time": 24})


#     open_dataset= xr.open_dataset
#     # wrap the open_dataset, getattr, and preprocess with delayed
#     open_ = dask.delayed(open_dataset)
#     getattr_ = dask.delayed(getattr)

#     datasets = [open_(p, **open_kwargs) for p in paths]
#     closers = [getattr_(ds, "_close") for ds in datasets]

#     datasets, closers = dask.compute(datasets, closers)
#     concat_dim= ""
#     concat_dim=["Time"]
#     compat = "no_conflicts"
#     data_vars = "all"
#     coords="minimal"
#     combined = xr.core.combine._nested_combine(
#                     datasets,
#                     concat_dims=concat_dim,
#                     compat=compat,
#                     data_vars=data_vars,
#                     coords=coords,
#                     ids=ids,
#                     # join=join,
#                     # combine_attrs=combine_attrs,
#                 )

# def open_with_mfdataset():
#     ds = xr.open_mfdataset(
#         "examples/net_cdf/single_file_range/inputs/*",
#         # "examples/net_cdf/multi_file_range/inputs/demo_wrf_2017-10_*",
#         engine="netcdf4",
#         concat_dim="Time",
#         #   concat_dim=["X", "Y", "Time"],
#         combine="nested",
#         parallel=True,
#         # chunks={'x': 1, 'y': 1, "Time": 24},
#         coords="minimal",
#     )

# from timeit import repeat
# t1 = repeat(open_manually, number=1, repeat=10)
# t2 = repeat(open_with_mfdataset, number=1, repeat=10)
# plt.plot(t1,label="1")
# plt.plot(t2,label="2")
# print(t1, t2)
# plt.legend()
# # %%

# dsa = xr.open_mfdataset(
#   "examples/net_cdf/single_file_hour/inputs/*",
#   # "examples/net_cdf/multi_file_range/inputs/demo_wrf_2017-10_*",
#   engine="netcdf4",
#   concat_dim="Time",
# #   concat_dim=["X", "Y", "Time"],
#   combine="nested",
#   parallel=True,
# )

# dsb = xr.open_mfdataset(
#   "examples/net_cdf/single_file_range/inputs/*",
#   # "examples/net_cdf/multi_file_range/inputs/demo_wrf_2017-10_*",
#   engine="netcdf4",
#   concat_dim="Time",
# #   concat_dim=["X", "Y", "Time"],
#   combine="nested",
#   parallel=True,
# )
# # %%
# dsa.td_2m[:,0,0].plot()
# dsb.td_2m[:,0,0].plot()
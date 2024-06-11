# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr


# %%

DX = 3
DY = 2
DT = 24 * 4
data_out_loc = 'examples/net_cdf/multi_file_range/inputs'

# %%
demo_data_SWDOWN = np.ones((DT, DX, DY)) * 100
demo_data_HFX_FORCE = np.zeros((DT, DX, DY))
demo_data_td_2m = np.ones((DT, DX, DY)) * 10
demo_data_rh = np.zeros((DT, DX, DY)) + 0.3
demo_data_o3 = np.ones((DT, DX, DY)) *10
demo_data_wspeed = np.ones((DT, DX, DY)) * 1.4
demo_data_pres = np.ones((DT, DX, DY)) * 101
demo_data_RAINNC = np.ones((DT, DX, DY)) * 4
demo_data_SNOWH = np.zeros((DT, DX, DY))


# %%
XTIME = pd.date_range("2017-12-27", periods=DT, freq="1H")
reference_time = pd.Timestamp("2017-01-01")
lon = np.full((DX,DY), np.arange(DY))
lat = np.full((DY,DX), np.arange(DX)).transpose()

dims = ['Time', 'x', 'y']
coords = dict(
    lon=(['x', 'y'], lon),
    lat=(['x', 'y'], lat),
    XTIME=(['Time'], XTIME),
    reference_time=reference_time,
)
SWDOWN_da = xr.DataArray(demo_data_SWDOWN, name="SWDOWN", dims=dims, coords=coords)
HFX_FORCE_da = xr.DataArray(demo_data_HFX_FORCE, name="HFX_FORCE", dims=dims, coords=coords)
td_2m_da = xr.DataArray(demo_data_td_2m, name="td_2m", dims=dims, coords=coords)
rh_da = xr.DataArray(demo_data_rh, name="rh", dims=dims, coords=coords)
o3_da = xr.DataArray(demo_data_o3, name="o3", dims=dims, coords=coords)
wspeed_da = xr.DataArray(demo_data_wspeed, name="wspeed", dims=dims, coords=coords)
pres_da = xr.DataArray(demo_data_pres, name="pres", dims=dims, coords=coords)
RAINNC_da = xr.DataArray(demo_data_RAINNC, name="RAINNC", dims=dims, coords=coords)
SNOWH_da = xr.DataArray(demo_data_SNOWH, name="SNOWH", dims=dims, coords=coords)
# %%
SWDOWN_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_SWDOWN")
HFX_FORCE_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_HFX_FORCE")
td_2m_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_td_2m")
rh_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_rh")
o3_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_o3")
wspeed_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_wspeed")
pres_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_pres")
RAINNC_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_RAINNC")
SNOWH_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-10_SNOWH")

# %%

XTIMEB = XTIME + DT * XTIME.freq
reference_time = pd.Timestamp("2017-01-01")
lon = np.full((DX,DY), np.arange(DY))
lat = np.full((DY,DX), np.arange(DX)).transpose()

dims = ['Time', 'x', 'y']
coords = dict(
    lon=(['x', 'y'], lon),
    lat=(['x', 'y'], lat),
    XTIME=(['Time'], XTIMEB),
    reference_time=reference_time,
)

SWDOWN_da = xr.DataArray(demo_data_SWDOWN, name="SWDOWN", dims=dims, coords=coords)
HFX_FORCE_da = xr.DataArray(demo_data_HFX_FORCE, name="HFX_FORCE", dims=dims, coords=coords)
td_2m_da = xr.DataArray(demo_data_td_2m, name="td_2m", dims=dims, coords=coords)
rh_da = xr.DataArray(demo_data_rh, name="rh", dims=dims, coords=coords)
o3_da = xr.DataArray(demo_data_o3, name="o3", dims=dims, coords=coords)
wspeed_da = xr.DataArray(demo_data_wspeed, name="wspeed", dims=dims, coords=coords)
pres_da = xr.DataArray(demo_data_pres, name="pres", dims=dims, coords=coords)
RAINNC_da = xr.DataArray(demo_data_RAINNC, name="RAINNC", dims=dims, coords=coords)
SNOWH_da = xr.DataArray(demo_data_SNOWH, name="SNOWH", dims=dims, coords=coords)
# %%
SWDOWN_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_SWDOWN")
HFX_FORCE_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_HFX_FORCE")
td_2m_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_td_2m")
rh_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_rh")
o3_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_o3")
wspeed_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_wspeed")
pres_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_pres")
RAINNC_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_RAINNC")
SNOWH_da.to_netcdf(f"{data_out_loc}/demo_wrf_2017-11_SNOWH")

# # %%
# SNOWH_da
# # %%
# ds_loaded = xr.open_mfdataset(f"{data_out_loc}/demo_wrf_2017-10*")
# ds_loaded

# # %%
# ds_loaded['XTIME'].dims[0] == ds_loaded['td_2m'].dims[0]

# # %%

# data_processed = ds_loaded
# # TODO: Failing on multi run
# time_key = "XTIME"
# zero_year = 2017
# year = int(data_processed[time_key].values[0].astype(str)[0:4])
# day_offset = (year - zero_year) * 365 if zero_year is not None else 0

# data_processed = data_processed.assign(hr=lambda d: getattr(d, time_key).dt.hour)
# data_processed = data_processed.assign(
#     dd=lambda d: getattr(d, time_key).dt.strftime('%j').astype(int) + day_offset)
# T = data_processed[time_key].shape[0]
# data_processed.dd[0].values


# %%
t = XTIME + DT * XTIME.freq
t[0]

# %%
XTIME[0].day_of_year,XTIMEB[-1].day_of_year

# %%
# %%


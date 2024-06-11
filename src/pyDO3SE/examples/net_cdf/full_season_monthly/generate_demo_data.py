# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr


# %%
# Demo data
demo_data = pd.read_csv("examples/bangor_2015/inputs/bangor_2015_hb_ww.csv")
demo_data = pd.read_csv("examples/xiaoji_2008/outputs/A_03/winter_wheat/latest/external_data.csv")
demo_data
# %%
demo_data['R'].plot()

# %%

DX = 3
DY = 2
DT = len(demo_data)
data_out_loc = 'examples/net_cdf/full_season_monthly/inputs'

# %%
demo_data_SWDOWN = np.ones((DT, DX, DY))
demo_data_HFX_FORCE = np.zeros((DT, DX, DY))
demo_data_td_2m = np.ones((DT, DX, DY))
demo_data_rh = np.zeros((DT, DX, DY))
demo_data_o3 = np.ones((DT, DX, DY))
demo_data_wspeed = np.ones((DT, DX, DY))
demo_data_pres = np.ones((DT, DX, DY))
demo_data_RAINNC = np.ones((DT, DX, DY))
demo_data_SNOWH = np.zeros((DT, DX, DY))

# %%
for x in range(DX):
    for y in range(DY):
        demo_data_SWDOWN[:,x,y] = demo_data['R']
        demo_data_HFX_FORCE[:,x,y] = demo_data['Hd']
        demo_data_td_2m[:,x,y] = demo_data['Ts_C']
        demo_data_rh[:,x,y] = demo_data['RH']
        demo_data_o3[:,x,y] = demo_data['O3']
        demo_data_wspeed[:,x,y] = demo_data['u']
        demo_data_pres[:,x,y] = demo_data['P']
        demo_data_RAINNC[:,x,y] = demo_data['precip']
        # demo_data_SNOWH = demo_data['']


# %%
demo_data.dd.iloc[-1] - 365
# %%
time_data  = pd.date_range("2017-01-01", periods=365*24, freq="1H")
time_data[24*158]

# %%
# TODO: Set as range(DT)
time_data  = pd.date_range("2017-11-16", periods=DT, freq="1H")
time_data[0],time_data[-1]
# XTIME = time_data
# f"{data_out_loc}/demo_wrf_{XTIME[0].year}-{str(XTIME[0].month).zfill(2)}-{str(XTIME[0].day).zfill(2)}-{str(XTIME[0].time()).replace(':', '-')}"
# # %%
# t = 0
# XTIME = time_data[t:t + 24]
# XTIME
# # %%
# reference_time = pd.Timestamp("2017-01-01")
# lon = np.full((DX,DY), np.arange(DY))
# lat = np.full((DY,DX), np.arange(DX)).transpose()

# dims = ['Time', 'x', 'y']
# coords = dict(
#     lon=(['x', 'y'], lon),
#     lat=(['x', 'y'], lat),
#     XTIME=(['Time'], XTIME),
#     reference_time=reference_time,
# )
# SWDOWN_da = xr.DataArray(demo_data_SWDOWN[t:t+24], name="SWDOWN", dims=dims, coords=coords)
# HFX_FORCE_da = xr.DataArray(demo_data_HFX_FORCE[t:t+24], name="HFX_FORCE", dims=dims, coords=coords)
# td_2m_da = xr.DataArray(demo_data_td_2m[t:t+24], name="td_2m", dims=dims, coords=coords)
# rh_da = xr.DataArray(demo_data_rh[t:t+24], name="rh", dims=dims, coords=coords)
# o3_da = xr.DataArray(demo_data_o3[t:t+24], name="o3", dims=dims, coords=coords)
# wspeed_da = xr.DataArray(demo_data_wspeed[t:t+24], name="wspeed", dims=dims, coords=coords)
# pres_da = xr.DataArray(demo_data_pres[t:t+24], name="pres", dims=dims, coords=coords)
# RAINNC_da = xr.DataArray(demo_data_RAINNC[t:t+24], name="RAINNC", dims=dims, coords=coords)
# SNOWH_da = xr.DataArray(demo_data_SNOWH[t:t+24], name="SNOWH", dims=dims, coords=coords)

# ds = xr.merge(dict(
#     SWDOWN=SWDOWN_da,
#     HFX_FORCE=HFX_FORCE_da,
#     td_2m=td_2m_da,
#     rh=rh_da,
#     o3=o3_da,
#     wspeed=wspeed_da,
#     pres=pres_da,
#     RAINNC=RAINNC_da,
#     SNOWH=SNOWH_da,
# ).values())
# ds
# %%
demo_data_HFX_FORCE.shape
# %%
# iterate over months in time_data
print(time_data[0], time_data[-1])
print(time_data[0].month, time_data[-1].month)
for mmT in range(time_data[0].month, time_data[-1].month + 12):
    mm = mmT % 12
    print(mm)
    start_index = np.where(time_data.month == mm+1)[0][0]
    end_index = np.where(time_data.month == mm+1)[0][-1]
    print(start_index, end_index)
    XTIME = time_data[start_index:end_index]
    print(XTIME[0], XTIME[-1])
    print(f"{data_out_loc}/demo_wrf_{XTIME[0].year}-{str(XTIME[0].month).zfill(2)}.nc")
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
    SWDOWN_da = xr.DataArray(demo_data_SWDOWN[start_index: end_index], name="SWDOWN", dims=dims, coords=coords)
    HFX_FORCE_da = xr.DataArray(demo_data_HFX_FORCE[start_index: end_index], name="HFX_FORCE", dims=dims, coords=coords)
    td_2m_da = xr.DataArray(demo_data_td_2m[start_index: end_index], name="td_2m", dims=dims, coords=coords)
    rh_da = xr.DataArray(demo_data_rh[start_index: end_index], name="rh", dims=dims, coords=coords)
    o3_da = xr.DataArray(demo_data_o3[start_index: end_index], name="o3", dims=dims, coords=coords)
    wspeed_da = xr.DataArray(demo_data_wspeed[start_index: end_index], name="wspeed", dims=dims, coords=coords)
    pres_da = xr.DataArray(demo_data_pres[start_index: end_index], name="pres", dims=dims, coords=coords)
    RAINNC_da = xr.DataArray(demo_data_RAINNC[start_index: end_index], name="RAINNC", dims=dims, coords=coords)
    SNOWH_da = xr.DataArray(demo_data_SNOWH[start_index: end_index], name="SNOWH", dims=dims, coords=coords)

    ds = xr.merge(dict(
        SWDOWN=SWDOWN_da,
        HFX_FORCE=HFX_FORCE_da,
        td_2m=td_2m_da,
        rh=rh_da,
        o3=o3_da,
        wspeed=wspeed_da,
        pres=pres_da,
        RAINNC=RAINNC_da,
        SNOWH=SNOWH_da,
    ).values())
    ds.to_netcdf(f"{data_out_loc}/demo_wrf_{XTIME[0].year}-{str(XTIME[0].month).zfill(2)}.nc")

# %%
# %%
# %%
# Check remerge

print(f"{DX}, {DY}, {DT}")
ds_loaded = xr.open_mfdataset(f"{data_out_loc}/*.nc", concat_dim="Time", combine="nested")
ds_loaded

# %%
ds_loaded.td_2m[:,0,0].plot()

# %%

# ==== Create e_State_overrides.nc ====
# BELOW IS A WIP!
# %%
# ext_overrides = xr.open_dataset('examples/net_cdf/full_season_monthly/e_state_overrides.nc')
# ext_overrides
# # %%
# ext_overrides['lat'].values
# # %%
# ext_overrides['x'].shape

# # %%
# x = np.arange(0, 5, 1)
# y = np.arange(0, 4, 1)
# # %%
# ext_overrides['lat'].plot()
# # %%

# lat = np.array([[i for i in range(len(x))] for _ in range(len(y))])
# lat.shape, ext_overrides['lat'].shape, lat
# # %%
# lon = np.array([[i for _ in range(len(x))] for i in range(len(y))])

# # %%
# ext_overrides['lat'] = ext_overrides['lat'] + lat
# # %%
# ext_overrides['lat'].plot()
# # %%
# ext_overrides.to_netcdf('examples/net_cdf/full_season_monthly/e_state_overrides_new.nc')
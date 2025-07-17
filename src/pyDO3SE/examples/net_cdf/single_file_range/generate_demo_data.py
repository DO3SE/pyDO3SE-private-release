# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr


# %%

DX = 3
DY = 2
DT = 24 * 8
data_out_loc = 'examples/net_cdf/single_file_range/inputs'

# %%
# daily_sqdown = [0,1,2,3,4,5,6,7,8,9,9.5, 9.8,10,9.8,9.5,9,8,7,6,5,4,3,2,1]
# NOTE: Par is offset by 3 hours to demostrate timezone issues
daily_sqdown = [3,4,5,6,7,8,9,9.5, 9.8,10,9.8,9.5,9,8,7,6,5,4,3,2,1,0,1,2]
demo_data_SWDOWN = np.array([[[daily_sqdown[i] for _ in range(DT) for i in range(24)] for _ in range(DX)] for _ in range(DY)]).transpose()
demo_data_HFX_FORCE = np.zeros((DT, DX, DY))
demo_data_td_2m = np.ones((DT, DX, DY)) * 10 + np.random.rand(DT, DX, DY)
demo_data_rh = np.zeros((DT, DX, DY)) + 0.3
demo_data_o3 = np.ones((DT, DX, DY)) *10
demo_data_wspeed = np.ones((DT, DX, DY)) * 1.4
demo_data_pres = np.ones((DT, DX, DY)) * 101
demo_data_RAINNC = np.ones((DT, DX, DY)) * 4
demo_data_SNOWH = np.zeros((DT, DX, DY))


# %%
# TODO: Set as range(DT)
time_data  = pd.date_range("2017-12-27", periods=DT, freq="1H")
XTIME = time_data
f"{data_out_loc}/demo_wrf_{XTIME[0].year}-{str(XTIME[0].month).zfill(2)}-{str(XTIME[0].day).zfill(2)}-{str(XTIME[0].time()).replace(':', '-')}"
# %%
t = 0
XTIME = time_data[t:t + 24]
XTIME
# %%
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
SWDOWN_da = xr.DataArray(demo_data_SWDOWN[t:t+24], name="SWDOWN", dims=dims, coords=coords)
HFX_FORCE_da = xr.DataArray(demo_data_HFX_FORCE[t:t+24], name="HFX_FORCE", dims=dims, coords=coords)
td_2m_da = xr.DataArray(demo_data_td_2m[t:t+24], name="td_2m", dims=dims, coords=coords)
rh_da = xr.DataArray(demo_data_rh[t:t+24], name="rh", dims=dims, coords=coords)
o3_da = xr.DataArray(demo_data_o3[t:t+24], name="o3", dims=dims, coords=coords)
wspeed_da = xr.DataArray(demo_data_wspeed[t:t+24], name="wspeed", dims=dims, coords=coords)
pres_da = xr.DataArray(demo_data_pres[t:t+24], name="pres", dims=dims, coords=coords)
RAINNC_da = xr.DataArray(demo_data_RAINNC[t:t+24], name="RAINNC", dims=dims, coords=coords)
SNOWH_da = xr.DataArray(demo_data_SNOWH[t:t+24], name="SNOWH", dims=dims, coords=coords)

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
ds
# %%
for t in range(0, DT, 24):
    XTIME = time_data[t:t + 24]
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
    SWDOWN_da = xr.DataArray(demo_data_SWDOWN[t:t+24], name="SWDOWN", dims=dims, coords=coords)
    HFX_FORCE_da = xr.DataArray(demo_data_HFX_FORCE[t:t+24], name="HFX_FORCE", dims=dims, coords=coords)
    td_2m_da = xr.DataArray(demo_data_td_2m[t:t+24], name="td_2m", dims=dims, coords=coords)
    rh_da = xr.DataArray(demo_data_rh[t:t+24], name="rh", dims=dims, coords=coords)
    o3_da = xr.DataArray(demo_data_o3[t:t+24], name="o3", dims=dims, coords=coords)
    wspeed_da = xr.DataArray(demo_data_wspeed[t:t+24], name="wspeed", dims=dims, coords=coords)
    pres_da = xr.DataArray(demo_data_pres[t:t+24], name="pres", dims=dims, coords=coords)
    RAINNC_da = xr.DataArray(demo_data_RAINNC[t:t+24], name="RAINNC", dims=dims, coords=coords)
    SNOWH_da = xr.DataArray(demo_data_SNOWH[t:t+24], name="SNOWH", dims=dims, coords=coords)

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
    ds.to_netcdf(f"{data_out_loc}/demo_wrf_{XTIME[0].year}-{str(XTIME[0].month).zfill(2)}-{str(XTIME[0].day).zfill(2)}.nc")
print(f"output saved to {data_out_loc}")
ds
# %%
time_data[-1].day_of_year,time_data[0].day_of_year
# %%
# %%
# Check remerge

print(f"{DX}, {DY}, {DT}")
ds_loaded = xr.open_mfdataset(f"{data_out_loc}/*.nc", concat_dim="Time", combine="nested")
ds_loaded

# %%
ds_loaded.td_2m[:,0,0].plot()
# %%
ds_loaded.SWDOWN[:,0,0].plot()

# %%

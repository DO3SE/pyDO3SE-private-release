# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr


# %%

DX = 3
DY = 2
DT = 24 * 8
data_out_loc = 'examples/net_cdf/single_file_hour/inputs'

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
# TODO: Set as range(DT)
time_data  = pd.date_range("2017-12-27", periods=DT, freq="1H")
XTIME = time_data
f"{data_out_loc}/demo_wrf_{XTIME[0].year}-{str(XTIME[0].month).zfill(2)}-{str(XTIME[0].day).zfill(2)}-{str(XTIME[0].time()).replace(':', '-')}"

# %%
for t in range(DT):
    XTIME = [time_data[t]]
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
    SWDOWN_da = xr.DataArray([demo_data_SWDOWN[t]], name="SWDOWN", dims=dims, coords=coords)
    HFX_FORCE_da = xr.DataArray([demo_data_HFX_FORCE[t]], name="HFX_FORCE", dims=dims, coords=coords)
    td_2m_da = xr.DataArray([demo_data_td_2m[t]], name="td_2m", dims=dims, coords=coords)
    rh_da = xr.DataArray([demo_data_rh[t]], name="rh", dims=dims, coords=coords)
    o3_da = xr.DataArray([demo_data_o3[t]], name="o3", dims=dims, coords=coords)
    wspeed_da = xr.DataArray([demo_data_wspeed[t]], name="wspeed", dims=dims, coords=coords)
    pres_da = xr.DataArray([demo_data_pres[t]], name="pres", dims=dims, coords=coords)
    RAINNC_da = xr.DataArray([demo_data_RAINNC[t]], name="RAINNC", dims=dims, coords=coords)
    SNOWH_da = xr.DataArray([demo_data_SNOWH[t]], name="SNOWH", dims=dims, coords=coords)

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
    ds.to_netcdf(f"{data_out_loc}/demo_wrf_{XTIME[0].year}-{str(XTIME[0].month).zfill(2)}-{str(XTIME[0].day).zfill(2)}-{str(XTIME[0].time()).replace(':', '-')}.nc")

# %%
time_data[-1].day_of_year,time_data[0].day_of_year
# %%

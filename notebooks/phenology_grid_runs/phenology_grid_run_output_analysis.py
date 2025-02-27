# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
from typing import List
import xarray as xr
import os

# %%
from pyDO3SE.Config.config_loader import config_loader_pickled, Config_Shape, config_loader
from pyDO3SE.Model_State.model_state_loader import model_state_loader_quick, Model_State_Shape
from pyDO3SE.Output.process_outputs import dump_config_to_file_json

# %%
ds = xr.open_mfdataset('examples/net_cdf/full_season_monthly/runs/phenology_grid_run/bangor_wheat/outputs_grid/*.nc')
ds

# %%
ds.dvi[:,:,-1].plot()
# ds.dvi[0,0,:].plot()
# %%

processed_configs_dir = "examples/net_cdf/full_season_monthly/runs/phenology_grid_run/bangor_wheat/processed_configs"
files = os.listdir(processed_configs_dir)

processed_configs: List[Config_Shape] = []
for f in files:
    processed_config: Config_Shape = config_loader_pickled(os.path.join(
        processed_configs_dir, f.replace('state', 'config')))
    processed_configs.append(processed_config)

# %%
sowing_dates = []
latitudes_phenology = []
latitudes = []
for i in range(len(processed_configs)):
    sowing_dates.append(processed_configs[i].Land_Cover.parameters[0].phenology.key_dates.sowing)
    latitudes_phenology.append(processed_configs[i].Land_Cover.phenology_options.latitude)
    latitudes.append(processed_configs[i].Location.lat)

sowing_dates, latitudes, latitudes_phenology
# %%

# %%
# %%
%load_ext autoreload
%autoreload 2
from do3se_met.wind import *

# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# %%
h = 20
w = 0.1 # LAI weighted Leaf width
SAI = 1
u_at_canopy_top = 1
canopy_layer_height = 1
z = 1
o_top_chamber = False
layer_depth = 1
wind_speed = calc_layer_windspeed(
    h=h,
    w=w,
    SAI=SAI,
    u_at_canopy_top=u_at_canopy_top,
    z=z,
    o_top_chamber=o_top_chamber,
    layer_depth=layer_depth,
)
wind_speed
# %%

for SAI in np.arange(1,10, 2):
    u_at_canopy_top = 20
    # SAI = 1
    wind_speeds = []
    layer_heights = np.arange(5, h, 2)
    layer_depth = 0
    for z in layer_heights:
        wind_speed = calc_layer_windspeed(
            h=h,
            w=w,
            SAI=SAI,
            u_at_canopy_top=u_at_canopy_top,
            z=z,
            o_top_chamber=o_top_chamber,
            layer_depth=layer_depth,
        )
        wind_speeds.append(wind_speed)
    wind_speeds.append(u_at_canopy_top)
    x = np.append(layer_heights, h)
    plt.scatter(x, wind_speeds, label=SAI)
plt.xlabel("Height (m)")
plt.ylabel("Wind speed (m/s)")
plt.legend(title="SAI")
# %%
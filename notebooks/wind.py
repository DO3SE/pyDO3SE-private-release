# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# %%
from do3se_met.wind import *

# %%
u_i_list = []
x = []
c = []
for ustar in range(0,3):
    for h in range(1,20):
        try:
            h_u = h
            izr = 20
            z_u=h
            canopy_d = 1 # 0.78
            canopy_z0 = 1# 0.07
            # Forest specific
            u_z0=min(1.0, h_u * canopy_z0)
            d=h * canopy_d
            u_d = h_u * canopy_d
            z0=min(1.0, h * canopy_z0)
            u = 1.234
            assert z0 == 1.0
            assert u_z0 == 1.0
            ustar_ref=ustar
            min_windspeed = 0.01

            z_u_a = z_u
            u_z0_a = u_z0
            u_d_a = u_d

            # Find windspeed at izR, over reference canopy
            u_i = velocity_from_ustar(ustar_ref, (izr - u_d_a), u_z0_a)
            u_i_lim =  max(min_windspeed,u_i)
            u_i_list.append(u_i_lim)
            x.append(h)
            c.append(ustar)
        except Exception as e:
            print(e)
            print(f"Failed with {h}")
            raise e
c = np.array(c)

cs = plt.scatter(x, u_i_list, c=c/3)
plt.colorbar()
# %%
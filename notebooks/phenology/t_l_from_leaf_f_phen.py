"""This is a python notebook to get the key leaf phenology stages(t_l...) from leaf_f_fphen data.

"""

# %%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from pyDO3SE.plugins.met.thermal_time import *

data = pd.read_csv('../icp_runs/08_multiplicative/inputs/SE87_NF_thermal.csv')

data.head()

# %%

data.columns

# %%
t_b = 0
t_o = 20
t_m = 30

# %%

data['td'] = calc_thermal_time(data['Ts_C'], 0, int(len(data) / 24), t_b)
data.head()
data['td'].plot()

# %%
leaf_f_phen = data['Leaf_f_phen'].groupby(data.index // 24).mean()
# np.repeat(leaf_f_phen.diff().values, 24)
data['leaf_f_phen_diff'] = np.repeat(leaf_f_phen.diff().values, 24)
data['leaf_f_phen_diff'] = data['leaf_f_phen_diff'].fillna(0)
data.head()

# %%
data['leaf_f_phen_diff'].max(), data['leaf_f_phen_diff'].min()

# %%
t_lse_start = None
t_l_end = None
for i, (
    leaf_f_phen,
    leaf_f_phen_diff,
    td,
) in enumerate(zip(
        data['Leaf_f_phen'],
        data['Leaf_f_phen'].diff(),
        data['td'],
    )
):
    if leaf_f_phen_diff < -0.03 and i > 5 and t_lse_start is None:
        t_lse_start = td
        print(i, leaf_f_phen_diff, td)

    if t_lse_start and leaf_f_phen_diff == 0:
        t_l_end = td

t_lse = t_l_end - t_lse_start
t_lse

# %%
t_l = 4.55844 * t_lse
sowing_date_td = t_l_end - t_l

sowing_date = None

for i, (
    dd,
    td,
) in enumerate(zip(
        data['dd'],
        data['td'],
    )
):
    if td >= sowing_date_td:
        sowing_date = dd
        break

sowing_date
# %%wa

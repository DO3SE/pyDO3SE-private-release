# %%
from pprint import pprint
from pyDO3SE.plugins.resistance.helpers import calc_Rsur_multilayer
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pyDO3SE.plugins.O3.helpers import calc_fst_leaf
# %%


data_1 = pd.read_csv(
    'examples/bangor_2015_multilayer/outputs/bangor_2015_hb_ww/bangor_wheat/4.16.3b/bangor_wheat-bangor_2015_hb_ww_out.csv')
data_2 = pd.read_csv(
    'examples/bangor_2015_multilayer/outputs/bangor_2015_hb_ww/bangor_wheat_single_layer/latest/bangor_wheat_single_layer-bangor_2015_hb_ww_out.csv')
data_3 = pd.read_csv(
    'examples/bangor_2015_multilayer/outputs/bangor_2015_hb_ww/bangor_wheat_no_ozone/latest/bangor_wheat_no_ozone-bangor_2015_hb_ww_out.csv')


data_b_1 = pd.read_csv(
    'examples/bangor_2015_multilayer/outputs/bangor_2015_hp_ww_skyfall/bangor_wheat/4.16.3b/bangor_wheat-bangor_2015_hp_ww_skyfall_out.csv')
data_b_2 = pd.read_csv(
    'examples/bangor_2015_multilayer/outputs/bangor_2015_hp_ww_skyfall/bangor_wheat_single_layer/latest/bangor_wheat_single_layer-bangor_2015_hp_ww_skyfall_out.csv')
data_b_3 = pd.read_csv(
    'examples/bangor_2015_multilayer/outputs/bangor_2015_hb_ww/bangor_wheat_no_ozone/latest/bangor_wheat_no_ozone-bangor_2015_hb_ww_out.csv')

# print(list(data_1.columns))

# %%

data_1['O3up_2'].iloc[2000:2800].plot()
data_1['O3up_1'].iloc[2000:2800].plot()
data_1['O3up_0'].iloc[2000:2800].plot()
# %%
data_2['O3up_2'].iloc[2000:2800].plot()
data_2['O3up_1'].iloc[2000:2800].plot()
data_2['O3up_0'].iloc[2000:2800].plot()

# %%
# for iP in range(3):
for iL in range(4):
    data_1[f'leaf_pop_distribution_{iL}_{2}'].iloc[2000:2800].plot()
for iL in range(1):
    data_2[f'leaf_pop_distribution_{iL}_{2}'].iloc[2000:2800].plot()
plt.legend()

# %%
# for iL in range(4):
    # data_1[f'leaf_pop_distribution_{iL}_{1}'].iloc[2000:2800].plot()
for iL in range(1):
    data_2[f'leaf_pop_distribution_{iL}_{1}'].iloc[2000:2800].plot()
plt.legend()
# %%
data_2['O3up_2'].iloc[1000:2800].plot()
data_1['O3up_2'].iloc[1000:2800].plot()

# %%
data_2['O3up_1'].iloc[1000:2800].plot()
data_2['O3up_2'].iloc[1000:2800].plot()
# %%
data_2['micro_O3_0'].iloc[1000:2800].plot()
data_1['micro_O3_0'].iloc[1000:2800].plot()
data_1['micro_O3_1'].iloc[1000:2800].plot()
data_1['micro_O3_2'].iloc[1000:2800].plot()

# %%
# These are very different!
data_2['o3_ppb_i'].iloc[2000:2300].plot()
data_1['o3_ppb_i'].iloc[2000:2300].plot()

# %%
data_2['canopy_height'].iloc[2000:2300].plot()
data_1['canopy_height'].iloc[2000:2300].plot()
# %%
data_2['rsto_c'].iloc[2000:2300].plot()

# %%
data_1['fO3_d'].iloc[1000:].plot(label="multi")
data_2['fO3_d'].iloc[1000:].plot(label="single")
data_3['fO3_d'].iloc[1000:].plot(label="noozone")
plt.legend()
# %%
data_1['f_LS'].iloc[1000:].plot(label="multi")
data_2['f_LS'].iloc[1000:].plot(label="single")
data_3['f_LS'].iloc[1000:].plot(label="noozone")
plt.legend()

# %%
data_1['fO3_l'].iloc[1000:].plot(label="multi")
data_2['fO3_l'].iloc[1000:].plot(label="single")
data_3['fO3_l'].iloc[1000:].plot(label="noozone")
plt.legend()
# %%

# ================ COMPARE Low and high ozone =========


data_1['O3up_2'].iloc[2000:2200].plot(label="low")
data_b_1['O3up_2'].iloc[2000:2200].plot(label="high")
plt.legend()
# %%
data_1['fO3_d'].iloc[2000:2200].plot(label="low")
data_b_1['fO3_d'].iloc[2000:2200].plot(label="high")
plt.legend()

# %%
data_1['f_LS'].iloc[2500:2800].plot(label="low")
data_b_1['f_LS'].iloc[2500:2800].plot(label="high")
plt.legend()

# %%
data_1['f_LA'].iloc[1000:3800].plot(label="low")
data_b_1['f_LA'].iloc[1000:3800].plot(label="high")
plt.legend()

# %%
data_1['fst_acc'].iloc[2000:3000].plot(label="low")
data_b_1['fst_acc'].iloc[2000:3000].plot(label="high")
plt.legend()


# %%
data_1['leaf_f_phen'].iloc[2000:3000].plot(label="low")
data_b_1['leaf_f_phen'].iloc[2000:3000].plot(label="high")
plt.legend()


# %%
data_1['fst_canopy'].iloc[1000:3800].plot(label="low")
data_b_1['fst_canopy'].iloc[1000:3800].plot(label="high")
plt.legend()
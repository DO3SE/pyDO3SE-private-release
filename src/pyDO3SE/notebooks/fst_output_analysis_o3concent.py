# %%
!python pyDO3SE_cli.py multi-run --project_directory=examples/bangor_2015_multilayer_simple --verbose=True --base-config-file=examples/bangor_2015_multilayer_simple/base_config.json  # type: ignore noqa E402
# %%
from pprint import pprint
from pyDO3SE.plugins.resistance.helpers import calc_Rsur_multilayer
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pyDO3SE.plugins.O3.helpers import calc_fst_leaf
# %%
# %%
# ========================COMPARE HIGH AND LOW OZONE =================


# %%
data_1 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/outputs/bangor_wheat_1_layer/bangor_2015_hp_ww_skyfall.csv')

data_1_b = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/outputs/bangor_wheat_1_layer/bangor_2015_hb_ww.csv')
# %%
(data_1['leaf_f_phen']*100).plot()
data_1['fst_acc'].plot()
# %%

data_1['fst_acc'].plot(label="high")
data_1_b['fst_acc'].plot(label="low")
plt.legend()
# %%

data_1['f_LS'].plot(label="high")
data_1_b['f_LS'].plot(label="low")
plt.legend()

# %%

data_1['f_LA'].plot(label="high")
data_1_b['f_LA'].plot(label="low")
plt.legend()
# %%

data_1['fO3_l'].plot(label="high")
data_1_b['fO3_l'].plot(label="low")
plt.legend()

# %%
data_1['fst'].plot(label="high")
data_1_b['fst'].plot(label="low")
plt.legend()


# %%
data_1['leaf_f_phen'].plot()
data_1['f_phen'].plot()

# %%
data_1['leaf_pop_distribution_0_0'].plot()
data_1['leaf_pop_distribution_0_0'].plot()


# %%
# external data
# %%
ext_data_1 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/inputs/bangor_2015_hp_ww_skyfall.csv')

ext_data_1 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/inputs/bangor_2015_hb_ww.csv')
# %%
ext_data_1['O3, ppb'].plot()
ext_data_1['O3, ppb'].plot()
# %%
data_1['leaf_pop_distribution_0_0'].plot()
data_1['leaf_pop_distribution_0_1'].plot()
data_1['leaf_pop_distribution_0_2'].plot()


# %%
data_1['O3up_0'].plot()
data_1['O3up_1'].plot()
# data_1['fst_acc'].plot()
# %%
data_1['fst_canopy'].plot()

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


data_1 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/outputs/bangor_wheat_1_layer/bangor_2015_hp_ww_skyfall.csv')
data_2 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/outputs/bangor_wheat_2_layer/bangor_2015_hp_ww_skyfall.csv')
data_3 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/outputs/bangor_wheat_3_layer/bangor_2015_hp_ww_skyfall.csv')
data_4 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/outputs/bangor_wheat_4_layer/bangor_2015_hp_ww_skyfall.csv')
data_5 = pd.read_csv(
    'examples/bangor_2015_multilayer_simple/outputs/bangor_wheat_5_layer/bangor_2015_hp_ww_skyfall.csv')
data = [
    data_1, data_2, data_3, data_4, data_5
]
pprint(list(data_2.columns))

# %%
# Check fst output
plt.figure(figsize=(12, 4))
for i, d in enumerate(data):
    d[f'O3up_{2}'].iloc[200:360].plot(label=i)

plt.legend()
# %%
# Check fst output
plt.figure(figsize=(12, 4))
for i, d in enumerate(data):
    d[f'pody'].iloc[320:360].plot(label=i)

plt.legend()

# %%
"""Here we can see that single layer model has same fst for all populations."""
row_index = 320
for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'O3up_{iP}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=iP)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("fst")
plt.title(f"fst at row {row_index}")
plt.show()

# %%
row_index = 320
for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'rb_l_{i}_{iP}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=iP)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("rb")
plt.title(f"rb at row {row_index} top layer")
plt.show()
# %%
row_index = 320
for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'rb_l_{0}_{iP}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=iP)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("Rb")
plt.title(f"Rb at row {row_index} at bottom layer")
plt.show()

# %%
"""Here we can see that single layer model has same Rsto for all populations."""
row_index = 320
for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'rsto_l_{i}_{iP}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=iP)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("Rsto")
plt.title(f"Rsto at row {row_index} at top layer")
plt.show()
# %%
row_index = 320
for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'gsto_l_{iP}_{0}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=iP)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("Gsto")
plt.title(f"Gsto at row {row_index} at bottom layer")
plt.show()

# %%
row_index = 320
for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'gsto_l_{iP}_{i}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=iP)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("Gsto")
plt.title(f"Gsto at row {row_index} at top layer")
plt.show()

# %%
# row_index = 320
# for iP in range(3):
#     x = [i + 1 for i, _ in enumerate(data)]
#     y = [d[f'rsto_l_{i}_{iP}'].iloc[row_index] for i, d in enumerate(data)]
#     plt.scatter(x, y, label=iP)
#     plt.plot(x, y)

for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'rsto_l_{0}_{iP}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=f"{iP}_bottom")
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("Gsto")
plt.title(f"Gsto at row {row_index} at top layer")
plt.show()

# %%
row_index = 320
iP = 2  # Flag leaf
for i, d in enumerate(data):
    x = [j+0.1*i for j in range(i + 1)]
    y = [d[f'gsto_l_{iP}_{j}'].iloc[row_index] for j in range(i + 1)]
    plt.scatter(x, y, label=i)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("Gsto")
plt.title(f"Gsto at row {row_index} for number of layers")
plt.show()

# %%
row_index = 320
iP = 2  # Flag leaf
for i, d in enumerate(data):
    x = [j for j in range(i + 1)]
    y = [d[f'O3up_{iP}'].iloc[row_index] for j in range(i + 1)]
    plt.scatter(x, y, label=i)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("fst")
plt.title(f"fst at row {row_index} for number of layers")
plt.show()

# %%
row_index = 320
for iP in range(3):
    x = [i + 1 for i, _ in enumerate(data)]
    y = [d[f'leaf_pop_distribution_{i}_{iP}'].iloc[row_index] for i, d in enumerate(data)]
    plt.scatter(x, y, label=iP)
    plt.plot(x, y)
plt.legend()
plt.xlabel("model layer count")
plt.ylabel("leaf_pop_distribution_0_0")
plt.title(f"leaf_pop_distribution_0_0 at row {row_index}")
plt.show()
# %%
for i, d in enumerate(data):
    plt.plot([d[f'leaf_pop_distribution_{i}_{iP}'].iloc[row_index] for iP in range(3)] ,label=i)
plt.legend()

# %%
for i, d in enumerate(data):
    plt.plot([d[f'fLAI_{iP}_{i}'].iloc[row_index] for iP in range(3)] ,label=i)
plt.legend()

# %%
for i, d in enumerate(data):
    plt.plot([d[f'rb_l_{i}_{iP}'].iloc[row_index]+i*0.00005 for iP in range(3)], label=i)
plt.legend()
# %%
plt.plot([d[f'micro_O3_{i}'].iloc[row_index] for i, d in enumerate(data)])


# %%
plt.plot([d[f'canopy_height'].iloc[row_index] for i, d in enumerate(data)])
# %%
plt.plot([d[f'rsur'].iloc[row_index] for i, d in enumerate(data)])
# %%
plt.plot([d[f'rsto_c'].iloc[row_index] for i, d in enumerate(data)])
# %%
row_index=140
plt.plot([d[f'layer_lai_{i}'].iloc[row_index] for i, d in enumerate(data)])

# %%
for i, d in enumerate(data):
    plt.plot(d['rsur'].iloc[138:155], label=i)
plt.legend()
# %%
for i, d in enumerate(data):
    plt.plot(d[f'layer_lai_{i}'].iloc[138:155], label=i)
plt.legend()
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

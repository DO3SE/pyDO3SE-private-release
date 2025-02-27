"""Compare outputs in examples/outputs directory against DO3SE-UI."""
# %%
from matplotlib import pyplot as plt
import pandas as pd
from pprint import pprint
import os
from pyDO3SE.version import version
from pyDO3SE.util.Objects import Field
from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.Analysis.charts import annual_graph, calc_moving_average, multi_series_annual_graph


# %%
# CONFIG
output_data_dir = os.path.normpath('examples/spanish_wheat/output')
# Set output_dir to save output
output_dir = None  # os.path.normpath(f'./examples/spanish_wheat/comparisons/{version}')
fields_to_graph = [  # Must match fields in Output_Shape.py
    'gsto_l',
    'gsto',
    'rsto',
    'rsto_l',
    'dvi',
    'PARsun',
    'PARshade',
    'fst',
    'A_n',
    'lai',
]

# make sure outputdir exists
output_dir and os.makedirs(output_dir, exist_ok=True)

# %%
# Find all csv files in outputs datafiles directory
data_files = ((
    dr.split('\\')[-1],
    dr,
    next((fl for fl in f if fl[-4:] == '.csv'), None),
) for dr, d, f in os.walk(output_data_dir) if len(f) > 0)

# %%
# Setup data
series_titles = []
data_sets = {f: [] for f in fields_to_graph}
for name, dir, file in data_files:
    print(name)
    series_titles.append(name)

    data = pd.read_csv(os.path.join(dir, file))
    for f in fields_to_graph:
        # Note if data does not exist we set to 0
        if f not in data.columns:
            print(f'WARNING {f} not in {name} data')
        field_data = data[f].values if f in data.columns else [0 for i in range(len(data.index))]
        data_sets[f].append(field_data)

# %%
# Create a graph for each field
for f in fields_to_graph:
    multi_series_annual_graph(
        f,
        data_sets[f],
        series_titles,
        output_fields_map[f],
        output_dir=output_dir,
        label_x_days=30,
        average_step=24,
        chart_id=version,
        linestyles=[{
            'style': 'dashed' if s == 'DO3SE_ui' else 'solid',
            'width': 0.8,
            'zorder': 100 if s == 'DO3SE_ui' else i,
        } for i, s in enumerate(series_titles)],
        figsize=(8, 3),
        dpi=300,
    )
plt.show()

# %%

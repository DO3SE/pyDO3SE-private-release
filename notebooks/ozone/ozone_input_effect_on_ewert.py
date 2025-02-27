"""Experimenting with the model outputs when we vary the input ozone."""

# %%
import pandas as pd
from pyDO3SE.Analysis.charts import annual_graph, multi_series_annual_graph
from pyDO3SE.main import Main_Overrides, run_model, setup_config, setup_external_state
from pyDO3SE.Analysis.util import output_log_to_field_data
from pyDO3SE.Output.Output_Shape import output_fields_map

# %%
config_file = "./examples/spanish_wheat/spanish_wheat_config.json"
data_file = "./examples/spanish_wheat/spanish_wheat_data.csv"
output_directory = './notebooks/ozone/comparisons_output'
output_fields = ['gsto', 'gsto_l', 'A_n', 'fst', 'fO3_d', 'td_dd', 'rsto_l']

# variation_count = 6

fractions = [0, 0.1, 0.6, 1.0, 2]
variation_count = len(fractions)


def get_ozone_fac(i, count):
    return fractions[i]
    # return 10 ** (i - 1)
    # return 2 * i / (count - 1)
    # return i * 0.333
# %%

# TODO: create a run with ozone overrides


overrides = Main_Overrides(0, 365)

# data_all = {f: [] for f in output_fields}

for i in range(variation_count):
    # Modify O3nmol
    config = setup_config(config_file)
    external_state = setup_external_state(config, data_file)
    ozone_fraction = get_ozone_fac(i, variation_count)
    print(ozone_fraction)
    external_state.O3_nmol = [o * ozone_fraction for o in external_state.O3_nmol]
    # external_state.O3_nmol = [(o * float(i) / float(variation_count)) for o in external_state.O3_nmol]
    # # Run model
    final_state, output_logs = run_model(config, external_state, overrides)
    output_filename = f'run_{i}.csv'
    full_logs = pd.DataFrame(output_logs)
    full_logs.to_csv(f"{output_directory}/{output_filename}")
    for f in output_fields:
        field = output_fields_map[f]
        # data_all[f].append(output_logs[f])
        data = output_log_to_field_data(output_logs, field)
        annual_graph(
            data,
            field,
            label_x_days=30,
            average_step=24,
            output_dir=output_directory,
            chart_id=f'{field.id}-{i}',
        )

# %%
start = 147
end = 230
data_all = {f: [] for f in output_fields}
for i in range(variation_count):
    output_filename = f'run_{i}.csv'
    loaded_data = pd.read_csv(f"{output_directory}/{output_filename}")
    # print(loaded_data.head())
    for f in output_fields:
        data_all[f].append(loaded_data[f].values)


for f in output_fields:
    series_titles = ['{:.0f}%'.format(100 * get_ozone_fac(i, variation_count))
                     for i in range(variation_count)]

    multi_series_annual_graph(
        f,
        [x[start * 24:end * 24] for x in data_all[f]],
        series_titles,
        output_fields_map[f],
        output_dir=output_directory,
        label_x_days=30,
        average_step=24,
        chart_id=f'{f}-comparison',
        linestyles=[{
            'style': 'dashed' if s == 'DO3SE_ui' else 'solid',
            'width': 0.8,
            'zorder': 100 if s == 'DO3SE_ui' else i,
        } for i, s in enumerate(series_titles)],
        figsize=(8, 3),
        dpi=300,
        output_format='svg',
        start_day=start,
        end_day=end,
    )


# %%
# CHECKING EWERT
# from pyDO3SE.plugins.gsto.ewert.ewert import ewert


# ewert(

# )

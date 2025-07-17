# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %%'
# %%
"""Test the full model."""

from proflow.Objects.Process import Process
from proflow.Objects.Interface import I
from pyDO3SE.plugins.gsto.ewert.ewert import ewert
import numpy as np
from functools import reduce
from matplotlib import pyplot as plt
from dataclasses import asdict
import json
import pandas as pd
from timeit import timeit
from helpers.list_helpers import flatten_list

from helpers.named_tuple_helpers import unpack
from proflow.ProcessRunnerCls import ProcessRunner


from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.External_State.external_state_loader import load_external_state
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Defaults.state_init_processes import state_init_processes
from pyDO3SE.Defaults.es_init_processes import external_state_init_processes
from pyDO3SE.Defaults.default_processes import full_model_processes, daily_process_list, hourly_processes, daily_start_processes, daily_end_processes
from pyDO3SE.Defaults.config_init_processes import config_init_processes


# %%
"""Time log helpers"""


def accumulator(acc, v):
    acc[v[0]] = acc[v[0]] + v[1] if v[0] in acc else 0
    return acc


def grouper(acc, v):
    if v[0] not in acc:
        acc[v[0]] = [v[1]]
    else:
        acc[v[0]].append(v[1])
    return acc


def averager(acc, v):
    acc[v[0]] = sum(v[1]) / len(v[1])
    return acc


def counter(acc, v):
    acc[v[0]] = len(v[1])
    return acc


# %%
"""Initial Setup"""

DEMO_START_DAY = 0
DEMO_END_DAY = 40
HOURS_IN_DAY = 24
DAYS_IN_YEAR = 365
DAYS_TO_RUN = DEMO_END_DAY - DEMO_START_DAY
HOURS_IN_YEAR = HOURS_IN_DAY * DAYS_IN_YEAR
HOURS_TO_RUN = HOURS_IN_DAY * DAYS_TO_RUN
EXT_DATA_COLS = [
    # TODO: This should be based on the config
    'PAR',
    'VPD',
    'Ts_C',
    'u',
    'P',
    'O3',
    'dd',
    'hr',
    'precip',
]


# %%
start_day = DEMO_START_DAY
end_day = DEMO_END_DAY
config_location = 'examples/spanish_wheat/spanish_wheat_config_for_short_test.json'
data_location = 'examples/spanish_wheat/spanish_wheat_data.csv'

# %%
"""SETUP CONFIG"""

config = config_loader(config_location, 'json')
process_runner = ProcessRunner(config, DEBUG_MODE=True)

config_amended = process_runner.run_processes(
    config_init_processes(config),
    config)
process_runner.config = config_amended
unpack(config_amended)

# %%
"""Setup External state"""

external_state_data = load_external_state(data_location, 'csv', EXT_DATA_COLS)
process_runner.external_state = external_state_data
external_state = process_runner.run_processes(
    external_state_init_processes(start_day, end_day, config),
    external_state_data)
process_runner.external_state = external_state

# json.dumps(prep_data_for_snapshot(asdict(process_runner.external_state)),
#                                      cls=NumpyEncoder, indent=4, sort_keys=True)


# %%
assert external_state.sinB[0] is not None
assert process_runner.external_state.sinB[0] is not None

# %%
"""Setup initial state"""

empty_state = Model_State_Shape()
state_init = process_runner.run_processes(
    state_init_processes(config, start_day, end_day),
    empty_state,
)
unpack(state_init)


# %%
assert state_init.canopy_component[0].season_Astart_td is not None
assert process_runner.external_state.sinB[0] is not None

# %%
"""Setup daily start process"""
day_start_state = process_runner.run_processes(
    flatten_list(daily_start_processes(config, 0)),
    state_init,
)

# %%
"""Run for 24 hours"""
process_runner.time_logs = []


# %%
# Check time for hour run
def run_hours(process_runner):
    current_state = day_start_state
    for hr in range(0, 24):
        current_state = process_runner.run_processes(
            flatten_list(hourly_processes(config, 0, hr)),
            current_state,
        )
    return current_state


state_after_run_24_hours = run_hours(process_runner)

# %%
tlgs = process_runner.time_logs

# %%
"""Get time info"""
grouped_times = reduce(grouper, tlgs, {})
average_times = reduce(averager, grouped_times.items(), {})
total_times = reduce(accumulator, tlgs, {})
process_counts = reduce(counter, grouped_times.items(), {})


# %%
"""Chart - average times"""

items = sorted(average_times.items(), key=lambda tup: tup[1])
len(items)
x = np.arange(len(items))
plt.bar(x, [l[1] for l in items])
plt.xticks(x, [l[0] for l in items], rotation='vertical')
plt.bar


# %%
"""Chart - process_counts"""
items = sorted(process_counts.items(), key=lambda tup: tup[1])
len(items)
x = np.arange(len(items))
plt.bar(x, [l[1] for l in items])
plt.xticks(x, [l[0] for l in items], rotation='vertical')
plt.bar


# %%
"""Chart - total_times"""
items = sorted(total_times.items(), key=lambda tup: tup[1])
len(items)
x = np.arange(len(items))
plt.bar(x, [l[1] for l in items])
plt.xticks(x, [l[0] for l in items], rotation='vertical')
plt.bar


# %%
sorted(list(total_times.items()), key=lambda tup: -tup[1])[0:10]


# %%
"""Test ewert on its own at day 1"""
iLC = 0
iL = 0
ewert_process = Process(
    func=ewert,
    config_inputs=lambda config: [
        # TODO: t_l etc should come from phenology output
        I(config.Land_Cover.parameters[iLC].pn_gsto.t_lse_constant,
          as_='t_lse_constant'),
        I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_1, as_='gamma_1'),
        I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_2, as_='gamma_2'),
        I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_3, as_='gamma_3'),

        I(config.Land_Cover.parameters[iLC].pn_gsto.g_sto_0, as_='g_sto_0'),
        I(config.Land_Cover.parameters[iLC].pn_gsto.m, as_='m'),

        I(config.Land_Cover.parameters[iLC].Lm, as_='Lm'),
    ],
    state_inputs=lambda state: [
        I(state.canopy_component[iLC].t_l_estimate, as_='t_l_estimate'),
        I(state.canopy_component[iLC].t_lem, as_='t_lem'),
        I(state.canopy_component[iLC].t_lma, as_='t_lma'),
        I(state.canopy_component[iLC].t_lep, as_='t_lep'),
        I(state.canopy_component[iLC].t_lse, as_='t_lse'),

        I(state.canopy_layer_component[iL][iLC].V_cmax_25, as_='V_cmax_25'),
        I(state.canopy_layer_component[iL][iLC].J_max_25, as_='J_max_25'),

        I(state.canopy_layers[iL].micro_met.PARsun, as_='PARsun'),
        I(state.canopy_layers[iL].micro_met.PARshade, as_='PARshade'),
        I(state.canopy_layer_component[iL][iLC].LAI, as_='LAI'),
        I(state.canopy_layer_component[iL][iLC].D_0, as_='D_0'),
        I(state.canopy_layer_component[iL][iLC].g_bv, as_='g_bv'),
        I(state.canopy_layer_component[iL][iLC].td_dd, as_='td_dd'),
        I(state.prev_hour.canopy_layer_component[iL][iLC].td_dd, as_='td_dd_prev'),
        I(state.prev_hour.canopy_layer_component[iL][iLC].fO3_d, as_='fO3_d_prev'),
        I(state.prev_hour.canopy_layer_component[iL][iLC].O3up, as_='O3up_prev'),
        I(state.prev_hour.canopy_layer_component[iL][iLC].O3up_acc,
          as_='O3up_acc_prev'),
        I(state.canopy_layers[iL].micro_met.micro_u, as_='uh'),
    ],
    external_state_inputs=lambda e_state, row_index: [
        I(e_state.Ts_C[row_index], as_='Tleaf_C'),
        I(e_state.O3_nmol[row_index], as_='O3_nmol'),
        I(e_state.eact[row_index], as_='eact'),
        I(e_state.sinB[row_index], as_='sinB'),
        I(e_state.CO2[row_index], as_='c_a'),
        I(e_state.is_daylight[row_index], as_='is_daylight'),
    ],
    additional_inputs=lambda: [
        I(2500, as_='Rext'),  # TODO: Check this is const
    ],
    state_outputs=[
        I('Tleaf_C', as_=f'canopy_layer_component.{iL}.{iLC}.Tleaf_C_estimate'),
        I('g_sv', as_=f'canopy_layer_component.{iL}.{iLC}.g_sv'),
        I('A_n', as_=f'canopy_layer_component.{iL}.{iLC}.A_n'),
        I('A_c', as_=f'canopy_layer_component.{iL}.{iLC}.A_c'),
        I('A_j', as_=f'canopy_layer_component.{iL}.{iLC}.A_j'),
        I('A_p', as_=f'canopy_layer_component.{iL}.{iLC}.A_p'),
        I('A_n_limit_factor',
          as_=f'canopy_layer_component.{iL}.{iLC}.A_n_limit_factor'),
        I('R_d', as_=f'canopy_layer_component.{iL}.{iLC}.R_d'),
        I('O3up_out', as_=f'canopy_layer_component.{iL}.{iLC}.O3up'),
        I('O3up_acc_out', as_=f'canopy_layer_component.{iL}.{iLC}.O3up_acc'),
        I('fO3_h_out', as_=f'canopy_layer_component.{iL}.{iLC}.fO3_h'),
        I('fO3_d_out', as_=f'canopy_layer_component.{iL}.{iLC}.fO3_d'),
        I('c_i', as_=f'canopy_layer_component.{iL}.{iLC}.c_i'),
        I('Canopy_A_n', as_=f'canopy_layer_component.{iL}.{iLC}.Canopy_A_n'),
    ]
)


def run_ewert(process_runner):
    return process_runner.run_processes(
        [ewert_process],
        state_after_run_24_hours,
    )


process_runner.time_logs = []
run_ewert(process_runner)
None
# %%

process_runner.time_logs

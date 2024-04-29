# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
"""Test the full model."""

from pstats import SortKey
import pstats
import cProfile
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

# %% [markdown]
#  # Initial Setup

# %%

DEMO_START_DAY = 0
DEMO_END_DAY = 365
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
config_location = 'examples/spanish_wheat/spanish_wheat_config.json'
data_location = 'examples/spanish_wheat/spanish_wheat_data.csv'

# %% [markdown]
# ## Setup config

# %%
# == 1. SETUP CONFIG

config = config_loader(config_location, 'json')
process_runner = ProcessRunner(config, DEBUG_MODE=False)

config_amended = process_runner.run_processes(
    config_init_processes(config),
    config)
process_runner.config = config_amended
unpack(config_amended)

# %% [markdown]
# ## Setup External state

# %%
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

# %% [markdown]
# ## Setup initial state

# %%

# == SETUP INITIAL STATE
initial_state = Model_State_Shape()
state_init = process_runner.run_processes(
    state_init_processes(config, start_day, end_day),
    initial_state,
)
unpack(state_init)


# %%
assert state_init.canopy_component[0].season_Astart_td is not None
assert process_runner.external_state.sinB[0] is not None

# %% [markdown]
# # Test for single hour
#
# %% [markdown]
# ## Setup daily start process

# %%
day_start_state = process_runner.run_processes(
    flatten_list(daily_start_processes(config, 0)),
    state_init,
)

# %% [markdown]
# ## run for hour 0

# %%
state_post_hour_0 = process_runner.run_processes(
    flatten_list(hourly_processes(config, 0, 0)),
    day_start_state,
)


# %%
unpack(state_post_hour_0.canopy)


# %%
# Check time for hour run
def run_hour_0_for_timeit():
    process_runner.run_processes(
        flatten_list(hourly_processes(config, 0, 0)),
        day_start_state,
    )


t0 = timeit(lambda: run_hour_0_for_timeit(), number=1)
print('t0', t0)

# %% [markdown]
# ## Run for 24 hours

# %%
# Check time for day run


def run_hours():
    current_state = day_start_state
    for hr in range(0, 24):
        current_state = process_runner.run_processes(
            flatten_list(hourly_processes(config, 0, hr)),
            current_state,
        )
    return current_state


state_after_run_24_hours = run_hours()

t1 = timeit(lambda: run_hours(), number=1)
print('t1', t1)

# %% [markdown]
# ------
# %% [markdown]
# # Test time per of running full model

# %%


def run_model_full():
    process_runner.reset()
    current_state = state_init
    total_processes_ran = 0
    for dd in range(0, DEMO_END_DAY):
        day_start_processes = flatten_list(daily_start_processes(config, dd))
        total_processes_ran += len(day_start_processes)
        current_state = process_runner.run_processes(
            day_start_processes,
            current_state,
        )
        for hr in range(0, 24):
            hour_processes = flatten_list(hourly_processes(config, dd, hr))
            total_processes_ran += len(hour_processes)
            current_state = process_runner.run_processes(
                hour_processes,
                current_state,
            )
        day_end_processes = flatten_list(daily_end_processes(config, dd))
        total_processes_ran += len(day_end_processes)
        current_state = process_runner.run_processes(
            day_end_processes,
            current_state,
        )
    print(total_processes_ran)
    return current_state


# %%
run_model_full()
# %%
t2 = timeit(lambda: run_model_full(), number=1)
print('t2', t2)
# 12.23s for 40 days (0.30s per day) with 51920 processes (1298 per day)
# 382.41s for 365 (1.04s per day) with 464580 processes (1272 per day)
# 41.6s for 365 days removed immutability
# 23.58s for 365 days proflow 0.3.3
# %%

# %%
# store time per hour
time_per_hour = []
time_per_day = []
def day_start_processes(dd): return flatten_list(daily_start_processes(config, dd))
def hour_processes(dd, hr): return flatten_list(hourly_processes(config, dd, hr))


def run_day_start(dd, s): return process_runner.run_processes(day_start_processes(dd), s)
def run_hour(dd, hr, s): return process_runner.run_processes(hour_processes(dd, hr), s)
def day_end_processes(dd): return flatten_list(daily_end_processes(config, dd))
def run_day_end(dd, s): return process_runner.run_processes(day_end_processes(dd), s)


def run_model_full_simple():
    process_runner.reset()
    current_state = state_init
    for dd in range(0, 365):
        current_state = run_day_start(dd, current_state)
        for hr in range(0, 24):
            current_state = run_hour(dd, hr, current_state)
        current_state = run_day_end(dd, current_state)
    return current_state
# t3=timeit(lambda: run_model_full_simple(), number=1)
# t3


profile_output_file = 'notebooks/full_model/run_model_full_stats'
cProfile.run('run_model_full_simple()', profile_output_file)

# %%
profile_output_file = 'notebooks/full_model/run_model_full_stats'
p = pstats.Stats(profile_output_file)
p.strip_dirs().sort_stats(SortKey.CALLS).print_stats(30)

# %%
p = pstats.Stats(profile_output_file)
p.strip_dirs().sort_stats(SortKey.TIME).print_callers(30)


# # %% [markdown]
# # -----

# # %%
# final_state = process_runner.run_processes(
#     full_model_processes(config, start_day, end_day),
#     state_init,
# )


# # %%


# # %%
# def run_hourly_processes():
#     process_runner.run_processes(
#         flatten_list(hourly_processes(config, 0, 0)),
#         state_init,
#     )


# # %%
# t3=timeit(lambda: run_hourly_processes(), number=1)
# print('t3', t3)


# # %%
# final_state_post_setup = process_runner.run_processes(
#     flatten_list(daily_process_list(config, 0)),
#     state_init,
# )


# # %%
# for p in flatten_list(daily_process_list(config, 0)):
#     print(p.comment or f'func: {p.func.__name__}')


# # %%
# unpack(final_state_post_setup.temporal)


# # %%
# unpack(final_state_post_setup.canopy)


# # %%


# %%

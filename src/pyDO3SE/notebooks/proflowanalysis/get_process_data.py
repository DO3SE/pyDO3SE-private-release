"""Convert list of hourly processes into nodes and edges."""
# %%
from helpers.list_helpers import flatten_list

from helpers.named_tuple_helpers import unpack
from proflow.ProcessRunnerCls import ProcessRunner

from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Defaults.default_processes import hourly_processes
from pyDO3SE.Defaults.config_init_processes import config_init_processes


from proflow.analysis.network.extract_nodes_and_edges import processes_to_nodes_and_edges
from proflow.analysis.network.export import asjson
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

hourly_processes_hr0 = flatten_list(hourly_processes(config_amended, 0, 0))
len(hourly_processes_hr0)

# %%


nodes, edges = processes_to_nodes_and_edges(hourly_processes_hr0)
asjson(nodes, edges, 'data.json')


# %%
len(nodes), len(hourly_processes_hr0)

# %%
hourly_processes_hr0

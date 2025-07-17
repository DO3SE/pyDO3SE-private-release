"""Test the full model."""

from data_helpers.list_helpers import flatten_list
from data_helpers.cls_parsing import unpack
from proflow.ProcessRunnerCls import ProcessRunner
from pyDO3SE.Config.config_loader import config_loader
from pyDO3SE.Pipelines.config_init_processes import config_init_processes
from pyDO3SE.util.test_utils import process_snapshot


DEMO_START_DAY = 0
DEMO_END_DAY = 40
HOURS_IN_DAY = 24
DAYS_IN_YEAR = 365
DAYS_TO_RUN = DEMO_END_DAY - DEMO_START_DAY
HOURS_IN_YEAR = HOURS_IN_DAY * DAYS_IN_YEAR
HOURS_TO_RUN = HOURS_IN_DAY * DAYS_TO_RUN
EXT_DATA_COLS = [
    # TODO: This should be based on the config
    "PAR",
    "VPD",
    "Ts_C",
    "u",
    "P",
    "O3",
    "dd",
    "hr",
    "precip",
]

config_location = "examples/spanish_wheat/configs/spanish_wheat_config.json"
data_location = "examples/spanish_wheat/inputs/spanish_wheat_data.csv"


def test_config_init_processes(snapshot):
    """Test using the config init processes to setup config."""
    config = config_loader(config_location, config_type="json")
    process_runner = ProcessRunner(config, DEBUG_MODE=True)
    processes = config_init_processes(config)
    config_amended = process_runner.run_processes(processes, config)
    process_runner.config = config_amended
    snapshot.assert_match(process_snapshot(processes), "proccesses")
    flattened_process_comments = [p.comment or p.func.__name__ for p in flatten_list(processes)]
    snapshot.assert_match(process_snapshot(flattened_process_comments), "Process_comments")
    snapshot.assert_match(process_snapshot(config_amended), "config_amended")

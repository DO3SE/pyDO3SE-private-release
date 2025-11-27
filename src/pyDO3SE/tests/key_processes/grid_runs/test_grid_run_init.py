"""A set of tests for checking we can initialise grid runs."""
import os
import pytest
from pyDO3SE.Config.config_loader import config_loader_pickled

#
from .utils import _assertTestSetup, _setup, _run_initialization, TestSetup


@pytest.fixture(scope="class")
def before_all_init(request):
    _self = request.cls
    print("Before all")
    _assertTestSetup(_self)
    _setup(_self)
    _run_initialization(_self)


@pytest.mark.usefixtures("before_all_init")
class TestGridModelInit(TestSetup):
    multi_file_netcdf = False
    runid = "TestGridModelInit"
    project_dir = "examples/net_cdf/single_file_hour"

    def test_state_setup_correctly(self: TestSetup):
        assert os.path.exists(
            f"{self.run_paths.live_state_dir}/0_0.state"), f"State file not created in {self.run_paths.live_state_dir}/0_0.state"
        assert os.path.exists(
            f"{self.run_paths.processed_configs_dir}/0_0.config"), f"Config file not created in {self.run_paths.processed_configs_dir}/0_0.config"

        # loaded_state = model_state_loader_quick(f"{self.run_paths.live_state_dir}/0_0.state")
        # assert loaded_state.canopy.canopy_height == 0.0, f"Canopy height not set to 0.0, instead {loaded_state.canopy.canopy_height}"

    def test_should_set_sowing_date_from_lat_function(self):
        loaded_config = config_loader_pickled(f"{self.run_paths.processed_configs_dir}/0_0.config")
        assert loaded_config.Land_Cover.parameters[0].phenology.key_dates.sowing == 291
        assert loaded_config.Met.thermal_time_start == 291
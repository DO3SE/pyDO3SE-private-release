"""A set of tests for checking we can initialise grid runs for global data.


This involves checking that the current hour is linked to timezone and not constant UTC.
"""
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
class TestGlobalGridRunInit(TestSetup):
    multi_file_netcdf = False
    runid = "TestGlobalGridRunInit"
    project_dir = "examples/net_cdf/single_file_global"


    def test_should_set_hour_correctly(self):
        """This involves getting the Location.hr_offset value from the e_state_overrides file."""
        loaded_config_a = config_loader_pickled(f"{self.run_paths.processed_configs_dir}/0_0.config")
        loaded_config_b = config_loader_pickled(f"{self.run_paths.processed_configs_dir}/2_0.config")

        assert loaded_config_a.Location.hr_offset == 0
        assert loaded_config_b.Location.hr_offset == 12



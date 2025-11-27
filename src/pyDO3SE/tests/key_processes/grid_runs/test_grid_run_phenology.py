"""A set of tests to check the phenology works in grid runs."""
import os
import pytest
from pathlib import Path
from do3se_phenology.plots import plot_phenology_from_config
from do3se_phenology.units import TimeTypes
from pyDO3SE.Grid_Model.setup_grid_model import process_grid_config
from .utils import _assertTestSetup, _setup, _run_initialization, TestSetup


@pytest.fixture(scope="class")
def before_all_init(request):
    _self = request.cls
    print("Before all")
    _assertTestSetup(_self)
    _setup(_self)
    _run_initialization(_self)


@pytest.mark.usefixtures("before_all_init")
class TestGridModelPhenology(TestSetup):
    multi_file_netcdf = False
    runid = "TestGridModelPhenology"
    project_dir = "examples/net_cdf/full_season"
    config_id = "bangor_wheat"

    def test_can_produce_pre_run_phenology_plots(self):
        processed_config = process_grid_config(
            Path(self.run_paths.config_path),
            grid_overrides_file=Path(self.project_paths.e_state_overrides_file_path),
            grid_overrides_field_map_path=Path(self.run_paths.e_state_overrides_field_map_path),
            base_config_file=Path(self.project_paths.base_config_path),
            grid_coord=(0, 0),
        )

        day_count = (
            int(processed_config.Location.end_day or 0)
            - int(processed_config.Location.start_day or 0)
            if processed_config.Location.start_day is not None
            and processed_config.Location.end_day is not None
            else 365
        )
        output_directory =  Path("tests/key_processes/grid_runs/test_grid_run_phenology_outputs")
        os.makedirs(output_directory, exist_ok=True)
        plot_phenology_from_config(
            processed_config.Land_Cover.parameters[0].phenology,
            processed_config.Land_Cover.phenology_options,
            nP=processed_config.Land_Cover.nP,
            output_location=Path(f"{output_directory}/phenology_{self.run_paths.config_id}.png"),
            # TODO: Add external data input
            day_count=day_count,
            plot_dd=processed_config.Land_Cover.phenology_options.time_type
            == TimeTypes.JULIAN_DAY,
            plot_td=processed_config.Land_Cover.phenology_options.time_type
            == TimeTypes.THERMAL_TIME,
            plot_f_phen=True,
            plot_lengths=True,
            plot_carbon=processed_config.Land_Cover.phenology_options.time_type
            == TimeTypes.THERMAL_TIME,
            plot_growing=processed_config.Land_Cover.phenology_options.time_type
            == TimeTypes.THERMAL_TIME,
        )


@pytest.mark.usefixtures("before_all_init")
class TestGridModelPhenologyMultiplicative(TestGridModelPhenology):
    config_id = "bangor_wheat_multiplicative"

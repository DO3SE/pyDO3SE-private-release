from copy import deepcopy
from dataclasses import replace

from do3se_phenology.plots import plot_phenology_from_config


from do3se_phenology.switchboard import (
    PhenologyMethods,
)
from do3se_phenology.config import (
    ModelConfig,
)
from do3se_phenology.presets.wheat import (
    SpringWheat,
)
from do3se_phenology.units import *


class TestPlotPhenologyFromConfig:

    def test_should_run_without_errors(self):
        model_config = ModelConfig(
            phenology_method=PhenologyMethods.SEASON_FRACTION,
        )
        nP = 3
        species_config = replace(
            deepcopy(SpringWheat)
            # deepcopy(default_expected_species_output(nP))
        )
        species_config.key_dates.sowing = 20
        species_config.key_dates_td.sowing = 0
        species_config.key_lengths_td.sowing_to_end = 2000
        species_config.f_phen_min = 0.1

        plot_phenology_from_config(
            species_config,
            model_config,
            nP=nP,
            output_location="tests/outputs/phenology_spring_wheat.png",
            plot_dd=True,
        )

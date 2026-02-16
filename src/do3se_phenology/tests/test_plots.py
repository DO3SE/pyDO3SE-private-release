from copy import deepcopy
from dataclasses import replace
from pathlib import Path

from do3se_phenology.plots import plot_phenology_from_config
from do3se_phenology.switchboard import (
    PhenologyMethods,
)
from do3se_phenology.config import (
    ModelConfig,
    PlantEmergeMethod,
    DVIMethods,
)
from do3se_phenology.presets.wheat import (
    SpringWheat,
    SpringWheatMultiplicative,
)
from do3se_phenology.units import TimeTypes


class TestPlotPhenologyFromConfig:
    def test_plot_thermal_time_plots(self):
        model_config = ModelConfig(
            phenology_method=PhenologyMethods.SEASON_FRACTION,
            dvi_method=DVIMethods.JULES,
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
            output_location=Path("tests/outputs/phenology_spring_wheat_td.png"),
            plot_td=True,
            plot_dd=False,
        )

    def test_plot_julian_day_plots(self):
        model_config = ModelConfig(
            phenology_method=PhenologyMethods.LEGACY_DAY_PLF,
            time_type=TimeTypes.JULIAN_DAY,
            plant_emerge_method=PlantEmergeMethod.FPHEN,
            dvi_method=DVIMethods.DISABLED,
        )
        nP = 1
        species_config = replace(
            deepcopy(SpringWheatMultiplicative),
        )

        plot_phenology_from_config(
            species_config,
            model_config,
            nP=nP,
            output_location=Path("tests/outputs/phenology_spring_wheat_dd.png"),
            plot_td=False,
            plot_dd=True,
            plot_carbon=False,
            plot_growing=False,

        )


    def test_plot_julian_day_multi_year_plots(self):
        model_config = ModelConfig(
            phenology_method=PhenologyMethods.LEGACY_DAY_PLF,
            time_type=TimeTypes.JULIAN_DAY,
            plant_emerge_method=PlantEmergeMethod.FPHEN,
            dvi_method=DVIMethods.DISABLED,
            allow_cross_year_plf_phenology=True
        )
        nP = 1
        species_config = replace(
            deepcopy(SpringWheatMultiplicative),
            key_dates=replace(
                SpringWheatMultiplicative.key_dates,
                harvest=20,
                sowing=200,
            )
        )

        fig, axss = plot_phenology_from_config(
            species_config,
            model_config,
            nP=nP,
            output_location=Path("tests/outputs/phenology_spring_wheat_multi_year_dd.png"),
            plot_td=False,
            plot_dd=True,
            plot_carbon=False,
            plot_growing=False,
            day_count=365*2
        )

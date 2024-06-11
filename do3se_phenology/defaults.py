from dataclasses import replace
from copy import deepcopy

from .config import (
    ModelConfig,
    SpeciesConfig,
    SpeciesPresets,
    Wheat,
    PhenologyMethods,
)


example_species = SpeciesConfig(
    PRESET=SpeciesPresets.WHEAT_SPRING,
    f_tt_emr=0.05,
    f_tt_veg=700 / 2000,
    f_tt_rep=1200 / 2000,
    f_phen_min=0.1,
)


class ConfigDefault():
    model_config = ModelConfig()
    species_config = replace(Wheat)


class FPhenThermalTime(ConfigDefault):
    model_config = ModelConfig(
        phenology_method=PhenologyMethods.FPHEN_THERMAL_TIME,
    )
    species_config = deepcopy(example_species)
    species_config.leaf_f_phen_a=0.3,
    species_config.leaf_f_phen_b=0.7,

    species_config.key_lengths_td.sowing_to_emerge=100.0
    species_config.key_lengths_td.sowing_to_f_phen_b=400.0
    species_config.key_lengths_td.sowing_to_f_phen_c=1240.0
    species_config.key_lengths_td.sowing_to_end=2000.0

    species_config.key_dates.sowing = 20
    species_config.key_dates_td.Astart = 1080

    species_config.key_lengths_flag_leaf_td.leaf_f_phen_e=160.0,  # 0.08
    species_config.key_lengths_flag_leaf_td.leaf_f_phen_g=100.0,  # 0.05
    species_config.key_lengths_flag_leaf_td.leaf_f_phen_h=440.0,  # 0.22
    species_config.key_lengths_flag_leaf_td.leaf_f_phen_i=760.0,  # 0.38
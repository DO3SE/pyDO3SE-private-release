from do3se_phenology.config import SpeciesPresets
from do3se_phenology.presets.wheat import SpringWheat
from do3se_phenology.presets.forest import ForestGeneric, TemperateMixedForest

SpeciesPresetsParams = {
    SpeciesPresets.WHEAT_SPRING: SpringWheat,
    SpeciesPresets.FOREST_GENERIC: ForestGeneric,
    SpeciesPresets.TEMPERATE_MIXED_FOREST: TemperateMixedForest,
}

from do3se_phenology.config import (
    SpeciesConfig,
    ModelConfig,
    SAILCMethods,
    FPhenMethods,
    LeafFPhenMethods,
    DayFphenPlf,
)

ForestGeneric = SpeciesConfig(
    SAI_method=SAILCMethods.LAI_max,
    f_phen_method=FPhenMethods.SIMPLE_DAY_PLF,
    leaf_f_phen_method=LeafFPhenMethods.F_PHEN,

)


ForestModelConfig = ModelConfig(
    flag_leaf_only=False,
    phenology_method="day plf",
    dvi_method="disabled",
    LAI_method="estimate total",
    time_type="julian_day",
    sgs_time_type="julian_day",
    sgs_key_day="sowing_day",
    zero_day="sowing",
    plant_emerge_method="fphen",
    flag_leaf_emerge_method="constant",
    use_vernalisation=False,
    use_photoperiod_factor=False,
    sowing_day_method="LATITUDE_FOREST",
)


TemperateMixedForest = SpeciesConfig(
    SAI_method=SAILCMethods.LAI_max,
    f_phen_method=FPhenMethods.SIMPLE_DAY_PLF,
    leaf_f_phen_method=LeafFPhenMethods.F_PHEN,
    day_fphen_plf=DayFphenPlf(
        f_phen_a=0.0,
        f_phen_c=1.0,
        f_phen_e=0.0,
        f_phen_1=20.0,
        f_phen_4=20.0,
    ),
    LAI_a=0.4,
    LAI_b=1.0,
    LAI_c=1.0,
    LAI_d=1.0,
    LAI_1=1.0,
    LAI_2=1.0,
)

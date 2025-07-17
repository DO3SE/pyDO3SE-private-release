from do3se_phenology.config import SpeciesConfig

# NOTE: Preset values that override presets in the dataclass will be overriden

#: SpringWheat defaults
SpringWheat = SpeciesConfig(
    f_Astart=0.54,
    f_mid_anthesis=0.62,
    f_fphen_a=0.05,
    f_fphen_b=0.20,
    f_fphen_c=0.62,
    f_fphen_d=1.0,
    f_tt_emr=0.05,
    f_tt_veg=0.95 * 800 / 2000,
    f_tt_rep=0.95 * 1200 / 2000,
    leaf_f_phen_a=0.3,
    leaf_f_phen_b=0.7,
    f_fphen_1_ets=0.08,
    f_fphen_3_ets=0.05,
    f_fphen_4_ets=0.22,
    f_fphen_5_ets=0.38,
    f_tt_fst_acc=0,
    f_t_lem=0.14,
    f_t_lma=0.46,
    f_t_lep=0.31,
    f_t_lse=0.15,
    f_t_lse_mature=0.33,
    f_phen_min=0,
    f_leaf_f_fphen=0.46,
)

#: Winter wheat defaults
Wheat = SpringWheat

"""Configuration for the Phenology Module.


Key sources for parameters
--------------------------
- Legacy Day plf runs: https://icpvegetation.ceh.ac.uk/sites/default/files/FinalnewChapter3v4Oct2017_000.pdf

Setting new Season fractions
----------------------------
The best way to set the phenology config is to work out the fractions and let the model work out everything else.
The fractions required are(with defaults)...
    {
        "f_phen_min": 0,
        "leaf_f_phen_a": 0.3,
        "leaf_f_phen_b": 0.7,
        "f_Astart": 0.54,
        "f_mid_anthesis": 0.62,
        "f_fphen_1_ets": 0.08,
        "f_fphen_3_ets": 0.05,
        "f_fphen_4_ets": 0.22,
        "f_fphen_5_ets": 0.38,
        "f_t_lem": 0.14,
        "f_t_lma": 0.46,
        "f_t_lep": 0.31,
        "f_t_lse": 0.15,
        "f_t_lse_mature": 0.33,
        "f_fphen_a": 0.05,
        "f_fphen_b": 0.2,
        "f_fphen_c": 0.62,
        "f_fphen_d": 1.0,
        "f_tt_emr": 0.05,
        "f_tt_veg": 0.38,
        "f_tt_rep": 0.57,
        "f_leaf_f_fphen": 0.46,
    }

These should all be fractions of growing season I.e from sowing to end of senescence(harvest).

As an example if you have "astart_td": 862.0,"harvest_td": 928.0
The f_Astart fraction would then be 862/928 = 0.92

A way to check it is then working is to set these fractions and run the pyDO3SE model.
In the processed config "key_dates_td" should match what you have.

# TODO: Add plots that explain the fractions here
See :module:`do3se_phenology.fphen` and :module:`do3se_phenology.phyllochron_dvi`

"""


from data_helpers.fill_np_array import fill_np_array_with_cls
from typing import List, Optional
from dataclasses import dataclass, field
from enum import Enum

from do3se_phenology.units import Fraction, PiecewiseFunction, TimeUnit, TimeTypes

DEFAULT_LC: int = 1  # default number of land covers


class SpeciesPresets(Enum):
    """Preset species configs.

    List of presets that can be used as a base.
    Each of these should have a matching preset.

    """

    WHEAT_SPRING = "WHEAT_SPRING"
    FOREST_GENERIC = "FOREST_GENERIC"
    TEMPERATE_MIXED_FOREST = "TEMPERATE_MIXED_FOREST"


# TODO: Implement accumulation period methods
# class AccumulationPeriodMethods(Enum):
#     """Methods of defining Astart and Aend.

#      - "constant":         Astart and Aend supplied
#      - "growing season":   Astart = SGS and Aend = EGS
#      - "wheat latitude":   Use wheat latitude model

#     """
#     CONSTANT = "constant"
#     GROWING_SEASON = "growing season"
#     WHEAT_LATITUDE = "wheat latitude"

class KeyDays(Enum):

    SGS = "SGS"
    EGS = "EGS"
    SOWING_DAY = "sowing_day"


class PhenologyMethods(Enum):
    """Methods of defining whole plant phenology.


    Options
    -------
    DISABLED = "disabled"
        Don't set any phenology values
    SEASON_FRACTION = "season_fraction"
        Define phenology by fractions of total growing season
    FLAG_LEAF_FRACTION = "flag_leaf_fraction"
        Define phenology by fractions of Astart to Aend
    FPHEN_THERMAL_TIME = "fphen_thermal_time"
        Define phenology using fphen thermal time intervals
    FPHEN_DATA = "fphen_data"
        Define phenology from existing fphen data
    LEAF_FPHEN_DATA = "leaf_fphen_data"
        Define phenology from existing leaf fphen data
    DEVELOPMENT_INDEX = "development index"
        Define phenology from DVI intervals
    LEGACY_DAY_PLF = "day plf"
        Use day plf functions to set phenology

    """

    DISABLED = "disabled"
    SEASON_FRACTION = "season_fraction"
    FLAG_LEAF_FRACTION = "flag_leaf_fraction"
    FPHEN_THERMAL_TIME = "fphen_thermal_time"
    FPHEN_DATA = "fphen_data"
    LEAF_FPHEN_DATA = "leaf_fphen_data"
    DEVELOPMENT_INDEX = "development index"
    LEGACY_DAY_PLF = "day plf"


class FPhenMethods(Enum):
    """Methods to defne f_phen.

    Options
    -------
    DISABLED = "disabled"
        f_phen supplied (or left at default value of 1.0)
    SIMPLE_DAY_PLF = "simple day PLF"
        "single hump" method using 3 values (_a, _c and
                           _e) and 2 slope periods (_1 and _4)
    COMPLEX_DAY_PLF = "complex day PLF"
        "double hump" method using 5 values (_a to _e)
                           and 4 slopes: _1 and _4 for start and end, _2
                           and _3 for the middle section, starting at _limA
                           and ending at _limB
    TT_DAY_PLF = "tt day PLF"
        Simple day PLF using thermal time values for intervals
    TT_GROWING_SEASON = "tt growing season"
        Define t_sgs and t_egs and use species fractions to define other values

    """
    DISABLED = "disabled"
    INPUT = "input"
    SIMPLE_DAY_PLF = "simple day PLF"
    COMPLEX_DAY_PLF = "complex day PLF"
    TT_DAY_PLF = "tt day PLF"
    TT_GROWING_SEASON = "tt growing season"


class SowingDateMethods(Enum):
    """Methods for defining the sowing date.

    Options
    -------
    INPUT = "INPUT"
        Sowing date is provided as input
    LATITUDE = "LATITUDE"
        DEPRECATED: Use LATITUDE_SPRING_EUROPE or LATITUDE_WINTER_CHINA instead
    LATITUDE_FOREST = "LATITUDE_FOREST"
        Sowing date set from latitude. Valid for European Forest
    LATITUDE_SPRING_EUROPE = "LATITUDE_SPRING_EUROPE"
        Sowing date set from latitude. Valid for European Spring Wheat
    LATITUDE_WINTER_CHINA = "LATITUDE_WINTER_CHINA"
        Sowing date set from latitude. Valid for Chinese Winter Wheat
    SKIP = "SKIP"
        Sowing day not required(E.g. Astart of flag only)

    """
    INPUT = "INPUT"
    LATITUDE = "LATITUDE"
    LATITUDE_FOREST = "LATITUDE_FOREST"
    LATITUDE_SPRING_EUROPE = "LATITUDE_SPRING_EUROPE"
    LATITUDE_WINTER_CHINA = "LATITUDE_WINTER_CHINA"
    SKIP = "SKIP"


class LeafFPhenMethods(Enum):
    """Methods to define leaf_f_phen in Multiplicative model.

    Options
    -------
    DISABLED = "disabled"
        leaf_f_phen supplied (or left at default value of 1.0)
    INPUT = "input"
        leaf_f_phen = input
    F_PHEN = "f_phen"
        leaf_f_phen = f_phen
    DAY_PLF = "day PLF"
        "single hump" PLF between Astart and Aend
    TT_DAY_PLF = "tt day PLF"
        Simple day PLF using thermal time values for intervals
    TT_GROWING_SEASON = "tt growing season"
        Define t_sgs and t_egs and use species fractions to define other values

    """
    DISABLED = "disabled"
    INPUT = "input"
    F_PHEN = "f_phen"
    DAY_PLF = "day PLF"
    TT_DAY_PLF = "tt day PLF"
    TT_GROWING_SEASON = "tt growing season"


class LifeSpanMethods(Enum):
    """Life span method for calculating t_l etc.

    Options
    -------
    CONSTANT = "constant"
        Use config parameters
    JULES = "JULES"
        Use JULES dvi phyllochron phenology method
    LEAF_F_PHEN = "leaf_f_phen"
        Estimate the life stages using input leaf_f_phen data
    T_LEAF_F_PHEN = "t_leaf_f_phen"
        Estimate the life stages using input leaf_f_phen thermal time intervals
    GROWING_SEASON = "growing_season"
        Estimate the life stages when given growing season

    """
    CONSTANT = "constant"
    JULES = "JULES"
    LEAF_F_PHEN = "leaf_f_phen"
    T_LEAF_F_PHEN = "t_leaf_f_phen"
    GROWING_SEASON = "growing_season"


class PlantEmergeMethod(Enum):
    """Methods of defining if a plant has emerged.

    Options
    -------
    CONSTANT = "constant"
        t_emerg or emerge_day set in config
    SGS = "SGS"
        emerge_day set to equal SGS
    DVI = "dvi"
        Plant emerges at DVI > 0
    FPHEN = "fphen"
        Plant emerges at fphen > 0
    LEAF_F_PHEN = "leaf_f_phen"
        Calculate emerge_t from fphen data

    """

    CONSTANT = "constant"
    SGS = "SGS"
    DVI = "dvi"
    FPHEN = "fphen"
    LEAF_F_PHEN = "leaf_f_phen"

# @deprecated
# class GrowingSeasonMethod(Enum):
#     """Methods of defining start of growing season.

#     Growing season is the time span which we run calculations.

#     Options
#     -------
#     CONSTANT = "constant"
#         t_sgs or SGS set in config
#     ASTART = "Astart"
#         Reverse calculate from Astart date flag leaf emergence.
#     FPHEN_DATA = "fphen_data"
#         Start growing season at fphen > 0
#     LEAF_F_PHEN_DATA = "leaf_f_phen_data"
#         Reverse calculate growing season from leaf_fphen > 0 == Astart

#     """

#     CONSTANT = "constant"
#     ASTART = "ASTART"
#     FPHEN_DATA = "NOT_IMPLEMENTED"  # "fphen_data"
#     LEAF_F_PHEN_DATA = "NOT_IMPLEMENTED"  # "leaf_f_phen_data"


@dataclass
class PhenologyKeyDates:
    """Key phenology dates.

    Days are relative to zero day.

    Parameters
    ----------

    sowing: TimeUnit
        Sowing day or SGS
    emergence: TimeUnit
        The start of the growing season/emergence date
    harvest: TimeUnit
        The end of the growing season/harvest date or EGS
    Astart: TimeUnit
        The start of leaf accumulation period(leaf f phen)
    Aend: TimeUnit
        The end of leaf accumulation period(leaf f phen)
    mid_anthesis: TimeUnit
        The mid anthesis

    """

    sowing: Optional[TimeUnit] = None
    emergence: Optional[TimeUnit] = None
    harvest: Optional[TimeUnit] = None
    Astart: Optional[TimeUnit] = None
    Aend: Optional[TimeUnit] = None
    mid_anthesis: Optional[TimeUnit] = None


@dataclass
class PhenologyKeyLengths:
    """Key phenology time periods.


    Parameters
    ----------
    sowing_to_emerge: TimeUnit
        sowing to canopy emergence(f_phen_a/tt_emr/t_f_phen_a)
    sowing_to_fst_acc: TimeUnit = None
        sowing to first accumulation of ozone
    sowing_to_f_phen_b: TimeUnit
        sowing to fphen b (t_f_phen_b)
    sowing_to_f_phen_c: TimeUnit
        sowing to fphen c (t_f_phen_c)
    sowing_to_end: TimeUnit
        sowing to fphen d (t_f_phen_d)
    emerg_to_end: TimeUnit
        Canopy emergence to end of growing season
    emerg_to_veg: TimeUnit
        Canopy emergence to vegitation/flowering stage(tt_veg)
    veg_to_harvest: TimeUnit
        Canopy flowering to maturity/harvest(tt_rep)

    """

    sowing_to_emerge: Optional[TimeUnit] = None
    sowing_to_f_phen_b: Optional[TimeUnit] = None
    sowing_to_f_phen_c: Optional[TimeUnit] = None
    sowing_to_astart: Optional[TimeUnit] = None
    sowing_to_end: Optional[TimeUnit] = None
    emerg_to_astart: Optional[TimeUnit] = None
    emerg_to_end: Optional[TimeUnit] = None
    emerg_to_veg: Optional[TimeUnit] = None
    veg_to_harvest: Optional[TimeUnit] = None


@dataclass
class PhenologyLeafKeyLengths:
    """Phenology key leaf lengths.

    Parameters
    ----------

    """
    tl: Optional[TimeUnit] = None
    tl_em: Optional[TimeUnit] = None
    tl_ma: Optional[TimeUnit] = None
    tl_ep: Optional[TimeUnit] = None
    tl_se: Optional[TimeUnit] = None

    #: Thermal time between Astart and leaf_f_phen_e(Mid anthesis) [degC Days]
    leaf_f_phen_e: Optional[TimeUnit] = None
    #: Thermal time between Mid Anthesis and start of seed setting [degC Days]
    leaf_f_phen_g: Optional[TimeUnit] = None
    #: Thermal time between Mid anthesis and start of senesence [degC Days]
    leaf_f_phen_h: Optional[TimeUnit] = None
    #: Thermal time between Mid anthesis and end of senescence [degC Days]
    leaf_f_phen_i: Optional[TimeUnit] = None

    #: Thermal time between plant emergence and leaf emergence
    plant_emerg_to_leaf_emerg: Optional[TimeUnit] = None
    #: Thermal time between plant emergence and leaf fst acc
    leaf_emerg_to_leaf_fst_acc: Optional[TimeUnit] = None

    # Only used in Julian Day mode
    leaf_emerg_to_fully_grown: Optional[TimeUnit] = None
    fully_grown_to_senescence: Optional[TimeUnit] = None


@dataclass
class PhenologyPLFs:
    fphen_intervals: Optional[PiecewiseFunction] = None
    leaf_fphen_intervals: Optional[PiecewiseFunction] = None
    dvi_interval: Optional[PiecewiseFunction] = None


class DVIMethods(Enum):
    """DVI methods

     - "disabled": Not calculated
     - "JULES": Use JULES thermal time method
     - "INPUT": Use external input (NOT IMPLEMENTED)

    """

    DISABLED = "disabled"
    JULES = "JULES"
    INPUT = "INPUT"


class SAILCMethods(Enum):
    """Canopy Land cover estimate height methods.

    Options
    -------
    LAI = "LAI"
        Make SAI = LAI
    LAI_max = "LAI_max"
        Make SAI = LAI_max
    LAI_BROWN_GREEN = "LAI_BROWN_GREEN"
        Make SAI = LAI brown and green leaf
    FOREST="forest"
        LEGACY Forest method: SAI = LAI + 1
    WHEAT="wheat"
        LEGACY Based on wheat lifecycle, requires LAI PLF parameters

    """

    LAI = "LAI"
    LAI_max = "LAI_max"
    LAI_BROWN_GREEN = "LAI_BROWN_GREEN"
    FOREST = "forest"
    WHEAT = "wheat"


class LAIMethods(Enum):
    """Canopy LAI methods.

     - "constant":              LAI provided as a constant value
     - "input":                 LAI provided in external inputs
     - "input total":           total LAI provided and split between layers
     - "estimate input":        estimate according to land cover growing season
     - "dvi limited constant":  Constant LAI but set to 0 outside DVI growing season
     - "carbon":                LAI calculated from carbon pools

    Additional Required Parameters
    ------------------------------

    constant
    ^^^^^^^^
        - config.config_land_cover.LAI
        - config.Land_Cover.fLAI[iL][iLC]

    input
    ^^^^^
        - external_data.LAI

    input total
    ^^^^^^^^^^
        - external_data.LAI
        - config.Land_Cover.fLAI[iL][iLC]

    estimate total
    ^^^^^^^^^^^^^^
        - config.Land_Cover.parameters[0].season.EGS
        - config.Land_Cover.parameters[0].season.LAI_1
        - config.Land_Cover.parameters[0].season.LAI_2
        - config.Land_Cover.parameters[0].season.LAI_a
        - config.Land_Cover.parameters[0].season.LAI_b
        - config.Land_Cover.parameters[0].season.LAI_c
        - config.Land_Cover.parameters[0].season.LAI_1
        - config.Land_Cover.parameters[0].season.LAI_d
        - config.Land_Cover.fLAI
        - config.Land_Cover.nL
        - config.Land_Cover.nLC

    dvi limited constant
    ^^^^^^^^^^^^^^^^^^^^
        - config.Land_Cover.fLAI[iL][iLC]
        - config.Land_Cover.LAI

    carbon
    ^^^^^^
        config.config_land_cover.LAI


    """

    CONSTANT = "constant"
    INPUT = "input"
    INPUT_TOTAL = "input total"
    ESTIMATE_TOTAL = "estimate total"
    DVI_LIMITED_CONSTANT = "dvi limited constant"
    CARBON = "carbon"


class ZeroDayOptions(Enum):
    """Options for the thermal time zero day.

    This is the day that we set thermal time to 0.
    Data before this day will either have thermal time of 0 or negative value.

    Options
    -------
    SOWING = "sowing"
        [description]
    EMERGENCE = "NOT_IMPLEMENTED"  # "emergence"
        [description]
    DATA_START = "NOT_IMPLEMENTED"  # "data_start"
        [description]
    ASTART = "Astart"
        [description]

    """

    SOWING = "sowing"
    EMERGENCE = "NOT_IMPLEMENTED"  # "emergence"
    DATA_START = "DATA_START"  # "data_start"
    ASTART = "Astart"


@dataclass
class DayFphenPlf:
    """Parameters for using fphen day plf functions.

    For documentation on how these parameters are used go to the
    :module:`do3se_phenology.f_phen`

    """

    # f_phen
    f_phen_limA: Optional[int] = None   #: Start of soil water limitation
    f_phen_limB: Optional[int] = None   #: End of soil water limitation
    f_phen_a: Optional[float] = None          #: f_phen at SGS
    f_phen_b: Optional[float] = None
    f_phen_c: Optional[float] = None
    f_phen_d: Optional[float] = None
    f_phen_e: Optional[float] = None          #: f_phen at EGS
    f_phen_1: Optional[int] = None      #: Time from f_phen_a to f_phen_b [days]
    f_phen_2: Optional[int] = None      #: Time from f_phen_b to f_phen_c [days]
    f_phen_3: Optional[int] = None      #: Time from f_phen_c to f_phen_d [days]
    f_phen_4: Optional[int] = None      #: Time from f_phen_d to f_phen_e [days]

    # leaf_f_phen
    leaf_f_phen_a: Optional[float] = None       #: f_phen at Astart
    leaf_f_phen_b: Optional[float] = None       #: f_phen at mid-season peak
    leaf_f_phen_c: Optional[float] = None       #: f_phen at Aend
    leaf_f_phen_1: Optional[int] = None   #: Time from _a to _b [days]
    leaf_f_phen_2: Optional[int] = None   #: Time from _b to _c [days]


@dataclass
class SpeciesConfig:
    """Species specific parameters.


    Piecewise intervals are assumed to start at sowing day.

    """
    #: Set PRESET to use preset values
    PRESET: Optional[SpeciesPresets] = None

    fphen_intervals: Optional[PiecewiseFunction] = None
    leaf_fphen_intervals: Optional[PiecewiseFunction] = None
    dvi_interval: Optional[PiecewiseFunction] = None

    key_dates: PhenologyKeyDates = field(default_factory=lambda: PhenologyKeyDates())
    key_dates_td: PhenologyKeyDates = field(default_factory=lambda: PhenologyKeyDates())
    key_lengths: PhenologyKeyLengths = field(default_factory=lambda: PhenologyKeyLengths())
    key_lengths_td: PhenologyKeyLengths = field(default_factory=lambda: PhenologyKeyLengths())

    key_lengths_leaf: PhenologyLeafKeyLengths = field(
        default_factory=lambda: PhenologyLeafKeyLengths())
    key_lengths_leaf_td: PhenologyLeafKeyLengths = field(
        default_factory=lambda: PhenologyLeafKeyLengths())

    key_lengths_flag_leaf: PhenologyLeafKeyLengths = field(
        default_factory=lambda: PhenologyLeafKeyLengths())
    key_lengths_flag_leaf_td: PhenologyLeafKeyLengths = field(
        default_factory=lambda: PhenologyLeafKeyLengths())

    # LAI piecewise linear function parameters
    # TODO: Merge this with LAI constant
    LAI_a: Optional[Fraction] = None       #: LAI value at SGS [m2 m-2]
    LAI_b: Optional[Fraction] = None       #: LAI value at SGS + LAI_1 [m2 m-2]
    LAI_c: Optional[Fraction] = None       #: LAI value at EGS - LAI_2 [m2 m-2]
    LAI_d: Optional[Fraction] = None       #: LAI value at EGS [m2 m-2]
    LAI_1: Optional[TimeUnit] = None   #: Time from LAI_a to LAI_b [days]
    LAI_2: Optional[TimeUnit] = None   #: Time from LAI_c to LAI_d [days]

    #: SAI method:
    SAI_method: SAILCMethods = SAILCMethods.LAI

    # FPHEN
    f_phen_method: FPhenMethods = FPhenMethods.DISABLED
    leaf_f_phen_method: LeafFPhenMethods = LeafFPhenMethods.DISABLED
    day_fphen_plf: DayFphenPlf = field(default_factory=lambda: DayFphenPlf())

    f_phen_min: Optional[TimeUnit] = None  #: thermal time between sowing day and phen_min []

    # Thermal Time leaf_f_phen
    leaf_f_phen_a: Optional[TimeUnit] = None  #: Gradient of descent during seed setting [degC Days]
    leaf_f_phen_b: Optional[TimeUnit] = None  #: Gradient of descent during senescence [degC Days]

    # Season fractions
    # Should all be fractions of growing season I.e from sowing to end of senescence(harvest)

    f_Astart: Optional[Fraction] = None
    f_mid_anthesis: Optional[Fraction] = None
    f_fphen_1_ets: Optional[Fraction] = None
    f_fphen_3_ets: Optional[Fraction] = None
    f_fphen_4_ets: Optional[Fraction] = None
    f_fphen_5_ets: Optional[Fraction] = None
    f_t_lem: Optional[Fraction] = None
    f_t_lma: Optional[Fraction] = None
    f_t_lep: Optional[Fraction] = None
    f_t_lse: Optional[Fraction] = None
    f_t_lse_mature: Optional[Fraction] = None  # fraction of mature leaf
    f_fphen_a: Optional[Fraction] = None
    f_fphen_b: Optional[Fraction] = None
    f_fphen_c: Optional[Fraction] = None
    f_fphen_d: Optional[Fraction] = None
    f_tt_emr: Optional[Fraction] = None
    #: Fraction of thermal time from sowing to first accumulation of ozone
    f_tt_fst_acc: Optional[Fraction] = None
    f_tt_veg: Optional[Fraction] = None
    f_tt_rep: Optional[Fraction] = None
    f_leaf_f_fphen: Optional[Fraction] = None

    # Vernalisation
    v_T_max: float = 30  #: Maximum temperature for vernalisation
    v_T_min: float = 15  #: Minimum temperature for vernalisation
    PIV: float = 1.5  #: Sensitivity to vernalisation

    # Latitude_function
    lat_f_k: Optional[float] = None  #: Latitude function parameter
    lat_f_b: Optional[float] = None  #: Latitude function parameter
    lat_f_c: Optional[float] = None  #: Latitude function parameter


@dataclass
class ModelConfig:
    """ Config options related to the model as a whole.

    Parameters
    ----------


    """
    flag_leaf_only: bool = False
    phenology_method: Optional[PhenologyMethods] = None
    dvi_method: DVIMethods = DVIMethods.DISABLED
    LAI_method: LAIMethods = LAIMethods.ESTIMATE_TOTAL
    time_type: TimeTypes = TimeTypes.THERMAL_TIME
    sgs_time_type: TimeTypes = TimeTypes.JULIAN_DAY
    sgs_key_day: KeyDays = KeyDays.SOWING_DAY
    zero_day: ZeroDayOptions = ZeroDayOptions.SOWING
    plant_emerge_method: PlantEmergeMethod = PlantEmergeMethod.CONSTANT
    flag_leaf_emerge_method: PlantEmergeMethod = PlantEmergeMethod.CONSTANT
    use_vernalisation: bool = False
    use_photoperiod_factor: bool = False
    sowing_day_method: SowingDateMethods = SowingDateMethods.INPUT
    latitude: Optional[float] = None  #: Only required for sowing_day_methods.LATITUDE


@dataclass
class PhenologyConfig:
    model: ModelConfig = field(default_factory=lambda: ModelConfig())
    species: List[SpeciesConfig] = \
        field(default_factory=lambda: fill_np_array_with_cls(
            DEFAULT_LC, SpeciesConfig)) # type: ignore

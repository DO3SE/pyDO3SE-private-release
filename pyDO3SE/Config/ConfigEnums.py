from enum import Enum
from deprecated import deprecated


class LandCoverType(Enum):
    CROP = "CROP"
    FOREST = "FOREST"


@deprecated(reason="Use senescence method instead")
class LeafFPhenAnetInfluence(Enum):
    """Methods by which we can influence Anet using leaf f phen.

    Options
    -------
    DISABLED = "disabled"
        No effect
    V_C_MAX = "v_c_max"
        Apply multiplicative leaf_f_phen to V_cmax_25 and J_max_25
    G_STO = "g_sto"
        Apply multiplicative leaf_f_phen to leaf_gsto

    """

    DISABLED = "disabled"
    V_C_MAX = "v_c_max"
    G_STO = "g_sto"


class LifeSpanMethods(Enum):
    """Life span method for calculating t_l etc.

    Options
    -------
    CONSTANT = "constant"
        Use config parameters
    JULES = "JULES"
        Use JULES phyllochron phenology method
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


class FVPDMethods(Enum):
    """f_VPD Methods.

    Options
    -------
    DISABLED = "disabled"
        f_VPD supplied (or left at default value of 1.0)
    LINEAR = "linear"
        Linear f_VPD relationship between VPD_max and VPD_min
    LOG = "log"
        Simple, unparameterised, logarithmic relationship
    LEUNING = "leuning"
        Calculate f_VPD using Leauning method inside Anet model.
    DANIELSSON = "danielsson"
        Calculate f_VPD using Danielsson et al.(2003) method inside Anet model.

    """

    DISABLED = "disabled"
    LINEAR = "linear"
    LOG = "log"
    LEUNING = "leuning"
    DANIELSSON = "danielsson"


class TLeafMethods(Enum):
    """Leaf temperature methods.

    Options
    -------
    INPUT = "input"
        Tleaf_C supplied
    AMBIENT = "ambient"
        Use ambient air temperature
    NIKOLOV = "Nikolov"
        NOT IMPLEMENTED! Use estimation method based on from Nikolov (1995)
    EB = "EB"
        NOT IMPLEMENTED! Use estimation method based on "An Introduction to
                     Environmental Biophysics" (Campbell & Norman, 1998)
    DE_BOECK = "de Boeck"
        Use estimation method based on "Leaf temperatures in
                     glasshouses and open-top chambers" (de Boeck, 2012)

    """

    INPUT = "input"
    AMBIENT = "ambient"
    NIKOLOV = "Nikolov"
    EB = "EB"
    DE_BOECK = "de Boeck"


class CanopyHeightMethods(Enum):
    """Canopy height methods.

    Options
    -------
    CONSTANT = "constant"
        canopy height provided as a constant value
    INPUT = "input"
        canopy height provided in external inputs
    CARBON = "carbon"
        canopy height calculated from carbon pools


    Additional Required Parameters
    ------------------------------


    """

    CONSTANT = "constant"
    INPUT = "input"
    CARBON = "carbon"


class LAIMethods(Enum):
    """Canopy height methods.

    Options
    -------
    CONSTANT = "constant"
        LAI provided as a constant value
    INPUT = "input"
        LAI provided in external inputs
    INPUT_HOURLY = "input hourly"
        LAI provided in external inputs hourly
    INPUT_TOTAL = "input total"
        total LAI provided and split between layers
    ESTIMATE_TOTAL = "estimate total"
        estimate according to land cover growing season (Using day_plf)
    DVI_LIMITED_CONSTANT = "dvi limited constant"
        Constant LAI but set to 0 outside DVI growing season
    LEAF_F_PHEN_LIMITED_CONSTANT = "leaf_f_phen limited constant"
        Constant LAI when leaf_f_phen > 0 else 0
    CARBON = "carbon"
        LAI calculated from carbon pools

    Additional Required Parameters
    ------------------------------

    constant
    ^^^^^^^^
        - config.config_land_cover.LAI
        - config.Land_Cover.fLAI[iL][iLC]

    input
    ^^^^^
        - external_data.LAI

    input hourly
    ^^^^^^^^^^^^

    If true then we run processes that rely on lai hourly instead of daily
    Note this is not fully tested and does not work with Penman Monteith calculations

        - external_data.LAI provided hourly

    input total
    ^^^^^^^^^^^
        - external_data.LAI
        - config.Land_Cover.fLAI[iL][iLC]

    estimate total
    ^^^^^^^^^^^^^^
        - config.Land_Cover.parameters[0].phenology.key_dates.harvest
        - config.Land_Cover.parameters[0].phenology.LAI_1
        - config.Land_Cover.parameters[0].phenology.LAI_2
        - config.Land_Cover.parameters[0].phenology.LAI_a
        - config.Land_Cover.parameters[0].phenology.LAI_b
        - config.Land_Cover.parameters[0].phenology.LAI_c
        - config.Land_Cover.parameters[0].phenology.LAI_1
        - config.Land_Cover.parameters[0].phenology.LAI_d
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
    INPUT_HOURLY = "input hourly"
    INPUT_TOTAL = "input total"
    ESTIMATE_TOTAL = "estimate total"
    DVI_LIMITED_CONSTANT = "dvi limited constant"
    LEAF_F_PHEN_LIMITED_CONSTANT = "leaf_f_phen limited constant"
    CARBON = "carbon"


class DVIMethods(Enum):
    """DVI methods

    DISABLED = "disabled"
       Not calculated
    JULES = "JULES"
       Use JULES thermal time method
    THERMAL TIME = "THERMAL TIME"
       Use THERMAL TIME PLF method
    INPUT = "INPUT"
       Use external input (NOT IMPLEMENTED)

    """

    DISABLED = "disabled"
    JULES = "JULES"
    THERMAL_TIME = "THERMAL_TIME"
    INPUT = "INPUT"


class LayerLAIDistributionMethods(Enum):
    """Methods of distributing LAI between layers.

    Options
    -------
    FRACTION = "fraction"
        Use a constant fraction per layer
    MAX_LAI_PER_LAYER = "max_lai_per_layer"
        Use a maximum value per layer

    """

    SKIP = "skip"
    FRACTION = "fraction"
    MAX_LAI_PER_LAYER = "max_lai_per_layer"


class SenescenceFunctionMethods(Enum):
    """Methods for defining leaf senesence.

    Only valid for Anet and Ewert photosynthesis models.

    Options
    -------
    DISABLED = "disabled"
        f_LS is always 1
    ANET = "anet"
        f_LS is set to equal leaf f phen and fO3_d is set to fO3 Multiplicative
    EWERT = "ewert"
        f_LS is calculated from f_LS phenology and fO3_d is calculated from fst

    """

    DISABLED = "disabled"
    ANET = "anet"
    EWERT = "ewert"


class OzoneDepositionMethods(Enum):
    """Methods for calculating ozone at top of canopy.

    Options
    -------
    SINGLE_LAYER = "single layer"
        Use big leaf method
    MULTI_LAYER = "multi layer"
        Use multilayer method

    """

    SINGLE_LAYER = "single layer"
    MULTI_LAYER = "multi layer"


class GstoMethods(Enum):
    MULTIPLICATIVE = "multiplicative"
    PHOTOSYNTHESIS = "photosynthesis"

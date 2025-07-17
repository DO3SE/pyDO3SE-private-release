from dataclasses import dataclass, field
from typing import List
import numpy as np

from data_helpers.fill_np_array import fill_np_array_with_cls
from do3se_met.soil_moisture.enums import FSW_Methods

from do3se_phenology.config import ModelConfig as PhenologyModelConfig
from do3se_phenology.config import SpeciesConfig as PhenologySpeciesConfig

from pyDO3SE.constants.model_constants import DEFAULT_LAYERS, DEFAULT_LC
from pyDO3SE.plugins.gsto.ewert.enums import EwertLoopMethods, AdjustNegativeAnMethods

from .ConfigEnums import (
    GstoMethods,
    TLeafMethods,
    FVPDMethods,
    CanopyHeightMethods,
    DVIMethods,
    LAIMethods,
    LayerLAIDistributionMethods,
    LandCoverType,
    EnabledOrDisabled,
    FO3_methods,
    FTempMethods,
    SenescenceFunctionMethods,
    LeafFPhenAnetInfluence,
)


@dataclass(frozen=False)
class Config_Gsto:
    """Configuration used in all gsto methods."""
    #: Stomatal conductance method:
    #:   - "multiplicative":   Use DO3SE multiplicative model
    #:   - "photosynthesis":   Use Farquar-based photosynthesis model (hybrid
    #:                         with some multiplicative components)
    method: GstoMethods = GstoMethods.MULTIPLICATIVE

    #: Tleaf method:
    Tleaf_method: TLeafMethods = TLeafMethods.AMBIENT

    fmin: float = None    #: Minimum stomatal conductance [fraction]

    f_VPD_method: FVPDMethods = FVPDMethods.DISABLED

    #: Critical daily VPD threshold above which stomatal conductance will stop
    #: increasing [kPa].
    VPD_crit: float = 1000.0

    VPD_min: float = None     #: VPD for minimum gsto [kPa]
    VPD_max: float = None     #: VPD for maximum gsto [kPa]

    #: == SOIL MOISTURE
    #: f_SW method:
    #:   - "disabled":     f_SW supplied (or left at default value of 1.0)
    #:   - "fSWP exp":     Use fSWP exponential curve (see fSWP_exp_curve)
    #:   - "fSWP linear":  Use linear fSWP function (see SWP_min and SWP_max)
    #:   - "fLWP exp":     Use fSWP exponential curve, but with LWP instead of SWP
    #:   - "fPAW":         Use fPAW relationship
    #: This was previously fxwp_method
    f_SW_method: FSW_Methods = FSW_Methods.DISABLED

    #: fSWP linear parameters:
    SWP_min: float = None   #: SWP for minimum gsto [MPa]
    SWP_max: float = None   #: SWP for maximum gsto [MPa]

    #: fSWP exponential curve:
    #:   - "custom":         fSWP_exp_a and fSWP_exp_b supplied
    #:   - "temperate":      a = 0.355, b = -0.706
    #:   - "mediterranean":  a = 0.619, b = -1.024
    fSWP_exp_curve: str = "temperate"
    fSWP_exp_a: float = None
    fSWP_exp_b: float = None

    # > ASW (available soil water) for minimum gsto (percent of ASW at field capacity)
    ASW_min: float = 0.0
    # > ASW (available soil water) for maximum gsto (percent of ASW at field capacity)
    ASW_max: float = 50.0


@dataclass(frozen=False)
class Config_Multip_Gsto:
    """Configuration for multiplicative stomatal conductance functions."""

    #: f_light method:
    #:   - "disabled":   f_light and leaf_f_light supplied (or left at default
    #:                   value of 1.0)
    #:   - "enabled":    f_light and leaf_f_light calculated
    f_light_method: EnabledOrDisabled = EnabledOrDisabled.DISABLED
    f_lightfac: float = 0.006  #: Single leaf f_light coefficient

    #: f_temp method:
    #:   - "disabled":       f_temp supplied (or left at default value of 1.0)
    #:   - "default":        Normal bell-shaped function over T_min -> T_opt -> T_max
    #:   - "square high":    Same as "default", but straight lines from
    #:                       (T_opt, 1.0) -> (T_max, 1.0) -> (T_max, 0.0)
    f_temp_method: FTempMethods = FTempMethods.DISABLED
    T_min: float = None     #: Minimum temperature [degrees C]
    T_opt: float = None     #: Optimum temperature, for max. gsto [degrees C]
    T_max: float = None     #: Maximum temperature [degrees C]

    gmax: float = None    #: Maximum stomatal conductance [mmol O3 m-2 PLA s-1]
    gmorph: float = 1.0    #: Sun/shade morphology factor [fraction]

    #: f_O3 method:
    #:   - "disabled":   f_O3 supplied (or left at default value of 1.0)
    #:   - "wheat":      Wheat f_O3 method
    #:   - "potato":     Potato f_O3 method
    f_O3_method: FO3_methods = FO3_methods.DISABLED


@dataclass(frozen=False)
class Config_PnGsto:
    """Parameters for photosynthesis-based stomatal conductance."""

    g_sto_0: float = None       #: Closed stomata conductance [umol m-2 s-1]
    m: float = None             #: Species-specific sensitivity to An [dimensionless]
    use_O3_damage: bool = True  #: use fO3 damage factors(fO3_d and fO3_l) in ewert gsto model.

    #: Method for defining f_LS and fO3d
    senescence_method: SenescenceFunctionMethods = SenescenceFunctionMethods.DISABLED

    #: Ewert loop method
    #:
    #: "iterative"
    #: "cubic"
    ewert_loop_method: EwertLoopMethods = EwertLoopMethods.CUBIC

    #: V/J max method:
    #:   - "input":      V_cmax_25 and J_max_25 supplied as hourly inputs
    #:   - "constant":   Use constant V_cmax_25 and J_max_25 values (below)
    V_J_method: str = "constant"

    V_cmax_25: float = None     #: Maximum catalytic rate at 25 degrees [umol m-2 s-1]
    J_max_25: float = None      #: Maximum rate of electron transport at 25 degrees [umol m-2 s-1]

    V_cmax_25_kN: float = 0.78     #: Nitrogen profile coefficient [-]

    R_d_coeff: float = None    #: Dark respiration coefficient [Fraction] (Clark et al 2011)
    r_g: float = 0.25  # : fraction of gross primary productivity less maintaenance respiration that
    #:  is assigned to growth respiration [fraction]

    #: D_0 method:
    #:   - "constant":   Use constant D_0 value
    #:   - "f_VPD":      Determine D_0 from multiplicative f_VPD method (at f_VPD = 0.5)
    D_0_method: str = "f_VPD"
    #: "The VPD at which g_sto is reduced by a factor of 2" [kPa] (Leuning et al. 1998)
    D_0: float = None

    #: Ewert Tleaf iteration parameters
    #: Threshold (from 0) to consider co2 concentration equation as "balanced"
    co2_concentration_balance_threshold: float = 0.001
    #: Maximum number of iterations to find co2 concentration solution
    co2_concentration_max_iterations: int = 50

    #:     If True then allow negative A_n values, else return NaN
    #:     If "last_resort" then allow negative A_n values if no other solution is found
    #:     If "clip" then clip negative A_n values to -R_d
    adjust_negative_A_n: AdjustNegativeAnMethods = AdjustNegativeAnMethods.FALSE


    # #: TODO: This is not currently implemented
    # #: Threshold (from 0) to consider leaf energy balance equation as "balanced"
    # Tleaf_balance_threshold: float = 0.001
    # #: Maximum number of iterations to find Tleaf solution
    # Tleaf_max_iterations: int = 50
    # #: Factor to apply to energy balance when adjusting leaf temperature
    # Tleaf_adjustment_factor: float = 0.02

    #: Parameters for martin2000 O3 effect method
    #: DEPRECATED
    K_z: float = 24.0      #: Coefficient for ozone damage [dimensionless]
    F_0: float = 1.0       #: Threshold flux for ozone damage [nmol O3 m-2 PLA s-1]

    #: Use leaf f phen to limit anet:
    #: @deprecate(reason="Use senescence method instead")
    leaf_f_phen_Anet_influence: LeafFPhenAnetInfluence = LeafFPhenAnetInfluence.DISABLED  #: formally phenology_method

    #: temperature sensitivity parameters for phenology (Osborne et al 2015)
    t_b: float = 0  #: base temperature [degC]
    t_o: float = 20  #: optimum temperature [degC]
    t_m: float = 30  #: maximum temperature [degC]

    #: Photoperiod parameters
    p_crit: float = 24  #: critical photoperiod
    p_sens: float = 0  #: sensitivity of development rate to photoperiod

    #: Thermal time interval parameters
    tt_emr: float = None  #: thermal time interval between sowing and emergence [deg days]
    tt_veg: float = None  #: thermal time interval between emergence and flowering [deg days]
    tt_rep: float = None  #: thermal time interval between flowering and maturity/harvest [deg days]

    #: Coefficients from Ewert Paper for impact of ozone on photosynthesis
    #: "low ozone concentrations are detoxi®ed without direct effects on the photosynthetic system"
    gamma_1: float = 0.06
    #: "describes the decrease in A c perunit of ozone influx" [nmol m-2 s-1] -1
    gamma_2: float = 0.0045
    gamma_3: float = 0.0005  # : multiplier for ozone impact fO3_l [umol m-2 O3]
    #: ozone impact on start of senesence (See Ewert for docs)[0-1]
    gamma_4_senes: float = 1
    #: ozone impact on end of senesense (See Ewert for docs)[0-1]
    gamma_5_harvest: float = 0

    #: critical cumulative stomatal ozone flux [umol m-2 O3]. I.e the point at which
    #: ozone accumulation affects the onset of senescence
    cL3: float = 20000


@dataclass(frozen=False)
class Config_Land_Cover_Parameters:
    """Land cover properties. Formally LandCover_t."""

    name: str = ""

    phenology: PhenologySpeciesConfig = field(default_factory=lambda: PhenologySpeciesConfig())
    gsto: Config_Gsto = field(default_factory=lambda: Config_Gsto())
    multip_gsto: Config_Multip_Gsto = field(default_factory=lambda: Config_Multip_Gsto())
    pn_gsto: Config_PnGsto = field(default_factory=lambda: Config_PnGsto())

    height: float = None  #: Constant canopy height[m] # TODO: Is this copy of height in
    cosA: float = 0.5          #: cos(A), A = mean leaf inclination (0.5 = 60 degrees)
    Lm: float = None      #: Leaf dimension[m]
    # TODO: Check ok to set Y to default 0
    Y: float = 0       #: POD threshold [nmol O3 m-2 s-1]
    PID: float = 40  # : photoperiod factor parameter


@dataclass(frozen=False)
class Config_Land_Cover:
    """Land cover configuration. Formally LandCoverConfig_t."""

    nL: int = 1           #: Number of layers
    nLC: int = 1          #: Number of land covers configured
    nP: int = 1           #: Number of leaf populations
    primary_LC: int = 0   #: Primary land cover for height, LAI, etc.
    land_cover_type: LandCoverType = LandCoverType.CROP

    # TODO: Check land cover props
    #: Land cover parameters
    parameters: List[Config_Land_Cover_Parameters] = \
        field(default_factory=lambda: fill_np_array_with_cls(
            DEFAULT_LC, Config_Land_Cover_Parameters))

    #: Height method:
    height_method: CanopyHeightMethods = CanopyHeightMethods.CONSTANT

    #: DVI method:
    dvi_method: DVIMethods = DVIMethods.DISABLED

    #: Distribution of layers heights within canopy; height of top of layer as
    #: a fraction of the canopy height.
    layer_height_frac: List[float] = field(default_factory=lambda: np.full((DEFAULT_LAYERS), 1.0))

    #: LAI method
    LAI_method: LAIMethods = LAIMethods.ESTIMATE_TOTAL

    LAI: float = None  #: Used if LAI_method is constant

    #: SAI method:
    #:   - "input":          Input all SAI values
    #:   - "input total":    Input total SAI, split according to fLAI
    #:   - "estimate total": Estimate total SAI from primary land cover and
    #:                       total LAI, split according to fLAI
    SAI_method: str = "estimate total"
    SAI: float = None  #: Used if LAI_method is constant

    LAI_distribution_method: LayerLAIDistributionMethods = LayerLAIDistributionMethods.FRACTION

    #: Distribution of LAI/SAI in canopy (default: uniform distribution)
    #: Only used for LayerLAIDistributionMethods.FRACTION
    #: REAL, dimension(MAX_LAYERS,MAX_LC) :: fLAI = 1
    fLAI: List[List[float]] = field(
        default_factory=lambda: np.full((DEFAULT_LAYERS, DEFAULT_LC), 1.0))

    max_lai_per_layer: float = None  #: Max leaf area index per layer
    leaf_emergence_multiplier: float = 1  #: rate of leaf emergence multiplier

    #: Contains global settings related to the processing of phenology
    phenology_options: PhenologyModelConfig = field(
        default_factory=lambda: PhenologyModelConfig(),
    )


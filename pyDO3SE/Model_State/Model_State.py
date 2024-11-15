"""Model state is the state within the mode that can be effected by model processes

It consists of named tuples that can be accessed via dot notation
i.e.
canopy_layer_component[0][0].leaf_gsto

State can be split into components, layers and leaf populations.

#: TODO: We should reorder to component, layer, population
Array shapes should be in this order. (layer, component, population)

"""
import numpy as np
from dataclasses import dataclass, field
from typing import List

from data_helpers.fill_np_array import fill_np_array_with_cls
from do3se_phenology.state import LeafPhenologyState, PhenologyState
from do3se_met.resistance.model import Leaf_Resistance_Model, Resistance_Model
from pyDO3SE.plugins.soil_moisture.state import PM_State_t, SMDData_t
from pyDO3SE import settings


@dataclass(frozen=False)
class Temporal_State:
    """Current time state."""
    dd: int = None  #: Days from day 0. Usually day of year (0 -> 364)
    hr: int = None  #: Hour of day (0 -> 23)
    row_index: int = None  #: id of current data row (Uses helpers.get_row)
    td: float = None  #: Thermal time
    dd_offset: int = 0  # : Number of days to offset dd value. Useful for multi year runs


@dataclass(frozen=False)
class External_Meteorological_State:
    """Meteorological state external to the model.

    Data brought here from external state shape
    """

    photoperiod: float = 24  #: The day length (Weir et al 1984)
    t_acc: float = 0  #: Accumulated daily temperature


@dataclass(frozen=False)
class Internal_Meteorological_State:
    """Meteorological state within the model.

    variable data from MetData_t in DO3SE fortran model
    """

    #: TODO: Can these be moved to pre model setup?
    u_i: float = None  #: Decoupled windspeed at 50m [m/s]
    ustar: float = None  #: Friction velocity over target canopy [m/s]
    ustar_ref: float = None  #: Friction velocity over reference canopy [m/s]
    O3_i: float = None  #: Decoupled O3 concentration at 50m [ppb]
    L: float = 9999999  #: Monin-Obukhov length [m]


@dataclass(frozen=False)
class layer_Meteorological_State:
    """Meteorological state at each layer of the model.

    "Local" meteorological data, i.e. dependent on the canopy.

    Represents the value at the top of the canopy - should be used as an array
    for multiple layers, representing the value at the top of each layer.

    MicroMetData_t in DO3SE fortran model
    """

    PARsun: float = None    #: PAR received by sunlit leaves [W m-2]
    PARshade: float = None  #: PAR received by shaded leaves [W m-2]

    #: umet
    micro_u: float = None         #: Windspeed at layer canopy [m/s]
    micro_O3_top: float = None        #: O3 concentration at layer canopy top [ppb]
    micro_O3: float = None        #: O3 concentration at layer canopy middle [ppb]
    micro_O3_nmol_m3: float = None        #: O3 concentration at layer canopy [nmol_m3]

    #: TODO: Remove this as now in Canopy_Layer_Component
    #: Tleaf_C: float = None   #: Leaf temperature [degrees C]


@dataclass(frozen=False)
class Multip_Gsto_params:
    """Multiplicative stomatal conductance parameters.

    ADDITIONAL PARAMS BROUGHT FROM Gsto_ml.F90
    #: TODO: Check for any duplication or redundancy
    GstoParams_t
    """

    #: Duplicate of config params
    #: fmin: float = None         #: Minimum stomatal conductance [fraction]
    #: gmax: float = None         #: Maximum stomatal conductance [mmol O3 m-2 PLA s-1]
    #: gmorph: float = 1.0        #: Sun/shade morphology factor [fraction]

    f_phen: float = 1.0        #: Phenology-related effect on gsto [fraction]
    leaf_f_phen: float = 1.0   #: Phenology-related effect on leaf gsto [fraction]
    f_light: float = 1.0       #: Irradiance effect on gsto [fraction]
    leaf_f_light: float = 1.0  #: Irradiance effect on leaf gsto [fraction]
    f_temp: float = 1.0        #: Temperature effect on gsto [fraction]
    f_VPD: float = 1.0         #: VPD effect on gsto [fraction]
    f_SW: float = 1.0          #: Soil water effect on gsto [fraction]
    f_O3: float = 1.0          #: O3 effect on gsto [fraction]


@dataclass(frozen=False)
class Whole_Canopy_State:
    """Variables that do not vary throughout the canopy. Formally V_t"""
    #: Phenology and Canopy Structure
    canopy_height: float = None           #: Overall canopy height
    SAI_total: float = None  #: sum of layer component SAI
    LAI_total: float = None  #: sum of layer component LAI
    LAI_brown_total: float = None  #: sum of layer component LAI_brown
    Lm_LAI: float = None  #: LAI weighted mean leaf width

    #: Soil moisture data
    SMD: SMDData_t = field(default_factory=lambda: SMDData_t())    #: Results of SMD calculations)
    #: Penman-Monteith results and accumulators
    PM: PM_State_t = field(default_factory=lambda: PM_State_t())

    #: TODO: RModels are storing data for multiple layers
    #:       Should we split the model for 1 per layer??

    #: Multi-layer canopy O3 resistance model for modelled canopy
    rmodel_O3: Resistance_Model = field(default_factory=lambda: Resistance_Model(
        settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS))  #: Must be initialised
    #: Multi-layer canopy O3 resistance model reference canopy
    rmodel_O3_ref: Resistance_Model = field(
        default_factory=lambda: Resistance_Model(1))  #: Must be initialised
    Vd: float = None  #: Velocity of O3 deposition to top of canopy [m/s]
    canopy_top_o3: float = None  #: O3 concentration at top of canopy [ppb]

    #: Single-layer canopy water vapour resistance model
    rmodel_H2O: Resistance_Model = field(default_factory=lambda: Resistance_Model(
        settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS))   #: Must be initialised
    Es_blocked: bool = False  #: Should soil evaporation be blocked?

    #: met state at top of canopy
    met: layer_Meteorological_State = field(default_factory=lambda: layer_Meteorological_State())

    d: float = None  #: Canopy displacement height [m]
    z0: float = None  #: Canopy roughness length [m]


@dataclass
class Canopy_Population_State:
    """Variables that vary by canopy leaf population."""

    #: Below are upscaled to total leaf values and multiplied by actual LAI
    A_n: float = None           #: Net CO2 assimilation rate [umol CO2 m-2 PLA s-1]
    #: Net CO2 assimilation rate at top sunlit part of canopy [umol CO2 m-2 PLA s-1]
    A_n_sunlit: float = None
    A_n_limit_factor: float = None           #: Net CO2 assimilation rate [umol CO2 m-2 PLA s-1]
    A_c: float = None  #: Rubisco activity limited rate of photosynthesis [umol m-2 s-1 CO2]
    A_j: float = None  #: limited assimilation rate (A_q in paper) [umol m-2 s-1 CO2]
    A_p: float = None  #: Triose phosphate utilisation limited assimilation rate [umol m-2 s-1 CO2]
    R_d: float = None  #: day respiration rate [micro mol/(m^2*s) CO2]
    c_i: float = None  #: CO2 concentration inside stomata   [micromol/mol]

    #: g_sv and leaf_gsto are average per m^2 (Not upscaled)
    g_sv_per_layer: List[float] = None  #: gsto for CO2 from Ewert model
    g_sv_sunlit: float = None  # : gsto for sunlit part of top layer

    #: mean Leaf stomatal conductance[mmol O3 m-2 PLA s-1]
    mean_gsto_per_layer: List[float] = field(default_factory=lambda: [None for _ in range(
        settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS)])        #: Leaf stomatal conductance [mmol O3 m-2 PLA s-1]
    #: mean Leaf stomatal conductance for sunlit top layer[mmol O3 m-2 PLA s-1]
    mean_gsto_sunlit: float = None
    #: mean Leaf stomatal conductance * LAI [mmol O3 PLA s-1]
    bulk_gsto_per_layer: List[float] = field(
        default_factory=lambda: [None for _ in range(settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS)])

    #: TODO: This is from gsto_params
    f_VPD: float = 1.0         #: VPD effect on gsto [fraction]

    #: Photosynthetic stomatal conductance parameters/results
    #: Moved to Canopy_Component_population
    V_cmax_25_per_layer: List[float] = field(default_factory=lambda: [None for _ in range(
        settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS)])  #: Maximum catalytic rate at 25 degrees [umol m-2 s-1]
    J_max_25: float = None      #: Maximum rate of electron transport at 25 degrees [umol m-2 s-1]

    V_cmax: float = None     #: catalytic rate [umol m-2 s-1]
    J_max: float = None      #: rate of electron transport [umol m-2 s-1]

    O3up: float = 0  #: Stomatal ozone flux (fst) [nmol O3 m-2 PLA s-1]
    O3up_sunlit: float = 0  #: Stomatal ozone flux (fst) [nmol O3 m-2 PLA s-1]
    #: Accumulated Stomatal ozone flux for current day [nmol O3 m-2 PLA s-1]
    O3up_acc_day: float = 0.0
    O3up_acc_day_prev: float = 0.0

    O3up_acc: float = 0.0  #: Accumulated Stomatal ozone flux [nmol O3 m-2 PLA h-1]

    POD_0: float = 0.0          #: Phytotoxic Ozone Dose, no threshold (mmol m-2 PLA)
    #: Phytotoxic Ozone Dose, no threshold for sunlit leaves (mmol m-2 PLA)
    POD_0_sunlit: float = 0.0
    POD_Y: float = 0.0          #: Phytotoxic Ozone Dose above threshold Y (mmol m-2 PLA)
    #: Phytotoxic Ozone Dose above threshold Y for sunlit leaves (mmol m-2 PLA)
    POD_Y_sunlit: float = 0.0
    OT_0: float = None         #: [ppm]
    OT_40: float = None        #: [ppm]
    AOT_0: float = 0.0          #: [ppm h]
    AOT_40: float = 0.0         #: [ppm h]

    #: Ozone damage factors
    #:   Hourly ozone impace factor [dimensionless][0-1]
    fO3_h: float = 1
    #:   Cumulative ozone effect from previous hour [dimensionless][0-1]
    fO3_d: float = 1
    #: Ozone exposure factor [umol/m^2 O3]
    fO3_l: float = 1
    #: fraction used to simulate the recovery from ozone damage as a factor of leaf age [Dimensionless][0-1]
    f_LA: float = 1
    #: factor effect of leaf senescence on A_c [Dimensionless][0-1]
    f_LS: float = 1
    #: Leaf capacity to recover from Ozone damage [dimensionless][0-1]
    rO3: float = 1

    #: limited life span lengths (ozone) [Thermal Time]
    t_lep_limited: float = None
    t_lse_limited: float = None

    #: Leaf-level O3 resistance model
    leaf_rmodel_O3: Leaf_Resistance_Model = field(default_factory=lambda: Leaf_Resistance_Model())
    leaf_sunlit_rmodel_O3: Leaf_Resistance_Model = field(
        default_factory=lambda: Leaf_Resistance_Model())

    #: Phenology
    phenology: LeafPhenologyState = field(default_factory=lambda: LeafPhenologyState())

    #: Leaf population total LAI
    LAI: float = 0
    #: Fraction of each layer that is this leaf population (Sums up to 1.0)
    fLAI_layer: List[float] = field(default_factory=lambda: np.zeros(
        settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS))


@dataclass(frozen=False)
class Canopy_Layer_State:
    """Variables that vary by canopy layer. Formally ML_t."""

    layer_height: float = None            #: Height of top of layer[m]
    layer_depth: float = None            #: depth of layer[m]

    #: Meteorological data
    micro_met: layer_Meteorological_State = field(default_factory=lambda:
                                                  layer_Meteorological_State())


@dataclass(frozen=False)
class Canopy_Component_State:
    """Variables that vary by canopy component (land cover). Formally MC_t."""

    LC_dist: float = None  #: LAI land cover distribution

    #: LAI in each leaf population per layer (nL, nP)
    leaf_pop_distribution: List[List[float]] = field(default_factory=lambda: np.zeros(
        (settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS, settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS)))

    td: float = None  #: Thermal time

    t_eff: float = None  #: Effective Thermal Temperature [degC] (Osborne 2015)
    t_eff_acc: float = None  #: accumulated effective temperature [degC]

    #: Vernalised Td
    td_v: float = 0  #: thermal time since plant emergence
    #: Age of plant(deg C day). difference in td between current and emergence
    td_dd: float = 0  #: thermal time since plant emergence

    #: Plant and Leaf emergence
    phenology: PhenologyState = field(default_factory=lambda: PhenologyState())
    #: has_emerged: bool = False  #: Plant
    flag_has_emerged: bool = False  #: Flag leaf
    emergence_rate: float = 0.0  #: leaf emergence rate
    total_emerged_leaf_populations: int = 0
    growing_populations: List[bool] = field(default_factory=lambda: np.zeros(
        settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS))

    #: Thermal time difference between leaf pop emergence and current td
    td_dd_leaf_pops: List[float] = field(default_factory=lambda: np.zeros(
        settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS))

    #: Plant life stages in thermal time for the Ewert model
    #: Full life span of leaf in thermal time as estimated in phenology module [thermal time]
    t_l_estimate: float = None
    #: time from seed to end of emerging leaf/start of mature leaf [thermal time]
    t_lem: float = None
    t_lma: float = None  #: Full time of mature leaf (t_lep + t_lse) [thermal time]
    t_lep: float = None  #: Time during which leaf is expanding [thermal time]
    t_lse: float = None  #: Time that the leaf is senescing [thermal time]

    rpe: float = None  #: relative photoperiod effect

    dvi: float = -1.001  #: Development Index [Dimensionless] (Osborne 2015)
    phyllochron: float = None  #: Phyllochron (McMaster and Wilhelm 1995)

    #: Canopy gsto
    canopy_gsto: float = 0  #: Canopy total gsto
    canopy_anet: float = 0  #: Canopy total net carbon assimilation
    canopy_an_gross: float = 0  #: Canopy total gross carbon assimilation

    #: D_0 "The VPD at which g_sto is reduced by a factor of 2" [kPa] (Leuning et al. 1998)
    D_0: float = None
    g_bv_H2O: float = None  #: boundary layer conductance for forced convection [umol m-2 s-1 H2O]

    #: Carbon state ==
    #: Carbon Pools
    c_root: float = 0
    c_stem: float = 0
    c_leaf: float = 0.00001
    c_lbrn: float = 0.00001
    c_harv: float = 0
    c_resv: float = 0
    #: Carbon distribution fractions
    p_root: float = 0
    p_stem: float = 0
    p_leaf: float = 0
    p_harv: float = 0
    #: gross carbon assimilated per hour  [kg C m^-2 hour^-1]
    gpp: float = 0
    #: Net carbon assimilated per hour  [kg C m^-2 hour^-1]
    npp: float = 0
    #: net carbon assimilated aggregated to end of day [kg C m^-2 day^-1]
    npp_acc: float = 0

    R_pg: float = 0  # : maintanance respiration [kg C m^-2 hour^-1]
    R_pm: float = 0  # : growth respiration [kg C m^-2 hour^-1]

    #: Yield
    #: stem_dm [?]
    stem_dm: float = 0
    #: leaf_dm [?]
    leaf_dm: float = 0
    #: lbrn_dm [?]
    lbrn_dm: float = 0
    #: total_leaf_dm [?]
    total_leaf_dm: float = 0
    #: straw_dm [?]
    straw_dm: float = 0
    #: ear_dm [?]
    ear_dm: float = 0
    #: aboveground_dm [?]
    aboveground_dm: float = 0
    #: belowground_dm [?]
    belowground_dm: float = 0
    #: grain_dm [?]
    grain_dm: float = 0
    #: harvest_index [?]
    harvest_index: float = 0
    #: yield_ha [?]
    yield_ha: float = 0

    plant_height: float = 0  #: Specific to component
    R_dc: float = 0  #: canopy leaf dark respiration  [umol/(m^2*s)

    #: component total LAI and SAI
    LAI: float = 0                     #: Leaf Area Index [m2 m-2]
    LAI_brown: float = 0                     #: brown Leaf Area Index [m2 m-2]
    SAI: float = 0                     #: Stand Area Index [m2 m-2]

    photoperiod_factor: float = 1


@dataclass(frozen=False)
class Canopy_Layer_Component_Pop:
    """Variables that vary by canopy layer, component and leaf pop.

    #: TODO: Move other variables here.
    """
    Tleaf_C_estimate: float = 0  #: estimated leaf temperature from photosynthesis


@dataclass(frozen=False)
class Canopy_Layer_Component_State:
    """Variables that vary by canopy layer and component. Formally MLMC_t."""

    LAI: float = None                     #: Leaf area index [m2 m-2]
    LAI_brown: float = None                     #: Leaf area index [m2 m-2]
    SAI: float = None                     #: Stand area index [m2 m-2]
    LAIsunfrac: float = None              #: Fraction of LAI that is sunlit

    #:   Multiplicative stomatal conductance parameters
    #: TODO: Move these to canopy component and population
    gsto_params: Multip_Gsto_params = field(default_factory=lambda: Multip_Gsto_params())

    FO3_eff: float = 0.0           #: (Accumulated) effective ozone dose [nmol O3 m-2 PLA]

    #: Gsto values are stored centrally as O3 then converted where needed
    mean_gsto: float = None        #: layer mean stomatal conductance [mmol O3 m-2 PLA s-1]
    bulk_gsto: float = None        #: layer total stomatal conductance [mmol O3 m-2 PLA s-1]

    #: Fraction of each population that makes up the layer
    fLAI_layer: List[float] = field(default_factory=lambda: np.zeros(
        settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS))


@dataclass(frozen=False)
class DebugState:
    """Debug state for storing intermediate values for debugging."""

    ewert_loop_iterations: int = 0  #: Max number of iterations in the Ewert loop for this hour


@dataclass(frozen=False)
class Model_State_Shape:
    """The data shape for the internal model state.

    Notes:
    V_t(Variables that do not vary throughout the canopy)
    split into temporal, met and canopy data
    """

    #: Formally part of V_t
    temporal: Temporal_State = field(default_factory=lambda: Temporal_State())

    external_met: External_Meteorological_State = field(
        default_factory=lambda: External_Meteorological_State())

    met: Internal_Meteorological_State = field(
        default_factory=lambda: Internal_Meteorological_State())

    canopy: Whole_Canopy_State = field(default_factory=lambda: Whole_Canopy_State())

    canopy_component_population: List[List[Canopy_Population_State]] = \
        field(default_factory=lambda:
              fill_np_array_with_cls(shape=(settings.global_settings.MAX_NUM_OF_CANOPY_COMPONENTS, settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS),
                                     cls=Canopy_Population_State))

    #: Formally ML_t
    canopy_layers: List[Canopy_Layer_State] = \
        field(default_factory=lambda:
              fill_np_array_with_cls(settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS, Canopy_Layer_State))


    #: Formally MC_t
    canopy_component: List[Canopy_Component_State] = \
        field(default_factory=lambda:
              fill_np_array_with_cls(settings.global_settings.MAX_NUM_OF_CANOPY_COMPONENTS, Canopy_Component_State))

    #: Formally MLMC_t
    canopy_layer_component: List[List[Canopy_Layer_Component_State]] = \
        field(default_factory=lambda: fill_np_array_with_cls(
            shape=(settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS,
                   settings.global_settings.MAX_NUM_OF_CANOPY_COMPONENTS),
            cls=Canopy_Layer_Component_State,
        ))

    canopy_layer_component_pop: List[List[List[Canopy_Layer_Component_Pop]]] = \
        field(default_factory=lambda: fill_np_array_with_cls(
            shape=(settings.global_settings.MAX_NUM_OF_CANOPY_LAYERS,
                   settings.global_settings.MAX_NUM_OF_CANOPY_COMPONENTS, settings.global_settings.MAX_NUM_OF_LEAF_POPULATIONS),
            cls=Canopy_Layer_Component_Pop,
        ))

    prev_hour: 'Model_State_Shape' = None  #: Stores the state from the previous hour

    debug: DebugState = field(default_factory=lambda: DebugState())

    #: def init():
    #:     #: TODO: Validate input
    #:     #: TODO: Initialize PM & SMD? below is from Fortran model
    #:     #: if (this%SMD_conf%source == "P-M") then
    #:     #:     this%V%PM = PM_state_t()
    #:     #:     this%V%SMD = soil_moisture_from_SWC(this%SMD_conf, this%SMD_conf%initial_SWC)
    #:     pass

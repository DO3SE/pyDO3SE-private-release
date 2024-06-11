""" External state config

The external state config contains configuration for state outside the plant.
These values should be calculated prior to the model run (Unless specifically stated).


Setting thermal time
--------------------
There are multiple ways to set thermal time.
- Input data
    To set thermal time to come from the external data file set the following
    - Set `td_method` to `input`
    - Set `thermal_time_method` to `ThermalTimeMethods.EXTERNAL`
- Calculated from first row of data
    To set the thermal time to be calculated from the thermal time data set the following
    - Set `td_method` to `calculated`
    - Set `td_base_temperature`
    - Set `thermal_time_method` to `ThermalTimeMethods.EXTERNAL`
- Calculate with custom start date and offset
    You can add an initial offset to the thermal time calculation using:
    - As Calculated from first row of data +
    - Set `thermal_time_start` to the `dd` value to start thermal time accumulation
    - Set `thermal_time_offset` to the `td` value to set at `thermal_time_start`
- Calculate Hourly
    When doing hourly runs we only have an hour of data so cannot pre-calculate thermal time.
    Thermal time needs to be calculated in the model state instead of the external state.
    To do this set:
    - Set `td_method` to `skip`
    - Set `thermal_time_method` to `ThermalTimeMethods.HOURLY`

"""
from enum import Enum
from dataclasses import dataclass, field
from typing import List, NamedTuple

from do3se_phenology.units import JulianDay, ThermalTime


class FileTypes(Enum):
    CSV = "csv"
    JSON = "json"
    NETCDF = "nc"


class EStateOptions(NamedTuple):
    """Options for loading external state form input data file.

    Parameters
    ----------
    file_type : FileTypes
        The input data file type. One of [csv, netcdf]
    has_header_row : bool, optional
        csv has header row, by default True
    row_indexes : List[int], optional
        [description], by default None
    variable_map : dict, optional
        Map of variables in input data file to DO3SE input variable names
    multi_file_data: bool, optional
        If true then NetCDF is split with file per variable
    preprocess_map: dict, optional
        A map of variable to preprocess function
    zero_year: int, optional
        DOY (where year == zero_year + 1) == DOY + 365
    data_filter: str, optional
        When using multi file mode can use filename filter.
    file_suffix: str, optional
        Allow overrideing input file suffix. This is after chosen filetype
    netcdf_loader_kwargs: dict, optional
        kwargs to pass to xr.open_dataset or xr.open_mfdataset

    """
    file_type: FileTypes = FileTypes.CSV
    has_header_row: bool = True
    row_indexes: List[int] = None
    variable_map: dict = None
    multi_file_data: bool = False
    preprocess_map: dict = {}
    zero_year: int = None
    data_filter: str = ''
    file_suffix: str = None
    netcdf_loader_kwargs: dict = {}
    allow_sparse_inputs: bool = False


class InputMethod(Enum):
    CONSTANT = "constant"
    INPUT = "input"
    OFFSET = "offset"
    MULTIPLY = "multiply"
    CALCULATED = "calculated"
    SKIP = "skip"


class ParSunShadeMethods(Enum):
    """Methods of calculating PARsun shade from indirect and diffuse light data.

    Options
    -------
    SIMPLE = "SIMPLE"
        Use the simple
    FARQUHAR1997 = "FARQUHAR1997"
        Use the Farquhar 1997 method
    DO3SEUI = "DO3SEUI"
        Use the method used in the DO3SE UI

    """

    SIMPLE = "SIMPLE"
    FARQUHAR1997 = "FARQUHAR1997"
    DO3SEUI = "DO3SEUI"


class ThermalTimeMethods(Enum):
    """Methods of calculating thermal time.

    These are seperate from the thermal time methods in the external data.

    Options
    -------
    EXTERNAL: "External"
        Thermal time is calculated before run and contained in external data state
    HOURLY: "HOURLY"
        Thermal time is calculated hourly. We need to store the accumulated then average temperature
        then update the thermal time at hour 0

    """
    EXTERNAL = "EXTERNAL"
    HOURLY = "HOURLY"
    SKIP = "SKIP"


class DayCalcMethods(Enum):
    """Methods for defining the dd value.


    """

    #: Generate the dd value from start day in config
    CALC_FROM_START = "CALC_FROM_START"

    #: Convert input day month year to julian day
    CALC_FROM_INPUT = "CALC_FROM_INPUT"


@dataclass(frozen=False)
class Config_Met_Inputs:
    """Configuration for external meteorological input data.

    This is set in the "Met" section of the config file. It defines which inputs
    should be provided in an external csv file and which should be calculated internally.

    Name must use naming convention:
    [External_State_Shape_var_name]_method
    [External_State_Shape_var_name]_constant
    [External_State_Shape_var_name]_fillna
    External_State_Shape_var_name]_required: bool = False
    [External_State_Shape_var_name]_required
    E.g. CO2_method = "constant"
    CO2_constant = 391.0

    All method options can be:
    - "constant": Use constant value set in config
    - "input": Use input data value
    - "offset": Offset input data by value set in config
    - "multiply": Multiply input data by value set in config
    - "calculated": Calculated using other values before model run
    - "skip": Not needed

    Note: If it is calculated we will initially set it to None then validate it before
    running the model.

    fill_na
    -------

    If _fillna is set then any empty inputs will be replaced with this value.
    If _required: bool = False
    I.e. if CO2_fillna = 100 and input CO2 data is [20,30,40,,60] (Note missing val) then
    CO2_required: bool = False
    CO2 fixed will be [20,30,40,100,60]

    required
    -------

    If _required is true then all values must be valid after external data processing
    If an input is supplied but required is set to False then any missing values could
    be passed to the model.

    """

    time_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    row_index_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    dom_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    mm_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    year_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    dd_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    dd_required: bool = True

    hr_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    hr_required: bool = True

    CO2_method: InputMethod = field(default_factory=lambda: InputMethod.CONSTANT)
    #: Constant CO2 concentration [ppm]
    CO2_constant: float = 391.0
    CO2_fillna: float = None
    CO2_required: bool = True

    O3_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    O3_constant: float = None
    O3_fillna: float = None
    O3_required: bool = True

    Ts_C_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    Ts_C_constant: float = None
    Ts_C_fillna: float = None
    Ts_C_required: bool = True

    Ts_K_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    Ts_K_constant: float = None
    Ts_K_fillna: float = None
    Ts_K_required: bool = False

    P_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    P_constant: float = None
    P_fillna: float = None
    P_required: bool = True

    precip_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    precip_constant: float = None
    precip_fillna: float = None
    precip_required: bool = True

    u_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    u_constant: float = None
    u_fillna: float = None
    u_required: bool = True

    uh_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    uh_constant: float = None
    uh_fillna: float = None
    uh_required: bool = False

    O3_nmol_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    O3_nmol_constant: float = None
    O3_nmol_fillna: float = None
    O3_nmol_required: bool = False

    Tleaf_C_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    Tleaf_C_constant: float = None
    Tleaf_C_fillna: float = None
    Tleaf_C_required: bool = False

    u__method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    u__constant: float = None
    u__fillna: float = None
    u__required: bool = False

    Rn_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    Rn_constant: float = None
    Rn_fillna: float = None
    Rn_required: bool = False

    R_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    R_constant: float = None
    R_fillna: float = None
    R_required: bool = False

    PAR_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    PAR_constant: float = None
    PAR_fillna: float = None
    PAR_required: bool = True

    PPFD_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    PPFD_constant: float = None
    PPFD_fillna: float = None
    PPFD_required: bool = False

    Idrctt_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    Idrctt_constant: float = None
    Idrctt_fillna: float = None
    Idrctt_required: bool = False

    Idfuse_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    Idfuse_constant: float = None
    Idfuse_fillna: float = None
    Idfuse_required: bool = False

    VPD_method: InputMethod = field(default_factory=lambda: InputMethod.INPUT)
    VPD_constant: float = None
    VPD_fillna: float = None
    VPD_required: bool = True

    RH_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    RH_constant: float = None
    RH_fillna: float = None
    RH_required: bool = False

    h_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    h_constant: float = None
    h_fillna: float = None
    h_required: bool = False

    SWP_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    SWP_constant: float = None
    SWP_fillna: float = None
    SWP_required: bool = False

    SWC_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    SWC_constant: float = None
    SWC_fillna: float = None
    SWC_required: bool = False

    VPD_dd_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    VPD_dd_constant: float = None
    VPD_dd_fillna: float = None
    VPD_dd_required: bool = False

    esat_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    esat_constant: float = None
    esat_fillna: float = None
    esat_required: bool = False

    eact_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    eact_constant: float = None
    eact_fillna: float = None
    eact_required: bool = False

    td_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    td_constant: float = None
    td_fillna: float = None
    td_required: bool = False

    is_daylight_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    is_daylight_constant: bool = None
    is_daylight_fillna: float = None
    is_daylight_required: bool = False

    sinB_method: InputMethod = field(default_factory=lambda: InputMethod.CALCULATED)
    sinB_constant: float = None
    sinB_fillna: float = None
    sinB_required: bool = True

    Hd_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    Hd_constant: float = None
    Hd_fillna: float = None
    Hd_required: bool = False

    leaf_fphen_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    leaf_fphen_constant: float = None
    leaf_fphen_fillna: float = None
    leaf_fphen_required: bool = False

    fphen_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    fphen_constant: float = None
    fphen_fillna: float = None
    fphen_required: bool = False

    LAI_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    LAI_constant: float = None
    LAI_fillna: float = None
    LAI_required: bool = False

    V_cmax_25_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    V_cmax_25_constant: float = None
    V_cmax_25_fillna: float = None
    V_cmax_25_required: bool = False

    J_max_25_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    J_max_25_constant: float = None
    J_max_25_fillna: float = None
    J_max_25_required: bool = False

    snow_depth_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    snow_depth_constant: float = None
    snow_depth_fillna: float = None
    snow_depth_required: bool = False

    cloudfrac_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    cloudfrac_constant: float = None
    cloudfrac_fillna: float = None
    cloudfrac_required: bool = False

    ustar_ref_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    ustar_ref_constant: float = None
    ustar_ref_fillna: float = None
    ustar_ref_required: bool = False

    ustar_method: InputMethod = field(default_factory=lambda: InputMethod.SKIP)
    ustar_constant: float = None
    ustar_fillna: float = None
    ustar_required: bool = False


@dataclass
class Config_Met:
    """Configuration for External Met data.

    This configuration defines how the model deals with external data inputs.
    """

    td_base_temperature: float = 0.0  #: Base temperature used for calculating thermal time

    #: The day at which thermal time is 0
    #: Can be set by get_thermal_time_base_and_offset
    thermal_time_start: JulianDay = 0
    #: Offset thermal time at thermal_time_start
    thermal_time_offset: ThermalTime = 0

    #: Options for how we calculate thermal time
    thermal_time_method: ThermalTimeMethods = field(
        default_factory=lambda: ThermalTimeMethods.HOURLY)

    #: Method to calculate dd if dd_method == calculated
    dd_calc_method: DayCalcMethods = None

    #: If true then we run processes that rely on lai hourly instead of daily
    PARsunshade_method: ParSunShadeMethods = ParSunShadeMethods.FARQUHAR1997

    #: Model input fields
    inputs: Config_Met_Inputs = field(default_factory=lambda: Config_Met_Inputs())

    #: Sparse data mode
    #: This should be set if there are hours missing in the data or it doesnt start
    #: and finish at the end of the day (Hour 0 or 23)
    sparse_data: bool = False #: run in sparse data mode
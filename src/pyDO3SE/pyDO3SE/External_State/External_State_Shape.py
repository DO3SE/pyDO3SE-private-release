"""External state shape is the data model for state outside of the plant state.

It consists of imported data and calculated data based on the config and imported data.

Solar Radiation
---------------

The following parameters contribute to solar radiation:

 - Rn: Optional[List[float]] = None  #: 'MJ m-2'         Net radiation
 - cloudfrac: Optional[List[float]] = None  #: cloud fraction [Fraction]
 - R: Optional[List[float]] = None  #: 'W m-2'          Global radiation
 - PAR: Optional[List[float]] = None  #: 'W m-2'          Photosynthetically active radiation
 - PPFD: Optional[List[float]] = None  #: 'umol m-2 s-     Photosynthetic photon flux density
 - Idrctt: Optional[List[float]] = None  #: 'W m-2'          Direct PAR irradiance
 - Idfuse: Optional[List[float]] = None  #: 'W m-2'          Diffuse PAR irradiance

"""

from typing import List, Optional
from dataclasses import dataclass


# TODO: Implement not optional type
# NOT_SET = Literal["NOT_SET"]

# type NotOptional = Union[NOT_SET, List[float]]


@dataclass(frozen=False)
class External_State_Shape:
    """Shape of data from external data file (.csv).

    Make sure to add config properties to External_State_Config.py

    """

    row_index: Optional[List[int]] = None
    time: Optional[List[str]] = None  #: date time str E.g. 2017-10-01T00:00:00
    year: Optional[List[str]] = None  #: year (YYYY)
    dom: Optional[List[str]] = None  #: day of month base 1
    mm: Optional[List[str]] = None  #: month of year 1-12 (MM)
    dd: Optional[List[int]] = None  #: Day of year
    hr: Optional[List[int]] = None  #: Hour of day
    Ts_C: Optional[List[float]] = None  #: 'degrees C'      Temperature
    Ts_K: Optional[List[float]] = None  #: 'degrees Kelvin'      Temperature
    P: Optional[List[float]] = None  #: 'kPa'            Pressure
    precip: Optional[List[float]] = None  #: 'mm'             Precipitation
    u: Optional[List[float]] = None  #: 'm s-1'          Measured wind speed
    O3: Optional[List[float]] = None  #: 'ppb'            Measured ozone concentration
    uh: Optional[List[float]] = None  #: windspeed at top of target canopy
    O3_nmol: Optional[List[float]] = None  #: O3 data converted to nmol/m^3
    td: Optional[List[float]] = None  #: Thermal Time
    Tleaf_C: Optional[List[float]] = None  #: 'degrees C'      Leaf temperature
    u_: Optional[List[float]] = None  #: 'm s-1'          Friction velocity
    Rn: Optional[List[float]] = None  #: 'MJ m-2'         Net radiation
    cloudfrac: Optional[List[float]] = None  #: cloud fraction [Fraction]
    R: Optional[List[float]] = None  #: 'W m-2'          Global radiation
    PAR: Optional[List[float]] = (
        None  #: 'W m-2'          Photosynthetically active radiation
    )
    PPFD: Optional[List[float]] = (
        None  #: 'umol m-2 s-     Photosynthetic photon flux density
    )
    Idrctt: Optional[List[float]] = None  #: 'W m-2'          Direct PAR irradiance
    Idfuse: Optional[List[float]] = None  #: 'W m-2'          Diffuse PAR irradiance
    VPD: Optional[List[float]] = None  #: 'kPa'            Vapour pressure deficit
    RH: Optional[List[float]] = None  #: '1'              Relative humidity
    CO2: Optional[List[float]] = None  #: 'ppm'            Ambient CO2 concentration
    h: Optional[List[float]] = None  #: 'm'              Canopy height
    SWP: Optional[List[float]] = None  #: 'MPa'            Soil water potential
    SWC: Optional[List[float]] = None  #: '%'            Soil water content
    ASW: Optional[List[float]] = None  #: '%'            Available Soil water content
    #: 'W/m^2'            Surface Heat Flux (Should be +ve in middle of summer day)
    Hd: Optional[List[float]] = None

    VPD_dd: Optional[List[float]] = None  #: < Daily VPD sum during daylight hours [kPa]
    esat: Optional[List[float]] = None  #: < Saturated vapour pressure [kPa]
    eact: Optional[List[float]] = None  #: < Actual vapour pressure [kPa]
    is_daylight: Optional[List[bool]] = None  #: True if is daylight
    sinB: Optional[List[float]] = None  #: sin() of solar elevation angle
    leaf_fphen: Optional[List[float]] = None  #: leaf f_phen input [fraction]

    #:    TODO: Implement below with shape(DIM_L, DIM_LC)
    LAI: Optional[List[float]] = None  #: Leaf area index  [m2 m-2]
    V_cmax_25: Optional[List[float]] = (
        None  #: 'umol m-2 s-1    Maximum catalytic rate at 25 degrees
    )
    #: 'umol m-2 s-1    Maximum rate of electron transport at 25 degrees
    J_max_25: Optional[List[float]] = None

    snow_depth: Optional[List[float]] = None  #: Snow depth [m]
    ustar_ref: Optional[List[float]] = None  #: Observed friction velocity [m/s]
    ustar: Optional[List[float]] = None  #: modelled friction velocity [m/s]


#: TODO: Implement description
#: @property
#: def description(self):
#:     out = self.long_name + ' (' + self.name
#:     if self.units is not None:
#:         out += ', ' + self.units
#:     out += ')'
#:     return out
@dataclass
class InputField:
    id: str
    label: str
    dtype: type | str
    units: str | None
    alt_names: Optional[set] = None

    def __post_init__(self):
        self.alt_names = self.alt_names or set()
        self.alt_names |= {self.label, self.label.lower()}

    def to_dict(self):
        return {
            "id": self.id,
            "label": self.label,
            "dtype": self.dtype,
            "units": self.units,
            "alt_names": self.alt_names,
        }


INPUT_FIELDS = [
    InputField("row_index", "Row Index", "intc", None, {"i"}),
    InputField("time", "Datetime", "str", None, {"date"}),
    InputField("year", "Year", "intc", None, {"yr"}),
    InputField("dom", "Day of month", "intc", None, {"day_of_month"}),
    InputField("mm", "Month", "intc", None, {"month"}),
    InputField("dd", "Day of year", "intc", None, {"jd", "day"}),
    InputField("td", "Thermal Time", "double", None),
    InputField("hr", "Hour", "intc", None, {"hours", "hour"}),
    InputField(
        "Ts_C",
        "Temperature",
        "double",
        "degrees C",
        {"tsc", "air_temperature", "t2m", "temperature"},
    ),
    InputField("Ts_K", "Temperature", "double", "degrees K", {"tsk"}),
    InputField("Tleaf_C", "Leaf temperature", "double", "degrees C", {"tleaf"}),
    InputField("P", "Pressure", "double", "kPa", {"ambient_pressure", "pressure"}),
    InputField(
        "precip", "Precipitation", "double", "mm", {"prec", "precipitation", "rain"}
    ),
    InputField(
        "ustar", "Friction velocity", "double", "m s-1", {"ustar", "friction_velocity"}
    ),
    InputField("u", "Measured wind speed", "double", "m s-1", {"wind_speed", "uh_zr"}),
    InputField("O3", "Measured ozone concentration", "double", "ppb", {"o3_ppb_zr"}),
    InputField("Rn", "Net radiation", "double", "MJ m-2"),
    InputField("R", "Global radiation", "double", "W m-2"),
    InputField("cloudfrac", "Cloud fraction", "double", "fraction", {"cloud_amt"}),
    InputField(
        "PAR",
        "Photosynthetically active radiation",
        "double",
        "W m-2",
        {"photosynthetically_active_radiation", "par"},
    ),
    InputField("PPFD", "Photosynthetic photon flux density", "double", "umol m-2 s-1"),
    InputField("Idrctt", "Direct PAR irradiance", "double", "W m-2"),
    InputField("Idfuse", "Diffuse PAR irradiance", "double", "W m-2"),
    InputField("VPD", "Vapour pressure deficit", "double", "kPa"),
    InputField("RH", "Relative humidity", "double", "1", {"relative_humidity"}),
    InputField(
        "CO2", "Ambient CO2 concentration", "double", "ppm", {"co2", "carbon_dioxide"}
    ),
    InputField("h", "Canopy height", "double", "m"),
    InputField("SWP", "Soil water potential", "double", "MPa"),
    InputField("SWC", "Soil water content", "double", "%", {"soil_moisture"}),
    InputField("ASW", "Available Soil water content", "double", "%", {"asw"}),
    InputField("leaf_fphen", "Leaf fphen", "double", "fraction", {"leaf_f_phen"}),
    InputField("fphen", "Fphen", "double", "fraction", {"f_phen"}),
    InputField("LAI", "Leaf area index", "double", "m2 m-2", {"leaf_area_index"}),
    InputField(
        "ustar_ref",
        "Grid Friction velocity",
        "double",
        "fraction",
        {"ustar", "friction_velocity"},
    ),
    InputField("Hd", "Sensible Heat Flux", "double", "W m-2", {"sensible_heat_flux"}),
]


def ifield(id: str):
    """Simple utility fnc to help merging input fields with detailed fields"""
    return {"id": id, **next(f for f in INPUT_FIELDS if f.id == id).to_dict()}


@dataclass
class InputFieldDetailed:
    id: str
    label: str
    dtype: type
    units: str
    alt_names: set
    required: bool = False
    shape: Optional[tuple] = None
    long_name: Optional[str] = None
    description: str = ""

    def __post_init__(self):
        self.alt_names = self.alt_names or set()
        self.alt_names |= {self.label, self.label.lower()}
        self.long_name = self.long_name or self.label

    def __str__(self) -> str:
        return self.label


INPUT_FIELDS_DETAILED = [
    InputFieldDetailed(
        **ifield("row_index"),
        required=False,
        shape=(1,),
        description="Row Index",
    ),
    InputFieldDetailed(
        **ifield("time"),
        required=False,
        shape=(1,),
        description="Time in UTC",
    ),
    InputFieldDetailed(
        **ifield("year"),
        required=False,
        shape=(1,),
        description="Year",
    ),
    InputFieldDetailed(
        **ifield("dom"),
        required=False,
        shape=(1,),
        description="Day of month",
    ),
    InputFieldDetailed(
        **ifield("mm"),
        required=False,
        shape=(1,),
        description="Month of the year 1-12",
    ),
    InputFieldDetailed(
        **ifield("dd"),
        required=True,
        shape=(1,),
        description="Day of year",
    ),
    InputFieldDetailed(
        **ifield("td"),
        required=False,
        shape=(1,),
        description="Thermal Time",
    ),
    InputFieldDetailed(
        **ifield("hr"),
        required=True,
        shape=(1,),
        description="Hour of day",
    ),
    InputFieldDetailed(
        **ifield("Ts_C"),
        required=True,
        shape=(1,),
        description="Temperature in degrees C",
    ),
    InputFieldDetailed(
        **ifield("Ts_K"),
        required=False,
        shape=(1,),
        description="Temperature in degrees K",
    ),
    InputFieldDetailed(
        **ifield("Tleaf_C"),
        required=False,
        shape=(1,),
        description="Leaf temperature in degrees C",
    ),
    InputFieldDetailed(
        **ifield("P"),
        required=True,
        shape=(1,),
        description="Pressure in kPa",
    ),
    InputFieldDetailed(
        **ifield("precip"),
        required=True,
        shape=(1,),
        description="Precipitation in mm",
    ),
    InputFieldDetailed(
        **ifield("ustar"),
        required=False,
        shape=(1,),
        description="Friction velocity in m s-1",
    ),
    InputFieldDetailed(
        **ifield("u"),
        required=True,
        shape=(1,),
        description="Measured wind speed in m s-1",
    ),
    InputFieldDetailed(
        **ifield("O3"),
        required=True,
        shape=(1,),
        description="Measured ozone concentration in ppb",
    ),
    InputFieldDetailed(
        **ifield("Rn"),
        required=False,
        shape=(1,),
        description="Net radiation in MJ m-2",
    ),
    InputFieldDetailed(
        **ifield("R"),
        required=False,
        shape=(1,),
        description="Global radiation in W m-2",
    ),
    InputFieldDetailed(
        **ifield("cloudfrac"),
        required=False,
        shape=(1,),
        description="Cloud fraction",
    ),
    InputFieldDetailed(
        **ifield("PAR"),
        required=False,
        shape=(1,),
        description="Photosynthetically active radiation in W m-2",
    ),
    InputFieldDetailed(
        **ifield("PPFD"),
        required=False,
        shape=(1,),
        description="Photosynthetic photon flux density in umol m-2 s-1",
    ),
    InputFieldDetailed(
        **ifield("Idrctt"),
        required=False,
        shape=(1,),
        description="Direct PAR irradiance in W m-2",
    ),
    InputFieldDetailed(
        **ifield("Idfuse"),
        required=False,
        shape=(1,),
        description="Diffuse PAR irradiance in W m-2",
    ),
    InputFieldDetailed(
        **ifield("VPD"),
        required=False,
        shape=(1,),
        description="Vapour pressure deficit in kPa",
    ),
    InputFieldDetailed(
        **ifield("RH"),
        required=False,
        shape=(1,),
        description="Relative humidity",
    ),
    InputFieldDetailed(
        **ifield("CO2"),
        required=False,
        shape=(1,),
        description="Ambient CO2 concentration in ppm",
    ),
    InputFieldDetailed(
        **ifield("h"),
        required=False,
        shape=(1,),
        description="Canopy height in m",
    ),
    InputFieldDetailed(
        **ifield("SWP"),
        required=False,
        shape=(1,),
        description="Soil water potential in MPa",
    ),
    InputFieldDetailed(
        **ifield("SWC"),
        required=False,
        shape=(1,),
        description="Soil water content in %",
    ),
    InputFieldDetailed(
        **ifield("ASW"),
        required=False,
        shape=(1,),
        description="Available Soil water content in %",
    ),
    InputFieldDetailed(
        **ifield("leaf_fphen"),
        required=False,
        shape=(1,),
        description="Leaf fphen",
    ),
    InputFieldDetailed(
        **ifield("fphen"),
        required=False,
        shape=(1,),
        description="Fphen",
    ),
    InputFieldDetailed(
        **ifield("LAI"),
        required=False,
        shape=(1,),
        description="Leaf area index in m2 m-2",
    ),
    InputFieldDetailed(
        **ifield("ustar_ref"),
        required=False,
        shape=(1,),
        description="Grid Friction velocity",
    ),
    InputFieldDetailed(
        **ifield("Hd"),
        required=False,
        shape=(1,),
        description="Sensible Heat Flux in W m-2",
    ),
]


# #: #: Groups of inputs where at least one input (or set of inputs) is required
# #: INPUT_REQUIRED_GROUPS = (
# #:     ('Irradiance', ({'R'}, {'PAR'}, {'PPFD'}, {'Idrctt', 'Idfuse'})),
# #:     ('Humidity', ({'VPD'}, {'RH'})),
# #: )

# #: def check_input_fields(fields):
# #:     """Check that required input fields are present.

# #:     Raises a :exc:`ValueError` if required fields are missing.
# #:     """
# #:     missing_required = INPUT_REQUIRED - fields
# #:     if missing_required:
# #:         raise ValueError('Required inputs missing', missing_required)
# #:     for name, group in INPUT_REQUIRED_GROUPS:
# #:         if not any(s.issubset(fields) for s in group):
# #:             raise ValueError('At least one input required', name, group)

# This is a check that we have detailed fields for all the fields in the input fields
assert all(f.id == fd.id for f, fd in zip(INPUT_FIELDS, INPUT_FIELDS_DETAILED))

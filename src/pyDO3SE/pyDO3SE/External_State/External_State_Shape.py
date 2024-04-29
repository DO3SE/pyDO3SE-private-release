"""External state shape is the data model for state outside of the plant state.

It consists of imported data and calculated data based on the config and imported data.

Solar Radiation
---------------

The following parameters contribute to solar radiation:

 - Rn: List[float] = None  #: 'MJ m-2'         Net radiation
 - cloudfrac: List[float] = None  #: cloud fraction [Fraction]
 - R: List[float] = None  #: 'W m-2'          Global radiation
 - PAR: List[float] = None  #: 'W m-2'          Photosynthetically active radiation
 - PPFD: List[float] = None  #: 'umol m-2 s-     Photosynthetic photon flux density
 - Idrctt: List[float] = None  #: 'W m-2'          Direct PAR irradiance
 - Idfuse: List[float] = None  #: 'W m-2'          Diffuse PAR irradiance

"""

from typing import List
from collections import namedtuple
from dataclasses import dataclass


@dataclass(frozen=False)
class External_State_Shape:
    """Shape of data from external data file (.csv).

    Make sure to add config properties to External_State_Config.py

    """
    row_index: List[int] = None
    time: List[str] = None  #: date time str E.g. 2017-10-01T00:00:00
    year: List[str] = None  #: year (YYYY)
    dom: List[str] = None  #: day of month base 1
    mm: List[str] = None  #: month of year 1-12 (MM)
    dd: List[int] = None  #: Day of year
    hr: List[int] = None  #: Hour of day
    Ts_C: List[float] = None  #: 'degrees C'      Temperature
    Ts_K: List[float] = None  #: 'degrees Kelvin'      Temperature
    P: List[float] = None  #: 'kPa'            Pressure
    precip: List[float] = None  #: 'mm'             Precipitation
    u: List[float] = None  #: 'm s-1'          Measured wind speed
    O3: List[float] = None  #: 'ppb'            Measured ozone concentration
    uh: List[float] = None  #: windspeed at top of target canopy
    O3_nmol: List[float] = None  #: O3 data converted to nmol/m^3
    td: List[float] = None  #: Thermal Time
    Tleaf_C: List[float] = None  #: 'degrees C'      Leaf temperature
    u_: List[float] = None  #: 'm s-1'          Friction velocity
    Rn: List[float] = None  #: 'MJ m-2'         Net radiation
    cloudfrac: List[float] = None  #: cloud fraction [Fraction]
    R: List[float] = None  #: 'W m-2'          Global radiation
    PAR: List[float] = None  #: 'W m-2'          Photosynthetically active radiation
    PPFD: List[float] = None  #: 'umol m-2 s-     Photosynthetic photon flux density
    Idrctt: List[float] = None  #: 'W m-2'          Direct PAR irradiance
    Idfuse: List[float] = None  #: 'W m-2'          Diffuse PAR irradiance
    VPD: List[float] = None  #: 'kPa'            Vapour pressure deficit
    RH: List[float] = None  #: '1'              Relative humidity
    CO2: List[float] = None  #: 'ppm'            Ambient CO2 concentration
    h: List[float] = None  #: 'm'              Canopy height
    SWP: List[float] = None  #: 'MPa'            Soil water potential
    SWC: List[float] = None  #: '%'            Soil water content
    ASW: List[float] = None  #: '%'            Available Soil water content
    #: 'W/m^2'            Surface Heat Flux (Should be +ve in middle of summer day)
    Hd: List[float] = None

    VPD_dd: List[float] = None   #: < Daily VPD sum during daylight hours [kPa]
    esat: List[float] = None  #: < Saturated vapour pressure [kPa]
    eact: List[float] = None  #: < Actual vapour pressure [kPa]
    is_daylight: List[bool] = None  #: True if is daylight
    sinB: List[float] = None  #: sin() of solar elevation angle
    leaf_fphen: List[float] = None  #: leaf f_phen input [fraction]

#:    TODO: Implement below with shape(DIM_L, DIM_LC)
    LAI: List[float] = None  #: Leaf area index  [m2 m-2]
    V_cmax_25: List[float] = None  #: 'umol m-2 s-1    Maximum catalytic rate at 25 degrees
    #: 'umol m-2 s-1    Maximum rate of electron transport at 25 degrees
    J_max_25: List[float] = None

    snow_depth: List[float] = None  #: Snow depth [m]
    ustar_ref: List[float] = None  #: Observed friction velocity [m/s]
    ustar: List[float] = None  #: modelled friction velocity [m/s]


class InputField(namedtuple('_InputField', 'name required dtype shape units long_name, alt_names')):
    """Input field types to map inputs data against"""
    def __new__(cls, name, required, dtype, shape, units, long_name, alt_names=None):
        alt_names = alt_names or set()
        alt_names |= {name, name.lower()}
        return super(InputField, cls).__new__(
            cls, name, required, dtype, shape, units, long_name, alt_names)

    #: TODO: Implement description
    #: @property
    #: def description(self):
    #:     out = self.long_name + ' (' + self.name
    #:     if self.units is not None:
    #:         out += ', ' + self.units
    #:     out += ')'
    #:     return out


#: All available input fields
#: Imported external data will be compared against this to determine which data to import
#: Check here first for import data typos
#: TODO: Update this
INPUT_FIELDS = (
    #: key | Required | type | Dimensions | Unit | Description | alt_names
    InputField('row_index', False, 'intc', (), None, 'Row Index', {}),  # noqa: E241 E501
    InputField('time', True, 'str', (), None, 'Datetime', {'date'}),  # noqa: E241 E501
    InputField('year', False, 'intc', (), None, 'Year', {'yr'}),  # noqa: E241 E501
    InputField('dom', False, 'intc', (), None, 'Day of month', {'day_of_month'}),  # noqa: E241 E501
    InputField('mm', False, 'intc', (), None, 'month of the year 1-12', {'month'}),  # noqa: E241 E501
    InputField('dd', True, 'intc', (), None, 'Day of year', {'jd', 'day'}),  # noqa: E241 E501
    InputField('td', False, 'double', (), None, 'Thermal Time'),  # noqa: E241 E501
    InputField('hr', True, 'intc', (), None, 'Hour of day', {'hours', 'hour'}),  # noqa: E241 E501
    InputField('Ts_C', True, 'double', (), 'degrees C', 'Temperature', {'tsc', 'air_temperature', 't2m', 'temperature'}),  # noqa: E241 E501
    InputField('Ts_K', False, 'double', (), 'degrees K', 'Temperature', {'tsk'}),  # noqa: E241 E501
    InputField('Tleaf_C', False, 'double', (), 'degrees C', 'Leaf temperature', {'tleaf'}),  # noqa: E241 E501
    InputField('P', True, 'double', (), 'kPa', 'Pressure', {'ambient_pressure', 'pressure'}),  # noqa: E241 E501
    InputField('precip', True, 'double', (), 'mm', 'Precipitation', {'prec', 'precipitation', 'rain'}),  # noqa: E241 E501
    InputField('ustar', False, 'double', (), 'm s-1', 'Friction velocity', {'ustar', 'friction_velocity'}),  # noqa: E241 E501
    InputField('u', True, 'double', (), 'm s-1', 'Measured wind speed', {'wind_speed', 'uh_zr'}),  # noqa: E241 E501
    InputField('O3', True, 'double', (), 'ppb', 'Measured ozone concentration', {'o3_ppb_zr'}),  # noqa: E241 E501
    #: NOTE: Ensure correct PAR input data
    InputField('Rn', False, 'double', (), 'MJ m-2', 'Net radiation'),  # noqa: E241 E501
    InputField('R', False, 'double', (), 'W m-2', 'Global radiation'),  # noqa: E241 E501
    InputField('cloudfrac', False, 'double', (), 'fraction', 'Cloud fraction', {'cloud_amt'}),  # noqa: E241 E501
    InputField('PAR', False, 'double', (), 'W m-2', 'Photosynthetically active radiation', {'photosynthetically_active_radiation','par'}),  # noqa: E241 E501
    InputField('PPFD', False, 'double', (), 'umol m-2 s-1', 'Photosynthetic photon flux density'),  # noqa: E241 E501
    InputField('Idrctt', False, 'double', (), 'W m-2', 'Direct PAR irradiance'),  # noqa: E241 E501
    InputField('Idfuse', False, 'double', (), 'W m-2', 'Diffuse PAR irradiance'),  # noqa: E241 E501
    InputField('VPD', False, 'double', (), 'kPa', 'Vapour pressure deficit'),  # noqa: E241 E501
    InputField('RH', False, 'double', (), '1', 'Relative humidity', {'relative_humidity'}),  # noqa: E241 E501
    InputField('CO2', False, 'double', (), 'ppm', 'Ambient CO2 concentration', {'co2', 'carbon_dioxide'}),  # noqa: E241 E501
    InputField('h', False, 'double', (), 'm', 'Canopy height'),  # noqa: E241 E501
    InputField('SWP', False, 'double', (), 'MPa', 'Soil water potential'),  # noqa: E241 E501
    InputField('SWC', False, 'double', (), '%', 'Soil water content', {'soil_moisture'}),  # noqa: E241 E501
    InputField('ASW', False, 'double', (), '%', 'Available Soil water content', {'asw'}),  # noqa: E241 E501
    InputField('leaf_fphen', False, 'double', (), 'fraction', 'leaf fphen', {'leaf_f_phen'}),  # noqa: E241 E501
    InputField('fphen', False, 'double', (), 'fraction', 'fphen', {'f_phen'}),  # noqa: E241 E501
    InputField('LAI', False, 'double', (), 'm2 m-2', 'Leaf area index', {'leaf_area_index'}),  # noqa: E241 E501
    InputField('ustar_ref', False, 'double', (), 'fraction', 'Grid Friction velocity', {'ustar', 'friction_velocity'}),  # noqa: E241 E501
    InputField('Hd', False, 'double', (), 'W m-2', 'Sensible Heat Flux', {'sensible_heat_flux'}),  # noqa: E241 E501

)


#: #: Inputs that are always required
#: INPUT_REQUIRED = {f.name for f in INPUT_FIELDS if f.required}


#: #: Groups of inputs where at least one input (or set of inputs) is required
#: INPUT_REQUIRED_GROUPS = (
#:     ('Irradiance', ({'R'}, {'PAR'}, {'PPFD'}, {'Idrctt', 'Idfuse'})),
#:     ('Humidity', ({'VPD'}, {'RH'})),
#: )

#: def check_input_fields(fields):
#:     """Check that required input fields are present.

#:     Raises a :exc:`ValueError` if required fields are missing.
#:     """
#:     missing_required = INPUT_REQUIRED - fields
#:     if missing_required:
#:         raise ValueError('Required inputs missing', missing_required)
#:     for name, group in INPUT_REQUIRED_GROUPS:
#:         if not any(s.issubset(fields) for s in group):
#:             raise ValueError('At least one input required', name, group)

from enum import Enum
from typing import TypeVar, Union, List, Tuple

UX = TypeVar('UX')
UY = TypeVar('UY')

JulianDay = int
ThermalTime = float
TimeUnit = Union[JulianDay, ThermalTime]
# PiecewiseFunction = List[Tuple[UX, UY]]
PiecewiseFunction = Tuple[List[UX], List[UY]]

ThermalTime = float
Fraction = float


class TimeTypes(Enum):

    JULIAN_DAY = "julian_day"
    THERMAL_TIME = "thermal_time"

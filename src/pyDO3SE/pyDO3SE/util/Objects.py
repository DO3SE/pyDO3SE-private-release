from typing import NamedTuple, Type, Union
from enum import Enum


class Unit(Enum):
    _unknown = 0
    degc = "Degrees Celcius"


class Field(NamedTuple):
    id: str
    type: Type
    short: str = ""
    unit: Union[Unit, str] = Unit._unknown
    long: str = ""
    required: bool = False
    constant: bool = True
    per_iC: bool = False
    per_iL: bool = False
    per_iP: bool = False
    per_iCH: bool = False

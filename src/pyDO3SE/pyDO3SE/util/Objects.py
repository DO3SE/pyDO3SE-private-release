from typing import NamedTuple, Type
from enum import Enum


class Unit(Enum):
    _unknown = 0
    degc = 'Degrees Celcius'


class Field(NamedTuple):
    id: str
    type: Type
    short: str = ''
    unit: Unit = Unit._unknown
    long: str = ''
    required: bool = False
    constant: bool = True
    per_component: bool = False
    per_layer: bool = False
    per_population: bool = False

"""Phenology specific state."""

from dataclasses import dataclass, field
from enum import Enum
from functools import total_ordering

def enum_ordering(cls):
    def __lt__(self, other):
        if type(other) == type(self):
            return self.value < other.value

        raise ValueError("Cannot compare different Enums")

    setattr(cls, '__lt__', __lt__)
    return total_ordering(cls)

def enum_addition(cls):
    def __add__(self, other):
        index = int(self.value.split("_")[0])
        new_index = index + other
        new_stage = next(v for v in cls if int(v.value.split('_')[0]) == new_index)
        print(new_stage)
        return new_stage
    setattr(cls, '__add__', __add__)
    return cls


@enum_addition
@enum_ordering
class PhenologyStage(Enum):
    """Phenology stages.

    The stage refers to the last key milestone they passed.

    """

    NOT_SOWN = "0_NOT_SOWN"
    SOWN = "1_SOWN"
    EMERGED = "2_EMERGED"
    ASTART = "3_ASTART"
    HARVEST = "4_HARVEST"

    # def __add__(self, other):
    #     index = int(self.value.split("_")[0])
    #     new_index = index + other
    #     new_stage = next(v for v in PhenologyStage if int(v.value.split('_')[0]) == new_index)
    #     return new_stage

@enum_addition
@enum_ordering
class LeafPhenologyStage(Enum):
    """Leaf phenology stages.

    The stage refers to the last key milestone they passed.

    """

    NOT_EMERGED = "0_NOT_EMERGED"
    GROWING = "1_GROWING"
    MATURE = "2_MATURE"
    SENESCENCE = "3_SENESCENCE"
    FULLY_SENESED = "4_FULLY_SENESED"

@dataclass
class PhenologyState:

    phenology_stage: PhenologyStage = field(default_factory=lambda: PhenologyStage.NOT_SOWN)

    # Vernalisation
    V_tot: float = 0
    V_acc: float = 0
    Vf: float = 1

@dataclass
class LeafPhenologyState:

    phenology_stage: LeafPhenologyStage = field(default_factory=lambda: LeafPhenologyStage.NOT_EMERGED)

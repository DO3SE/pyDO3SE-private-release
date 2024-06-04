from dataclasses import dataclass, field
from typing import List
import numpy as np

from pyDO3SE.constants.model_constants import MAX_LAYERS


@dataclass(frozen=False)
class Resistance_Model:
    """Resistance model.  All LAI-related resistances are bulk.
    F90: ResistanceModel_t
    """

    nL: int     # Number of layers used
    Ra_c: float = None  # Aerodynamic resistance [s m-1] between 50m and
    # inside the canopy.
    Ra: float = None    # Aerodynamic resistance [s m-1] between 50m and the
    # reference height for this model (e.g. canopy
    # height or measurement height for O3).
    Rb: float = None    # Quasi-laminar boundary layer resistance [s m-1]
    # In-canopy aerodynamic resistance [s m-1] per layer
    Rinc: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))
    # External plant cuticle resistance [s m-1] per layer
    Rext: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))
    # Bulk Stomatal resistance [s m-1]) per layer
    Rsto: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))
    # Ground surface resistance [s m-1]
    Rgs: float = None
    # Combined surface resistance [s m-1] per layer
    Rsur: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))
    # Total resistance for each layer downwards [s m-1]
    Rtotal: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))


@dataclass(frozen=False)
class Leaf_Resistance_Model:
    """Leaf-level resistance model.  All LAI-related resistances are mean.
    F90: LeafResistanceModel_t
    """
    Rb: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))  # Leaf boundary layer resistances [s m-1]
    Rext: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))  # Leaf external plant cuticle resistance [s m-1]
    Rsto: List[float] = field(default_factory=lambda: np.full((MAX_LAYERS,), None))  # Leaf stomatal resistance [s m-1]

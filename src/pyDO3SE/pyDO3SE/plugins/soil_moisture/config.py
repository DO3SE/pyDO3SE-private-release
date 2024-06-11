"""Soil moisture parameters.

"""
from dataclasses import dataclass, field


@dataclass(frozen=False)
class Soil_t:
    """SMD soil texture parameters."""
    b: float = None  #: Texture dependent soil conductivity parameter
    FC: float = None  #: Field capacity (m3 m-3)
    SWP_AE: float = None  #: Water potential at air entry [MPa]
    Ksat: float = None  #: Saturated soil conductance (s-2 MPa-1 mm-1)


#: Commonly used soil textures:
#: Sandy loam
SOIL_SANDY_LOAM: Soil_t = \
    Soil_t(b=3.31, FC=0.16, SWP_AE=-0.00091, Ksat=0.0009576)
#: Commonly used soil textures:
#: Silt loam
SOIL_SILT_LOAM: Soil_t = \
    Soil_t(b=4.38, FC=0.26, SWP_AE=-0.00158, Ksat=0.0002178)
#: Commonly used soil textures:
#: Loam
SOIL_LOAM: Soil_t = \
    Soil_t(b=6.58, FC=0.29, SWP_AE=-0.00188, Ksat=0.0002286)
#: Commonly used soil textures:
#: Clay loam (Ksat estimated)
SOIL_CLAY_LOAM: Soil_t = \
    Soil_t(b=7.00, FC=0.37, SWP_AE=-0.00588, Ksat=0.00016)


@dataclass(frozen=False)
class Soil_Moisture_Config:
    """Soil moisture properties
    SMDConfig_t"""
    #: Soil texture:
    #:    - "sandy loam"
    #:    - "silt loam"
    #:    - "loam"
    #:    - "clay loam"
    #:    - "" or "custom":   Soil texture parameters set individually
    soil_texture: str = "loam"

    #: Soil texture parameters
    soil_config: Soil_t = field(default_factory=lambda: Soil_t())

    root: float = 1.2         #: Root depth[m]
    # TODO: seed PWP from SWP_min if available?
    PWP: float = -4.0         #: "Permanent wilting point", minimum level for SWP [MPa]
    ASW_FC: float = None     #: Available soil water at field capacity[m]
    #: (calculated by check_SMDConfig)

    #: Soil water data source:
    #:    - "disabled":   No soil water data
    #:    - "input SWP":  Input soil water potential [MPa]
    #:    - "input SWC":  Input soil water volumetric content (m3 m-3)
    #:    - "P-M":        Use Penman-Monteith method to track soil water content
    source: str = "disabled"
    #: Penman-Monteith method parameters
    initial_SWC: float = None      #: Initial soil water content, defaults to soil%FC
    run_off_fraction: float = 0.0   #: Irrigation inefficiency lost to run-off [fraction]

    #: (WEAP) "Management Allowable Depletion"[m], maximum SMD before irrigation is requested.
    MAD: float = 0.0

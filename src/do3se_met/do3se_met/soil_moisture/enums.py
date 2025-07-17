from enum import Enum


class SoilMoistureSource(Enum):
    """Methods to input soild moisture data.


    #: Soil water data source:
    #:    - "disabled":   No soil water data
    #:    - "input SWP":  Input soil water potential [MPa]
    #:    - "input SWC":  Input soil water volumetric content (m3 m-3)
    #:    - "P-M":        Use Penman-Monteith method to track soil water content
    source: str = "disabled"

    """
    DISABLED = "disabled"
    INPUT_SWP = "input SWP"
    INPUT_SWC_EXTERNAL = "external input SWC"
    INPUT_SWC = "input SWC"
    INPUT_ASW_EXTERNAL = "external input ASW"
    P_M = "P-M"


class FSW_Methods(Enum):
    DISABLED = "disabled"
    FSWP_EXP = "fSWP exp"
    FLWP_EXP = "fLWP exp"
    FSWP_LINEAR = "fSWP linear"
    FPAW = "fPAW"

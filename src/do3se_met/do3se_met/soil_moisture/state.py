from dataclasses import dataclass


@dataclass(frozen=False)
class SMDData_t:
    """Soil Moisture model data."""

    Sn: float = None        # Soil water content [m3 m-3]
    SWP: float = None       # Soil water potential [MPa]
    ASW: float = None       # Available soil water, above PWP [m]
    SMD: float = None       # Soil moisture deficit, below FC [m]


@dataclass(frozen=False)
class PM_State_t:
    """Penman monteith state."""

    PEt_rsto: float = None  # Potential rsto [m/s]

    # Latest evapotranspiration values
    Ei_hr: float = 0.0          # Evaporation from canopy [m]
    Et_hr: float = 0.0          # Plant transpiration [m]
    PEt_hr: float = 0.0          # Potential transpiration from canopy [m]
    Es_hr: float = 0.0          # Soil surface evaporation [m]
    Eat_hr: float = 0.0         # Evapotranspiration [m]

    # Accumulated values (processed and cleared daily)
    Ei_acc: float = 0.0      # Accumulated evaporation from canopy [m]
    Et_acc: float = 0.0      # Accumulated plant transpiration [m]
    PEt_acc: float = 0.0      # Accumulated potential transpiration [m]
    Es_acc: float = 0.0      # Accumulated soil surface evaporation [m]
    Eat_acc: float = 0.0     # Accumulated evapotranspiration [m]
    precip_acc_prev_day: float = 0.0  # Prev day final Accumulated precipitation for current day [m]
    precip_acc_dd: float = 0.0  # Current day Accumulated precipitation [m]

    # Soil water content tracking
    Sn_diff: float = 0.0     # Latest change in Sn [m3 m-3]

    # Water destination tracking (updated daily)
    rain_input: float = 0.0           # Input from rainfall + irrigation [m]
    run_off: float = 0.0         # Loss to run-off [m]
    effective_irrig: float = 0.0       # Effective irrigation [m]
    intercepted_evaporated: float = 0.0  # Loss to evaporation of intercepted [m]
    evapotranspiration: float = 0.0      # Loss to evapotranspiration [m]
    percolated: float = 0.0      # Loss to deep percolation [m]

    run_off_acc: float = 0.0     # Accumulated run-off [m]
    percolated_acc: float = 0.0  # Accumulated deep percolation [m]

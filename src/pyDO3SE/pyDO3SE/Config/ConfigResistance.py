from dataclasses import dataclass


@dataclass
class ResistanceConfig:
    """Config relating to resistance module.

    Attributes
    ----------
    ra_calc_method:{'simple', 'heat_flux'}
        Ra calculation method
            - simple: use ra_simple
            - heat_flux: use ra_heat_flux

    """

    #: Ra calculation method
    #:  - "simple" use ra_simple
    #:  - "heat_flux" use ra_simple
    ra_calc_method: str = 'simple'  #: simple or heat_flux

    #: Rsur calculation method
    #:  - "single_layer"
    #:  - "multi_layer"
    rsur_calc_method: str = 'single_layer'

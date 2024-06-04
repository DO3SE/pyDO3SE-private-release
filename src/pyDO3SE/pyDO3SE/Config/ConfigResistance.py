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

    # Ra calculation method
    # #  - "simple" use ra_simple
    # #  - "heat_flux" use ra_simple
    ra_calc_method: str = 'simple'  # simple or heat_flux

    # Rsur calculation method
    # #  - "single_layer"
    # #  - "multilayer"
    rsur_calc_method: str = 'single_layer'

    # Method of calculating Rext. Options:
    # - "const" - use Rext_const value
    # - "calculate_SAI" - calculate using Rext base * SAI
    rext_calc_method: str = 'const'

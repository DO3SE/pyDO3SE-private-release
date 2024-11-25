"""Carbon allocation module configuration parameters.

Initial parameters taken from Osborne et al 2015 - Jules Crop A parameterization.

"""
from dataclasses import dataclass


@dataclass(frozen=False)
class CarbonAllocationConfig:
    """Carbon allocation module config.

    References
    ----------
    Osborne et al 2015 - Jules Crop A parameterization

    """
    #: If true then allocate net primary photosynthesis to carbon pools
    use_carbon_allocation: bool = False
    #: Parameters for partition fractions - a_root
    a_root: float = None
    #: Parameters for partition fractions - a_stem
    a_stem: float = None
    #: Parameters for partition fractions - a_leaf
    a_leaf: float = None
    #: Parameters for partition fractions - b_root
    b_root: float = None
    #: Parameters for partition fractions - b_stem
    b_stem: float = None
    #: Parameters for partition fractions - b_leaf
    b_leaf: float = None

    #: additional carbon allocation parameter
    gamma: float = None
    #: additional carbon allocation parameter
    delta: float = None
    #: additional carbon allocation parameter
    theta: float = None  #: pr_resv =
    #: additional carbon allocation parameter
    k: float = None
    #: additional carbon allocation parameter
    lambdav: float = None

    #: fraction of leaf carbon that is retained in the green leaf per day
    #: during senescence [fraction] (default is 0.95)
    f_green_leaf: float = 0.95
    #: fraction of remaining leaf carbon that is retained in the brown leaf per day
    #: during senescence. Remainder goes to harvest [fraction] (default is 0.95)
    f_brown_leaf: float = 0.85

    #: Yield parameters

    #: https://edis.ifas.ufl.edu/publication/AG442, by default 1/0.84
    dry_to_wet_biomass: float = 1 / 0.84

    #: ratio of total grain mass in an ear to the ear mass
    #: avg ratio in Bangor_2015 dataset is 0.7, use this as a first guess then calibrate?,
    #: by default 0.75
    grain_to_ear: float = 0.75


#: Presets
#: Osborne et al 2015 - Jules Crop A parameterization
Wheat = CarbonAllocationConfig(
    a_root=18.5,
    a_stem=16.0,
    a_leaf=18.0,
    b_root=-20.0,
    b_stem=-15.0,
    b_leaf=-18.5,
    gamma=27.3,
    delta=-0.0507,
    theta=0.4,
    k=1.4,
    lambdav=0.4,
)

#: Presets
#: Osborne et al 2015 - Jules Crop A parameterization
Soybean = CarbonAllocationConfig(
    a_root=20.0,
    a_stem=18.5,
    a_leaf=19.5,
    b_root=-16.5,
    b_stem=-14.5,
    b_leaf=-15.0,
    gamma=25.9,
    delta=-0.1451,
    theta=0.18,
    k=1.6,
    lambdav=0.4,
)

#: Presets
#: Osborne et al 2015 - Jules Crop A parameterization
Maize = CarbonAllocationConfig(
    a_root=13.5,
    a_stem=12.5,
    a_leaf=13.0,
    b_root=-15.5,
    b_stem=-12.5,
    b_leaf=-14.0,
    gamma=22.5,
    delta=-0.2587,
    theta=0.35,
    k=3.5,
    lambdav=0.4,
)

#: Presets
#: Osborne et al 2015 - Jules Crop A parameterization
Rice = CarbonAllocationConfig(
    a_root=18.5,
    a_stem=19.0,
    a_leaf=19.5,
    b_root=-19.0,
    b_stem=-17.0,
    b_leaf=-18.5,
    gamma=20.9,
    delta=-0.2724,
    theta=0.25,
    k=1.4,
    lambdav=0.4,
)

#: List of cultivars with parameters
cultivars = ["Wheat", "Soybean", "Maize", "Rice"]

"""Carbon allocation module configuration parameters.

Initial parameters taken from Osborne et al 2015 - Jules Crop A parameterization.

"""
from dataclasses import dataclass


@dataclass(frozen=False)
class CarbonAllocationConfig:
    """Carbon allocation module config.


    """
    use_carbon_allocation: bool = False
    # Parameters for partition fractions
    a_root: float = None
    a_stem: float = None
    a_leaf: float = None
    b_root: float = None
    b_stem: float = None
    b_leaf: float = None

    # additional parameters
    gamma: float = None
    delta: float = None
    theta: float = None  # pr_resv =
    k: float = None
    lambdav: float = None


# Presets
# Osborne et al 2015 - Jules Crop A parameterization
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

cultivars = ["Wheat", "Soybean", "Maize", "Rice"]

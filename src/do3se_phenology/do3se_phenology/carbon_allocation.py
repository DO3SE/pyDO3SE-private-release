from collections import namedtuple
from math import exp

CarbonPartitionCoefficients = namedtuple(
    'CarbonPartitionCoefficients', 'p_root p_leaf p_stem p_harv')


def calc_partition_coefficients(
    DVI: float,  # -1 -> 2
    a_root: float,
    a_leaf: float,
    a_stem: float,
    b_root: float,
    b_leaf: float,
    b_stem: float,
) -> CarbonPartitionCoefficients:
    """Define partition coefficients based on DVI and cultivar parameters


    Here we define the partition coefficients as a function of
    thermal time using six parameters to describe continuously
    varying partition coefficients over the duration of the crop
    cycle" Osborne et al 2015


    References
    ----------

    - Eq6. Osborne et al 2015


    Parameters
    ----------
    DVI : float
        Development Index
    a_leaf : float
        partition coefficient
    a_stem : float
        partition coefficient
    b_root : float
        partition coefficient
    b_leaf : float
        partition coefficient
    b_stem : float
        partition coefficient

    Returns
    -------
    Tuple[float, float, float, float]
        Partition fractions
        p_root p_leaf p_stem p_harv

    """
    bottom_of_eq = exp(a_root + (b_root * DVI)) + exp(a_stem +
                                                      (b_stem * DVI)) + exp(a_leaf + (b_leaf * DVI)) + 1
    p_root = (exp(a_root + (b_root * DVI))) / bottom_of_eq
    p_leaf = (exp(a_leaf + (b_leaf * DVI))) / bottom_of_eq
    p_stem = (exp(a_stem + (b_stem * DVI))) / bottom_of_eq
    p_harv = 1 / bottom_of_eq
    return CarbonPartitionCoefficients(p_root, p_leaf, p_stem, p_harv)

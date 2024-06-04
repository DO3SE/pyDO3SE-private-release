"""Functions relating to canopy structure (Height, LAI and SAI)."""

from data_helpers.list_helpers import flatten_list
from do3se_phenology.config import SpeciesConfig
from math import ceil
from typing import List, Tuple
import numpy as np

from .utils import PLF_value
from .utils import offset
from do3se_phenology.phyllochron_dvi import calc_emergence_rate, get_growing_populations


def height_method_input():
    """Do nothing if using input."""
    pass


def height_method_constant(
    primary_land_cover_height: float,
) -> float:
    """Use height of primary land cover if height method is constant."""
    return primary_land_cover_height


def LAI_method_estimate_total(
    dd: int,
    nL: int,
    nLC: int,
    SGS: int,
    EGS: int,
    LAI_1: float,
    LAI_2: float,
    LAI_a: float,
    LAI_b: float,
    LAI_c: float,
    LAI_d: float,
    fLAI: List[List[float]]
) -> List[List[float]]:
    """Estimation method for calculating LAI.

    Use primary land cover's estimate of total LAI
    and spread over layers and LCs.

    Uses :func:`LAI_day_PLF`

    Parameters
    ----------
    dd : int
        day of year
    nL : int
        number of layers
    nLC : int
        number of components
    SGS : int
        Start of growing season
    EGS : int
        End of growing season
    LAI_1 : float
        Time from LAI_a to LAI_b [days]
    LAI_2 : float
        Time from LAI_c to LAI_d [days]
    LAI_a : float
        LAI value at SGS [m2 m-2]
    LAI_b : float
        LAI value at SGS + LAI_1 [m2 m-2]
    LAI_c : float
        LAI value at EGS - LAI_2 [m2 m-2]
    LAI_d : float
        LAI value at EGS [m2 m-2]
    fLAI : List[List[float]]
        Distribution of LAI/SAI in canopy (default: uniform distribution)

    Returns
    -------
    List[List[float]] of shape(nL, nLC)
        LAI per component per layer

    """
    # Use primary land cover's estimate of total
    t = LAI_day_PLF(dd, SGS, EGS, LAI_1, LAI_2, LAI_a, LAI_b, LAI_c, LAI_d)

    layer_LAI = [[t * fLAI[iL][iLC] for iLC in range(nLC)] for iL in range(nL)]
    return layer_LAI


def LAI_method_estimate_canopy_total(
    dd: int,
    SGS: int,
    EGS: int,
    LAI_1: float,
    LAI_2: float,
    LAI_a: float,
    LAI_b: float,
    LAI_c: float,
    LAI_d: float,
) -> float:
    """Estimation method for calculating LAI.

    Use primary land cover's estimate of total LAI

    Uses :func:`LAI_day_PLF`

    Parameters
    ----------
    dd : int
        day of year
    nLC : int
        number of components
    SGS : int
        Start of growing season
    EGS : int
        End of growing season
    LAI_1 : float
        Time from LAI_a to LAI_b [days]
    LAI_2 : float
        Time from LAI_c to LAI_d [days]
    LAI_a : float
        LAI value at SGS [m2 m-2]
    LAI_b : float
        LAI value at SGS + LAI_1 [m2 m-2]
    LAI_c : float
        LAI value at EGS - LAI_2 [m2 m-2]
    LAI_d : float
        LAI value at EGS [m2 m-2]
    fLAI : List[List[float]]
        Distribution of LAI/SAI in canopy (default: uniform distribution)

    Returns
    -------
    float
        Canopy LAI

    """
    # Use primary land cover's estimate of total
    assert SGS is not None, "SGS is None"
    assert EGS is not None, "EGS is None"
    t = LAI_day_PLF(dd, SGS, EGS, LAI_1, LAI_2, LAI_a, LAI_b, LAI_c, LAI_d)
    return t


def LAI_day_PLF(
    dd: int,
    SGS: int,
    EGS: int,
    LAI_1: float,
    LAI_2: float,
    LAI_a: float,
    LAI_b: float,
    LAI_c: float,
    LAI_d: float,
) -> float:
    r"""Estimate LAI based on a piecewise linear function of the day of year.

    To handle all situations, including winter growing seasons, everything is
    re-indexed to be relative to SGS=0.

    Interpolation wraps around between SGS and EGS to handle LAI seasons that
    aren't a simple "bump".

    TODO: The code in DO3SE UI is a much clearer implementation of this.


    .. code-block:: python

        ==========================================================================
         LAI calculation based on growing season - uses polygon calculation

             |        LAI_b                        LAI_c
             |          ______________________________
          L  |         /                              \
          A  |        /                                \
          I  |       /                                  \
             |      /                                    \
             |_____/            <1> = LAI_1               \_____
             |  LAI_a                                   LAI_d
             +--------------------------------------------------
                  |<1> |                              | <2>|
                 SGS                                      EGS
                               Day of year
        ==========================================================================


    Parameters
    ----------
    dd : int
        day of year
    SGS : int
        Start of growing season
    EGS : int
        End of growing season
    LAI_1 : float
        Time from LAI_a to LAI_b [days]
    LAI_2 : float
        Time from LAI_c to LAI_d [days]
    LAI_a : float
        LAI value at SGS [m2 m-2]
    LAI_b : float
        LAI value at SGS + LAI_1 [m2 m-2]
    LAI_c : float
        LAI value at EGS - LAI_2 [m2 m-2]
    LAI_d : float
        LAI value at EGS [m2 m-2]

    Returns
    -------
    lai: float
        Leaf area index
    """
    # x values
    # NOTE: Below assumes EGS is < SGS + 366
    GS_values = [SGS, (SGS + LAI_1), (EGS - LAI_2), EGS, (SGS + 366)]
    # Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    GS_offset = offset(GS_values, float(SGS), 365.0)

    # y values
    LAI_values = [LAI_a, LAI_b, LAI_c, LAI_d, LAI_a]

    gs_values_in_size_order = all([a < b for a, b in zip(GS_offset[0: 4], GS_offset[1: 5])])
    if not gs_values_in_size_order:
        raise ValueError(f"LAI_day_PLF: points not in order, {GS_offset}")

    dd_adj = dd - SGS if SGS - dd <= 0 else dd - SGS + 365

    func = [GS_offset, LAI_values]
    # Lookup value in PLF
    result = PLF_value(func, float(dd_adj))
    return result


def calc_SAI_wheat(
    dd: int,
    LAI: float,
    SGS: int,
    EGS: int,
    LAI_1: float,
) -> float:
    """Wheat stand area index based on the growing season.

    Parameters
    ----------
    dd : int
        day of year
    LAI: float
        Leaf area index [m2 m-2]
    SGS : int
        Start of growing season
    EGS : int
        End of growing season
    LAI_1 : float
        Time from LAI_a to LAI_b [days]

    Returns
    -------
    SAI: float
        Stand area index
    """
    # x values
    GS_values = [SGS, (SGS + LAI_1), (EGS + 1), (SGS + 365)]
    GS_offset = offset(GS_values, float(SGS), 365.0)
    # y_values
    LAI_values = [LAI, LAI + ((5.0 / 3.5) - 1) * LAI, LAI + 1.5, LAI]

    dd_adj = dd - SGS if SGS - dd <= 0 else dd - SGS + 365

    # TODO: Check this range (should it start at 0)
    result = None
    for i in range(0, 4):
        if (dd_adj < GS_offset[i]):
            result = LAI_values[i]
            break

    assert result is not None  # "SAI_wheat: no result"
    return result


def SAI_wheat_and_LAI(
    dd: int,
    SGS: int,
    EGS: int,
    LAI_1: float,
    LAI_2: float,
    LAI_a: float,
    LAI_b: float,
    LAI_c: float,
    LAI_d: float,
) -> float:
    """Calculate wheat SAI.

    Parameters
    ----------
    dd : int
        day of year
    SGS : int
        Start of growing season
    EGS : int
        End of growing season
    LAI_1 : float
        Time from LAI_a to LAI_b [days]
    LAI_2 : float
        Time from LAI_c to LAI_d [days]
    LAI_a : float
        LAI value at SGS [m2 m-2]
    LAI_b : float
        LAI value at SGS + LAI_1 [m2 m-2]
    LAI_c : float
        LAI value at EGS - LAI_2 [m2 m-2]
    LAI_d : float
        LAI value at EGS [m2 m-2]


    Returns
    -------
    SAI: float
        Stand area index
    """
    # TODO: Why do we recalculate LAI here?
    LAI = LAI_day_PLF(dd, SGS, EGS, LAI_1, LAI_2, LAI_a, LAI_b, LAI_c, LAI_d)
    SAI_wheat = calc_SAI_wheat(dd, LAI, SGS, EGS, LAI_1)
    return SAI_wheat


def calc_distribution_of_LAI_between_lcs(
    LAI_values: List[float],
    nL: int,
    nLC: int,
) -> List[float]:
    """Calculate the distribution of LAI between components.

    Parameters
    ----------
    LAI_values : List[List[float]]
        LAI per layer per component, shape (nL, nLC)
    nL : int
        Number of layers
    nLC : int
        Number of components

    Returns
    -------
    List[float]
        fraction of LAI for each component, shape (nLC)

    """
    LAI_values_np = np.array(LAI_values)
    LAI_total = np.sum(LAI_values_np)
    assert len(LAI_values) == nL, "LAI_values length not equal to nL"
    assert len(LAI_values[0]) == nLC, "LAI_values[0] length not equal to nLC"
    if LAI_total <= 0.0:
        LC_dist = [1.0 / (nL * nLC) for i in range(nLC)]
    else:
        # fraction of component total LAI over total LAI
        LAI_total_component = np.sum(LAI_values_np, 0)
        LC_dist = [lai / LAI_total for lai in LAI_total_component]
    return LC_dist


def distribute_lai_per_layers(
    total_lai: float,
    nL: int,
    max_lai_per_layer: float,
):  # TODO: Returns np array
    """Distribute lai per layers based on max lai per layer.


    """
    layers_lai = np.zeros(nL)
    if total_lai > nL * max_lai_per_layer:
        raise ValueError(
            "Not enough layers for LAI. Try increasing the layer count or LAI per layer.")
    for iL in range(nL):
        if iL == 0:
            layers_lai[int(iL)] = min(total_lai, max_lai_per_layer)
        else:
            sum_prev_layers = sum(layers_lai[0:int(iL)])
            layers_lai[int(iL)] = min(max(0, total_lai - sum_prev_layers), max_lai_per_layer)

    return layers_lai


# def get_growing_populations(
#     total_emerged_leaves: int,
#     max_leaf_lai: float,
#     prev_leaf_lais: List[float],
#     nP: int,
# ) -> List[bool]:
#     """Get a boolean list of growing leaf populations.

#     Parameters
#     ----------
#     total_emerged_leaves : int
#         then number of leaf populations that have emerged
#     max_leaf_lai : float
#         The maximum lai per leaf population
#     prev_leaf_lais : List[float]
#         previous hour leaf population lai
#     nP : int
#         number of leaf populations

#     Returns
#     -------
#     List[bool]
#         Boolean list of growing leaf populations (True if growing)

#     """
#     total_growing_populations = 0
#     growing_populations = [False for _ in range(nP)]
#     for i, iP in enumerate(range(nP)):
#         prev_leaf_lai = prev_leaf_lais[iP]
#         is_growing = True
#         if prev_leaf_lai >= max_leaf_lai:
#             is_growing = False
#         if i >= total_emerged_leaves:
#             is_growing = False

#         if is_growing:
#             total_growing_populations += 1
#             growing_populations[iP] = True

#     return growing_populations


LeafPopShape = List[List[float]]


def calc_leaf_pops_per_layer(
    prev_leaf_pops_per_lai: LeafPopShape,
    canopy_lai: float,
    layers_lai: List[float],
    nP: int,
    nL: int,
    growing_populations: List[bool],
) -> LeafPopShape:
    """Calculate the lai per leaf population per layer.

    Parameters
    ----------
    prev_leaf_pops_per_lai : LeafPopShape
        lai per leaf pop per layer
    canopy_lai : float
        total canopy lai [m2/m2]
    prev_canopy_lai : float
        total canopy lai from previous hour [m2/m2]
    layers_lai : List[float]
        layer current lai
    nP : int
        Number of leaf populations
    nL : int
        Number of layers
    growing_populations: List[bool]
        List of growing population (True is growing)

    Returns
    -------
    LeafPopShape
        lai per leaf pop per layer.
        Numpy array with shape (nL,nP)

    """
    prev_leaf_lais = [sum([prev_leaf_pops_per_lai[iL][iP] for iL in range(nL)])
                      for iP in range(nP)]  # sum of population over all layers
    total_growing_populations = sum(growing_populations)
    if total_growing_populations == 0:
        return prev_leaf_pops_per_lai

    # Limiting to min 0 stops population lai reducing
    # increase_in_canopy_lai = canopy_lai - prev_canopy_lai
    increase_in_canopy_lai = max(0, canopy_lai - sum(prev_leaf_lais))

    if increase_in_canopy_lai == 0:
        return prev_leaf_pops_per_lai
    new_lai_per_pop = increase_in_canopy_lai / \
        total_growing_populations if total_growing_populations > 0 else 0

    remaining_layer_lai = layers_lai.copy()
    lai_per_leaf_pop_per_layer = np.zeros((nL, nP))

    # TODO: Should we only distribute to last growing leaf?
    # growing_leaf = next(nP - i - 1 for i, g in enumerate(reversed(growing_populations)) if g)

    for iP in range(nP):
        prev_leaf_lai = prev_leaf_lais[iP]
        is_growing = growing_populations[iP]
        new_leaf_lai = prev_leaf_lai + new_lai_per_pop if is_growing else prev_leaf_lai

        # Split between layers
        remaining_leaf_lai = new_leaf_lai
        for iL in range(nL):
            layer_lai = remaining_layer_lai[iL]

            if layer_lai > 0:
                f_p_in_layer = min(1, max(0, remaining_leaf_lai / layer_lai))
                p_in_layer = layer_lai * f_p_in_layer
                remaining_layer_lai[iL] -= p_in_layer
                remaining_leaf_lai -= p_in_layer
                lai_per_leaf_pop_per_layer[iL, iP] = p_in_layer
    # TODO: Below has been disabled as LAI needs fixing
    # try:
    #     assert isclose(np.sum(lai_per_leaf_pop_per_layer, axis=None), sum(layers_lai), abs_tol=1e-3)
    # except AssertionError as e:
    #     print(
    #         f"lai_per_leaf_pop_per_layer: {np.sum(lai_per_leaf_pop_per_layer, axis=None), sum(layers_lai)}")
    #     raise e
    return lai_per_leaf_pop_per_layer


# def get_growing_populations_range(
#         nP: int,
#         td: List[float],
#         leaf_population_t_emerg: List[float],
#         leaf_population_t_lems: List[float],
# ) -> List[bool]:
#     """Calculate which leaf populations are still growing"""
#     growing_populations = [False for _ in range(nP)]
#     for iP in range(nP):
#         t_emerg = leaf_population_t_emerg[iP]
#         t_lem = leaf_population_t_lems[iP]
#         growing_populations[iP] = t_emerg < td <= t_emerg + t_lem
#     return growing_populations


def get_growing_populations_range_from_config(
    species_config: SpeciesConfig,
    nP: int,
    td: List[float],
) -> Tuple[np.ndarray, np.ndarray]:
    """Get an array representing which leaf populations are growing.

    Parameters
    ----------
    species_config : SpeciesConfig
        Species Config
    nP : int
        Number of populations
    td : List[float]
        Thermal Time data from sowing

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        0: Growing populations with shape (nrows, nP)
        1: Emerged leaf populations shape (nrows)

    """
    plant_emerge_to_flag_emerg = species_config.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg
    sowing_to_emerg_td = species_config.key_lengths_td.sowing_to_emerge
    td_at_sowing = species_config.key_dates_td.sowing
    # last_t_emerg = plant_emerge_to_flag_emerg + td_at_sowing + sowing_to_emerg_td
    t_lem = species_config.key_lengths_leaf_td.tl_em
    flag_t_lem = species_config.key_lengths_flag_leaf_td.tl_em

    td_dd_plant_emerge = td - td_at_sowing - sowing_to_emerg_td
    # td_dd_plant_sowing = td - td_at_sowing

    emergence_rate = calc_emergence_rate(nP, plant_emerge_to_flag_emerg)
    emerged_leaf_populations_count = [
        max(0, min(nP, ceil(t * emergence_rate)))
        for t in td_dd_plant_emerge]

    # assumes each leaf population emerges after the previous
    leaf_population_t_emerg = [i * t_lem for i in range(nP)]
    # Alt method not working
    # leaf_population_t_emerg = np.array(flatten_list([
    #     [t for _ in range(c - emerged_leaf_populations_count[i - 1])] for (i, c), t in zip(enumerate(emerged_leaf_populations_count), td_dd_plant_emerge)
    #     if c >= 1 and c > emerged_leaf_populations_count[i - 1]
    # ])[0:nP - 1] + [last_t_emerg])

    td_emerge_leaf_pops = [[tdd - t_emerge for t_emerge in leaf_population_t_emerg]
                           for tdd in td_dd_plant_emerge]

    leaf_population_t_lems = [t_lem for _ in range(nP - 1)] + [flag_t_lem]

    growing_populations_all = [
        get_growing_populations(td_dd_emerg, leaf_population_t_lems, flag_t_lem) for td_dd_emerg in td_emerge_leaf_pops
    ]

    return growing_populations_all, emerged_leaf_populations_count


def distribute_canopy_lai_to_leaf_pops(
    nL: int,
    nP: int,
    no_emerged_pops: int,
    canopy_lai: float,
) -> List[List[float]]:
    """Distribute lai fractions between leaf populations.

    This is an alternative to distributing using carbon allocation and leaf growth stages.

    Parameters
    ----------
    nL : int
        Total number of layers
    nP : int
        Total number of populations
    no_emerged_pops : int
        The number of populations that have emerged
    canopy_lai: float
        total canopy lai to be distributed between populations

    Returns
    -------
    List[List[float]]
        The lai per leaf population per layer with shape (nL, nP)

    """
    fLAI = np.zeros((nL, nP))
    if no_emerged_pops == 0:
        return fLAI.tolist()
    for iP in range(nP):
        for iL in range(nL):
            if no_emerged_pops >= iP + 1:
                layer_tot = sum(fLAI[iL, 0:iP])
                pop_tot = sum(fLAI[0:nL, iP])
                cap = no_emerged_pops / nL
                flai = min(1 - pop_tot, cap - layer_tot)
                fLAI[iL, iP] = flai
    return (fLAI * canopy_lai / no_emerged_pops).tolist()

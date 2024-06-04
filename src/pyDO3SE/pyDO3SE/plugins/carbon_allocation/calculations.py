"""Carbon Allocation Module.

# References
- Osborne et al 2015 - Jules Crop A parameterization
- Williams et al 2016 - Evaluation of JULES-crop performance against site observations of
irrigated maize from Mead, Nebraska


"""


from math import exp
from collections import namedtuple


def calc_net_prod(
    canopy_An: float,  # net canopy photosynthesis
    c_root: float,  # root carbon pool
    c_stem: float,  # stem carbon pool
    c_leaf: float,  # leaf carbon pool
    R_dc: float,
    r_g: float = 0.25,
    beta: float = 1,
) -> float:
    """Calculating hourly accumulation of Net Primary Productivity.

    "To simulate crop growth, net primary productivity (5) is accumulated over a day
    and then partitioned between five carbon pools: root (Croot), structural stem (Cstem),
    stem reserves (Cresv), leaves (Cleaf), and harvested organs (Charv). The original
    formulation for 5 in JULES includes assumptions about the sizes of the leaf, stem
    and root carbon pools in order to" *Osborne et al 2015*

    References
    ----------

    - Eq4. Osborne et al 2015
    - D. B. Clark et al.: JULES: carbon fluxes and vegetation dynamics


    Alt names
    ---------

    (NPP) (net_prod)

    Parameters
    ----------
    canopy_An : float
        net canopy photosynthesis [kg C m^-2]
    c_root: float
        root carbon pool [kg C m^-2]
    c_stem: float
        stem carbon pool [kg C m^-2]
    c_leaf: float
        leaf carbon pool [kg C m^-2]
    R_dc: float
        rate of non-moisture-stessed canopy dark respiration [kg C m^-2]
    r_g: float
        fraction of gross primary productivity less maintaenance respiration that
        is assigned to growth respiration [fraction]
    beta: float
        soil moisture dark respiration modifier

    Returns
    -------
    float
        Net primary productivity [kg C m^-2 hour^-1]


    Notes
    -----

    ""
    The gross primary productivity (5G) is:
    5G = Ac +βRdc (38)

    where βRdc is the soil moisture-modified canopy dark respi-
    ration. The net primary productivity (5) is:

    5 = 5g −Rp (39)
    where Rp is plant respiration. Rp is split into maintenance
    and growth respiration (Cox et al., 1999):
    Rp = Rpm +Rpg (40)
    Growth respiration is assumed to be a fixed fraction of the
    net primary productivity:
    Rpg = rg
    5G −Rpm
    . (41)
    The growth respiration coefficient (rg) is set to 0.25 for all
    PFTs. Leaf maintenance respiration is equivalent to the soil
    moisture-modified canopy dark respiration, βRdc, while root
    and stem respiration are assumed to be independent of soil

    moisture, but to have the same dependences on nitrogen con-
    tent and temperature. Thus total maintenance respiration is

    given by:eq 42
    "" Clark et al 2011


    Using Osbore method causes issues owing to discrepency in units and naming convention.
    Instead we use above Clark et al 2011 method.

    If $ R_{dc} * \frac{c_{root}+c_{stem}}{c_{leaf}} > A_{canopy} $ then $A_{netprod} $ < 0

    rg=0.25 for all plant functional types in the original documentation for the JULES
    Vegetation Model part 2, it's stated just under eqn 41

    Day Cycle
    ---------

    "Net primary productivity is accumulated over a day and then divided into five
    crop components according to paritition coeff" Osborne et al 2015

    """
    if c_leaf <= 0:
        return 0
    GPP = canopy_An + beta * R_dc

    # TODO: Below carbon ratio should be Nitrogen
    R_pm = R_dc * (beta + ((c_root + c_stem) / c_leaf))
    R_pg = r_g * (GPP - R_pm)
    R_p = R_pm + R_pg

    NPP = GPP - R_p
    return NPP

    # Below Osborne method has unknown units
    # NPP = (1 - r_g) * (canopy_An - R_dc * ((c_root + c_stem) / c_leaf))
    # return NPP


CarbonPoolChange = namedtuple(
    'CarbonPools', 'c_root_diff c_leaf_diff c_stem_diff c_harv_diff c_resv_diff')


def calc_carbon_pool_change(
    net_prod_acc: float,
    p_root: float,
    p_leaf: float,
    p_stem: float,
    p_harv: float,
    theta: float = 0.4,
) -> CarbonPoolChange:
    """Distribution of carbon between pools

    At the end of the day we distribute the carbon between the pools based on partition
    coefficients. These coefficients vary across the growing season as a function of
    thermal time since emergence.

    References
    ----------

    - Eq5. Osborne et al 2015


    Parameters
    ----------
    net_prod_acc: float
        acumulated net productivity [kg C m^-2]
    p_root: float
        partition coefficient for root [fraction]
    p_leaf: float
        partition coefficient for leaf [fraction]
    p_stem: float
        partition coefficient for stem [fraction]
    p_harv: float
        partition coefficient for harv [fraction]
    theta: float
        fraction of stem carbon that is partitioned into reserve pool - crop specific parameter

    Returns
    -------
    CarbonPoolChange
        Change in carbon for each pool
        c_root_diff c_leaf_diff c_stem_diff c_harv_diff c_resv_diff

    """
    c_root_diff = p_root * net_prod_acc
    c_leaf_diff = p_leaf * net_prod_acc
    c_stem_diff = p_stem * net_prod_acc * (1 - theta)
    c_harv_diff = p_harv * net_prod_acc
    c_resv_diff = p_stem * net_prod_acc * theta
    return CarbonPoolChange(
        c_root_diff,
        c_leaf_diff,
        c_stem_diff,
        c_harv_diff,
        c_resv_diff,
    )


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


AdjustedCarbonPools = namedtuple('AdjustedCarbonPools', 'c_leaf c_harv c_resv')


def adjust_c_pools_at_eol(
    DVI: float,
    c_leaf: float,
    p_stem: float,
    c_harv: float,
    c_resv: float,
) -> AdjustedCarbonPools:
    """Adjust the carbon pools at the end of the plants life.

    "...once carbon is no longer partitioned to stems, carbon from the
    stem reserve pool is mobilised to the harvest pool at a rate
    of 10 % a day" Osborne et al 2015

    "Leaf senescence is treated simplistically by mobilising carbon from the leaf to
    the harvest pool at a rate of 0.05 d−1 once DVI has reached 1." Osborne et al 2015


    References
    ----------

    - Eq7. Osborne et al 2015
    - Eq8. Osborne et al 2015



    Returns
    -------
    AdjustedCarbonPools[NamedTuple]
        Adjusted carbon pools
    """
    c_harv_out = c_harv
    c_resv_out = c_resv
    c_leaf_out = c_leaf
    # EQ 7
    if p_stem < 0.01:
        c_harv_out = c_harv + 0.1 * c_resv
        c_resv_out = 0.9 * c_resv

    # EQ 8
    if DVI > 1.5:
        c_harv_out = c_harv_out + 0.05 * c_leaf
        c_leaf_out = 0.95 * c_leaf

    return AdjustedCarbonPools(
        c_leaf=c_leaf_out,
        c_harv=c_harv_out,
        c_resv=c_resv_out,
    )


def calc_LAI_from_DVI_and_carbon(
    DVI: float,
    c_leaf: float,
    gamma: float,
    delta: float,
    emerged_leaf_count: int,
    f_c: float = 0.5,  # Const
) -> float:
    """Carbon pool to Leaf Area Index.

    At the end of each growth time step (24 h), the amount of
    carbon in the leaves is related to leaf area index" Osborne et al 2015.


    References
    ----------

    - Eq9. Osborne et al 2015
    - Eq10. Osborne et al 2015


    Parameters
    ----------
    DVI : float
        Development Index [unitless][-1 -> 2]
    c_leaf : float
        Leaf carbon pool [kg C m^-2]
    gamma : float
        crop specific parameter
    delta : float
        crop specific parameter
    emerged_leaf_count: int,
        number of emerged leaves
    f_c : float, optional
        constant, by default 0.5

    Returns
    -------
    float
        Leaf area index(LAI)[m^2/m^2]

    """
    # if DVI <= -0.06:
    #     return 0
    if DVI <= 0.0:
        return 0
    # TODO: Below should not occur when c_leaf > 0
    if emerged_leaf_count <= 0:
        return 0
    SLA = gamma * (DVI + 0.06) ** delta
    LAI = (c_leaf / f_c) * SLA
    return LAI


def get_plant_height_from_carbon(
    c_stem: float,
    k: float,
    lambdav: float = 0.4,
    f_c: float = 0.5,  # Const
) -> float:
    """Carbon to plant height.

    "The amount of carbon in the stem is related to the crop height by (Hunt, 1990)"


    References
    ----------

    - Eq11. Osborne et al 2015


    Parameters
    ----------
    c_stem : float
        stem carbon pool [kg C m^-2]
    k : float
        Crop specific coefficient
    lambdav : float
        Crop specific coefficient
    f_c: float
        Constant value, default = 0.5

    Returns
    -------
    float
        plant height[m]
    """
    height = k * (c_stem / f_c) ** lambdav
    return height


def calc_root_fraction_from_carbon(
    c_root: float,
    d_r: float = 0.5,  # Param - for all crop types
    r_dir: float = 0.0,  # Param - Crop specific table 4(Actually 0.0 for all crops!)
    z: float = 0.5,  # Root depth
    f_c: float = 0.5,
) -> float:
    """Carbon to root fraction.

    NOTE: Not used in DO3SE.

    "Because root biomass increases during the crop growing
    season the fraction of roots in each JULES soil layer varies
    according to the equation of Arora and Boer (2003) which
    defines the fraction of roots at depth z a" Osborne et al


    # TODO: Check r_dir is meant to be 0.0 for all params

    References
    ----------

    - Eq12. Osborne et al 2015
    - Eq13. Osborne et al 2015


    Parameters
    ----------
    c_root : float
        root carbon pool [kg C m^-2]
    d_r : float, optional
        [description], by default 0.5
    r_dir: float = 0.0,
        Param - Crop specific table 4(Actually 0.0 for all crops!)
    z: float,
        Root depth
    f_c: float,
        const value, default = 0.5
    Returns
    -------
    float
        Fraction of root??
    """
    a = d_r * (c_root / f_c) ** r_dir
    f = 1 - exp(-z / a)
    return f


CarbonPools = namedtuple('CarbonPools', 'c_root c_stem c_leaf c_harv c_resv')


def daily_carbon_allocation(
    net_prod_acc: float,
    DVI: float,
    c_root: float,
    c_stem: float,
    c_leaf: float,
    c_harv: float,
    c_resv: float,
    a_root: float,
    a_leaf: float,
    a_stem: float,
    b_root: float,
    b_leaf: float,
    b_stem: float,
    theta: float,
    c_init: float = 8e-4,
) -> CarbonPools:
    """Calculate the carbon pools at the end of a days accumulation.

    "The crop carbon pools are initialised at DVIinit, which is
    at or just after emergence. At initialisation, the crops are
    given a certain amount of carbon Cinit, which is distributed
    between the carbon pools according to the values of pi at
    DVI = DVIinit."  Williams et al 2016

    Parameters
    ----------
    net_prod_acc : float
        acumulated net productivity [kg C m^-2]
    DVI : float
        Development Index [unitless][-1 -> 2]
    c_root : float
        root carbon pool [kg C m^-2]
    c_stem : float
        stem carbon pool [kg C m^-2]
    c_leaf : float
        leaf carbon pool [kg C m^-2]
    c_harv : float
        harv carbon pool [kg C m^-2]
    c_resv : float
        resv carbon pool [kg C m^-2]
    a_root : float
        partition fraction parameter [unitless]
    a_leaf : float
        partition fraction parameter [unitless]
    a_stem : float
        partition fraction parameter [unitless]
    b_root : float
        partition fraction parameter [unitless]
    b_leaf : float
        partition fraction parameter [unitless]
    b_stem : float
        partition fraction parameter [unitless]
    theta: float
        fraction of stem carbon that is partitioned into reserve pool [fraction]
        crop specific parameter
    c_init: float
        amount of carbon at DVI_init(Emergence)

    Returns
    -------
    CarbonPools [namedtuple]
        c_root c_stem c_leaf c_harv c_resv

    """
    # Note we set DVI to be min 0 here to ensure correct emergence day values
    (p_root, p_leaf, p_stem, p_harv) = calc_partition_coefficients(
        DVI if DVI > 0 else 0,  # -1 -> 2
        a_root,
        a_leaf,
        a_stem,
        b_root,
        b_leaf,
        b_stem,
    )

    (
        c_root_diff,
        c_leaf_diff,
        c_stem_diff,
        c_harv_diff,
        c_resv_diff
    ) = calc_carbon_pool_change(
        net_prod_acc if DVI > 0 else c_init,  # acumulated net productivity
        p_root,  # partition coefficient for root
        p_leaf,  # partition coefficient for leaf
        p_stem,  # partition coefficient for stem
        p_harv,  # partition coefficient for harv
        theta,
    )

    if DVI > 2:
        # Stop accumulation after harvest
        c_root_out = c_root
        c_stem_out = c_stem
        c_leaf_out = c_leaf
        c_harv_out = c_harv
        c_resv_out = c_resv
    if DVI > 0:
        # After emergence increase carbon
        c_root_out = max(0, c_root + c_root_diff)
        c_stem_out = max(0, c_stem + c_stem_diff)
        c_leaf_out = max(0, c_leaf + c_leaf_diff)
        c_harv_out = max(0, c_harv + c_harv_diff)
        c_resv_out = max(0, c_resv + c_resv_diff)
    elif DVI > -1:
        # Before emergence set init carbon
        c_root_out = c_root_diff
        c_stem_out = c_stem_diff
        c_leaf_out = c_leaf_diff
        c_harv_out = c_harv_diff
        c_resv_out = c_resv_diff
    else:
        # Before sowing set to 0
        c_root_out = 0
        c_stem_out = 0
        c_leaf_out = 0
        c_harv_out = 0
        c_resv_out = 0

    if p_stem < 0.01:
        c_harv_out = c_harv_out + 0.1 * c_resv_out
        c_resv_out = 0.9 * c_resv_out

    if DVI > 1.5:
        c_harv_out = c_harv_out + 0.05 * c_leaf_out
        c_leaf_out = 0.95 * c_leaf_out

    return CarbonPools(
        c_root=c_root_out,
        c_stem=c_stem_out,
        c_leaf=c_leaf_out,
        c_harv=c_harv_out,
        c_resv=c_resv_out,
    )

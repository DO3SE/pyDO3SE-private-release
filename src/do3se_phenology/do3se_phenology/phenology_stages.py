"""Functions related to the leaf and plant phenology stages."""

from do3se_phenology.state import LeafPhenologyStage, PhenologyStage
from do3se_phenology.config import PhenologyKeyLengths, PhenologyKeyDates
from do3se_phenology.utils import wrap_day_of_year


def get_plant_phenology_stage_td(
    td_dd: float,
    key_lengths: PhenologyKeyLengths,
    phenology_stage: PhenologyStage,
) -> PhenologyStage:
    """Get the plant phenology stage.

    Parameters
    ----------
    td_dd : float
        thermal time since sowing
    key_lengths : PhenologyKeyLengths
        key phenology lengths
    phenology_stage: PhenologyStage
        Current phenology stage

    Returns
    -------
    PhenologyStage
        Output phenology stage

    """
    if phenology_stage < PhenologyStage.EMERGED:
        return phenology_stage
    elif td_dd <= key_lengths.emerg_to_astart:
        return PhenologyStage.EMERGED
    elif td_dd <= key_lengths.emerg_to_end:
        return PhenologyStage.ASTART
    else:
        return PhenologyStage.HARVEST


def get_plant_phenology_stage(
    dd: int,
    key_dates_sowing: PhenologyKeyDates.sowing,
    key_dates_emergence: PhenologyKeyDates.emergence,
    key_dates_harvest: PhenologyKeyDates.harvest,
    phenology_stage: PhenologyStage,
) -> PhenologyStage:
    """Get the plant phenology stage.

    Parameters
    ----------
    dd : int
        day of year
    key_dates_sowing: PhenologyKeyDates.sowing,
    key_dates_emergence: PhenologyKeyDates.emergence,
    key_dates_harvest: PhenologyKeyDates.harvest,
    phenology_stage: PhenologyStage
        Current phenology stage

    Returns
    -------
    PhenologyStage
        Output phenology stage

    """
    dd_wrapped = wrap_day_of_year(dd - key_dates_sowing)
    if phenology_stage < PhenologyStage.EMERGED:
        # Emergence handled elsewhere
        return phenology_stage
    # NOTE Currently disabled as we don't always have emerge
    # elif dd <= key_lengths.emerg_to_astart:
    #     return PhenologyStage.EMERGED
    # elif dd <= key_lengths.emerg_to_end:
    #     return PhenologyStage.ASTART
    elif key_dates_emergence is not None and dd_wrapped <= (
        key_dates_emergence - key_dates_sowing
    ):
        return PhenologyStage.SOWN
    elif dd_wrapped <= wrap_day_of_year(key_dates_harvest - key_dates_sowing):
        # This doesn't really do anything as we have already checked for emergence
        return PhenologyStage.EMERGED
    else:
        return PhenologyStage.HARVEST


def get_leaf_phenology_stage_td(
    td_dd: float,
    t_lem: float,
    t_lep: float,
    t_lse: float,
) -> LeafPhenologyStage:
    if td_dd <= 0:
        return LeafPhenologyStage.NOT_EMERGED
    elif td_dd <= t_lem:
        return LeafPhenologyStage.GROWING
    elif td_dd <= t_lem + t_lep:
        return LeafPhenologyStage.MATURE
    elif td_dd <= t_lem + t_lep + t_lse:
        return LeafPhenologyStage.SENESCENCE
    else:
        return LeafPhenologyStage.FULLY_SENESED


def get_leaf_phenology_stage(
    dd: int,
    dd_at_sowing: int,
    sowing_to_plant_emerg: int,
    plant_emerg_to_leaf_emerg: int,
    leaf_emerg_to_fully_grown: int,
    fully_grown_to_senescence: int,
    sowing_to_end: int,
) -> LeafPhenologyStage:
    dd_diff = wrap_day_of_year(dd - dd_at_sowing)
    if dd_diff <= sowing_to_plant_emerg + plant_emerg_to_leaf_emerg:
        return LeafPhenologyStage.NOT_EMERGED
    elif (
        dd_diff
        <= sowing_to_plant_emerg + plant_emerg_to_leaf_emerg + leaf_emerg_to_fully_grown
    ):
        return LeafPhenologyStage.GROWING
    elif (
        dd_diff
        <= sowing_to_plant_emerg
        + plant_emerg_to_leaf_emerg
        + leaf_emerg_to_fully_grown
        + fully_grown_to_senescence
    ):
        return LeafPhenologyStage.MATURE
    elif dd_diff <= sowing_to_end:
        return LeafPhenologyStage.SENESCENCE
    else:
        return LeafPhenologyStage.FULLY_SENESED

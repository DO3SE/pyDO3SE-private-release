"""A module for processing various phenology configurations into a consistent setup.

This is in response to the many methods of describing phenology that results in additional confusion in
the DO3SE model.


    JulianDay = int
    ThermalTime = float
    XUnit = Union[JulianDay, ThermalTime]
    YUnit = float
    PiecewiseFunction = List[Tuple[XUnit, YUnit]]

You will need different parameter setups depending on the type of input data you have. Below are examples of possible input data:

- SEASON_FRACTION - Sowing Day[Julian day] and season length(f_phen_d)[Thermal time]
    - All phenology values are set as a fraction of the season length with thermal time starting at sowing day
- FLAG_LEAF_FRACTION - Astart[Julian day] and season length(f_phen_d)[Thermal time]
    - All phenology values are set as a fraction of the season length with thermal time starting at Astart
- FPHEN_THERMAL_TIME - t_fphen values are all specified
    - All phenology values are set as input constants. Missing values are calculated as a % of season length
- LEAF_FPHEN_DATA - leaf_f_phen data supplied
    - Astart, Senesence start and Aend are calculated from data. All other phenology values are calculated as a %
...


"""

from copy import deepcopy
from typing import Callable, NamedTuple, Tuple, List
from dataclasses import replace
import warnings
import numpy as np

from data_helpers.dictionary_helpers import merge_dataclasses
from thermal_time.calcs import calc_thermal_time_range
from do3se_phenology.latitude_function import (
    lat_function,
    lat_function_spring_wheat_europe,
    lat_function_winter_wheat_china,
    lat_function_forest_europe,
)

from do3se_phenology.units import *
from do3se_phenology.phyllochron_dvi import get_dvi_PLF
from do3se_phenology.td_percent_definition import (
    calculate_growing_season_from_leaf_f_phen_data,
    get_canopy_td_intervals_f,
    get_dvi_from_season_length,
    get_leaf_td_intervals_f,
)
from do3se_phenology.f_phen import get_fphen_PLF, get_leaf_fphen_PLF
from do3se_phenology.config import (
    ModelConfig,
    PhenologyLeafKeyLengths,
    PhenologyMethods,
    LeafFPhenMethods,
    FPhenMethods,
    SowingDateMethods,
    SpeciesConfig,
    PhenologyKeyDates,
    SpeciesPresetsParams,
    ZeroDayOptions,
)
from do3se_phenology.utils import get_day_from_td, get_td_from_day

from .error_handling import ConfigError
from .units import ThermalTime, JulianDay


def phenology_from_legacy_day_plf(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
) -> Tuple[ModelConfig, SpeciesConfig]:
    sowing_day = get_sowing_day_from_config(
        species_config,
        model_config,
    )
    # 2. Get season length
    season_length = get_season_length_from_config(
        species_config,
        model_config,
        sowing_day,
    )
    assert season_length, "season_length could not be defined!"
    species_config.key_lengths.sowing_to_end = season_length

    assert model_config.time_type == TimeTypes.JULIAN_DAY, "Must set phenology_options.time_type to TimeTypes.JULIAN_DAY to use PhenologyMethods.LEGACY_DAY_PLF"

    sowing_to_emerg = species_config.key_lengths.sowing_to_emerge or 0
    egs = sowing_day + season_length

    if species_config.leaf_f_phen_method in [
            LeafFPhenMethods.DAY_PLF, LeafFPhenMethods.TT_DAY_PLF, LeafFPhenMethods.TT_GROWING_SEASON]:

        sowing_to_astart = species_config.key_lengths.sowing_to_astart
        assert sowing_to_astart, "sowing_to_astart could not be defined!"
        Astart = sowing_day + sowing_to_astart
        emerg_to_astart = Astart - sowing_to_emerg

        flag_emerg_to_astart = 0  # NOTE: This method assumes flag emergence time is 0
        plant_emerg_to_flag_emerg = Astart - sowing_to_emerg - flag_emerg_to_astart
        astart_to_senescence = egs - sowing_to_astart

        species_config.key_lengths_flag_leaf.plant_emerg_to_leaf_emerg = plant_emerg_to_flag_emerg
        species_config.key_lengths_flag_leaf.leaf_emerg_to_astart = flag_emerg_to_astart
        species_config.key_lengths_flag_leaf.astart_to_senescence = astart_to_senescence
        species_config.key_lengths.emerg_to_astart = emerg_to_astart
        species_config.key_dates.Astart = Astart
        species_config.key_dates.Aend = egs

        assert species_config.key_lengths_flag_leaf.plant_emerg_to_leaf_emerg is not None, "key_lengths_flag_leaf.plant_emerg_to_leaf_emerg could not be defined!"
        assert species_config.key_lengths_flag_leaf.leaf_emerg_to_astart is not None, "key_lengths_flag_leaf.leaf_emerg_to_astart could not be defined!"
        assert species_config.key_lengths_flag_leaf.astart_to_senescence is not None, "key_lengths_flag_leaf.astart_to_senescence could not be defined!"
        assert species_config.key_lengths.emerg_to_astart is not None, "key_lengths.emerg_to_astart could not be defined!"
        assert species_config.key_dates.Astart, "AStart day could not be defined!"
        assert species_config.key_dates.Aend, "AEnd day could not be defined!"
    else:
        species_config.key_lengths_flag_leaf.plant_emerg_to_leaf_emerg = 0

    emerg_to_end = season_length - sowing_to_emerg

    species_config.key_dates.sowing = sowing_day
    species_config.key_dates.harvest = egs

    species_config.key_lengths.emerg_to_end = emerg_to_end

    species_config.key_dates.sowing = sowing_day
    species_config.key_dates.harvest = egs

    assert species_config.key_lengths.emerg_to_end is not None, "key_lengths.emerg_to_end could not be defined!"

    assert species_config.key_dates.harvest is not None, "harvest day could not be defined!"
    assert species_config.key_dates.sowing is not None, "Sowing day could not be defined!"

    if species_config.f_phen_method == FPhenMethods.SIMPLE_DAY_PLF:
        assert species_config.day_fphen_plf.f_phen_1 is not None, "day_fphen_plf.f_phen_1 must be set in phenology species config."
        assert species_config.day_fphen_plf.f_phen_4 is not None, "day_fphen_plf.f_phen_4 must be set in phenology species config."
        assert species_config.day_fphen_plf.f_phen_a is not None, "day_fphen_plf.f_phen_a must be set in phenology species config."
        assert species_config.day_fphen_plf.f_phen_c is not None, "day_fphen_plf.f_phen_c must be set in phenology species config."
        assert species_config.day_fphen_plf.f_phen_e is not None, "day_fphen_plf.f_phen_e must be set in phenology species config."
        assert species_config.key_dates.sowing is not None, "key_dates.sowing must be set in phenology species config."
        assert species_config.key_dates.harvest is not None, "key_dates.harvest must be set in phenology species config."

    return model_config, species_config


def phenology_from_leaf_f_phen_data(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
    leaf_f_phen_data: List[float],
    td_data: List[ThermalTime],
    dd_data: List[JulianDay],
) -> Tuple[ModelConfig, SpeciesConfig]:
    """Generate all phenology parameters from leaf fphen data."""

    sowing_offset = 0
    # Thermal Time values
    t_sgs, t_Astart, t_egs, t_mature_leaf, t_senes_start = calculate_growing_season_from_leaf_f_phen_data(
        leaf_f_phen_data,
        td_data,
        f_leaf_f_fphen=species_config.f_fphen_1_ets + species_config.f_fphen_5_ets,
    )

    # Julian day values
    SGS, Astart, EGS, mature_leaf, senes_start = calculate_growing_season_from_leaf_f_phen_data(
        leaf_f_phen_data,
        dd_data,
        f_leaf_f_fphen=species_config.f_fphen_1_ets + species_config.f_fphen_5_ets,
    )

    canopy_intervals = get_canopy_td_intervals_f(
        t_sgs=t_sgs,
        t_egs=t_egs,
        f_Astart=species_config.f_Astart,
        f_mid_anthesis=species_config.f_mid_anthesis,
        f_fphen_a=species_config.f_fphen_a,
        f_fphen_b=species_config.f_fphen_b,
        f_fphen_c=species_config.f_fphen_c,
        f_fphen_d=species_config.f_fphen_d,
        f_tt_emr=species_config.f_tt_emr,
        f_tt_veg=species_config.f_tt_veg,
        f_tt_rep=species_config.f_tt_rep,
    )

    leaf_intervals = get_leaf_td_intervals_f(
        t_sgs=t_sgs,
        t_egs=t_egs,
        f_Astart=species_config.f_Astart,
        f_mid_anthesis=species_config.f_mid_anthesis,
        f_fphen_1_ets=species_config.f_fphen_1_ets,
        f_fphen_3_ets=species_config.f_fphen_3_ets,
        f_fphen_4_ets=species_config.f_fphen_4_ets,
        f_fphen_5_ets=species_config.f_fphen_5_ets,
        f_t_lem=species_config.f_t_lem,
        f_t_lma=species_config.f_t_lma,
        f_t_lep=species_config.f_t_lep,
        f_t_lse=species_config.f_t_lse,
    )

    tt_emr, tt_veg, tt_rep = get_dvi_from_season_length(
        t_sgs,
        t_egs,
        species_config.f_tt_emr,
        species_config.f_tt_veg,
        species_config.f_tt_rep,
    )
    dvi_interval = list(
        zip(*get_dvi_PLF(tt_emr, tt_veg, tt_rep, x_offset=-canopy_intervals.t_Astart)))

    # 1. Get sowing offset
    if model_config.zero_day == ZeroDayOptions.SOWING:
        raise NotImplementedError("Cannot use SOWING zero day for leaf f phen data method")
    elif model_config.zero_day == ZeroDayOptions.ASTART:
        sowing_offset = -canopy_intervals.t_Astart
    else:
        raise NotImplementedError(f"zero day method {model_config.zero_day} not implemented")

    key_dates_td = PhenologyKeyDates(
        sowing=sowing_offset,  # All dates relative to sowing day
        emergence=sowing_offset + canopy_intervals.t_fphen_a,
        harvest=sowing_offset + canopy_intervals.t_fphen_d,
        Astart=sowing_offset + canopy_intervals.t_Astart,
        Aend=sowing_offset + canopy_intervals.t_fphen_d,
        mid_anthesis=sowing_offset + leaf_intervals.t_Astart + leaf_intervals.t_fphen_1_ets,
    )
    key_dates = PhenologyKeyDates(
        Astart=Astart,
    )

    # We need to get leaf_f_phen_g and leaf_f_phen_h
    t_mid_anthesis = t_Astart + leaf_intervals.t_fphen_1_ets
    mid_anthesis_to_mature = t_mature_leaf - t_mid_anthesis
    mid_anthesis_to_senes = t_senes_start - t_mid_anthesis
    emerg_to_astart = canopy_intervals.t_Astart - canopy_intervals.t_fphen_a
    emerg_to_flag_emerg = emerg_to_astart - leaf_intervals.t_lem

    key_lengths_flag_leaf_td = PhenologyLeafKeyLengths(
        leaf_f_phen_e=leaf_intervals.t_fphen_1_ets,  # 0.08
        leaf_f_phen_g=mid_anthesis_to_mature,
        leaf_f_phen_h=mid_anthesis_to_senes,
        leaf_f_phen_i=leaf_intervals.t_fphen_5_ets,  # 0.38
        tl=leaf_intervals.t_l,
        tl_em=leaf_intervals.t_lem,
        tl_ma=leaf_intervals.t_lma,
        tl_ep=leaf_intervals.t_lep,
        tl_se=leaf_intervals.t_lse,
        plant_emerg_to_leaf_emerg=emerg_to_flag_emerg,
    )
    key_lengths_td = replace(
        species_config.key_lengths_td,
        sowing_to_emerge=canopy_intervals.t_fphen_a,
        sowing_to_f_phen_b=canopy_intervals.t_fphen_b,
        sowing_to_f_phen_c=canopy_intervals.t_fphen_c,
        sowing_to_astart=canopy_intervals.t_Astart,
        sowing_to_end=canopy_intervals.t_fphen_d,
        emerg_to_astart=canopy_intervals.t_Astart - canopy_intervals.t_fphen_a,
        emerg_to_end=canopy_intervals.t_fphen_d - canopy_intervals.t_fphen_a,
        emerg_to_veg=tt_veg,
        veg_to_harvest=tt_rep,
    )

    model_config_out = replace(
        model_config,
    )

    species_config_out = replace(
        species_config,
        fphen_intervals=None,
        leaf_fphen_intervals=None,
        dvi_interval=dvi_interval,
        key_lengths_td=key_lengths_td,
        key_lengths_flag_leaf_td=key_lengths_flag_leaf_td,
        key_dates_td=key_dates_td,
        key_dates=key_dates,
    )

    return model_config_out, species_config_out


def get_sowing_day_from_config(
    species_config: SpeciesConfig,
    model_config: ModelConfig,
) -> int:
    sowing_day = None
    # Get sowing day
    try:
        if model_config.sowing_day_method == SowingDateMethods.LATITUDE:
            if model_config.latitude is None:
                raise ConfigError(
                    "Must supply latitude in phenology model config to use SowingDateMethods.LATITUDE")
            if species_config.lat_f_k is None:
                raise ConfigError(
                    "Must supply species_config.lat_f_k to use SowingDateMethods.LATITUDE")
            if species_config.lat_f_b is None:
                raise ConfigError(
                    "Must supply species_config.lat_f_b to use SowingDateMethods.LATITUDE")
            if species_config.lat_f_c is None:
                raise ConfigError(
                    "Must supply species_config.lat_f_c to use SowingDateMethods.LATITUDE")
            sowing_day = int(lat_function(model_config.latitude, species_config.lat_f_k, species_config.lat_f_b, species_config.lat_f_c))
        elif model_config.sowing_day_method == SowingDateMethods.LATITUDE_FOREST:
            if model_config.latitude is None:
                raise ConfigError(
                    "Must supply latitude in phenology model config to use SowingDateMethods.LATITUDE_FOREST")
            sowing_day = int(lat_function_forest_europe(model_config.latitude))
        elif model_config.sowing_day_method == SowingDateMethods.LATITUDE_SPRING_EUROPE:
            if model_config.latitude is None:
                raise ConfigError(
                    "Must supply latitude in phenology model config to use SowingDateMethods.LATITUDE_SPRING_EUROPE")
            sowing_day = int(lat_function_spring_wheat_europe(model_config.latitude))
        elif model_config.sowing_day_method == SowingDateMethods.LATITUDE_WINTER_CHINA:
            if model_config.latitude is None:
                raise ConfigError(
                    "Must supply latitude in phenology model config to use SowingDateMethods.LATITUDE_WINTER_CHINA")
            sowing_day = int(lat_function_winter_wheat_china(model_config.latitude))
        elif model_config.sowing_day_method == SowingDateMethods.INPUT:
            if species_config.key_dates.sowing is None:
                raise ConfigError(
                    "Must supply sowing day in species phenology config to use SowingDateMethods.INPUT or set to SowingDateMethods.SKIP")
            sowing_day = int(species_config.key_dates.sowing)
        elif model_config.sowing_day_method == SowingDateMethods.SKIP:
            pass
        else:
            raise ConfigError(f"Unknown SowingDateMethod: {model_config.sowing_day_method}")
    except ConfigError as e:
        raise ConfigError(f"""Could not get sowing day from config. Check previous error for details:
        {str(e)}
        If using flag only runs set phenology_config.sowing_day_method to "SKIP".
        """)

    return sowing_day


def get_td_at_sowing_from_config(
    species_config: SpeciesConfig,
    model_config: ModelConfig,
) -> float:
    # 1. Get sowing offset
    sowing_offset = 0  # offset between sowing and zero day
    if model_config.zero_day == ZeroDayOptions.SOWING:
        sowing_offset = 0
    elif model_config.zero_day == ZeroDayOptions.ASTART:
        # Get astart sowing offset
        if species_config.key_dates_td.sowing is not None:
            sowing_offset = species_config.key_dates_td.sowing
        elif species_config.key_lengths_td.sowing_to_astart is not None:
            sowing_offset = -species_config.key_lengths_td.sowing_to_astart
        else:
            raise NotImplementedError("Must set species_config.key_dates_td.Astart")
    elif model_config.zero_day == ZeroDayOptions.DATA_START:
        warnings.warn("Using untested zero day method: DATA_START")
        sowing_offset = 0
    else:
        raise NotImplementedError(f"zero day method {model_config.zero_day} not implemented")
    return sowing_offset


def get_season_length_from_config(
    species_config: SpeciesConfig,
    model_config: ModelConfig,
    sowing_day: int,
    external_data: dict = None,
    td_base_temperature: float = 0,
) -> float:
    season_length = None
    # get season length
    if model_config.sowing_day_method == SowingDateMethods.LATITUDE_FOREST:
        sowing_day = int(lat_function_forest_europe(model_config.latitude))
        egs = lat_function_forest_europe(model_config.latitude, k=-2, c=297)
        season_length = egs - sowing_day

    elif model_config.time_type == TimeTypes.JULIAN_DAY:
        if species_config.key_lengths.sowing_to_end:
            season_length = species_config.key_lengths.sowing_to_end
        elif species_config.key_lengths.sowing_to_end:
            season_length = species_config.key_lengths.sowing_to_end
        elif species_config.key_lengths.sowing_to_emerge and species_config.key_lengths.emerg_to_end:
            season_length = species_config.key_lengths.sowing_to_emerge + \
                species_config.key_lengths.emerg_to_end
        elif species_config.key_dates.harvest:
            season_length = species_config.key_dates.harvest
    else:
        if species_config.key_lengths_td.sowing_to_end:
            season_length = species_config.key_lengths_td.sowing_to_end
        elif species_config.key_lengths_td.sowing_to_end:
            season_length = species_config.key_lengths_td.sowing_to_end
        elif species_config.key_lengths_td.sowing_to_emerge and species_config.key_lengths_td.emerg_to_end:
            season_length = species_config.key_lengths_td.sowing_to_emerge + \
                species_config.key_lengths_td.emerg_to_end
        elif species_config.key_dates_td.harvest:
            season_length = species_config.key_dates_td.harvest
        elif species_config.key_dates.harvest:
            td_data = calc_thermal_time_range(external_data.get('Ts_C'), t_b=td_base_temperature)
            EGS = species_config.key_dates.harvest
            try:
                t_sgs = next(t for t, d in zip(td_data, external_data.get('dd')) if d == sowing_day)
                t_egs = next(t for t, d in zip(td_data, external_data.get('dd')) if d == EGS)
            except StopIteration:
                raise StopIteration(
                    f'Sowing day or harvest date are outside range of input data. Sowing day: {sowing_day}, Harvest day: {EGS}')
            season_length = t_egs - t_sgs
    if season_length is None:
        raise ConfigError("Could not establish season length from current configuration.")
    return season_length


def get_astart_from_config(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
) -> float:
    Astart = None
    if species_config.key_dates.Astart is not None:
        Astart = species_config.key_dates.Astart
    elif species_config.key_dates_td.Astart is not None:
        raise NotImplementedError("Cannot currently set astart day as thermal time")
    else:
        raise ValueError("Astart must be set in config for flag_leaf_fraction method")
    assert Astart is not None, "Astart could not be defined. Check phenology config."
    return Astart


def phenology_from_season_fraction(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
    external_data: dict,
    td_base_temperature: float,
) -> Tuple[ModelConfig, SpeciesConfig]:
    """Generate the phenology parameters from growing season length fractions.

    To use this method you must supply parameters that can be used to derive the
    length of the growing season in thermal time.

    Parameters
    ----------
    model_config : ModelConfig
        phenology model config
    species_config : SpeciesConfig
        species config
    external_data: dict
        dict of external data
    td_base_temperature: float
        Base temp for thermal time. Only required for julian day season definition

    Returns
    -------
    Tuple[ModelConfig, SpeciesConfig]
        Updated comprehensive model config and species config

    """
    model_config_out = deepcopy(model_config)
    species_config_out = deepcopy(species_config)

    sowing_day = int(get_sowing_day_from_config(
        species_config,
        model_config,
    ))

    # 1. Get sowing offset
    sowing_offset_td = get_td_at_sowing_from_config(
        species_config, model_config)  # offset between sowing and zero day

    if model_config.zero_day == ZeroDayOptions.ASTART:
        Astart = get_astart_from_config(model_config, species_config)

    # 2. Get season length
    season_length = get_season_length_from_config(
        species_config,
        model_config,
        sowing_day,
        external_data,
        td_base_temperature,
    )

    # 3. Get other values
    canopy_td = get_canopy_td_intervals_f(
        0,
        season_length,
        species_config.f_Astart,
        species_config.f_mid_anthesis,
        species_config.f_fphen_a,
        species_config.f_fphen_b,
        species_config.f_fphen_c,
        species_config.f_fphen_d,
        species_config.f_tt_emr,
        species_config.f_tt_veg,
        species_config.f_tt_rep,
    )

    leaf_intervals = get_leaf_td_intervals_f(
        t_sgs=0,
        t_egs=season_length,
        f_Astart=species_config.f_Astart,
        f_mid_anthesis=species_config.f_mid_anthesis,
        f_fphen_1_ets=species_config.f_fphen_1_ets,
        f_fphen_3_ets=species_config.f_fphen_3_ets,
        f_fphen_4_ets=species_config.f_fphen_4_ets,
        f_fphen_5_ets=species_config.f_fphen_5_ets,
        f_t_lem=species_config.f_t_lem,
        f_t_lma=species_config.f_t_lma,
        f_t_lep=species_config.f_t_lep,
        f_t_lse=species_config.f_t_lse,
    )
    # SET PARAMS
    species_config_out.key_dates.sowing = sowing_day

    species_config_out.key_dates_td.sowing = sowing_offset_td
    species_config_out.key_dates_td.emergence = sowing_offset_td + canopy_td.t_fphen_a
    species_config_out.key_dates_td.Astart = sowing_offset_td + canopy_td.t_Astart
    species_config_out.key_dates_td.mid_anthesis = sowing_offset_td + canopy_td.t_mid_anthesis
    species_config_out.key_dates_td.harvest = sowing_offset_td + canopy_td.t_fphen_d
    species_config_out.key_dates_td.Aend = sowing_offset_td + canopy_td.t_fphen_d

    species_config_out.key_lengths_td.sowing_to_emerge = canopy_td.t_fphen_a
    species_config_out.key_lengths_td.sowing_to_f_phen_b = canopy_td.t_fphen_b
    species_config_out.key_lengths_td.sowing_to_f_phen_c = canopy_td.t_fphen_c
    species_config_out.key_lengths_td.sowing_to_astart = canopy_td.t_Astart
    species_config_out.key_lengths_td.sowing_to_end = canopy_td.t_fphen_d
    species_config_out.key_lengths_td.emerg_to_astart = canopy_td.t_Astart - canopy_td.t_fphen_a
    species_config_out.key_lengths_td.emerg_to_end = canopy_td.t_fphen_d - canopy_td.t_fphen_a
    species_config_out.key_lengths_td.emerg_to_veg = canopy_td.tt_veg
    species_config_out.key_lengths_td.veg_to_harvest = canopy_td.tt_rep

    species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_e = leaf_intervals.t_fphen_1_ets
    species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_g = leaf_intervals.t_fphen_3_ets
    species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_h = leaf_intervals.t_fphen_4_ets
    species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_i = leaf_intervals.t_fphen_5_ets
    species_config_out.key_lengths_flag_leaf_td.tl = leaf_intervals.t_l
    species_config_out.key_lengths_flag_leaf_td.tl_em = leaf_intervals.t_lem
    species_config_out.key_lengths_flag_leaf_td.tl_ma = leaf_intervals.t_lma
    species_config_out.key_lengths_flag_leaf_td.tl_ep = leaf_intervals.t_lep
    species_config_out.key_lengths_flag_leaf_td.tl_se = leaf_intervals.t_lse
    species_config_out.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg = species_config_out.key_lengths_td.emerg_to_astart - leaf_intervals.t_lem

    # 4. get intervals
    fphen_intervals = list(zip(*get_fphen_PLF(
        species_config_out.key_lengths_td.sowing_to_emerge,
        species_config_out.key_lengths_td.sowing_to_f_phen_b,
        species_config_out.key_lengths_td.sowing_to_f_phen_c,
        species_config_out.key_lengths_td.sowing_to_end,
        species_config_out.f_phen_min,
        x_offset=sowing_offset_td,
    )))

    leaf_fphen_intervals = list(zip(*get_leaf_fphen_PLF(
        species_config_out.leaf_f_phen_a,
        species_config_out.leaf_f_phen_b,
        species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_e,
        species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_g,
        species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_h,
        species_config_out.key_lengths_flag_leaf_td.leaf_f_phen_i,
        species_config_out.key_dates_td.Astart,
    )))

    dvi_interval = list(zip(*get_dvi_PLF(
        tt_emr=species_config_out.key_lengths_td.sowing_to_emerge,
        tt_veg=species_config_out.key_lengths_td.emerg_to_veg,
        tt_rep=species_config_out.key_lengths_td.veg_to_harvest,
        x_offset=sowing_offset_td,
    )))

    species_config_out.fphen_intervals = fphen_intervals
    species_config_out.leaf_fphen_intervals = leaf_fphen_intervals
    species_config_out.dvi_interval = dvi_interval

    return model_config_out, species_config_out


def phenology_from_flag_leaf_fraction(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
    external_data: NamedTuple,
    td_base_temperature: float,
) -> Tuple[ModelConfig, SpeciesConfig]:
    """Generate the phenology parameters from flag leaf length fractions.

    To use this method you must supply parameters that can be used to derive the
    length of the growing season in thermal time.

    Parameters
    ----------
    model_config : ModelConfig
        phenology model config
    species_config : SpeciesConfig
        species config
    external_data: NamedTuple
        Named tuple of external data
    td_base_temperature: float
        Base temp for thermal time. Only required for julian day season definition


    Returns
    -------
    Tuple[ModelConfig, SpeciesConfig]
        Updated comprehensive model config and species config

    """

    model_config_out = deepcopy(model_config)
    species_config_out = deepcopy(species_config)

    # 1. get Astart day
    Astart = None
    # TODO: Set astart
    if species_config.key_dates.Astart is not None:
        Astart = species_config.key_dates.Astart
    elif species_config.key_dates_td.Astart is not None:
        raise NotImplementedError("Cannot currently set Astart as thermal time")
    else:
        raise ValueError("Astart must be set in config for flag_leaf_fraction method")

    if species_config.key_dates.Aend:
        raise NotImplementedError("Must currently set sowing length in thermal time")
        # Aend = species_config.key_dates.Aend

        # growing_season_length = estimate_growing_season_from_flag_leaf_period(
        #     Astart,
        #     Aend,
        #     species_config.f_leaf_f_fphen,
        # )
        # SGS = Aend - growing_season_length
        # species_config_out.key_dates.sowing = int(SGS)

    elif species_config.key_lengths_td.sowing_to_end:
        growing_season_length_td = species_config.key_lengths_td.sowing_to_end
        t_astart = growing_season_length_td * species_config.f_Astart
        t_sgs = - t_astart
        t_egs = growing_season_length_td - t_astart
        species_config_out.key_dates_td.sowing = t_sgs
        species_config_out.key_dates_td.harvest = t_egs
        species_config_out.key_lengths_td.sowing_to_astart = t_astart
    else:
        raise ValueError("Aend must be set in config for flag_leaf_fraction method")
    assert t_sgs, "t_SGS missing"
    return phenology_from_season_fraction(
        model_config_out,
        species_config_out,
        external_data,
        td_base_temperature,
    )


def phenology_from_fphen_td_int(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
    output_time_type: TimeTypes,
) -> Tuple[ModelConfig, SpeciesConfig]:
    """Generate all phenology parameters from fphen td intervals.

    """
    growing_season_length = species_config.key_lengths_td.sowing_to_end
    # TODO: Implement Complex fphen
    fphen_intervals = list(zip(*get_fphen_PLF(
        species_config.key_lengths_td.sowing_to_emerge,
        species_config.key_lengths_td.sowing_to_f_phen_b,
        species_config.key_lengths_td.sowing_to_f_phen_c,
        species_config.key_lengths_td.sowing_to_end,
        species_config.f_phen_min,
    )))

    leaf_fphen_intervals = list(zip(*get_leaf_fphen_PLF(
        species_config.leaf_f_phen_a,
        species_config.leaf_f_phen_b,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_e,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_g,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_h,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_i,
        species_config.key_dates_td.Astart,
    )))
    tt_emr, tt_veg, tt_rep = get_dvi_from_season_length(
        0,
        growing_season_length,
        species_config.f_tt_emr,
        species_config.f_tt_veg,
        species_config.f_tt_rep,
    )
    dvi_interval = list(zip(*get_dvi_PLF(tt_emr, tt_veg, tt_rep)))

    canopy_td = get_canopy_td_intervals_f(
        0,
        growing_season_length,
        species_config.f_Astart,
        species_config.f_mid_anthesis,
        species_config.f_fphen_a,
        species_config.f_fphen_b,
        species_config.f_fphen_c,
        species_config.f_fphen_d,
        species_config.f_tt_emr,
        species_config.f_tt_veg,
        species_config.f_tt_rep,
    )
    leaf_intervals = get_leaf_td_intervals_f(
        t_sgs=0,
        t_egs=growing_season_length,
        f_Astart=species_config.f_Astart,
        f_mid_anthesis=species_config.f_mid_anthesis,
        f_fphen_1_ets=species_config.f_fphen_1_ets,
        f_fphen_3_ets=species_config.f_fphen_3_ets,
        f_fphen_4_ets=species_config.f_fphen_4_ets,
        f_fphen_5_ets=species_config.f_fphen_5_ets,
        f_t_lem=species_config.f_t_lem,
        f_t_lma=species_config.f_t_lma,
        f_t_lep=species_config.f_t_lep,
        f_t_lse=species_config.f_t_lse,
    )

    emergence = species_config.key_lengths_td.sowing_to_emerge
    harvest = growing_season_length
    Aend = growing_season_length
    mid_anthesis = species_config.key_dates_td.Astart + \
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_e
    emerg_to_astart = canopy_td.t_Astart - canopy_td.t_fphen_a
    emerg_to_flag_emerg = emerg_to_astart - leaf_intervals.t_lem

    key_dates_td = PhenologyKeyDates(
        sowing=0,  # All dates relative to sowing day
        emergence=emergence,
        harvest=harvest,
        Astart=species_config.key_dates_td.Astart,
        Aend=Aend,
        mid_anthesis=mid_anthesis,
    )

    key_lengths_td = replace(
        species_config.key_lengths_td,
        sowing_to_astart=canopy_td.t_Astart,
        emerg_to_astart=canopy_td.t_Astart - canopy_td.t_fphen_a,
        emerg_to_end=growing_season_length - emergence,
        emerg_to_veg=tt_veg,
        veg_to_harvest=tt_rep,
    )

    key_lengths_flag_leaf = replace(
        species_config.key_lengths_flag_leaf_td,
        tl=leaf_intervals.t_l,
        tl_em=leaf_intervals.t_lem,
        tl_ma=leaf_intervals.t_lma,
        tl_ep=leaf_intervals.t_lep,
        tl_se=leaf_intervals.t_lse,
        plant_emerg_to_leaf_emerg=emerg_to_flag_emerg,
    )

    model_config_out = replace(model_config,
                               time_type=output_time_type,
                               )

    species_config_out = replace(species_config,
                                 fphen_intervals=fphen_intervals,
                                 leaf_fphen_intervals=leaf_fphen_intervals,
                                 dvi_interval=dvi_interval,
                                 key_dates_td=key_dates_td,
                                 key_lengths_td=key_lengths_td,
                                 key_lengths_flag_leaf_td=key_lengths_flag_leaf,
                                 )

    return model_config_out, species_config_out


def calc_key_dates(
    species_config: SpeciesConfig,
    td_data: np.ndarray,
    dd_data: np.ndarray,
) -> SpeciesConfig:
    species_config_out = deepcopy(species_config)
    sowing_dd = int(species_config_out.key_dates.sowing)
    assert sowing_dd is not None
    sowing_td = species_config_out.key_dates_td.sowing or get_td_from_day(
        species_config_out.key_dates.sowing, dd_data, td_data)[1]
    td_dd_data = td_data - sowing_td
    _, harvest_dd = get_day_from_td(species_config_out.key_dates_td.harvest, dd_data, td_dd_data)
    _, emergence_dd = get_day_from_td(
        species_config_out.key_dates_td.emergence, dd_data, td_dd_data)
    _, f_phen_b_dd = get_day_from_td(
        species_config_out.key_lengths_td.sowing_to_f_phen_b, dd_data, td_dd_data)
    _, f_phen_c_dd = get_day_from_td(
        species_config_out.key_lengths_td.sowing_to_f_phen_c, dd_data, td_dd_data)
    # species_config_out.key_dates_td.sowing = sowing_td
    species_config_out.key_dates.harvest = harvest_dd
    species_config_out.key_dates.emergence = emergence_dd - sowing_dd
    species_config_out.key_lengths.sowing_to_emerge = emergence_dd - sowing_dd
    species_config_out.key_lengths.sowing_to_f_phen_b = f_phen_b_dd - sowing_dd
    species_config_out.key_lengths.sowing_to_f_phen_c = f_phen_c_dd - sowing_dd
    species_config_out.key_lengths.sowing_to_end = harvest_dd - sowing_dd
    return species_config_out


def validate_config(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
) -> bool:
    # TODO: Assert values add up correctly
    # TODO: Should these be warnings?
    try:
        if type(species_config.key_dates.sowing) != int and type(species_config.key_dates.sowing) != float:
            # TODO: Ideally we would validate sowing as int but can be float if set from netcdf
            raise ConfigError(
                f"Sowing date must be an integer but got '{species_config.key_dates.sowing}'. Check it has been set correctly")
        if model_config.time_type == TimeTypes.THERMAL_TIME:
            if species_config.f_Astart + species_config.f_fphen_1_ets + species_config.f_fphen_5_ets != 1.0:
                val = species_config.f_Astart + species_config.f_fphen_1_ets + species_config.f_fphen_5_ets
                raise ConfigError(
                    f"f_Astart + f_fphen_1_ets + f_fphen_5_ets != 1.0 Value is: {val}")
            if species_config.f_tt_emr + species_config.f_tt_rep + species_config.f_tt_veg != 1.0:
                val = species_config.f_tt_emr + species_config.f_tt_rep + species_config.f_tt_veg
                raise ConfigError(f"f_tt_emr + f_tt_rep + f_tt_veg != 1.0 Value is: {val}")
            if species_config.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg < 0:
                raise ConfigError(
                    f"flag emerges before plant emerges. Diff is ({species_config.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg}). Make sure f_t_lem({species_config.f_t_lem}) > f_Astart({species_config.f_Astart})")

        # CHECK ZERO DAY
        if model_config.zero_day == ZeroDayOptions.SOWING:
            pass
        elif model_config.zero_day == ZeroDayOptions.ASTART:
            pass
        elif model_config.zero_day == ZeroDayOptions.DATA_START:
            pass
        else:
            raise NotImplementedError(f"Zero day type {model_config.zero_day} not implemented")
        return True
    except (ConfigError, NotImplementedError) as e:
        raise e
    except Exception as e:
        raise ConfigError("Failed to validate phenology config")


def get_leaf_pop_phenology(
    species_config: SpeciesConfig,
    nP: int,
):
    """Get the phenology for none flag leaf populations.

    # TODO: This is currently arbitary and just fits the other leaf populations between
    plant emerge and flag leaf emerg

    """
    emerge_to_flag_emerg = species_config.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg
    flag_emer_to_end = species_config.key_lengths_td.emerg_to_end - emerge_to_flag_emerg
    # We split emerge_to_flag_emerge between all remaining populations
    # If we only have flag leaf(nP = 1) then we set a dummy population to grow between plant
    # emerge and flag emerg
    t_lem = emerge_to_flag_emerg / (nP - 1) if nP > 1 else emerge_to_flag_emerg
    # t_lma is a ratio of t_lem
    t_lma = t_lem * species_config.f_t_lma / species_config.f_t_lem

    # we cap the t_lma of sub leaf populations to the length between flag leaf emergence and end of canopy senescence.
    t_lma = min(t_lma, flag_emer_to_end)

    t_lep = t_lma * species_config.f_t_lep / species_config.f_t_lma
    t_lse = t_lma * species_config.f_t_lse / species_config.f_t_lma

    tl = t_lma + t_lem

    key_lengths_leaf_td = PhenologyLeafKeyLengths(
        tl=tl,
        tl_em=t_lem,
        tl_ma=t_lma,
        tl_ep=round(t_lep, 4),
        tl_se=round(t_lse, 4),
        leaf_f_phen_e=None,
        leaf_f_phen_g=None,
        leaf_f_phen_h=None,
        leaf_f_phen_i=None,
        # TODO: plant emerg to leaf emerg should be different
        # for each population
        plant_emerg_to_leaf_emerg=None,
    )

    return key_lengths_leaf_td


def process_phenology_config(
    model_config: ModelConfig,
    species_config: SpeciesConfig,
    external_data: dict,
    td_base_temperature: float = None,
    nP: int = 1,
    calculate_key_dates: bool = False,
    logger: Callable[[str, str], None] = print,
) -> Tuple[ModelConfig, SpeciesConfig]:
    """Convert many phenology options into a single DO3SE phenology config.

    Note: All values are relative to SGS
    """
    phenology_method = model_config.phenology_method
    if phenology_method == PhenologyMethods.DISABLED:
        return model_config, species_config

    species_config_set = deepcopy(species_config)

    species_preset = SpeciesPresetsParams.get(species_config.PRESET, None)
    if species_preset is not None:
        logger(f"Using preset: {species_config.PRESET}")
        species_config_set = merge_dataclasses(species_preset, species_config)

    time_type = model_config.time_type

    # sgs_key_day = model_config.sgs_key_day
    # if sgs_key_day == KeyDays.SGS:
    #     # TODO: if sgs day is not sowing day then we must reverse calculate sowing day
    #     if output_time_type != TimeTypes.JULIAN_DAY:
    #         raise NotImplementedError(
    #             "Cannot currently reverse calculate sowing day from emergence day for thermal time phenology")
    #     raise NotImplementedError(f"sgs day not implemented: {sgs_key_day}")
    # elif sgs_key_day == KeyDays.SOWING_DAY:
    #     sowing_offset = 0
    # else:
    #     raise NotImplementedError(f"sgs day not implemented: {sgs_key_day}")

    # GET PHENOLOGY VARS
    generated_configs = None
    if phenology_method == PhenologyMethods.FPHEN_THERMAL_TIME:
        logger("Phenology: Using thermal time")
        generated_configs = phenology_from_fphen_td_int(
            model_config=model_config,
            species_config=species_config_set,
            output_time_type=time_type,
        )
    elif phenology_method == PhenologyMethods.LEAF_FPHEN_DATA:
        logger("Phenology: Using leaf fphen data")
        # TODO manage using dd here also
        # TODO: Can we do this without td?

        time_data = calc_thermal_time_range(external_data.get('Ts_C'), t_b=td_base_temperature)
        leaf_f_phen_data = external_data.get('leaf_fphen')

        generated_configs = phenology_from_leaf_f_phen_data(
            model_config=model_config,
            species_config=species_config_set,
            leaf_f_phen_data=leaf_f_phen_data,
            td_data=time_data,
            dd_data=external_data.get('dd'),
        )
    elif phenology_method == PhenologyMethods.SEASON_FRACTION:
        logger("Phenology: Using season fraction")
        generated_configs = phenology_from_season_fraction(
            model_config=model_config,
            species_config=species_config_set,
            external_data=external_data,
            td_base_temperature=td_base_temperature,
        )
    elif phenology_method == PhenologyMethods.FLAG_LEAF_FRACTION:
        logger("Phenology: Using flag fraction")
        generated_configs = phenology_from_flag_leaf_fraction(
            model_config=model_config,
            species_config=species_config_set,
            external_data=external_data,
            td_base_temperature=td_base_temperature,
        )
    elif phenology_method == PhenologyMethods.LEGACY_DAY_PLF:
        logger("Phenology: Using legacy fraction")
        generated_configs = phenology_from_legacy_day_plf(
            model_config=model_config,
            species_config=species_config_set,
        )
    else:
        raise NotImplementedError(f"Phenology method not implemented: {phenology_method}")

    generated_model_config, generated_species_config = generated_configs
    if time_type == TimeTypes.THERMAL_TIME:
        if nP > 1:
            leaf_pop_phenology_t = get_leaf_pop_phenology(generated_species_config, nP)
            generated_species_config.key_lengths_leaf_td = leaf_pop_phenology_t
        else:
            # NOTE: Do we need to include leaf phenology if only using flag leaf?
            leaf_pop_phenology_t = get_leaf_pop_phenology(generated_species_config, nP)
            generated_species_config.key_lengths_leaf_td = leaf_pop_phenology_t
    elif time_type == TimeTypes.JULIAN_DAY:
        # TODO: Check if we need to set any leaf phenology if using JULIAN_DAY method
        pass
    else:
        raise NotImplementedError(f"TimeType: {time_type} is not implemented")

    # We override the generated parameters with parameters provided in initial config
    # This allows the user to override some generated parameters
    final_species_config = merge_dataclasses(generated_species_config, species_config)

    if calculate_key_dates:
        final_species_config = calc_key_dates(
            final_species_config, external_data.get('td'), external_data.get('dd'))

    if not validate_config(generated_model_config, final_species_config):
        raise ConfigError("Phenology Config Validation Failed")

    return generated_model_config, final_species_config

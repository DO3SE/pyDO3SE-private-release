# %%
import warnings
import functools
from copy import deepcopy
from data_helpers.cls_parsing import rgetattr, rsetattr, rdelattr

from pyDO3SE.version import config_version


class MigrationError(Exception):

    def __init__(self, exception, message=None):
        self.exception = exception
        self.message = message

    def __str__(self):
        if self.message:
            return 'Migration Error, {0} '.format(self.message)
        else:
            return 'Migration Error has been raised'


def set_version(input_config: dict, version) -> dict:
    output_config = deepcopy(input_config)
    output_config["VERSION"] = version
    return output_config


class Migrations():
    migrations = []

    def __init__(self):
        pass

    @classmethod
    def add_migration(cls, version: int):
        """Decorator that defines the version this migration is applied to."""
        def wrapper(func):
            @functools.wraps(func)
            def _inner(config, verbose=True, *args, **kwargs):
                if verbose:
                    print(f"Running migration: {func.__name__}")
                return func(config, *args, **kwargs)
            cls.migrations.append((version, _inner))
            return _inner
        return wrapper

    @classmethod
    def run_migrations(cls, input_config, input_version, verbose=True):
        if input_version == config_version:
            return input_config
        if verbose:
            print(f"Running migrations between version {input_version} and {config_version}")
        output_config = deepcopy(input_config)
        if not input_version:
            raise ValueError("Must define VERSION in the input config file")
        for (v, func) in cls.migrations:
            if int(v) > int(input_version):
                try:
                    output_config = func(output_config, verbose=verbose)
                except Exception as e:
                    raise MigrationError(e, f"Failed to run migration: {func.__name__}")
        output_config = set_version(output_config, config_version)
        return output_config


def _move_value(dict_in, path, target_path):
    """Simple helper to move dictionary fields."""
    initial_value = rgetattr(dict_in, path, None)
    new_dict = rsetattr(dict_in, target_path,
                        initial_value) if initial_value is not None else dict_in
    new_dict = rdelattr(new_dict, path)
    return new_dict


def _copy_value(dict_in, path, target_path):
    """Simple helper to copy dictionary fields."""
    initial_value = rgetattr(dict_in, path, None)
    new_dict = rsetattr(dict_in, target_path, initial_value)
    return new_dict


def _remove_value(dict_in, path):
    """Simple helper to remove dictionary fields."""
    new_dict = rdelattr(dict_in, path)
    return new_dict


def _replace_none(dict_in, path):
    if rgetattr(dict_in, path) is None:
        return rsetattr(dict_in, path, {})
    return dict_in


@Migrations.add_migration(1)
def update_version(config: dict) -> dict:
    """Update the config version."""
    config_out = deepcopy(config)
    config_out['VERSION'] = config_version
    return config_out


@Migrations.add_migration(1)
def update_met_nested(config: dict) -> dict:
    """Met inputs config is now nested to include additional met config."""
    met_inputs = config['Met']
    new_met = {"inputs": met_inputs}
    config_out = {**config, "Met": new_met}
    return config_out


@Migrations.add_migration(1)
def move_gsto_params(config: dict) -> dict:
    """Move the shared gsto params from multip and pn_gsto to gsto."""
    config_out = deepcopy(config)
    new_params = config_out['Land_Cover']['parameters']
    for params in new_params:
        params = _move_value(params, 'multip_gsto.f_SW_method', 'gsto.f_SW_method')
        params = _move_value(params, 'multip_gsto.SWP_min', 'gsto.SWP_min')
        params = _move_value(params, 'multip_gsto.SWP_max', 'gsto.SWP_max')
        params = _move_value(params, 'multip_gsto.fSWP_exp_curve', 'gsto.fSWP_exp_curve')
        params = _move_value(params, 'multip_gsto.fSWP_exp_a', 'gsto.fSWP_exp_a')
        params = _move_value(params, 'multip_gsto.fSWP_exp_b', 'gsto.fSWP_exp_b')
        params = _move_value(params, 'multip_gsto.f_VPD_method', 'gsto.f_VPD_method')
        params = _move_value(params, 'multip_gsto.VPD_min', 'gsto.VPD_min')
        params = _move_value(params, 'multip_gsto.VPD_max', 'gsto.VPD_max')
        params = _move_value(params, 'multip_gsto.VPD_crit', 'gsto.VPD_crit')
    return config_out


@Migrations.add_migration(1)
def rename_leaf_phenology_method(config: dict) -> dict:
    """Rename leaf_phenology_method to leaf_f_phen_Anet_influence."""
    config_out = deepcopy(config)
    new_params = config_out['Land_Cover']['parameters']
    for params in new_params:
        params = _move_value(params, 'pn_gsto.leaf_phenology_method',
                             'pn_gsto.leaf_f_phen_Anet_influence')
    return config_out


@Migrations.add_migration(1)
def rename_leaf_phenology_method(config: dict) -> dict:
    """Rename leaf_phenology_method to leaf_f_phen_Anet_influence."""
    config_out = deepcopy(config)
    new_params = config_out['Land_Cover']['parameters']
    for params in new_params:
        params = _move_value(params, 'pn_gsto.leaf_phenology_method',
                             'pn_gsto.leaf_f_phen_Anet_influence')
    return config_out


@Migrations.add_migration(5)
def rename_fvpd_methods(config: dict) -> dict:
    """Rename leaf_phenology_method to leaf_f_phen_Anet_influence."""
    config_out = deepcopy(config)
    # nL = rgetattr(config, 'Land_Cover.nL', 1) or 1
    nLC = rgetattr(config, 'Land_Cover.nLC', 1) or 1
    assert nLC is not None
    for iLC in range(nLC):
        if rgetattr(config, f'Land_Cover.parameters.{iLC}.gsto.f_VPD_method') == "photosynthesis":
            config_out = rsetattr(
                config_out, f'Land_Cover.parameters.{iLC}.gsto.f_VPD_method', "leuning")
    return config_out


@Migrations.add_migration(6)
def set_emerg_day(config: dict) -> dict:
    """Set emerg day if life span method is constant."""
    config_out = deepcopy(config)

    gsto_method = rgetattr(config, 'Land_Cover.parameters.0.gsto.method', None)
    life_span_method = rgetattr(config, "Land_Cover.parameters.0.pn_gsto.life_span_method")
    if gsto_method == "photosynthesis" and life_span_method == "constant":
        warnings.warn("Check t_emerg in Land_Cover.parameters.0.season.t_emerg")
        nLC = rgetattr(config, 'Land_Cover.nLC', 1) or 1
        assert nLC is not None
        for iLC in range(nLC):
            config_out = rsetattr(
                config_out, f'Land_Cover.parameters.{iLC}.season.t_emerg', 0)

    return config_out


@Migrations.add_migration(7)
def set_dvi_method(config: dict) -> dict:
    """Set the dvi method to JULES if life span method was jules."""
    config_out = deepcopy(config)

    life_span_method = rgetattr(config, "Land_Cover.parameters.0.pn_gsto.life_span_method")
    if life_span_method == "JULES":
        config_out = rsetattr(config_out, 'Land_Cover.dvi_method', "JULES")
    return config_out


@Migrations.add_migration(8)
def set_r_dc_coeff(config: dict) -> dict:
    """Set R_dc coefficient for dark repiration."""
    config_out = deepcopy(config)
    default_R_d_coeff = 0.015

    gsto_method = rgetattr(config, 'Land_Cover.parameters.0.gsto.method', None)
    if gsto_method == "photosynthesis":
        warnings.warn("Dark respiration coefficient (R_dc) has been set to wheat default.")
        nLC = rgetattr(config, 'Land_Cover.nLC', 1) or 1
        assert nLC is not None
        for iLC in range(nLC):
            config_out = rsetattr(
                config_out, f'Land_Cover.parameters.{iLC}.pn_gsto.R_d_coeff', default_R_d_coeff)

    return config_out


@Migrations.add_migration(9)
def set_emergence_method(config: dict) -> dict:
    """Set the method of defining if plant has emerged."""
    config_out = deepcopy(config)
    config_out = rsetattr(config_out, 'phenology', {})
    config_out = rsetattr(config_out, 'phenology.model', {})
    config_out = rsetattr(config_out, 'phenology.species', [{"key_dates": {}, "key_dates_td": {}}])

    t_emerge = rgetattr(config_out, 'Land_Cover.parameters.0.season.t_emerg')
    emerge_day = rgetattr(config_out, 'Land_Cover.parameters.0.season.emerge_day')
    if t_emerge or emerge_day:
        warnings.warn("Only 1st land cover updated")
        config_out = rsetattr(config_out, 'phenology.model.plant_emerge_method', 'constant')
        config_out = rsetattr(config_out, 'phenology.species.0.key_dates.emergence', emerge_day)
        config_out = rsetattr(config_out, 'phenology.species.0.key_dates_td.emergence', t_emerge)
    life_span_method = rgetattr(config_out, 'Land_Cover.parameters.0.pn_gsto.life_span_method')
    if life_span_method == "constant":
        if not (t_emerge or emerge_day):
            warnings.warn("Must set t_emerge or emerge_day")
        config_out = rsetattr(config_out, 'phenology.species.0.key_dates.emergence', 0)
        config_out = rsetattr(config_out, 'phenology.species.0.key_dates_td.emergence', 0)
        config_out = rsetattr(config_out, 'phenology.model.plant_emerge_method', 'constant')
    elif life_span_method == "JULES":
        config_out = rsetattr(config_out, 'phenology.model.plant_emerge_method', 'dvi')
        # TODO: Check dvi method
    elif life_span_method == "leaf_f_phen":
        config_out = rsetattr(config_out, 'phenology.model.plant_emerge_method', 'leaf_f_phen')
    return config_out
    # TODO: Set other plant emerge methods


@Migrations.add_migration(10)
def move_phenology_options(config: dict) -> dict:
    """Move all phenology options to a single location."""

    config_out = deepcopy(config)
    config_out = _move_value(config_out, 'phenology.model', 'Land_Cover.phenology_options')
    config_out = _move_value(config_out, 'phenology.species.0', 'Land_Cover.parameters.0.phenology')

    t_sgs = rgetattr(config_out, 'Land_Cover.parameters.0.season.t_sgs')
    t_egs = rgetattr(config_out, 'Land_Cover.parameters.0.season.t_egs')
    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.season.t_sgs')
    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.season.t_egs')
    sowing_to_end = (t_sgs is not None and t_egs is not None) and t_egs - t_sgs

    config_out = _replace_none(config_out, 'Land_Cover.parameters.0.phenology.key_lengths_leaf_td')
    config_out = _replace_none(config_out, 'Land_Cover.parameters.0.phenology.key_lengths_td')
    config_out = _replace_none(config_out, 'Land_Cover.phenology_options')

    if sowing_to_end:
        config_out = rsetattr(
            config_out, 'Land_Cover.parameters.0.phenology.key_lengths_td.sowing_to_end', sowing_to_end)
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.sowing_day',
                             'Land_Cover.parameters.0.phenology.key_dates.sowing')

    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.SGS',
                             'Land_Cover.parameters.0.phenology.key_dates.emergence')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.EGS',
                             'Land_Cover.parameters.0.phenology.key_dates.EGS')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.Astart',
                             'Land_Cover.parameters.0.phenology.key_dates.Astart')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.Aend',
                             'Land_Cover.parameters.0.phenology.key_dates.Aend')

    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.LAI_a',
                             'Land_Cover.parameters.0.phenology.LAI_a')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.LAI_b',
                             'Land_Cover.parameters.0.phenology.LAI_b')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.LAI_c',
                             'Land_Cover.parameters.0.phenology.LAI_c')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.LAI_d',
                             'Land_Cover.parameters.0.phenology.LAI_d')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.LAI_1',
                             'Land_Cover.parameters.0.phenology.LAI_1')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.LAI_2',
                             'Land_Cover.parameters.0.phenology.LAI_2')

    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.SAI_method',
                             'Land_Cover.parameters.0.phenology.SAI_method')

    life_span_method = rgetattr(config_out, 'Land_Cover.parameters.0.pn_gsto.life_span_method')
    if life_span_method == "constant":
        warnings.warn(
            "Old life_span_method \"constant\" is invalid. Check \"Land_Cover.parameters.0.phenology.key_lengths_leaf_td\" values")
        config_out = _copy_value(config_out, 'Land_Cover.parameters.0.pn_gsto.t_l_estimate',
                                 'Land_Cover.parameters.0.phenology.key_lengths_td.sowing_to_end')
        config_out = _move_value(config_out, 'Land_Cover.parameters.0.pn_gsto.t_l_estimate',
                                 'Land_Cover.parameters.0.phenology.key_lengths_leaf_td.tl')
        config_out = _move_value(config_out, 'Land_Cover.parameters.0.pn_gsto.t_lem',
                                 'Land_Cover.parameters.0.phenology.key_lengths_leaf_td.tl_em')
        config_out = _move_value(config_out, 'Land_Cover.parameters.0.pn_gsto.t_lma',
                                 'Land_Cover.parameters.0.phenology.key_lengths_leaf_td.tl_ma')
        config_out = _move_value(config_out, 'Land_Cover.parameters.0.pn_gsto.t_lse',
                                 'Land_Cover.parameters.0.phenology.key_lengths_leaf_td.tl_se')
        config_out = _move_value(config_out, 'Land_Cover.parameters.0.pn_gsto.t_lep',
                                 'Land_Cover.parameters.0.phenology.key_lengths_leaf_td.tl_ep')
        config_out = rsetattr(
            config_out, 'Land_Cover.phenology_options.phenology_method', 'season_fraction')

        sowing_to_end = rgetattr(
            config_out, 'Land_Cover.parameters.0.phenology.key_lengths_td.sowing_to_end')
        assert sowing_to_end is not None and sowing_to_end > 0

    elif life_span_method == "JULES":
        config_out = rsetattr(
            config_out, 'Land_Cover.phenology_options.phenology_method', 'development index')
    elif life_span_method == "leaf_f_phen":
        config_out = rsetattr(
            config_out, 'Land_Cover.phenology_options.phenology_method', 'leaf_fphen_data')
        config_out = rsetattr(
            config_out, 'Land_Cover.phenology_options.zero_day', 'Astart')
    elif life_span_method == "t_leaf_f_phen":
        config_out = rsetattr(
            config_out, 'Land_Cover.phenology_options.phenology_method', 'fphen_thermal_time')
    elif life_span_method == "growing_season":
        config_out = rsetattr(
            config_out, 'Land_Cover.phenology_options.phenology_method', 'season_fraction')
        assert sowing_to_end is not None and sowing_to_end > 0

    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.pn_gsto.life_span_method')
    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.pn_gsto.tt_emr')
    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.pn_gsto.tt_veg')
    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.pn_gsto.tt_rep')

    config_out = _move_value(config_out, 'Land_Cover.parameters.0.season.EGS',
                             'Land_Cover.parameters.0.phenology.key_dates.EGS')

    if rgetattr(config_out, 'Land_Cover.parameters.0.phenology.PRESET') is None:
        config_out = rsetattr(
            config_out, 'Land_Cover.parameters.0.phenology.PRESET', "WHEAT_SPRING")
        warnings.warn("Using Spring Wheat preset")

    return config_out


@Migrations.add_migration(10)
def move_fphen_day_params(config: dict) -> dict:
    config_out = deepcopy(config)
    photosynthesis_method = rgetattr(config, 'Land_Cover.parameters.0.gsto.method')
    config_out = rsetattr(config_out, 'Land_Cover.parameters.0.phenology.day_fphen_plf', {})

    if photosynthesis_method == "multiplicative":
        f_phen_method = rgetattr(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_method')
        if f_phen_method in ["simple day PLF", "complex day PLF"]:

            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_limA',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_limA')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_limB',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_limB')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_a',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_a')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_b',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_b')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_c',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_c')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_d',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_d')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_e',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_e')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_1',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_1')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_2',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_2')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_3',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_3')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_4',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.f_phen_4')

            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.leaf_f_phen_a',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_f_phen_a')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.leaf_f_phen_b',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_f_phen_b')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.leaf_f_phen_c',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_f_phen_c')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.leaf_f_phen_1',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_f_phen_1')
            config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.leaf_f_phen_2',
                                     'Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_f_phen_2')
            config_out = rsetattr(
                config_out, 'Land_Cover.phenology_options.phenology_method', "fphen_julian_days")

        else:
            raise ValueError(f"Multiplicative Fphen method {f_phen_method} migration not setup")
            # TODO: Update phenology method etc

    config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.f_phen_method',
                             'Land_Cover.parameters.0.phenology.f_phen_method')
    config_out = _move_value(config_out, 'Land_Cover.parameters.0.multip_gsto.leaf_f_phen_method',
                             'Land_Cover.parameters.0.phenology.leaf_f_phen_method')

    return config_out


@Migrations.add_migration(10)
def remove_redundant_options(config: dict) -> dict:
    """Remove all options that are deprecated."""
    config_out = deepcopy(config)
    if rgetattr(config_out, 'Land_Cover.parameters.0.pn_gsto.O3_method') is not None:
        warnings.warn("pn_gsto O3_method is deprecated")
    if rgetattr(config_out, 'Land_Cover.parameters.0.season.growing_season_method') is not None:
        warnings.warn("growing_season_method is deprecated")
    if rgetattr(config_out, 'Land_Cover.parameters.0.season.accumulation_period_method') is not None:
        warnings.warn("accumulation_period_method is deprecated")

    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.pn_gsto.O3_method')
    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.season.growing_season_method')
    config_out = _remove_value(
        config_out, 'Land_Cover.parameters.0.season.accumulation_period_method')

    config_out = _remove_value(
        config_out, 'Land_Cover.parameters.0.season.t_emerg')

    config_out = _remove_value(config_out, 'phenology')

    old_season_config = rgetattr(config_out, 'Land_Cover.parameters.0.season')
    if len(old_season_config.keys()) > 0:
        print(old_season_config.keys())
        raise ValueError("Old season config not fully migrated")

    config_out = _remove_value(config_out, 'Land_Cover.parameters.0.season')

    return config_out


@Migrations.add_migration(11)
def change_sai_method_default(config: dict) -> dict:
    """Change default SAI version."""
    config_out = deepcopy(config)
    try:
        if rgetattr(config_out, 'Land_Cover.parameters.0.phenology.SAI_method') == "LAI":
            config_out = rsetattr(
                config_out, 'Land_Cover.parameters.0.phenology.SAI_method', "LAI_max")
    except AttributeError as e:
        if "'NoneType' object has no attribute" not in str(e):
            raise e
    return config_out


@Migrations.add_migration(12)
def set_senescence_method_default(config: dict) -> dict:
    """Default senescence_method has been set to disabled."""
    config_out = deepcopy(config)
    try:
        if not rgetattr(config_out, 'Land_Cover.parameters.0.pn_gsto.senescence_method', {}):
            gsto_method = rgetattr(config_out, 'Land_Cover.parameters.0.gsto.method', None)
            if gsto_method == "photosynthesis":
                config_out = rsetattr(
                    config_out, 'Land_Cover.parameters.0.pn_gsto.senescence_method', "ewert", create_missing_dicts=True)
            elif gsto_method == "multiplicative":
                config_out = rsetattr(
                    config_out, 'Land_Cover.parameters.0.pn_gsto.senescence_method', "anet", create_missing_dicts=True)
    except AttributeError as e:
        if "'NoneType' object has no attribute" not in str(e):
            raise e
    return config_out

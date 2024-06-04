"""Config parameters are validated against the lists of required fields below."""


from data_helpers.list_helpers import flatten_list
from data_helpers.cls_parsing import rgetattr
from pyDO3SE.version import config_version
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Config.ConfigEnums import CanopyHeightMethods


class ConfigValidationError(Exception):
    pass


def required_config_fields_multip_gsto(config): return list(filter(lambda v: v, flatten_list([
    [[
        f"Land_Cover.parameters.{iLC}.phenology.f_phen_method",
        f"Land_Cover.parameters.{iLC}.phenology.leaf_f_phen_method",
        f"Land_Cover.parameters.{iLC}.multip_gsto.f_light_method",
        f"Land_Cover.parameters.{iLC}.multip_gsto.f_temp_method",
        [
            f"Land_Cover.parameters.{iLC}.multip_gsto.T_min",
            f"Land_Cover.parameters.{iLC}.multip_gsto.T_opt",
            f"Land_Cover.parameters.{iLC}.multip_gsto.T_max",
        ] if config.Land_Cover.parameters[iLC].multip_gsto.f_temp_method != "disabled" else None,

        f"Land_Cover.parameters.{iLC}.multip_gsto.gmax",
        f"Land_Cover.parameters.{iLC}.multip_gsto.gmorph",
        f"Land_Cover.parameters.{iLC}.multip_gsto.gmorph",
    ] for iLC in range(config.Land_Cover.nLC or 0)],
])))


def required_config_fields_anet_gsto(config): return list(filter(lambda v: v, flatten_list([
    [[
        f"Land_Cover.parameters.{iLC}.pn_gsto.g_sto_0",
        f"Land_Cover.parameters.{iLC}.pn_gsto.m",
        config.Land_Cover.parameters[iLC].pn_gsto.D_0_method == "constant" and [
            f"Land_Cover.parameters.{iLC}.pn_gsto.D_0",
        ],

        config.Land_Cover.parameters[iLC].pn_gsto.D_0_method == "f_VPD" and [
            f"Land_Cover.parameters.{iLC}.gsto.f_VPD_method",  # Should not be disabled!
        ],

        f"Land_Cover.parameters.{iLC}.pn_gsto.V_cmax_25",
        f"Land_Cover.parameters.{iLC}.pn_gsto.J_max_25",
    ] for iLC in range(config.Land_Cover.nLC or 0)],
])))


def required_config_fields(config): return list(filter(lambda v: v, flatten_list([
    "Location.lat",
    "Location.lon",
    "Location.elev",
    not config.Location.OTC and "Location.z_u",
    not config.Location.OTC and "Location.z_O3",
    "Land_Cover.nL",
    "Land_Cover.nLC",
    # [f"Land_Cover.parameters.{iLC}.season.SGS" or f"Land_Cover.parameters.{iLC}.season.t_sgs" for iLC in range(
    #     config.Land_Cover.nLC or 0)],
    [f"Land_Cover.parameters.{iLC}.height" for iLC in range(
        config.Land_Cover.nLC or 0) if config.Land_Cover.height_method == CanopyHeightMethods.CONSTANT],
    [f"Land_Cover.parameters.{iLC}.Lm" for iLC in range(config.Land_Cover.nLC or 0)],
    [f"Land_Cover.parameters.{iLC}.Y" for iLC in range(config.Land_Cover.nLC or 0)],
    # === GSTO ===
    config.Land_Cover.parameters[0].gsto.method == "multiplicative"
    and required_config_fields_multip_gsto(config),
    config.Land_Cover.parameters[0].gsto.method == "photosynthesis"
    and required_config_fields_anet_gsto(config),

    # f_VPD
    [config.Land_Cover.parameters[iLC].gsto.f_VPD_method == "linear" and [
        f"Land_Cover.parameters.{iLC}.gsto.VPD_min",
        f"Land_Cover.parameters.{iLC}.gsto.VPD_max",
    ] for iLC in range(config.Land_Cover.nLC or 0)],
    # fSWP
    [config.Land_Cover.parameters[iLC].gsto.fSWP_exp_curve in ["temperate", "mediterranean"] and [
        f"Land_Cover.parameters.{iLC}.gsto.fSWP_exp_a",
        f"Land_Cover.parameters.{iLC}.gsto.fSWP_exp_b",
    ] for iLC in range(config.Land_Cover.nLC or 0)],
    # f_SW
    [config.Land_Cover.parameters[iLC].gsto.f_SW_method in ["fSWP exp", "fLWP exp"] and [
        f"Land_Cover.parameters.{iLC}.gsto.fSWP_exp_a",
        f"Land_Cover.parameters.{iLC}.gsto.fSWP_exp_b",
    ] for iLC in range(config.Land_Cover.nLC or 0)],
    [config.Land_Cover.parameters[iLC].gsto.f_SW_method in ["linear"] and [
        f"Land_Cover.parameters.{iLC}.gsto.SWP_min",
        f"Land_Cover.parameters.{iLC}.gsto.SWP_max",
    ] for iLC in range(config.Land_Cover.nLC or 0)],
    # Phenology
])))


def validate_config_fields(config: Config_Shape):
    for path in required_config_fields(config):
        val = rgetattr(config, path)
        if val is None:
            raise ConfigValidationError(
                f"Config value missing: {path}")


def validate_config_version(config: Config_Shape):
    if config.VERSION != config_version:
        raise ConfigValidationError(
            f"Config version is incorrect. Run config migration to update.")

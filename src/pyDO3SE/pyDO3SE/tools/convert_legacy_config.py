"""Load old DO3SE UI 3.0 configs and convert to pyDO3SE 4.0"""

import json
from data_helpers.cls_parsing import rsetattr
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.Output.process_outputs import dump_config_to_file_json

conversion_map = {
    "input_fields": None,
    "input_trim": None,
    "lat": "Location.lat",
    "lon": "Location.lon",
    "elev": "Location.elev",
    "albedo": "Location.albedo",
    # "co2_constant": "Met.inputs.CO2_constant", // TODO: This is an object as input
    "soil_tex": "soil_moisture.soil_texture",
    "rsoil": "Location.Rsoil",
    "o3zr": "Location.z_O3",
    "o3_h": "Location.h_O3",
    "uzr": "Location.z_u",
    "u_h": "Location.z_u",
    # "ustar_method": None, # TODO: Handle ustar_method
    # "d_meas": "d_meas", # Not currently used in pyDO3SE
    "h": "Land_Cover.parameters.0.height",
    "root": "soil_moisture.root",
    "lm": "Land_Cover.parameters.0.Lm",
    "gmax": "Land_Cover.parameters.0.multip_gsto.gmax",
    "gmorph": "Land_Cover.parameters.0.multip_gsto.gmorph",
    "fmin": "Land_Cover.parameters.0.gsto.fmin",
    # "rext": "rext", # TODO: Not currently used in pyDO3SE
    "y": "Land_Cover.parameters.0.Y",
    "f_lightfac": "Land_Cover.parameters.0.multip_gsto.f_lightfac",
    "t_min": "Land_Cover.parameters.0.multip_gsto.T_min",
    "t_opt": "Land_Cover.parameters.0.multip_gsto.T_opt",
    "t_max": "Land_Cover.parameters.0.multip_gsto.T_max",
    "vpd_max": "Land_Cover.parameters.0.gsto.VPD_max",
    "vpd_min": "Land_Cover.parameters.0.gsto.VPD_min",
    "vpd_crit": "Land_Cover.parameters.0.gsto.VPD_crit",
    "swp_min": "Land_Cover.parameters.0.gsto.SWP_min",
    "swp_max": "Land_Cover.parameters.0.gsto.SWP_max",
    # "ra_method": "ra_method", # TODO: Match up Ramethod
    "gsto": "Land_Cover.parameters.0.gsto.method",
    # "tleaf": "Land_Cover.parameters.0.gsto.Tleaf_method", # TODO: "estimate" is not in new methods
    "fo3": "Land_Cover.parameters.0.multip_gsto.f_O3_method",  # TODO: Match up
    # "fxwp": "fxwp",# NOTE: managed seperately
    # "lwp": "lwp",# NOTE: managed seperately
    # "asw": "asw",# NOTE: managed seperately
    "asw_fc_override": "soil_moisture.ASW_FC",  # TODO: Match up
    "asw_min": "Land_Cover.parameters.0.gsto.ASW_min",
    "asw_max": "Land_Cover.parameters.0.gsto.ASW_max",
    "sgs_egs_calc": "Land_Cover.phenology_options.sowing_day_method",  # Also needs to capitalize!
    "sgs": "Land_Cover.parameters.0.phenology.key_dates.sowing",
    "egs": "Land_Cover.parameters.0.phenology.key_dates.harvest",
    "mid_anthesis": "Land_Cover.parameters.0.phenology.key_dates.mid_anthesis",
    "lai_a": "Land_Cover.parameters.0.phenology.LAI_a",
    "lai_b": "Land_Cover.parameters.0.phenology.LAI_b",
    "lai_c": "Land_Cover.parameters.0.phenology.LAI_c",
    "lai_d": "Land_Cover.parameters.0.phenology.LAI_d",
    "lai_1": "Land_Cover.parameters.0.phenology.LAI_1",
    "lai_2": "Land_Cover.parameters.0.phenology.LAI_2",
    "sai": "Land_Cover.parameters.0.phenology.SAI_method",
    "fphen_a": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_a",
    "fphen_b": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_b",
    "fphen_c": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_c",
    "fphen_d": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_d",
    "fphen_e": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_e",
    "fphen_1": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_1",
    "fphen_lima": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_lima",
    "fphen_2": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_2",
    "fphen_3": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_3",
    "fphen_limb": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_limb",
    "fphen_4": "Land_Cover.parameters.0.phenology.day_fphen_plf.fphen_4",
    # "leaf_fphen": "Land_Cover.parameters.0.phenology.leaf_f_phen_method", # values do not match
    "astart": "Land_Cover.parameters.0.phenology.key_dates.Astart",
    "aend": "Land_Cover.parameters.0.phenology.key_dates.Aend",
    "leaf_fphen_a": "Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_fphen_a",
    "leaf_fphen_b": "Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_fphen_b",
    "leaf_fphen_c": "Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_fphen_c",
    "leaf_fphen_1": "Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_fphen_1",
    "leaf_fphen_2": "Land_Cover.parameters.0.phenology.day_fphen_plf.leaf_fphen_2",
}

defaults_map = {
    "Land_Cover.parameters.0.gsto.f_VPD_method": "linear",
    "Land_Cover.parameters.0.multip_gsto.f_temp_method": "default",
    "Land_Cover.parameters.0.multip_gsto.f_light_method": "enabled",
    "Land_Cover.phenology_options.phenology_method": "day plf",  # TODO: Check this
    "Land_Cover.phenology_options.sgs_key_day": "SGS",  # TODO: Check this
    "Land_Cover.phenology_options.time_type": "julian_day",  # TODO: Check this
    "Land_Cover.phenology_options.plant_emerge_method": "SGS"
}


def convert_legacy_config(config_path: str):
    with open(config_path) as config_file_data:
        legacy_config = json.load(config_file_data)
        new_config = Config_Shape()
        # Convert config here
        # Save to output path
        for k in conversion_map.keys():
            if conversion_map[k] is not None:
                rsetattr(new_config, conversion_map[k], legacy_config.get(k, None))
            else:
                print(f"Key {k} not found in conversion map")
        # Apply defaults
        for k in defaults_map.keys():
            rsetattr(new_config, k, defaults_map[k])

        # Fix soil moisture methods
        if legacy_config.get('fxwp', None) == "disabled":
            rsetattr(new_config, "Land_Cover.parameters.0.gsto.f_SW_method", "disabled")
        if legacy_config.get('fxwp', None) == "fswp":
            if legacy_config.get('fswp', None) == "input":
                # rsetattr(new_config, "Land_Cover.parameters.0.gsto.f_SW_method","disabled")
                # TODO: Check what to do with this
                raise NotImplementedError("fswp input not implemented")
            elif legacy_config.get('fswp', None) == "exp":
                rsetattr(new_config, "Land_Cover.parameters.0.gsto.f_SW_method", "fSWP exp")
            elif legacy_config.get('fswp', None) == "linear":
                rsetattr(new_config, "Land_Cover.parameters.0.gsto.f_SW_method", "fSWP linear")
        if legacy_config.get('fxwp', None) == "flwp":
            # NOTE: flwp is not implemented in pyDO3SE
            raise NotImplementedError("flwp not implemented")
        if legacy_config.get('fxwp', None) == "fpaw":
            rsetattr(new_config, "Land_Cover.parameters.0.gsto.f_SW_method", "fPAW")

        # Fix nested legacy params
        if legacy_config.get('co2_constant', None):
            if legacy_config.get('co2_constant', None)['disabled'] == False:
                rsetattr(new_config, "Met.inputs.CO2_constant",
                         legacy_config.get('co2_constant', None)['value'])
            else:
                raise NotImplementedError("CO2 constant disabled not implemented")
        if legacy_config.get('o3_h', None):
            #
            if legacy_config.get('o3_h', None)['disabled'] == False:
                rsetattr(new_config, "Location.h_O3", legacy_config.get('o3_h', None)['value'])
            else:
                rsetattr(new_config, "Location.h_O3", legacy_config.get('h', None))
        if legacy_config.get('u_h', None):
            #
            if legacy_config.get('u_h', None)['disabled'] == False:
                rsetattr(new_config, "Location.z_u", legacy_config.get('u_h', None)['value'])
            else:
                rsetattr(new_config, "Location.h_O3", legacy_config.get('h', None))

        # Fix ASW
        if legacy_config.get('asw', None) == "input":
            rsetattr(new_config, "Met.inputs.ASW_method", "input")
            rsetattr(new_config, "soil_moisture.source", "external input ASW")
        # Fix ustar
        if legacy_config.get('ustar_method', None) == "input_ref":
            rsetattr(new_config, "Met.inputs.ustar_ref_method", "input")
        if legacy_config.get('ustar_method', None) == "input":
            rsetattr(new_config, "Met.inputs.ustar_method", "input")

        # Fix leaf f phen method
        if legacy_config.get('leaf_fphen', None):
            if legacy_config.get('leaf_fphen', None) == "copy":
                rsetattr(new_config, "Land_Cover.parameters.0.phenology.leaf_f_phen_method", "f_phen")
            if legacy_config.get('leaf_fphen', None) == "fixedday":
                rsetattr(new_config, "Land_Cover.parameters.0.phenology.leaf_f_phen_method", "day PLF")
            if legacy_config.get('leaf_fphen', None) == "input":
                rsetattr(new_config, "Land_Cover.parameters.0.phenology.leaf_f_phen_method", "input")

        # Fix fO3 method
        if legacy_config.get('fo3', "none") == "none":
            rsetattr(new_config, "Land_Cover.parameters.0.multip_gsto.f_O3_method", "disabled")

        # Fix sgs_method
        # ('inputs',       options.sgs_egs_use_inputs,     'Use inputs below'),
        # ('latitude',     options.sgs_egs_latitude,       'Latitude function'),
        # ('thermal_time', options.sgs_egs_tt,
        # 'ETS ((Bread) Wheat (atlantic, boreal, continental))'),
        # ('thermal_time_mb', options.sgs_egs_tt_mb,
        # 'ETS ((Bread) Wheat (Mediterranean))'),
        # ('thermal_time_md', options.sgs_egs_tt_md,
        # 'ETS ((Durum) Wheat (Mediterranean))'),
        # ('thermal_time_pot', options.sgs_egs_tt_pot,      'ETS (Potato)'),
        # ('thermal_time_tom', options.sgs_egs_tt_tom,      'ETS (Tomato)'),
        # TO:
        #  ['INPUT', 'LATITUDE', 'LATITUDE_FOREST', 'LATITUDE_SPRING_EUROPE', 'LATITUDE_WINTER_CHINA', 'SKIP']

        if legacy_config.get('sgs_egs_calc', "none") == "inputs":
            rsetattr(new_config, "Land_Cover.phenology_options.sowing_day_method", "INPUT")
        elif legacy_config.get('sgs_egs_calc', "none") == "latitude":
            rsetattr(new_config, "Land_Cover.phenology_options.sowing_day_method", "LATITUDE_FOREST")
            rsetattr(new_config, "Land_Cover.phenology_options.latitude",
                     legacy_config.get('lat', None))
        else:
            raise NotImplementedError(
                f"SGS EGS method not implemented: {legacy_config.get('sgs_egs_calc', 'none')}")

        print("Make sure to update `Met.inputs` and `output.fields` to match the new config structure")

        return new_config


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python convert_legacy_config.py <legacy_config_path> <output_path>")
        sys.exit(1)
    new_config = convert_legacy_config(sys.argv[1])
    output_path = sys.argv[2]
    dump_config_to_file_json(new_config, output_path)

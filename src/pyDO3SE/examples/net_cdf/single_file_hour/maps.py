
# Superceded by individual files but kept for comments
variable_map = {
    'time': 'XTIME',
    'Ts_C': 'td_2m',
    'P': 'pres',
    'PAR': 'SWDOWN',  # Check is this PAR or PPFD
    'precip': 'RAINNC',  # Should be sum of this and RAINC
    'RH': 'rh',
    'u': 'wspeed',
    'O3': 'o3',
    'Hd': 'HFX_FORCE',
    # 'SWC': 'SMOISREL',
    'snow_depth': 'SNOWH',

}

preprocess_map = {}

e_state_overrides_field_map = {
    'Location.lat': 'lat',
    'Location.lon': 'lon',
    # TODO: Add additional values here
}

variable_map = dict([
    ('time', 'time'),
    ('O3', 'O3_45m'),
    ('u', 'u_45m'),
    ('Ts_C', 't2m'), # Needs transform
    ('Hd', 'SH_Wm2'),
    ('precip', 'Precip'),
    ('cloudfrac', 'CloudFrac'),
    ('P', 'Psurf'),
    ('RH', 'rh2m'),
    # ('ustar', 'ustar'),
])


preprocess_map = dict([
    ('Ts_C', lambda row: row - 273.15),
])

# %%

from math import floor
from matplotlib import pyplot as plt
t_l = 30
t_lem = t_l * 0.3
t_lma = t_l * 0.6
# %%

fO3_h_prev = 1
fO3_d_prev = 1
rO3_h = 1
fO3_d_values = []
for i in range(24 * t_l + 2):
    hr = i % 24
    # print(hr)
    dd = floor(i / 24)
    fO3_h = 0.9 if hr in [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] else 1
    f_LA = 1 if dd < t_lem else 1 - (dd - t_lem) / t_lma if dd < t_l - 2 else 0
    rO3_h = fO3_d_prev + (1 - fO3_d_prev) * f_LA
    fO3_d_base = fO3_d_prev * fO3_h  # eq 11 fO3_d degrades during the day
    fO3_d_night = rO3_h * fO3_h  # eq 12  fO3_d recovers during the night
    if hr == 0:
        fO3_d = fO3_d_night
    else:
        fO3_d = fO3_d_base

    fO3_h_prev = fO3_h
    fO3_d_prev = fO3_d
    fO3_d_values.append(fO3_d)

plt.plot(fO3_d_values)


# %%

fO3_h_prev = 1
fO3_d_prev = 1
rO3_h = 1
fO3_d_values_alt = []
f_LA_values = []
for i in range(24 * t_l + 2):
    hr = i % 24
    # print(hr)
    dd = floor(i / 24)
    fO3_h = 0.9
    f_LA = 1 if dd < t_lem else 1 - (dd - t_lem) / t_lma if dd < t_l - 2 else 0

    rO3_h = fO3_d_prev + (1 - fO3_d_prev) * f_LA
    is_daylight = hr in [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    fO3_d_base = fO3_d_prev * fO3_h  # eq 11 fO3_d degrades during the day
    fO3_d_night = rO3_h * fO3_h  # eq 12  fO3_d recovers during the night
    fO3_d = fO3_d_base if is_daylight == 1 else fO3_d_night

    fO3_h_prev = fO3_h
    fO3_d_prev = fO3_d

    fO3_d_values_alt.append(fO3_d)

plt.plot(fO3_d_values_alt)


# %%
plt.plot(fO3_d_values, label="fO3_d_values")
plt.plot(fO3_d_values_alt, label="fO3_d_values_alt")
plt.legend()
# %%

plt.plot([1 if dd < t_lem else 1 - (dd - t_lem) / t_lma if dd < t_l else 0 for dd in range(100)])

# %%

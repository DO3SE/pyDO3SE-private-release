# %%
from do3se_met.wind import calc_ustar_and_L, ustar_from_velocity

u_h = 25

calc_ustar_and_L(
    1.6,
    102,
    102.3,
    2.08 + 271.15,
    3,
    u_h * 0.7,
    u_h * 0.1,
)

# %%

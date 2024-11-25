# %%
from math import sqrt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

# O3up = calc_fst(
#     Gsto_l=19.33,
#     Lm=0.01,  #?
#     uh=0.51,
#     Rsto_l=2120.15,
#     Rext=2500,
#     O3_nmol_m3=2291,
# )
# O3up

# %%


@interact(
    Gsto_l=widgets.FloatSlider(min=0.1, max=20000, value=19.33),
    Lm=widgets.FloatSlider(min=0.01, max=1, value=0.01),
    uh=widgets.FloatSlider(min=0.1, max=2, value=0.51),
    Rsto_l=widgets.FloatSlider(min=0.1, max=5000, value=2120),
    Rext=widgets.FloatSlider(min=0.1, max=5000, value=2500),
    O3_nmol_m3=widgets.FloatSlider(min=0.1, max=4000, value=2291),
)
def calc_fst(
    Gsto_l: float,
    Lm: float,
    uh: float,
    Rsto_l: float,
    Rext: float,
    O3_nmol_m3: float,
) -> float:
    if (Gsto_l > 0):
        leaf_rb = 1.3 * 150 * sqrt(Lm / uh)   # leaf boundary layer resistance (s/m)
        leaf_r = 1.0 / ((1.0 / Rsto_l) + (1.0 / Rext))  # leaf resistance in s/m
        Fst = O3_nmol_m3 * (1 / Rsto_l) * (leaf_r / (leaf_rb + leaf_r))
    else:
        Fst = 0
    return Fst

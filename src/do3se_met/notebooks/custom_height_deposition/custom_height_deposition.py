# %%
"""This notebook explores methods to get the ozone at a height that is not the canopy height."""

# %%
%load_ext autoreload
%autoreload 2
from do3se_met.deposition import (
    calc_canopy_ozone_concentration,
    calc_multi_layer_O3_ozone_concentration,
)
from do3se_met.resistance.resistance import calc_resistance_model

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
# %%
# Calculate single layer canopy top ozone
top_layer_index = 0
resistances = calc_resistance_model(
    nL=1,
    nLC=1,
    ustar_above_canopy=1.5,
    canopy_height=1.5,
    SAI_values=[[0]],
    LAI_values=[[0]],
    mean_gsto_values=[0],
    Rsoil=200,
    measured_height=0
)
print(resistances)
calc_canopy_ozone_concentration(
    O3_ppb_zR=100,
    Ra_ref_canopy=3,
    Ra_ref_measured=50,
    Ra_tar_canopy=resistances.Ra_canopy_to_izr,
    Ra_tar_canopy_top=resistances.Ra_canopy_top_to_izr,
    Rsur_ref=100001.998102272,
    Rsur_top_layer=resistances.Rsur[top_layer_index],
    Rb_ref=3.99,
    Rb_top_layer=resistances.Rb[top_layer_index],
)

# %%
nL = 20

y = calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [3 for _ in range(nL)],
    rm_Rsur = [3 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
)
x = np.arange(nL)
plt.scatter(
    x,
    y,
)
plt.xlabel("Layer")
plt.ylabel("Ozone concentration (ppb)")

# %%
nL = 8
plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [10,5,5,3,2,1,1,1],
    rm_Rsur = [10,5,5,3,2,1,1,1],
    rm_Rgs = 100  # Soil resistance is constant
), label="High resistance at top")
plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [1,1,1,1,2,3,4,10],
    rm_Rsur = [1,1,1,1,2,3,4,10],
    rm_Rgs = 100  # Soil resistance is constant
), label="High resistance at bottom")
plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [9999 for _ in range(nL)],
    rm_Rsur = [9999 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
), label="Uniform resistance")

plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rsur = [0.0000000001 for _ in range(nL)],
    rm_Rinc = [0.0000000001 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
), label="Uniform low resistance")

plt.legend()
plt.xlabel("Layer")
plt.ylabel("Ozone concentration (ppb)")
# Plot line at 0
plt.axhline(0, color='black', lw=1)
plt.title('Ozone concentration through canopy for different resistance profiles')
# %%
nL = 8
plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [1,0.5,0.5,0.3,0.2,0.1,0.1,0.1],
    rm_Rsur = [1,0.5,0.5,0.3,0.2,0.1,0.1,0.1],
    rm_Rgs = 100  # Soil resistance is constant
), label="High resistance at top * 0.1")


plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [10,5,5,3,2,1,1,1],
    rm_Rsur = [10,5,5,3,2,1,1,1],
    rm_Rgs = 100  # Soil resistance is constant
), label="High resistance at top")

plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [100,50,50,30,20,10,10,10],
    rm_Rsur = [100,50,50,30,20,10,10,10],
    rm_Rgs = 100  # Soil resistance is constant
), label="High resistance at top * 10")


plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [1000,500,500,300,200,100,100,100],
    rm_Rsur = [1000,500,500,300,200,100,100,100],
    rm_Rgs = 100  # Soil resistance is constant
), label="High resistance at top * 100")

plt.legend()
plt.xlabel("Layer")
plt.ylabel("Ozone concentration (ppb)")
# Plot line at 0
plt.axhline(0, color='black', lw=1)
plt.title('Ozone concentration through canopy for different resistances scaled')

# %%


nL = 8
values = [
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [10,5,5,3,2,1,1,1],
    rm_Rsur = [10,5,5,3,2,1,1,1],
    rm_Rgs = 100  # Soil resistance is constant
)[-1],
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [1,1,1,1,2,3,4,10],
    rm_Rsur = [1,1,1,1,2,3,4,10],
    rm_Rgs = 100  # Soil resistance is constant
)[-1],
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [9999 for _ in range(nL)],
    rm_Rsur = [9999 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
)[-1],


calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rsur = [0.0000000001 for _ in range(nL)],
    rm_Rinc = [0.0000000001 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
)[-1]]

labels = [
    "High resistance at top",
"High resistance at bottom",
"Uniform high resistance",
"Uniform low resistance",
]

# rotate x labels 90 degrees
plt.scatter(
    labels,
    values,
)
plt.xticks(rotation=90)
# plot h line at 0
plt.axhline(0, color='black', lw=1)
plt.title('Ozone concentration at bottom of canopy for different resistance profiles')

# %%


nL = 8
r_vals = [
    0.00001,0.1,1,10,100,1000,10000,100000
]
values = [
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=80,
    rm_Ra=0,
    rm_Rsur = [r for _ in range(nL)],
    rm_Rinc = [r for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
)[-1] for r in r_vals]

# labels = [
#     "High resistance at top",
# "High resistance at bottom",
# "Uniform resistance",
# "Uniform low resistance",
# ]

# rotate x labels 90 degrees
plt.scatter(
    [str(r) for r in r_vals],
    values,
)
plt.xticks(rotation=90)
# plot h line at 0
plt.axhline(0, color='black', lw=1)
plt.title('Ozone concentration at bottom of canopy for different constant resistances')
plt.xlabel("Resistance")
plt.ylabel('Ozone concentration (ppb)')
# %%
from scipy.linalg.lapack import sgesv


def calc_multi_layer_O3_ozone_concentration_update(
    nL: int,
    O3_in: float,
    rm_Ra: float,
    rm_Rinc,
    rm_Rsur,
    rm_Rgs: float,
):
    """Calculate O3 concentration for all layers.

    Requires that the value for the top layer (umet(1)%O3) is already known.

    # TODO: Is this absorbed O3 at each layer?

    Assumes layer 0 is top layer.

    This uses SGESV which is an old Fortran function. The documentation is vague on this.

    Parameters
    ----------
    nL: float
        number of model layers
    O3_in: float
        O3 for top layer micromet[ppb]
    rm_Ra: float
        rmodel_O3 Ra
    rm_Rinc: List[float]
        rmodel_O3 Rinc per layer
    rm_Rsur: List[float]
        rmodel_O3 Rsur per layer
    rm_Rgs: float
        rmodel_O3 Rgs

    Output
        O3 per layer[ppb]

    """
    C = np.full((nL + 1), 0, dtype=float)
    X = np.full((nL + 1, nL + 1), 0, dtype=float)

    bigR = np.array([rm_Ra] + rm_Rinc)
    assert bigR.shape == (nL + 1,)

    # TODO: per-layer Rb
    smallR = np.array(rm_Rsur + [rm_Rgs])
    assert smallR.shape == (nL + 1,)
    # Iterate over columns
    for j in range(0, nL + 1):
        X[0:j, j] = bigR[0:j]
        X[j, j] = X[j, j] + smallR[j]
        if j < nL:
            X[j + 1, j] = -smallR[j]


    C[0] = O3_in
    print(C)
    out = sgesv(X, C)
    lu, IPIV, C_out, info = sgesv(X, C)

    if info != 0:
        raise Exception('SGESV Failed')
    C_final = smallR * C_out
    O3_out = C_final[0: nL]
    return O3_out

plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration_update(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rsur = [1000 for _ in range(nL)],
    rm_Rinc = [1000 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
), label="High resistance at top * 100")

plt.legend()
plt.xlabel("Layer")
plt.ylabel("Ozone concentration (ppb)")
# Plot line at 0
plt.axhline(0, color='black', lw=1)
plt.title('Ozone concentration through canopy for different resistances scaled')
# %%
# import ipywidgets as widgets
import ipywidgets as widgets
from ipywidgets import interact, interact_manual

@interact
def plot_multi_layer_ozone(
    O3_in=widgets.FloatSlider(min=0, max=100, step=1, value=20),
    rm_Ra=widgets.FloatSlider(min=0, max=100, step=1, value=0),
    rm_Rinc=widgets.FloatSlider(min=1e-10, max=10000, step=0.01, value=600),
    rm_Rsur=widgets.FloatSlider(min=1e-10, max=50000, step=0.01, value=10000),
    rm_Rgs=widgets.FloatSlider(min=1e-10, max=2000, step=0.1, value=100),
):
    plt.scatter(
    np.arange(nL),
        calc_multi_layer_O3_ozone_concentration(
            nL=nL,
            O3_in=20,
            rm_Ra=0,
            rm_Rinc = [10,5,5,3,2,1,1,1],
            rm_Rsur = [10,5,5,3,2,1,1,1],
            rm_Rgs = 100  # Soil resistance is constant
        ), label="baseline high at top")
    plt.scatter(
        np.arange(nL),
    calc_multi_layer_O3_ozone_concentration(
        nL=nL,
        O3_in=O3_in,
        rm_Ra=rm_Ra,
        rm_Rinc = [rm_Rinc for _ in range(nL)],
        rm_Rsur = [rm_Rsur for _ in range(nL)],
        rm_Rgs = rm_Rgs  # Soil resistance is constant
    ), label="Ozone concentration")

    plt.legend()
    plt.xlabel("Layer")
    plt.ylabel("Ozone concentration (ppb)")
    # Plot line at 0
    plt.axhline(0, color='black', lw=1)
    plt.title('Ozone concentration through canopy for different resistances scaled')



# %%

# %%
nL = 8
plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [999999 for _ in range(nL)],
    rm_Rsur = [100000 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
), label="100")


plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [0.1 for _ in range(nL)],
    rm_Rsur = [0.1 for _ in range(nL)],
    rm_Rgs = 100  # Soil resistance is constant
), label="100")


plt.scatter(
    np.arange(nL),
calc_multi_layer_O3_ozone_concentration(
    nL=nL,
    O3_in=20,
    rm_Ra=0,
    rm_Rinc = [999999 for _ in range(nL)],
    rm_Rsur = [999999 for _ in range(nL)],
    rm_Rgs = 800  # Soil resistance is constant
), label="800")



plt.legend()
plt.xlabel("Layer")
plt.ylabel("Ozone concentration (ppb)")
# Plot line at 0
plt.axhline(0, color='black', lw=1)
plt.title('Ozone concentration through canopy for different resistance profiles')

# %%

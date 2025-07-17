# %%
from typing import List


def calc_Rtotal(
    nL: int,
    Rsur: List[float],
    Rinc: List[float],
    Rgs: float,
) -> List[float]:
    """Calculate multi-layer Rtotal.

    The total resistance for each layer and everything below that layer.

    The equation below will tend towards the lowest of Rsur or Rinc+prev_layer_total.
    If the previous layer is 0 then the resistance is equal to Rsur + Rinc

    Parameters
    ----------
    nL: int
        Number of model layers
    Rsur: List[float]
        Combined surface resistance [s m-1] per layer
    Rinc: List[float]
        In-canopy aerodynamic resistance [s m-1] per layer
    Rgs: float
        Ground surface resistance [s m-1]

    """
    tmp = [None for _ in range(nL + 1)]
    tmp[-1] = Rgs
    for i in range(nL-1, -1, -1):
        print(i)
        print(tmp[i+1])
        tmp[i] = 1 / (1 / Rsur[i] + 1 / (Rinc[i] + tmp[i+1]))
    Rtotal = tmp[0:nL]
    return Rtotal

nL= 4
Rsur = [100 for _ in range(nL)]
Rinc = [3 for _ in range(nL)]
Rgs = 200
calc_Rtotal(nL, Rsur, Rinc, Rgs)



# %%
def calc_Rtotal_reversed(
    nL: int,
    Rsur: List[float],
    Rinc: List[float],
    Rgs: float,
) -> List[float]:
    """Calculate multi-layer Rtotal.

    The total resistance for each layer and everything below that layer.

    The equation below will tend towards the lowest of Rsur or Rinc+prev_layer_total.
    If the previous layer is 0 then the resistance is equal to Rsur + Rinc

    NOTE: Layer 0 is top

    Parameters
    ----------
    nL: int
        Number of model layers
    Rsur: List[float]
        Combined surface resistance [s m-1] per layer
    Rinc: List[float]
        In-canopy aerodynamic resistance [s m-1] per layer
    Rgs: float
        Ground surface resistance [s m-1]

    """
    tmp = [None for _ in range(nL + 1)]
    tmp[0] = Rgs
    for i in range(1, nL+1):
        print(i)
        print(tmp[i-1])
        tmp[i] = 1 / (1 / Rsur[i-1] + 1 / (Rinc[i-1] + tmp[i-1]))
    Rtotal = tmp[1:nL+1]
    return Rtotal

nL= 4
Rsur = [100 for _ in range(nL)]
Rinc = [3 for _ in range(nL)]
Rgs = 200
calc_Rtotal_reversed(nL, Rsur, Rinc, Rgs)


# %%


Rsur = [100 for _ in range(nL)]
Rinc = [3 for _ in range(nL)]
Rgs = 200
calc_Rtotal_reversed(
    nL=3,
    Rsur=[100, 100, 100],
    Rinc=[0.1, 0.1, 0.1],
    Rgs=333,
)


# %%

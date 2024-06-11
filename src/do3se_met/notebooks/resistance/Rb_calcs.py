# %%
from do3se_met.resistance import calc_Rb
from matplotlib import pyplot as plt

from do3se_met.wind import calc_monin_obukhov_length
# %%

CANOPY_D = 0.7
CANOPY_Z0 = 0.1
canopy_height = 1
ustar_list = [i/10 for i in range(1,10)]
z1 = (canopy_height * (CANOPY_D + CANOPY_Z0))
z2 = 4
d = canopy_height * CANOPY_D
diff = 1
Ra_simple = [calc_Rb(ustar, diff) for ustar in ustar_list]

plt.plot(Ra_simple, label='Ra_simple')

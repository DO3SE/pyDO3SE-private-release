# %%
from do3se_met.resistance import calc_Ra_simple, calc_Ra_with_heat_flux
from matplotlib import pyplot as plt

from do3se_met.wind import calc_monin_obukhov_length
# %%
CANOPY_D = 0.7
CANOPY_Z0 = 0.1
canopy_height = 1
ustar = 1.1
z1 = (canopy_height * (CANOPY_D + CANOPY_Z0))
z2_list = list(range(0,10))
d = canopy_height * CANOPY_D

Ra_simple = [calc_Ra_simple(ustar, z1, z2, d) for z2 in z2_list]

plt.plot(Ra_simple, label='Ra_simple')

Tk = 0 + 273.15
Hd = -11
P = 101
L = calc_monin_obukhov_length(Tk, ustar, Hd, P)
Ra_with_heat_flux = [calc_Ra_with_heat_flux(ustar, z1, z2, L) for z2 in z2_list]

plt.plot(Ra_with_heat_flux, label='Ra_with_heat_flux')
plt.legend()


# %%
from functools import partial
CANOPY_D = 0.7
CANOPY_Z0 = 0.1
canopy_height = 1
ustar_list = [i/10 for i in range(1,10)]
z1 = (canopy_height * (CANOPY_D + CANOPY_Z0))
z2 = 4
d = canopy_height * CANOPY_D

Ra_simple = [calc_Ra_simple(ustar, z1, z2, d) for ustar in ustar_list]

plt.plot(Ra_simple, label='Ra_simple')

Tk = 0 + 273.15
Hd = -11
P = 101
L = partial(calc_monin_obukhov_length, Tk=Tk, Hd=Hd,P=P)
Ra_with_heat_flux = [calc_Ra_with_heat_flux(ustar, z1, z2, L(ustar=ustar)) for ustar in ustar_list]

plt.plot(Ra_with_heat_flux, label='Ra_with_heat_flux')
plt.legend()
# RMODEL H2O CALCS
# %%
from math import inf
Rinc = [1]
nL = 1
Rinc_sum = sum([(1.0 / r) if r > 0 else inf for r in Rinc[0:nL]])
Rinc = 1.0 / Rinc_sum if Rinc_sum > 0 else inf
Rinc
# %%
# We need to work out these values
# I(state.canopy.rmodel_H2O.Rb, as_='rm_Rb'),   => Rb_H2O
# I(state.canopy.rmodel_H2O.Rinc[0], as_='rm_Rinc_l0'), => Rinc
# I(state.canopy.rmodel_H2O.Rsto[0], as_='rm_Rsto_l0'), => Rsto_c
# I(state.canopy.rmodel_H2O.Rgs, as_='rm_Rgs'), Rsoil

# Python Calcs take rmodel_O3 as input

# rb_H2O  === OK
# FORTRAN UI
Rb_H2O = rb(ustar, DH2O) = (2.0 / (K * ustar)) * (((V / d) / PR)**(2.0 / 3.0))
# Python
Rb_H2O = calc_Rb(ustar, DIFF_H2O)

# Rinc == OK
# FORTRAN UI
Rinc_b = 14
Rinc = Calc_Rinc() = Rinc_b * SAI * h / ustar
# Python
Rinc = Rinc_b * SAI * h / ustar if ustar > 0 else MAX_RINC
Rinc_sum = sum([(1.0 / r) if r > 0 else inf for r in rmodel.Rinc[0:rmodel.nL]])
Rinc = 1.0 / Rinc_sum if Rinc_sum > 0 else inf

# Rsto_c
# FORTRAN UI
Gsto_c = Gsto * LAI
Rsto_c = rsto_from_gsto(Gsto_c) = min(MAX_RSTO, 41000.0 / gsto)
# PYTHON
# TODO: Check DRATIO conversion needed
Rsto = DRATIO * 1.0 / sum([(1.0 / r) for r in rmodel.Rsto[0:rmodel.nL]])

# Rsoil === OK
Rsoil = 100

# %%
# CHECK O3 RESISTANCE MODEL SETUP
# We need to work out these values
# rmodel.Rb, # OK
# rmodel.Rinc[0], #OK
# rmodel.Rsto[0],
# rmodel.Rgs, # OK

# Calc_O3_Concentration
# FORTRAN
# %%

# ustar over reference canopy
ustar_ref = estimate_ustar(uh_50, izR - O3_d, O3_zo)  # OK
# Ra between reference canopy and izR
Ra_ref_50 = ra_simple(ustar_ref, O3_zo + O3_d, izR, O3_d)  # OK
# Rb for reference canopy
Rb_ref_50 = rb_func(ustar_ref, DO3)  # OK
# Deposition velocity at izR over reference canopy
# (assuming that Rsur_ref = Rsur)
# TODO: Check Ra_c = Ra_ref_50
# TODO: Check rmodel_total_top_layer = rb + rsur
Vd_50 = 1.0 / (Ra_ref_50 + Rb_ref_50 + Rsur)
# Ra between measurement height and izR
Ra_O3zR_50 = ra_simple(ustar_ref, O3zR, izR, O3_d)  # OK
# O3 concentration at izR
O3_50 = O3_ppb_zR / (1.0 - (Ra_O3zR_50 * Vd_50))  # OK
# Ra between target canopy and izR
# (ustar already calculated for target canopy)
Ra_tar_50 = ra_simple(ustar, zo + d, izR, d)
# Deposition velocity at izR over target canopy
Vd = 1.0 / (Ra_tar_50 + Rb + Rsur)
# O3 concentration at target canopy
# (Ra already calculated between canopy height and izR)
O3_ppb = O3_50 * (1.0 - (Ra * Vd))

# Specific molar volume of an ideal gas at current temp + pressure
Vn = 8.314510 * ((Ts_C + Ts_K) / P)
# Convert to nmol/m^3
O3_nmol_m3 = (1.0 / Vn) * O3_ppb * M_O3 * 20.833  # 1 microgram O3 = 20.833 nmol/m^3

######################################
# FORTRAN IMPORT
# DO3SE_PhysicalConstants_ml.F90
######################################

# TODO: Use math library instead
# PI: float = 3.141592653589793238  # real, parameter
# DEG2RAD: float = 0.017453292519943295  # real, parameter

# > Atmospheric pressure at sea level [kPa]
seaP: float = 101.325  # real, parameter

# > 0 degrees Celsius in Kelvin
T0: float = 273.15  # real, parameter
# > Stefan-Boltzmann constant (W m-2 K-4)
SBC: float = 5.670373e-8  # real, parameter

# > Universal gas constant (J K-1 mol-1)
# TODO: update this value to 8.3144621
R: float = 8.314472  # real, parameter
Rmass: float = 287.0   # mass gas constant for dry air (J/Kg/K)
# > Approximate fraction of global radiation in PAR waveband
PARfrac: float = 0.45  # real, parameter
# > PAR conversion from W m-2 to umol photons m-2 s-1
PAR_Wm2_to_photons: float = 4.57  # real, parameter
# > Net radiation conversion from MJ m-2 h-1 to W m-2
Rn_MJ_to_W: float = 277.8  # real, parameter

# > Molecular diffusivity of O3 in air (m2 s-1)
DIFF_O3: float = 0.000015  # real, parameter
# > Molecular diffusivity of H2O (m2 s-1)
DIFF_H2O: float = 0.000025  # real, parameter

# Ratios from Massman paper (1998)
# TODO: Check for still air & semi turb air
# Still for gsto and turb for g_bl
# > Ratio between molecular diffusivity of O3 and H2O at 20degc
DRATIO_O3_H20: float = 0.663  # real, parameter DIFF_O3 / DIFF_H2O
DRATIO_O3_CO2: float = 1.043  # real, parameter DIFF_O3 / DIFF_CO2
DRATIO_H20_O3: float = 1.508  # real, parameter DIFF_H2O / DIFF_O3


VON_KAR: float = 0.41  # von Karman's constant

GSC: float = 0.082  # Solar constant (MJ/m^2/min)
# TODO: This is wrong value  \/\/
# SBC: float = 4.903e-9 / 24  # Stephan Boltzman constant

g: float = 9.8          # gravitational acceleration
cp: float = 1005        # specific head of air at constant

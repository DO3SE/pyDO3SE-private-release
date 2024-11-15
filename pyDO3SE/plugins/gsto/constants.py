# Defined in Farquhar model! parameters considered (or defined) to be constant for all species
O_i = 210.0           # O2 concentration                   [mmol/mol]
E_K_C = 79430.0       # activation energy of K_C           [J/mol]            Medlyn2002
E_K_O = 36380.0       # activation energy of K_O           [J/mol]            Medlyn2002
E_Gamma_star = 37830.0  # activation energy for C-comp-point [J/mol]            Medlyn2002
K_C_25 = 404.9        # K.C at reference temperature 25    [micro mol/mol]    Medlyn2002
K_O_25 = 278.4        # K.O at reference temperature 25    [mmol/mol]         Medlyn2002
Gamma_star_25 = 42.75  # CO2 compensation point at T= 25    [micro mol/mol]    Medlyn2002
A_j_a = 4.0           # electron requirement for NADPH formation
A_j_b = 8.0           # electron requirement for ATP formation

#  species spedific model parameters (that don't tend to have species specific
#  values, others are supplied as arguments)
alpha = 0.3           # efficiency light energy conversion [mol electrons/mol photons]
Teta = 0.90           # shape of J~Q determining factor    []
H_a_jmax = 50300      # activation energy for J_max        [J/mol]
H_d_jmax = 152044     # deactivation energy for J_max      [J/mol]
H_a_vcmax = 73637     # activation energy for V_cmax       [J/mol]
H_d_vcmax = 149252    # deactivation energy for V_cmax     [J/mol]
S_V_vcmax = 486       # entropy terms                      [J/(mol*K)]
S_V_jmax = 495        # entropy terms                      [J/(mol*K)
fDO3 = 0.93
Gamma_3 = 0.5

# parameters considered (or defined) to be constant for all species
A_j_a = 4.0  # electron requirement for NADPH formation
A_j_b = 8.0  # electron requirement for ATP formation

# species spedific model paarameters (that don't tend to have species specific
# values, others are supplied as arguments)
fDO3 = 0.93
Gamma_3 = 0.5

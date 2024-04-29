# %%
# IMPORTS
from math import exp, pi
from scipy.integrate import quad as integrate

from pyDO3SE.plugins.met.helpers import sunlit_LAI

from pyDO3SE.plugins.gsto.ewert.ewert_helpers import Ewert_Input_Factors, calc_input_factors

from pyDO3SE.plugins.gsto.ewert.ewert import CO2_Constant_Loop_Inputs, CO2_loop_State, \
    co2_concentration_in_stomata_loop, Output_Shape

# %%
# Model coefficients
sigma = 0.15  # leaf scattering coefficient of PAR (p_i + T_i)
k_d_alt = 0.719  # diffuse and scattered diffuse PAR extinction coefficient, 0.719
seaP = 101.325  # real, parameter

# %%
# Equations ======================================


def calc_Idrctt_Idfuse(
    PAR: float,
    sinB: float,
    P: float,
) -> float:
    """Estimate diffuse and direct PAR components.

    Args:
    PAR: float       !< Photosynthetically active radiation [W m-2]
    sinB: float      !< sin() of solar elevation angle
    P: float         !< Atmospheric pressure [kPa]

    Returns:
    Idrctt: float   !< Direct PAR irradiance [W m-2]
    Idfuse: float   !< Diffuse PAR irradiance [W m-2]

    real :: m, pPARdir, pPARdif, pPARtotal, ST, fPARdir, fPARdif
    """

    if sinB > 0.0:
        m = 1.0 / sinB

        # Potential direct PAR
        pPARdir = 600 * exp(-0.185 * (P / seaP) * m) * sinB
        # Potential diffuse PAR
        pPARdif = 0.4 * (600 - pPARdir) * sinB
        # Potential total PAR
        pPARtotal = pPARdir + pPARdif

        # Sky transmissivity
        ST = max(0.21, min(0.9, PAR / pPARtotal))

        # Direct and diffuse fractions
        fPARdir = (pPARdir / pPARtotal) * (1.0 - ((0.9 - ST) / 0.7)**(2.0/3.0))
        fPARdif = 1 - fPARdir

        # Apply calculated direct and diffuse fractions to PARtotal
        Idrctt = fPARdir * PAR
        Idfuse = fPARdif * PAR
    else:
        Idrctt = 0.0
        Idfuse = 0.0

    return Idrctt, Idfuse


def calc_diffuse_irradiance_refl(
    # P_cb: float,  # beam irradiance for uniform leaf angle distribution
    sinB: float,
    P_h: float,
    # N_d: float,  #
    Ir_dfuse_0: float,  # Diffuse PAR per unit ground area at top of canopy[umol m^-2 s^-2]
) -> float:
    """Calculate diffuse irradiance reflection coefficient EQ A21."""
    # f
    def f(alpha):
        """Calculate integral per rotation unit of sun[radian]"""
        # TODO: Check this is correct
        # We should be
        k_b = 0.5 / sinB
        # TODO: This probably varies based on leaf normal
        # k_b = 0.5 / (sinB * cos(alpha))
        # Diffuse photon radiance of the sky (per radian??)[umol m^-2 sr ^-1]
        N_d = Ir_dfuse_0 / (2 * pi)
        P_cb = 1 - exp(2 * P_h * k_b / (1 + k_b))
        return N_d * P_cb

    integ, err = integrate(f, 0, pi/2)
    P_cd = (1 / Ir_dfuse_0) * integ
    return P_cd


def calc_PAR_sun_shade(
    Ir_beam_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    Ir_dfuse_0: float,  # PAR per unit ground area at top of canopy[umol m^-2 s^-2]
    sinB: float,
    cosA: float,
    LAI_c: float,  # cumulative leaf-area index from top of canopy (L = 0 at top)
):
    """Calculate the sun and shade PAR values.

    Uses equations from Farquhar 1997 (Simple scaling of photosynthesis from leaves to canopies
    without the errors of big-leaf models)
    """
    k_b = 0.5 / sinB  # Beam radiation extinction coefficient
    k_b_alt = 0.46 / sinB  # beam and scattered beam PAR extinction coefficient
    P_h = (1 - (1 - sigma) ** 0.5) / (1 + (1 - sigma) ** 0.5)
    P_cb = 1 - exp(2 * P_h * k_b / (1 + k_b))
    P_cd = calc_diffuse_irradiance_refl(sinB, P_h, Ir_dfuse_0)
    Ir_diffuse = (1 - P_cd) * k_d_alt * Ir_dfuse_0 * exp(-k_d_alt * LAI_c)
    Ir_beam_sun = (1 - sigma) * Ir_beam_0 * cosA / sinB
    Ir_scattered_b = Ir_beam_0 * (
        ((1 - P_cb)*k_b_alt * exp(-k_b * LAI_c)) /
        -(1 - sigma) * k_b * exp(-k_b * LAI_c)
    )
    PAR_shade = Ir_diffuse + Ir_scattered_b
    PAR_sun = PAR_shade + Ir_beam_sun
    return PAR_sun, PAR_shade


# %%
# Ewert Model Inputs
t_lem_constant = 0.15
t_lse_constant = 0.33

t_l = 800
t_lem = t_l * t_lem_constant
t_lma = t_l - t_lem
t_lse = t_l - (t_lma * t_lse_constant)
t_lep = t_l - (t_lem + t_lse)  # check this


Tleaf_C = 20.0
PARsun = 800
PARshade = 400
V_cmax_25 = 180.0
J_max_25 = 400.0
LAI = 0.01
sinB = 0.3  # TODO: Check in radians

t_lse_constant = t_lse_constant
t_l_estimate = t_l
t_lem = t_lem
t_lep = t_lep
t_lse = t_lse
t_lma = t_lma

c_a = 391.0
eact = 1.0
g_bv = 1469999.0
g_sto_0 = 20000
m = 8.12
D_0 = 2.27
O3_nmol = 14179
O3up_prev = 80
O3up_acc_prev = 300
Lm = 0.02
uh = 0.4
# Rsto_l=82.7398147583  # We now calculate this internally
Rext = 2500  # TODO: Check this value!
fO3_d_prev = 0.89
td_dd = 24.1
td_dd_prev = 11.3
gamma_1 = 0.06
gamma_2 = 0.0045
gamma_3 = 0.5
is_daylight = True

# %%
# ===================RUN MODEL OLD SETUP==================== #
# ========================================================== #
# ========================================================== #
# ========================================================== #
# ========================================================== #

e_a = eact * 1e3

# %%
# Setup farquhar for PARSun

inputs_1 = calc_input_factors(
    Tleaf_C=Tleaf_C,
    Q=PARsun * 4.57,
    V_cmax_25=V_cmax_25,
    J_max_25=J_max_25,
)

loop_inputs_1 = CO2_Constant_Loop_Inputs(
    c_a=c_a,
    e_a=e_a,
    g_bl=g_bv,
    g_sto_0=g_sto_0,
    m=m,
    D_0=D_0 * 1e3,
    O3_nmol=O3_nmol,
    O3up_prev=O3up_prev,
    O3up_acc_prev=O3up_acc_prev,
    Lm=Lm,
    uh=uh,
    # Rsto_l=Rsto_l,  # We now calculate this internally
    Rext=Rext,
    fO3_d_prev=fO3_d_prev,
    td_dd=td_dd,
    td_dd_prev=td_dd_prev,
    gamma_1=gamma_1,
    gamma_2=gamma_2,
    gamma_3=gamma_3,
    is_daylight=is_daylight,
    t_lse_constant=t_lse_constant,
    t_l_estimate=t_l_estimate,
    t_lem=t_lem,
    t_lep=t_lep,
    t_lse=t_lse,
    t_lma=t_lma,
    Gamma=inputs_1.Gamma,
    Gamma_star=inputs_1.Gamma_star,
    V_cmax=inputs_1.V_cmax,
    K_C=inputs_1.K_C,
    K_O=inputs_1.K_O,
    J=inputs_1.J,
    R_d=inputs_1.R_d,
    e_sat_i=inputs_1.e_sat_i,
)

# %%
# Run farquhar for PARSun
state_out_1 = co2_concentration_in_stomata_loop(loop_inputs_1)
state_out_1
# %%
# Setup farquhar for PARshade

inputs_2: Ewert_Input_Factors = calc_input_factors(
    Tleaf_C=Tleaf_C,
    Q=PARshade * 4.57,
    V_cmax_25=V_cmax_25,
    J_max_25=J_max_25,
)

loop_inputs_2: CO2_Constant_Loop_Inputs = loop_inputs_1._replace(
    Gamma=inputs_2.Gamma,
    Gamma_star=inputs_2.Gamma_star,
    V_cmax=inputs_2.V_cmax,
    K_C=inputs_2.K_C,
    K_O=inputs_2.K_O,
    J=inputs_2.J,
    R_d=inputs_2.R_d,
    e_sat_i=inputs_2.e_sat_i,
)

# %%
# Run farquhar for PARshade
state_out_2 = co2_concentration_in_stomata_loop(loop_inputs_2)
state_out_2
# %%

# setup farquhar for PAR

inputs_3: Ewert_Input_Factors = calc_input_factors(
    Tleaf_C=Tleaf_C,
    Q=PARshade * 4.57,
    V_cmax_25=V_cmax_25,
    J_max_25=J_max_25,
)

loop_inputs_3: CO2_Constant_Loop_Inputs = loop_inputs_2._replace(
    Gamma=inputs_3.Gamma,
    Gamma_star=inputs_3.Gamma_star,
    V_cmax=inputs_3.V_cmax,
    K_C=inputs_3.K_C,
    K_O=inputs_3.K_O,
    J=inputs_3.J,
    R_d=inputs_3.R_d,
    e_sat_i=inputs_3.e_sat_i,
)

# %%
# run farquhar for PAR
state_out_3 = co2_concentration_in_stomata_loop(loop_inputs_3)
state_out_3
# %%
# Ewert output

# Get canopy A_n
A_n_sun = state_out_1.A_n
A_n_shade = state_out_2.A_n

# Farquhar 1997
LAIsun = sunlit_LAI(LAI, sinB)
LAIshade = LAI - LAIsun

Canopy_A_n = LAIsun * A_n_sun + LAIshade * A_n_shade

result = Output_Shape(
    Tleaf_C=Tleaf_C,
    g_sv=state_out_3.g_sto,
    f_LS=state_out_3.f_LS,
    A_n=state_out_3.A_n,
    A_c=state_out_3.A_c,
    A_j=state_out_3.A_j,
    A_p=state_out_3.A_p,
    A_n_limit_factor=state_out_3.A_n_limit_factor,
    R_d=inputs_3.R_d,
    O3up_out=state_out_3.O3up,
    O3up_acc_out=state_out_3.O3up_acc,
    fO3_h_out=state_out_3.fO3_h,
    fO3_d_out=state_out_3.fO3_d,
    c_i=state_out_3.c_i,
    Canopy_A_n=Canopy_A_n,
    Rsto_l=state_out_3.Rsto_l,
)
result

# %%

# ================= Run model (NEW SETUP) ================== #
# ========================================================== #
# ========================================================== #
# ========================================================== #

# TODO: Set these as sliders
PAR = 800 * 4.57  # PAR converted to [umol/m^2/s]
P = 98.1
cosA = 0.5
LAI_c = LAI
# New PAR sun shade equations
Ir_beam_0, Ir_dfuse_0 = calc_Idrctt_Idfuse(PAR, sinB, P)
PARsun, PARshade = calc_PAR_sun_shade(Ir_beam_0, Ir_dfuse_0, sinB, cosA, LAI_c)
k_b = 0.5 / sinB
# f_sun = exp(-k_b * LAI_c)
LAIsun = ((1 - exp(-k_b * LAI)) / k_b)
LAIshade = LAI - LAIsun
# ======
e_a = eact * 1e3

# %%
# Setup farquhar for PARSun
inputs_1: Ewert_Input_Factors = calc_input_factors(
    Tleaf_C=Tleaf_C,
    Q=PARsun,  # * 4.57,  already converted to umol above
    V_cmax_25=V_cmax_25,
    J_max_25=J_max_25,
)

loop_inputs_1: CO2_Constant_Loop_Inputs = CO2_Constant_Loop_Inputs(
    c_a=c_a,
    e_a=e_a,
    g_bl=g_bv,
    g_sto_0=g_sto_0,
    m=m,
    D_0=D_0 * 1e3,
    O3_nmol=O3_nmol,
    O3up_prev=O3up_prev,
    O3up_acc_prev=O3up_acc_prev,
    Lm=Lm,
    uh=uh,
    # Rsto_l=Rsto_l,
    Rext=Rext,
    fO3_d_prev=fO3_d_prev,
    td_dd=td_dd,
    td_dd_prev=td_dd_prev,
    gamma_1=gamma_1,
    gamma_2=gamma_2,
    gamma_3=gamma_3,
    is_daylight=is_daylight,
    t_lse_constant=t_lse_constant,
    t_l_estimate=t_l_estimate,
    t_lem=t_lem,
    t_lep=t_lep,
    t_lse=t_lse,
    t_lma=t_lma,
    Gamma=inputs_1.Gamma,
    Gamma_star=inputs_1.Gamma_star,
    V_cmax=inputs_1.V_cmax,
    K_C=inputs_1.K_C,
    K_O=inputs_1.K_O,
    J=inputs_1.J,
    R_d=inputs_1.R_d,
    e_sat_i=inputs_1.e_sat_i,
)
loop_inputs_1
# %%
# Run farquhar for PARSun
state_out_1 = co2_concentration_in_stomata_loop(loop_inputs_1)
state_out_1
# %%
# Setup farquhar for PARshade

inputs_2: Ewert_Input_Factors = calc_input_factors(
    Tleaf_C=Tleaf_C,
    Q=PARshade,  # * 4.57,  already converted to umol above
    V_cmax_25=V_cmax_25,
    J_max_25=J_max_25,
)

loop_inputs_2: CO2_Constant_Loop_Inputs = loop_inputs_1._replace(
    Gamma=inputs_2.Gamma,
    Gamma_star=inputs_2.Gamma_star,
    V_cmax=inputs_2.V_cmax,
    K_C=inputs_2.K_C,
    K_O=inputs_2.K_O,
    J=inputs_2.J,
    R_d=inputs_2.R_d,
    e_sat_i=inputs_2.e_sat_i,
)
loop_inputs_2
# %%
# Run farquhar for PARShade

state_out_2 = co2_concentration_in_stomata_loop(loop_inputs_2)
state_out_2

# %%
# Setup farquhar for PAR

inputs_3: Ewert_Input_Factors = calc_input_factors(
    Tleaf_C=Tleaf_C,
    Q=PAR,  # Use unsplit PAR here
    V_cmax_25=V_cmax_25,
    J_max_25=J_max_25,
)

loop_inputs_3: CO2_Constant_Loop_Inputs = loop_inputs_2._replace(
    Gamma=inputs_3.Gamma,
    Gamma_star=inputs_3.Gamma_star,
    V_cmax=inputs_3.V_cmax,
    K_C=inputs_3.K_C,
    K_O=inputs_3.K_O,
    J=inputs_3.J,
    R_d=inputs_3.R_d,
    e_sat_i=inputs_3.e_sat_i,
)

# %%
# Run farquhar for PAR

state_out_3 = co2_concentration_in_stomata_loop(loop_inputs_3)
state_out_3

# %%
# Get ewert output

# Get canopy A_n
A_n_sun = state_out_1.A_n
A_n_shade = state_out_2.A_n

Canopy_A_n = LAIsun * A_n_sun + LAIshade * A_n_shade

result = Output_Shape(
    Tleaf_C=Tleaf_C,
    g_sv=state_out_3.g_sto,
    f_LS=state_out_3.f_LS,
    A_n=state_out_3.A_n,
    A_c=state_out_3.A_c,
    A_j=state_out_3.A_j,
    A_p=state_out_3.A_p,
    A_n_limit_factor=state_out_3.A_n_limit_factor,
    R_d=inputs_3.R_d,
    O3up_out=state_out_3.O3up,
    O3up_acc_out=state_out_3.O3up_acc,
    fO3_h_out=state_out_3.fO3_h,
    fO3_d_out=state_out_3.fO3_d,
    c_i=state_out_3.c_i,
    Canopy_A_n=Canopy_A_n,
    Rsto_l=state_out_3.Rsto_l,
)
result

# %%

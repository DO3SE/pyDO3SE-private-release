# %%
from math import isclose
from pyDO3SE.constants.physical_constants import DRATIO
from math import exp


def penman_monteith_daily(
    LAI: float,
    root_depth: float,
    run_off_fraction: float,
    ASW: float,
    SMD: float,
    pm_state_precip_acc: float,
    pm_state_run_off_acc: float,
    pm_state_Ei_acc: float,
    pm_state_Eat_acc: float,
    pm_state_percolated_acc: float,
):

    # real :: max_ET, delta_SM

    # Start with full amount of precipitation
    # Estimate loss to run-off
    rain_input = pm_state_precip_acc
    run_off = run_off_fraction * rain_input
    run_off_acc = pm_state_run_off_acc + run_off

    # Calculate "effective irrigation"
    effective_irrig = rain_input - run_off

    # Estimate loss of intercepted precipitation to evaporation.  Intercepted
    # precipitation is estimated as 0.0001*LAI, which is therefore a limit on
    # how much can be evaporated.
    intercepted_evaporated = min(effective_irrig, 0.0001 * LAI, pm_state_Ei_acc)

    # Can't lose water below PWP, so constrain evapotranspiration to ASW by
    # restricting evapotranspiration.
    max_ET = ASW + effective_irrig - intercepted_evaporated
    evapotranspiration = min(max_ET, pm_state_Eat_acc)

    # Total balance = input - run_off - evaporated - evapotranspiration
    # precip_acc - Ei - AEt
    # effective_irrig is ok
    # TODO: Ei is different
    # TODO: Aet is different
    delta_SM = effective_irrig - intercepted_evaporated - evapotranspiration
    # Converted to volumetric change using root depth.
    Sn_diff = delta_SM / root_depth

    # Amount that will go to deep percolation = remainder if water balance
    # refills soil water, i.e. if it is greater than SMD.
    percolated = max(0.0, delta_SM - SMD)
    percolated_acc = pm_state_percolated_acc + percolated
    return Sn_diff


out = penman_monteith_daily(
    LAI=3.0,
    root_depth=0.75,
    run_off_fraction=0.0,
    ASW=0.12388163805,
    SMD=0.0,
    pm_state_precip_acc=0.000410999957239,
    pm_state_run_off_acc=99.9,
    pm_state_Ei_acc=4.26064944E-03,
    pm_state_Eat_acc=7.86058372E-04,
    pm_state_percolated_acc=99.9,
)
print(out)
# %%
-5.00077906E-04

# %%


def Fortran_calc_swp(
    precip_acc,
    LAI,
    Ei,
    AEt,
    root,
):

    if (precip_acc > 0):
        P_input = (precip_acc - (0.0001 * LAI)) + ((0.0001 * LAI) - min(Ei, 0.0001 * LAI))
    else:
        P_input = 0
    # Can't lose water through Ei
    P_input = max(0.0, P_input)
    Sn_diff = (P_input - AEt) / root
    return Sn_diff
    # # Calculate new Sn, with field capacity as a maximum
    # Sn = min(Fc_m, Sn + Sn_diff)
    # Sn = max(PWP_vol, Sn)
    # per_vol = Sn * 100

    # # Calculate ASW and SWP for new water content
    # ASW = (Sn - PWP_vol) * root
    # SWP = SWP_AE * ((SWC_sat / Sn)**soil_b)

    # # Calculate SMD for new water content
    # SMD = (Fc_m - Sn) * root


Fortran_calc_swp(
    precip_acc=0.000410999957239,
    LAI=3.0,
    Ei=4.26064944E-03,
    AEt=7.86058372E-04,
    root=0.75,
)


# =========================================
# %%
Ts_K = 273.15


def FORtran_pen_mon(
    Rn_MJ,
    VPD,
    LAI,
    P,
    esat_kPa,
    eact_kPa,
    Ts_c,
    Rb_H2O,
    Rinc,
    Rsoil,
    Rsto_c,
    Es_blocked=False,
    Et_hr=0.0,
    Ei_dd=0.0,
    AEt_dd=0.0,
):
    Rn = Rn_MJ * 1000000.0

    VPD_Pa = VPD * 1000
    P_Pa = P * 1000

    esat = 1000 * esat_kPa
    eact = 1000 * eact_kPa

    Tvir = (Ts_c + Ts_K) / (1 - (0.378 * (eact / P_Pa)))
    delta = ((4098 * esat) / ((Ts_c + 237.3)**2))
    lmbd = (2501000 - (2361 * Ts_c))
    psychro = 1628.6 * (P_Pa / lmbd)
    Pair = (0.003486 * (P_Pa / Tvir))
    Cair = (0.622 * ((lmbd * psychro) / P_Pa))

    G = 0.1 * Rn

    # OK  ============
    Et_1 = (delta * (Rn - G)) / lmbd
    Et_2 = 3600 * Pair * Cair * VPD_Pa / Rb_H2O / lmbd

    Ei_3 = delta + psychro
    Ei_hr = (Et_1 + Et_2) / Ei_3 / 1000

    # PEt_3 = delta + psychro * (1 + (Rsto_PEt * Dratio) / Rb_H2O)
    # PEt_hr = (Et_1 + Et_2) / PEt_3 / 1000

    Et_3 = delta + psychro * (1 + (Rsto_c * DRATIO) / Rb_H2O)
    Et_hr_prev = Et_hr
    Et_hr = (Et_1 + Et_2) / Et_3 / 1000

    if (Es_blocked):
        Es_hr = 0
    else:
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lmbd
        Es_2 = ((3600 * Pair * Cair * VPD_Pa) - (delta * Rinc *
                                                 ((Rn - G) - (Es_Rn - Es_G)))) / (Rinc + Rb_H2O) / lmbd
        Es_3 = delta + (psychro * (1.0 + (Rsoil / (Rb_H2O + Rinc))))
        Es_hr = (Es_1 + Es_2) / Es_3 / 1000

    # Calculate AEt from Et and Es (after Shuttleworth and Wallace, 1985)
    SW_a = (delta + psychro) * Rb_H2O
    SW_s = (delta + psychro) * Rinc + (psychro * Rsoil)
    SW_c = psychro * (Rsto_c * DRATIO)  # Boundary layer resistance = 0
    C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
    C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
    if (Es_hr <= 0):
        AEt_hr = Et_hr
    else:
        AEt_hr = (C_canopy * Et_hr) + (C_soil * Es_hr)

    Ei_dd = Ei_dd + Ei_hr
    # PEt_dd = PEt_dd + PEt_hr
    # Et_dd = Et_dd + Et_hr
    Es_dd = 0.0 + Es_hr
    AEt_dd = AEt_dd + AEt_hr
    return (
        Ei_dd,
        # PEt_dd,
        # Et_dd,
        Es_dd,
        AEt_dd,
    )


(
    Ei_dd,
    # PEt_dd,
    # Et_dd,
    Es_dd,
    AEt_dd,
) = FORtran_pen_mon(
    Rn_MJ=0.00000000E+00,
    VPD=5.69900014E-02,
    LAI=3.00000000E+00,
    P=91.8599976,
    esat_kPa=0.773210466,
    eact_kPa=0.716220438,
    Ts_c=3.27999997,
    Rb_H2O=22.2595272,
    Rinc=216.42455968133493,
    Rsoil=200.000000,
    Rsto_c=100000.000,
    Es_blocked=False,

)
print(Ei_dd, AEt_dd, Es_dd)
isclose(Ei_dd, 3.769305726470758e-05, abs_tol=1e-8) \
    and isclose(AEt_dd, 2.465408666674146e-06, abs_tol=1e-8) \
    and isclose(Es_dd, 2.4442570459627592e-06, abs_tol=1e-8)


# %%

def penman_monteith_hourly(
    Rn_MJ: float,
    P_kPa: float,
    Ts_C: float,
    esat_kPa: float,
    eact_kPa: float,
    VPD_kPa: float,
    rm_nL: int,
    rm_Rb: float,
    rm_Rinc_l0: float,
    rm_Rsto_l0: float,
    rm_Rgs: float,
    LAI: float,
    Es_blocked: bool,
    pm_state_Ei_acc: float,
    pm_state_Et_acc: float,
    pm_state_Es_acc: float,
    pm_state_Eat_acc: float,
):

    # real :: Tvir, delta, lambda, psychro, Pair, Cair, G

    # real :: Et_1, Et_2, Ei_3, Et_3
    # real :: t, Es_Rn, Es_G, Es_1, Es_2, Es_3
    # real :: SW_a, SW_s, SW_c, C_canopy, C_soil

    # This model (probably) makes some one-layer assumptions, so don't allow
    # multi-layer resistance model.
    assert rm_nL == 1

    # Convert units and associate values
    Rn = Rn_MJ * 1000000
    P = P_kPa * 1000
    esat = esat_kPa * 1000
    eact = eact_kPa * 1000
    VPD = VPD_kPa * 1000

    Rb_H2O = rm_Rb
    Rinc = rm_Rinc_l0
    Rsto_H2O = rm_Rsto_l0
    Rsoil = rm_Rgs
    Ei_acc = pm_state_Ei_acc
    Et_acc = pm_state_Et_acc
    Es_acc = pm_state_Es_acc
    Eat_acc = pm_state_Eat_acc

    # OK  ============
    Tvir = (Ts_C + Ts_K) / (1 - (0.378 * (eact / P)))
    delta = ((4098 * esat) / ((Ts_C + 237.3)**2))
    lmd = (2501000 - (2361 * Ts_C))
    psychro = 1628.6 * (P / lmd)
    Pair = (0.003486 * (P / Tvir))
    Cair = (0.622 * ((lmd * psychro) / P))

    G = 0.1 * Rn

    # OK  ============
    Et_1 = (delta * (Rn - G)) / lmd
    Et_2 = 3600 * Pair * Cair * VPD / Rb_H2O / lmd

    Ei_3 = delta + psychro
    Ei = (Et_1 + Et_2) / Ei_3 / 1000

    # TODO: Make sure Rsto_H20 being passed in here
    Et_3 = delta + psychro * (1 + Rsto_H2O / Rb_H2O)
    Et = (Et_1 + Et_2) / Et_3 / 1000

    if Es_blocked:
        Es = 0
    else:
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lmd
        Es_2 = ((3600 * Pair * Cair * VPD) - (delta * Rinc *  # noqa: W504
                                              ((Rn - G) - (Es_Rn - Es_G)))) / (Rinc + Rb_H2O) / lmd
        Es_3 = delta + (psychro * (1.0 + (Rsoil / (Rb_H2O + Rinc))))
        Es = (Es_1 + Es_2) / Es_3 / 1000
    # Calculate Eat from Et and Es (after Shuttleworth and Wallace, 1985)
    SW_a = (delta + psychro) * Rb_H2O
    SW_s = (delta + psychro) * Rinc + (psychro * Rsoil)
    SW_c = psychro * Rsto_H2O  # Boundary layer resistance = 0
    C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
    C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
    if Es <= 0:
        Eat = Et
    else:
        Eat = (C_canopy * Et) + (C_soil * Es)

    # Accumulate values
    Ei_acc = pm_state_Ei_acc + Ei
    Et_acc = pm_state_Et_acc + Et
    Es_acc = pm_state_Es_acc + Es
    Eat_acc = pm_state_Eat_acc + Eat
    return(
        Ei_acc,
        Eat_acc,
        Es_acc,
    )


(Ei_acc, Eat_acc, Es_dd) = penman_monteith_hourly(
    Rn_MJ=0.0,  # OK
    P_kPa=91.86,  # CHECK!
    Ts_C=3.28,  # CHECK(Slightly different accuracy)
    esat_kPa=0.7732104619053026,  # OK
    eact_kPa=0.7162204619053026,  # OK
    VPD_kPa=0.05699,  # OK
    rm_nL=1.0,  # OK
    rm_Rb=22.25952506067983,  # OK
    rm_Rinc_l0=216.42455968133493,  # CHECK
    rm_Rsto_l0=99999.99999999999 * DRATIO,  # OK
    rm_Rgs=200,  # OK
    LAI=3.0,  # OK
    Es_blocked=False,  # OK
    pm_state_Ei_acc=0.0,
    pm_state_Et_acc=0.0,
    pm_state_Es_acc=0.0,
    pm_state_Eat_acc=0.0,
)
print(Ei_acc, Eat_acc, Es_dd)
# isclose(out.Ei_acc, 6.62323582E-05, abs_tol=1e-5) and isclose(out.Eat_acc, 1.29596865E-05, abs_tol=1e-5)
isclose(Ei_acc, 3.769305726470758e-05, abs_tol=1e-8) \
    and isclose(Eat_acc, 2.465408666674146e-06, abs_tol=1e-8) \
    and isclose(Es_dd, 2.4442570459627592e-06, abs_tol=1e-8)
# isclose(Ei_acc, 3.7693060526526357e-05, abs_tol=1e-5) and isclose(Eat_acc,
#                                                                   1.603862256714372e-08, abs_tol=1e-5)

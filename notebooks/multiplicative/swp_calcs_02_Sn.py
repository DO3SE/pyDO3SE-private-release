# %%
from math import isclose

from pyDO3SE.plugins.soil_moisture.helpers import (
    SWC_to_SWP,
    SWP_to_SWC,
    get_soil_config,
    soil_moisture_from_SWC,
    soil_moisture_from_SWP,
)
from pyDO3SE.plugins.soil_moisture.config import SOIL_LOAM, Soil_t
# ROW 26

# %%

# TODO: Check AEt
# TODO: Ckeck Ei
# TODO: Check AEt is the same as pm_state_Eat_acc
# TODO: Check Ei is the same as pm_state_Ei_accsb2457


def Calc_SWP_FORTRAN(
    precip_acc=0.000410999957239,
    AEt=0.00,  # CHECK THIS
    root=0.8,
    LAI=3,
    Ei=0.0,  # CHECK THIS
):
    if (precip_acc > 0):
        P_input = (precip_acc - (0.0001 * LAI)) + ((0.0001 * LAI) - min(Ei, 0.0001 * LAI))
    else:
        P_input = 0
    P_input = max(0.0, P_input)
    Sn_diff = (P_input - AEt) / root
    return Sn_diff


Calc_SWP_FORTRAN(
    precip_acc=0.000410999957239,
    AEt=0.00,  # CHECK THIS
    root=0.8,
    LAI=3,
    Ei=0.0,  # CHECK THIS
)

# %%


def calc_sn_diff(
    pm_state_precip_acc=0.000410999957239,  # OK
    run_off_fraction=0.0,  # CHECK THIS
    pm_state_run_off_acc=0.0,
    LAI=3,  # OK
    root_depth=0.8,  # OK
    pm_state_Eat_acc=0.0,
    pm_state_Ei_acc=0.0,  # OK
    ASW=0.13214041097790075,  # WRONG
):
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
    # precip_acc- Ei - AEt
    delta_SM = effective_irrig - intercepted_evaporated - evapotranspiration
    # Converted to volumetric change using root depth.
    Sn_diff = delta_SM / root_depth
    return Sn_diff


calc_sn_diff()


# ====== Day 0

Calc_SWP_FORTRAN(
    precip_acc=0.000410999957239,
    AEt=0.00,  # CHECK THIS
    root=0.8,
    LAI=3,
    Ei=0.0,  # CHECK THIS
) == calc_sn_diff(
    pm_state_precip_acc=0.000410999957239,
    run_off_fraction=0.0,  # CHECK THIS
    pm_state_run_off_acc=0.0,
    LAI=3,
    root_depth=0.8,
    pm_state_Eat_acc=0.0,
    pm_state_Ei_acc=0.0,
    ASW=0.13214041097790075,
)
# %%
# ======= Day 1
Calc_SWP_FORTRAN(
    precip_acc=1.000410999957239,
    AEt=1.00,  # CHECK THIS
    root=0.8,
    LAI=3,
    Ei=9.9,  # CHECK THIS
) == calc_sn_diff(
    pm_state_precip_acc=1.000410999957239,
    run_off_fraction=0.0,  # CHECK THIS
    pm_state_run_off_acc=0.0,
    LAI=3,
    root_depth=0.8,
    pm_state_Eat_acc=1.0,
    pm_state_Ei_acc=9.9,
    ASW=0.13214041097790075,
)

# %%
from math import isclose

from do3se_met.soil_moisture.helpers import (
    SWC_to_SWP,
    SWP_to_SWC,
    get_soil_config,
    soil_moisture_from_SWC,
    soil_moisture_from_SWP,
)
from pyDO3SE.plugins.soil_moisture.config import SOIL_LOAM, Soil_t


# %%
def Calc_SWP_FORTRAN(
    precip_acc,
    AEt,
    root=0.8,
    Fc_m=0.29,
    Sn=0.29,
    SWP_AE=-0.00188,
    SWC_sat=0.4,
    soil_b=6.58,
    PWP=-4,
):
    # FROM Init_SoilWater
    PWP_vol = 1.0 / (((PWP / SWP_AE)**(1.0 / soil_b)) / SWC_sat)
    print(PWP_vol)
    # PWP_vol = 0.12482448627762405
    print(PWP_vol)
    # Calc_SWP
    if (precip_acc > 0):
        P_input = (precip_acc - (0.0001 * LAI)) + ((0.0001 * LAI) - min(Ei, 0.0001 * LAI))
    else:
        P_input = 0
    # Can't lose water through Ei
    P_input = max(0.0, P_input)
    Sn_diff = (P_input - AEt) / root

    # Calculate new Sn, with field capacity as a maximum
    Sn = min(Fc_m, Sn + Sn_diff)
    Sn = max(PWP_vol, Sn)
    #Sn = 0.29
    # Calculate ASW and SWP for new water content
    ASW = (Sn - PWP_vol) * root

    # Calculate SMD for new water content
    return ASW


Calc_SWP_FORTRAN(
    0.0,
    0.0
)

# %%


def soil_moisture_from_SWP(soil_config: Soil_t, PWP: float, root_depth: float, SWP: float):
    """Fill soil moisture data from a soil water potential value.

    Args:
        real, intent(in) :: SWP     !< Soil water potential [MPa]
    """
    Output = namedtuple('Output', 'Sn SWP ASW SMD')
    output: Output = soil_moisture_from_SWC(
        soil_config,
        PWP,
        root_depth,
        Sn_in=SWP_to_SWC(soil_config.SWP_AE, soil_config.b, SWP))
    return output


soil_config = Soil_t(6.58, 0.29, -0.00188, 0.0002286)
soil_moisture = soil_moisture_from_SWP(
    soil_config,
    PWP=-4,
    root_depth=0.8,
    SWP=-0.0156003293092588,
)
assert isclose(soil_moisture.Sn, 0.289270281792, abs_tol=1e-3)
assert isclose(soil_moisture.SWP, -0.0158611088991, abs_tol=1e-3)
assert isclose(soil_moisture.ASW, 0.123334348202, abs_tol=1e-3)
assert isclose(soil_moisture.SMD, 0.000547282397747, abs_tol=1e-3)

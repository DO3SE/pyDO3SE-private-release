"""Helper functions for calculating atomspheric deposition."""

import numpy as np
from typing import List, Tuple
from collections import namedtuple
from scipy.linalg.lapack import sgesv

from do3se_met.resistance import (
    calc_Ra_simple,
    calc_Ra_with_heat_flux,
    calc_deposition_velocity,
)
from do3se_met.model_constants import izR


def O3_transfer_up(
    Ra: float,
    O3: float,
    Vd: float
) -> float:
    """Scale O3 concentration up from the resistance model's reference height.

    to 50m.

    Parameters
    ----------
    Ra: float
        resistance model Ra [unit]
    O3: float
        Ozone value [unit]
    Vd: float
        deposition velocity [unit]

    Returns
    -------
    O3_i: float
        O3 at 50m [unit]
    """
    O3_transfered = O3 / (1.0 - (Ra * Vd))
    return O3_transfered


def O3_transfer_down(Ra: float, O3_i: float, Vd: float) -> float:
    """Scale O3 concentration down from 50m to the resistance model's reference height.

    Ra: float
        resistance model Ra [unit]
    O3_i: float
        Ozone value at 50m [unit]
    Vd: float
        deposition velocity [unit]
    """
    O3_transfered = O3_i * (1.0 - (Ra * Vd))
    return O3_transfered


OzoneConcentrationOutput = namedtuple('OzoneConcentrationOutput', 'O3_i micro_O3 Vd_i Vd')


def calc_canopy_ozone_concentration(
    O3_ppb_zR: float,
    Ra_ref_canopy: float,
    Ra_ref_measured: float,
    Ra_tar_canopy: float,
    Ra_tar_canopy_top: float,
    Rsur_ref: float,
    Rsur_top_layer: float,
    Rb_ref: float,
    Rb_top_layer: float,
) -> OzoneConcentrationOutput:
    r"""Calculate the ozone concentration at top of canopy.

    This method is taken from the DO3SE-UI model and should be
    replaced with calc_canopy_ozone_concentration.

    This procedure results in the calculation of deposition velocity (Vd) and
    the ozone concentration at the canopy in both parts-per-billion and
    nmol/m^3

    We translate the ozone from the measured height up to a decoupled height
    then back down to the target canopy. As the measured data may have had a canopy
    of a different height we need to calculate the deposition velocity for a
    reference canopy.


    --------------------------------------------------- Decoupled height
              ^ O3_ppb_i                |
              |                         |
              | Ra_O3zR_i * vd_i        | Ra * Vd
              |                         |/_\ Vd =1/(Ra + Rb + Rsur)
    O3_ppb_zR |                         |
    ------------------------------------|------------------------- Measured height
                                        |
                                        |
    ------------------------------------|--------------------------Measured Canopy height
                                        |
                                        |
    ------------------------------------|----------------Canopy height
    ------------------------------------|---------------
                                        |               ^
                                        |               zo
                                        |               v
    ------------------------------------v---------------
                                                        ^
                                                        d
                                                        v
    ---------------------------------------------------- Soil



    uh   |=>|ustar_ref | => |Ra_O3zR_i| ===========> |
    invL |             |                             |
                    | => |Ra_ref_i | => |vd_i |=> |O3_ppb_i|
                    | => |Rb_ref   |
                                    |
                        ||Rsur     |

    Ra_O3zR_i = calc_ra_with_heat_flux(ustar_ref, O3zR, izR, invL)
    vd_i = 1.0 / (Ra_ref_i + Rb_ref + Rsur)


    Parameters
    ----------
    O3_ppb_zR: float
        input O3 data[ppb]
    Ra_ref_canopy: float
        Ra for reference canopy to izr
    Ra_ref_measured: float
        Ra for measured height to izr
    Ra_tar_canopy: float
        Ra for target canopy to izr
    Ra_tar_canopy_top: float
        Ra for top of target canopy to izr
    Rsur_ref: float
        Rsur for ref canopy
    Rsur_top_layer: float
        Rsur for top layer
    Rb_top_layer: float
        Rb for top layer
    Rb_ref: float
        Rb reference canopy
    Ra_top_layer: float
        Ra for top layer

    Returns
    -------
    O3_i: float
        Ozone concentration at izR[ppb]
    micro_O3
        Ozone concentration at canopy top[ppb]

    """
    # z1 = canopy_ref_z0 and z2 = izr
    Vd_i = 1.0 / (Ra_ref_canopy + Rb_ref + Rsur_ref)

    # z1 = measure height - canopy_ref_d and z2 = izr
    O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_ref_measured * Vd_i))

    # z1 = canopy_zo and z2 = izr
    Vd = 1.0 / (Ra_tar_canopy + Rb_top_layer + Rsur_top_layer)

    # z1 = canopy_height - canopy_d and z2 = izr
    O3_ppb = O3_ppb_i * (1.0 - (Ra_tar_canopy_top * Vd))

    return OzoneConcentrationOutput(
        O3_i=O3_ppb_i,
        micro_O3=O3_ppb,
        Vd_i=Vd_i,
        Vd=Vd,
    )


def calc_canopy_ozone_concentration_alt(
    O3: float,
    ustar_ref: float,
    canopy_height: float,
    u_i: float,
    z_O3: float,
    Rtotal_top_layer: float,
    O3_d: float,
    O3_z0: float,
    L: float,
    h_O3_in: float = None,
    izr=izR,
    # TODO: Check ra method strings
    ra_method: str = "simple"
) -> Tuple[float, float]:
    """Calculate the ozone concentration at the canopy.

    To do this we must assume the measured ozone is potentially influenced by the vegitation
    at the measured location. To resolve this issue we must estimate the ustar at a reference height
    where the atmosphere is assumed to be decoupled from the underlying vegitation.
    Once we have this we can then estimate the concentration at our canopy height by
    transfering down from the reference height.

    # TODO: Check this against the calculations in Calc_O3_concentration in DO3SE UI

    Parameters
    ----------
    O3: float
        input O3 data
    ustar_ref : float
        Friction velocity for ozone at reference height [m/s]
    canopy_height: float
        The overall height of the canopy [m]
    u_i: float
        The windspeed at 50m [m^s]
    h_O3_in: float
        target height [m]
    z_O3: float
        measured height [m]
    Rtotal_top_layer: float,
        Total resistance at top of canopy
    O3_d: float
        Ozone measured canopy displacement height [m]
    O3_z0: float
        Ozone measured canopy roughness length [m]

    Returns
    -------
    O3_i: float
        [description] [unit]
    micro_O3
        [description] [unit]
    """
    assert O3 is not None
    assert canopy_height is not None
    assert u_i is not None
    assert z_O3 is not None
    assert Rtotal_top_layer is not None
    Output = namedtuple('Output', 'O3_i micro_O3')
    # h_o3 is the assumed canopy height
    h_O3 = h_O3_in if h_O3_in is not None else canopy_height
    assert h_O3 is not None

    # O3_d = h_O3 * canopy_d
    # O3_z0 = h_O3 * canopy_z0

    # 1. transfer the ozone measurement from measured height to 50m
    # Calculate the resistances at the O3 ref height

    # ustar_ref = ustar_from_velocity(u_i, 50.0 - O3_d, O3_z0)

    # TODO: Should input ra method option here
    # 2. Get ra between canopy and 50m
    rmodel_Ra_c_i = calc_Ra_simple(ustar_ref, O3_z0, izr, O3_d) if ra_method == "simple" \
        else calc_Ra_with_heat_flux(ustar_ref, O3_z0 + O3_d, izr, L)

    # 3. Get ra between measured height and 50m
    rmodel_Ra_i = calc_Ra_simple(ustar_ref, z_O3, izr, O3_d) if ra_method == "simple" \
        else calc_Ra_with_heat_flux(ustar_ref, z_O3, izr, L)

    # rmodel_Rb_i = calc_Rb(ustar_ref, DIFF_O3)

    # TODO: Check inputs and calcs for calc_deposition_velocity
    # 4. Calculate the average deposition velocity between 50m and canopy
    Vd_c = calc_deposition_velocity(rmodel_Ra_c_i, Rtotal_top_layer)
    Vd = calc_deposition_velocity(rmodel_Ra_i, Rtotal_top_layer)

    # 5. use Ra between measured and 50m to get O3 at 50m
    # TODO: Should we use Vd between measured and 50m?
    O3_i = O3_transfer_up(rmodel_Ra_i, O3, Vd)
    # 6. use Ra between 50m and canopy to get O3 at canopy height
    # R1 = Rsto_top_layer + Rext_top_layer
    # R_total = Ra + Rb + R1
    # O3 = O3_i * (R1 / R_total)
    O3_at_target_height = O3_transfer_down(rmodel_Ra_c_i, O3_i, Vd_c)

    return Output(
        O3_i=O3_i,
        micro_O3=O3_at_target_height,
    )


def calc_multi_layer_O3_ozone_concentration(
    nL: int,
    O3_in: float,
    rm_Ra: float,
    rm_Rinc: List[float],
    rm_Rsur: List[float],
    rm_Rgs: float,
) -> List[float]:
    """Calculate O3 concentration for all layers.

    Requires that the value for the top layer (umet(1)%O3) is already known.

    # TODO: Is this absorbed O3 at each layer?

    Assumes layer 0 is top layer.

    This uses SGESV which is an old Fortran function. The documentation is vague on this.

    Parameters
    ----------
    nL: float
        number of model layers
    O3_in: float
        O3 for top layer micromet[ppb]
    rm_Ra: float
        rmodel_O3 Ra
    rm_Rinc: List[float]
        rmodel_O3 Rinc per layer
    rm_Rsur: List[float]
        rmodel_O3 Rsur per layer
    rm_Rgs: float
        rmodel_O3 Rgs

    Output
        O3 per layer[ppb]

    """
    # real, dimension(nL+1) :: bigR, smallR, C
    # real, dimension(nL+1,nL+1) :: X
    # bigR = np.full((nL + 1), None, dtype=float)
    # smallR = np.full((nL + 1), None, dtype=float)
    C = np.full((nL + 1), 0, dtype=float)
    X = np.full((nL + 1, nL + 1), 0, dtype=float)

    bigR = np.array([rm_Ra] + rm_Rinc)
    assert bigR.shape == (nL + 1,)

    # TODO: per-layer Rb
    smallR = np.array(rm_Rsur + [rm_Rgs])
    assert smallR.shape == (nL + 1,)
    # Iterate over columns
    for j in range(0, nL + 1):
        X[0:j, j] = bigR[0:j]
        X[j, j] = X[j, j] + smallR[j]
        if j < nL:  # TODO: Check this
            X[j + 1, j] = -smallR[j]

    # NOTE: Below was first attempt at copying fortran.
    # Can probably just remove this once happy with above.
    # for i in range(0, nL + 1):
    #     # Fill in "R" values in upper triangle
    #     # TODO: Check this works in python
    #     for j in range(0, i + 1):
    #         X[j, i] = bigR[j]
    #     # Fill in "r" values on diagonal
    #     X[i, i] = X[i, i] + smallR[i]
    #     # Fill in "r" values below diagonal
    #     if i <= nL - 1:
    #         X[i + 1, i] = -smallR[i]

    C[0] = O3_in

    # NOTE: Very little documentation on SGESV
    # N = nL + 1
    # NRHS = 1,
    # A = X
    # LDA = nL + 1
    # IPIV = ipiv_
    # B = C
    # LDB = nL + 1

    # call SGESV(nL + 1, 1, X, nL + 1, ipiv_, C, nL + 1, info)
    out = sgesv(X, C)
    lu, IPIV, C_out, info = sgesv(X, C)

    if info != 0:
        raise Exception('SGESV Failed')
    C_final = smallR * C_out
    O3_out = C_final[0: nL]
    return O3_out


def calc_multi_layer_O3_ozone_concentration_simple(
    nL: int,
    O3_in: float,
    rm_Ra: float,
    rm_Rb: float,
    rm_Rinc: List[float],
    rm_Rsur: List[float],
    rm_Rgs: float,
) -> List[float]:
    """Calculate O3 concentration for all layers.

    Requires that the value for the top layer (umet(1)%O3) is already known.

    Assumes layer 0 is top layer.

    Parameters
    ----------
    nL: float
        number of model layers
    O3_in: float
        O3 for top layer micromet
    rm_Ra: float
        rmodel_O3 Ra
    rm_Ra: float
        rmodel_O3 Rb
    rm_Rinc: List[float]
        rmodel_O3 Rinc per layer
    rm_Rsur: List[float]
        rmodel_O3 Rsur per layer
    rm_Rgs: float
        rmodel_O3 Rgs

    Output
        O3 per layer
    """
    total_resistance_per_layer = [sum([rm_Ra, rm_Rb] + rm_Rinc[0:iL]) for iL in range(nL)]

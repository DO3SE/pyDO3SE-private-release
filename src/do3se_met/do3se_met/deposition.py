"""Helper functions for calculating atomspheric deposition."""

import numpy as np
from typing import List
from collections import namedtuple
from scipy.linalg.lapack import sgesv  # type: ignore
from do3se_met.resistance import calc_deposition_velocity_multilayer

def O3_transfer_up(Ra: float, O3: float, Vd: float) -> float:
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


OzoneConcentrationOutput = namedtuple(
    "OzoneConcentrationOutput", "O3_i micro_O3 Vd_i Vd"
)


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
    vd_i: float
        deposition velocity at izR [m/s]
    vd: float
        deposition velocity at canopy top [m/s]

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


def calc_canopy_ozone_concentration_multilayer(
    O3_ppb_zR: float,
    Ra_ref_canopy: float,
    Ra_ref_measured: float,
    Ra_tar_canopy: float,
    Ra_tar_canopy_top: float,
    Rsur_ref: float,
    Rb_ref: float,
    Rtotal: float,
):
    r"""Calculate the ozone concentration at top of canopy.

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
    Rb_ref: float
        Rb reference canopy
    Rtotal: float
        Total resistance for canopy

    Returns
    -------
    O3_i: float
        Ozone concentration at izR[ppb]
    micro_O3
        Ozone concentration at canopy top[ppb]
    vd_i: float
        deposition velocity at izR [m/s]
    vd: float
        deposition velocity at canopy top [m/s]

    """
    # TODO: Compare this with single layer method
    # z1 = canopy_ref_z0 and z2 = izr
    Vd_i = 1.0 / (Ra_ref_canopy + Rb_ref + Rsur_ref)

    # z1 = measure height - canopy_ref_d and z2 = izr
    O3_ppb_i = O3_ppb_zR / (1.0 - (Ra_ref_measured * Vd_i))

    # z1 = canopy_zo and z2 = izr
    Vd = calc_deposition_velocity_multilayer(
        Ra_tar_canopy, Rtotal,
    )

    # z1 = canopy_height - canopy_d and z2 = izr
    O3_ppb = O3_ppb_i * (1.0 - (Ra_tar_canopy_top * Vd))

    return OzoneConcentrationOutput(
        O3_i=O3_ppb_i,
        micro_O3=O3_ppb,
        Vd_i=Vd_i,
        Vd=Vd,
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
    if sum(bigR) == 0:
        return [O3_in] * (nL + 1)
    assert bigR.shape == (nL + 1,)

    # TODO: per-layer Rb
    smallR = np.array(rm_Rsur + [rm_Rgs])
    assert smallR.shape == (nL + 1,)
    # Iterate over columns
    for j in range(0, nL + 1):
        X[0 : j + 1, j] = bigR[0 : j + 1]
        if j < nL + 1:
            X[j, j] = X[j, j] + smallR[j]
        if j < nL:
            X[j + 1, j] = -smallR[j]
    C[0] = O3_in

    # NOTE: Very little documentation on SGESV
    # N = nL + 1
    # NRHS = 1,
    # A = X
    # LDA = nL + 1
    # IPIV = ipiv_
    # B = C
    # LDB = nL + 1

    # call SGESV(nL + 1, 1, X, nL + 1, ipiv_, C, nL + 1, info) # Eq fortran code
    lu, IPIV, C_out, info = sgesv(X, C)

    if info != 0:
        raise Exception("SGESV Failed")
    C_final = smallR * C_out
    return C_final


def calc_ozone_at_custom_height_linear(
    ozone_at_layers: List[float] | np.ndarray,
    layer_heights: List[float] | np.ndarray,
    custom_height: float,
    izr_height: float | None = None,
    ozone_at_izr_height: float | None = None,
) -> float:
    """Calculate the ozone concentration at a custom height.

    This method uses linear interpolation to calculate the ozone at at the custom height.

    Parameters
    ----------
    ozone_at_layers: List[float]
        Ozone concentration at each layer
    layer_heights: List[float]
        Height of each layer
    custom_height: float
        Custom height to calculate ozone concentration
    izr_height: float
        Decoupled height
    ozone_at_izr_height: float
        Ozone concentration at decoupled height

    """
    layer_heights = (
        layer_heights if izr_height is None else layer_heights + [izr_height]
    )
    ozone_at_layers = (
        ozone_at_layers
        if ozone_at_izr_height is None
        else ozone_at_layers + [ozone_at_izr_height]
    )
    return np.interp(custom_height, layer_heights, ozone_at_layers)


def calc_ozone_at_custom_height_polynomial(
    ozone_at_layers: List[float] | np.ndarray,
    layer_heights: List[float] | np.ndarray,
    custom_height: float,
    izr_height: float | None = None,
    ozone_at_izr_height: float | None = None,
) -> float:
    """Calculate the ozone concentration at a custom height.

    This method fits a polynomial to the ozone concentration at each layer and
    then evaluates the polynomial at the custom height.

    Parameters
    ----------
    ozone_at_layers: List[float]
        Ozone concentration at each layer
    layer_heights: List[float]
        Height of each layer
    custom_height: float
        Custom height to calculate ozone concentration
    izr_height: float
        Decoupled height
    ozone_at_izr_height: float
        Ozone concentration at decoupled height

    """
    layer_heights = (
        layer_heights if izr_height is None else layer_heights + [izr_height]
    )
    ozone_at_layers = (
        ozone_at_layers
        if ozone_at_izr_height is None
        else ozone_at_layers + [ozone_at_izr_height]
    )
    f = np.polynomial.polynomial.Polynomial.fit(layer_heights, ozone_at_layers, 5)
    return float(f(custom_height))

"""Photosynthesis gsto helpers."""


from functools import partial
from do3se_met.f_functions import inverse_f_VPD_linear, inverse_f_VPD_log
from pyDO3SE.Config.ConfigEnums import FVPDMethods


def calc_D_0_f_VPD_Method(
    f_VPD_method: FVPDMethods,
    fmin: float,
    VPD_max: float = None,
    VPD_min: float = None,
) -> float:
    """Get D_0 when D_0_method is "f_VPD".

    Note use of partial to front load functions
    """
    vpd_linear_fn = partial(inverse_f_VPD_linear, 0.5, VPD_max, VPD_min, fmin)
    vpd_log_fn = partial(inverse_f_VPD_log, .5, fmin)

    D_0 = vpd_linear_fn() if f_VPD_method == FVPDMethods.LINEAR else None
    D_0 = vpd_log_fn() if f_VPD_method == FVPDMethods.LOG else D_0
    return D_0


def calc_D_0(
    D_0_method: str,
    f_VPD_method: FVPDMethods = None,
    constant_D_0: float = None,
    VPD_max: float = None,
    VPD_min: float = None,
    fmin: float = None,
) -> float:
    """Get D_0 based on the D_0_methods specified in the config.

    D_0 "The VPD at which g_sto is reduced by a factor of 2" [kPa] (Leuning et al. 1998)

    Parameters
    ----------
    D_0_method : str
        One of ['constant', 'f_VPD']
    f_VPD_method : FVPDMethods, optional
        One of FVPDMethods, by default None
    constant_D_0 : float, optional
        Constant D_0 value, by default None
    VPD_max : float, optional
        [description], by default None
    VPD_min : float, optional
        [description], by default None
    fmin : float, optional
        [description], by default None

    Returns
    -------
    D_0: float
        [description]

    Raises
    ------
    ValueError
        [description]
    Exception
        [description]
    """
    if D_0_method == "constant":
        assert constant_D_0 is not None
        D_0 = constant_D_0
    elif D_0_method == "f_VPD":
        D_0 = calc_D_0_f_VPD_Method(f_VPD_method, fmin, VPD_max, VPD_min)
    else:
        raise ValueError(f'Invalid D_0_method {D_0_method}')
    if D_0 is None:
        raise Exception('Invalid D_0 method or missing arguments')
    return D_0

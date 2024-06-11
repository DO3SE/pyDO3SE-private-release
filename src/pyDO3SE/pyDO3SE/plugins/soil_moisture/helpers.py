"""Helper functions related to soil moisture."""
from collections import namedtuple
from typing import NamedTuple

from pyDO3SE.plugins.soil_moisture.config import (
    SOIL_CLAY_LOAM,
    SOIL_LOAM,
    SOIL_SANDY_LOAM,
    SOIL_SILT_LOAM,
    Soil_t,
)

SWC_sat = 0.4  # Saturated soil water content for soil water release curve


def init_soilWater_state(
    SWP_min: float,
    SWP_AE: float,
    soil_b: float,
    Fc_m: float,
    root_depth: float,
    SWC_sat: float = SWC_sat,
) -> NamedTuple:
    """Initialize soil moisture variables.

    TODO: Copied directly from DO3SE UI. Check all calcs are needed.
    Unused equations are commented out
    """
    # use Constants, only: SWC_sat
    # use Parameters, only: Fc_m, soil_b, SWP_AE, D_meas
    # use Parameters, only: SWP_min, SWP_max, fmin, root
    # use Variables, only: Sn_star, Sn, per_vol, ASW, SWP, &
    #                         fSWP, SMD, AEt, Et, Es, PEt, Ei, fLWP, &
    #                         Sn_meas, SWP_meas, SMD_meas
    Output = namedtuple("Output", [
        "Sn", "ASW", "SWP", "SMD", "ASW_FC"
    ])
    # Ei_dd = 0
    # PEt_dd = 0
    # Et_dd = 0
    # Es_dd = 0
    # AEt_dd = 0
    # Ei = 0
    # PEt = 0
    # Et = 0
    # Es = 0
    # AEt = 0

    # PWP can't be any higher than SWP_min
    PWP = min(-4.0, SWP_min)
    # Convert PWP to volumetric (MPa -> m^3/m^3)
    PWP_vol = 1.0 / (((PWP / SWP_AE)**(1.0 / soil_b)) / SWC_sat)

    # Volumetric water content, initially at field capacity
    Sn_star = Fc_m
    Sn = Sn_star

    # As a percentage
    # per_vol = Sn * 100

    # ASW and SWP for initial volumetric water content
    ASW = (Sn - PWP_vol) * root_depth
    ASW_FC = (Fc_m - PWP_vol) * root_depth
    SWP = SWP_AE * ((SWC_sat / Sn)**soil_b)

    # Calculate fSWP and SMD for initial water content
    # fSWP = fSWP_exp_curve(SWP, fmin)
    SMD = (Fc_m - Sn) * root_depth

    # Initial fLWP = 1
    # fLWP = 1

    # Initialised SWP_meas
    # r_meas = (1 - (0.97**(D_meas * 100)))
    # Sn_meas = Sn_star
    # SWP_meas = SWP_AE * ((SWC_sat / Sn_meas)**soil_b)

    return Output(
        Sn=Sn,
        ASW=ASW,
        SWP=SWP,
        SMD=SMD,
        ASW_FC=ASW_FC,
    )


def get_soil_config(name=None, soil_in=None) -> Soil_t:
    """Check that a Soil_t is valid.  The optional *name* argument can be used to
    select one of the preset values, passing name="" or name="custom" has the
    same effect as omitting the argument.

    Returns
    -------
        type(Soil_t), intent(inout): : soil
        character(len=*), intent(in ), optional : : name
    """
    if name == "" or name == "custom" or name is None:
        # Parameters should have been set manually
        soil = soil_in
    elif name == "sandy loam":
        soil = SOIL_SANDY_LOAM
    elif name == "silt loam":
        soil = SOIL_SILT_LOAM
    elif name == "loam":
        soil = SOIL_LOAM
    elif name == "clay loam":
        soil = SOIL_CLAY_LOAM
    else:
        raise ValueError(f"unknown soil texture: {name}")

    # Check that parameters were set by the name or manually
    assert soil.b is not None
    assert soil.FC is not None
    assert soil.SWP_AE is not None
    assert soil.Ksat is not None

    return soil


def calc_ASW(soil_config: Soil_t, PWP: float, root_depth: float) -> float:
    """Calculate available soil water at field capacity"""
    tmp_data = soil_moisture_from_SWC(
        soil_config,
        PWP,
        root_depth,
        Sn_in=soil_config.FC
    )
    return tmp_data.ASW


# TODO: Below has been split into default processes. Can be removed
# def check_SMDConfig(config: Soil_Moisture_Config):
#     """  # Check that a SMDConfig_t is valid.

#     - Checks that the Soil_t is valid
#     - Applies named soil texture if necessary
#     - Calculates various soil-dependent constants

#     Inputs
#     ======
#     type(SMDConfig_t), intent(inout):: config

#     """
#     # TODO: Test Before implementing
#     raise NotImplementedError()

#     # Check/set soil texture parameters
#     check_Soil(config.soil, config.soil_texture)

#     # Calculate available soil water at field capacity
#     tmp_data: SMDData_t = soil_moisture_from_SWC(config, config.soil.FC)
#     # TODO: Should not be setting config here
#     config.ASW_FC = tmp_data.ASW

#     # Set up soil water data source parameters
#     if config.source == "P-M":
#         if config.initial_SWC is None:
#             config.initial_SWC = config.soil.FC


def SWC_to_SWP(SWP_AE, b, SWC) -> float:
    """Use soil water release curve to convert from soil water content (m3 m-3) to
    soil water potential [MPa].

    Returns
    -------
        real, intent( in ) : : SWP_AE    !< Water potential at air entry [MPa]
        real, intent( in ) : : b         !< Texture dependent soil conductivity parameter
        real, intent( in ) : : SWC       !< Soil water content (m3 m-3)
    """

    SWC_sat = 0.4  # Saturated soil water content for soil water release curve

    SWP = SWP_AE * ((SWC_sat / SWC)**b)
    return SWP


def SWP_to_SWC(SWP_AE, b, SWP):
    """Use soil water release curve to convert from soil water potential [MPa] to
    soil water content (m3 m-3).

    Returns
    -------
        real, intent(in) :: SWP_AE    !< Water potential at air entry [MPa]
        real, intent(in) :: b         !< Texture dependent soil conductivity parameter
        real, intent(in) :: SWP       !< Soil water potential [MPa]
    """
    SWC = 1.0 / (((SWP / SWP_AE)**(1.0 / b)) / SWC_sat)
    return SWC


# TODO: Specify args instead of using soil_config
def soil_moisture_from_SWC(soil_config: Soil_t, PWP: float, root_depth: float, Sn_in: float):
    """# Fill soil moisture data from a soil water content value.

    Returns
    -------
        soil_config !< Soil texture parameters
        real, intent(in) :: PWP  !< # "Permanent wilting point", minimum level for SWP [MPa]
        real, intent(in) :: root_depth !< Root depth[m]
        real, intent(in) :: Sn   !< Soil water content (m3 m-3)
    """
    Output = namedtuple('Output', 'Sn SWP ASW SMD')
    # Convert PWP to volumetric content to use as a minimum soil water content
    PWP_vol = SWP_to_SWC(soil_config.SWP_AE, soil_config.b, PWP)
    # Constrain soil water content to be between field capacity and PWP
    Sn = max(PWP_vol, min(soil_config.FC, Sn_in))

    # Calculate soil water potential (SWP)
    SWP = SWC_to_SWP(soil_config.SWP_AE, soil_config.b, Sn)

    # Calculate available soil water (ASW)
    ASW = (Sn - PWP_vol) * root_depth

    # Calculate soil moisture deficit (SMD)
    SMD = (soil_config.FC - Sn) * root_depth

    return Output(
        Sn=Sn,
        SWP=SWP,
        ASW=ASW,
        SMD=SMD,
    )


def soil_moisture_from_SWP(soil_config: Soil_t, PWP: float, root_depth: float, SWP: float):
    """Fill soil moisture data from a soil water potential value.

    Returns
    -------
        real, intent(in) :: SWP     !< Soil water potential [MPa]
    """
    Output = namedtuple('Output', 'Sn SWP ASW SMD')
    output: Output = soil_moisture_from_SWC(
        soil_config,
        PWP,
        root_depth,
        Sn_in=SWP_to_SWC(soil_config.SWP_AE, soil_config.b, SWP))
    return output

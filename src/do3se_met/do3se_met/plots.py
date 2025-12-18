"""Plotting utils"""

from typing import Optional
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np

from do3se_met.f_functions import f_PAW, f_SWP_linear
from do3se_met.soil_moisture.enums import FSW_Methods


def plot_f_SW_curve(
    f_SW_method: FSW_Methods,
    ASW_FC: float = None,
    fmin: float = None,
    ASW_min: float = None,
    ASW_max: float = None,
    SWP_min: float = None,
    SWP_max: float = None,
    fig: Optional[Figure] = None,
    axes: Optional[Axes] = None,
):
    """Plot f_SW curve"""
    fig = fig or plt.figure()
    axes = axes or fig.add_subplot(1, 1, 1)

    if f_SW_method == FSW_Methods.FPAW:
        ASW_values = np.linspace(0, 1, 100)
        f_SW_curve = [
            f_PAW(
                ASW_FC=ASW_FC,
                fmin=fmin,
                ASW=ASW,
                ASW_min=ASW_min,
                ASW_max=ASW_max,
            )
            for ASW in ASW_values
        ]
        axes.plot(ASW_values, f_SW_curve)
        axes.set_xlabel("Available soil water (m3 m-3)")
        axes.set_ylabel("f_SW")
        axes.set_title("fPAW method")
    elif f_SW_method == FSW_Methods.FSWP_LINEAR:
        SWP_values = np.linspace(0, 1, 100)
        f_SW_curve = [
            f_SWP_linear(
                SWP_min=SWP_min,
                SWP_max=SWP_max,
                fmin=fmin,
                SWP=SWP,
            )
            for SWP in SWP_values
        ]
        axes.plot(SWP_values, f_SW_curve)
        axes.set_xlabel("Available soil water (m3 m-3)")
        axes.set_ylabel("f_SW")
        axes.set_title("Linear SWP method")
    else:
        raise NotImplementedError(f"f_SW method {f_SW_method} not implemented")

    return fig, axes

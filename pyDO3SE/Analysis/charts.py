"""Create charts of pypDO3SE output."""

# TODO: Should import output data

from math import ceil
from typing import Any, List, Optional
from pathlib import Path
import warnings
from matplotlib import pyplot as plt, colors
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from datetime import datetime
import pandas as pd
from pyDO3SE.Analysis.util import MONTHS

from pyDO3SE.util.Objects import Field


# plt.style.use('pyDO3SE/Analysis/style.mplstyle')


def calc_moving_average(data: List[float] | list[int], size: int) -> List[float] | List[int]:
    if len(data) == 0:
        raise ValueError("Data is empty")
    data_window = [data[0] for i in range(size)]

    moving_average = []
    for i in data:
        data_window = data_window[1:] + [i]
        moving_average.append(sum(data_window) / size)
    return moving_average


def log(log_level):
    def _inner(*args, **kwargs):
        if log_level:
            print(*args, **kwargs)

    return _inner


def annual_graph(
    data: List[float] | list[int],
    field: Field,
    start_day: int = 0,
    end_day: int = 365,
    output_dir: Optional[str | Path] = None,
    chart_id: Optional[str | int] = None,
    label_x_days: int = 1,
    average_step: int = 1,
    output_format: str = "png",
    log_level: int = 1,
    fig: Optional[Figure] = None,
    axes: Optional[Axes] = None,
    plot_kwargs: dict = {},
    *args,
    **kwargs,
) -> tuple[Figure, Axes]:
    """Create a graph of a single series between two dates.

    Parameters
    ----------
    data : List[float|int]
        Array of the data to plot
    field : Field
        Field object describing the data
    start_day : int, optional
        Start day dd, by default 0
    end_day : int, optional
        End day dd, by default 365
    output_dir : str, optional
        Directory to save plots, by default None
    chart_id : str, optional
        Unique identifier to save the file as, by default None
    label_x_days : int, optional
        Number of days between x-axis labels, by default 1
    moving_average: int, optional
        Size of the moving average window, by default 1
    fig : Figure, optional
        Matplotlib figure to use, by default None
    axes : Axes, optional
        Matplotlib axes to use, by default None
    plot_kwargs: dict, optional
        kwargs to pass to axes.plot()
    Returns
    -------
    Figure
        Matplotlib figure

    Raises
    ------
    ValueError
        Raise a value error if input is incorrect
    """
    logger = log(log_level)
    fig = fig or plt.figure(*args, **kwargs)
    title = ""

    # 9 days  with 3 day spacing
    # Get x ticks
    t_count = ceil((end_day - start_day) / label_x_days)
    t_spacing = (end_day - start_day) / t_count

    labels = [f"{int(i * t_spacing) + start_day}" for i in range(t_count + 1)]
    labelx = [i * 24 * t_spacing for i in range(t_count + 1)]

    if isinstance(field, list):
        # assume fields map to outputdata lists
        raise ValueError("Use multi_series_annual_graph for multi series graph")
    # Assume only a single data series
    title = field.short
    fig.suptitle(field.short)
    axes = axes or fig.add_subplot(1, 1, 1)
    plot_data = data if average_step == 1 else calc_moving_average(data, average_step)
    axes.plot(plot_data, **plot_kwargs)
    axes.set_xticks(labelx)
    axes.set_xticklabels(labels)
    axes.set_xlabel("Days")
    axes.set_ylabel(f"{field.short} ({field.unit})")
    axes.yaxis.grid(True)

    if output_dir:
        filename = f"{title}-{chart_id or str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S'))}"
        fig.savefig(
            f"{output_dir}/{filename}.{output_format}", format=output_format, transparent=True
        )
        logger(f"Output saved to {filename}")

    return fig, axes


def monthly_diurnal_graph(
    data: pd.DataFrame,
    observed_diurnal_data: Optional[pd.DataFrame],
    field: Field,
    month: int,
    ylim: tuple[float, float],
    output_dir: Optional[Path | str] = None,
    chart_id: Optional[str] = None,
    output_format: str = "png",
    log_level: int = 1,
    close: bool = True,
    fig: Optional[Figure] = None,
    axes: Optional[Axes] = None,
    fill_color: str = "#1f77b4",
    observed_fill_color: str = "#ff7f0e",
    fill_opacity: float = 0.1,
    *args,
    **kwargs,
) -> tuple[Figure, Axes]:
    """Create a graph for monthly average diurnal plot

    Parameters
    ----------
    data : DataFrame
        pandas dataframe with 'hr' column and field.id column
    observed_diurnal_data : DataFrame
        pandas dataframe with 'hr' column and 'observed' column
    field : Union[Field, List[Field]]
        Field object describing the data
    month : int
        Month number 0-11
    ylim : tuple[float,float]
        y-axis limits [min, max]
    output_dir : str, optional
        Directory to save plots, by default None
    chart_id : str, optional
        Unique identifier to save the file as, by default None
    output_format : str, optional
        File format to save the plot as, by default 'png'
    log_level : int, optional
        Log level, by default 1
    close : bool, optional
        Whether to close the figure after saving, by default True
    fig : Figure, optional
        Matplotlib figure to use, by default None
    axes : Axes, optional
        Matplotlib axes to use, by default None
    fill_color : str, optional
        Color to use for the fill, by default '#1f77b4'
    observed_fill_color : str, optional
        Color to use for the observed fill, by default '#ff7f0e'
    fill_opacity : float, optional
        Opacity to use for the fill, by default 0.1

    Returns
    -------
    tuple[Figure, Axes]
        Matplotlib figure and axes

    Raises
    ------
    ValueError
        Raise a value error if input is incorrect
    """
    logger = log(log_level)
    fig = fig or plt.figure(*args, **kwargs)
    ax = axes or fig.add_subplot(1, 1, 1)
    if isinstance(field, list):
        # assume fields map to outputdata lists
        raise ValueError("Multi series not supported")
    # Assume only a single data series
    title = field.short
    # fig.suptitle(f'{field.short} - {MONTHS[month]}')
    ax.set_title(f"{field.short} - {MONTHS[month]}")
    _df = data.groupby("hr")[field.id].agg(["min", "max", "mean"])

    _df["mean"].plot(grid=True, color=fill_color, ax=ax)

    if observed_diurnal_data is not None:
        _odf = observed_diurnal_data.groupby("hr")["observed"].agg(["min", "max", "mean"])
        _odf["mean"].plot(ax=ax, color=observed_fill_color, grid=True)
        ax.fill_between(
            _odf.index,
            _odf["min"],
            _odf["max"],
            color=colors.to_rgba(observed_fill_color, fill_opacity),
        )
        ax.legend([field.id, "observed"], loc="upper right")

    ax.set_xlim((0.0, 23.0))
    ax.set_ylim(ylim)
    ax.fill_between(
        _df.index, _df["min"], _df["max"], color=colors.to_rgba(fill_color, fill_opacity)
    )
    ax.axvspan(0, 8, color="lightgray", alpha=0.2, lw=0)
    ax.axvspan(20, 23, color="lightgray", alpha=0.2, lw=0)
    ax.set_xlabel("Hours")
    ax.set_ylabel(f"{field.short} ({field.unit})")

    if output_dir:
        filename = (
            f"{title}-{chart_id or str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S'))}-{month + 1}"
        )
        fig.savefig(
            f"{output_dir}/{filename}.{output_format}", format=output_format, transparent=True
        )
        logger(f"Output saved to {filename}")
    if close:
        plt.close(fig)
    return fig, ax


def multi_series_annual_graph(
    title: str,
    output_data: List[List[Any]],
    series_titles: List[str],
    field: Field,
    start_day: int = 0,
    end_day: int = 365,
    output_dir: Optional[str | Path] = None,
    chart_id: Optional[str | int] = None,
    label_x_days: int = 1,
    average_step: int = 1,
    linestyles: List[dict] = [],
    output_format: str = "png",
    top_y_lim: Optional[float] = None,
    bottom_y_lim: float = 0,
    log_level: int = 1,
    FONTSIZE: int = 8,
    *args,
    **kwargs,
) -> Figure:
    """Create a multi series graph between two dates.


    # assume fields map to outputdata

    Parameters
    ----------
    title : str
        [description]
    output_data : List[List[Any]]
        [description]
    series_titles : List[str]
        [description]
    field: Field
        [description]
    start_day : int, optional
        [description], by default 0
    end_day : int, optional
        [description], by default 365
    output_dir : str, optional
        [description], by default None
    chart_id : str, optional
        [description], by default None
    label_x_days : int, optional
        [description], by default 1
    average_step: int, optional
        the average step for moving average graph, default 1 (i.e ignore)
    *args,**kwargs
        Passed to Figure class


    Returns
    -------
    Figure
        Matplotlib figure

    Raises
    ------
    ValueError
        If input is invalid
    """
    logger = log(log_level)
    if not isinstance(series_titles, list):
        raise ValueError("series_titles should be a list of fields")
    plt.rc("font", size=FONTSIZE)

    fig = plt.figure(*args, **kwargs)

    fig.suptitle(title)

    # Get x ticks
    t_count = ceil((end_day - start_day) / label_x_days)
    t_spacing = (end_day - start_day) / t_count
    labels = [f"{int(i * t_spacing) + start_day}" for i in range(t_count + 1)]
    labelxpos = [i * 24 * t_spacing for i in range(t_count + 1)]

    axes = fig.add_subplot(1, 4, (1, 3))

    for i, series in enumerate(series_titles):
        try:
            linestyle = (
                linestyles[i]["style"] if linestyles and "style" in linestyles[i] else "solid"
            )
            linewidth = (
                linestyles[i]["width"]
                if linestyles and "width" in linestyles[i]
                else 0.5 + (i * 0.1) / len(series_titles)
            )
            zorder = linestyles[i]["zorder"] if linestyles and "zorder" in linestyles[i] else i
            plot_data = (
                output_data[i]
                if average_step == 1
                else calc_moving_average(output_data[i], average_step)
            )

            axes.plot(
                plot_data,
                label=series,
                linestyle=linestyle,
                zorder=zorder,
                linewidth=linewidth,
            )
        except Exception as e:
            print(e)
            warnings.warn(f"Failed to plot {series}")

    axes.set_xticks(labelxpos)
    axes.set_xticklabels(labels, fontsize=FONTSIZE)
    axes.set_xlabel("Days", fontsize=FONTSIZE)
    axes.set_ylabel("y" if not field else f"{field.short} ({field.unit})", fontsize=FONTSIZE)
    axes.yaxis.grid(b=True)
    if top_y_lim and bottom_y_lim:
        axes.set_ylim(top=top_y_lim, bottom=bottom_y_lim)
    elif top_y_lim:
        axes.set_ylim(top=top_y_lim, bottom=0)
    else:
        bottom_y_lim = bottom_y_lim if -9999999 < bottom_y_lim < 99999999 else 0
        axes.set_ylim(bottom=bottom_y_lim or 0)

    if len(series_titles) < 10:
        axes.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize="xx-small")
    else:
        axes.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize="xx-small")
    plt.tight_layout()

    if output_dir:
        filename = f"{title}-{chart_id if chart_id is not None else str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S'))}"
        fig.savefig(
            f"{output_dir}/{filename}.{output_format}", format=output_format, transparent=True
        )
        logger(f"Output saved to {filename}")

    return fig


if __name__ == "__main__":
    import os
    from pyDO3SE.Output.Output_Shape import output_fields_map

    data = [0, 1, 2, 3, 4, 5, 29, 7, 8, 9, 10, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 2]
    fig = annual_graph(
        [data[i % 24] for i in list(range(24 * 6))] + [0],
        output_fields_map["gsto"],
        output_dir=os.getcwd(),
        chart_id=1,
        start_day=0,
        end_day=12,
        label_x_days=3,
        average_step=12,
    )

    fig = annual_graph(
        [data[i % 24] for i in list(range(24 * 6))] + [0],
        output_fields_map["gsto"],
        output_dir=os.getcwd(),
        chart_id=1,
        start_day=0,
        end_day=12,
        label_x_days=3,
        average_step=1,
    )

    series = ["gsto_l", "A_n", "PARsun", "PARsun", "PARsun", "PARsun"]
    demo_data = [[j * (i + 1) for j in list(range(0, 24 * 6))] + [0] for i in range(6)]
    multi_series_annual_graph(
        "demo",
        demo_data,
        series,
        output_fields_map["gsto_l"],
        output_dir=os.getcwd(),
        chart_id=2,
        start_day=0,
        end_day=6,
        label_x_days=3,
        linestyles=[
            {
                "style": "solid" if f != "A_n" else "dashed",
                "width": i + 1 if i < 2 else 0.2,
            }
            for i, f in enumerate(series)
        ],
    )

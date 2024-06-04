"""Create charts of pypDO3SE output."""

# TODO: Should import output data

from math import ceil
from typing import Any, List
import warnings
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from datetime import datetime


from pyDO3SE.util.Objects import Field


# plt.style.use('pyDO3SE/Analysis/style.mplstyle')


def calc_moving_average(data: List[float], size: int) -> List[float]:
    if len(data) == 0:
        raise ValueError(f"Data is empty")
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
    data: List[Any],
    field: Field,
    start_day: int = 0,
    end_day: int = 365,
    output_dir: str = None,
    chart_id: str = None,
    label_x_days: int = 1,
    average_step: int = 1,
    output_format: str = 'png',
    log_level: int = 1,
    *args,
    **kwargs,
) -> Figure:
    """Create a graph of a single series between two dates.

    Parameters
    ----------
    data : List[Any]
        [description]
    field : Union[Field, List[Field]]
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
    moving_average: int, optional
        [description], by default 1

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
    fig = plt.figure(*args, **kwargs)
    title = ''

    # 9 days  with 3 day spacing
    # Get x ticks
    t_count = ceil((end_day - start_day) / label_x_days)
    t_spacing = (end_day - start_day) / t_count

    labels = [f'{int(i * t_spacing) + start_day}' for i in range(t_count + 1)]
    labelx = [i * 24 * t_spacing for i in range(t_count + 1)]

    if isinstance(field, list):
        # assume fields map to outputdata lists
        raise ValueError('Use multi_series_annual_graph for multi series graph')
    # Assume only a single data series
    title = field.short
    fig.suptitle(field.short)
    axes = fig.add_subplot(1, 1, 1)
    plot_data = data if average_step == 1 else calc_moving_average(data, average_step)
    axes.plot(plot_data)
    axes.set_xticks(labelx)
    axes.set_xticklabels(labels)
    axes.set_xlabel('Days')
    axes.set_ylabel(f'{field.short} ({field.unit})')
    axes.yaxis.grid(b=True)

    if output_dir:
        filename = f'{title}-{chart_id or str(datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))}'
        fig.savefig(f'{output_dir}/{filename}.{output_format}',
                    format=output_format, transparent=True)
        logger(f'Output saved to {filename}')

    return fig


def multi_series_annual_graph(
    title: str,
    output_data: List[List[Any]],
    series_titles: List[str],
    field: Field,
    start_day: int = 0,
    end_day: int = 365,
    output_dir: str = None,
    chart_id: str = None,
    label_x_days: int = 1,
    average_step: int = 1,
    linestyles: List[dict] = [],
    output_format: str = 'png',
    top_y_lim: float = None,
    bottom_y_lim: float = 0,
    log_level: int = 1,
    FONTSIZE: int = 8 ,
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
        raise ValueError('series_titles should be a list of fields')
    plt.rc('font', size=FONTSIZE)

    fig = plt.figure(*args, **kwargs)

    fig.suptitle(title)

    # Get x ticks
    t_count = ceil((end_day - start_day) / label_x_days)
    t_spacing = (end_day - start_day) / t_count
    labels = [f'{int(i * t_spacing) + start_day}' for i in range(t_count + 1)]
    labelxpos = [i * 24 * t_spacing for i in range(t_count + 1)]

    axes = fig.add_subplot(1, 4, (1, 3))

    for i, series in enumerate(series_titles):
        try:
            linestyle = linestyles[i]['style'] if linestyles and 'style' in linestyles[i] \
                else 'solid'
            linewidth = linestyles[i]['width'] if linestyles and 'width' in linestyles[i] \
                else 0.5 + (i * 0.1) / len(series_titles)
            zorder = linestyles[i]['zorder'] if linestyles and 'zorder' in linestyles[i] \
                else i
            plot_data = output_data[i] if average_step == 1 else calc_moving_average(
                output_data[i], average_step)

            axes.plot(
                plot_data,
                label=series,
                linestyle=linestyle,
                zorder=zorder,
                linewidth=linewidth,
            )
        except Exception as e:
            warnings.warn(f"Failed to plot {series}")

    axes.set_xticks(labelxpos)
    axes.set_xticklabels(labels, fontsize=FONTSIZE)
    axes.set_xlabel('Days', fontsize=FONTSIZE)
    axes.set_ylabel('y' if not field else f'{field.short} ({field.unit})', fontsize=FONTSIZE)
    axes.yaxis.grid(b=True)
    if top_y_lim and bottom_y_lim:
        axes.set_ylim(top=top_y_lim, bottom=bottom_y_lim)
    elif top_y_lim:
        axes.set_ylim(top=top_y_lim, bottom=0)
    else:
        bottom_y_lim = bottom_y_lim if -9999999 < bottom_y_lim < 99999999 else 0
        axes.set_ylim(bottom=bottom_y_lim or 0)

    if len(series_titles) < 10:
        axes.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='xx-small')
    else:
        axes.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='xx-small')
    plt.tight_layout()

    if output_dir:
        filename = f'{title}-{chart_id if chart_id is not None else str(datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))}'
        fig.savefig(f'{output_dir}/{filename}.{output_format}',
                    format=output_format, transparent=True)
        logger(f'Output saved to {filename}')

    return fig


if __name__ == "__main__":
    import os
    from pyDO3SE.Output.Output_Shape import output_fields_map
    data = [0, 1, 2, 3, 4, 5, 29, 7, 8, 9, 10, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 2]
    fig = annual_graph(
        [data[i % 24] for i in list(range(24 * 6))] + [0],
        output_fields_map['gsto'],
        output_dir=os.getcwd(),
        chart_id=1,
        start_day=0,
        end_day=12,
        label_x_days=3,
        average_step=12
    )

    fig = annual_graph(
        [data[i % 24] for i in list(range(24 * 6))] + [0],
        output_fields_map['gsto'],
        output_dir=os.getcwd(),
        chart_id=1,
        start_day=0,
        end_day=12,
        label_x_days=3,
        average_step=1
    )

    series = ['gsto_l', 'A_n', 'PARsun', 'PARsun', 'PARsun', 'PARsun']
    demo_data = [
        [j * (i + 1) for j in list(range(0, 24 * 6))] + [0]
        for i in range(6)]
    multi_series_annual_graph(
        'demo',
        demo_data,
        series,
        output_fields_map['gsto_l'],
        output_dir=os.getcwd(),
        chart_id=2,
        start_day=0,
        end_day=6,
        label_x_days=3,
        linestyles=[{
            'style': 'solid' if f != 'A_n' else 'dashed',
            'width': i + 1 if i < 2 else 0.2,
        } for i, f in enumerate(series)],
    )

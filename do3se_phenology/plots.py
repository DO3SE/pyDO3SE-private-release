from itertools import count
from pathlib import Path
from typing import List, Optional
import numpy as np

from do3se_phenology.carbon_allocation import calc_partition_coefficients
from do3se_phenology.f_phen import f_phen_simple_PLF_range, f_phen_simple_PLF_value
from do3se_phenology.switchboard import process_phenology_config
from do3se_phenology.canopy_structure import get_growing_populations_range_from_config
from do3se_phenology.phyllochron_dvi import get_dvi_range_from_species_config
from do3se_phenology.config import (
    ModelConfig,
    SpeciesConfig,
    DVIMethods,
    LeafFPhenMethods,
)
from do3se_phenology.utils import generate_days_data, generate_example_td_data

try:
    from matplotlib import pyplot as plt
    from matplotlib import colormaps as mpcm
    import matplotlib.patches as patches
    import matplotlib.figure
    import matplotlib.axes
except:
    plt = lambda *args, **kwargs: Exception("matplotlib required to draw graphs!")


def plot_f_phen_dd_box_data(
    f_phen_1: int,
    f_phen_4: int,
    f_phen_a: float,
    f_phen_c: float,
    f_phen_e: float,
    SGS: int,
    EGS: int,
    ax: matplotlib.axes.Axes,
    box_y_start=0,
    **kwargs,
):
    EGS_adj = EGS if EGS > SGS else EGS + 365
    COLOR = "blue"
    ax.axvline(SGS, linestyle="dotted", color=COLOR)

    y = count(box_y_start, 2)
    text_offset = (EGS_adj - SGS) / 100
    draw_range(
        ax,
        SGS,
        length=f_phen_1,
        label="f_phen_1",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        EGS_adj - f_phen_4,
        length=f_phen_4,
        label="f_phen_4",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    ax.scatter([], [], label="f_phen")
    return ax


def plot_f_phen_tt_box_data(
    t_sgs,
    t_egs,
    t_f_phen_a,
    t_f_phen_b,
    t_f_phen_c,
    t_f_phen_d,
    ax: matplotlib.axes.Axes,
    box_y_start=0,
    **kwargs,
):
    COLOR = "blue"
    ax.axvline(t_sgs, linestyle="dotted", color=COLOR)

    y = count(box_y_start, 2)
    text_offset = (t_egs - t_sgs) / 100
    draw_range(
        ax,
        t_sgs,
        length=t_f_phen_a,
        label="t_f_phen_a",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs,
        length=t_f_phen_b,
        label="t_f_phen_b",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs,
        length=t_f_phen_c,
        label="t_f_phen_c",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    # draw_range(ax, t_sgs, length=t_f_phen_d,
    #            label="t_f_phen_d", COLOR=COLOR, text_offset=text_offset, y=next(y))
    ax.scatter([], [], label="f_phen")
    return ax


def plot_f_phen_td_data(
    accumulated_temperatures,
    f_phen_data,
    t_sgs,
    t_egs,
    t_f_phen_a,
    t_f_phen_b,
    t_f_phen_c,
    t_f_phen_d,
    ax: matplotlib.axes.Axes,
    **kwargs,
):
    TEXT_HEIGHT = 2.5
    TEXT_ANGLE = 90
    COLOR = "blue"

    ax.plot(
        accumulated_temperatures,
        f_phen_data,
        label="f_phen_data",
        color=COLOR,
        **kwargs,
    )
    ax.axvline(t_sgs, linestyle="dotted", color=COLOR)
    ax.axvline(t_sgs + t_f_phen_a, linestyle="dotted", color=COLOR)
    ax.axvline(t_sgs + t_f_phen_b, linestyle="dotted", color=COLOR)
    ax.axvline(t_sgs + t_f_phen_c, linestyle="dotted", color=COLOR)
    ax.axvline(t_sgs + t_f_phen_d, linestyle="dotted", color=COLOR)
    ax.text(t_sgs, TEXT_HEIGHT, "SGS", rotation=TEXT_ANGLE, color=COLOR)
    ax.text(t_egs, TEXT_HEIGHT, "EGS", rotation=TEXT_ANGLE, color=COLOR)
    ax.text(
        t_sgs + t_f_phen_a,
        TEXT_HEIGHT,
        "Plant Emergence",
        rotation=TEXT_ANGLE,
        color=COLOR,
    )
    ax.text(
        t_sgs + t_f_phen_b, TEXT_HEIGHT, "f_phen_b", rotation=TEXT_ANGLE, color=COLOR
    )

    return ax


def draw_range(
    ax,
    start,
    end=None,
    label=None,
    length=None,
    COLOR=None,
    y=0.0,
    end_height=0.2,
    text_offset=1.0,
):
    end = end or start + length
    ax.axvline(start, linestyle="dotted", color=COLOR)
    ax.axvline(end, linestyle="dotted", color=COLOR)
    ax.plot(
        [start, start],
        [y - end_height, y + end_height],
        color=COLOR,
        linestyle="dotted",
    )
    ax.plot(
        [end, end], [y - end_height, y + end_height], color=COLOR, linestyle="dotted"
    )
    ax.text(end + text_offset, y, label, color=COLOR)
    ax.plot([start, end], [y, y], color=COLOR)


def plot_leaf_f_phen_dd_box_plot(
    leaf_f_phen_1,
    leaf_f_phen_2,
    astart,
    sgs,
    egs,
    ax: matplotlib.axes.Axes,
    box_y_start=0.2,
    **kwargs,
):
    egs_adg = egs if egs > sgs else egs + 365
    COLOR = "red"
    ax.axvline(sgs, linestyle="dotted", color=COLOR)
    y = count(box_y_start, 2)
    text_offset = (egs_adg - sgs) / 100
    draw_range(
        ax,
        astart,
        length=leaf_f_phen_1,
        label="leaf_f_phen_1",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
        **kwargs,
    )
    draw_range(
        ax,
        egs_adg - leaf_f_phen_2,
        length=leaf_f_phen_2,
        label="leaf_f_phen_2",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
        **kwargs,
    )
    ax.scatter([], [], label="leaf_f_phen")
    return ax


def plot_leaf_f_phen_td_box_plot(
    t_leaf_f_phen_a,
    t_leaf_f_phen_b,
    t_leaf_f_phen_e,
    t_leaf_f_phen_g,
    t_leaf_f_phen_h,
    t_leaf_f_phen_i,
    t_astart,
    t_sgs,
    t_egs,
    ax: matplotlib.axes.Axes,
    box_y_start=0.2,
    **kwargs,
):
    COLOR = "red"
    ax.axvline(t_sgs, linestyle="dotted", color=COLOR)

    y = count(box_y_start, 2)
    text_offset = (t_egs - t_sgs) / 100
    draw_range(
        ax,
        t_sgs,
        length=t_astart,
        label="sgs -> astart",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs + t_astart,
        length=t_leaf_f_phen_e,
        label="leaf_f_phen_e",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs + t_astart + t_leaf_f_phen_e,
        length=t_leaf_f_phen_g,
        label="leaf_f_phen_g",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs + t_astart + t_leaf_f_phen_e,
        length=t_leaf_f_phen_h,
        label="leaf_f_phen_h",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs + t_astart + t_leaf_f_phen_e,
        length=t_leaf_f_phen_i,
        label="leaf_f_phen_i",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    ax.scatter([], [], label="leaf_f_phen")
    return ax


def plot_leaf_f_phen_data(
    x_data,
    leaf_f_phen_data,
    t_sgs,
    ax: matplotlib.axes.Axes,
    **kwargs,
):
    COLOR = "red"
    ax.plot(x_data, leaf_f_phen_data, label="leaf_f_phen_data", color=COLOR, **kwargs)
    ax.axvline(t_sgs, linestyle="dotted", color=COLOR)
    return ax


def plot_f_phen_dd_vlines(
    sgs,
    egs,
    astart,
    ax: matplotlib.axes.Axes,
    TEXT_HEIGHT=1.0,
    TEXT_ANGLE=90,
    COLOR="blue",
    plot_labels=True,
    **kwargs,
):

    egs_adj = egs if egs > sgs else egs + 365

    if plot_labels:
        ax.annotate("SGS", (sgs, TEXT_HEIGHT), rotation=TEXT_ANGLE, color=COLOR)
        ax.annotate("EGS", (egs_adj, TEXT_HEIGHT), rotation=TEXT_ANGLE, color=COLOR)
        ax.annotate("Astart", (astart, TEXT_HEIGHT), rotation=TEXT_ANGLE, color=COLOR) if astart else None
    ax.axvline(sgs, linestyle="dotted", color=COLOR)
    ax.axvline(egs_adj, linestyle="dotted", color=COLOR)
    ax.axvline(astart, linestyle="dotted", color=COLOR) if astart else None


def plot_f_phen_dd_data(
    x_data,
    f_phen_data,
    sgs,
    ax: matplotlib.axes.Axes,
    **kwargs,
):
    COLOR = "blue"
    ax.plot(x_data, f_phen_data, label="f_phen_data", color=COLOR, **kwargs)
    ax.axvline(sgs, linestyle="dotted", color=COLOR)
    return ax


def plot_ewert_phenology_data_box_plot(
    t_lem,
    t_lma,
    t_lep,
    t_lse,
    t_sgs,
    t_egs,
    t_Astart,
    ax: matplotlib.axes.Axes,
    box_y_start=0,
    **kwargs,
):
    COLOR = "green"

    ax.axvline(t_sgs, linestyle="dotted", color=COLOR)

    y = count(box_y_start, 2)
    text_offset = (t_egs - t_sgs) / 100
    draw_range(
        ax,
        t_sgs + t_Astart - t_lem,
        length=t_lem,
        label="t_lem",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs + t_Astart,
        length=t_lma,
        label="t_lma",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs + t_Astart,
        length=t_lep,
        label="t_lep",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    draw_range(
        ax,
        t_sgs + t_Astart + t_lep,
        length=t_lse,
        label="t_lse",
        COLOR=COLOR,
        text_offset=text_offset,
        y=next(y),
    )
    ax.scatter([], [], label="ewert_phenology")

    return ax


def plot_key_dates(
    accumulated_temperatures,
    f_LA_data,
    f_LS_data,
    t_lem,
    t_lma,
    t_lep,
    t_lse,
    t_sgs,
    t_egs,
    t_Astart,
    ax: matplotlib.axes.Axes,
    **kwargs,
):
    return ax


def plot_ewert_phenology_data(
    accumulated_temperatures,
    f_LA_data,
    f_LS_data,
    t_lem,
    t_lma,
    t_lep,
    t_lse,
    t_sgs,
    t_egs,
    t_Astart,
    ax: matplotlib.axes.Axes,
    **kwargs,
):
    TEXT_ANGLE = 90
    COLOR = "green"
    ax.plot(accumulated_temperatures, f_LA_data, label="f_LA", color=COLOR, **kwargs)

    ax.plot(
        accumulated_temperatures, f_LS_data, "--", label="f_LS", color=COLOR, **kwargs
    )

    ax.axvline(t_sgs, linestyle="dotted", color=COLOR)
    ax.axvline(t_sgs + t_Astart, linestyle="dotted", color=COLOR)
    ax.axvline(t_sgs + t_Astart - t_lem, linestyle="dotted", color=COLOR)

    ax.axvline(t_sgs + t_Astart + t_lma, linestyle="dotted", color=COLOR)
    ax.axvline(t_sgs + t_Astart + t_lep, linestyle="dotted", color=COLOR)

    ax.text(t_sgs, 2.5, "SGS", rotation=TEXT_ANGLE)
    ax.text(t_sgs + t_Astart - t_lem, 2.5, "Flag Emergence", rotation=TEXT_ANGLE)
    ax.text(t_sgs + t_Astart + t_lep, 2.5, "Flag Onset senes", rotation=TEXT_ANGLE)
    ax.text(t_sgs + t_Astart, 2.5, "AStart", rotation=TEXT_ANGLE)
    ax.text(t_egs, 2.5, "EGS", rotation=TEXT_ANGLE)

    return ax


def plot_growing_fractions_from_config(
    species_config,
    td,
    dvi,
    nP,
    a_root=18.5,
    a_leaf=16.0,
    a_stem=18.0,
    b_root=-20.0,
    b_leaf=-15.0,
    b_stem=-18.5,
    ax: matplotlib.axes.Axes = plt.axes(),
):
    ax.set_title("Growing populations")
    growing_populations_all, emerged_leaf_populations_count = (
        get_growing_populations_range_from_config(species_config, nP, td)
    )

    p_root, p_leaf, p_stem, p_harv = list(
        zip(
            *[
                calc_partition_coefficients(
                    d,
                    a_root,
                    a_leaf,
                    a_stem,
                    b_root,
                    b_leaf,
                    b_stem,
                )
                for d in dvi
            ]
        )
    )
    leaf_carbon_fracs = [
        [
            sum(g[i + 1 :]) * p / sum(g)
            if ff == 0 or sum(g) == 0
            else p * (ff + sum(g[i + 1 :])) / sum(g)
            if i < len(g)
            else (p * ff) / sum(g)
            for i, ff in enumerate(g)
        ]
        for g, p in zip(growing_populations_all, p_leaf)
    ]

    colors = [i / nP for i in range(nP)]
    # cmap=plt.cm.bgrcmyk
    cmap = mpcm.get_cmap("Spectral")
    c = cmap(colors)

    for iP in range(nP):
        data = [d[iP] for d in leaf_carbon_fracs]
        ax.plot(td, data, label=f"{iP} leaf_carbon_fractions", c=c[iP])
    return ax


def plot_carbon_fractions_from_species_config(
    species_config,
    td,
    dvi,
    nP,
    a_root=18.5,
    a_leaf=16.0,
    a_stem=18.0,
    b_root=-20.0,
    b_leaf=-15.0,
    b_stem=-18.5,
    ax: matplotlib.axes.Axes = plt.axes(),
):
    p_root, p_leaf, p_stem, p_harv = list(
        zip(
            *[
                calc_partition_coefficients(
                    d,
                    a_root,
                    a_leaf,
                    a_stem,
                    b_root,
                    b_leaf,
                    b_stem,
                )
                for d in dvi
            ]
        )
    )

    ax.set_title("Carbon distribution fractions")
    ax.plot(td, p_root, label="p_root")
    ax.plot(td, p_leaf, label="p_leaf")
    ax.plot(td, p_stem, label="p_stem")
    ax.plot(td, p_harv, label="p_harv")
    return ax


class _PhenologyPlotClass:
    speciesConfig: SpeciesConfig
    td_data: List[float]
    dd_data: List[int]

    def __init__(
        self,
        speciesConfig: SpeciesConfig,
        td_data: List[float],
        dd_data: List[int],
        td_plot_args: dict = {},
        dd_plot_args: dict = {},
    ) -> None:
        self.speciesConfig = speciesConfig
        self.td_data = td_data
        self.dd_data = dd_data
        self.td_plot_args = td_plot_args
        self.dd_plot_args = dd_plot_args


class PlotLeafFPhenDataFromConfig(_PhenologyPlotClass):
    def plot_td(self, ax: matplotlib.axes.Axes):
        species_config = self.speciesConfig
        td_data = self.td_data
        assert species_config.leaf_fphen_intervals is not None, (
            "Leaf f_phen intervals not defined in species config"
        )
        leaf_f_phen_data_td = np.interp(
            td_data,
            [i[0] for i in species_config.leaf_fphen_intervals],
            [i[1] for i in species_config.leaf_fphen_intervals],
        )

        plot_leaf_f_phen_data(
            td_data,
            leaf_f_phen_data_td,
            species_config.key_dates_td.sowing,
            ax=ax,
        )

    def plot_dd(self, ax: matplotlib.axes.Axes):
        # TODO: Get lead_f_phen data from dd intervals instead
        leaf_f_phen_intervals = self.speciesConfig.leaf_fphen_intervals
        assert leaf_f_phen_intervals is not None, (
            "Leaf f_phen intervals not defined in species config"
        )
        leaf_f_phen_data_dd = np.interp(
            self.dd_data - self.speciesConfig.key_dates.Astart +1,
            leaf_f_phen_intervals[0],
            leaf_f_phen_intervals[1],
        )

        plot_leaf_f_phen_data(
            self.dd_data,
            leaf_f_phen_data_dd,
            self.speciesConfig.key_dates.sowing,
            ax=ax,
        )


class PlotLeafFPhenDataBoxPlotFromConfig(_PhenologyPlotClass):
    def plot_td(self, ax: matplotlib.axes.Axes):
        species_config = self.speciesConfig

        plot_leaf_f_phen_td_data_box_plot_from_config(
            species_config,
            ax=ax,
            **self.td_plot_args,
        )

    def plot_dd(self, ax: matplotlib.axes.Axes):
        species_config = self.speciesConfig
        plot_leaf_f_phen_dd_data_box_plot_from_config(
            species_config,
            ax=ax,
            **self.dd_plot_args,
        )


class PlotFPhenDataFromConfig(_PhenologyPlotClass):
    def plot_td(self, ax: matplotlib.axes.Axes):
        species_config = self.speciesConfig
        td_data = self.td_data

        plot_f_phen_td_data_from_config(
            species_config,
            td_data,
            ax=ax,
            **self.td_plot_args,
        )

    def plot_dd(self, ax: matplotlib.axes.Axes):
        species_config = self.speciesConfig
        dd_data = self.dd_data

        plot_f_phen_dd_data_from_config(
            species_config,
            dd_data,
            ax=ax,
            **self.dd_plot_args,
        )


class PlotFPhenDataBoxPlotFromConfig(_PhenologyPlotClass):
    def plot_td(self, ax: matplotlib.axes.Axes, **kwargs):
        species_config = self.speciesConfig

        plot_f_phen_tt_data_box_plot_from_config(
            species_config, ax=ax, **kwargs, **self.td_plot_args
        )

    def plot_dd(self, ax: matplotlib.axes.Axes, **kwargs):
        species_config = self.speciesConfig
        plot_f_phen_dd_data_box_plot_from_config(
            species_config, ax=ax, **kwargs, **self.dd_plot_args
        )


class PlotEwertPhenologyDataFromConfig(_PhenologyPlotClass):
    def plot_td(self, ax: matplotlib.axes.Axes):
        plot_ewert_phenology_data_from_config(
            self.speciesConfig,
            self.td_data,
            ax=ax,
            **self.td_plot_args,
        )

    def plot_dd(self, ax: matplotlib.axes.Axes):
        raise NotImplementedError(
            "Ewert phenology box plot not implemented for DD data yet."
        )


class PlotEwertPhenologyDataBoxPlotFromConfig(_PhenologyPlotClass):
    def plot_td(self, ax: matplotlib.axes.Axes):
        species_config = self.speciesConfig

        plot_ewert_phenology_data_box_plot_from_config(
            species_config,
            ax=ax,
        )

    def plot_dd(self, ax: matplotlib.axes.Axes):
        raise NotImplementedError(
            "Ewert phenology box plot not implemented for DD data yet."
        )


def plot_leaf_f_phen_td_data_box_plot_from_config(
    species_config: SpeciesConfig,
    ax: matplotlib.axes.Axes,
    box_y_start=0.8,
):
    plot_leaf_f_phen_td_box_plot(
        species_config.leaf_f_phen_a,
        species_config.leaf_f_phen_b,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_e,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_g,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_h,
        species_config.key_lengths_flag_leaf_td.leaf_f_phen_i,
        species_config.key_dates_td.Astart,
        species_config.key_dates_td.sowing,
        species_config.key_dates_td.harvest,
        ax=ax,
        box_y_start=box_y_start,
    )


def plot_leaf_f_phen_dd_data_box_plot_from_config(
    species_config: SpeciesConfig,
    ax: matplotlib.axes.Axes,
    box_y_start=0.2,
):
    plot_leaf_f_phen_dd_box_plot(
        species_config.day_fphen_plf.leaf_f_phen_1,
        species_config.day_fphen_plf.leaf_f_phen_2,
        species_config.key_dates.Astart,
        species_config.key_dates.sowing,
        species_config.key_dates.harvest,
        ax=ax,
        box_y_start=box_y_start,
    )


def plot_f_phen_dd_data_from_config(
    species_config: SpeciesConfig, dd_data, ax: matplotlib.axes.Axes
):
    assert species_config.day_fphen_plf.f_phen_1 is not None
    assert species_config.day_fphen_plf.f_phen_4 is not None
    assert species_config.day_fphen_plf.f_phen_a is not None
    assert species_config.day_fphen_plf.f_phen_c is not None
    assert species_config.day_fphen_plf.f_phen_e is not None
    assert species_config.key_dates.sowing is not None
    assert species_config.key_dates.harvest is not None
    f_phen_values = f_phen_simple_PLF_range(
        dd_data,
        f_phen_1=species_config.day_fphen_plf.f_phen_1,
        f_phen_4=species_config.day_fphen_plf.f_phen_4,
        f_phen_a=species_config.day_fphen_plf.f_phen_a,
        f_phen_c=species_config.day_fphen_plf.f_phen_c,
        f_phen_e=species_config.day_fphen_plf.f_phen_e,
        SGS=int(species_config.key_dates.sowing),
        EGS=int(species_config.key_dates.harvest),
        f_phen_min=species_config.f_phen_min or 0.0,
    )
    plot_f_phen_dd_data(
        dd_data,
        f_phen_values,
        sgs=species_config.key_dates.sowing,
        ax=ax,
    )
    plot_f_phen_dd_vlines(
        species_config.key_dates.sowing,
        species_config.key_dates.harvest,
        species_config.key_dates.Astart,
        ax=ax,
        TEXT_HEIGHT=1.0,
    )


def plot_f_phen_td_data_from_config(
    species_config: SpeciesConfig, td_data, ax: matplotlib.axes.Axes, **kwargs
):
    assert species_config.fphen_intervals is not None, (
        "fphen_intervals not defined in species config"
    )
    f_phen_data = np.interp(
        td_data,
        [i[0] for i in species_config.fphen_intervals],
        [i[1] for i in species_config.fphen_intervals],
    )

    plot_f_phen_td_data(
        td_data,
        f_phen_data,
        species_config.key_dates_td.sowing,
        species_config.key_dates_td.harvest,
        species_config.key_lengths_td.sowing_to_emerge,
        species_config.key_lengths_td.sowing_to_f_phen_b,
        species_config.key_lengths_td.sowing_to_f_phen_c,
        species_config.key_lengths_td.sowing_to_end,
        ax=ax,
        **kwargs,
    )


def plot_f_phen_tt_data_box_plot_from_config(
    species_config: SpeciesConfig, ax: matplotlib.axes.Axes
):
    plot_f_phen_tt_box_data(
        species_config.key_dates_td.sowing,
        species_config.key_dates_td.harvest,
        species_config.key_lengths_td.sowing_to_emerge,
        species_config.key_lengths_td.sowing_to_f_phen_b,
        species_config.key_lengths_td.sowing_to_f_phen_c,
        species_config.key_lengths_td.sowing_to_end,
        ax=ax,
    )


def plot_f_phen_dd_data_box_plot_from_config(
    species_config: SpeciesConfig, ax: matplotlib.axes.Axes, **kwargs
):
    assert species_config.day_fphen_plf.f_phen_1 is not None, (
        "species_config.day_fphen_plf.f_phen_1 is None!"
    )
    assert species_config.day_fphen_plf.f_phen_4 is not None, (
        "species_config.day_fphen_plf.f_phen_4 is None!"
    )
    assert species_config.day_fphen_plf.f_phen_a is not None, (
        "species_config.day_fphen_plf.f_phen_a is None!"
    )
    assert species_config.day_fphen_plf.f_phen_c is not None, (
        "species_config.day_fphen_plf.f_phen_c is None!"
    )
    assert species_config.day_fphen_plf.f_phen_e is not None, (
        "species_config.day_fphen_plf.f_phen_e is None!"
    )
    assert species_config.key_dates.sowing is not None, (
        "species_config.key_dates.sowing is None!"
    )
    assert species_config.key_dates.harvest is not None, (
        "species_config.key_dates.harvest is None!"
    )

    plot_f_phen_dd_box_data(
        species_config.day_fphen_plf.f_phen_1,
        species_config.day_fphen_plf.f_phen_4,
        species_config.day_fphen_plf.f_phen_a,
        species_config.day_fphen_plf.f_phen_c,
        species_config.day_fphen_plf.f_phen_e,
        int(species_config.key_dates.sowing),
        int(species_config.key_dates.harvest),
        ax=ax,
        **kwargs,
    )
    plot_f_phen_dd_vlines(
        species_config.key_dates.sowing,
        species_config.key_dates.harvest,
        species_config.key_dates.Astart,
        ax=ax,
        plot_labels=False,
    )


def plot_growing_populations_from_config(
    species_config: SpeciesConfig,
    nP: int,
    td: List[float],
    row_height=0.15,
    base_offset=1,
    ax: matplotlib.axes.Axes = plt.axes(),
):
    growing_populations_all, emerged_leaf_populations_count = (
        get_growing_populations_range_from_config(species_config, nP, td)
    )

    colors = [i / nP for i in range(nP)]
    # cmap=plt.cm.bgrcmyk
    cmap = mpcm.get_cmap("Spectral")
    c = cmap(colors)
    for iP in range(nP):
        # growing = np.array([g[iP] for g in growing_populations_all])
        growing_boundary = [i for i, g in zip(td, growing_populations_all) if g[iP]]
        base_point = (growing_boundary[0], base_offset + (iP * row_height))
        width = growing_boundary[-1] - growing_boundary[0]
        L = ax.plot([0], label=f"{iP} is growing", c=c[iP])
        r = patches.Rectangle(base_point, width, row_height, facecolor=L[0].get_color())
        ax.add_patch(r)

    ax.yaxis.set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    ax.spines["left"].set_visible(False)

    # ax.legend(rects, labels)


def plot_ewert_phenology_data_from_config(
    species_config: SpeciesConfig, td, ax: matplotlib.axes.Axes
):
    td_dd_list = [t - species_config.key_dates_td.sowing for t in td]
    t_astart = species_config.key_lengths_td.sowing_to_astart
    assert t_astart is not None, "t_astart is None!"
    t_lma = species_config.key_lengths_flag_leaf_td.tl_ma
    assert t_lma is not None, "t_lma is None!"
    t_lep = species_config.key_lengths_flag_leaf_td.tl_ep
    assert t_lep is not None, "t_lep is None!"
    t_lem = species_config.key_lengths_flag_leaf_td.tl_em
    assert t_lem is not None, "t_lem is None!"
    f_LA = [max(0, min(1, 1 - (td_dd - t_astart) / (t_lma))) for td_dd in td_dd_list]
    growing = [
        max(0, min(1, (td_dd - (t_astart - t_lem)) / (t_lem))) for td_dd in td_dd_list
    ]
    fO3_l = 1
    f_LS = [
        max(0, min(1, 1 - ((td_dd - t_astart - t_lep) / (t_lma / fO3_l - t_lep))))
        for td_dd in td_dd_list
    ]

    ax.plot(td, growing, label="flag_growing")

    plot_ewert_phenology_data(
        td,
        f_LA,
        f_LS,
        species_config.key_lengths_flag_leaf_td.tl_em,
        species_config.key_lengths_flag_leaf_td.tl_ma,
        species_config.key_lengths_flag_leaf_td.tl_ep,
        species_config.key_lengths_flag_leaf_td.tl_se,
        species_config.key_dates_td.sowing,
        species_config.key_dates_td.harvest,
        species_config.key_dates_td.Astart,
        ax=ax,
    )


def plot_ewert_phenology_data_box_plot_from_config(
    species_config: SpeciesConfig, ax: matplotlib.axes.Axes
):
    plot_ewert_phenology_data_box_plot(
        species_config.key_lengths_flag_leaf_td.tl_em,
        species_config.key_lengths_flag_leaf_td.tl_ma,
        species_config.key_lengths_flag_leaf_td.tl_ep,
        species_config.key_lengths_flag_leaf_td.tl_se,
        species_config.key_dates_td.sowing,
        species_config.key_dates_td.harvest,
        species_config.key_dates_td.Astart,
        ax=ax,
        box_y_start=17,
    )


def plot_phenology_from_config(
    species_config: SpeciesConfig,
    model_config: ModelConfig,
    nP: int,
    output_location: Optional[Path] = None,
    external_data=None,
    day_count=365,
    plot_f_phen: bool = True,
    plot_lengths: bool = True,
    plot_carbon: bool = True,
    plot_growing: bool = True,
    plot_td: bool = True,
    plot_dd: bool = False,
    logger=print,
):
    """Output plots related to the phenology from input config.

    Parameters
    ----------
    species_config : SpeciesConfig
        Species config (unprocessed)
    model_config : ModelConfig
        Model config (unprocessed)
    nP : int
        Number of leaf populations
    output_location : Path, optional
        location to save plots, by default None
    external_data : [type], optional
        external data to use for config processing, by default None
    day_count : int, optional
        Total Number of days, by default 365
    plot_f_phen : bool, optional
        If true plot f_phen data, by default True
    plot_lengths : bool, optional
        If true plot the key phenology lengths, by default True
    plot_carbon : bool, optional
        If true plot the carbon plots, by default True
    plot_growing : bool, optional
        If true plot the leaf population plots, by default True
    plot_dd: bool, optional
        If true plots additional plots with day count on x axis
    plot_dd: bool, optional
        If true plots additional plots with day count on x axis

    """
    nrows = plot_f_phen + plot_lengths + plot_carbon + plot_growing
    ncols = (plot_dd and 1 or 0) + (plot_td and 1 or 0)
    assert nrows > 0, (
        "At least one plot type must be selected. Set plot_dd or plot_td to True."
    )
    external_data = {} if external_data is None else external_data
    # Get plot indexes
    f_phen_ax_i = 0 if plot_f_phen else -1
    lengths_ax_i = 0 + plot_f_phen if plot_lengths else -1
    carbon_ax_i = 0 + plot_f_phen + plot_lengths if plot_carbon else -1
    growing_ax_i = 0 + plot_f_phen + plot_lengths + plot_carbon if plot_growing else -1

    td_data = None
    dd_data = None
    # Generate mock data if not provided #
    if plot_td:
        T_b, T_o, T_m = [0, 20, 50]

        td_data = (
            external_data.get("td", None)
            if "td" in external_data
            else generate_example_td_data(day_count, T_b)
        )
        external_data["td"] = td_data

    if plot_dd:
        dd_data = (
            external_data.get("dd", None)
            if "dd" in external_data
            else generate_days_data(day_count)
        )
        external_data["dd"] = dd_data

    processed_model_config, processed_species_config = process_phenology_config(
        model_config,
        species_config,
        external_data,
        calculate_key_dates=False,
        td_base_temperature=0,
        nP=nP,
        logger=logger,
    )

    dvi = (
        get_dvi_range_from_species_config(processed_species_config, td_data)
        if model_config.dvi_method != DVIMethods.DISABLED
        else None
    )

    # ===== PLOTS ===== #
    fig, _axss = plt.subplots(
        ncols=ncols, nrows=nrows, figsize=(10 * 2, nrows * 5), dpi=120
    )
    # Handle axss dimensions vary with ncols and nrows. Check this is correct
    axss: list[list[matplotlib.axes.Axes]] = (
        _axss
        if ncols > 1 and nrows > 1
        else [_axss]
        if ncols > 1
        else [[a] for a in _axss]
        if nrows > 1
        else [[_axss]]
    )

    td_col_index = 0
    dd_col_index = 1 if plot_td and plot_dd else 0
    axs_td = plot_td and [
        ax[td_col_index] for ax in axss
    ]  #  if nrows > 1 else [axss[0]]
    axs_dd = plot_dd and [
        ax[dd_col_index] for ax in axss
    ]  #  if nrows > 1 else [axss[1]]

    if plot_td and axs_td:
        f_phen_ax_td = axs_td[f_phen_ax_i] if plot_f_phen else None
        lengths_ax_td = axs_td[lengths_ax_i] if plot_lengths else None
        carbon_ax_td = axs_td[carbon_ax_i] if plot_carbon else None
        growing_ax_td = axs_td[growing_ax_i] if plot_growing else None

    if plot_dd and axs_dd:
        f_phen_ax_dd = axs_dd[f_phen_ax_i] if plot_f_phen else None
        lengths_ax_dd = axs_dd[lengths_ax_i] if plot_lengths else None
        # carbon_ax_dd = axs_dd[carbon_ax_i] if plot_carbon else None
        # growing_ax_dd = axs_dd[growing_ax_i] if plot_growing else None

    phenology_plot_factory = PhenologyPlotClass(
        processed_species_config, td_data, dd_data
    )
    # Leaf f_phen plots
    if processed_species_config.leaf_f_phen_method not in [
        LeafFPhenMethods.DISABLED,
        LeafFPhenMethods.F_PHEN,
    ]:
        phenology_plot_factory.leaf_f_phen.plot_td(
            ax=f_phen_ax_td
        ) if plot_td and plot_f_phen and f_phen_ax_td else None
        phenology_plot_factory.leaf_f_phen.plot_dd(
            ax=f_phen_ax_dd
        ) if plot_dd and plot_f_phen and f_phen_ax_dd else None
        phenology_plot_factory.leaf_f_phen_box.plot_td(
            ax=lengths_ax_td
        ) if plot_td and plot_lengths and lengths_ax_td else None
        phenology_plot_factory.leaf_f_phen_box.plot_dd(
            ax=lengths_ax_dd
        ) if plot_dd and plot_lengths and lengths_ax_dd else None
    # f_phen plots
    phenology_plot_factory.f_phen.plot_td(
        ax=f_phen_ax_td
    ) if plot_td and f_phen_ax_td else None
    phenology_plot_factory.f_phen.plot_dd(
        ax=f_phen_ax_dd
    ) if plot_dd and f_phen_ax_dd else None
    phenology_plot_factory.f_phen_box.plot_td(
        ax=lengths_ax_td
    ) if plot_td and lengths_ax_td else None
    phenology_plot_factory.f_phen_box.plot_dd(
        ax=lengths_ax_dd
    ) if plot_dd and lengths_ax_dd else None

    # Ewert phenology plots
    phenology_plot_factory.ewert.plot_td(
        ax=f_phen_ax_td
    ) if plot_td and plot_f_phen and f_phen_ax_td else None
    # phenology_plot_factory.ewert.plot_dd(ax=f_phen_ax_dd) if plot_dd and plot_f_phen and f_phen_ax_dd else None
    phenology_plot_factory.ewert_box.plot_td(
        ax=lengths_ax_td
    ) if plot_td and plot_lengths and lengths_ax_td else None
    # phenology_plot_factory.ewert_box.plot_dd(ax=lengths_ax_dd) if plot_dd and plot_lengths and lengths_ax_dd else None

    if plot_td:
        assert td_data is not None
        # We only plot carbon and growing for TD data as DVI is not calculated for DD data
        plot_carbon_fractions_from_species_config(
            processed_species_config, td_data, dvi, nP, ax=carbon_ax_td
        ) if plot_carbon and carbon_ax_td and dvi is not None else None

        plot_growing_fractions_from_config(
            processed_species_config, td_data, dvi, nP, ax=growing_ax_td
        ) if plot_growing and growing_ax_td and dvi is not None else None
        plot_growing_populations_from_config(
            processed_species_config, nP, td_data, ax=growing_ax_td
        ) if plot_growing and growing_ax_td else None

        f_phen_ax_td.plot(
            td_data, dvi, label="DVI"
        ) if plot_f_phen and f_phen_ax_td and dvi is not None else None
        growing_ax_td.plot(
            td_data, dvi, label="DVI"
        ) if plot_f_phen and growing_ax_td and dvi is not None else None

    if plot_td and axs_td:
        t_sgs = processed_species_config.key_dates_td.sowing
        t_egs = processed_species_config.key_dates_td.harvest
        assert t_sgs is not None and t_egs is not None, (
            "SGS or EGS not defined in species config"
        )
        # Setup axes
        for ax in axs_td:
            ax.set_xlim((t_sgs - 100, t_egs + 300))
            ax.legend()

    if plot_dd and axs_dd:
        sgs = processed_species_config.key_dates.sowing
        egs = processed_species_config.key_dates.harvest
        assert sgs is not None and egs is not None, (
            "SGS or EGS not defined in species config"
        )

        # Setup axes
        for ax in axs_dd:
            end_lim = egs + 5 if egs > sgs else egs + 365 + 5
            ax.set_xlim((sgs - 5, end_lim))
            ax.legend()

    if plot_dd and plot_td:
        f_phen_ax_dd.set_title("Phenology") if f_phen_ax_dd else None
        f_phen_ax_dd.set_xlabel("days") if f_phen_ax_dd else None
        f_phen_ax_td.set_title("Phenology") if f_phen_ax_td else None
        lengths_ax_td.set_title("Phenology lengths") if lengths_ax_td else None
        f_phen_ax_td.set_xlabel("thermal time") if f_phen_ax_td else None
        f_phen_ax_td.set_ylabel("fraction") if f_phen_ax_td else None
    elif plot_td:
        if f_phen_ax_td:
            f_phen_ax_td.set_ylabel("fraction")
        if carbon_ax_td:
            carbon_ax_td.set_ylabel("fraction")
    elif plot_dd:
        if f_phen_ax_dd:
            f_phen_ax_dd.set_ylabel("fraction")

    axss[0][0].set_title("phenology")
    axss[-1][0].set_xlabel("days" if plot_dd else "thermal time")

    if output_location:
        plt.savefig(output_location)
    return fig, axss


class PhenologyPlotClass(_PhenologyPlotClass):
    leaf_f_phen: PlotLeafFPhenDataFromConfig
    leaf_f_phen_box: PlotLeafFPhenDataBoxPlotFromConfig
    f_phen: PlotFPhenDataFromConfig
    f_phen_box: PlotFPhenDataBoxPlotFromConfig
    ewert: PlotEwertPhenologyDataFromConfig
    ewert_box: PlotEwertPhenologyDataBoxPlotFromConfig

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        plot_args = [
            self.speciesConfig,
            self.td_data,
            self.dd_data,
        ]
        self.leaf_f_phen = PlotLeafFPhenDataFromConfig(*plot_args)
        self.leaf_f_phen_box = PlotLeafFPhenDataBoxPlotFromConfig(*plot_args)
        self.f_phen = PlotFPhenDataFromConfig(*plot_args)
        self.f_phen_box = PlotFPhenDataBoxPlotFromConfig(*plot_args)
        self.ewert = PlotEwertPhenologyDataFromConfig(*plot_args)
        self.ewert_box = PlotEwertPhenologyDataBoxPlotFromConfig(*plot_args)

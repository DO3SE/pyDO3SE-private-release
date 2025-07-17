from itertools import count
from pathlib import Path
from typing import List
import numpy as np
import matplotlib.patches as patches

from do3se_phenology.carbon_allocation import calc_partition_coefficients
from do3se_phenology.switchboard import process_phenology_config
from do3se_phenology.canopy_structure import get_growing_populations_range_from_config
from do3se_phenology.phyllochron_dvi import get_dvi_range_from_species_config
from do3se_phenology.config import ModelConfig, SpeciesConfig
from do3se_phenology.utils import generate_days_data, generate_example_td_data

try:
    from matplotlib import pyplot as plt
    from matplotlib import colormaps as mpcm
except:
    plt = lambda *args, **kwargs: Exception("matplotlib required to draw graphs!")


def plot_f_phen_box_data(
    t_sgs,
    t_egs,
    t_f_phen_a,
    t_f_phen_b,
    t_f_phen_c,
    t_f_phen_d,
    ax=plt,
    box_y_start=0,
    **kwargs,
):
    TEXT_HEIGHT = 1.1
    TEXT_ANGLE = 90
    COLOR = "blue"
    ax.axvline(t_sgs, linestyle='dotted', color=COLOR)

    y = count(box_y_start, 2)
    text_offset = (t_egs - t_sgs) / 100
    draw_range(ax, t_sgs, length=t_f_phen_a,
               label="t_f_phen_a", COLOR=COLOR, text_offset=text_offset, y=next(y))
    draw_range(ax, t_sgs, length=t_f_phen_b,
               label="t_f_phen_b", COLOR=COLOR, text_offset=text_offset, y=next(y))
    draw_range(ax, t_sgs, length=t_f_phen_c,
               label="t_f_phen_c", COLOR=COLOR, text_offset=text_offset, y=next(y))
    # draw_range(ax, t_sgs, length=t_f_phen_d,
    #            label="t_f_phen_d", COLOR=COLOR, text_offset=text_offset, y=next(y))
    ax.scatter([], [], label="f_phen")
    return ax


def plot_f_phen_data(
    accumulated_temperatures,
    f_phen_data,
    t_sgs,
    t_egs,
    t_f_phen_a,
    t_f_phen_b,
    t_f_phen_c,
    t_f_phen_d,
    ax=plt,
    **kwargs,
):
    TEXT_HEIGHT = 2.5
    TEXT_ANGLE = 90
    COLOR = "blue"

    ax.plot(accumulated_temperatures, f_phen_data, label="f_phen_data", color=COLOR, **kwargs)
    ax.axvline(t_sgs, linestyle='dotted', color=COLOR)
    ax.axvline(t_sgs + t_f_phen_a, linestyle='dotted', color=COLOR)
    ax.axvline(t_sgs + t_f_phen_b, linestyle='dotted', color=COLOR)
    ax.axvline(t_sgs + t_f_phen_c, linestyle='dotted', color=COLOR)
    ax.axvline(t_sgs + t_f_phen_d, linestyle='dotted', color=COLOR)
    ax.text(t_sgs, TEXT_HEIGHT, 'SGS', rotation=TEXT_ANGLE, color=COLOR)
    ax.text(t_egs, TEXT_HEIGHT, 'EGS', rotation=TEXT_ANGLE, color=COLOR)
    ax.text(t_sgs + t_f_phen_a, TEXT_HEIGHT, 'Plant Emergence', rotation=TEXT_ANGLE, color=COLOR)
    ax.text(t_sgs + t_f_phen_b, TEXT_HEIGHT, 'f_phen_b', rotation=TEXT_ANGLE, color=COLOR)

    return ax


def draw_range(ax, start, end=None, label=None, length=None, COLOR=None, y=0, end_height=0.2, text_offset=1):
    end = end or start + length
    ax.axvline(start, linestyle='dotted', color=COLOR)
    ax.axvline(end, linestyle='dotted', color=COLOR)
    ax.plot([start, start], [y - end_height, y + end_height], color=COLOR, linestyle='dotted')
    ax.plot([end, end], [y - end_height, y + end_height], color=COLOR, linestyle='dotted')
    ax.text(end + text_offset, y, label, color=COLOR)
    ax.plot([start, end], [y, y], color=COLOR)


def plot_leaf_f_phen_box_plot(
    t_leaf_f_phen_a,
    t_leaf_f_phen_b,
    t_leaf_f_phen_e,
    t_leaf_f_phen_g,
    t_leaf_f_phen_h,
    t_leaf_f_phen_i,
    t_astart,
    t_sgs,
    t_egs,
    ax=None,
    box_y_start=0,
    **kwargs,
):
    TEXT_HEIGHT = 1.1
    TEXT_ANGLE = 90
    COLOR = "red"
    offset = t_sgs + t_astart
    if ax:
        ax.axvline(t_sgs, linestyle='dotted', color=COLOR)

        y = count(box_y_start, 2)
        text_offset = (t_egs - t_sgs) / 100
        draw_range(ax, t_sgs, length=t_astart,
                   label="sgs -> astart", COLOR=COLOR, text_offset=text_offset, y=next(y))
        draw_range(ax, t_sgs + t_astart, length=t_leaf_f_phen_e,
                   label="leaf_f_phen_e", COLOR=COLOR, text_offset=text_offset, y=next(y))
        draw_range(ax, t_sgs + t_astart + t_leaf_f_phen_e, length=t_leaf_f_phen_g,
                   label="leaf_f_phen_g", COLOR=COLOR, text_offset=text_offset, y=next(y))
        draw_range(ax, t_sgs + t_astart + t_leaf_f_phen_e, length=t_leaf_f_phen_h,
                   label="leaf_f_phen_h", COLOR=COLOR, text_offset=text_offset, y=next(y))
        draw_range(ax, t_sgs + t_astart + t_leaf_f_phen_e, length=t_leaf_f_phen_i,
                   label="leaf_f_phen_i", COLOR=COLOR, text_offset=text_offset, y=next(y))
    ax.scatter([], [], label="leaf_f_phen")
    return ax


def plot_leaf_f_phen_data(
    accumulated_temperatures,
    leaf_f_phen_data,
    t_sgs,
    ax=plt,
    **kwargs,
):
    COLOR = "red"
    ax.plot(accumulated_temperatures, leaf_f_phen_data,
            label="leaf_f_phen_data", color=COLOR, **kwargs)
    ax.axvline(t_sgs, linestyle='dotted', color=COLOR)
    return ax


def plot_ewert_phenology_data_box_plot(
    t_lem,
    t_lma,
    t_lep,
    t_lse,
    t_sgs,
    t_egs,
    t_Astart,
    ax=plt,
    box_y_start=0,
    **kwargs,
):
    TEXT_HEIGHT = 1.1
    TEXT_ANGLE = 90
    COLOR = "green"

    ax.axvline(t_sgs, linestyle='dotted', color=COLOR)

    y = count(box_y_start, 2)
    text_offset = (t_egs - t_sgs) / 100
    draw_range(ax, t_sgs + t_Astart - t_lem, length=t_lem,
               label="t_lem", COLOR=COLOR, text_offset=text_offset, y=next(y))
    draw_range(ax, t_sgs + t_Astart, length=t_lma,
               label="t_lma", COLOR=COLOR, text_offset=text_offset, y=next(y))
    draw_range(ax, t_sgs + t_Astart, length=t_lep,
               label="t_lep", COLOR=COLOR, text_offset=text_offset, y=next(y))
    draw_range(ax, t_sgs + t_Astart + t_lep, length=t_lse,
               label="t_lse", COLOR=COLOR, text_offset=text_offset, y=next(y))
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
    ax=plt,
    **kwargs,
):
    TEXT_HEIGHT = 1.1
    TEXT_ANGLE = 90
    COLOR = "green"

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
    ax=plt,
    **kwargs,
):
    TEXT_HEIGHT = 1.1
    TEXT_ANGLE = 90
    COLOR = "green"
    ax.plot(accumulated_temperatures, f_LA_data,
            label="f_LA", color=COLOR, **kwargs)

    ax.plot(accumulated_temperatures, f_LS_data, '--',
            label="f_LS", color=COLOR, **kwargs)

    ax.axvline(t_sgs, linestyle='dotted', color=COLOR)
    ax.axvline(t_sgs + t_Astart, linestyle='dotted', color=COLOR)
    ax.axvline(t_sgs + t_Astart - t_lem, linestyle='dotted', color=COLOR)

    ax.axvline(t_sgs + t_Astart + t_lma, linestyle='dotted', color=COLOR)
    ax.axvline(t_sgs + t_Astart + t_lep, linestyle='dotted', color=COLOR)


    ax.text(t_sgs, 2.5, "SGS", rotation=TEXT_ANGLE)
    ax.text(t_sgs + t_Astart - t_lem, 2.5, "Flag Emergence", rotation=TEXT_ANGLE)
    ax.text(t_sgs + t_Astart + t_lep, 2.5, "Flag Onset senes", rotation=TEXT_ANGLE)
    ax.text(t_sgs + t_Astart, 2.5, "AStart", rotation=TEXT_ANGLE)
    ax.text(t_egs, 2.5, "EGS", rotation=TEXT_ANGLE)

    return ax


def plot_growing_fractions_from_config(
    species_config, td, dvi, nP,
    a_root=18.5,
    a_leaf=16.0,
    a_stem=18.0,
    b_root=-20.0,
    b_leaf=-15.0,
    b_stem=-18.5, ax=plt,
):
    ax.set_title("Growing populations")
    growing_populations_all, emerged_leaf_populations_count = get_growing_populations_range_from_config(
        species_config, nP, td)

    p_root, p_leaf, p_stem, p_harv = list(zip(*[calc_partition_coefficients(d, a_root,
                                                                            a_leaf,
                                                                            a_stem,
                                                                            b_root,
                                                                            b_leaf,
                                                                            b_stem,) for d in dvi]))
    leaf_carbon_fracs = [[
        sum(g[i + 1:]) * p / sum(g) if ff == 0 or sum(g) == 0 else p *
        (ff + sum(g[i + 1:])) / sum(g)  if i < len(g) else (p * ff) / sum(g)
        for i, ff in enumerate(g)] for g, p in zip(growing_populations_all, p_leaf)]

    colors = [i / nP for i in range(nP)]
    # cmap=plt.cm.bgrcmyk
    cmap = mpcm.get_cmap('Spectral')
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
    ax=plt,
):
    p_root, p_leaf, p_stem, p_harv = list(zip(*[calc_partition_coefficients(d, a_root,
                                                                            a_leaf,
                                                                            a_stem,
                                                                            b_root,
                                                                            b_leaf,
                                                                            b_stem,) for d in dvi]))

    ax.set_title("Carbon distribution fractions")
    ax.plot(td, p_root, label="p_root")
    ax.plot(td, p_leaf, label="p_leaf")
    ax.plot(td, p_stem, label="p_stem")
    ax.plot(td, p_harv, label="p_harv")
    return ax


def plot_leaf_f_phen_data_from_config(species_config, td_data, dd_data, ax=plt, axb=None):
    leaf_f_phen_data_td = np.interp(
        td_data, [i[0] for i in species_config.leaf_fphen_intervals],
        [i[1]for i in species_config.leaf_fphen_intervals])
    plot_leaf_f_phen_data(
        td_data,
        leaf_f_phen_data_td,
        species_config.key_dates_td.sowing,
        ax=ax,
    )
    if axb:
        leaf_f_phen_data_dd = np.concatenate(
        [np.zeros(species_config.key_dates.sowing * 24), leaf_f_phen_data_td])[0:len(leaf_f_phen_data_td)]

        plot_leaf_f_phen_data(
            dd_data,
            leaf_f_phen_data_dd,
            species_config.key_dates.sowing,
            ax=axb,
        )


def plot_leaf_f_phen_data_box_plot_from_config(species_config, ax=plt):
    plot_leaf_f_phen_box_plot(
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
    )


def plot_f_phen_data_from_config(species_config, td_data, dd_data, ax=plt, axb=None):
    f_phen_data = np.interp(td_data, [i[0] for i in species_config.fphen_intervals],
                            [i[1]for i in species_config.fphen_intervals])

    plot_f_phen_data(
        td_data,
        f_phen_data,
        species_config.key_dates_td.sowing,
        species_config.key_dates_td.harvest,
        species_config.key_lengths_td.sowing_to_emerge,
        species_config.key_lengths_td.sowing_to_f_phen_b,
        species_config.key_lengths_td.sowing_to_f_phen_c,
        species_config.key_lengths_td.sowing_to_end,
        ax=ax,
    )
    if axb:
        f_phen_data_dd = np.concatenate(
        [np.zeros(species_config.key_dates.sowing * 24), f_phen_data])[0:len(f_phen_data)]

        plot_f_phen_data(
            dd_data,
            f_phen_data_dd,
            species_config.key_dates.sowing,
            species_config.key_dates.harvest,
            species_config.key_lengths.sowing_to_emerge,
            species_config.key_lengths.sowing_to_f_phen_b,
            species_config.key_lengths.sowing_to_f_phen_c,
            species_config.key_lengths.sowing_to_end,
            ax=axb,
        )


def plot_f_phen_data_box_plt_from_config(species_config, ax=plt):
    plot_f_phen_box_data(
        species_config.key_dates_td.sowing,
        species_config.key_dates_td.harvest,
        species_config.key_lengths_td.sowing_to_emerge,
        species_config.key_lengths_td.sowing_to_f_phen_b,
        species_config.key_lengths_td.sowing_to_f_phen_c,
        species_config.key_lengths_td.sowing_to_end,
        ax=ax,
    )


def plot_growing_populations_from_config(species_config, nP, td, row_height=0.15, base_offset=1, ax=plt):
    growing_populations_all, emerged_leaf_populations_count = get_growing_populations_range_from_config(
        species_config, nP, td)

    colors = [i / nP for i in range(nP)]
    # cmap=plt.cm.bgrcmyk
    cmap = mpcm.get_cmap('Spectral')
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
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # ax.legend(rects, labels)


def plot_ewert_phenology_data_from_config(species_config, td, ax=plt):
    td_dd_list = [t - species_config.key_dates_td.sowing for t in td]
    t_astart = species_config.key_lengths_td.sowing_to_astart
    t_lma = species_config.key_lengths_flag_leaf_td.tl_ma
    t_lep = species_config.key_lengths_flag_leaf_td.tl_ep
    t_lem = species_config.key_lengths_flag_leaf_td.tl_em
    f_LA = [max(0, min(1, 1 - (td_dd - t_astart) / (t_lma))) for td_dd in td_dd_list]
    growing = [max(0, min(1, (td_dd - (t_astart - t_lem)) / (t_lem))) for td_dd in td_dd_list]
    fO3_l = 1
    f_LS = [max(0, min(1, 1 - ((td_dd - t_astart - t_lep) / (t_lma / fO3_l - t_lep))))
            for td_dd in td_dd_list]

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


def plot_ewert_phenology_data_box_plot_from_config(species_config, ax=plt):

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
    output_location: Path = None,
    external_data={},
    day_count=365,
    plot_f_phen: bool = True,
    plot_lengths: bool = True,
    plot_carbon: bool = True,
    plot_growing: bool = True,
    plot_dd: bool = False,
    logger = print,
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

    """
    nrows = plot_f_phen + plot_lengths + plot_carbon + plot_growing
    ncols = 2 if plot_dd else 1
    f_phen_ax_i = 0 if plot_f_phen else None
    lengths_ax_i = 0 + plot_f_phen if plot_f_phen else -1
    carbon_ax_i = 0 + plot_f_phen + plot_lengths if plot_carbon else -1
    growing_ax_i = 0 + plot_f_phen + plot_lengths + plot_carbon if plot_growing else -1

    # Generate missing data #
    T_b, T_o, T_m = [0, 20, 50]

    td_data = external_data.get(
        'td', None) if 'td' in external_data else generate_example_td_data(day_count, T_b)
    dd_data = external_data.get(
        'dd', None) if 'dd' in external_data else generate_days_data(day_count)

    external_data['dd'] = dd_data
    external_data['td'] = td_data

    processed_model_config, processed_species_config = process_phenology_config(
        model_config,
        species_config,
        external_data,
        calculate_key_dates=True,
        td_base_temperature=0,
        nP=nP,
        logger=logger
    )

    dvi = get_dvi_range_from_species_config(processed_species_config, td_data)

    # ===== PLOTS ===== #
    fig, axss = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10 * 2, nrows * 5), dpi=120)
    # Handle axss dimensions vary with ncols and nrows. Check this is correct
    axss = axss if ncols > 1 and nrows > 1 else \
        [axss] if ncols > 1 else \
        [[a] for a in axss] if nrows > 1 else \
            [[axss]]

    axs_td = [ax[0] for ax in axss]#  if nrows > 1 else [axss[0]]
    axs_dd = plot_dd and [ax[1] for ax in axss]#  if nrows > 1 else [axss[1]]

    f_phen_ax = axs_td[f_phen_ax_i] if plot_f_phen else None
    lengths_ax = axs_td[lengths_ax_i] if plot_f_phen else None
    carbon_ax = axs_td[carbon_ax_i] if plot_carbon else None
    growing_ax = axs_td[growing_ax_i] if plot_growing else None

    f_phen_ax_dd = axs_dd[f_phen_ax_i] if plot_dd and plot_f_phen else None
    lengths_ax_dd = axs_dd[lengths_ax_i] if plot_dd and plot_f_phen else None
    carbon_ax_dd = axs_dd[carbon_ax_i] if plot_dd and plot_carbon else None
    growing_ax_dd = axs_dd[growing_ax_i] if plot_dd and plot_growing else None

    plot_leaf_f_phen_data_from_config(processed_species_config, td_data, dd_data,
                                      ax=f_phen_ax, axb=f_phen_ax_dd)if plot_f_phen else None
    plot_leaf_f_phen_data_box_plot_from_config(
        processed_species_config, ax=lengths_ax) if plot_lengths else None

    plot_f_phen_data_from_config(processed_species_config, td_data, dd_data,
                                 ax=f_phen_ax, axb=f_phen_ax_dd) if plot_f_phen else None
    plot_f_phen_data_box_plt_from_config(
        processed_species_config, ax=lengths_ax) if plot_lengths else None

    plot_ewert_phenology_data_from_config(
        processed_species_config, td_data, ax=f_phen_ax) if plot_f_phen else None
    plot_ewert_phenology_data_box_plot_from_config(
        processed_species_config, ax=lengths_ax) if plot_lengths else None

    plot_carbon_fractions_from_species_config(
        processed_species_config, td_data, dvi, nP, ax=carbon_ax) if plot_carbon else None

    plot_growing_fractions_from_config(processed_species_config, td_data, dvi, nP,
                                       ax=growing_ax) if plot_growing else None

    plot_growing_populations_from_config(
        processed_species_config, nP, td_data, ax=growing_ax) if plot_growing else None

    f_phen_ax.plot(td_data, dvi, label="DVI") if plot_f_phen else None
    growing_ax.plot(td_data, dvi, label="DVI") if plot_f_phen else None

    t_sgs = processed_species_config.key_dates_td.sowing
    t_egs = processed_species_config.key_dates_td.harvest

    # Setup axes
    for ax in axs_td:
        ax.set_xlim((t_sgs - 100, t_egs + 300))
        ax.legend()

    f_phen_ax.set_title("Phenology") if plot_f_phen else None
    f_phen_ax_dd.set_title("Phenology days") if plot_dd and plot_f_phen  else None
    f_phen_ax_dd.set_xlabel("days") if plot_dd and plot_f_phen  else None
    lengths_ax.set_title("Phenology lengths") if plot_lengths else None
    f_phen_ax.set_xlabel("thermal time") if plot_f_phen else None
    f_phen_ax.set_ylabel("%") if plot_f_phen else None

    if output_location:
        plt.savefig(output_location)
    return fig, axss

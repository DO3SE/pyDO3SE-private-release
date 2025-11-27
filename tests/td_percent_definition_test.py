import pytest
from math import isclose
import numpy as np

from do3se_phenology.presets.wheat import Wheat
from do3se_phenology.td_percent_definition import (
    calculate_growing_season_from_leaf_f_phen_data,
    get_canopy_td_intervals,
    get_current_f_phen_from_t_sgs_t_egs,
    get_current_leaf_f_phen_from_t_sgs_t_egs,
    get_leaf_td_intervals,
)
from thermal_time.calcs import calc_thermal_time_range
from do3se_phenology.f_phen import tt_f_phen_simple_PLF_range, tt_leaf_f_phen_PLF_range
from do3se_phenology.plots import (
    plot_ewert_phenology_data,
    plot_ewert_phenology_data_box_plot,
    plot_f_phen_tt_box_data,
    plot_f_phen_td_data,
    plot_leaf_f_phen_td_box_plot,
    plot_leaf_f_phen_data,
)


def test_get_canopy_td_intervals():
    td = np.arange(0, 2000, 2000 / 365)
    dd = np.arange(365)

    t_sgs = 200
    t_egs = 1775 + 200

    (
        SGS,
        EGS,
        Astart,
        mid_anthesis,
        fphen_a,
        fphen_b,
        fphen_c,
        fphen_d,
        dd_emr,
        dd_veg,
        dd_rep,
        t_Astart,
        t_mid_anthesis,
        t_fphen_a,
        t_fphen_b,
        t_fphen_c,
        t_fphen_d,
        tt_emr,
        tt_veg,
        tt_rep,
    ) = get_canopy_td_intervals(
        td,
        dd,
        t_sgs=t_sgs,
        t_egs=t_egs,
    )

    assert fphen_a < fphen_b < fphen_c < fphen_d
    assert Astart < mid_anthesis

    assert mid_anthesis == fphen_c
    assert tt_emr + tt_veg + tt_rep == t_egs - t_sgs

    assert SGS == 37
    assert EGS == 361
    assert Astart == 212
    assert mid_anthesis == 238
    assert fphen_a == 53
    assert fphen_b == 102
    assert fphen_c == 238
    assert fphen_d == 361
    assert dd_emr == 53
    assert dd_veg == 160
    assert dd_rep == 222

    assert isclose(t_Astart, 958.5, abs_tol=1e-3)
    assert isclose(t_mid_anthesis, 1100.5, abs_tol=1e-3)
    assert isclose(t_fphen_a, 88.75, abs_tol=1e-3)
    assert isclose(t_fphen_b, 355.0, abs_tol=1e-3)
    assert isclose(t_fphen_c, 1100.5, abs_tol=1e-3)
    assert isclose(t_fphen_d, 1775.0, abs_tol=1e-3)


def test_get_leaf_td_intervals():
    td = np.arange(0, 2000, 2000 / 365)
    dd = np.arange(365)
    t_season_length = 1775
    (
        SGS,
        EGS,
        Astart,
        mid_anthesis,
        fphen_1_ets,
        fphen_3_ets,
        fphen_4_ets,
        fphen_5_ets,
        t_Astart,
        t_mid_anthesis,
        t_fphen_1_ets,
        t_fphen_3_ets,
        t_fphen_4_ets,
        t_fphen_5_ets,
        t_lem,
        t_lse,
        t_lma,
        t_lep,
        lem,
        lse,
        lma,
        lep,
    ) = get_leaf_td_intervals(
        td,
        dd,
        t_sgs=200,
        t_egs=t_season_length + 200,
        f_Astart=Wheat.f_Astart,
        f_mid_anthesis=Wheat.f_mid_anthesis,
        f_fphen_1_ets=Wheat.f_fphen_1_ets,
        f_fphen_3_ets=Wheat.f_fphen_3_ets,
        f_fphen_4_ets=Wheat.f_fphen_4_ets,
        f_fphen_5_ets=Wheat.f_fphen_5_ets,
        f_t_lem=Wheat.f_t_lem,
        f_t_lma=Wheat.f_t_lma,
        f_t_lep=Wheat.f_t_lep,
        f_t_lse=Wheat.f_t_lse,
    )
    season_length = EGS - SGS
    assert Astart + fphen_1_ets + fphen_5_ets == EGS - SGS
    assert isclose(t_Astart, t_season_length - t_lma)
    assert lma == season_length - Astart
    # assert lem + lma == season_length - Astart
    assert lep + lse == lma

    # TODO: We don't know what lem size should be

    assert SGS == 37
    assert EGS == 361
    assert Astart == 175
    assert mid_anthesis == 201
    assert fphen_1_ets == 26
    assert fphen_3_ets == 16
    assert fphen_4_ets == 71
    assert fphen_5_ets == 123

    assert lem == 45
    assert lse == 48
    assert lma == 149
    assert lep == 101

    assert isclose(t_Astart, 958.50, abs_tol=1e-3)
    assert isclose(t_mid_anthesis, 1100.5, abs_tol=1e-3)
    assert isclose(t_fphen_1_ets, 142.0, abs_tol=1e-3)
    assert isclose(t_fphen_3_ets, 88.75, abs_tol=1e-3)
    assert isclose(t_fphen_4_ets, 390.5, abs_tol=1e-3)
    assert isclose(t_fphen_5_ets, 674.5, abs_tol=1e-3)
    assert isclose(t_lem, 248.5, abs_tol=1e-3)
    assert isclose(t_lse, 266.25, abs_tol=1e-3)
    assert isclose(t_lma, 816.5, abs_tol=1e-3)
    assert isclose(t_lep, 550.25, abs_tol=1e-3)


@pytest.mark.parametrize(
    ["td", "t_sgs", "t_egs", "expected_f_phen"],
    [
        (1, 0, 100, 0.0),
        (10, 0, 100, 0.333),
        (40, 0, 100, 1),
        (80, 0, 100, 0.526),
        (100, 0, 100, 0.0),
        (101, 0, 100, 0.0),
    ],
)
def test_get_fphen_from_t_sgs_t_egs(td, t_sgs, t_egs, expected_f_phen):
    out = get_current_f_phen_from_t_sgs_t_egs(
        td,
        t_sgs,
        t_egs,
        f_phen_min=Wheat.f_phen_min,
        f_Astart=Wheat.f_Astart,
        f_mid_anthesis=Wheat.f_mid_anthesis,
        f_fphen_a=Wheat.f_fphen_a,
        f_fphen_b=Wheat.f_fphen_b,
        f_fphen_c=Wheat.f_fphen_c,
        f_fphen_d=Wheat.f_fphen_d,
    )

    assert isclose(out, expected_f_phen, abs_tol=1e-3)


@pytest.mark.parametrize(
    ["td", "t_sgs", "t_egs", "expected_f_phen"],
    [
        (1, 0, 100, 0.0),
        (10, 0, 100, 0.0),
        (40, 0, 100, 0),
        (80, 0, 100, 0.7705),
        (100, 0, 100, 0.0),
        (101, 0, 100, 0.0),
    ],
)
def test_get_leaf_fphen_from_t_sgs_t_egs(td, t_sgs, t_egs, expected_f_phen):
    out = get_current_leaf_f_phen_from_t_sgs_t_egs(
        td,
        t_sgs,
        t_egs,
        f_Astart=Wheat.f_Astart,
        f_mid_anthesis=Wheat.f_mid_anthesis,
        f_fphen_1_ets=Wheat.f_fphen_1_ets,
        f_fphen_3_ets=Wheat.f_fphen_3_ets,
        f_fphen_4_ets=Wheat.f_fphen_4_ets,
        f_fphen_5_ets=Wheat.f_fphen_5_ets,
        f_leaf_f_phen_a=Wheat.leaf_f_phen_a,
        f_leaf_f_phen_b=Wheat.leaf_f_phen_b,
    )

    assert isclose(out, expected_f_phen, abs_tol=1e-3)


def test_calculate_growing_season_from_leaf_f_phen_data():
    example_SGS_t = 100

    leaf_f_phen_data = np.interp(
        np.arange(0, 2200 + 20, 20 / 24),
        np.array([0, 1080, 1081, 1240, 1340, 1680, 2000]) - example_SGS_t,
        [0.0, 0.0, 1.0, 1.0, 1.0, 0.8, 0.1],
    ).tolist()

    dd = [j for j in range(111) for _ in range(24)]
    td = [j * 20 for j in range(111) for _ in range(24)]

    # %%
    SGS, Astart, Aend, mature_start, senes_start = calculate_growing_season_from_leaf_f_phen_data(
        leaf_f_phen_data,
        td,
    )
    assert isclose(SGS, -452.1739, abs_tol=1e-3)
    assert isclose(Astart, 980, abs_tol=1e-3)
    assert isclose(Aend, 2200, abs_tol=1e-3)
    assert isclose(mature_start, 1240, abs_tol=1e-3)
    assert isclose(senes_start, 1580, abs_tol=1e-3)


def test_plot_overlays():
    from matplotlib import pyplot as plt

    t_sgs = 100
    t_egs = 1000
    T_b, T_o, T_m = [0, 20, 50]
    day_count = 365
    dd_data = np.array([[d for i in range(24)] for d in range(day_count)]).reshape(day_count * 24)
    hrs_data = np.array([[i for i in range(24)] for d in range(day_count)]).reshape(day_count * 24)
    demo_temp_data = [24 - abs(hr - 12) for hr in hrs_data]
    tsc = demo_temp_data
    td = calc_thermal_time_range(tsc, t_b=T_b)

    fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(16, 10))
    (
        SGS,
        EGS,
        Astart,
        mid_anthesis,
        fphen_a,
        fphen_b,
        fphen_c,
        fphen_d,
        dd_emr,
        dd_veg,
        dd_rep,
        t_Astart,
        t_mid_anthesis,
        t_f_phen_a,
        t_f_phen_b,
        t_f_phen_c,
        t_f_phen_d,
        tt_emr,
        tt_veg,
        tt_rep,
    ) = get_canopy_td_intervals(
        td,
        dd_data,
        t_sgs=t_sgs,
        t_egs=t_egs,
    )

    f_phen_data = tt_f_phen_simple_PLF_range(
        td,
        t_f_phen_a,
        t_f_phen_b,
        t_f_phen_c,
        t_f_phen_d,
        f_phen_min=0.2,
        td_at_sgs=t_sgs,
    )

    (
        SGS,
        EGS,
        Astart,
        mid_anthesis,
        fphen_1_ets,
        fphen_3_ets,
        fphen_4_ets,
        fphen_5_ets,
        t_Astart,
        t_mid_anthesis,
        t_fphen_1_ets,
        t_fphen_3_ets,
        t_fphen_4_ets,
        t_fphen_5_ets,
        t_lem,
        t_lse,
        t_lma,
        t_lep,
        lem,
        lse,
        lma,
        lep,
    ) = get_leaf_td_intervals(
        td,
        dd_data,
        t_sgs=t_sgs,
        t_egs=t_egs,
    )

    t_leaf_f_phen_a = 0.3
    t_leaf_f_phen_b = 0.7
    t_leaf_f_phen_e = t_fphen_1_ets
    t_leaf_f_phen_g = t_fphen_3_ets
    t_leaf_f_phen_h = t_fphen_4_ets
    t_leaf_f_phen_i = t_fphen_5_ets
    t_astart = t_Astart
    td_at_sgs = t_sgs

    leaf_f_phen_data = tt_leaf_f_phen_PLF_range(
        td,
        t_leaf_f_phen_a,
        t_leaf_f_phen_b,
        t_leaf_f_phen_e,
        t_leaf_f_phen_g,
        t_leaf_f_phen_h,
        t_leaf_f_phen_i,
        t_astart,
        td_at_sgs,
    )

    # ===== PLOTS ===== #
    # Leaf f phen
    # Plot thermal time
    plot_leaf_f_phen_data(
        td,
        leaf_f_phen_data,
        t_egs,
        ax=axs[0][0],
    )

    plot_leaf_f_phen_td_box_plot(
        t_leaf_f_phen_a,
        t_leaf_f_phen_b,
        t_leaf_f_phen_e,
        t_leaf_f_phen_g,
        t_leaf_f_phen_h,
        t_leaf_f_phen_i,
        t_astart,
        t_sgs,
        t_egs,
        ax=axs[1][0],
    )

    # Plot days
    plot_leaf_f_phen_data(
        dd_data,
        leaf_f_phen_data,
        SGS,
        ax=axs[0][1],
    )

    plot_leaf_f_phen_td_box_plot(
        t_leaf_f_phen_a,
        t_leaf_f_phen_b,
        fphen_1_ets,
        fphen_3_ets,
        fphen_4_ets,
        fphen_5_ets,
        Astart,
        SGS,
        EGS,
        ax=axs[1][1],
    )

    # f phen
    plot_f_phen_td_data(
        td,
        f_phen_data,
        t_sgs,
        t_egs,
        t_f_phen_a,
        t_f_phen_b,
        t_f_phen_c,
        t_f_phen_d,
        ax=axs[0][0],
    )

    plot_f_phen_td_data(
        dd_data,
        f_phen_data,
        SGS,
        EGS,
        fphen_a - SGS,
        fphen_b - SGS,
        fphen_c - SGS,
        fphen_d - SGS,
        ax=axs[0][1],
    )

    # f phen
    plot_f_phen_tt_box_data(
        t_sgs,
        t_egs,
        t_f_phen_a,
        t_f_phen_b,
        t_f_phen_c,
        t_f_phen_d,
        ax=axs[1][0],
        box_y_start=10,
    )

    plot_f_phen_tt_box_data(
        SGS,
        EGS,
        fphen_a - SGS,
        fphen_b - SGS,
        fphen_c - SGS,
        fphen_d - SGS,
        ax=axs[1][1],
        box_y_start=10,
    )

    # Plot fla fls
    td_dd_list = [t - t_sgs for t in td]
    f_LA = [max(0, min(1, 1 - (td_dd - t_astart) / (t_lma))) for td_dd in td_dd_list]
    # f_LA = [max(0,min(1,1 - (td_dd - t_lem) / (t_lma))) for td_dd in td_dd_list]
    fO3_l = 1
    f_LS = [
        max(0, min(1, 1 - ((td_dd - t_astart - t_lep) / (t_lma / fO3_l - t_lep))))
        for td_dd in td_dd_list
    ]
    # f_LS = [max(0,min(1,1 - ((td_dd - t_lem - t_lep) / (t_lma / fO3_l - t_lep)))) for td_dd in td_dd_list]
    t_l = t_lma + t_lem

    plot_ewert_phenology_data(
        td,
        f_LA,
        f_LS,
        t_lem,
        t_lma,
        t_lep,
        t_lse,
        t_sgs,
        t_egs,
        t_Astart,
        ax=axs[0][0],
    )
    plot_ewert_phenology_data(
        dd_data,
        f_LA,
        f_LS,
        lem,
        lma,
        lep,
        lse,
        SGS,
        EGS,
        Astart,
        ax=axs[0][1],
    )

    plot_ewert_phenology_data_box_plot(
        t_lem,
        t_lma,
        t_lep,
        t_lse,
        t_sgs,
        t_egs,
        t_Astart,
        ax=axs[1][0],
        box_y_start=17,
    )
    plot_ewert_phenology_data_box_plot(
        lem,
        lma,
        lep,
        lse,
        SGS,
        EGS,
        Astart,
        ax=axs[1][1],
        box_y_start=17,
    )

    dvi_x = [
        td_at_sgs,
        td_at_sgs + tt_emr,
        td_at_sgs + tt_emr + tt_veg,
        td_at_sgs + tt_emr + tt_veg + tt_rep,
    ]
    dvi_y = [-1, 0, 1, 2]

    axs[0][0].plot(dvi_x, dvi_y, label="DVI", color="yellow")

    # Setup axes
    axs[0][0].set_xlim((t_sgs - 100, t_egs + 300))
    axs[1][0].set_xlim((t_sgs - 100, t_egs + 300))
    axs[0][0].set_xlabel("thermal time")
    axs[0][0].set_ylabel("leaf_f_phen")
    axs[0][0].legend()

    axs[0][1].set_xlim(((SGS - 10), (EGS + 10)))
    axs[1][1].set_xlim(((SGS - 10), (EGS + 10)))

    axs[0][1].set_xlabel("days")
    axs[0][1].set_ylabel("leaf_f_phen")
    axs[0][1].legend()
    # plt.show()
    plt.savefig("tests/outputs/td_percent_overlay.png")

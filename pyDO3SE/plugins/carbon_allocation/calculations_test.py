import numpy as np
from math import isclose
from collections import namedtuple
import pytest
from .calculations import (
    calc_root_fraction_from_carbon,
    get_plant_height_from_carbon,
    calc_LAI_from_DVI_and_carbon,
    adjust_c_pools_at_eol,
    calc_partition_coefficients,
    calc_carbon_pool_change,
    daily_carbon_allocation,
    calc_net_prod,
)
from pyDO3SE.util.test_utils import process_snapshot


class TestCalcNetProd:

    def test_output_values(self):
        out = calc_net_prod(
            canopy_An=20,
            c_root=0.2,
            c_stem=0.2,
            c_leaf=0.2,
            R_dc=0.32,
            r_g=0.25,
        ).NPP

        assert isclose(out, 14.52, abs_tol=1e-3)

    def test_should_be_0_when_stem_root_leaf_ratio_high(self):
        c_root = 3
        c_stem = 3
        c_leaf = 0.01

        out = calc_net_prod(
            canopy_An=20,
            c_root=c_root,
            c_stem=c_stem,
            c_leaf=c_leaf,
            R_dc=0.32,
            r_g=0.25,
        ).NPP

        assert out == 0

    def test_should_be_less_near_acan_in_when_stem_root_leaf_ratio_low(self):
        c_root = 0.1
        c_stem = 0.1
        c_leaf = 10
        out = calc_net_prod(
            canopy_An=20.0,
            c_root=c_root,
            c_stem=c_stem,
            c_leaf=c_leaf,
            R_dc=0.32,
            r_g=0.25,
        )

        assert out.NPP > 0
        assert isclose(out.NPP, 14.9952, abs_tol=1e-3)


class TestCalcCarbonPoolChange:

    def test_output_values(self):
        out = calc_carbon_pool_change(
            net_prod_acc=3.81,
            p_root=0.1,
            p_leaf=0.2,
            p_stem=0.3,
            p_harv=0.4,
            theta=0.4,
        )
        assert isclose(out.c_root_diff, 0.381, abs_tol=1e-3)
        assert isclose(out.c_leaf_diff, 0.762, abs_tol=1e-3)
        assert isclose(out.c_stem_diff, 0.6858, abs_tol=1e-3)
        assert isclose(out.c_harv_diff, 1.524, abs_tol=1e-3)
        assert isclose(out.c_resv_diff, 0.4572, abs_tol=1e-3)


class TestCalcPartitionCoefficients:

    def test_output_values(self, snapshot):
        out = [calc_partition_coefficients(
            DVI=DVI,
            a_root=18.5,
            a_stem=16.0,
            a_leaf=18.0,
            b_root=-20.0,
            b_stem=-15.0,
            b_leaf=-18.5,
        ) for DVI in [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]]

        snapshot.assert_match(process_snapshot([o.p_root for o in out]), 'p_root')
        snapshot.assert_match(process_snapshot([o.p_leaf for o in out]), 'p_leaf')
        snapshot.assert_match(process_snapshot([o.p_stem for o in out]), 'p_stem')
        snapshot.assert_match(process_snapshot([o.p_harv for o in out]), 'p_harv')

    @pytest.mark.parametrize('DVI', [-1, 0, 1, 2])
    def test_fractions_add_to_1(self, DVI):
        out = calc_partition_coefficients(
            DVI=DVI,
            a_root=18.5,
            a_stem=16.0,
            a_leaf=18.0,
            b_root=-20.0,
            b_stem=-15.0,
            b_leaf=-18.5,
        )
        assert isclose(sum(out), 1, abs_tol=1e-7)


class TestAdjustCPoolsAtEol:

    @pytest.mark.parametrize('DVI', [-1, 0, 1])
    @pytest.mark.parametrize('p_stem', [0.3, 0.5, 0.8])
    def test_no_change_at_early_stage(self, DVI, p_stem):
        c_leaf = 10.0
        c_harv = 15.0
        c_resv = 20.0
        out = adjust_c_pools_at_eol(
            DVI=DVI,
            c_leaf=c_leaf,
            p_stem=p_stem,
            c_harv=c_harv,
            c_resv=c_resv,
        )
        assert out.c_leaf == c_leaf
        assert out.c_harv == c_harv
        assert out.c_resv == c_resv

    @pytest.mark.parametrize('DVI', [1.55, 1.8, 2])
    @pytest.mark.parametrize('p_stem', [0.008, 0.0002, 0])
    def test_change_in_pools_at_end_stage(self, DVI, p_stem):
        c_leaf = 10.0
        c_harv = 15.0
        c_resv = 20.0
        out = adjust_c_pools_at_eol(
            DVI=DVI,
            c_leaf=c_leaf,
            p_stem=p_stem,
            c_harv=c_harv,
            c_resv=c_resv,
        )
        assert out.c_leaf < c_leaf
        assert out.c_harv > c_harv
        assert out.c_resv < c_resv

    @pytest.mark.parametrize('DVI', [0, 0.5, 1.3, 1.8, 2])
    @pytest.mark.parametrize('p_stem', [0, 0.0004, 0.001, 0.1, 0.8, 0.99])
    @pytest.mark.parametrize('c_leaf', [0.1, 1, 10, 100])
    @pytest.mark.parametrize('c_harv', [0.1, 1, 10, 100])
    @pytest.mark.parametrize('c_resv', [0.1, 1, 10, 100])
    def test_total_carbon_does_not_change(self, DVI, p_stem, c_leaf, c_harv, c_resv):
        out = adjust_c_pools_at_eol(
            DVI=DVI,
            c_leaf=c_leaf,
            p_stem=p_stem,
            c_harv=c_harv,
            c_resv=c_resv,
        )
        total_in = sum([c_leaf, c_harv, c_resv])
        total_out = sum([out.c_leaf, out.c_harv, out.c_resv])
        assert isclose(total_in, total_out, abs_tol=1e-6)


class TestCalcLAIFromDVIAndCarbon:

    @pytest.mark.parametrize('DVI', [0, 0.5, 1.3, 1.8, 2])
    @pytest.mark.parametrize('c_leaf', [0, 0.001, 0.01, 0.1, 1])
    def test_output_values(self, snapshot, DVI, c_leaf):
        out = calc_LAI_from_DVI_and_carbon(
            DVI=DVI,
            c_leaf=c_leaf,
            gamma=27.3,
            delta=-0.0507,
            f_c=0.5,
            emerged_leaf_count=2,
        )

        snapshot.assert_match(process_snapshot(out), f"DVI {DVI} c_leaf {c_leaf}")


class TestGetPlantHeightFromCarbon:

    @pytest.mark.parametrize('c_stem', [0, 0.001, 0.01, 0.1, 1])
    def test_output_values(self, c_stem, snapshot):
        out = get_plant_height_from_carbon(
            c_stem=c_stem,
            k=1.4,
            lambdav=0.4,
            f_c=0.5,
        )
        snapshot.assert_match(process_snapshot(out), f"c_stem {c_stem}")


class TestCalcRootFractionFromCarbon:

    @pytest.mark.parametrize('c_root', [0, 0.001, 0.01, 0.1, 1])
    def test_output_values(self, snapshot, c_root):
        out = calc_root_fraction_from_carbon(
            c_root=c_root,
            d_r=0.5,
            r_dir=0.0,
            z=0.5,
            f_c=0.5,
        )
        snapshot.assert_match(process_snapshot(out), f"c_root {c_root}")


class TestDailyCarbonAllocation:

    def test_single_day(self):

        c_root, c_stem, c_leaf, c_harv, c_resv, c_lbrn, *fractions = daily_carbon_allocation(
            net_prod_acc=20,
            DVI=1.3,
            c_root=0.2,
            c_stem=0.2,
            c_leaf=0.2,
            c_harv=0.2,
            c_resv=0.2,
            c_lbrn=0.2,
            a_root=18.5,
            a_stem=16.0,
            a_leaf=18.0,
            b_root=-20.0,
            b_stem=-15.0,
            b_leaf=-18.5,
            theta=0.4,
        )

        assert isclose(c_root, 0.2107, abs_tol=1e-3)
        assert isclose(c_stem, 0.5507, abs_tol=1e-3)
        assert isclose(c_leaf, 0.2456, abs_tol=1e-3)
        assert isclose(c_harv, 19.559, abs_tol=1e-3)
        assert isclose(c_resv, 0.4338, abs_tol=1e-3)

    def test_growing_season(self, snapshot):
        cparams = {
            "a_root": 18.5,  # Param
            "a_stem": 16.0,  # Param
            "a_leaf": 18.0,  # Param
            "b_root": -20.0,  # Param
            "b_stem": -15.0,  # Param
            "b_leaf": -18.5,  # Param
            "gamma": 27.3,
            "delta": -0.0507,
            "r_g": 0.1,
            "R_dc": 0.32,
            "theta": 0.4,
            "k": 1.4,
            "lambdav": 0.4,
        }
        # State
        ModelState = namedtuple('ModelState', [
            "c_root",
            "c_stem",
            "c_leaf",
            "c_harv",
            "c_resv",
            "c_lbrn",
            "net_prod_acc",
            "lai",
            "plant_height",
        ])
        # Parameters

        # Initial state
        model_state = ModelState(
            c_root=0,
            c_stem=0,
            c_leaf=0.01,
            c_harv=0,
            c_resv=0,
            c_lbrn=0,
            net_prod_acc=0,
            lai=0,
            plant_height=0,
        )
        # Setup demo data
        hrs = [hr for dd in range(365) for hr in range(24)]
        dd_hourly = [dd for dd in range(365) for hr in range(24)]
        dd_daily = [dd for dd in range(365)]
        x = [0, 50, 100, 150, 200, 250, 270, 280, 290, 300, 320, 350]
        y = [0, 0.5, 1.5, 3, 8, 20, 38, 37, 23, 7, 1, 0]

        UMOL_TO_G = 0.0432
        a = 0.5
        hourly_ratio = np.concatenate(
            (np.arange(0, 1, 1 / 12), np.flip(np.arange(0, 1, 1 / 12))), axis=None)
        # convert from umol to grams and seconds to hours
        anet_daily_max = np.interp(dd_daily, x, y) * UMOL_TO_G * a
        anet_hourly = [a * hm for a in anet_daily_max for hm in hourly_ratio]

        # Generate DVI data
        a = -50
        b = 0.8
        x = b * np.array([a, 100, 130, 140, 190, 210, 230, 270, 340, 341, 342, 350]) - a
        y = np.array([-1, -1, -0.99, 0.0, 0.4, 0.5, 0.7, 1.1, 2, 2, 2, 2])
        dvi_hourly = np.interp(dd_hourly, x, y)

        # RUN
        for dd in range(365):
            net_prod_acc = 0
            for hr in range(24):
                row_index = (dd * 24) + hr
                An_canopy = anet_hourly[row_index]
                net_prod_acc += calc_net_prod(
                    An_canopy,
                    model_state.c_root,
                    model_state.c_stem,
                    model_state.c_leaf,
                    cparams["r_g"],
                    cparams["R_dc"],
                ).NPP
            net_prod_acc = max(0, net_prod_acc)

            DVI = dvi_hourly[row_index]

            c_root, c_stem, c_leaf, c_harv, c_resv, c_lbrn, *fractions = daily_carbon_allocation(
                net_prod_acc=net_prod_acc,
                DVI=DVI,
                c_root=model_state.c_root,
                c_stem=model_state.c_stem,
                c_leaf=model_state.c_leaf,
                c_harv=model_state.c_harv,
                c_resv=model_state.c_resv,
                c_lbrn=model_state.c_lbrn,
                a_root=cparams["a_root"],
                a_leaf=cparams["a_leaf"],
                a_stem=cparams["a_stem"],
                b_root=cparams["b_root"],
                b_leaf=cparams["b_leaf"],
                b_stem=cparams["b_stem"],
                theta=cparams["theta"],
            )
            lai = calc_LAI_from_DVI_and_carbon(
                DVI,
                c_leaf,
                cparams["gamma"],
                cparams["delta"],
                emerged_leaf_count=3,
            )
            plant_height = get_plant_height_from_carbon(
                c_stem,
                cparams["k"],
                cparams["lambdav"],
            )
            model_state = ModelState(
                c_root=c_root,
                c_stem=c_stem,
                c_leaf=c_leaf,
                c_harv=c_harv,
                c_resv=c_resv,
                c_lbrn=c_lbrn,
                net_prod_acc=net_prod_acc,
                lai=lai,
                plant_height=plant_height,
            )

        snapshot.assert_match(process_snapshot(model_state), "growing season carbon allocation")

    def test_reducing_green_leaf_fraction(self):
        common_args = dict(
            net_prod_acc=20,
            DVI=1.8,
            c_root=0.2,
            c_stem=0.2,
            c_leaf=0.2,
            c_harv=0.2,
            c_resv=0.2,
            c_lbrn=0.2,
            a_root=18.5,
            a_stem=16.0,
            a_leaf=18.0,
            b_root=-20.0,
            b_stem=-15.0,
            b_leaf=-18.5,
            theta=0.4,
            plant_is_senescing=True,
        )
        out_default = daily_carbon_allocation(
            **common_args,
            f_green_leaf = 0.95, # Defaults
            f_brown_leaf = 0.85, # Defaults
        )

        assert out_default.c_harv > 0.2

        out_decrease_green_f = daily_carbon_allocation(
            **common_args,
            f_green_leaf=0.5,
        )

        assert out_decrease_green_f.c_harv > out_default.c_harv

        out_decrease_brown_f = daily_carbon_allocation(
            **common_args,
            f_brown_leaf = 0.55, # Defaults
        )

        assert out_decrease_brown_f.c_harv > out_default.c_harv

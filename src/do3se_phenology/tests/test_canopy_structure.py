"""helpers tests

example values taken from spanish wheat config
"""

import pytest
from math import isclose
import numpy as np
from numpy.testing import assert_array_equal

from do3se_phenology.switchboard import process_phenology_config
from do3se_phenology.config import ModelConfig, PhenologyMethods, SpeciesConfig, SpeciesPresets
from do3se_phenology.canopy_structure import (
    LAI_day_PLF,
    LAI_method_estimate_canopy_total,
    LAI_method_estimate_total, SAI_wheat_and_LAI,
    calc_SAI_wheat,
    calc_distribution_of_LAI_between_lcs,
    calc_leaf_pops_per_layer,
    distribute_canopy_lai_to_leaf_pops,
    distribute_lai_per_layers,
    get_growing_populations_range_from_config,
)


def test_LAI_day_PLF():
    lai_out = LAI_day_PLF(
        dd=140,
        SGS=118,
        EGS=210,
        LAI_1=21,
        LAI_2=21,
        LAI_a=0.0,
        LAI_b=3.5,
        LAI_c=3.5,
        LAI_d=0.0,
    )
    assert lai_out == 3.5


def test_SAI_wheat_and_LAI():
    """Test that SAI_wheat matches UI output."""
    sai_out = SAI_wheat_and_LAI(
        dd=0,
        SGS=118,
        EGS=210,
        LAI_1=21,
        LAI_2=21,
        LAI_a=3.0,
        LAI_b=3.0,
        LAI_c=3.0,
        LAI_d=3.0,
    )
    assert sai_out == 3.0


def test_SAI_wheat():
    sai_out = calc_SAI_wheat(
        dd=140,
        LAI=3.5,
        SGS=118,
        EGS=210,
        LAI_1=21,
    )
    assert sai_out == 5.0
    sai_out = calc_SAI_wheat(
        dd=10,
        LAI=3.5,
        SGS=10,
        EGS=210,
        LAI_1=21,
    )
    assert sai_out == 5.0


def test_calc_distribution_of_LAI_between_lcs():
    nL = 3
    nLC = 2
    lai_values = [[1, 2, 3], [4, 5, 6]]
    assert np.array(lai_values).shape == (nLC, nL)
    LC_dist = calc_distribution_of_LAI_between_lcs(lai_values, nL, nLC)
    assert LC_dist == [6 / 21, 15 / 21]
    assert len(LC_dist) == nLC

    nL = 1
    nLC = 1
    lai_values = [[3.16]]
    LC_dist = calc_distribution_of_LAI_between_lcs(lai_values, nL, nLC)
    assert LC_dist == [1.0]
    assert len(LC_dist) == nLC



def test_LAI_method_estimate_total():
    """Test that LAI_method_estimate_total returns lai per layer."""
    out = LAI_method_estimate_total(
        dd=80,
        nL=2,
        nLC=3,
        SGS=20,
        EGS=200,
        LAI_1=21,
        LAI_2=21,
        LAI_a=0.0,
        LAI_b=3.5,
        LAI_c=3.5,
        LAI_d=0.0,
        fLAI=[[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]],
    )
    assert isclose(out[0][0], 0.35, abs_tol=1e-3)
    assert isclose(out[0][1], 0.7, abs_tol=1e-3)
    assert isclose(out[0][2], 1.05, abs_tol=1e-3)


def test_LAI_method_estimate_canopy_total():
    """Test that LAI_method_estimate_canopy_total returns lai per layer."""
    out = LAI_method_estimate_canopy_total(
        dd=80,
        SGS=20,
        EGS=200,
        LAI_1=21,
        LAI_2=21,
        LAI_a=0.0,
        LAI_b=3.5,
        LAI_c=3.5,
        LAI_d=0.0,
    )
    assert isclose(out, 3.5, abs_tol=1e-3)


@pytest.mark.parametrize('total_lai', [0, 1, 4, 8])
@pytest.mark.parametrize('nL', [1, 2, 3])
@pytest.mark.parametrize('max_lai_per_layer', [0.1, 1.5, 5])
class TestDistributeLaiPerLayers():

    def test_should_split_lai_between_all_layers(self, total_lai, nL, max_lai_per_layer):
        if nL * max_lai_per_layer >= total_lai:
            out = distribute_lai_per_layers(total_lai, nL, max_lai_per_layer)
            assert len(out) == nL

    def test_sum_of_layer_lai_should_equal_total_lai(self, total_lai, nL, max_lai_per_layer):
        if nL * max_lai_per_layer >= total_lai:
            out = distribute_lai_per_layers(total_lai, nL, max_lai_per_layer)
            assert sum(out) == total_lai

    def test_no_layer_should_exceed_max_lai(self, total_lai, nL, max_lai_per_layer):
        if nL * max_lai_per_layer >= total_lai:
            out = distribute_lai_per_layers(total_lai, nL, max_lai_per_layer)
            assert all([i <= max_lai_per_layer for i in out])

    def test_should_fail_if_not_enough_layers(self, total_lai, nL, max_lai_per_layer):
        if nL * max_lai_per_layer < total_lai:
            with pytest.raises(ValueError) as e:
                distribute_lai_per_layers(total_lai, nL, max_lai_per_layer)
            assert "Not enough layers for LAI" in str(e.value)


# @pytest.mark.parametrize('f_total_emerged_leaves', [0, 0.5, 1])
# @pytest.mark.parametrize('nL', [1, 2, 3])
# @pytest.mark.parametrize('nP', [1, 2, 3])
# @pytest.mark.parametrize('max_leaf_lai', [0.1, 1.5, 5])
# class TestGetGrowingLeafPopulations:

#     def test_should_return_no_growing_populations_before_emerged(self, f_total_emerged_leaves, nP, nL, max_leaf_lai):
#         prev_leaf_lais = [0 for _ in range(nP)]
#         out = get_growing_populations(
#             total_emerged_leaves=0,
#             max_leaf_lai=max_leaf_lai,
#             prev_leaf_lais=prev_leaf_lais,
#             nP=nP,
#         )
#         assert out == [False for _ in range(nP)]

#     def test_should_return_all_growing_if_all_emerged_and_not_above_max_lai(self, f_total_emerged_leaves, nP, nL, max_leaf_lai):
#         prev_leaf_lais = [0 for _ in range(nP)]
#         out = get_growing_populations(
#             total_emerged_leaves=nP,
#             max_leaf_lai=max_leaf_lai,
#             prev_leaf_lais=prev_leaf_lais,
#             nP=nP,
#         )
#         assert out == [True for _ in range(nP)]

#     def test_should_return_none_growing_if_all_emerged_and_at_max_lai(self, f_total_emerged_leaves, nP, nL, max_leaf_lai):
#         prev_leaf_lais = [max_leaf_lai for _ in range(nP)]
#         total_emerged_leaves = int(f_total_emerged_leaves * nP)
#         out = get_growing_populations(
#             total_emerged_leaves=total_emerged_leaves,
#             max_leaf_lai=max_leaf_lai,
#             prev_leaf_lais=prev_leaf_lais,
#             nP=nP,
#         )
#         assert out == [False for _ in range(nP)]

#     def test_growing_count_should_equal_emerged_count_if_none_above_max_lai(self, f_total_emerged_leaves, nP, nL, max_leaf_lai):
#         prev_leaf_lais = [0 for _ in range(nP)]
#         total_emerged_leaves = int(f_total_emerged_leaves * nP)
#         out = get_growing_populations(
#             total_emerged_leaves=total_emerged_leaves,
#             max_leaf_lai=max_leaf_lai,
#             prev_leaf_lais=prev_leaf_lais,
#             nP=nP,
#         )
#         assert sum(out) == total_emerged_leaves


@pytest.mark.parametrize('growing_populations', ['none', 'some', 'all'])
@pytest.mark.parametrize('nL', [1, 2, 3])
@pytest.mark.parametrize('nP', [1, 2, 3])
class TestCalcLeafPopsPerlayer:

    def test_should_all_be_0_before_any_emerged(self, nL, nP, growing_populations):
        prev_leaf_pops_per_lai = np.zeros((nL, nP))
        layers_lai = [0 for _ in range(nL)]
        out = calc_leaf_pops_per_layer(
            prev_leaf_pops_per_lai,
            canopy_lai=0,
            layers_lai=layers_lai,
            nP=nP,
            nL=nL,
            growing_populations=[False for _ in range(nP)],
        )
        assert_array_equal(out, np.zeros((nL, nP)))

    def test_distribute_lai_to_emerged_leaf_pops(self, nL, nP, growing_populations):
        prev_leaf_pops_per_lai = np.zeros((nL, nP))
        increase_in_canopy_lai = 99 * nL
        growing_populations_val = [False for _ in range(nP)] if growing_populations == "none" \
            else [True] + [False for _ in range(nP - 1)] if growing_populations == "some" \
            else [True for _ in range(nP)] if growing_populations == "all" else None

        layers_lai = [99 for _ in range(nL)] if growing_populations != "none" else [
            0 for _ in range(nL)]
        out = calc_leaf_pops_per_layer(
            prev_leaf_pops_per_lai.tolist(),
            canopy_lai=increase_in_canopy_lai,
            layers_lai=layers_lai,
            nP=nP,
            nL=nL,
            growing_populations=growing_populations_val,
        )
        assert len(growing_populations_val)
        for iP in range(nP):
            if growing_populations_val[iP]:
                try:
                    assert np.sum(out[:, iP]) > 0
                except AssertionError:
                    print(out)
                    raise AssertionError()


example_species = SpeciesConfig(
    PRESET=SpeciesPresets.WHEAT_SPRING,
    f_tt_emr=0.05,
    f_tt_veg=700 / 2000,
    f_tt_rep=1200 / 2000,
    f_phen_min=0.1,
)


class TestGetGrowingPopulationsRangeFromConfig:
    species_config = example_species
    species_config.key_dates.sowing = 20
    species_config.key_lengths_td.sowing_to_end = 2000
    species_config.f_phen_min = 0.1
    model_config = ModelConfig(
        phenology_method=PhenologyMethods.SEASON_FRACTION,
    )

    def processed_config(self, nP): return process_phenology_config(
        model_config=self.model_config,
        species_config=self.species_config,
        external_data=None,
        td_base_temperature=0,
        nP=nP,
    )[1]

    def test_should_get_growing_populations_single_population(self):
        nP = 1
        species_config = self.processed_config(nP)
        td = np.arange(0, 2200 + 20, 20)
        growing_populations, emerged_leaf_populations_count = get_growing_populations_range_from_config(
            species_config,
            nP,
            td,
        )
        assert len(growing_populations) == len(td)
        assert len(growing_populations[0]) == nP

    def test_should_get_growing_populations_multi_population(self):
        nP = 3
        species_config = self.processed_config(nP)
        td = np.arange(0, 2200 + 20, 20)
        growing_populations, emerged_leaf_populations_count = get_growing_populations_range_from_config(
            species_config,
            nP,
            td,
        )

        assert len(growing_populations) == len(td)
        assert len(growing_populations[0]) == nP

    def test_should_not_have_overlapping_growing_populations(self):
        nP = 3
        species_config = self.processed_config(nP)
        td = np.arange(0, 2200 + 20, 20)
        growing_populations, emerged_leaf_populations_count = get_growing_populations_range_from_config(
            species_config,
            nP,
            td,
        )
        total_growing_populations = np.sum(growing_populations, axis=1)
        assert len(total_growing_populations) == len(td)
        assert max(total_growing_populations) == 1

@pytest.mark.parametrize([
    'nL',
    'nP',
    'no_emerged_pops',
    'canopy_lai',
], [
    (4, 3, 0, 1),
    (4, 3, 1, 1),
    (4, 3, 2, 1),
    (1, 3, 2, 1),
    (1, 1, 1, 1),
    (4, 1, 1, 1),
    (1, 1, 0, 1),
    (4, 1, 0, 1),
])
class TestDistributeCanopyLaiToLeafPops:

    def test_works(self, nL, nP, no_emerged_pops, canopy_lai):
        distribute_canopy_lai_to_leaf_pops(nL, nP, no_emerged_pops, canopy_lai)

    def test_has_correct_shape(self, nL, nP, no_emerged_pops, canopy_lai):
        out = distribute_canopy_lai_to_leaf_pops(nL, nP, no_emerged_pops, canopy_lai)
        assert np.array(out).shape == (nL, nP)

    def test_matches_canopy_lai(self, nL, nP, no_emerged_pops, canopy_lai):
        out = distribute_canopy_lai_to_leaf_pops(nL, nP, no_emerged_pops, canopy_lai)
        if no_emerged_pops > 0:
            assert np.array(out).sum() == canopy_lai

    def test_is_zero_when_no_growing(self, nL, nP, no_emerged_pops, canopy_lai):
        out = distribute_canopy_lai_to_leaf_pops(nL, nP, no_emerged_pops, canopy_lai)
        if no_emerged_pops == 0:
            assert np.array(out).sum() == 0

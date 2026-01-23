"""Test the phenology stage functions."""

from do3se_phenology.config import PhenologyKeyLengths
from do3se_phenology.state import LeafPhenologyStage, PhenologyStage
from do3se_phenology.phenology_stages import (
    get_leaf_phenology_stage_td,
    get_plant_phenology_stage_td,
)



class TestPlantPhenologyStage:
    def test_gets_correct_phenology_stages(self):
        tds = [6,20]
        key_lengths = PhenologyKeyLengths(
            emerg_to_astart=5,
            emerg_to_end=19,
        )
        phenology_stage = PhenologyStage.EMERGED
        for td in tds:
            next_phenology_stage = get_plant_phenology_stage_td(
                td,key_lengths,phenology_stage,
            )

            assert next_phenology_stage == phenology_stage + 1
            phenology_stage = next_phenology_stage

    def test_returns_input_phenology_stage_before_emerged(self):
        # Sowing and emergence are handled seperately.
        tds = [6,20]
        key_lengths = PhenologyKeyLengths(
            emerg_to_astart=5,
            emerg_to_end=19,
        )
        phenology_stage = PhenologyStage.NOT_SOWN
        for td in tds:
            next_phenology_stage = get_plant_phenology_stage_td(
                td,key_lengths,phenology_stage,
            )
            assert next_phenology_stage == phenology_stage


        phenology_stage = PhenologyStage.SOWN
        for td in tds:
            next_phenology_stage = get_plant_phenology_stage_td(
                td,key_lengths,phenology_stage,
            )
            assert next_phenology_stage == phenology_stage



class TestLeafPhenologyStage:
    def test_gets_correct_phenology_stages(self):
        tds = [6,101, 121, 200]
        t_lem=100
        t_lep=20
        t_lse=20
        phenology_stage = get_leaf_phenology_stage_td(
                0,
                t_lem=t_lem,
                t_lep=t_lep,
                t_lse=t_lse,
            )
        assert phenology_stage == LeafPhenologyStage.NOT_EMERGED
        for td in tds:
            next_phenology_stage = get_leaf_phenology_stage_td(
                td,
                t_lem=t_lem,
                t_lep=t_lep,
                t_lse=t_lse,
            )

            assert next_phenology_stage == phenology_stage + 1
            phenology_stage = next_phenology_stage

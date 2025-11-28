from math import isclose, floor
import numpy as np
from copy import deepcopy
from dataclasses import replace

from do3se_phenology.switchboard import (
    PhenologyMethods,
    TimeTypes,
    calc_key_dates,
    process_phenology_config,
)
from do3se_phenology.config import (
    LeafFPhenMethods,
    PlantEmergeMethod,
    SowingDateMethods,
    DVIMethods,
    KeyDays,
    LAIMethods,
    ModelConfig,
    PhenologyKeyDates,
    PhenologyKeyLengths,
    PhenologyLeafKeyLengths,
    SpeciesConfig,
    SpeciesPresets,
    ZeroDayOptions,
)

from do3se_phenology.presets.wheat import SpringWheatMultiplicative, Wheat

# TODO: Tidy up example data
example_SGS = 20
t_astart = 1080
t_egs = 2000

dd = np.array([j for j in range(140) for _ in range(24)])
tdb = np.array([i for i in np.arange(0, 2200 + 20, 20) for _ in range(24)])
day_count = floor(len(tdb) / 24)
td = np.concatenate([tdb, tdb + tdb[-1]])[0 : 140 * 24]
example_SGS_t = next(((d, t) for d, t in zip(dd, td) if d == example_SGS))[1]
example_EGS = next((d for d, t in zip(dd, td) if t >= example_SGS_t + t_egs))

example_data = dict(
    td=td,
    Ts_C=np.array([20 for _ in range(140) for _ in range(24)]),
    dd=dd,
    leaf_fphen=np.interp(
        np.arange(0, 2200 + 20, 20 / 24),
        np.array([0, 1080, 1081, 1240, 1340, 1680, 2000, 2000]) - (t_astart - example_SGS_t),
        [0.0, 0.0, 1.0, 1.0, 1.0, 0.8, 0.0, 0.0],
    ),
)

example_species = SpeciesConfig(
    PRESET=SpeciesPresets.WHEAT_SPRING,
    f_tt_emr=0.05,
    f_tt_veg=700 / 2000,
    f_tt_rep=1200 / 2000,
    f_tt_fst_acc=0,
    f_phen_min=0.1,
)


def default_expected_species_output(nP):
    return replace(
        Wheat,
        PRESET=SpeciesPresets.WHEAT_SPRING,
        f_tt_emr=0.05,
        f_tt_veg=700 / 2000,
        f_tt_rep=1200 / 2000,
        f_phen_min=0.1,
        fphen_intervals=[
            (0.0, 0.1),
            (100.0, 0.1),
            (400.0, 1.0),
            (1240.0, 1.0),
            (2000.0, 0.0),
        ],
        leaf_fphen_intervals=[
            (0.0, 0.0),
            (1080.0, 0.0),
            (1080.0, 1.0),
            (1240.0, 1.0),
            (1340.0, 1.0),
            (1680.0, 0.7),
            (2000.0, 0.0),
            (2000.0, 0.0),
        ],
        dvi_interval=[
            (0, -1.000000001),
            (100, 0.0),
            (800.0, 1.0),
            (2000.0, 2.0),
        ],
        key_dates=PhenologyKeyDates(
            sowing=20,
            emergence=None,
            harvest=None,
            Astart=None,
            Aend=None,
            mid_anthesis=None,
        ),
        key_dates_td=PhenologyKeyDates(
            sowing=0.0,
            emergence=100.0,
            harvest=2000.0,
            Astart=1080.0,
            Aend=2000.0,
            mid_anthesis=1240.0,
        ),
        key_lengths_td=PhenologyKeyLengths(
            sowing_to_emerge=100.0,
            sowing_to_f_phen_b=400.0,
            sowing_to_f_phen_c=1240.0,
            sowing_to_astart=1080,
            emerg_to_astart=980,
            sowing_to_end=2000.0,
            emerg_to_end=1900.0,
            emerg_to_veg=700.0,
            veg_to_harvest=1200.0,
        ),
        key_lengths_leaf_td=PhenologyLeafKeyLengths(
            tl=1900,
            # TODO: Should equal the td between plant emerge and flag emerge divided by no of pops
            tl_em=(980 - 280) / nP,
            tl_ma=1200,
            tl_ep=round(1200 * (0.31 / 0.46), 4),
            tl_se=round(1200 * (0.15 / 0.46), 4),
            leaf_f_phen_e=None,
            leaf_f_phen_g=None,
            leaf_f_phen_h=None,
            leaf_f_phen_i=None,
            # TODO: plant emerg to leaf emer should be different
            # for each population
            plant_emerg_to_leaf_emerg=None,
            leaf_emerg_to_leaf_fst_acc=None,
        ),
        key_lengths_flag_leaf_td=PhenologyLeafKeyLengths(
            tl=920.0 + 280.0,
            tl_em=280.0,
            tl_ma=920.0,
            tl_ep=620.0,
            tl_se=300.0,
            leaf_f_phen_e=160.0,  # 0.08
            leaf_f_phen_g=100.0,  # 0.05
            leaf_f_phen_h=440.0,  # 0.22
            leaf_f_phen_i=760.0,  # 0.38
            plant_emerg_to_leaf_emerg=980 - 280,
            leaf_emerg_to_leaf_fst_acc=0,
        ),
    )


default_expected_species_output_astart = replace(
    Wheat,
    PRESET=SpeciesPresets.WHEAT_SPRING,
    f_tt_emr=0.05,
    f_tt_veg=700 / 2000,
    f_tt_rep=1200 / 2000,
    f_phen_min=0.1,
    fphen_intervals=[
        (-1080.0, 0.1),
        (-980.0, 0.1),
        (-680.0, 1.0),
        (160.0, 1.0),
        (920.0, 0.0),
    ],
    leaf_fphen_intervals=[
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 1.0),
        (160.0, 1.0),
        (260.0, 1.0),
        (600.0, 0.7),
        (920.0, 0.0),
        (920.0, 0.0),
    ],
    dvi_interval=[
        (-1080, -1.000000001),
        (-980, 0.0),
        (-280.0, 1.0),
        (920.0, 2.0),
    ],
    key_dates=PhenologyKeyDates(
        sowing=None,
        emergence=None,
        harvest=None,
        Astart=20,
        Aend=None,
        mid_anthesis=None,
    ),
    key_dates_td=PhenologyKeyDates(
        sowing=-1080,
        emergence=-980.0,
        harvest=920.0,
        Astart=0.0,
        Aend=920.0,
        mid_anthesis=160.0,
    ),
    key_lengths_td=PhenologyKeyLengths(
        sowing_to_emerge=100.0,
        sowing_to_f_phen_b=400.0,
        sowing_to_f_phen_c=1240.0,
        sowing_to_astart=1080,
        emerg_to_astart=980,
        sowing_to_end=2000.0,
        emerg_to_end=1900.0,
        emerg_to_veg=700.0,
        veg_to_harvest=1200.0,
    ),
    key_lengths_leaf_td=PhenologyLeafKeyLengths(
        tl=1900,
        tl_em=980
        - 280,  # TODO: Should equal the td between plant emerge and flag emerge divided by no of pops
        tl_ma=1200,
        tl_ep=round(1200 * (0.31 / 0.46), 4),
        tl_se=round(1200 * (0.15 / 0.46), 4),
        leaf_f_phen_e=None,
        leaf_f_phen_g=None,
        leaf_f_phen_h=None,
        leaf_f_phen_i=None,
        # TODO: plant emerg to leaf emer should be different
        # for each population
        plant_emerg_to_leaf_emerg=None,
    ),
    key_lengths_flag_leaf_td=PhenologyLeafKeyLengths(
        tl=920.0 + 280.0,
        tl_em=280.0,
        tl_ma=920.0,
        tl_ep=620.0,
        tl_se=300.0,
        leaf_f_phen_e=160.0,  # 0.08
        leaf_f_phen_g=100.0,  # 0.05
        leaf_f_phen_h=440.0,  # 0.22
        leaf_f_phen_i=760.0,  # 0.38
        # TODO: plant emerg to leaf emer should be different
        # for each population
        plant_emerg_to_leaf_emerg=980 - 280,
    ),
)

default_expected_model_output = ModelConfig(
    phenology_method=PhenologyMethods.FPHEN_THERMAL_TIME,
    dvi_method=DVIMethods.DISABLED,
    LAI_method=LAIMethods.ESTIMATE_TOTAL,
    time_type=TimeTypes.THERMAL_TIME,
    sgs_time_type=TimeTypes.JULIAN_DAY,
    sgs_key_day=KeyDays.SOWING_DAY,
)


class SwitchBoardTestBase:
    model_config = ModelConfig()
    species_config = SpeciesConfig(PRESET=SpeciesPresets.WHEAT_SPRING)
    external_data = None
    nP = 1
    expected_model_config = deepcopy(default_expected_model_output)
    expected_species_config = deepcopy(default_expected_species_output(nP))

    def _run(
        self, model_config=None, species_config=None, external_data: dict | None = None, nP=nP
    ):
        external_data_in = external_data or example_data
        return process_phenology_config(
            model_config=model_config or self.model_config,
            species_config=species_config or self.species_config,
            external_data=external_data_in,
            td_base_temperature=0,
            nP=nP,
        )

    def test_returns_expected_output(self):
        out_model_config, out_species_config = self._run()
        assert out_model_config == self.expected_model_config
        assert out_species_config == self.expected_species_config

    def test_defines_all_key_dates(self):
        out_model_config, out_species_config = self._run()
        assert out_species_config.key_lengths_td == self.expected_species_config.key_lengths_td
        assert out_species_config.key_lengths == self.expected_species_config.key_lengths

        assert out_species_config.key_dates_td == self.expected_species_config.key_dates_td
        assert out_species_config.key_dates == self.expected_species_config.key_dates

    def test_defines_all_key_dates_leaf(self):
        out_model_config, out_species_config = self._run()
        assert out_species_config.key_lengths_leaf == self.expected_species_config.key_lengths_leaf
        assert (
            out_species_config.key_lengths_leaf_td
            == self.expected_species_config.key_lengths_leaf_td
        )

    def test_defines_all_key_intervals(self):
        out_model_config, out_species_config = self._run()
        assert out_species_config.fphen_intervals == self.expected_species_config.fphen_intervals
        assert (
            out_species_config.leaf_fphen_intervals
            == self.expected_species_config.leaf_fphen_intervals
        )
        assert out_species_config.dvi_interval == self.expected_species_config.dvi_interval

    def test_should_define_initial_thermal_time(self):
        pass

    def test_should_set_sgs_and_egs(self):
        # TODO: We need to know the date of EGS?
        pass

    def test_should_have_set_zero_day(self):
        out_model_config, out_species_config = self._run()
        assert out_model_config.zero_day is not None

    def test_values_should_be_sensible(self):
        out_model_config, out_species_config = self._run()
        if self.model_config.time_type == TimeTypes.THERMAL_TIME:
            assert out_species_config.key_lengths_leaf_td.tl is not None
            assert out_species_config.key_lengths_td.emerg_to_end is not None
            assert (
                out_species_config.key_lengths_leaf_td.tl
                <= out_species_config.key_lengths_td.emerg_to_end
            )

        elif self.model_config.time_type == TimeTypes.JULIAN_DAY:
            # TODO: Check julian day values
            pass

    def test_should_work_with_multiple_leaf_pops(self):
        nP = 3
        if self.model_config.time_type != TimeTypes.THERMAL_TIME:
            return
        emerge_to_flag_emerg = (
            self.expected_species_config.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg
        )
        assert emerge_to_flag_emerg is not None
        assert emerge_to_flag_emerg > 0
        assert self.expected_species_config.f_t_lma is not None
        assert self.expected_species_config.f_t_lep is not None
        assert self.expected_species_config.f_t_lse is not None
        assert self.expected_species_config.f_t_lem is not None
        assert self.expected_species_config.f_t_lem is not None
        assert self.expected_species_config.f_t_lem is not None
        t_lem = emerge_to_flag_emerg / (nP - 1)
        t_lma = t_lem * self.expected_species_config.f_t_lma / self.expected_species_config.f_t_lem
        t_lep = t_lem * self.expected_species_config.f_t_lep / self.expected_species_config.f_t_lem
        t_lse = t_lem * self.expected_species_config.f_t_lse / self.expected_species_config.f_t_lem
        tl = t_lma + t_lem
        assert isclose(t_lep + t_lse, t_lma, abs_tol=1e-3)
        key_lengths_leaf_td = PhenologyLeafKeyLengths(
            tl=tl,
            tl_em=t_lem,
            tl_ma=t_lma,
            tl_ep=round(t_lep, 4),
            tl_se=round(t_lse, 4),
            leaf_f_phen_e=None,  # 0.08
            leaf_f_phen_g=None,  # 0.05
            leaf_f_phen_h=None,  # 0.22
            leaf_f_phen_i=None,  # 0.38
            # TODO: plant emerg to leaf emerg should be different
            # for each population
            plant_emerg_to_leaf_emerg=None,
        )

        out_model_config, out_species_config = self._run(nP=3)

        expected_species_config = deepcopy(self.expected_species_config)
        expected_species_config.key_lengths_leaf_td = key_lengths_leaf_td

        assert out_model_config == self.expected_model_config
        assert out_species_config == expected_species_config

    def test_should_use_preset_if_set(self):
        """Should fill missing values with preset if set."""
        # NOTE: All tests using preset so will fail if not working
        self._run()


class TestSwitchBoardFphenTd(SwitchBoardTestBase):
    """Test phenology from fphen thermal time intervals."""

    model_config = ModelConfig(
        phenology_method=PhenologyMethods.FPHEN_THERMAL_TIME,
    )
    species_config = replace(
        deepcopy(example_species),
        leaf_f_phen_a=0.3,
        leaf_f_phen_b=0.7,
    )
    species_config.key_lengths_td = replace(
        species_config.key_lengths_td,
        sowing_to_emerge=100.0,
        sowing_to_f_phen_b=400.0,
        sowing_to_f_phen_c=1240.0,
        sowing_to_end=2000.0,
    )

    species_config.key_dates.sowing = 20
    species_config.key_dates_td.Astart = 1080
    species_config.key_lengths_flag_leaf_td = replace(
        species_config.key_lengths_flag_leaf_td,
        leaf_f_phen_e=160.0,  # 0.08
        leaf_f_phen_g=100.0,  # 0.05
        leaf_f_phen_h=440.0,  # 0.22
        leaf_f_phen_i=760.0,  # 0.38
    )

    expected_model_config = replace(
        deepcopy(default_expected_model_output),
        phenology_method=PhenologyMethods.FPHEN_THERMAL_TIME,
    )
    expected_species_config = replace(
        default_expected_species_output(1),
        f_phen_min=0.1,
    )
    expected_species_config.key_lengths_flag_leaf_td = replace(
        expected_species_config.key_lengths_flag_leaf_td,
        leaf_f_phen_e=160.0,  # 0.08
        leaf_f_phen_g=100.0,  # 0.05
        leaf_f_phen_h=440.0,  # 0.22
        leaf_f_phen_i=760.0,  # 0.38
    )

    def test_uses_thermal_time_type(self):
        model_config, series_config = self._run()
        assert model_config.time_type == TimeTypes.THERMAL_TIME


class TestSwitchBoardDayPLF(SwitchBoardTestBase):
    """Test phenology from fphen legacy Julian Day values."""

    model_config = ModelConfig(
        phenology_method=PhenologyMethods.LEGACY_DAY_PLF,
        plant_emerge_method=PlantEmergeMethod.SGS,
        sowing_day_method=SowingDateMethods.INPUT,
        time_type=TimeTypes.JULIAN_DAY,
    )
    species_config = replace(
        deepcopy(SpringWheatMultiplicative),
    )
    expected_model_config = replace(
        deepcopy(model_config),
        phenology_method=PhenologyMethods.LEGACY_DAY_PLF,
        plant_emerge_method=PlantEmergeMethod.SGS,
        time_type=TimeTypes.JULIAN_DAY,
    )
    expected_species_config = replace(
        deepcopy(species_config),
        key_dates=replace(
            deepcopy(species_config.key_dates),
            Astart=145,  # sowing + sowing_to_astart
            Aend=species_config.key_dates.harvest,  # Same as harvest
        ),
        key_lengths=replace(
            deepcopy(species_config.key_lengths),
            sowing_to_emerge=0,
            sowing_to_end=197 - 105,  # harvest - sowing
            emerg_to_end=197 - 105 - 0,  # harvest - sowing - sowing_to_emerge
            emerg_to_astart=145,  # Astart - sowing
        ),
        key_lengths_flag_leaf=replace(
            deepcopy(species_config.key_lengths_flag_leaf),
            plant_emerg_to_leaf_emerg=145,  # Same as sowing +
            leaf_emerg_to_fully_grown=20,  # f_phen_1
            fully_grown_to_senescence=2.0,  # 197 - 167
        ),
        leaf_fphen_intervals=([145.0, 145.0, 165.0, 167.0, 197.0, 197.0], [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]),
    )

    def test_sets_sgs_and_egs(self):
        model_config, species_config = self._run()
        assert species_config.key_dates.sowing is not None
        assert species_config.key_dates.harvest is not None

        # These values come from the mapping manual and
        # Simpson et al., (2003). Transboundary Acidification, Eutrophication and Ground Level Ozone in Europe PART I. Unified EMEP Model Description. EMEP Report 1/2003
        # Table 5.1
        assert species_config.key_dates.sowing == 105
        assert species_config.key_dates.harvest == 197

    def test_sets_day_fphen_plf(self):
        model_config, species_config = self._run()
        assert species_config.day_fphen_plf == self.expected_species_config.day_fphen_plf
        assert species_config.key_dates.Astart is not None
        assert species_config.key_dates.Aend is not None
        assert species_config.key_lengths_flag_leaf.plant_emerg_to_leaf_emerg is not None
        assert species_config.key_lengths_flag_leaf.leaf_emerg_to_fully_grown is not None
        assert species_config.key_lengths_flag_leaf.fully_grown_to_senescence is not None

        assert species_config.key_lengths.emerg_to_astart is not None
        assert species_config.key_lengths.emerg_to_end is not None


class TestSwitchBoardDayPLFPlantOnly(TestSwitchBoardDayPLF):
    model_config = ModelConfig(
        phenology_method=PhenologyMethods.LEGACY_DAY_PLF,
        plant_emerge_method=PlantEmergeMethod.SGS,
        sowing_day_method=SowingDateMethods.INPUT,
        time_type=TimeTypes.JULIAN_DAY,
    )
    species_config = replace(
        deepcopy(SpringWheatMultiplicative),
        leaf_f_phen_method=LeafFPhenMethods.F_PHEN,
    )

    expected_model_config = replace(
        deepcopy(TestSwitchBoardDayPLF.expected_model_config),
    )
    expected_species_config = replace(
        deepcopy(TestSwitchBoardDayPLF.expected_species_config),
        leaf_f_phen_method=LeafFPhenMethods.F_PHEN,
        key_lengths=replace(
            deepcopy(TestSwitchBoardDayPLF.species_config.key_lengths),
            sowing_to_astart=40.0,
            sowing_to_end=197 - 105,
            emerg_to_end=197 - 105,  # harvest - sowing
        ),
        key_lengths_flag_leaf=replace(
            deepcopy(species_config.key_lengths_flag_leaf),
            plant_emerg_to_leaf_emerg=0,
        ),
        key_dates=replace(
            deepcopy(TestSwitchBoardDayPLF.species_config.key_dates),
            Astart=None,
            Aend=None,
        ),
        leaf_fphen_intervals=None,
    )

    def test_sets_day_fphen_plf(self):
        model_config, species_config = self._run()
        assert species_config.day_fphen_plf == self.expected_species_config.day_fphen_plf
        assert species_config.key_lengths.emerg_to_end is not None


class TestSwitchBoardLeafFPhenData(SwitchBoardTestBase):
    """Test phenology from leaf f phen data."""

    model_config = ModelConfig(
        phenology_method=PhenologyMethods.LEAF_FPHEN_DATA,
        zero_day=ZeroDayOptions.ASTART,
    )

    species_config = replace(
        deepcopy(example_species),
    )

    external_data = example_data

    expected_model_config = replace(
        deepcopy(default_expected_model_output),
        phenology_method=PhenologyMethods.LEAF_FPHEN_DATA,
        zero_day=ZeroDayOptions.ASTART,
    )

    expected_species_config = replace(
        deepcopy(default_expected_species_output_astart),
        fphen_intervals=None,
        leaf_fphen_intervals=None,
    )


class TestSwitchBoardGrowingSeasonFraction(SwitchBoardTestBase):
    """Test phenology from fraction of growing season."""

    model_config = ModelConfig(
        phenology_method=PhenologyMethods.SEASON_FRACTION,
    )
    species_config = replace(
        deepcopy(example_species),
    )
    species_config.key_dates.sowing = 20
    species_config.key_lengths_td.sowing_to_end = 2000
    species_config.f_phen_min = 0.1

    expected_model_config = replace(
        deepcopy(default_expected_model_output),
        phenology_method=PhenologyMethods.SEASON_FRACTION,
    )
    expected_species_config = replace(
        deepcopy(
            default_expected_species_output(1),
        )
    )


class TestSwitchBoardGrowingSeasonFractionHarvestDay(SwitchBoardTestBase):
    """Test phenology from fraction of growing season."""

    model_config = ModelConfig(
        phenology_method=PhenologyMethods.SEASON_FRACTION,
    )
    species_config = replace(
        deepcopy(example_species),
    )
    species_config.key_dates.sowing = 20
    species_config.key_dates.harvest = 120
    # species_config.key_lengths_td.sowing_to_end = 2000
    species_config.f_phen_min = 0.1

    external_data = example_data

    expected_model_config = replace(
        deepcopy(default_expected_model_output),
        phenology_method=PhenologyMethods.SEASON_FRACTION,
    )
    expected_species_config = replace(
        deepcopy(
            default_expected_species_output(1),
        )
    )
    expected_species_config.key_dates.harvest = 120


class TestSwitchBoardOverrides(SwitchBoardTestBase):
    """Test switchboard does not override set params.

    We need to know that parameters set in the input config will
    not be overriden by the calculated values.
    """

    model_config = ModelConfig(
        phenology_method=PhenologyMethods.SEASON_FRACTION,
    )
    species_config = replace(
        # We set input config to expected config.
        # Even if we change the input season length it should not override these values.
        deepcopy(default_expected_species_output(1))
    )
    species_config.key_dates.sowing = 20
    # We have set sowing to end to 4000 to show this does not override input
    species_config.key_lengths_td.sowing_to_end = 4000
    species_config.f_phen_min = 0.1
    species_config.key_lengths_leaf_td = deepcopy(
        default_expected_species_output(1).key_lengths_leaf_td
    )

    expected_model_config = replace(
        deepcopy(default_expected_model_output),
        phenology_method=PhenologyMethods.SEASON_FRACTION,
    )
    expected_species_config = replace(deepcopy(default_expected_species_output(1)))
    expected_species_config.key_lengths_leaf_td = deepcopy(
        default_expected_species_output(1).key_lengths_leaf_td
    )

    expected_species_config.key_lengths_td.sowing_to_end = 4000

    def test_should_work_with_multiple_leaf_pops(self):
        # Cannot run this test for this test class as it overrides expected config
        pass
        # return super().test_should_work_with_multiple_leaf_pops()


# TODO: Implement Astart to Aend none thermal time method
# class TestSwitchBoardFractionFlagLeaf(SwitchBoardTestBase):
#     """Test phenology from fraction of flag leaf."""

#     model_config = ModelConfig(
#         phenology_method=PhenologyMethods.FLAG_LEAF_FRACTION,
#     )
#     model_config.zero_day = ZeroDayOptions.ASTART
#     species_config = replace(
#         deepcopy(example_species),
#     )
#     species_config.key_dates.Astart = 200
#     species_config.key_dates.Aend = 353.333
#     species_config.f_phen_min = 0.1

#     expected_model_config = replace(
#         deepcopy(default_expected_model_output),
#         phenology_method=PhenologyMethods.FLAG_LEAF_FRACTION,
#     )
#     expected_species_config = replace(
#         deepcopy(default_expected_species_output),
#     )
#     expected_species_config.key_dates.Astart = 200
#     expected_species_config.key_dates.Aend = 353.333


class TestSwitchBoardFractionFlagLeafTd(SwitchBoardTestBase):
    """Test phenology from fraction of flag leaf using season thermal time length."""

    model_config = ModelConfig(
        phenology_method=PhenologyMethods.FLAG_LEAF_FRACTION,
        sowing_day_method=SowingDateMethods.SKIP,
    )
    model_config.zero_day = ZeroDayOptions.ASTART

    species_config = replace(
        deepcopy(example_species),
    )
    species_config.key_dates.Astart = 200
    species_config.key_lengths_td.sowing_to_end = 2000
    species_config.f_phen_min = 0.1

    expected_model_config = replace(
        deepcopy(default_expected_model_output),
        phenology_method=PhenologyMethods.FLAG_LEAF_FRACTION,
        sowing_day_method=SowingDateMethods.SKIP,
    )
    expected_model_config.zero_day = ZeroDayOptions.ASTART
    expected_species_config = replace(
        deepcopy(default_expected_species_output_astart),
    )
    expected_species_config.key_dates.Astart = 200
    expected_species_config.key_lengths_td.sowing_to_end = 2000
    # We don't need leaf_emerg_to_leaf_fst_acc but it is being set
    expected_species_config.key_lengths_flag_leaf_td.leaf_emerg_to_leaf_fst_acc = 0


def vert_line(ax, x, ylim, **kwargs):
    ax.plot([x, x], ylim, **kwargs)


class TestCalculateKeyDates:
    def test_can_calculate_key_dates_from_ext_data(self):
        nP = 3
        td_data = td
        dd_data = dd

        output_config: SpeciesConfig = calc_key_dates(
            deepcopy(default_expected_species_output(nP)),
            td_data,
            dd_data,
        )

        print(output_config.key_dates)
        print(output_config.key_lengths)
        assert output_config.key_dates.sowing == example_SGS
        assert output_config.key_dates.harvest == example_EGS
        assert output_config.key_lengths.sowing_to_emerge is not None
        assert output_config.key_lengths.sowing_to_f_phen_b is not None
        assert output_config.key_lengths.sowing_to_f_phen_c is not None
        assert output_config.key_lengths.sowing_to_end is not None

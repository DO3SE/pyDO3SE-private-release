"""Default DO3SE processes.

This should be the only place the model state is modified

NOTE: Always set iLC inside iL
i.e
[[iLC for iLC in range(nLC)]
for iL in range(nL)]

"""

import os
from pathlib import Path
from typing import Callable, List, Optional
from copy import deepcopy
from data_helpers.list_helpers import flatten_list

from proflow.ProcessRunnerCls import advance_time_step_process
from proflow.Objects.Process import Process, ProcessType
from proflow.Objects.Interface import I
from proflow.Switch import switch
from proflow.helpers import NOT_IMPLEMENTED_PROCESS, set_value, skip
from proflow.logger import log_values, merge_logs

from do3se_phenology.config import (
    PlantEmergeMethod,
    FPhenMethods,
    LeafFPhenMethods,
    PhenologyMethods,
    TimeTypes,
    SAILCMethods,
)
from do3se_phenology.state import LeafPhenologyStage, PhenologyStage
from pyDO3SE.Model_State.model_state_loader import dump_state_to_file

from pyDO3SE.util.error_handling import ConfigError
from pyDO3SE.constants.enums import GAS
from pyDO3SE.constants.physical_constants import DRATIO_H20_O3, DRATIO_O3_CO2, T0, DIFF_O3
from pyDO3SE.Config.Config_Shape import Config_Shape, LeafFPhenAnetInfluence
from pyDO3SE.Model_State.Model_State import Model_State_Shape
from pyDO3SE.External_State.External_State_Config import (
    InputMethod, ParSunShadeMethods, ThermalTimeMethods,
)
from pyDO3SE.Config.ConfigEnums import (
    CanopyHeightMethods,
    DVIMethods,
    FVPDMethods,
    LAIMethods,
    LayerLAIDistributionMethods,
    OzoneDepositionMethods,
    SenescenceFunctionMethods,
    TLeafMethods,
    LandCoverType,
)
# Processes
from do3se_met import irradiance as met_irrad_helpers
from do3se_met import wind as met_wind_helpers
from do3se_met import deposition as met_deposition_helpers
from do3se_met import resistance as resistance_helpers
from do3se_phenology import phyllochron_dvi
from do3se_phenology import canopy_structure
from do3se_phenology import f_phen as f_phen_helpers
from do3se_phenology import vernalisation as vernalisation_helpers

from pyDO3SE.constants import model_constants
from pyDO3SE.plugins.gsto import photosynthesis_helpers as pn_helpers
from pyDO3SE.plugins.gsto import helpers as gsto_helpers
from pyDO3SE.plugins.gsto.multiplicative import multiplicative
from pyDO3SE.plugins.gsto.ewert.ewert import ewert_leaf_pop, ewert_leaf_pop_cubic
from pyDO3SE.plugins.gsto.ewert.enums import EwertLoopMethods
from pyDO3SE.plugins.gsto.ewert.ewert_helpers import calc_all_ozone_damage_factors, calc_mean_gsto
from pyDO3SE.plugins.leaf_temperature.de_boeck import get_leaf_temp_de_boeck
from pyDO3SE.plugins.nitrogen.v_cmax import multilayer_vcmax25
from pyDO3SE.plugins.soil_moisture import helpers as SMD_helpers
from pyDO3SE.plugins.soil_moisture.enums import SoilMoistureSource
from pyDO3SE.plugins.soil_moisture import penman_monteith as SMD_PM_helpers
from pyDO3SE.plugins.O3 import helpers as O3_helpers
from pyDO3SE.plugins.carbon_allocation import calculations as carbon_calcs
from pyDO3SE.plugins.carbon_allocation.conversions import umol_c_to_kg_c

from pyDO3SE.Pipelines.validation_processes import setup_validation_processes
from .es_hour_processes import (
    accumulate_hourly_temperature_process,
    accumulate_precipitation_process,
    calc_effective_temperature_process,
    calc_photoperiod_process,
    calc_photoperiod_factor_process,
    calculate_daily_thermal_time_process,
    calculate_relative_photoperiod_process,
    store_accumulate_precipitation_process,
    set_thermal_time_process,
)


def lget(v: Optional[List[any]], i: int, fallback_value=None) -> any:
    try:
        return v[i]
    except (TypeError, IndexError):
        return fallback_value


def nullop():
    return None


def tag_process(comment):
    return Process(
        func=nullop,
        comment=comment,
    )


def perNlc(nLC, fn, *args, **kwargs):
    return [fn(iLC, *args, **kwargs) for iLC in range(nLC)],


def raise_error(e):
    """Allow raising an error in a lambda function."""
    raise e


def accumulate(**kwargs):
    out = sum(kwargs.values())
    return out


# ======================== PHENOLOGY ================= #

def calc_DVI_process(iLC, dvi_method: DVIMethods) -> Process:
    return switch(
        gate=dvi_method,
        comment="Choose DVI method",
        options={
            DVIMethods.DISABLED: [],
            DVIMethods.JULES: [
                Process(
                    func=phyllochron_dvi.calc_dvi,
                    comment="Calculate the development index",
                    config_inputs=lambda config, iLC=iLC: [
                        I(config.Land_Cover.parameters[iLC].phenology.key_dates.sowing,
                            as_='sowing_day'),
                        I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.sowing_to_emerge, as_='tt_emr'),
                        I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.emerg_to_veg, as_='tt_veg'),
                        I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.veg_to_harvest, as_='tt_rep'),
                    ],
                    state_inputs=lambda state, iLC=iLC: [
                        I(state.temporal.dd, as_='dd'),
                        I(state.prev_hour.canopy_component[iLC].dvi, as_='prev_dvi'),
                        I(state.canopy_component[iLC].rpe, as_='rpe'),
                        I(state.canopy_component[iLC].t_eff, as_='t_eff'),
                    ],
                    state_outputs=lambda result, iLC=iLC: [
                        (result, f'canopy_component.{iLC}.dvi'),
                    ],
                ),
            ],
            DVIMethods.THERMAL_TIME: [
                Process(
                    func=phyllochron_dvi.calc_dvi_tt_PLF,
                    comment="Calculate the development index",
                    config_inputs=lambda config, iLC=iLC: [
                        I(config.Land_Cover.parameters[iLC].phenology.dvi_interval,
                            as_='dvi_interval'),
                    ],
                    state_inputs=lambda state, iLC=iLC: [
                        I(state.canopy_component[iLC].td_v, as_='td'),
                    ],
                    state_outputs=lambda result, iLC=iLC: [
                        (result, f'canopy_component.{iLC}.dvi'),
                    ],
                ),
            ],
        },
    )


def calc_phyllochron_process(iLC) -> Process:
    return Process(
        func=phyllochron_dvi.calc_phyllochron,
        comment="Calculate the phyllochron",
        state_inputs=lambda state: [
            I(state.prev_hour.external_met.photoperiod - state.external_met.photoperiod, as_='dl'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.phyllochron'),
        ],
    )


def calc_if_plant_is_sown_process(
    iLC: int,
) -> Process:
    def constant_plant_sown_check(dd, td, sowing_day, t_sowing, phenology_stage):
        """Phenology stage must be SOWN for plant to emerge."""
        return phenology_stage if phenology_stage >= PhenologyStage.SOWN \
            else PhenologyStage.SOWN if (
                dd > sowing_day if sowing_day is not None else
                td > t_sowing if t_sowing is not None
                else raise_error(ValueError("Must set sowing day"))) else phenology_stage

    return Process(
        func=constant_plant_sown_check,
        comment="Set plant is sown when dd > sowing_day",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].phenology.key_dates.sowing,
                as_='sowing_day'),
            I(config.Land_Cover.parameters[iLC].phenology.key_dates_td.sowing,
                as_='t_sowing'),
        ],
        state_inputs=lambda state: [
            I(state.temporal.dd, as_='dd'),
            I(state.canopy_component[iLC].td, as_='td'),
            I(state.canopy_component[iLC].phenology.phenology_stage, as_="phenology_stage"),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
        ],
    )


def calc_if_plant_has_emerged_process(
    iLC: int,
    plant_emerge_method: PlantEmergeMethod,
) -> Process:
    def constant_plant_emerged_method(td, t_emerge, dd, emerge_day, phenology_stage):
        """Phenology stage must be SOWN for plant to emerge."""
        return phenology_stage if PhenologyStage.SOWN != phenology_stage \
            else PhenologyStage.EMERGED if (td > t_emerge if t_emerge is not None
                                            else dd > emerge_day if emerge_day is not None
                                            else raise_error(ValueError("Must set t_emerge or emerge_day"))) else phenology_stage
    return switch(
        gate=plant_emerge_method,
        comment="Choose method for defining if plant has emerged",
        options={
            PlantEmergeMethod.CONSTANT: Process(
                func=constant_plant_emerged_method,
                comment="Set plant has emerged when td > t_emerge or dd > emerge_day",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_dates_td.emergence, as_='t_emerge'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_dates.emergence,
                      as_='emerge_day'),
                ],
                state_inputs=lambda state: [
                    I(state.temporal.dd, as_='dd'),
                    I(state.canopy_component[iLC].td_v, as_='td'),
                    I(state.canopy_component[iLC].phenology.phenology_stage, as_="phenology_stage"),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
                ],
            ),
            PlantEmergeMethod.SGS: Process(
                func=lambda SGS, dd, phenology_stage: PhenologyStage.EMERGED if (
                    phenology_stage == PhenologyStage.SOWN and dd > SGS) else phenology_stage,
                comment="Set plant has emerged when dd > SGS",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].phenology.key_dates.sowing, as_='SGS'),
                ],
                state_inputs=lambda state: [
                    I(state.temporal.dd, as_='dd'),
                    I(state.canopy_component[iLC].phenology.phenology_stage, as_="phenology_stage"),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
                ],
            ),
            PlantEmergeMethod.DVI: Process(
                func=lambda dvi, phenology_stage: PhenologyStage.EMERGED if (
                    phenology_stage == PhenologyStage.SOWN and dvi > 0) else phenology_stage,
                comment="Set plant has emerged when dvi > 0",
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_component[iLC].dvi, as_='dvi'),
                    I(state.canopy_component[iLC].phenology.phenology_stage, as_="phenology_stage"),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
                ],
            ),
            PlantEmergeMethod.FPHEN: Process(
                func=lambda f_phen, phenology_stage: PhenologyStage.EMERGED if (
                    phenology_stage == PhenologyStage.SOWN and f_phen > 0) else phenology_stage,
                comment="Set plant has emerged when f_phen > 0",
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_layer_component[0][iLC].gsto_params.f_phen, as_='f_phen'),
                    I(state.canopy_component[iLC].phenology.phenology_stage, as_="phenology_stage"),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
                ],
            ),
            PlantEmergeMethod.LEAF_F_PHEN: Process(
                # TODO: This should actually be true when td_dd > 0 as set in state init
                func=lambda leaf_f_phen, phenology_stage: PhenologyStage.EMERGED if (
                    phenology_stage == PhenologyStage.SOWN and leaf_f_phen > 0) else phenology_stage,
                comment="Set plant has emerged when leaf_fphen > 0",
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_layer_component[0]
                        [iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
                    I(state.canopy_component[iLC].phenology.phenology_stage, as_="phenology_stage"),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
                ],
            ),
        },
    )


def calc_if_flag_leaf_has_emerged_process(
    iLC: int,
    flag_leaf_emerge_method: PlantEmergeMethod,
) -> Process:
    return switch(
        gate=flag_leaf_emerge_method,
        comment="Choose method for defining if flag leaf has emerged",
        options={
            PlantEmergeMethod.CONSTANT: Process(
                func=lambda td, t_emerg, dd, emerg: td > t_emerg if t_emerg is not None
                else dd > emerg if emerg is not None
                else raise_error(ValueError("Must set key_lengths_flag_leaf.plant_emerg_to_leaf_emerg or key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg")),
                comment="Set flag leaf has emerged when td or dd pass flag leaf emergence date",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_lengths_flag_leaf.plant_emerg_to_leaf_emerg, as_='emerg'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg, as_='t_emerg'),
                ],
                state_inputs=lambda state, iLC=iLC: [
                    I(state.temporal.dd, as_='dd'),
                    I(state.canopy_component[iLC].td_dd, as_='td'),

                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.flag_has_emerged'),
                ],
            ),
            PlantEmergeMethod.SGS: NOT_IMPLEMENTED_PROCESS("SGS method not implemented for flag leaf emergence"),
            PlantEmergeMethod.DVI: NOT_IMPLEMENTED_PROCESS("DVI method not implemented for flag leaf emergence"),
            PlantEmergeMethod.FPHEN: Process(
                # Note we use leaf f phen
                func=lambda leaf_f_phen: leaf_f_phen > 0,
                comment="Set flag leaf has emerged when leaf_f_phen > 0",
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_layer_component[0]
                      [iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.flag_has_emerged'),
                ],
            ),
            PlantEmergeMethod.LEAF_F_PHEN: Process(
                # NOTE: This method only defines fully emerged flag leaf
                func=lambda leaf_f_phen: leaf_f_phen > 0,
                comment="Set flag leaf has emerged when leaf_f_phen > 0",
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_layer_component[0]
                      [iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f'canopy_component.{iLC}.flag_has_emerged'),
                ],
            ),
        },
    )


def calc_emergence_rate_process(nP: int, iLC: int) -> Process:
    return Process(
        func=phyllochron_dvi.calc_emergence_rate,
        comment="Calculate the emergence rate",
        gate=nP > 1,
        # TODO Add phyllochron version here
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.nP, as_='nP'),
            I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg,
              as_='t_emerg_to_flag_emerg'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.emergence_rate'),
        ],
    )


def calc_emerged_leaf_count_process(nP: int, iLC: int, time_type: LeafPhenologyStage) -> Process:
    return [
        Process(
            func=phyllochron_dvi.calc_emerged_leaf_count,
            comment="Calculate the number of leaf populations that have emerged",
            gate=time_type == TimeTypes.THERMAL_TIME,
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.nP, as_='nP'),
                # Only required if nP == 1
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg,
                  as_="t_emerge_flag"),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].td_dd, as_='td_dd'),
                I(state.canopy_component[iLC].emergence_rate, as_='emergence_rate'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.total_emerged_leaf_populations'),
                # TODO: This should set phenology stage
                # *[(True, f'canopy_component_population.{iLC}.{iP}.phenology.phenology_stage' >= LeafPhenologyStage.GROWING)
                # for iP in range(result)],
            ],
        ),
        Process(
            func=phyllochron_dvi.calc_emerged_leaf_count,
            comment="Calculate the number of leaf populations that have emerged",
            gate=time_type == TimeTypes.JULIAN_DAY,
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.nP, as_='nP'),
                # Only required if nP == 1
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf.plant_emerg_to_leaf_emerg,
                  as_="t_emerge_flag"),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].emergence_rate, as_='emergence_rate'),
                I(state.temporal.dd, as_='td_dd'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.total_emerged_leaf_populations'),
            ],
        ),
        Process(
            func=lambda emerged_populations, phenology_stages: [
                ps if ps > LeafPhenologyStage.GROWING \
                else LeafPhenologyStage.GROWING if i < emerged_populations \
                else LeafPhenologyStage.NOT_EMERGED \
                for i, ps in enumerate(phenology_stages)],
            comment="Set leaf phenology stages",
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].total_emerged_leaf_populations,
                  as_='emerged_populations'),
                I([state.canopy_component_population[iLC]
                   [iP].phenology.phenology_stage for iP in range(nP)], as_='phenology_stages'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                *[(stage, f'canopy_component_population.{iLC}.{iP}.phenology.phenology_stage')
                  for iP, stage in enumerate(result)],
            ],
        )
    ]
    # TODO: Add phyllochron method
    # Process(
    #             func=lambda acc, phyllochron, td, t_emerg, multiplier: acc +
    #             multiplier / phyllochron if td > t_emerg else 0,
    #             gate=nP > 1,
    #             comment="Calculate number of emerged leaf populations",
    #             config_inputs=lambda config: [
    #                 I(config.Land_Cover.leaf_emergence_multiplier, as_="multiplier"),
    #             ],
    #             external_state_inputs=lambda e_state, row_index: [
    #                 I(lget(e_state.td, row_index), as_="td"),
    #             ],
    #             state_inputs=lambda state, iLC=iLC: [
    #                 I(state.canopy_component[iLC].total_emerged_leaf_populations, as_="acc"),
    #                 I(state.canopy_component[iLC].phyllochron, as_='phyllochron'),
    #                 I(state.canopy_component[iLC].t_emerg, as_="t_emerg"),
    #             ],
    #             state_outputs=lambda result, iLC=iLC: [
    #                 (result, f"canopy_component.{iLC}.total_emerged_leaf_populations"),
    #             ],
    #         )


def get_phenology_stage_process_td(nP: int, iLC: int) -> List[Process]:
    return [
        [Process(
            func=phyllochron_dvi.get_leaf_phenology_stage_td,
            comment="Get leaf phenology stage for population {iP}",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_leaf_td.tl_em, as_="t_lem"),
            ],
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I(state.canopy_component[iLC].td_dd_leaf_pops[iP], as_='td_dd'),
                I(state.canopy_component_population[iLC][iP].t_lep_limited, as_='t_lep'),
                I(state.canopy_component_population[iLC][iP].t_lse_limited, as_='t_lse'),
            ],
            state_outputs=lambda result, iLC=iLC, iP=iP: [
                (result, f'canopy_component_population.{iLC}.{iP}.phenology.phenology_stage'),
            ],
        ) for iP in range(nP - 1)],
        Process(
            func=phyllochron_dvi.get_leaf_phenology_stage_td,
            comment="Get leaf phenology stage for flag leaf(td)",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_em, as_="t_lem"),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].td_dd_leaf_pops[-1], as_='td_dd'),
                I(state.canopy_component_population[iLC][-1].t_lep_limited, as_='t_lep'),
                I(state.canopy_component_population[iLC][-1].t_lse_limited, as_='t_lse'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component_population.{iLC}.{-1}.phenology.phenology_stage'),
            ],
        ),
        Process(
            func=phyllochron_dvi.get_plant_phenology_stage_td,
            comment="Get phenology stage for plant",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td, as_="key_lengths"),
            ],
            state_inputs=lambda state, iLC=iLC: [
                # Thermal time since plant has been sown
                I(state.canopy_component[iLC].td_dd, as_='td_dd'),
                I(state.canopy_component[iLC].phenology.phenology_stage, as_='phenology_stage'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
            ],
        ),
    ]


def get_phenology_stage_process(nP: int, iLC: int, leaf_f_phen_method: LeafFPhenMethods) -> List[Process]:
    assert nP == 1, "get_phenology_stage_process NOT IMPLEMENTED for nP > 1"
    return [
        Process(
            func=phyllochron_dvi.get_leaf_phenology_stage,
            comment="Get leaf phenology stage for flag leaf(dd)",
            gate=leaf_f_phen_method not in [LeafFPhenMethods.DISABLED, LeafFPhenMethods.F_PHEN],
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].phenology.key_dates.sowing, as_="dd_at_sowing"),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths.sowing_to_emerge,
                  as_="sowing_to_plant_emerg"),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf.plant_emerg_to_leaf_emerg,
                  as_="plant_emerg_to_leaf_emerg"),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf.leaf_emerg_to_astart,
                    as_="leaf_emerg_to_astart"),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf.astart_to_senescence,
                    as_="Astart_to_senescence"),
            ],
            state_inputs=lambda state: [
                I(state.temporal.dd, as_='dd'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component_population.{iLC}.{-1}.phenology.phenology_stage'),
            ],
        ),
        Process(
            func=phyllochron_dvi.get_plant_phenology_stage,
            comment="Get phenology stage for plant",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths, as_="key_lengths"),
            ],
            state_inputs=lambda state: [
                I(state.temporal.dd, as_='dd'),
                I(state.canopy_component[iLC].phenology.phenology_stage, as_='phenology_stage'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.phenology.phenology_stage'),
            ],
        ),
    ]


def get_growing_populations_process(nP: int, iLC: int) -> Process:
    return [
        # TODO: Can probably optimize this by getting from process above,
        Process(
            func=phyllochron_dvi.get_growing_populations,
            comment="Define which populations are growing",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_em,
                  as_="flag_leaf_t_lem"),
                I([config.Land_Cover.parameters[iLC].phenology.key_lengths_leaf_td.tl_em for _ in range(
                    nP - 1)], as_="leaf_population_t_lems"),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].td_dd_leaf_pops, as_='td_dd_emerg'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.growing_populations'),
                # *[(result[iP], f'canopy_component_population.{iLC}.{iP}.is_growing') for iP in range(nP)],
            ],
        ),
        # Process(
        #     func=lambda growing_populations, phenology_stages: [
        #         LeafPhenologyStage.GROWING if growing else ps
        #                     for growing, ps in zip(growing_populations, phenology_stages)],
        #     comment="Define phenology state of growing populations",
        #     state_inputs=lambda state, iLC=iLC: [
        #         I(state.canopy_component[iLC].growing_populations, as_='growing_populations'),
        #         I([state.canopy_component_population[iLC][iP].phenology.phenology_stage for iP in range(nP)], as_='phenology_stages'),
        #     ],
        #     state_outputs=lambda result, iLC=iLC: [
        #         *[(stage, f'canopy_component_population.{iLC}.{iP}.phenology.phenology_stage')
        #         for iP, stage in enumerate(result)],
        #     ],
        # ),
    ]

# =========================== CANOPY STRUCTURE =========================== #


def calc_canopy_height_process(nL: int, nLC: int, height_method: str) -> List[Process]:
    """Calculate the height of each layer.

    Canopy height is either an input or taken from the "primary land cover".
    """
    return [
        switch(
            gate=height_method,
            comment="Get canopy height",
            options={
                CanopyHeightMethods.CONSTANT: Process(
                    func=canopy_structure.height_method_constant,
                    comment="Use height of primary land cover if height method is constant",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.parameters[config.Land_Cover.primary_LC].height,
                          as_='primary_land_cover_height')
                    ],
                    state_outputs=lambda result: [(result, 'canopy.canopy_height')],
                ),
                CanopyHeightMethods.INPUT: Process(
                    func=set_value,
                    comment="Set canopy height from input",
                    external_state_inputs=lambda e_state, row_index: [
                        I(lget(e_state.h, row_index), as_='canopy_height'),
                    ],
                    state_outputs=lambda result: [
                        (result['canopy_height'], 'canopy.canopy_height')],
                ),
                CanopyHeightMethods.CARBON: Process(
                    func=lambda plant_heights: max(0.001, sum(plant_heights) / len(plant_heights)),
                    comment="Set canopy height to max plant height",
                    state_inputs=lambda state: [
                        I([state.canopy_component[iLC].plant_height for iLC in range(nLC)], as_="plant_heights"),
                    ],
                    state_outputs=lambda result: [
                        (result, f'canopy.canopy_height'),
                    ],
                ),
            },
        ),
        [
            Process(
                func=lambda canopy_height, layer_frac: canopy_height * layer_frac,
                comment="Calculate height of top of each layer as a fraction of total canopy height",
                config_inputs=lambda config, iL=iL: [
                    I(config.Land_Cover.layer_height_frac[iL], as_='layer_frac')],
                state_inputs=lambda state: [
                    I(state.canopy.canopy_height, as_='canopy_height'),
                ],
                state_outputs=lambda result, iL=iL: [
                    (result, f'canopy_layers.{iL}.layer_depth'),
                ],
            ) for iL in range(nL)
        ],
        [
            Process(
                func=lambda layer_depths: sum(layer_depths),
                comment="Get layer height as sum of layer depths below",
                state_inputs=lambda state, iL=iL: [
                    I([state.canopy_layers[jL].layer_depth for jL in range(iL + 1)], as_='layer_depths'),
                ],
                state_outputs=lambda result, iL=iL: [
                    (result, f'canopy_layers.{iL}.layer_height'),
                ],
            ) for iL in range(nL)
        ],
    ]


def calc_canopy_LAI_processes(nLC: int, LAI_method: LAIMethods) -> List[Process]:
    """Calculate LAI change in canopy.

    Supplying LAI and SAI for multi-layer multi-component models can be
    awkward.  To generate a nL*nLC array for LAI/SAI, multiple schemes are
    possible:

    - All LAI/SAI values are known and supplied as input.
    - Total LAI/SAI is known and supplied as input in the top-left cell
        (e.g. LAI(1,1)).  Values are divided among layers and land covers
        according to a (normalised) nL*nLC fLAI array.
    - dvi limited constant - constant LAI value but set to 0 when 0 < dvi < 2
        i.e. when not emerged or after harvest
    - Estimate total - Total LAI/SAI is estimated from the properties of the "primary land
        cover" and divided according to fLAI.

    """
    return [
        # TODO: Set canopy LAI instead of per layer LAI
        switch(
            gate=LAI_method,
            comment="Get LAI",
            options={
                LAIMethods.INPUT: Process(
                    func=set_value,
                    comment="Distribute input LAI into layers",
                    external_state_inputs=lambda e_state, row_index: [
                        I(lget(e_state.LAI, row_index), as_='LAI'),
                    ],
                    state_outputs=lambda result: [
                        (result["LAI"], 'canopy.LAI_total'),
                        (result["LAI"], 'canopy_component.0.LAI'),
                    ],
                ),
                LAIMethods.INPUT_HOURLY: Process(
                    func=set_value,
                    comment="Distribute input LAI into layers",
                    external_state_inputs=lambda e_state, row_index: [
                        I(lget(e_state.LAI, row_index), as_='LAI'),
                    ],
                    state_outputs=lambda result: [
                        (result["LAI"], 'canopy.LAI_total'),
                        (result["LAI"], 'canopy_component.0.LAI'),
                    ],
                ),
                LAIMethods.CONSTANT: Process(
                    func=set_value,
                    comment="Set LAI to constant value",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.LAI, as_='LAI'),
                    ],
                    state_outputs=lambda result: [
                        (result["LAI"], 'canopy.LAI_total'),
                        (result["LAI"], 'canopy_component.0.LAI'),
                    ],
                ),
                # TODO: Check this method works
                LAIMethods.INPUT_TOTAL: Process(
                    func=set_value,
                    comment="Spread single LAI value to layers and LCs",
                    state_inputs=lambda state: [
                        # TODO: Use primary land cover instead of 0.0
                        I(state.canopy_layer_component[0][0].LAI, as_='LAI'),
                    ],
                    state_outputs=lambda result: [
                        (result["LAI"], 'canopy.LAI_total'),
                        (result["LAI"], 'canopy_component.0.LAI'),
                    ],
                ),
                LAIMethods.DVI_LIMITED_CONSTANT: Process(
                    func=lambda LAI, dvi: LAI if 0 < dvi < 2 else 0,
                    comment="Set LAI to constant value limited to dvi range",
                    config_inputs=lambda config: [
                            I(config.Land_Cover.LAI, as_='LAI'),
                    ],
                    state_inputs=lambda state: [
                        I(state.canopy_component[0].dvi, as_='dvi'),
                    ],
                    state_outputs=lambda result: [
                        (result, 'canopy.LAI_total'),
                        (result, 'canopy_component.0.LAI'),
                    ],
                ),
                LAIMethods.LEAF_F_PHEN_LIMITED_CONSTANT: Process(
                    func=lambda LAI, leaf_f_phen: LAI if 0 < leaf_f_phen else 0,
                    comment="Set LAI to constant value limited to leaf_f_phen range",
                    config_inputs=lambda config: [
                            I(config.Land_Cover.LAI, as_='LAI'),
                    ],
                    state_inputs=lambda state: [
                        I(state.canopy_layer_component[0]
                          [0].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
                    ],
                    state_outputs=lambda result: [
                        (result, 'canopy.LAI_total'),
                        (result, 'canopy_component.0.LAI'),
                    ],
                ),
                LAIMethods.ESTIMATE_TOTAL: Process(
                    func=canopy_structure.LAI_method_estimate_canopy_total,
                    comment="""Use primary land cover's estimate of total LAI and spread over
                        layers and LCs""",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.parameters[0].phenology.LAI_1, as_='LAI_1'),
                        I(config.Land_Cover.parameters[0].phenology.LAI_2, as_='LAI_2'),
                        I(config.Land_Cover.parameters[0].phenology.LAI_a, as_='LAI_a'),
                        I(config.Land_Cover.parameters[0].phenology.LAI_b, as_='LAI_b'),
                        I(config.Land_Cover.parameters[0].phenology.LAI_c, as_='LAI_c'),
                        I(config.Land_Cover.parameters[0].phenology.LAI_1, as_='LAI_1'),
                        I(config.Land_Cover.parameters[0].phenology.LAI_d, as_='LAI_d'),
                        #   TODO: Check that sowing date is SGS
                        I(config.Land_Cover.parameters[0]
                          .phenology.key_dates.sowing, as_='SGS'),
                        I(config.Land_Cover.parameters[0]
                          .phenology.key_dates.harvest, as_='EGS'),
                    ],
                    state_inputs=lambda state: [
                        I(state.temporal.dd, as_='dd'),
                    ],
                    state_outputs=lambda result: [
                        (result, 'canopy.LAI_total'),
                        (result, 'canopy_component.0.LAI'),
                    ],
                ),
                LAIMethods.CARBON: [
                    Process(
                        func=lambda LAIs: sum(LAIs),
                        comment="set canopy LAI to total plant LAI",
                        state_inputs=lambda state: [
                            I([state.canopy_component[iLC].LAI for iLC in range(nLC)],
                                as_="LAIs")
                        ],
                        state_outputs=lambda result: [
                            (result, 'canopy.LAI_total'),
                        ],
                    ),
                    Process(
                        func=lambda LAIs: sum(LAIs),
                        comment="set canopy brown leaf LAI to total plant LAI",
                        state_inputs=lambda state: [
                            I([state.canopy_component[iLC].LAI_brown for iLC in range(nLC)],
                                as_="LAIs")
                        ],
                        state_outputs=lambda result: [
                            (result, 'canopy.LAI_brown_total'),
                        ],
                    ),
                ],
            }
        ),

    ]


def distribute_lai_per_layer_processes(
    nL: int,
    nLC: int,
    LAI_distribution_method: LayerLAIDistributionMethods,
    LAI_method: LAIMethods,
):
    return [
        switch(
            gate=LAI_distribution_method,
            comment="Distribute the LAI between layers",
            options={
                LayerLAIDistributionMethods.SKIP: Process(
                    func=skip,
                    comment="LAI distribution disabled",
                ),
                LayerLAIDistributionMethods.FRACTION: Process(
                    func=lambda LAI, fLAIs: [[LAI * fLAI for fLAI in fLAIc] for fLAIc in fLAIs],
                    comment="Distribute LAI using constant fractions from config",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.fLAI, as_="fLAIs"),
                    ],
                    state_inputs=lambda state: [
                        I(state.canopy.LAI_total, as_="LAI"),
                    ],
                    state_outputs=lambda result: [
                        (result[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.LAI')
                        for iL in range(nL) for iLC in range(nLC)
                    ],
                ),
                LayerLAIDistributionMethods.MAX_LAI_PER_LAYER: Process(
                    func=canopy_structure.distribute_lai_per_layers,
                    comment="Distribute lai from bottom layer up with max per layer",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.max_lai_per_layer, as_="max_lai_per_layer"),
                        I(config.Land_Cover.nL, as_="nL"),
                    ],
                    state_inputs=lambda state: [
                        I(state.canopy.LAI_total, as_='total_lai'),
                    ],
                    state_outputs=lambda result: [
                        (result[iL], f'canopy_layer_component.{iL}.{iLC}.LAI')
                        for iL in range(nL) for iLC in range(nLC)
                    ],
                ),
            },
        ),
        switch(
            gate=LAI_distribution_method,
            comment="Distribute the LAI brown between layers",
            options={
                LayerLAIDistributionMethods.SKIP: Process(
                    func=skip,
                    comment="LAI distribution disabled",
                ),
                LayerLAIDistributionMethods.FRACTION: Process(
                    func=lambda LAI, fLAIs: [[LAI * fLAI for fLAI in fLAIc] for fLAIc in fLAIs],
                    comment="Distribute LAI using constant fractions from config",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.fLAI, as_="fLAIs"),
                    ],
                    state_inputs=lambda state: [
                        I(state.canopy.LAI_brown_total, as_="LAI"),
                    ],
                    state_outputs=lambda result: [
                        (result[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.LAI_brown')
                        for iL in range(nL) for iLC in range(nLC)
                    ],
                ),
                LayerLAIDistributionMethods.MAX_LAI_PER_LAYER: Process(
                    func=canopy_structure.distribute_lai_per_layers,
                    comment="Distribute lai from bottom layer up with max per layer",
                    config_inputs=lambda config: [
                        I(config.Land_Cover.max_lai_per_layer, as_="max_lai_per_layer"),
                        I(config.Land_Cover.nL, as_="nL"),
                    ],
                    state_inputs=lambda state: [
                        I(state.canopy.LAI_brown_total, as_='total_lai'),
                    ],
                    state_outputs=lambda result: [
                        (result[iL], f'canopy_layer_component.{iL}.{iLC}.LAI_brown')
                        for iL in range(nL) for iLC in range(nLC)
                    ],
                ),
            },
        ) if LAI_method == LAIMethods.CARBON else [],  # Only run brown leaf calcs for carbon runs
        [Process(
            func=canopy_structure.calc_distribution_of_LAI_between_lcs,
            comment="Calculate the distribution of LAI between land covers",
            gate=LAI_distribution_method != LayerLAIDistributionMethods.SKIP,
            config_inputs=lambda config: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.nLC, as_='nLC'),
            ],
            state_inputs=lambda state: [
                I([[state.canopy_layer_component[jL][jLC].LAI for jL in range(nL)]
                for jLC in range(nLC)], as_='LAI_values')
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result[iLC], f'canopy_component.{iLC}.LC_dist')
            ]
        ) for iLC in range(nLC)],
        Process(
            # sum(this%MC(:)%LC_dist * this%LCs(:)%Lm)
            func=lambda MCs, LCPs: sum([(LC_dist * Lm) for Lm in LCPs for LC_dist in MCs]),
            gate=LAI_distribution_method != LayerLAIDistributionMethods.SKIP,
            comment="get layer LAI-weighted mean leaf width(Lm_LAI)",
            config_inputs=lambda config: [
                I([config.Land_Cover.parameters[iLC].Lm for iLC in range(nLC)], as_='LCPs'),
            ],
            state_inputs=lambda state: [
                I([state.canopy_component[iLC].LC_dist for iLC in range(nLC)], as_='MCs'),
            ],
            state_outputs=lambda result: [(result, 'canopy.Lm_LAI')],
        ),
    ]


def distribute_lai_to_leaf_populations_processes(
    nL: int,
    nP: int,
    nLC: int,
    LAI_method: LAIMethods,
):
    """Distribute lai between leaf populations.

    Parameters
    ----------

    nL: int
        Number of layers
    nP: int
        Number of leaf populations
    nLC: int
        Number of plant components
    """
    def get_multi_layer_leaf_pop_fractions(leaf_pop_distribution, nL, nP):
        # Fraction of each population in layer
        fLAI_layer = [[leaf_pop_distribution[iL][iP] / sum(leaf_pop_distribution[iL]) if sum(
            leaf_pop_distribution[iL]) > 0 else 0 for iP in range(nP)] for iL in range(nL)]

        # fraction of each layer that makes up population
        fLAI_pop = [[leaf_pop_distribution[iL][iP] / sum([leaf_pop_distribution[iL][iP] for iL in range(nL)]) if sum(
            [leaf_pop_distribution[iL][iP] for iL in range(nL)]) > 0 else 0 for iL in range(nL)] for iP in range(nP)]
        return fLAI_layer, fLAI_pop

    return [
        [
            # only using single population.
            Process(
                func=set_value,
                comment="Set leaf lai for single population",
                gate=not nP > 1,
                state_inputs=lambda state, iLC=iLC: [
                    I([state.canopy_layer_component[iL][iLC].LAI for iL in range(nL)], as_="layers_lai"),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    # NOTE: Need to reshape layer lai here to (nL, nP)
                    ([[result["layers_lai"][iL] for iP in range(nP)]
                      for iL in range(nL)], f"canopy_component.{iLC}.leaf_pop_distribution"),
                    (sum(result["layers_lai"]), f'canopy_component_population.{iLC}.0.LAI'),
                    # ([1 for _ in range(nL)],
                    #  f'canopy_component_population.{iLC}.0.fLAI_layer'),
                ],
            ) for iLC in range(nLC)],
        [
            # Using multiple populations and Carbon Allocation
            Process(
                func=canopy_structure.calc_leaf_pops_per_layer,
                comment="Distribute lai between leaf populations",
                gate=nP > 1 and LAI_method == LAIMethods.CARBON,
                config_inputs=lambda config: [
                    I(config.Land_Cover.nL, as_="nL"),
                    I(config.Land_Cover.nP, as_="nP"),
                ],
                state_inputs=lambda state, iLC=iLC: [
                    I(state.prev_hour.canopy_component[iLC].leaf_pop_distribution,
                      as_="prev_leaf_pops_per_lai"),
                    I(state.canopy_component[iLC].LAI, as_="canopy_lai"),
                    I([state.canopy_layer_component[iL][iLC].LAI for iL in range(nL)], as_="layers_lai"),
                    I(state.canopy_component[iLC].growing_populations, as_="growing_populations")
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f"canopy_component.{iLC}.leaf_pop_distribution"),
                    *[(sum([result[iL][iP] for iL in range(nL)]),
                       f"canopy_component_population.{iLC}.{iP}.LAI") for iP in range(nP)],
                ],
            ) for iLC in range(nLC)
        ],
        [
            # Multiple populations not using Carbon Allocation
            # If not using CARBON LAI then distribute evenly
            Process(
                func=canopy_structure.distribute_canopy_lai_to_leaf_pops,
                comment="Distribute lai between growing leaf populations using multiplicative method",
                gate=nP > 1 and LAI_method != LAIMethods.CARBON,
                config_inputs=lambda config: [
                    I(config.Land_Cover.nL, as_="nL"),
                    I(config.Land_Cover.nP, as_="nP"),
                ],
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_component[iLC].LAI, as_="canopy_lai"),
                    I(state.canopy_component[iLC].total_emerged_leaf_populations,
                      as_='no_emerged_pops')
                ],
                state_outputs=lambda result, iLC=iLC: [
                    (result, f"canopy_component.{iLC}.leaf_pop_distribution"),
                    *[(sum([result[iL][iP] for iL in range(nL)]),
                       f"canopy_component_population.{iLC}.{iP}.LAI") for iP in range(nP)],
                ],
            ) for iLC in range(nLC)
        ],
        [
            Process(
                func=get_multi_layer_leaf_pop_fractions,
                comment="Convert lai distribution to pop fractions",
                config_inputs=lambda config: [
                    I(config.Land_Cover.nL, as_="nL"),
                    I(config.Land_Cover.nP, as_="nP"),
                ],
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_component[iLC].leaf_pop_distribution,
                      as_="leaf_pop_distribution"),
                ],
                state_outputs=lambda result, iLC=iLC: [
                    *[(result[0][iL], f"canopy_layer_component.{iL}.{iLC}.fLAI_layer") for iL in range(nL)],
                    *[(result[1][iP], f"canopy_component_population.{iLC}.{iP}.fLAI_layer") for iP in range(nP)],
                ],
            ) for iLC in range(nLC)
        ]
    ]


def calc_canopy_SAI(nL: int, nLC: int, SAI_method: str, primary_SAI_method: str) -> List[Process]:
    """Calculate the SAI throughout the canopy.

    Supplying LAI and SAI for multi-layer multi-component models can be
    awkward.  To generate a nL*nLC array for LAI/SAI, three schemes are
    possible:

      - "input" - All LAI/SAI values are known and supplied as input.
      - "input total" - Total LAI/SAI is known and supplied as input in the top-left cell
        (e.g. LAI(1,1)).  Values are divided among layers and land covers
        according to a (normalised) nL*nLC fLAI array.
      - "estimate total" - Total LAI/SAI is estimated from the properties of the "primary land
        cover" and divided according to fLAI.
    """
    return [
        switch(
            gate=SAI_method,
            comment="calc_canopy_SAI",
            options={
                "input": Process(
                    func=skip,
                    comment="Do nothing if using input",
                    # TODO: Should get from ext data
                ),
                "constant": [Process(
                    func=lambda SAI, fLAI: SAI * fLAI,
                    comment="Set SAI to constant value",
                    config_inputs=lambda config, iL=iL, iLC=iLC: [
                        I(config.Land_Cover.fLAI[iL][iLC], as_='fLAI'),
                        I(config.Land_Cover.SAI, as_='SAI'),
                    ],
                    state_outputs=lambda result, iL=iL, iLC=iLC: [
                        (result, f'canopy_layer_component.{iL}.{iLC}.SAI'),
                    ],
                ) for iL in range(nL) for iLC in range(nLC)],
                "input total": [Process(
                    func=lambda SAI, fLAI: SAI * fLAI,
                    comment="Spread single SAI value to layers and LCs",
                    config_inputs=lambda config, iL=iL, iLC=iLC: [
                        I(config.Land_Cover.fLAI[iL][iLC], as_='fLAI'),
                    ],
                    state_inputs=lambda state: [
                        # TODO: Use primary land cover instead of 0.0
                        I(state.canopy_layer_component[0][0].SAI, as_='SAI'),
                    ],
                    state_outputs=lambda result, iL=iL, iLC=iLC: [
                        (result, f'canopy_layer_component.{iL}.{iLC}.SAI')
                    ]
                ) for iL in range(nL) for iLC in range(nLC)],
                "estimate total": [
                    switch(
                        gate=primary_SAI_method,
                        comment="Use primary land cover's estimate of total SAI",
                        options={
                            SAILCMethods.LAI: Process(
                                func=set_value,
                                comment="Calc_SAI Estimate_total - LAI",
                                state_inputs=lambda state: [
                                    I(state.canopy.LAI_total, as_='LAI_total'),
                                ],
                                state_outputs=lambda result: [
                                    (result['LAI_total'], 'canopy.SAI_total'),
                                ],
                            ),
                            SAILCMethods.LAI_BROWN_GREEN: Process(
                                func=set_value,
                                comment="SAI = sum(LAI and LAI_brown)",
                                state_inputs=lambda state: [
                                    I(state.canopy.LAI_total, as_='LAI_total'),
                                    I(state.canopy.LAI_brown_total, as_='LAI_brown_total'),
                                ],
                                state_outputs=lambda result: [
                                    (sum([result['LAI_total'], result['LAI_brown_total']]),
                                     'canopy.SAI_total'),
                                ],
                            ),
                            SAILCMethods.LAI_max: Process(
                                func=lambda LAI_total, SAI_total_prev: max(
                                    LAI_total or 0, SAI_total_prev or 0),
                                comment="Calc_SAI Estimate_total - max of prev SAI or current LAI",
                                state_inputs=lambda state: [
                                    I(state.canopy.LAI_total, as_='LAI_total'),
                                    I(state.canopy.SAI_total, as_='SAI_total_prev'),
                                ],
                                state_outputs=lambda result: [
                                    (result, 'canopy.SAI_total'),
                                ],
                            ),
                            SAILCMethods.FOREST: Process(
                                func=set_value,
                                comment="Calc_SAI Estimate_total - LAI",
                                state_inputs=lambda state: [
                                    I(state.canopy.LAI_total + 1, as_='LAI_total'),
                                ],
                                state_outputs=lambda result: [
                                    (result['LAI_total'], 'canopy.SAI_total'),
                                ],
                            ),
                            "wheat": NotImplementedError(),
                            # "wheat": Process(
                            #     func=canopy_structure.SAI_wheat_and_LAI,
                            #     comment="Calc_SAI_Estimate_total - wheat",
                            #     config_inputs=lambda config: [
                            #         I(config.Land_Cover.parameters[0].phenology.LAI_1, as_='LAI_1'),
                            #         I(config.Land_Cover.parameters[0].phenology.LAI_2, as_='LAI_2'),
                            #         I(config.Land_Cover.parameters[0].phenology.LAI_a, as_='LAI_a'),
                            #         I(config.Land_Cover.parameters[0].phenology.LAI_b, as_='LAI_b'),
                            #         I(config.Land_Cover.parameters[0].phenology.LAI_c, as_='LAI_c'),
                            #         I(config.Land_Cover.parameters[0].phenology.LAI_1, as_='LAI_1'),
                            #         I(config.Land_Cover.parameters[0].phenology.LAI_d, as_='LAI_d'),
                            #     ],
                            #     state_inputs=lambda state: [
                            #         I(state.temporal.dd, as_='dd'),
                            #         I(state.canopy_component[0].SGS, as_='SGS'),
                            #         I(state.canopy_component[0].EGS, as_='EGS'),
                            #     ],
                            #     state_outputs=lambda result: [
                            #         (result, 'canopy.SAI_total'),
                            #     ]
                            # ),
                        }
                    ),
                    Process(
                        # this%MLMC(:,:)%SAI = this%MLMC(1,1)%SAI * this%fLAI(:,:)
                        func=lambda fLAI, SAI, nL, nLC: [
                            [SAI * fLAI[iL][iLC] for iLC in range(nLC)] for iL in range(nL)],
                        comment="Spread single SAI value to layers and LCs",
                        config_inputs=lambda config: [
                            I(config.Land_Cover.fLAI, as_='fLAI'),
                            I(config.Land_Cover.nL, as_='nL'),
                            I(config.Land_Cover.nLC, as_='nLC'),
                        ],
                        state_inputs=lambda state: [
                            I(state.canopy.SAI_total, as_='SAI'),
                        ],
                        state_outputs=lambda result: [
                            (result[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.SAI')
                            for iL in range(nL) \
                            for iLC in range(nLC)
                        ]
                    ),
                ],
            },
        ),
        Process(
            func=lambda SAI_values: sum([SAI for LC in SAI_values for SAI in LC]),
            comment="get layer SAI_total",
            state_inputs=lambda state: [
                    I([[state.canopy_layer_component[iL][iLC].SAI for iL in range(nL)]
                       for iLC in range(nLC)], as_='SAI_values')
            ],
            state_outputs=lambda result: [(result, 'canopy.SAI_total')],
        ),
    ]


def calc_PAR_sun_shade_process(parsun_shade_method: ParSunShadeMethods, iL: int, iLC: int, nL: int, nLC: int) -> Process:
    return switch(
        gate=parsun_shade_method,
        options={
            ParSunShadeMethods.FARQUHAR1997: Process(
                func=met_irrad_helpers.calc_PAR_sun_shade_farq_b,
                comment="calculate PAR_sun_shade using Farquhar 1997 method",
                config_inputs=lambda config: [
                    # Average cosA
                    I(sum([config.Land_Cover.parameters[jLC].cosA for jLC in range(nLC)]) / nLC,
                      as_='cosA')
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.sinB, row_index), as_='sinB'),
                    I(lget(e_state.Idrctt, row_index), as_='Ir_beam_0'),
                    I(lget(e_state.Idfuse, row_index), as_='Ir_dfuse_0'),
                ],
                state_inputs=lambda state, iL=iL: [
                    # We calculate cumulative Lai here
                    I(sum([state.canopy_layer_component[jL][jLC].SAI for jLC in range(nLC)
                           for jL in range(iL + 1, nL)]), as_='LAI_c'),
                ],
                state_outputs=lambda result, iL=iL: [
                    (result.PARsun, f'canopy_layers.{iL}.micro_met.PARsun'),
                    (result.PARshade, f'canopy_layers.{iL}.micro_met.PARshade'),
                ],
            ),
            ParSunShadeMethods.DO3SEUI: Process(
                func=met_irrad_helpers.calc_PAR_sun_shade_UI,
                comment="Calc par sun shade using DO3SE UI method",
                config_inputs=lambda config, iLC=iLC: [
                    I(config.Land_Cover.parameters[iLC].cosA, as_='cosA'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.sinB, row_index), as_='sinB'),
                    I(lget(e_state.PAR, row_index), as_='PAR'),
                    I(lget(e_state.P, row_index), as_='P'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.LAI_total, as_='LAI'),
                ],

                state_outputs=lambda result, iL=iL: [
                    (result.PARsun, f'canopy_layers.{iL}.micro_met.PARsun'),
                    (result.PARshade, f'canopy_layers.{iL}.micro_met.PARshade'),
                ],

            ),
            ParSunShadeMethods.SIMPLE: Process(
                func=met_irrad_helpers.get_parSunShade_from_par,
                comment="Calc par sun shade using simple method",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].cosA, as_='cosA'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.sinB, row_index), as_='sinB'),
                    I(lget(e_state.PAR, row_index), as_='PAR'),
                    I(lget(e_state.P, row_index), as_='P'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SAI_total, as_='LAI'),
                ],

                state_outputs=lambda result: [
                    (result.PARsun, f'canopy_layers.{iL}.micro_met.PARsun'),
                    (result.PARshade, f'canopy_layers.{iL}.micro_met.PARshade'),
                ],

            ),
        }
    )


def calc_ustar_ref_process(have_Hd_Data: bool, have_ustar_ref_data: bool, have_ustar_data: bool, is_OTC: bool) -> Process:
    return [
        Process(
            func=set_value,
            comment="Set ustar and ustar_ref to input ustar",
            gate=have_ustar_data,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.ustar[row_index], as_='ustar_in'),
            ],
            state_inputs=lambda state: [
                I(state.met.u_i, as_='u'),
            ],
            state_outputs=lambda result:[
                (result['ustar_in'], 'met.ustar_ref'),
                (result['ustar_in'], 'met.ustar'),
            ],
        ),
        Process(
            func=met_wind_helpers.calc_ustar_and_L,
            comment="Calculate ustar_ref and monin obukhov length",
            gate=have_Hd_Data and have_ustar_ref_data,
            config_inputs=lambda config: [
                I(config.Location.z_u, as_='z_u'),
                I(config.Location.u_d, as_='u_d'),
                I(config.Location.u_z0, as_='u_z0'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(lget(e_state.Hd, row_index), as_='Hd'),
                I(lget(e_state.P, row_index), as_='P'),
                I(e_state.Ts_C[row_index] + T0, as_='Tk'),
                I(e_state.ustar_ref[row_index], as_='ustar_ref_in'),
            ],
            state_inputs=lambda state: [
                I(state.met.u_i, as_='u'),
            ],
            state_outputs=lambda result:[
                (result[0], 'met.ustar_ref'),
                (result[1], 'met.L'),
            ],
        ),
        Process(
            func=met_wind_helpers.calc_ustar_and_L,
            comment="Calculate ustar_ref and monin obukhov length",
            gate=have_Hd_Data and not have_ustar_ref_data,
            config_inputs=lambda config: [
                I(config.Location.z_u, as_='z_u'),
                I(config.Location.u_d, as_='u_d'),
                I(config.Location.u_z0, as_='u_z0'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(lget(e_state.Hd, row_index), as_='Hd'),
                I(lget(e_state.P, row_index), as_='P'),
                I(e_state.Ts_C[row_index] + T0, as_='Tk'),
            ],
            state_inputs=lambda state: [
                I(state.met.u_i, as_='u'),
                I(state.met.ustar_ref, as_='ustar_ref_in'),
            ],
            state_outputs=lambda result:[
                (result[0], 'met.ustar_ref'),
                (result[1], 'met.L'),
            ],
        ),
        Process(
            func=set_value,
            comment="Set ustar_ref to external input",
            gate=not have_Hd_Data and have_ustar_ref_data,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.ustar_ref[row_index], as_='ustar_ref_in'),
            ],
            state_outputs=lambda result:[
                (result['ustar_ref_in'], 'met.ustar_ref'),
                (None, 'met.L'),
            ],
        ),
        Process(
            func=met_wind_helpers.ustar_from_velocity_simple,
            comment="Calculate ustar_ref using simple method",
            gate=not have_Hd_Data and not have_ustar_ref_data and not is_OTC,
            config_inputs=lambda config: [
                I(config.Location.z_u - config.Location.u_d, as_='z'),
                I(config.Location.u_z0, as_='z0'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(lget(e_state.u, row_index), as_='u'),
            ],
            state_outputs=lambda result:[
                (result, 'met.ustar_ref'),
                (None, 'met.L'),
            ],
        ),
        Process(
            func=met_wind_helpers.ustar_from_velocity_simple,
            comment="Calculate ustar_ref using simple method and OTC",
            gate=not have_Hd_Data and not have_ustar_ref_data and is_OTC,
            state_inputs=lambda state: [
                I(state.canopy.canopy_height - state.canopy.d, as_="z"),
                I(state.canopy.z0, as_="z0"),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(lget(e_state.u, row_index), as_='u'),
            ],
            state_outputs=lambda result:[
                (result, 'met.ustar_ref'),
                (None, 'met.L'),
            ],
        ),
    ]


def calc_windspeed_parameters_process(using_ustar_input_data: bool) -> Process:
    """Calculate estimated windspeed at canopy."""
    if using_ustar_input_data:
        return Process(
            func=met_wind_helpers.calc_windspeed_ustar_in,
            comment="Calculate estimated windspeed at canopy with ustar input",
            config_inputs=lambda config: [
                I(config.Location.OTC, as_='o_top_chamber'),
                I(config.Location.z_u, as_='z_u'),
                I(config.Location.u_z0, as_='u_z0'),
                I(config.Location.u_d, as_='u_d'),
                I(config.Location.izr, as_='izr'),
            ],
            state_inputs=lambda state: [
                I(state.canopy.canopy_height, as_='h'),
                I(state.canopy.d, as_='d'),
                I(state.canopy.z0, as_='z0'),
                I(state.met.L, as_='L'),
                I(state.met.ustar, as_='ustar_in'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(lget(e_state.u, row_index), as_="u"),
            ],
            additional_inputs=lambda: [
                # TODO: Get min windspeed from config
                I(0.01, as_='min_windspeed'),
            ],
            state_outputs=lambda result: [
                (result.u_i, 'met.u_i'),
                (result.ustar, 'met.ustar'),
                (result.micro_u, 'canopy.met.micro_u'),  # micro_u at top of canopy
            ],
        )
    else:
        return Process(
            func=met_wind_helpers.calc_windspeed,
            comment="Calculate estimated windspeed at canopy",
            config_inputs=lambda config: [
                I(config.Location.OTC, as_='o_top_chamber'),
                I(config.Location.z_u, as_='z_u'),
                I(config.Location.u_z0, as_='u_z0'),
                I(config.Location.u_d, as_='u_d'),
                I(config.Location.izr, as_='izr'),
            ],
            state_inputs=lambda state: [
                I(state.canopy.canopy_height, as_='h'),
                I(state.canopy.d, as_='d'),
                I(state.canopy.z0, as_='z0'),
                I(state.met.L, as_='L'),
                I(state.met.ustar_ref, as_='ustar_ref'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(lget(e_state.u, row_index), as_="u"),
            ],
            additional_inputs=lambda: [
                # TODO: Get min windspeed from config
                I(0.01, as_='min_windspeed'),
            ],
            state_outputs=lambda result: [
                (result.u_i, 'met.u_i'),
                (result.ustar, 'met.ustar'),
                (result.micro_u, 'canopy.met.micro_u'),  # micro_u at top of canopy
            ],
        )


def calc_layer_windspeed_process(iL: int) -> Process:
    """Set constant wind speed at each layer."""
    return Process(
        func=met_wind_helpers.calc_layer_windspeed,
        comment="Set constant wind speed at each layer",
        config_inputs=lambda config: [
            I(config.Location.OTC, as_='o_top_chamber'),
        ],
        state_inputs=lambda state, iL=iL: [
            I(state.canopy.canopy_height, as_='h'),
            I(state.canopy.Lm_LAI, as_='w'),
            # TODO: Why is SAI different
            I(state.canopy.SAI_total, as_='SAI'),
            I(state.canopy.met.micro_u, as_='u_at_canopy_top'),
            # TODO: Check layer height not layer depth
            I(state.canopy_layers[iL].layer_height, as_='z'),
            # TODO: Calculate layer depth
            # I(state.canopy_layers[iL].layer_depth, as_='layer_depth'),
            I(0, as_='layer_depth'),
        ],
        state_outputs=lambda result, iL=iL: [
            (result, f'canopy_layers.{iL}.micro_met.micro_u'),
        ]
    )


def MLMC_sunlit_LAI_process(nL: int, nLC: int) -> Process:
    """Calculate the sunlit and shaded LAI fractions."""
    # NOTE: Reversed layers
    return Process(
        func=met_irrad_helpers.MLMC_sunlit_LAI,
        comment="Estimate sunlit LAI fractions",
        config_inputs=lambda config: [
            I(config.Land_Cover.nL, as_='nL'),
            I(config.Land_Cover.nLC, as_='nLC'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(lget(e_state.sinB, row_index), as_="sinB"),
        ],
        state_inputs=lambda state: [
            I([[state.canopy_layer_component[iL][iLC].SAI for iLC in range(nLC)]
               for iL in reversed(range(nL))], as_='LAI'),
        ],
        state_outputs=lambda result: [
            (list(reversed(result))[iL][iLC], f'canopy_layer_component.{iL}.{iLC}.LAIsunfrac')
            for iL in range(nL)
            for iLC in range(nLC)
        ],
    )


def calc_soil_moisture_hourly_process(soil_moisture_source: str) -> List[Process]:
    """Calculate soil moisture changes during hour.

    Sources:
    - "disabled" - Not set
    - "input SWP" -
    - "input SWC" -
    - "P-M" -
    """
    return switch(
        gate=soil_moisture_source,
        comment="Calculate SMD values from whichever soil moisture input is being used.",
        options={
            SoilMoistureSource.DISABLED: Process(
                func=skip,
                comment="Disabled - Nothing to do"
            ),
            SoilMoistureSource.INPUT_SWP: Process(
                func=SMD_helpers.soil_moisture_from_SWP,
                comment="Calculate SMD values from SWP",
                group="soil-moisture",
                config_inputs=lambda config: [
                    I(config.soil_moisture.soil_config, as_='soil_config'),
                    I(config.soil_moisture.PWP, as_='PWP'),
                    I(config.soil_moisture.root, as_='root_depth'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.SWP, as_='SWP'),
                ],
                state_outputs=lambda result: [
                    (result.Sn, 'canopy.SMD.Sn'),
                    (result.SWP, 'canopy.SMD.SWP'),
                    (result.ASW, 'canopy.SMD.ASW'),
                    (result.SMD, 'canopy.SMD.SMD'),
                ]
            ),
            SoilMoistureSource.INPUT_SWC_EXTERNAL: Process(
                func=SMD_helpers.soil_moisture_from_SWC,
                comment="Calculate SMD values from external SWC",
                group="soil-moisture",
                config_inputs=lambda config: [
                    I(config.soil_moisture.soil_config, as_='soil_config'),
                    I(config.soil_moisture.PWP, as_='PWP'),
                    I(config.soil_moisture.root, as_='root_depth'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.SWC, row_index), as_='Sn_in'),
                ],
                state_outputs=lambda result: [
                    (result.Sn, 'canopy.SMD.Sn'),
                    (result.SWP, 'canopy.SMD.SWP'),
                    (result.ASW, 'canopy.SMD.ASW'),
                    (result.SMD, 'canopy.SMD.SMD'),
                ]
            ),
            SoilMoistureSource.INPUT_ASW_EXTERNAL: Process(
                func=set_value,
                comment="Calculate SMD values from external SWC",
                group="soil-moisture",
                config_inputs=lambda config: [
                    I(config.soil_moisture.soil_config, as_='soil_config'),
                    I(config.soil_moisture.PWP, as_='PWP'),
                    I(config.soil_moisture.root, as_='root_depth'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.ASW, row_index), as_='ASW'),
                ],
                state_outputs=lambda result: [
                    (result['ASW'], 'canopy.SMD.ASW'),
                ]
            ),
            SoilMoistureSource.INPUT_SWC: Process(
                func=SMD_helpers.soil_moisture_from_SWC,
                comment="Calculate SMD values from calculated SWC",
                group="soil-moisture",
                config_inputs=lambda config: [
                    I(config.soil_moisture.soil_config, as_='soil_config'),
                    I(config.soil_moisture.PWP, as_='PWP'),
                    I(config.soil_moisture.root, as_='root_depth'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.Sn, as_='Sn_in'),
                ],
                state_outputs=lambda result: [
                    (result.Sn, 'canopy.SMD.Sn'),
                    (result.SWP, 'canopy.SMD.SWP'),
                    (result.ASW, 'canopy.SMD.ASW'),
                    (result.SMD, 'canopy.SMD.SMD'),
                ]
            ),
            SoilMoistureSource.P_M: Process(
                func=skip,
                group="soil-moisture",
                gate=False,
                comment="PM Soil moisture calculated daily",
            ),
        }
    )


def f_phen_method_process(iL: int, iLC: int, f_phen_method: str) -> Process:
    """Calculate f_phen using specified method.

    Method options:
    - "simple day PLF" -
    - "complex day PLF" -
    - "tt day PLF" -
    """
    # TODO: Can reduce this now phenology pre calculated
    return switch(
        gate=f_phen_method,
        comment="select a f_phen method",
        options={
            # TODO: Add input
            FPhenMethods.DISABLED: Process(
                func=skip,
                comment="f_phen supplied (or left at default value of 1.0)",
                state_outputs=lambda result: [
                    (1.0,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
                ]
            ),
            FPhenMethods.INPUT: Process(
                func=set_value,
                comment="f_phen being set from input data",
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.leaf_fphen[row_index], as_='fphen'),
                ],
                state_outputs=lambda result: [
                    (result['fphen'],
                     f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
                ],
            ),
            FPhenMethods.SIMPLE_DAY_PLF: Process(
                func=f_phen_helpers.f_phen_simple_PLF,
                comment="f_phen_method - f_phen_simple_PLF",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.day_fphen_plf.f_phen_1, as_='f_phen_1'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.day_fphen_plf.f_phen_4, as_='f_phen_4'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.day_fphen_plf.f_phen_a, as_='f_phen_a'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.day_fphen_plf.f_phen_c, as_='f_phen_c'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.day_fphen_plf.f_phen_e, as_='f_phen_e'),
                    #   TODO: Check that sowing date is SGS
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_dates.sowing, as_='SGS'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_dates.harvest, as_='EGS'),
                ],
                state_inputs=lambda state: [
                    I(state.temporal.dd, as_='dd'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
                ],
            ),
            # TODO: Not validated with UI
            FPhenMethods.COMPLEX_DAY_PLF: NOT_IMPLEMENTED_PROCESS('Use td method'),
            # FPhenMethods.COMPLEX_DAY_PLF: Process(
            #     func=f_phen_helpers.f_phen_complex_PLF,
            #     comment="f_phen_method - f_phen_complex_PLF",
            #     config_inputs=lambda config: [
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_1, as_='f_phen_1'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_2, as_='f_phen_2'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_3, as_='f_phen_3'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_4, as_='f_phen_4'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_a, as_='f_phen_a'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_b, as_='f_phen_b'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_c, as_='f_phen_c'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_d, as_='f_phen_d'),
            #         I(config.Land_Cover.parameters[iLC]
            #           .phenology.day_fphen_plf.f_phen_e, as_='f_phen_e'),
            #         I(config.Land_Cover.parameters[iLC].phenology.day_fphen_plf.f_phen_limA,
            #           as_='f_phen_limA'),
            #         I(config.Land_Cover.parameters[iLC].phenology.day_fphen_plf.f_phen_limB,
            #           as_='f_phen_limB'),
            #     ],
            #     state_inputs=lambda state: [
            #         I(state.temporal.dd, as_='dd'),
            #         I(state.canopy_component[iLC].SGS, as_='SGS'),
            #         I(state.canopy_component[0].EGS, as_='EGS'),
            #     ],
            #     state_outputs=lambda result: [
            #         (result,
            #             f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
            #     ]
            # ),
            FPhenMethods.TT_DAY_PLF: Process(
                func=f_phen_helpers.tt_f_phen_simple_PLF,
                comment="f_phen_method - tt_f_phen_simple_PLF",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.sowing_to_emerge, as_='t_f_phen_a'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.sowing_to_f_phen_b,
                      as_='t_f_phen_b'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.sowing_to_f_phen_c,
                      as_='t_f_phen_c'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.sowing_to_end,
                      as_='t_f_phen_d'),
                    I(config.Land_Cover.parameters[iLC].phenology.f_phen_min,
                      as_='f_phen_min'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_dates_td.sowing, as_="td_at_sgs")
                ],
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_component[iLC].td_v, as_='td'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
                ],
            ),
            # DEPRECATEDas calculated in setup
            # FPhenMethods.TT_GROWING_SEASON: Process(
            #     func=get_current_f_phen_from_t_sgs_t_egs,
            #     comment="f_phen_method - tt growing season",
            #     config_inputs=lambda config: [
            #         I(0, as_='t_sgs'),
            #         I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.sowing_to_emerge, as_='t_egs'),
            #         # TODO: Add fraction parameters
            #     ],
            #     external_state_inputs=lambda e_state, row_index: [
            #         # TODO: Check td is 0 at sgs
            #         I(lget(e_state.td, row_index), as_='td'),
            #     ],
            #     state_outputs=lambda result: [
            #         (result,
            #             f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_phen'),
            #     ],
            # ),
        }
    )


def leaf_f_phen_process(iL: int, iLC: int, leaf_f_phen_method) -> Process:
    """Calculate leaf_f_phen using specified method.

    Method options:
    - "disabled" - Not set
    - "f_phen" -
    - "day PLF" -
    - "leaf tt day PLF" -
    """
    # TODO: Can reduce this now phenology pre calculated

    return switch(
        gate=leaf_f_phen_method,
        comment="Choose leaf_f_phen_method",
        options={
            LeafFPhenMethods.DISABLED: Process(
                func=skip,
                comment="calculate leaf_f_phen_method - disabled",
            ),
            LeafFPhenMethods.INPUT: Process(
                func=set_value,
                comment="leaf_f_phen being set from input data",
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.leaf_fphen[row_index], as_='leaf_fphen'),
                ],
                state_outputs=lambda result: [
                    (result['leaf_fphen'],
                     f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
                ],
            ),
            LeafFPhenMethods.F_PHEN: Process(
                func=set_value,
                comment="calculate leaf_f_phen_method - f_phen",
                state_inputs=lambda state: [
                    I(state.canopy_layer_component[iL][iLC].gsto_params.f_phen,
                        as_='leaf_f_phen'),
                ],
                state_outputs=lambda result: [
                    (result['leaf_f_phen'],
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
                ]
            ),
            LeafFPhenMethods.DAY_PLF: Process(
                func=f_phen_helpers.leaf_f_phen_PLF,
                comment="calculate leaf_f_phen_method - day PLF",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].phenology.day_fphen_plf.leaf_f_phen_1,
                        as_='leaf_f_phen_1'),
                    I(config.Land_Cover.parameters[iLC].phenology.day_fphen_plf.leaf_f_phen_2,
                        as_='leaf_f_phen_2'),
                    I(config.Land_Cover.parameters[iLC].phenology.day_fphen_plf.leaf_f_phen_a,
                        as_='leaf_f_phen_a'),
                    I(config.Land_Cover.parameters[iLC].phenology.day_fphen_plf.leaf_f_phen_b,
                        as_='leaf_f_phen_b'),
                    I(config.Land_Cover.parameters[iLC].phenology.day_fphen_plf.leaf_f_phen_c,
                        as_='leaf_f_phen_c'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_dates.Astart, as_='Astart'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_dates.Aend, as_='Aend'),
                ],
                state_inputs=lambda state: [
                    I(state.temporal.dd, as_='dd'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
                ],
            ),
            LeafFPhenMethods.TT_DAY_PLF: Process(
                # TODO: td is time since zero day
                func=f_phen_helpers.tt_leaf_f_phen_PLF,
                comment="calculate leaf_f_phen_method - leaf tt day PLF",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].phenology.leaf_f_phen_a,
                      as_='t_leaf_f_phen_a'),
                    I(config.Land_Cover.parameters[iLC].phenology.leaf_f_phen_b,
                      as_='t_leaf_f_phen_b'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.leaf_f_phen_e,
                      as_='t_leaf_f_phen_e'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.leaf_f_phen_g,
                      as_='t_leaf_f_phen_g'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.leaf_f_phen_h,
                      as_='t_leaf_f_phen_h'),
                    I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.leaf_f_phen_i,
                      as_='t_leaf_f_phen_i'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_lengths_td.sowing_to_astart, as_='t_astart'),
                    I(config.Land_Cover.parameters[iLC]
                      .phenology.key_dates_td.sowing, as_="td_at_sgs")
                ],
                state_inputs=lambda state, iLC=iLC: [
                    I(state.canopy_component[iLC].td_v, as_='td'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
                ],
            ),
            # LeafFPhenMethods.TT_GROWING_SEASON: Process(
            #     func=get_current_leaf_f_phen_from_t_sgs_t_egs,
            #     comment="leaf_f_phen_method - tt growing season",
            #     config_inputs=lambda config: [
            #         I(0, as_='t_sgs'),
            #         I(config.Land_Cover.parameters[iLC].phenology.key_lengths_td.sowing_to_emerge, as_='t_egs'),
            #         # TODO: Add fraction parameters
            #     ],
            #     external_state_inputs=lambda e_state, row_index: [
            #         # TODO: Check td is 0 at sgs
            #         I(lget(e_state.td, row_index), as_='td'),
            #     ],
            #     state_outputs=lambda result: [
            #         (result,
            #             f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_phen'),
            #     ],
            # ),
        }
    )


def calc_f_light_process(iL: int, iLC: int) -> Process:
    """Calculate f_light."""
    return Process(
        func=gsto_helpers.calc_f_light_method,
        # gate=f_light_method == "enabled",
        comment="Calculate f_light",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].multip_gsto.f_lightfac, as_='f_lightfac'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(lget(e_state.sinB, row_index), as_='sinB'),
            I(lget(e_state.PAR, row_index), as_='PAR'),
        ],
        state_inputs=lambda state, iLC=iLC, iL=iL: [
            I(state.canopy.LAI_total, as_='LAI'),
            I(state.canopy_layers[iL].micro_met.PARsun, as_='PARsun'),
            I(state.canopy_layers[iL].micro_met.PARshade, as_='PARshade'),
            I(state.canopy_layer_component[iL][iLC].LAIsunfrac, as_='LAIsunfrac'),
        ],
        state_outputs=lambda result, iLC=iLC, iL=iL: [
            (result.f_light, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_light'),
            (result.leaf_f_light,
             f'canopy_layer_component.{iL}.{iLC}.gsto_params.leaf_f_light'),
        ],
    )


def calc_f_temp_process(iL: int, iLC: int, f_temp_method: str) -> Process:
    """Calculate f_temp using specified method.

    Method options:
    - "disabled" - Not set
    - "default" -
    - "square high" -
    """
    return switch(
        gate=f_temp_method,
        comment="Choose f_temp Method",
        options={
            "disabled": Process(
                func=skip,
                comment="f_temp - disabled",
            ),
            "default": Process(
                func=gsto_helpers.calc_f_temp,
                comment="f_temp - default",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_min, as_='T_min'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_opt, as_='T_opt'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_max, as_='T_max'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.Ts_C[row_index], as_='Ts_C'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_temp'),
                ],
            ),
            "square high": Process(
                func=gsto_helpers.calc_f_temp_square_high,
                comment="f_temp - square high",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_min, as_='T_min'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_opt, as_='T_opt'),
                    I(config.Land_Cover.parameters[iLC].multip_gsto.T_max, as_='T_max'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.Ts_C[row_index], as_='Ts_C'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_temp'),
                ],
            ),
        },
    )


def calc_f_VPD_process(iL: int, iLC: int, f_VPD_method: FVPDMethods) -> Process:
    """Calculate f_VPD using specified method.

    Method options:
     - "disabled" - Not set
     - "linear" -
     - "log" -
     - "leuning" - Calculated inside photosynthesis gsto An model
     - "danielsson" - Calculated inside photosynthesis gsto An model
    """
    return switch(
        gate=f_VPD_method,
        comment="Choose f_VPD method",
        options={
            FVPDMethods.DISABLED: Process(
                func=skip,
                comment="f_VPD method - disabled",
            ),
            FVPDMethods.LINEAR: Process(
                func=gsto_helpers.f_VPD_linear,
                comment="f_VPD_method - linear",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.VPD_max, as_='VPD_max'),
                    I(config.Land_Cover.parameters[iLC].gsto.VPD_min, as_='VPD_min'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.VPD, row_index), as_='VPD'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_VPD'),
                ]
            ),
            FVPDMethods.LOG: Process(
                func=gsto_helpers.f_VPD_log,
                comment="f_VPD_method - log",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.VPD, row_index), as_='VPD'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_VPD'),
                ]
            ),
            FVPDMethods.LEUNING: Process(
                # f_VPD calculated during gsto calculation
                func=skip,
                comment="f_VPD method - photosynthesis",
            ),

            FVPDMethods.DANIELSSON: Process(
                # f_VPD calculated during gsto calculation
                func=skip,
                comment="f_VPD method - Danielsson",
            ),
        },
    )


def f_SW_process(iL: int, iLC: int, f_SW_method: str) -> Process:
    """Calculate f_sw using specified method.

    Method options:
     - "disabled" - Not set
     - "fSWP exp" -
     - "fSWP linear" -
     - "fPAW" -
    """
    return switch(
        gate=f_SW_method,
        comment="Choose f_SW method",
        options={
            "disabled": Process(
                func=set_value,
                comment="f_SW_method - Disabled",
                state_outputs=lambda result: [
                    (1.0, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
            "fSWP exp": Process(
                func=gsto_helpers.f_SWP_exp,
                comment="f_SW_method - fSWP exp",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.fSWP_exp_a, as_='a'),
                    I(config.Land_Cover.parameters[iLC].gsto.fSWP_exp_b, as_='b'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.SWP, as_='SWP'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
            "fSWP linear": Process(
                func=gsto_helpers.f_SWP_linear,
                comment="f_SW_method - linear",
                config_inputs=lambda config: [
                    I(config.Land_Cover.parameters[iLC].gsto.SWP_min, as_='SWP_min'),
                    I(config.Land_Cover.parameters[iLC].gsto.SWP_max, as_='SWP_max'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.SWP, as_='SWP'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
            "fPAW": Process(
                func=gsto_helpers.f_PAW,
                comment="f_SW_method - fPAW",
                config_inputs=lambda config: [
                    I(config.soil_moisture.ASW_FC, as_='ASW_FC'),
                    I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
                    I(config.Land_Cover.parameters[iLC].gsto.ASW_min, as_='ASW_min'),
                    I(config.Land_Cover.parameters[iLC].gsto.ASW_max, as_='ASW_max'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.SMD.ASW, as_='ASW'),
                ],
                state_outputs=lambda result: [
                    (result,
                        f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_SW'),
                ],
            ),
        },
    )


def f_O3_process(iL: int, iLC: int, f_O3_method: str, nP: int = 0) -> Process:
    """Calculate f_O3 using specified method.

    NOTE: Only valid for single layer model

    Method options:
     - "disabled" - Not set
     - "wheat" - Using wheat calc
     - "potato" - Using Potato calc
     - "Steph wheat" - Using Steph Wheat calc

    """
    if f_O3_method != "disabled":
        if iL > 0 or nP > 1:
            raise ConfigError("Cannot use fO3 with multi layer or multi population models")
    iP = 0
    return switch(
        gate=f_O3_method,
        comment="Choose f_O3 method",
        options={
            "disabled": Process(
                func=skip,
                comment="f_O3 method - disabled",
            ),
            "wheat": Process(
                func=lambda POD_0: ((1 + (POD_0 / 11.5)**10)**(-1)),
                comment="f_O3_method - wheat",
                state_inputs=lambda state, iLC=iLC, iP=iP: [
                    I(state.canopy_component_population[iLC][iP].POD_0, as_='POD_0'),
                ],
                state_outputs=lambda result: [
                    (result, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_O3'),
                ],
            ),
            "potato": Process(
                func=lambda AOT_0: ((1 + (AOT_0 / 40)**5)**(-1)),
                comment="f_O3_method - potato",
                state_inputs=lambda state, iLC=iLC, iP=iP: [
                    I(state.canopy_component_population[iLC][iP].AOT_0, as_='AOT_0'),
                ],
                state_outputs=lambda result: [
                    (result, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_O3'),
                ],
            ),
            "Steph wheat": Process(
                func=lambda POD_0: ((1 + (POD_0 / 27)**10)**(-1)),
                comment="f_O3_method - Steph wheat",
                state_inputs=lambda state, iLC=iLC, iP=iP: [
                    I(state.canopy_component_population[iLC][iP].POD_0, as_='POD_0'),
                ],
                state_outputs=lambda result: [
                    (result, f'canopy_layer_component.{iL}.{iLC}.gsto_params.f_O3'),
                ],
            ),
        },
    )


def set_vcmax_and_jmax(nL: int, iLC: int, iP: int, V_J_method: str) -> Process:
    """Set V_cmax_25 and J_max_25 using the specified method.

    Method Options:
     - "constant" - "V_cmax_25 and J_max_25 supplied as constant.
     - "input" - "V_cmax_25 and J_max_25 supplied as hourly inputs.

    """
    return [switch(
        gate=V_J_method,
        comment="Choose f_O3 method",
        options={
            "constant": Process(
                func=set_value,
                comment="V_cmax_25 and J_max_25 supplied as constant",
                config_inputs=lambda config, iLC=iLC: [
                    I(config.Land_Cover.parameters[iLC].pn_gsto.V_cmax_25, as_='V_cmax_25'),
                    I(config.Land_Cover.parameters[iLC].pn_gsto.J_max_25, as_='J_max_25'),
                ],
                state_outputs=lambda result, iLC=iLC, iP=iP: [
                    (result['V_cmax_25'],
                     f'canopy_component_population.{iLC}.{iP}.V_cmax_25_per_layer.{0}'),
                    (result['J_max_25'], f'canopy_component_population.{iLC}.{iP}.J_max_25'),
                ],
            ),
            "input": Process(
                func=set_value,
                comment="V_cmax_25 and J_max_25 supplied as hourly inputs",
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.V_cmax_25[row_index], as_='V_cmax_25'),
                    I(e_state.J_max_25[row_index], as_='J_max_25'),
                ],
                state_outputs=lambda result, iLC=iLC, iP=iP: [
                    (result['V_cmax_25'],
                     f'canopy_component_population.{iLC}.{iP}.V_cmax_25_per_layer.{0}'),
                    (result['J_max_25'], f'canopy_component_population.{iLC}.{iP}.J_max_25'),
                ],
            )
        },
    ),
        Process(
        func=multilayer_vcmax25,
        gate=nL > 1,
        comment="Scale v_cmax_25 to layers below canopy",
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.parameters[iLC].pn_gsto.V_cmax_25_kN, as_='kN'),
        ],
        state_inputs=lambda state, iLC=iLC, iP=iP: [
            I([state.canopy_layer_component[iL][iLC].LAI for iL in range(nL)], as_="layer_LAI"),
            I(state.canopy_component_population[iLC][iP].V_cmax_25_per_layer[0], as_='V_cmax_25'),
        ],
        state_outputs=lambda result, iLC=iLC, iP=iP: [
            (result, f'canopy_component_population.{iLC}.{iP}.V_cmax_25_per_layer'),
        ],
    )]


def calc_leaf_f_phen_effect_on_V_cmax_25_process(
    iL: int,
    iLC: int,
    iP: int,
    leaf_f_phen_Anet_influence: str,
) -> Process:
    """Calculate gsto phenology effect on V_cmax_25 and J_max_25."""
    return Process(
        func=f_phen_helpers.calc_leaf_f_phen_effect_on_V_cmax_25,
        comment="Calculate gsto phenology effect on V_cmax_25 and J_max_25",
        gate=leaf_f_phen_Anet_influence == LeafFPhenAnetInfluence.V_C_MAX,
        state_inputs=lambda state, iL=iL, iLC=iLC, iP=iP: [
            I(state.canopy_layer_component[0][iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
            I(state.canopy_component_population[iLC][iP].V_cmax_25[iL], as_='V_cmax_25_in'),
            I(state.canopy_component_population[iLC][iP].J_max_25, as_='J_max_25_in'),
        ],
        state_outputs=lambda result, iL=iL, iLC=iLC, iP=iP: [
            (result.V_cmax_25, f'canopy_component_population.{iLC}.{iP}.V_cmax_25_per_layer.{iL}'),
            (result.J_max_25, f'canopy_component_population.{iLC}.{iP}.J_max_25'),
        ],
    )


def calc_leaf_temperature_process(
    iL: int,
    iLC: int,
    iP: int,
    t_leaf_method: TLeafMethods,
) -> Process:
    return switch(
        gate=t_leaf_method,
        # TODO: May need to be per leaf population
        comment="Choose t_leaf method",
        options={
            TLeafMethods.INPUT: Process(
                func=set_value,
                comment="Use input t_leaf_temperature",
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.Tleaf_C[row_index], as_='Tleaf_C'),
                ],
                state_outputs=lambda result, iL=iL, iLC=iLC, iP=iP: [
                    (result['Tleaf_C'],
                     f'canopy_layer_component_pop.{iL}.{iLC}.{iP}.Tleaf_C_estimate'),
                ],
            ),
            TLeafMethods.AMBIENT: Process(
                func=set_value,
                comment="Use air temperature",
                external_state_inputs=lambda e_state, row_index: [
                    I(e_state.Ts_C[row_index], as_='Tleaf_C'),
                ],
                state_outputs=lambda result, iL=iL, iLC=iLC, iP=iP: [
                    (result['Tleaf_C'],
                     f'canopy_layer_component_pop.{iL}.{iLC}.{iP}.Tleaf_C_estimate'),
                ],
            ),
            TLeafMethods.DE_BOECK: Process(
                func=get_leaf_temp_de_boeck,
                comment="Use de Boeck method",
                config_inputs=lambda config: [
                    I(config.Location.albedo, as_='albedo'),
                    I(True, as_='hypostomatous'),
                    I(config.Land_Cover.parameters[iLC].Lm, as_='d'),
                ],
                state_inputs=lambda state, iL=iL, iLC=iLC, iP=iP: [
                    I((state.canopy_component_population[iLC]
                       [iP].mean_gsto_per_layer[iL] or 0) * DRATIO_H20_O3, as_='g_vs'),
                    I(state.canopy_layers[iL].micro_met.micro_u, as_='u_speed'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.R, row_index), as_='R'),
                    I(lget(e_state.eact, row_index), as_='eact'),
                    I(e_state.Ts_C[row_index], as_='T_air'),
                    I(e_state.Ts_C[row_index], as_='initial_T_leaf'),
                    I(lget(e_state.P, row_index), as_='P'),
                    # TODO: If cloudfrac is None then use 1.0
                    I(lget(e_state.cloudfrac, row_index, 1.0), as_='cloud_cover'),
                    # I(e_state.cloudfrac and lget(e_state.cloudfrac, row_index) or 1.0, as_='cloud_cover'),
                ],
                state_outputs=lambda result, iL=iL, iLC=iLC, iP=iP: [
                    # TODO: Should be per leaf pop
                    (result, f'canopy_layer_component_pop.{iL}.{iLC}.{iP}.Tleaf_C_estimate'),
                ],
            ),
        },
    )


def gsto_multiplicative_process(iL: int, iLC: int, iP: int) -> Process:
    """Calculate Gsto using the multiplicative method."""
    return Process(
        # TODO: Could seperate this into leaf and canopy calculations
        func=multiplicative,
        comment="Calculate gsto - multiplicative",
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.parameters[iLC].gsto.VPD_crit, as_='VPD_crit'),
            I(config.Land_Cover.parameters[iLC].multip_gsto.gmax, as_='gmax'),
            I(config.Land_Cover.parameters[iLC].multip_gsto.gmorph, as_='gmorph'),
            I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.VPD_dd[row_index], as_='VPD_dd')
        ],
        state_inputs=lambda state, iL=iL, iLC=iLC, iP=iP: [
            # TODO: Make leaf_gsto multilayer
            I(state.canopy_component_population[iLC]
              [iP].mean_gsto_per_layer[iL], as_='initial_leaf_gsto'),
            I(state.canopy_layer_component[iL][iLC].mean_gsto, as_='initial_mean_gsto'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.f_phen, as_='f_phen'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.f_light, as_='f_light'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.leaf_f_light, as_='leaf_f_light'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.f_temp, as_='f_temp'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.f_VPD, as_='f_VPD'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.f_SW, as_='f_SW'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.f_O3, as_='f_O3'),
        ],
        state_outputs=lambda result, iL=iL, iLC=iLC, iP=iP: [
            # TODO: Make leaf_gsto multilayer
            # TODO: Storing gsto in multiple places!
            (result.new_leaf_gsto,
             f'canopy_component_population.{iLC}.{iP}.mean_gsto_per_layer.{iL}'),
            #  TODO: mean_gsto_sunlit may not be valid for multiplicative
            (result.new_leaf_gsto,
             f'canopy_component_population.{iLC}.{iP}.mean_gsto_sunlit'),
            # NOTE: Multiplicative mean gsto is for canopy not leaf population
            (result.new_mean_gsto, f'canopy_layer_component.{iL}.{iLC}.mean_gsto'),
        ]
    )


def calc_D_0_process(iLC: int) -> Process:
    """Calculate the D_0 for photosynthesis.

    # TODO: This is a constant so could be calculated at start of model
    """
    return Process(
        func=pn_helpers.calc_D_0,
        comment="Calculate the D_0 for photosynthesis",
        # TODO: Use a switch below for D_0 methods
        config_inputs=lambda config, iLC=iLC: [
            # TODO: Check correct approach for photosynthesis approach
            I(config.Land_Cover.parameters[iLC].pn_gsto.D_0_method, as_='D_0_method'),
            I(config.Land_Cover.parameters[iLC].gsto.f_VPD_method, as_='f_VPD_method'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.D_0, as_='constant_D_0'),
            I(config.Land_Cover.parameters[iLC].gsto.VPD_max, as_='VPD_max'),
            I(config.Land_Cover.parameters[iLC].gsto.VPD_min, as_='VPD_min'),
            I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.D_0'),
        ]
    )


def get_vernalisation_factor_process(
    iLC: int,
) -> Process:
    return [
        Process(
            func=vernalisation_helpers.calculate_vernalisation_factor,
            comment="Calculate vernalisation factor",
            config_inputs=lambda config: [
                I(config.Land_Cover.parameters[iLC].phenology.v_T_max, as_="v_T_max"),
                I(config.Land_Cover.parameters[iLC].phenology.v_T_min, as_="v_T_min"),
                I(config.Land_Cover.parameters[iLC].phenology.PIV, as_="PIV"),
            ],
            external_state_inputs=lambda e_state, row_index: [
                # TODO: Should get prev day here?
                I(max(e_state.Ts_C[row_index:row_index + 24]), as_='max_ambient_temp'),
                I(min(e_state.Ts_C[row_index:row_index + 24]), as_='min_ambient_temp'),
            ],
            state_inputs=lambda state: [
                I(state.canopy_component[iLC].phenology.phenology_stage, as_="phenology_stage"),
                I(state.prev_hour.canopy_component[iLC].phenology.V_acc, as_="V_acc_prev"),
            ],
            state_outputs=lambda result: [
                (result[0], f'canopy_component.{iLC}.phenology.Vf'),
                (result[1], f'canopy_component.{iLC}.phenology.V_acc'),
                (result[2], f'canopy_component.{iLC}.phenology.V_pos'),
                (result[3], f'canopy_component.{iLC}.phenology.V_neg'),
                (result[2] + result[3], f'canopy_component.{iLC}.phenology.V_tot'),
            ],
        ),
    ]


def get_td_dd_process(iLC: int) -> Process:
    """Get the difference between current td and td at plant emergence.

    Also takes into account Vernalisation.

    """

    def calc_td_dd(td_prev, td, accumulate, td_dd_prev, Vf, ppf):
        # if td_prev is None assume first row so ignore.
        f = Vf * ppf if Vf < 1 else 1
        # NOTE: Alt method uses min instead of *
        # f = min(Vf, ppf) if Vf < 1 else 1
        return td_dd_prev + f * (td - td_prev) if accumulate and td_prev is not None else td_dd_prev

    return [
        Process(
            func=calc_td_dd,
            comment="Get the difference between current td and td at Sowing(Vernalised if on)",
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].td, as_="td"),
                I(state.prev_hour.canopy_component[iLC].td, as_="td_prev"),
                I(state.canopy_component[iLC].phenology.phenology_stage >=
                  PhenologyStage.SOWN, as_='accumulate'),
                I(state.canopy_component[iLC].phenology.Vf, as_='Vf'),
                I(state.canopy_component[iLC].photoperiod_factor, as_='ppf'),
                I(state.prev_hour.canopy_component[iLC].td_v, as_='td_dd_prev'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.td_v'),
            ],
        ),
        Process(
            func=calc_td_dd,
            comment="Get the difference between current td and td at emergence",
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].td, as_="td"),
                I(state.prev_hour.canopy_component[iLC].td, as_="td_prev"),
                I(state.canopy_component[iLC].phenology.phenology_stage >=
                  PhenologyStage.EMERGED, as_='accumulate'),
                I(state.canopy_component[iLC].phenology.Vf, as_='Vf'),
                I(state.canopy_component[iLC].photoperiod_factor, as_='ppf'),
                I(state.prev_hour.canopy_component[iLC].td_dd, as_='td_dd_prev'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.td_dd'),
            ],
        ),
    ]


def calc_td_dd_per_leaf_pop_process(iLC: int) -> Process:
    return Process(
        func=phyllochron_dvi.calc_td_dd_per_leaf_pop,
        comment="Calculate the thermal time difference between leaf pop emergence and current td",
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.nP, as_='nP'),
            I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.plant_emerg_to_leaf_emerg,
              as_='t_emerg_to_flag_emerg'),
        ],
        state_inputs=lambda state, iLC=iLC: [
            # TODO: Flag leaf may have already emerged?
            I(state.canopy_component[iLC].td_dd, as_='td'),
            I(state.prev_hour.canopy_component[iLC].td_dd, as_='td_prev'),
            I(state.prev_hour.canopy_component[iLC].td_dd_leaf_pops, as_='td_dd_prev'),
            I(state.canopy_component[iLC].total_emerged_leaf_populations, as_='emerged_leaf_count'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.td_dd_leaf_pops')
        ],
    )


def calc_g_bv_process(iLC: int) -> Process:
    """Calculate g_bv_H2O for photosynthesis."""
    return Process(
        func=pn_helpers.calc_g_bv,
        comment="Calculate g_bv for photosynthesis",
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.parameters[iLC].Lm, as_='Lm'),
        ],
        # TODO: Should this be per layer?
        external_state_inputs=lambda e_state, row_index: [
            I(lget(e_state.u, row_index), as_='u'),
        ],
        additional_inputs=lambda: [
            I(GAS.H2O, as_='gas'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.g_bv_H2O'),
        ],
    )


def ozone_damage_processes(
    iLC: int,
    iP: int,
    flag_leaf: bool,
    senescence_method: SenescenceFunctionMethods,
) -> List[Process]:
    return [
        Process(
            func=calc_all_ozone_damage_factors,
            comment="Calculate the per population ozone damage Ewert Method",
            gate=senescence_method == SenescenceFunctionMethods.EWERT,
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_1, as_='gamma_1'),
                I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_2, as_='gamma_2'),
                I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_3, as_='gamma_3'),
                I(config.Land_Cover.parameters[iLC].pn_gsto.cL3, as_='cL3'),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_em
                    if flag_leaf else config.Land_Cover.parameters[iLC].phenology.key_lengths_leaf_td.tl_em,
                    as_='t_lem'),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_se
                    if flag_leaf else config.Land_Cover.parameters[iLC].phenology.key_lengths_leaf_td.tl_se,
                    as_='t_lse'),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_ep
                    if flag_leaf else config.Land_Cover.parameters[iLC].phenology.key_lengths_leaf_td.tl_ep,
                    as_='t_lep'),
                I(config.Land_Cover.parameters[iLC].pn_gsto.use_O3_damage, as_='use_O3_damage'),
                I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_4_senes, as_='gamma_4_senes'),
                I(config.Land_Cover.parameters[iLC].pn_gsto.gamma_5_harvest, as_='gamma_5_harvest'),
                # TODO: Link these to global config
                I(False, as_='opt_full_night_recovery'),
            ],
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I(state.canopy_component[iLC].td_dd_leaf_pops[iP], as_='td_dd'),
                I(state.canopy_component_population[iLC][iP].O3up_acc,
                  as_='O3up_acc'),
                I(state.canopy_component_population[iLC][iP].O3up, as_='O3up'),
                I(state.prev_hour.canopy_component_population[iLC][iP].fO3_d, as_='fO3_d_prev'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.is_daylight[row_index], as_='is_daylight'),
                I(lget(e_state.hr, row_index), as_='hr'),
            ],
            state_outputs=lambda result, iP=iP, iLC=iLC: [
                (result.fO3_h, f'canopy_component_population.{iLC}.{iP}.fO3_h'),
                (result.fO3_d, f'canopy_component_population.{iLC}.{iP}.fO3_d'),
                (result.fO3_l, f'canopy_component_population.{iLC}.{iP}.fO3_l'),
                (result.f_LA, f'canopy_component_population.{iLC}.{iP}.f_LA'),
                (result.f_LS, f'canopy_component_population.{iLC}.{iP}.f_LS'),
                (result.rO3, f'canopy_component_population.{iLC}.{iP}.rO3'),
                (result.t_lep_limited, f'canopy_component_population.{iLC}.{iP}.t_lep_limited'),
                (result.t_lse_limited, f'canopy_component_population.{iLC}.{iP}.t_lse_limited'),
            ],
        ),
        Process(
            func=set_value,
            comment="Override f_LS and fO3_d value(Anet run)",
            gate=senescence_method == SenescenceFunctionMethods.ANET,
            config_inputs=lambda config: [
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_se
                  if flag_leaf else config.Land_Cover.parameters[iLC].phenology.key_lengths_leaf_td.tl_se,
                  as_='t_lse'),
                I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.tl_ep
                    if flag_leaf else config.Land_Cover.parameters[iLC].phenology.key_lengths_leaf_td.tl_ep,
                    as_='t_lep'),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_layer_component[0][iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
                I(state.canopy_layer_component[0][iLC].gsto_params.f_O3, as_='f_O3'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result['f_O3'], f'canopy_component_population.{iLC}.{iP}.fO3_d'),
                (result['leaf_f_phen'], f'canopy_component_population.{iLC}.{iP}.f_LS'),
                (1, f'canopy_component_population.{iLC}.{iP}.f_LA'),
                (result['t_lep'], f'canopy_component_population.{iLC}.{iP}.t_lep_limited'),
                (result['t_lse'], f'canopy_component_population.{iLC}.{iP}.t_lse_limited'),
            ],
        ),
    ]


def ewert_leaf_process(iLC: int, iP: int, nL: int, ewert_loop_method: EwertLoopMethods) -> Process:
    """Ewert Photosynthesis model."""
    ewert_leaf_fn = ewert_leaf_pop if ewert_loop_method is EwertLoopMethods.ITERATIVE \
        else ewert_leaf_pop_cubic if ewert_loop_method is EwertLoopMethods.CUBIC else None

    return Process(
        func=ewert_leaf_fn,
        comment="Ewert Photosynthesis model",
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.nL, as_="nL"),
            I(config.Land_Cover.parameters[iLC].pn_gsto.g_sto_0, as_='g_sto_0'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.m, as_='m'),
            I(config.Land_Cover.parameters[iLC].gsto.f_VPD_method, as_='f_VPD_method'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.co2_concentration_balance_threshold,
              as_='co2_concentration_balance_threshold'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.co2_concentration_max_iterations,
              as_='co2_concentration_max_iterations'),
            I(config.Land_Cover.parameters[iLC].pn_gsto.R_d_coeff, as_="R_d_coeff"),
        ],
        state_inputs=lambda state, iLC=iLC, iP=iP: [
            I(state.canopy_component_population[iLC][iP].f_LS, as_='f_LS'),
            I(state.canopy_component_population[iLC][iP].fO3_d, as_='fO3_d'),
            I([state.canopy_component_population[iLC][iP].V_cmax_25_per_layer[iL]
               for iL in range(nL)], as_='V_cmax_25'),
            I([state.canopy_component_population[iLC][iP].J_max_25
               for _ in range(nL)], as_='J_max_25'),
            I([state.canopy_layers[iL].micro_met.PARsun for iL in range(nL)], as_='PARsun'),
            I([state.canopy_layers[iL].micro_met.PARshade for iL in range(nL)], as_='PARshade'),
            I([state.canopy_layer_component[iL][iLC].LAIsunfrac for iL in range(nL)], as_='LAIsunfrac'),
            I(state.canopy_component[iLC].D_0, as_='D_0'),
            I(state.canopy_component[iLC].g_bv_H2O, as_='g_bv'),
            I(state.canopy_component_population[iLC][iP].fLAI_layer, as_='layer_lai_frac'),
            I([state.canopy_component[iLC].leaf_pop_distribution[iL][iP]
               for iL in range(nL)], as_='layer_lai'),
            # NOTE: Below should all be same for all layers and populations
            I([state.canopy_layer_component_pop[iL][iLC]
               [iP].Tleaf_C_estimate for iL in range(nL)], as_='Tleaf_C'),
            I(state.canopy_layer_component[0][iLC].gsto_params.f_SW, as_='f_SW'),
            I(state.canopy_layer_component[0][iLC].gsto_params.f_VPD, as_='f_VPD'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(lget(e_state.eact, row_index), as_='eact'),
            I(lget(e_state.CO2, row_index), as_='c_a'),
        ],
        state_outputs=lambda result, iLC=iLC, iP=iP: [
            (result.g_sv_per_layer, f'canopy_component_population.{iLC}.{iP}.g_sv_per_layer'),
            (result.g_sv_sunlit, f'canopy_component_population.{iLC}.{iP}.g_sv_sunlit'),
            # NOTE: Below are already upscaled to actual leaf values
            (result.A_n, f'canopy_component_population.{iLC}.{iP}.A_n'),
            (result.A_n_sunlit, f'canopy_component_population.{iLC}.{iP}.A_n_sunlit'),
            (result.A_c, f'canopy_component_population.{iLC}.{iP}.A_c'),
            (result.A_j, f'canopy_component_population.{iLC}.{iP}.A_j'),
            (result.A_p, f'canopy_component_population.{iLC}.{iP}.A_p'),
            (result.A_n_limit_factor,
                f'canopy_component_population.{iLC}.{iP}.A_n_limit_factor'),
            (result.R_d, f'canopy_component_population.{iLC}.{iP}.R_d'),
            (result.c_i, f'canopy_component_population.{iLC}.{iP}.c_i'),
            # TODO: Output f_VPD
            # (result.f_VPD, f'canopy_component_population.{iLC}.{iP}.gsto_params.f_VPD'),
            (result.v_cmax, f'canopy_component_population.{iLC}.{iP}.V_cmax'),
            (result.j_max, f'canopy_component_population.{iLC}.{iP}.J_max'),
            (result.loop_iterations, 'debug.ewert_loop_iterations'),
        ]
    )


def convert_gsto_CO2umol_to_O3mmol_process(iLC: int, iP: int, nL: int) -> Process:
    """Convert g_sto from umol to mmol, and from CO2 to O3."""
    def gsv_to_gsto(g_sv: float) -> float:
        return DRATIO_O3_CO2 * (max(0.0, g_sv) / 1000)
    return [
        Process(
            func=lambda g_sv: [gsv_to_gsto(g_sv[iL]) for iL in range(nL)],
            comment="Convert g_sto per layer from umol to mmol, and from CO2 to O3",
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I(state.canopy_component_population[iLC][iP].g_sv_per_layer, as_='g_sv'),
            ],
            state_outputs=lambda result: [
                (result, f'canopy_component_population.{iLC}.{iP}.mean_gsto_per_layer'),
            ],
        ), Process(
            func=gsv_to_gsto,
            comment="Convert g_sto sunlit from umol to mmol, and from CO2 to O3",
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I(state.canopy_component_population[iLC][iP].g_sv_sunlit, as_='g_sv'),
            ],
            state_outputs=lambda result: [
                (result, f'canopy_component_population.{iLC}.{iP}.mean_gsto_sunlit'),
            ],
        ),
    ]


def limit_gsto_l_with_leaf_fphen_process(
    iL: int,
    iLC: int,
    iP: int,
    leaf_f_phen_Anet_influence: bool,
) -> Process:
    """Optionally limit the output leaf_gsto using the leaf_f_phen.

    After gsto calculateions.

    """
    if leaf_f_phen_Anet_influence != LeafFPhenAnetInfluence.DISABLED:
        raise ConfigError("Using leaf f phen anet influence is deprecated")
    return Process(
        func=lambda leaf_f_phen, leaf_gsto: leaf_f_phen * leaf_gsto,
        comment="Limit the leaf_gsto with leaf_f_phen",
        # TODO: Use per layer leaf_gsto
        gate=leaf_f_phen_Anet_influence == LeafFPhenAnetInfluence.G_STO,
        state_inputs=lambda state, iL=iL, iLC=iLC, iP=iP: [
            I(state.canopy_component_population[iLC][iP].mean_gsto_per_layer[iL], as_='leaf_gsto'),
            I(state.canopy_layer_component[iL][iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
        ],
        state_outputs=lambda result, iL=iL, iLC=iLC, iP=iP: [
            (result, f'canopy_component_population.{iLC}.{iP}.mean_gsto_per_layer.{iL}'),
        ],
    )


def calc_layer_mean_gsto_process(iLC: int, nP: int, nL: int) -> Process:
    """Set mean gsto to equal leaf gsto using weighted population average."""
    return Process(
        func=calc_mean_gsto,
        comment="set mean gsto",
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.nL, as_="nL"),
        ],
        state_inputs=lambda state, iLC=iLC: [
            I([[state.canopy_component_population[iLC]
                [iP].mean_gsto_per_layer[iL] for iP in range(nP)] for iL in range(nL)], as_='leaf_gsto_list'),
            I([state.canopy_layer_component[iL][iLC].fLAI_layer for iL in range(nL)], as_='leaf_fLAI_list'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result[iL], f'canopy_layer_component.{iL}.{iLC}.mean_gsto') for iL in range(nL)
        ],
    )


def scale_layer_mean_gsto_to_layer_bulk_gsto_process(iL: int, iLC: int) -> Process:
    """Scale mean gsto up to bulk gsto."""
    return Process(
        func=lambda mean_gsto, LAI: mean_gsto * LAI,
        comment="Scale mean gsto up to bulk gsto",
        state_inputs=lambda state, iL=iL, iLC=iLC: [
            I(state.canopy_layer_component[iL][iLC].mean_gsto, as_='mean_gsto'),
            I(state.canopy_layer_component[iL][iLC].LAI, as_='LAI'),
        ],
        state_outputs=lambda result, iL=iL, iLC=iLC: [
            (result, f'canopy_layer_component.{iL}.{iLC}.bulk_gsto')
        ],
    )


def scale_leaf_layer_mean_gsto_to_leaf_layer_bulk_gsto_process(iLC: int, iP: int, nL: int) -> Process:
    """Scale mean gsto up to bulk gsto."""
    return Process(
        func=lambda mean_gsto, LAI: [mean_gsto[iL] * LAI[iL] for iL in range(nL)],
        comment="Scale mean gsto up to bulk gsto",
        state_inputs=lambda state, iLC=iLC, iP=iP: [
            I(state.canopy_component_population[iLC][iP].mean_gsto_per_layer, as_='mean_gsto'),
            I([state.canopy_component[iLC].leaf_pop_distribution[iL][iP]
               for iL in range(nL)], as_='LAI'),
        ],
        state_outputs=lambda result, iLC=iLC, iP=iP: [
            (result, f'canopy_component_population.{iLC}.{iP}.bulk_gsto_per_layer')
        ],
    )


def scale_layer_bulk_gsto_to_canopy_gsto_process(nL: int, iLC: int) -> Process:
    """Scale bulk gsto to canopy gsto."""
    return Process(
        func=lambda bulk_gsto: sum(bulk_gsto),
        comment="Scale bulk gsto to canopy gsto",
        state_inputs=lambda state, iLC=iLC: [
            I([state.canopy_layer_component[iL][iLC].bulk_gsto for iL in range(nL)], as_='bulk_gsto'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.canopy_gsto')
        ],
    )


def scale_rd_to_canopy_rd_process(iLC: int, nP: int) -> Process:
    """Scale population Rd to canopy Rd."""
    return Process(
        # NOTE: Already upscaled to per leaf LAI
        # func=lambda R_d, LAI: sum([a * l for a, l in zip(R_d, LAI)]),
        func=lambda R_d: sum(R_d),
        comment="Scale population Rd to canopy Rd",
        state_inputs=lambda state, iLC=iLC: [
            I([state.canopy_component_population[iLC][iP].R_d for iP in range(nP)], as_='R_d'),
            # I([state.canopy_component_population[iLC][iP].LAI for iP in range(nP)], as_='LAI'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.R_dc')
        ],
    )


def scale_anet_to_canopy_anet_process(iLC: int, nP: int) -> Process:
    """Scale population anet to canopy anet."""
    return Process(
        # NOTE: Already upscaled to per leaf LAI
        # func=lambda An, LAI: sum([a * l for a, l in zip(An, LAI)]),
        func=lambda A_n: sum(A_n),
        comment="Scale layer anet to canopy anet",
        state_inputs=lambda state, iLC=iLC: [
            I([state.canopy_component_population[iLC][iP].A_n for iP in range(nP)], as_='A_n'),
            # I([state.canopy_component_population[iLC][iP].LAI for iP in range(nP)], as_='LAI'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.canopy_anet')
        ],
    )


def scale_anet_to_canopy_an_gross_process(iLC: int, nP: int) -> Process:
    """Scale layer anet to canopy gross anet."""
    return Process(
        # NOTE: Already upscaled to per leaf LAI
        # func=lambda A_n, R_d, LAI: sum([(a + r) * l for a, r, l in zip(A_n, R_d, LAI)]),
        func=lambda A_n, R_d: sum([(a + r) for a, r in zip(A_n, R_d)]),
        comment="Scale layer anet to canopy gross anet",
        state_inputs=lambda state, iLC=iLC: [
            I([state.canopy_component_population[iLC][iP].A_n for iP in range(nP)], as_='A_n'),
            # I([state.canopy_component_population[iLC][iP].LAI for iP in range(nP)], as_='LAI'),
            I([state.canopy_component_population[iLC][iP].R_d for iP in range(nP)], as_='R_d'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.canopy_an_gross')
        ],
    )


def calculate_canopy_npp_process(iLC: int) -> Process:
    """Calculate the net primary productivity(NPP)."""
    return Process(
        func=carbon_calcs.calc_net_prod,
        comment="Calculate the net primary productivity(NPP)",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[iLC].pn_gsto.r_g, as_="r_g"),
        ],
        state_inputs=lambda state, iLC=iLC: [
            I(umol_c_to_kg_c(state.canopy_component[iLC].canopy_anet), as_='canopy_An'),
            I(umol_c_to_kg_c(state.canopy_component[iLC].R_dc), as_='R_dc'),
            I(state.canopy_component[iLC].c_root, as_='c_root'),
            I(state.canopy_component[iLC].c_stem, as_='c_stem'),
            I(state.canopy_component[iLC].c_leaf, as_='c_leaf'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result.NPP, f'canopy_component.{iLC}.npp'),
            (result.GPP, f'canopy_component.{iLC}.gpp'),
            (result.R_pm, f'canopy_component.{iLC}.R_pm'),
            (result.R_pg, f'canopy_component.{iLC}.R_pg'),
        ],
    )


def accumulate_canopy_npp_process(iLC: int) -> Process:
    """Accumulate hourly net primary productivity."""
    return Process(
        func=accumulate,
        comment="Accumulate hourly net primary productivity",
        state_inputs=lambda state, iLC=iLC: [
            I(state.canopy_component[iLC].npp_acc, as_='npp_acc'),
            I(state.canopy_component[iLC].npp, as_='npp'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result, f'canopy_component.{iLC}.npp_acc')
        ],
    )


def calc_resistance_model_process(
    nL: int,
    nLC: int,
    ra_calc_method: str,
    rsur_calc_method: str,
    is_OTC: bool,
) -> Process:
    """Reset and Calculate the resistance model for O3 over the target canopy.

    Parameters
    ----------
    nL : int
        number of model layers
    nLC : int
        number of model components
    ra_calc_method : str
        ra calc method. Options:
        - "simple" -
        - "heat_flux" -

    """
    return [
        Process(
            func=resistance_helpers.calc_resistance_model,
            comment="Reset and Calculate the resistance model for O3 over the reference canopy",
            gate=not is_OTC,
            config_inputs=lambda config: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.nLC, as_='nLC'),
                I(config.Location.Rsoil, as_='Rsoil'),
                I(ra_calc_method, as_='ra_calc_method'),
                I(rsur_calc_method, as_='rsur_calc_method'),
                I(config.Location.h_O3, as_='canopy_height'),
                I(config.Location.z_O3, as_='measured_height'),
                I(config.Location.izr, as_='izr'),
                I(config.Location.canopy_d, as_='CANOPY_D'),
                I(config.Location.canopy_zo, as_='CANOPY_Z0'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.Ts_C[row_index], as_='Ts_C'),
                I(lget(e_state.snow_depth, row_index), as_='snow_depth'),
            ],
            state_inputs=lambda state: [
                I(state.met.ustar_ref, as_='ustar'),
                I(state.met.L, as_='L'),  # TODO: Only required for Ra heat flux
                I([[state.canopy_layer_component[iL][iLC].SAI for iLC in range(nLC)]
                   for iL in range(nL)], as_='SAI_values'),
                I([[state.canopy_layer_component[iL][iLC].LAI for iLC in range(nLC)]
                   for iL in range(nL)], as_='LAI_values'),
                I([[state.canopy_layer_component[iL][iLC].mean_gsto for iLC in range(nLC)]
                   for iL in range(nL)], as_='mean_gsto_values'),
                # I([[state.canopy_layer_component[iL][iLC].bulk_gsto for iLC in range(nLC)]
                #    for iL in range(nL)], as_='bulk_gsto_values'),
            ],
            additional_inputs=lambda:[
                I(DIFF_O3, as_='Rb_diff'),
            ],
            state_outputs=lambda result: [
                (result, 'canopy.rmodel_O3_ref'),
            ],
        ),
        Process(
            func=resistance_helpers.calc_resistance_model,
            comment="Reset and Calculate the resistance model for O3 over the reference canopy OTC",
            gate=is_OTC,
            # TODO: Make sure using reference canopy values here
            config_inputs=lambda config: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.nLC, as_='nLC'),
                I(config.Location.Rsoil, as_='Rsoil'),
                I(ra_calc_method, as_='ra_calc_method'),
                I(rsur_calc_method, as_='rsur_calc_method'),
                I(config.Location.izr, as_='izr'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.Ts_C[row_index], as_='Ts_C'),
                I(lget(e_state.snow_depth, row_index), as_='snow_depth'),
            ],
            state_inputs=lambda state: [
                I(state.canopy.canopy_height, as_='canopy_height'),
                I(state.canopy.canopy_height, as_='measured_height'),
                I(state.canopy.d, as_='CANOPY_D'),
                I(state.canopy.z0, as_='CANOPY_Z0'),
                I(state.met.ustar_ref, as_='ustar'),
                I(state.met.L, as_='L'),  # TODO: Only required for Ra heat flux
                I([[state.canopy_layer_component[iL][iLC].SAI for iLC in range(nLC)]
                   for iL in range(nL)], as_='SAI_values'),
                I([[state.canopy_layer_component[iL][iLC].LAI for iLC in range(nLC)]
                   for iL in range(nL)], as_='LAI_values'),
                I([[state.canopy_layer_component[iL][iLC].mean_gsto for iLC in range(nLC)]
                   for iL in range(nL)], as_='mean_gsto_values'),
                # I([[state.canopy_layer_component[iL][iLC].bulk_gsto for iLC in range(nLC)]
                #    for iL in range(nL)], as_='bulk_gsto_values'),
            ],
            additional_inputs=lambda:[
                I(DIFF_O3, as_='Rb_diff'),
                # Added to stop double counting of CANOPY_D
                I(1, as_="CANOPY_D"),
                I(1, as_="CANOPY_Z0"),
            ],
            state_outputs=lambda result: [
                (result, 'canopy.rmodel_O3_ref'),
            ],
        ),
        Process(
            func=resistance_helpers.calc_resistance_model,
            comment="Reset and Calculate the resistance model for O3 over the target canopy",
            config_inputs=lambda config: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.nLC, as_='nLC'),
                I(config.Location.Rsoil, as_='Rsoil'),
                I(ra_calc_method, as_='ra_calc_method'),
                I(rsur_calc_method, as_='rsur_calc_method'),
                I(config.Location.z_O3, as_='measured_height'),
                I(config.Location.izr, as_='izr'),
                I(config.Location.canopy_d, as_='CANOPY_D'),
                I(config.Location.canopy_zo, as_='CANOPY_Z0'),
            ],
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.Ts_C[row_index], as_='Ts_C'),
                I(lget(e_state.snow_depth, row_index), as_='snow_depth'),
            ],
            state_inputs=lambda state: [
                I(state.met.ustar, as_='ustar'),
                I(state.canopy.canopy_height, as_='canopy_height'),
                I(state.met.L, as_='L'),  # TODO: Only required for Ra heat flux
                I([[state.canopy_layer_component[iL][iLC].SAI for iLC in range(nLC)]
                   for iL in range(nL)], as_='SAI_values'),
                I([[state.canopy_layer_component[iL][iLC].LAI for iLC in range(nLC)]
                   for iL in range(nL)], as_='LAI_values'),
                I([[state.canopy_layer_component[iL][iLC].mean_gsto for iLC in range(nLC)]
                   for iL in range(nL)], as_='mean_gsto_values'),
                # I([[state.canopy_layer_component[iL][iLC].bulk_gsto for iLC in range(nLC)]
                #    for iL in range(nL)], as_='bulk_gsto_values'),
            ],
            additional_inputs=lambda:[
                I(DIFF_O3, as_='Rb_diff'),
            ],
            state_outputs=lambda result: [
                (result, 'canopy.rmodel_O3'),
            ],
        ),
    ]


def calc_deposition_velocity_process() -> Process:
    """Vd calculation duplicated here just for comparison purposes."""
    # TODO: This may now be redundant as calculated as part of ozone deposition
    return Process(
        func=resistance_helpers.calc_deposition_velocity,
        comment="Vd calculation duplicated here just for comparison purposes",
        # ! TODO: better way to expose Vd
        state_inputs=lambda state: [
            # TODO: Fix this
            I(state.canopy.rmodel_O3.Ra_canopy_to_izr, as_='rmodel_Ra_c'),
            I(state.canopy.rmodel_O3.Rtotal[0], as_='rmodel_Rtotal_top_layer'),
        ],
        state_outputs=lambda result: [
            (result, 'canopy.Vd'),
        ],
    )


def calc_O3_otc_process(nL) -> Process:
    """Set O3 values to O3 input in OTC. Uniform across layers."""
    return Process(
        func=set_value,
        comment="Set O3 values to O3 input in OTC. Uniform across layers",
        gate="OTC == True",
        external_state_inputs=lambda e_state, row_index: [
            I(lget(e_state.O3, row_index), as_='O3'),
        ],
        state_outputs=lambda result: [
            (result['O3'], 'met.O3_i'),
            (result['O3'], 'canopy.canopy_top_o3'),
            *[(result['O3'], f'canopy_layers.{iL}.micro_met.micro_O3') for iL in range(nL)],
        ],
    )


def calc_canopy_ozone_concentration_process(
    canopy_ozone_method: str,
    nL: int,
) -> Process:
    """Calculate Top Layer Canopy ozone."""
    top_layer_index = nL - 1
    return switch(
        gate=canopy_ozone_method,
        comment="Calculate canopy_ozone_method",
        options={
            OzoneDepositionMethods.SINGLE_LAYER: Process(
                func=met_deposition_helpers.calc_canopy_ozone_concentration,
                comment="Calculate Top Layer Canopy ozone single layer",
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.O3, row_index), as_='O3_ppb_zR'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.rmodel_O3.Rsur[top_layer_index], as_='Rsur_top_layer'),
                    I(state.canopy.rmodel_O3.Rb, as_='Rb_top_layer'),
                    I(state.canopy.rmodel_O3_ref.Rb, as_='Rb_ref'),
                    I(state.canopy.rmodel_O3_ref.Rsur[-1], as_='Rsur_ref'),
                    I(state.canopy.rmodel_O3.Ra_canopy_to_izr, as_='Ra_tar_canopy'),
                    I(state.canopy.rmodel_O3.Ra_canopy_top_to_izr, as_='Ra_tar_canopy_top'),
                    I(state.canopy.rmodel_O3_ref.Ra_canopy_to_izr, as_='Ra_ref_canopy'),
                    I(state.canopy.rmodel_O3_ref.Ra_measured_to_izr, as_='Ra_ref_measured'),
                ],
                state_outputs=lambda result: [
                    (result.O3_i, 'met.O3_i'),
                    (result.micro_O3, f'canopy.canopy_top_o3'),
                    # NOTE: Below will be overriden
                    (result.micro_O3, f'canopy_layers.{top_layer_index}.micro_met.micro_O3'),
                    (result.Vd, 'canopy.Vd'),
                ]),
            OzoneDepositionMethods.MULTI_LAYER: Process(
                # TODO: update this!!
                func=met_deposition_helpers.calc_canopy_ozone_concentration_alt,
                comment="Calculate Top Layer Canopy ozone multi layer",
                config_inputs=lambda config: [
                    I(config.Location.h_O3, as_='h_O3_in'),
                    I(config.Location.z_O3, as_='z_O3'),
                    I(config.Location.O3_d, as_='O3_d'),
                    I(config.Location.O3_z0, as_='O3_z0'),
                    I(config.resistance.ra_calc_method, as_='ra_method'),
                    I(config.Location.izr, as_='izr'),
                ],
                external_state_inputs=lambda e_state, row_index: [
                    I(lget(e_state.O3, row_index), as_='O3'),
                ],
                state_inputs=lambda state: [
                    I(state.canopy.canopy_height, as_='canopy_height'),
                    I(state.canopy.rmodel_O3.Rtotal[top_layer_index], as_='Rtotal_top_layer'),
                    I(state.met.u_i, as_='u_i'),
                    I(state.met.L, as_='L'),  # TODO: Should be ozone L
                    I(state.met.ustar_ref, as_='ustar_ref'),  # TODO: Should be ozone ustar_ref
                ],
                state_outputs=lambda result: [
                    (result.O3_i, 'met.O3_i'),
                    (result.micro_O3, f'canopy.canopy_top_o3'),
                    # NOTE: Below will be overriden
                    (result.micro_O3, f'canopy_layers.{top_layer_index}.micro_met.micro_O3'),
                ],
            ),
        },
    )


def calc_multi_layer_O3_ozone_concentration_process(nL: int) -> Process:
    """Calculate Ozone concentration for other layers.

    NOTE: We have to reverse layers.
    """
    return Process(
        func=met_deposition_helpers.calc_multi_layer_O3_ozone_concentration,
        comment="Calculate Ozone concentration for other layers",
        config_inputs=lambda config: [
            I(config.Land_Cover.nL, as_='nL'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.canopy_top_o3, as_='O3_in'),
            I(0, as_='rm_Ra'),  # 0 because we have already brought ozone to top of canopy
            # TODO: Check using correct resistances
            # I(state.canopy.rmodel_O3.Ra_canopy_to_izr, as_='rm_Ra'),
            I(list(reversed(state.canopy.rmodel_O3.Rinc)), as_='rm_Rinc'),
            I(list(reversed(state.canopy.rmodel_O3.Rsur)), as_='rm_Rsur'),
            I(state.canopy.rmodel_O3.Rgs, as_='rm_Rgs'),
        ],
        state_outputs=lambda result: [
            (list(reversed(result))[iL], f'canopy_layers.{iL}.micro_met.micro_O3')
            for iL in range(nL)
        ],
    )
    # return Process(
    #     func=set_value,
    #     comment="Calculate Ozone concentration for other layers",
    #     state_inputs=lambda state: [
    #         I(state.canopy_layers[top_layer_index].micro_met.micro_O3, as_='O3_in'),
    #     ],
    #     state_outputs=lambda result: [
    #         (result['O3_in'], f'canopy_layers.{iL}.micro_met.micro_O3')
    #         for iL in range(nL)
    #     ],
    # )


def O3_ppb_to_nmol_process(iL: int) -> Process:
    """Convert O3 ppb to O3 nmol."""
    return Process(
        func=O3_helpers.O3_ppb_to_nmol,
        comment="Convert O3 ppb to O3 nmol",
        state_inputs=lambda state, iL=iL: [
            I(state.canopy_layers[iL].micro_met.micro_O3, as_='O3_ppb'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.Ts_C[row_index], as_='Ts_C'),
            I(lget(e_state.P, row_index), as_='P'),
        ],
        state_outputs=lambda result, iL=iL: [
            (result, f'canopy_layers.{iL}.micro_met.micro_O3_nmol_m3')
        ],
    )


def calc_leaf_resistance_model_process(iLC: int, iP: int, nL: int) -> Process:
    """Reset and setup Leaf Ozone Resistance Model."""
    return [
        Process(
            func=resistance_helpers.calc_leaf_resistance_model,
            comment="Reset and setup Leaf Resistance Model",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.parameters[iLC].Lm, as_='Lm'),
            ],
            additional_inputs=lambda: [
                I(model_constants.LEAF_G_O3, as_="LEAF_G"),
            ],
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I([state.canopy_layers[iL].micro_met.micro_u for iL in range(nL)],
                  as_='u_per_layer'),
                I(state.canopy_component_population[iLC]
                  [iP].mean_gsto_per_layer, as_='leaf_gsto_per_layer'),
                #   [iP].bulk_gsto_per_layer, as_='leaf_gsto_per_layer'),
            ],
            state_outputs=lambda result, iLC=iLC, iP=iP: [
                (result, f'canopy_component_population.{iLC}.{iP}.leaf_rmodel_O3'),
            ],
        ),
        Process(
            func=resistance_helpers.calc_leaf_resistance_model,
            comment="Reset and setup Leaf Resistance Model sunlit",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.nL, as_='nL'),
                I(config.Land_Cover.parameters[iLC].Lm, as_='Lm'),
            ],
            additional_inputs=lambda: [
                I(model_constants.LEAF_G_O3, as_="LEAF_G"),
            ],
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I([state.canopy_layers[iL].micro_met.micro_u for iL in range(nL)],
                  as_='u_per_layer'),
                I([state.canopy_component_population[iLC]
                   [iP].mean_gsto_sunlit for _ in range(nL)], as_='leaf_gsto_per_layer'),
            ],
            state_outputs=lambda result, iLC=iLC, iP=iP: [
                (result, f'canopy_component_population.{iLC}.{iP}.leaf_sunlit_rmodel_O3'),
            ],
        )]


def calc_fst_leaf_process(iLC: int, iP: int, nL: int):
    """Calculate the O3uptake of each leaf population at each level and set average for leaf."""
    # TODO: Can we split this into shaded and sunlit
    return [
        Process(
            func=O3_helpers.calc_fst_leaf,
            comment=f"Calculate the O3up (fst) {iLC}_{iP}",
            config_inputs=lambda config: [
                I(config.Land_Cover.nL, as_='nL'),
            ],
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I([state.canopy_layers[iL].micro_met.micro_O3_nmol_m3 for iL in range(nL)], as_='O3_nmol_m3'),
                I(state.canopy_component_population[iLC][iP].leaf_rmodel_O3.Rb, as_='Rb_l'),
                I(state.canopy_component_population[iLC][iP].leaf_rmodel_O3.Rsto, as_='Rsto_l'),
                I(state.canopy_component_population[iLC][iP].leaf_rmodel_O3.Rext, as_='Rext'),
                I(state.canopy_component_population[iLC][iP].mean_gsto_per_layer, as_='Gsto_l'),
                # TODO: Should this be m^2 or upscaled with LAI?
                I(state.canopy_component_population[iLC][iP].fLAI_layer, as_='fLAI'),
                # I([state.canopy_component[iLC].leaf_pop_distribution[iL][iP]
                #    for iL in range(nL)], as_='LAI_per_layer'),
            ],
            state_outputs=lambda result, iLC=iLC, iP=iP: [
                (result, f'canopy_component_population.{iLC}.{iP}.O3up'),
            ],
        ),
        Process(
            func=O3_helpers.calc_fst_leaf,
            comment=f"Calculate the O3up (fst) for sunlit component {iLC}_{iP}",
            # NOTE: We only use top layer for O3up sunlit
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I([state.canopy_layers[nL - 1].micro_met.micro_O3_nmol_m3], as_='O3_nmol_m3'),
                I([state.canopy_component_population[iLC][iP].leaf_sunlit_rmodel_O3.Rb[nL - 1]], as_='Rb_l'),
                I([state.canopy_component_population[iLC]
                  [iP].leaf_sunlit_rmodel_O3.Rsto[nL - 1]], as_='Rsto_l'),
                I([state.canopy_component_population[iLC]
                  [iP].leaf_sunlit_rmodel_O3.Rext[nL - 1]], as_='Rext'),
                I([state.canopy_component_population[iLC]
                  [iP].mean_gsto_sunlit for _ in range(1)], as_='Gsto_l'),
                # TODO: Multiply by fLAI_sunlit?
                # I(state.canopy_component_population[iLC][iP].fLAI_layer, as_='fLAI'),
                # I([fLAI * state.canopy_layer_component[iL][iLC]
                #   .LAIsunfrac for iL, fLAI in enumerate(state.canopy_component_population[iLC][iP].fLAI_layer)], as_='fLAI'),
                # I([state.canopy_component[iLC].leaf_pop_distribution[iL][iP]
                #    for iL in range(nL)], as_='LAI_per_layer'),
            ],
            additional_inputs=lambda: [
                I(1, as_="nL"),
                I([1], as_='fLAI'),
            ],
            state_outputs=lambda result, iLC=iLC, iP=iP: [
                (result, f'canopy_component_population.{iLC}.{iP}.O3up_sunlit'),
            ],
        )]


def integrate(acc, a, b, ta, tb):
    # return acc + b
    return acc + ((tb - ta) * (a + ((b - a) / 2)))


def calc_fst_leaf_acc_hour_process(iLC: int, iP: int) -> Process:
    """Calculate the accumulated O3up(fst) for each leaf population.

    Note that the use of the leaf_emerg_to_leaf_fst_acc restricts uptake of O3 until
    the LAI is increasing at a rate greater than the knee value. This avoids high
    uptake at the start of the season. It was discussed during DO3SE meetings
    Jan 2024. By default the value will be 0 which disables this feature.

    This is still a bit of a hack and only tested for 2024 China runs.

    """

    def fst_acc(
        O3up, O3up_acc, has_emerged, leaf_emerg_to_leaf_fst_acc, td
    ):
        can_accumulate_fst = has_emerged and (leaf_emerg_to_leaf_fst_acc is None or td > leaf_emerg_to_leaf_fst_acc)
        return 0 if not can_accumulate_fst else  O3up + O3up_acc

    return Process(
        func=fst_acc,
        comment="Accumulate O3up(fst)",
        config_inputs=lambda config: [
            # TODO: Should not always use flag leaf here!
            I(config.Land_Cover.parameters[iLC].phenology.key_lengths_flag_leaf_td.leaf_emerg_to_leaf_fst_acc,
              as_='leaf_emerg_to_leaf_fst_acc'),
        ],
        state_inputs=lambda state, iLC=iLC, iP=iP: [
            I(state.canopy_component[iLC].td_dd_leaf_pops[iP], as_='td'),
            I(state.canopy_component_population[iLC][iP].phenology.phenology_stage >=
              LeafPhenologyStage.GROWING, as_='has_emerged'),
            I(state.canopy_component_population[iLC][iP].O3up, as_='O3up'),
            I(state.canopy_component_population[iLC][iP].O3up_acc_day, as_='O3up_acc'),
        ],
        state_outputs=lambda result, iLC=iLC, iP=iP: [
            (result, f'canopy_component_population.{iLC}.{iP}.O3up_acc_day'),
        ],
    )


def store_fst_acc_prev_day(iLC: int, iP: int) -> Process:
    return Process(
        func=set_value,
        comment="Store the O3up_acc_day from previous day and reset current day to 0.",
        state_inputs=lambda state, iLC=iLC, iP=iP: [
            I(state.canopy_component_population[iLC][iP].O3up_acc_day, as_='O3up'),
        ],
        state_outputs=lambda result, iLC=iLC, iP=iP: [
            (result['O3up'], f'canopy_component_population.{iLC}.{iP}.O3up_acc_day_prev'),
            (0, f'canopy_component_population.{iLC}.{iP}.O3up_acc_day'),
        ],
    )


def calc_fst_leaf_acc_day_process(iLC: int, iP: int) -> Process:
    """Calculate the accumulated O3up(fst) for each leaf population.

    Advanced method uses average daily value.
    """
    return Process(
        # NOTE: Conversion of nmol/s to nmol/hour
        func=lambda O3up, O3up_acc: O3up_acc + O3up * 60 * 60,
        comment="Accumulate O3up(fst)",
        state_inputs=lambda state, iLC=iLC, iP=iP: [
            I(state.canopy_component_population[iLC][iP].O3up_acc_day, as_='O3up'),
            # I(state.canopy_component_population[iLC][iP].O3up_acc_day_prev, as_='O3up_prev'),
            I(state.canopy_component_population[iLC][iP].O3up_acc, as_='O3up_acc'),
        ],
        state_outputs=lambda result, iLC=iLC, iP=iP: [
            (result, f'canopy_component_population.{iLC}.{iP}.O3up_acc'),
        ],
    )


# def calc_fst_leaf_acc_process(iLC: int, iP: int) -> Process:
#     """Calculate the accumulated O3up(fst) for each leaf population.

#     Simple method. Only uses ozone at last hour of day!
#     """
#     return Process(
#         # TODO: Check if we should be using td as Δt
#         # TODO: Should be 0 until leaf pop has emerged.
#         func=lambda O3up, O3up_prev, O3up_acc, td_dd, td_dd_prev: integrate(
#             O3up_acc, O3up_prev, O3up, td_dd_prev, td_dd),
#         comment="Accumulate O3up(fst)",
#         state_inputs=lambda state, iLC=iLC, iP=iP: [
#             I(state.canopy_component_population[iLC][iP].O3up, as_='O3up'),
#             I(state.prev_hour.canopy_component_population[iLC][iP].O3up, as_='O3up_prev'),
#             I(state.canopy_component_population[iLC][iP].O3up_acc, as_='O3up_acc'),
#             I(state.canopy_component[iLC].td_dd, as_='td_dd'),
#             I(state.prev_hour.canopy_component[iLC].td_dd, as_='td_dd_prev'),
#         ],
#         state_outputs=lambda result, iLC=iLC, iP=iP: [
#             (result, f'canopy_component_population.{iLC}.{iP}.O3up_acc'),
#         ],
#     )


def calc_POD_leaf_process(iLC: int, iP: int) -> Process:
    """Calculate the Phytotoxic Ozone Dose."""
    return [
        Process(
            func=O3_helpers.calc_POD,
            comment=f"Calculate the Phytotoxic Ozone Dose iLC:{iLC} iP:{iP}",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].Y, as_='Y', required=True),
            ],
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I(state.canopy_component_population[iLC][iP].O3up, as_='fst', required=True),
                I(state.canopy_component_population[iLC][iP].POD_0,
                  as_='POD_0_prev', required=True),
                I(state.prev_hour.canopy_component_population[iLC][iP].POD_Y,
                  as_='POD_Y_prev', required=True),
                I(state.canopy_component_population[iLC][iP].phenology.phenology_stage >=
                  LeafPhenologyStage.GROWING, as_='has_emerged'),
                I(state.canopy_component_population[iLC][iP].phenology.phenology_stage ==
                  LeafPhenologyStage.GROWING, as_='is_growing'),
            ],
            state_outputs=lambda result, iLC=iLC, iP=iP: [
                (result.POD_0, f'canopy_component_population.{iLC}.{iP}.POD_0'),
                (result.POD_Y, f'canopy_component_population.{iLC}.{iP}.POD_Y'),
            ],
        ),
        Process(
            func=O3_helpers.calc_POD,
            comment=f"Calculate the Phytotoxic Ozone Dose sunlit iLC:{iLC} iP:{iP}",
            config_inputs=lambda config, iLC=iLC: [
                I(config.Land_Cover.parameters[iLC].Y, as_='Y', required=True),
            ],
            state_inputs=lambda state, iLC=iLC, iP=iP: [
                I(state.canopy_component_population[iLC][iP].O3up_sunlit, as_='fst', required=True),
                I(state.canopy_component_population[iLC][iP].POD_0_sunlit,
                  as_='POD_0_prev', required=True),
                I(state.prev_hour.canopy_component_population[iLC][iP].POD_Y_sunlit,
                  as_='POD_Y_prev', required=True),
                I(state.canopy_component_population[iLC][iP].phenology.phenology_stage >=
                  LeafPhenologyStage.GROWING, as_='has_emerged'),
                I(state.canopy_component_population[iLC][iP].phenology.phenology_stage ==
                  LeafPhenologyStage.GROWING, as_='is_growing'),
            ],
            state_outputs=lambda result, iLC=iLC, iP=iP: [
                (result.POD_0, f'canopy_component_population.{iLC}.{iP}.POD_0_sunlit'),
                (result.POD_Y, f'canopy_component_population.{iLC}.{iP}.POD_Y_sunlit'),
            ],
        )]


def calc_OT_leaf_process(iLC: int, iP: int, nL: int, nP: int) -> Process:
    """Calculate the OT accumulation per leaf population."""
    return Process(
        func=O3_helpers.calc_OT_leaf,
        comment="Calculate the OT accumulation",
        gate=iP == nP - 1,  # only calculated for flag leaf
        config_inputs=lambda config, iLC=iLC: [
            I(config.Land_Cover.nL, as_='nL'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(e_state.is_daylight[row_index], as_='is_daylight'),
        ],
        state_inputs=lambda state, iLC=iLC, iP=iP: [
            I([state.canopy_layers[iL].micro_met.micro_O3 for iL in range(nL)], as_='micro_O3'),
            I(state.canopy_layer_component[0][iLC].gsto_params.f_phen, as_='f_phen'),
            I(state.canopy_layer_component[0][iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
            I(state.prev_hour.canopy_component_population[iLC][iP].AOT_0, as_='AOT_0_prev'),
            I(state.prev_hour.canopy_component_population[iLC][iP].AOT_40, as_='AOT_40_prev'),
        ],
        state_outputs=lambda result, iLC=iLC, iP=iP: [
            (result.OT_0, f'canopy_component_population.{iLC}.{iP}.OT_0'),
            (result.OT_40, f'canopy_component_population.{iLC}.{iP}.OT_40'),
            (result.AOT_0, f'canopy_component_population.{iLC}.{iP}.AOT_0'),
            (result.AOT_40, f'canopy_component_population.{iLC}.{iP}.AOT_40'),
        ],
    )


def calc_FO3_eff_process(iL, iLC) -> Process:
    """Calculate FO3_eff Effective ozone dose for effect on photosynthesis."""
    Warning('FO3_eff not implemented')
    return []
    # TODO: Implement FO3_eff calc
    # [Process(
    #     func=O3_helpers.calc_FO3_eff,
    #     comment="Effective ozone dose for effect on photosynthesis",
    #     gate=config.Land_Cover.parameters[iLC].pn_gsto.O3_method == 'martin2000',
    #     config_inputs=lambda config: [
    #         (result.',''),
    #     ],
    #     state_inputs=lambda state: [
    #         (result.',''),
    #     ],
    #     state_outputs=lambda result: [
    #         (result.',''),
    #     ]
    # ) for iL in range(nL) for iLC in range(nLC)],


def check_soil_evaporation_blocked_process() -> Process:
    return Process(
        # TODO: this assumes the first land cover is the only one that matters
        func=SMD_PM_helpers.check_soil_evaporation_blocked,
        comment="Is soil evaporation blocked?",
        # TODO: Split this into function per method
        group="soil-moisture",
        config_inputs=lambda config: [
            I(config.Land_Cover.parameters[0].gsto.f_SW_method, as_='f_SW_method'),
            I(config.Land_Cover.parameters[0].gsto.SWP_max, as_='SWP_max'),
            I(config.Land_Cover.parameters[0].gsto.fmin, as_='fmin'),
            I(config.Land_Cover.parameters[0].gsto.fSWP_exp_a, as_='fSWP_exp_a'),
            I(config.Land_Cover.parameters[0].gsto.fSWP_exp_b, as_='fSWP_exp_b'),
            I(config.soil_moisture.ASW_FC, as_='ASW_FC'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.SMD.SWP, as_='SWP'),
            I(state.canopy.SMD.ASW, as_='ASW'),
        ],
        state_outputs=lambda result: [
            (result, 'canopy.Es_blocked'),
        ],
    )


def penman_monteith_hourly_process(gsto_method) -> Process:
    return [
        Process(
            # NOTE: Not valid for multilayer or pn_gsto method
            func=SMD_PM_helpers.calc_PET_gsto,
            comment="P-M - Calculate PEt RSTO",
            group="soil-moisture",
            gate=gsto_method == "multiplicative",
            config_inputs=lambda config: [
                I(config.Land_Cover.parameters[0].multip_gsto.gmax, as_='gmax'),
            ],
            state_inputs=lambda state: [
                I(state.canopy_layer_component[0][0].gsto_params.f_phen, as_='f_phen'),
                I(state.canopy_layer_component[0][0].gsto_params.f_light, as_='f_light'),
                I(state.canopy_layer_component[0][0].gsto_params.f_temp, as_='f_temp'),
                I(state.canopy_layer_component[0][0].gsto_params.f_VPD, as_='f_VPD'),
                I(state.canopy_layer_component[0][0].LAI, as_='LAI'),
            ],
            state_outputs=lambda result: [
                (result, 'canopy.PM.PEt_rsto'),
            ]
        ),
        Process(
            # NOTE: Not valid for multilayer or pn_gsto method
            func=SMD_PM_helpers.calc_PET_gsto_An,
            comment="P-M - Calculate PEt RSTO",
            group="soil-moisture",
            gate=gsto_method == "photosynthesis",
            state_inputs=lambda state: [
                # NOTE: Assumes using first component
                I(state.canopy_component[0].canopy_gsto, as_='bulk_gsto'),
                # NOTE: Assumes using top layer
                I(state.canopy_layer_component[0][0].gsto_params.f_SW, as_='f_SW'),
            ],
            state_outputs=lambda result: [
                (result, 'canopy.PM.PEt_rsto'),
            ]
        ),
        Process(
            func=SMD_PM_helpers.penman_monteith_hourly,
            comment="P-M - Hourly Penman-Monteith calculations for evaporation and transpiration",
            group="soil-moisture",
            external_state_inputs=lambda e_state, row_index: [
                I(lget(e_state.Rn, row_index), as_='Rn_MJ'),
                I(lget(e_state.P, row_index), as_='P_kPa'),
                I(e_state.Ts_C[row_index], as_='Ts_C'),
                I(lget(e_state.esat, row_index), as_='esat_kPa'),
                I(lget(e_state.eact, row_index), as_='eact_kPa'),
                I(lget(e_state.VPD, row_index), as_='VPD_kPa'),
            ],
            state_inputs=lambda state: [
                I(state.canopy.rmodel_H2O.nL, as_='rm_h2o_nL'),
                I(state.canopy.rmodel_H2O.Rb, as_='rm_h2o_Rb'),
                # TODO: Check multi layer model
                I(state.canopy.rmodel_H2O.Rinc[0], as_='rm_h2o_Rinc_l0'),
                I(state.canopy.rmodel_H2O.Rsto[0], as_='rm_h2o_Rsto_l0'),
                I(state.canopy.rmodel_H2O.Rgs, as_='rm_h2o_Rgs'),
                I(state.canopy.LAI_total, as_='LAI'),
                I(state.canopy.Es_blocked, as_='Es_blocked'),
                I(state.canopy.PM.PEt_rsto, as_='rm_pet_Rsto_l0'),
                I(state.canopy.PM.Ei_acc, as_='pm_state_Ei_acc'),
                I(state.canopy.PM.PEt_acc, as_='pm_state_PEt_acc'),
                I(state.canopy.PM.Et_acc, as_='pm_state_Et_acc'),
                I(state.canopy.PM.Es_acc, as_='pm_state_Es_acc'),
                I(state.canopy.PM.Eat_acc, as_='pm_state_Eat_acc'),
            ],
            state_outputs=lambda result: [
                (result.Ei_acc, 'canopy.PM.Ei_acc'),
                (result.Et_acc, 'canopy.PM.Et_acc'),
                (result.PEt_acc, 'canopy.PM.PEt_acc'),
                (result.Es_acc, 'canopy.PM.Es_acc'),
                (result.Eat_acc, 'canopy.PM.Eat_acc'),
                (result.Ei_hr, 'canopy.PM.Ei_hr'),
                (result.Et_hr, 'canopy.PM.Et_hr'),
                (result.PEt_hr, 'canopy.PM.PEt_hr'),
                (result.Es_hr, 'canopy.PM.Es_hr'),
                (result.Eat_hr, 'canopy.PM.Eat_hr'),
            ],
        )]


def penman_monteith_daily_process():
    """Calculate the daily sn diff and accumulate run off and percolated."""
    return Process(
        func=SMD_PM_helpers.penman_monteith_daily,
        comment="P-M - Daily soil water content update from accumulated \
Penman-Monteith values",
        group="soil-moisture",
        config_inputs=lambda config: [
            I(config.soil_moisture.root, as_='root_depth'),
            I(config.soil_moisture.run_off_fraction, as_='run_off_fraction'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.LAI_total, as_='LAI'),
            I(state.canopy.SMD.ASW, as_='ASW'),
            I(state.canopy.SMD.SMD, as_='SMD'),
            I(state.canopy.PM.precip_acc_prev_day, as_='pm_state_precip_acc'),
            I(state.canopy.PM.run_off_acc, as_='pm_state_run_off_acc'),
            I(state.canopy.PM.Ei_acc, as_='pm_state_Ei_acc'),
            I(state.canopy.PM.Eat_acc, as_='pm_state_Eat_acc'),
            I(state.canopy.PM.percolated_acc, as_='pm_state_percolated_acc'),
        ],
        state_outputs=lambda result: [
            (result.rain_input, 'canopy.PM.rain_input'),
            (result.run_off, 'canopy.PM.run_off'),
            (result.run_off_acc, 'canopy.PM.run_off_acc'),
            (result.effective_irrig, 'canopy.PM.effective_irrig'),
            (result.intercepted_evaporated, 'canopy.PM.intercepted_evaporated'),
            (result.evapotranspiration, 'canopy.PM.evapotranspiration'),
            (result.Sn_diff, 'canopy.PM.Sn_diff'),
            (result.percolated, 'canopy.PM.percolated'),
            (result.percolated_acc, 'canopy.PM.percolated_acc'),
        ],
    )


def calc_h2o_r_model_process() -> List[Process]:
    return [
        Process(
            func=SMD_PM_helpers.multi_layer_r_model_to_single_H20,
            comment="Adapt multi-layer O3 resistance model to single-layer H2O resistance",
            group="soil-moisture",
            state_inputs=lambda state: [
                I(state.canopy.rmodel_O3, as_='rmodel'),
                I(state.met.ustar, as_='ustar'),
                I(state.canopy.LAI_total, as_='total_LAI'),
                I(state.canopy.SAI_total, as_='total_SAI'),
            ],
            state_outputs=lambda result: [(result, 'canopy.rmodel_H2O')],
        ),
        Process(
            func=resistance_helpers.calc_Rtotal,
            comment="Calculate the total resistance in H2O resistance model",
            group="soil-moisture",
            config_inputs=lambda config: [I(1, as_='nL')],
            # config_inputs=lambda config: [I(config.Land_Cover.nL, as_='nL')],
            state_inputs=lambda state: [
                I(state.canopy.rmodel_H2O.Rsur, as_='Rsur'),
                I(state.canopy.rmodel_H2O.Rinc, as_='Rinc'),
                I(state.canopy.rmodel_H2O.Rgs, as_='Rgs'),
            ],
            state_outputs=lambda result: [(result, 'canopy.rmodel_H2O.Rtotal')]
        ),
    ]


def soil_moisture_PM_process():
    return Process(
        func=SMD_PM_helpers.PM_soil_moisture_calc,
        comment="P-M - soil moisture calc",
        group="soil-moisture",
        config_inputs=lambda config: [
            I(config.soil_moisture.soil_config, as_='soil_config'),
            I(config.soil_moisture.PWP, as_='PWP'),
            I(config.soil_moisture.root, as_='root_depth'),
        ],
        state_inputs=lambda state: [
            I(state.canopy.SMD.Sn, as_='Sn_in'),
            I(state.canopy.PM.Sn_diff, as_='Sn_diff'),
        ],
        state_outputs=lambda result: [
            # ASW Output here is incorrect
            (result.Sn, 'canopy.SMD.Sn'),
            (result.SWP, 'canopy.SMD.SWP'),
            (result.ASW, 'canopy.SMD.ASW'),
            (result.SMD, 'canopy.SMD.SMD'),
        ]
    )


def penman_monteith_reset_process():
    return Process(
        func=SMD_PM_helpers.penman_monteith_reset,
        comment="reset penman monteith daily accumulators to 0",
        group="soil-moisture",
        state_outputs=lambda result: [
            (result.Ei_acc, 'canopy.PM.Ei_acc'),
            (result.Et_acc, 'canopy.PM.Et_acc'),
            (result.Es_acc, 'canopy.PM.Es_acc'),
            (result.Eat_acc, 'canopy.PM.Eat_acc'),
        ],
    )


def calc_daily_carbon_allocation_process(iLC) -> Process:
    return Process(
        func=carbon_calcs.daily_carbon_allocation,
        comment="Calculate the daily carbon allocation",
        config_inputs=lambda config: [
            I(config.carbon_allocation.a_root, as_='a_root'),
            I(config.carbon_allocation.a_stem, as_='a_stem'),
            I(config.carbon_allocation.a_leaf, as_='a_leaf'),
            I(config.carbon_allocation.b_root, as_='b_root'),
            I(config.carbon_allocation.b_stem, as_='b_stem'),
            I(config.carbon_allocation.b_leaf, as_='b_leaf'),
            I(config.carbon_allocation.theta, as_='theta'),
            I(config.carbon_allocation.f_green_leaf, as_='f_green_leaf'),
            I(config.carbon_allocation.f_brown_leaf, as_='f_brown_leaf'),
        ],
        state_inputs=lambda state, iLC=iLC: [
            I(state.canopy_component[iLC].npp_acc, as_='net_prod_acc'),
            I(state.canopy_component[iLC].dvi, as_='DVI'),
            I(
                state.canopy_component_population[iLC][-1].phenology.phenology_stage >= LeafPhenologyStage.SENESCENCE,
                as_='plant_is_senescing',
            ),
            I(state.canopy_component[iLC].c_root, as_='c_root'),
            I(state.canopy_component[iLC].c_stem, as_='c_stem'),
            I(state.canopy_component[iLC].c_leaf, as_='c_leaf'),
            I(state.canopy_component[iLC].c_harv, as_='c_harv'),
            I(state.canopy_component[iLC].c_resv, as_='c_resv'),
            I(state.canopy_component[iLC].c_lbrn, as_='c_lbrn'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result.c_root, f'canopy_component.{iLC}.c_root'),
            (result.c_stem, f'canopy_component.{iLC}.c_stem'),
            (result.c_leaf, f'canopy_component.{iLC}.c_leaf'),
            (result.c_lbrn, f'canopy_component.{iLC}.c_lbrn'),
            (result.c_harv, f'canopy_component.{iLC}.c_harv'),
            (result.c_resv, f'canopy_component.{iLC}.c_resv'),
            (result.p_root, f'canopy_component.{iLC}.p_root'),
            (result.p_stem, f'canopy_component.{iLC}.p_stem'),
            (result.p_leaf, f'canopy_component.{iLC}.p_leaf'),
            (result.p_harv, f'canopy_component.{iLC}.p_harv'),
        ],
    )


def calc_yield_process(iLC: int) -> Process:
    return Process(
        func=carbon_calcs.calculate_yield,
        comment="Calculate the daily yield",
        config_inputs=lambda config: [
            I(config.carbon_allocation.dry_to_wet_biomass, as_='dry_to_wet_biomass'),
            I(config.carbon_allocation.grain_to_ear, as_='grain_to_ear'),
        ],
        state_inputs=lambda state, iLC=iLC: [
            I(state.canopy_component[iLC].c_root, as_='c_root'),
            I(state.canopy_component[iLC].c_stem, as_='c_stem'),
            I(state.canopy_component[iLC].c_leaf, as_='c_leaf'),
            I(state.canopy_component[iLC].c_harv, as_='c_harv'),
            I(state.canopy_component[iLC].c_resv, as_='c_resv'),
            I(state.canopy_component[iLC].c_lbrn, as_='c_lbrn'),
        ],
        state_outputs=lambda result, iLC=iLC: [
            (result.stem_dm, f'canopy_component.{iLC}.stem_dm'),
            (result.leaf_dm, f'canopy_component.{iLC}.leaf_dm'),
            (result.lbrn_dm, f'canopy_component.{iLC}.lbrn_dm'),
            (result.total_leaf_dm, f'canopy_component.{iLC}.total_leaf_dm'),
            (result.straw_dm, f'canopy_component.{iLC}.straw_dm'),
            (result.ear_dm, f'canopy_component.{iLC}.ear_dm'),
            (result.aboveground_dm, f'canopy_component.{iLC}.aboveground_dm'),
            (result.belowground_dm, f'canopy_component.{iLC}.belowground_dm'),
            (result.grain_dm, f'canopy_component.{iLC}.grain_dm'),
            (result.harvest_index, f'canopy_component.{iLC}.harvest_index'),
            (result.yield_ha, f'canopy_component.{iLC}.yield_ha'),
        ],
    )


def calculate_LAI_from_DVI_and_carbon_process(iLC) -> Process:
    return [
        Process(
            func=carbon_calcs.calc_LAI_from_DVI_and_carbon,
            comment="Calculate the plant LAI from accumulated carbon",
            config_inputs=lambda config: [
                I(config.carbon_allocation.gamma, as_='gamma'),
                I(config.carbon_allocation.delta, as_='delta'),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].dvi, as_='DVI'),
                I(state.canopy_component[iLC].c_leaf, as_='c_leaf'),
                I(state.canopy_component[iLC].total_emerged_leaf_populations,
                  as_='emerged_leaf_count'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.LAI'),
            ],
        ),
        Process(
            func=carbon_calcs.calc_LAI_from_DVI_and_carbon,
            comment="Calculate the plant brownleaf LAI from accumulated carbon",
            config_inputs=lambda config: [
                I(config.carbon_allocation.gamma, as_='gamma'),
                I(config.carbon_allocation.delta, as_='delta'),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].dvi, as_='DVI'),
                I(state.canopy_component[iLC].c_lbrn, as_='c_leaf'),
                I(state.canopy_component[iLC].total_emerged_leaf_populations,
                  as_='emerged_leaf_count'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.LAI_brown'),
            ],
        ),
    ]


def calc_canopy_displacement_height_process() -> Process:
    return Process(
        func=resistance_helpers.calc_displacement_and_roughness_parameters,
        comment="calculate measured wind canopy displacement parameters",
        config_inputs=lambda config: [
            I(config.Land_Cover.land_cover_type == LandCoverType.FOREST, as_="is_forest"),
        ],
        state_inputs=lambda state: [
            I(state.canopy.canopy_height, as_="h"),
        ],
        state_outputs=lambda result: [
            (result[0], 'canopy.d'),
            (result[1], 'canopy.z0'),
        ],
    )


def get_plant_height_from_carbon_process(iLC) -> Process:
    return [
        Process(
            func=carbon_calcs.get_plant_height_from_carbon,
            comment="Calculate the plant height from accumulated carbon",
            config_inputs=lambda config: [
                I(config.carbon_allocation.k, as_='k'),
            ],
            state_inputs=lambda state, iLC=iLC: [
                I(state.canopy_component[iLC].c_stem, as_='c_stem'),
            ],
            state_outputs=lambda result, iLC=iLC: [
                (result, f'canopy_component.{iLC}.plant_height'),
            ]
        ),
    ]


def reset_carbon_accumulators_process(nLC):
    return Process(
        func=set_value,
        comment="reset daily carbon accumulators",
        state_outputs=lambda result: [
            (0, f'canopy_component.{iLC}.npp_acc') for iLC in range(nLC)
        ],
    )


# TODO: Get from external data?
def set_hour(hr: int) -> Process:
    """Set the state hour."""
    return Process(
        func=set_value,
        comment="Set Hour",
        additional_inputs=lambda: [I(hr, as_='hr')],
        state_outputs=lambda result: [(result['hr'], 'temporal.hr')]
    )


# def set_row_index() -> Process:
#     """Set the current data row index."""
#     return Process(
#         func=get_row_index,
#         comment="Set data row index",
#         state_inputs=lambda state: [
#             I(state.temporal.dd, as_='day'),
#             I(state.temporal.hr, as_='hour'),
#         ],
#         state_outputs=lambda result: [
#             (result, 'temporal.row_index')
#         ]
#     )

def sync_row_index() -> Process:
    """Sync the internal state row index with the external state row index."""
    return Process(
        func=set_value,
        comment="Set data row index",
        external_state_inputs=lambda e_state, row_index: [
            I(row_index, as_="row_index"),
        ],
        state_outputs=lambda result: [
            (result['row_index'], 'temporal.row_index')
        ]
    )


# def store_prev_state() -> List[Process]:
#     """Store a copy of the previous State."""
#     return [
#         Process(
#             func=lambda: Model_State_Shape(),
#             comment="Reset previous state",
#             state_outputs=lambda result: [
#                 (result, 'prev_hour'),
#             ]
#         ),
#         Process(
#             func=set_value,
#             comment="Copy State",
#             state_inputs=lambda state: [
#                 I(state.temporal, as_='temporal'),
#                 I(state.external_met, as_='external_met'),
#                 I(state.met, as_='met'),
#                 I(state.canopy, as_='canopy'),
#                 I(state.canopy_layers, as_='canopy_layers'),
#                 I(state.canopy_component, as_='canopy_component'),
#                 I(state.canopy_layer_component, as_='canopy_layer_component'),
#                 I(state.canopy_component_population, as_='canopy_component_population'),
#             ],
#             # Note we need to deep copy here to ensure we don't modify prev state
#             # when updating current state
#             state_outputs=lambda result: [
#                 (deepcopy(result['temporal']), 'prev_hour.temporal'),
#                 (deepcopy(result['external_met']), 'prev_hour.external_met'),
#                 (deepcopy(result['met']), 'prev_hour.met'),
#                 (deepcopy(result['canopy']), 'prev_hour.canopy'),
#                 (deepcopy(result['canopy_layers']), 'prev_hour.canopy_layers'),
#                 (deepcopy(result['canopy_component']), 'prev_hour.canopy_component'),
#                 (deepcopy(result['canopy_layer_component']), 'prev_hour.canopy_layer_component'),
#                 (deepcopy(result['canopy_component_population']),
#                  'prev_hour.canopy_component_population'),
#             ]
#         ),
#     ]
def store_prev_state() -> List[Process]:
    """Store a copy of the previous State."""
    return Process(
        func=set_value,
        comment="Copy State",
        state_inputs=lambda state: [
            I(state, as_='state'),
        ],
        # Note we need to deep copy here to ensure we don't modify prev state
        # when updating current state
        state_outputs=lambda result: [
            (
                Model_State_Shape(
                    temporal=deepcopy(result['state'].temporal),
                    external_met=deepcopy(result['state'].external_met),
                    met=deepcopy(result['state'].met),
                    canopy=deepcopy(result['state'].canopy),
                    canopy_layers=deepcopy(result['state'].canopy_layers),
                    canopy_component=deepcopy(result['state'].canopy_component),
                    canopy_layer_component=deepcopy(result['state'].canopy_layer_component),
                    canopy_component_population=deepcopy(
                        result['state'].canopy_component_population),
                ), 'prev_hour'
            ),
        ],
    )


def set_day_func(dd: int, dd_offset: int) -> int:
    return int(dd) + int(dd_offset)


def set_day_offset_func(current_dd: int, current_dd_offset: int, multi_season_day: int) -> int:
    if current_dd > multi_season_day:
        return current_dd_offset - 365
    else:
        return current_dd_offset


def set_day_offset(multi_season: bool):
    """Set the offset day if using multi season."""
    return [
        Process(
            func=set_day_offset_func,
            gate=multi_season,
            comment="set day offset",
            config_inputs=lambda config: [
                I(config.Location.multi_season_day, as_='multi_season_day'),
            ],
            state_inputs=lambda state: [
                I(state.temporal.dd, as_='current_dd'),
                I(state.temporal.dd_offset, as_='current_dd_offset'),
            ],
            state_outputs=lambda result: [
                (result, 'temporal.dd_offset'),
            ]
        ),
        # TODO: Setup if not using multi season
        # Process(
        #     func=set_day_func,
        #     comment="set day",
        #     gate=not multi_season,
        #     config_inputs=lambda config: [
        #         I(config.Met.dd_offset, as_='dd_offset'),
        #     ],
        #     additional_inputs=lambda: [I(dd, as_='dd')],
        #     state_outputs=lambda result: [(result, 'temporal.dd')]
        # ),
    ]


def reset_model():
    # TODO: We should reset the model at the end of the season
    return []


def set_day(dd: int = None) -> Process:
    """Set the state day.

    If dd is None we use dd from external state.
    """

    return Process(
        func=set_day_func,
        comment="set day",
        state_inputs=lambda state: [
            I(state.temporal.dd_offset, as_='dd_offset'),
        ],
        additional_inputs=lambda: [I(dd, as_='dd')],
        state_outputs=lambda result: [(result, 'temporal.dd')]
    ) if dd is not None else Process(
        func=set_day_func,
        comment="set day",
        state_inputs=lambda state: [
            I(state.temporal.dd_offset, as_='dd_offset'),
        ],
        external_state_inputs=lambda e_state, row_index: [
            I(lget(e_state.dd, row_index, f"e_state.dd is None at row_index: {row_index}!"), as_='dd')],
        state_outputs=lambda result: [
            (result, 'temporal.dd'),
        ]
    )


def save_state(run_dir: Path, dump_state_n: int) -> Process:
    if dump_state_n and not run_dir:
        raise ValueError("Must supply run_dir when using dump_state_n")
    os.makedirs(f"{run_dir}/running_state", exist_ok=True)

    def _save_process(state: Model_State_Shape, run_dir: Path, n: int):
        row_index = state.temporal.row_index
        if row_index >= int(n) and int(row_index) % int(n) == 0:
            out_path = f"{run_dir}/running_state/dump_state_{row_index}.json"
            print("Saving dumped state to", out_path)
            dump_state_to_file(state, out_path)

    return Process(
        func=_save_process,
        comment="Save out state",
        additional_inputs=lambda: [
            I(run_dir, as_="run_dir"),
            I(dump_state_n, as_="n"),
        ],
        state_inputs=lambda state: [
            I(state, as_="state"),
        ],
    )


def daily_start_processes(config: Config_Shape, run_dir: Path, dump_state_n: int) -> List[Process]:
    """Get processes to be ran at start of day."""
    nL = config.Land_Cover.nL
    nLC = config.Land_Cover.nLC
    nP = config.Land_Cover.nP
    thermal_time_method = config.Met.thermal_time_method
    multi_season = config.Location.multi_season
    time_type = config.Land_Cover.phenology_options.time_type
    height_method = config.Land_Cover.height_method
    LAI_method = config.Land_Cover.LAI_method
    LAI_distribution_method = config.Land_Cover.LAI_distribution_method
    SAI_method = config.Land_Cover.SAI_method
    primary_SAI_method = config.Land_Cover.parameters[0].phenology.SAI_method
    phenology_method = config.Land_Cover.phenology_options.phenology_method
    f_phen_method = [
        config.Land_Cover.parameters[iLC].phenology.f_phen_method for iLC in range(nLC)]
    leaf_f_phen_method = [
        config.Land_Cover.parameters[iLC].phenology.leaf_f_phen_method for iLC in range(nLC)]
    # f_SW_method = [config.Land_Cover.parameters[iLC].gsto.f_SW_method for iLC in range(nLC)]
    soil_moisture_source = config.soil_moisture.source
    use_carbon_allocation = config.carbon_allocation.use_carbon_allocation
    dvi_method = config.Land_Cover.dvi_method
    plant_emerge_method = config.Land_Cover.phenology_options.plant_emerge_method
    flag_leaf_emerge_method = config.Land_Cover.phenology_options.flag_leaf_emerge_method
    use_vernalisation = config.Land_Cover.phenology_options.use_vernalisation
    use_photoperiod_factor = config.Land_Cover.phenology_options.use_photoperiod_factor
    sparse_data = config.Met.sparse_data
    have_precip_data = config.Met.inputs.precip_method != InputMethod.SKIP

    return [
        tag_process("===== Start of Daily Processes ====="),
        set_day_offset(multi_season),
        reset_model() if multi_season else [],
        set_day(),
        save_state(run_dir, dump_state_n) if dump_state_n is not None else [],

        # MET
        tag_process("== Met Processes =="),
        # perNlc(nLC, calc_effective_temperature_process),
        [calc_effective_temperature_process(iLC) for iLC in range(nLC)] if not sparse_data else [],
        [calculate_daily_thermal_time_process(
            thermal_time_method, iLC) for iLC in range(nLC)] \
        if thermal_time_method != ThermalTimeMethods.EXTERNAL and not sparse_data else [],
        calc_photoperiod_process(),
        [calc_photoperiod_factor_process(iLC)
         for iLC in range(nLC)]if use_photoperiod_factor else [],

        # PHENOLOGY
        tag_process("== Phenology Processes =="),
        # f_phen
        [[
            f_phen_method_process(iL, iLC, f_phen_method[iLC]),
            leaf_f_phen_process(iL, iLC, leaf_f_phen_method[iLC]),
        ] for iL in range(nL) for iLC in range(nLC)] if phenology_method != PhenologyMethods.DISABLED else [],


        [
            [
                calc_DVI_process(iLC, dvi_method),
                calc_phyllochron_process(iLC),
                calculate_relative_photoperiod_process(iLC),
                calc_if_plant_is_sown_process(iLC),
                calc_if_plant_has_emerged_process(iLC, plant_emerge_method),
                get_vernalisation_factor_process(iLC) if use_vernalisation else [],
                get_td_dd_process(iLC) if time_type == TimeTypes.THERMAL_TIME else [],
                calc_if_flag_leaf_has_emerged_process(iLC, flag_leaf_emerge_method),
                calc_emergence_rate_process(nP, iLC),
                calc_emerged_leaf_count_process(nP, iLC, time_type),
                calc_td_dd_per_leaf_pop_process(iLC) if time_type == TimeTypes.THERMAL_TIME else [],
                get_phenology_stage_process_td(
                    nP, iLC) if time_type == TimeTypes.THERMAL_TIME else [],
                get_phenology_stage_process(
                    nP, iLC, leaf_f_phen_method[0]) if time_type == TimeTypes.JULIAN_DAY else [],
                get_growing_populations_process(nP, iLC),
            ]
            for iLC in range(nLC)
        ] if phenology_method != PhenologyMethods.DISABLED else [],

        # Carbon Allocation
        tag_process("== Carbon Allocation Processes =="),
        [
            [calc_daily_carbon_allocation_process(iLC) for iLC in range(nLC)],
            [calc_yield_process(iLC) for iLC in range(nLC)],
            [calculate_LAI_from_DVI_and_carbon_process(iLC) for iLC in range(
                nLC)] if LAI_method == LAIMethods.CARBON else [],
            [get_plant_height_from_carbon_process(iLC) for iLC in range(nLC)],
            reset_carbon_accumulators_process(nLC),
        ] if use_carbon_allocation else [],


        # CANOPY STRUCTURE
        tag_process("== Canopy Structure Processes =="),
        calc_canopy_height_process(nL, nLC, height_method),
        calc_canopy_LAI_processes(nLC, LAI_method),
        distribute_lai_per_layer_processes(nL, nLC, LAI_distribution_method, LAI_method),
        calc_canopy_SAI(nL, nLC, SAI_method, primary_SAI_method),
        distribute_lai_to_leaf_populations_processes(nL, nP, nLC, LAI_method),
        calc_canopy_displacement_height_process(),

        # Soil Moisture
        tag_process("== Soil Moisture Processes =="),
        store_accumulate_precipitation_process() if have_precip_data else [],

        [
            penman_monteith_daily_process(),
            soil_moisture_PM_process(),
            penman_monteith_reset_process(),
        ] if soil_moisture_source == SoilMoistureSource.P_M else [],

        # fst
        tag_process("== Fst Processes =="),
        [calc_fst_leaf_acc_day_process(iLC, iP) for iLC in range(nLC) for iP in range(nP)],
        [store_fst_acc_prev_day(iLC, iP) for iLC in range(nLC) for iP in range(nP)],
        tag_process("===== End of Daily Processes ====="),
    ]


def log_external_state_values(
    external_state_inputs: Callable[[Model_State_Shape], List[I]]
) -> Process:
    return Process(
        func=merge_logs,
        ptype=ProcessType.LOG,
        comment="Logging External Values",
        external_state_inputs=lambda e_state, row_index: [
            *external_state_inputs(e_state, row_index)
        ],
        format_output=True,
    )


# def dump_state_process():
#     # === DEBUG MODE ===
#     log_row_ids = [
#         "71_1",
#         "71_2",
#         "71_12",
#         "71_13",
#         "100_1",
#         "100_2",
#         "110_12",
#         "110_13",
#     ]
#     return Process(
#         func=lambda state, target_path, dd, hr: dump_state_to_file(
#             state, target_path) if f"{dd}_{hr}" in log_row_ids else None,
#         state_inputs=lambda state: [
#             I(state, 'state'),
#             I(state.temporal.dd, 'dd'),
#             I(state.temporal.hr, 'hr'),
#             I(f'tmp/{state.temporal.dd}_{state.temporal.hr}.json', 'target_path'),
#         ],
#     )


def log_processes(nL: int, nLC: int, nP: int, fields: List[str], log_multilayer: bool = False):
    """Log state to output data.

    Where a single layer is used should output top layer
    Where a single population is used should output flag leaf

    """
    flag_index = nP - 1
    top_layer_index = nL - 1
    component_index = 0
    def __SPACER__(name): return I(f'|{name}|', as_='|')
    # TODO: Can we optimize the below logging?
    return [
        # TODO: This should come from config or output_values_map
        log_external_state_values(lambda e_state, row_index: list(filter(lambda f: f, [
            __SPACER__("External Values >"),
            I(row_index, as_='row_index') if 'row_index' in fields else None,
            I(lget(e_state.dd, row_index), as_='dd_e') if 'dd_e' in fields else None,
            I(lget(e_state.hr, row_index), as_='hr') if 'hr' in fields else None,
            I(lget(e_state.Ts_C, row_index), as_='ts_c') if 'ts_c' in fields else None,
            I(lget(e_state.RH, row_index), as_='rh') if 'rh' in fields else None,
            I(lget(e_state.PAR, row_index), as_='par') if 'par' in fields else None,
            I(lget(e_state.P, row_index), as_='p') if 'p' in fields else None,
            I(lget(e_state.u, row_index), as_='uh_zr') if 'uh_zr' in fields else None,
            I(lget(e_state.precip, row_index), as_='precip') if 'precip' in fields else None,
            I(lget(e_state.O3, row_index), as_='o3_ppb_zr') if 'o3_ppb_zr' in fields else None,
            I(lget(e_state.RH, row_index), as_='rh') if 'rh' in fields else None,
            I(lget(e_state.CO2, row_index), as_='co2') if 'co2' in fields else None,
            I(lget(e_state.Hd, row_index), as_='hd') if 'hd' in fields else None,
            I(lget(e_state.R, row_index), as_='r') if 'r' in fields else None,
            I(lget(e_state.Rn, row_index), as_='rn') if 'rn' in fields else None,
            I(lget(e_state.VPD, row_index), as_='vpd') if 'vpd' in fields else None,
            I(lget(e_state.SWC, row_index), as_='swc') if 'swc' in fields else None,
        ]))),
        log_values(lambda state: list(filter(lambda f: f, [
            # Values here should match values in Output_Shape
            __SPACER__("Canopy Level >"),
            I(state.temporal.dd, as_='dd') if 'dd' in fields else None,
            I(state.temporal.dd_offset, as_='dd_offset') if 'dd_offset' in fields else None,
            I(state.canopy_component[component_index].td, as_='td') if 'td' in fields else None,
            I(state.canopy.LAI_total, as_='canopy_lai') if 'canopy_lai' in fields else None,
            I(state.canopy.LAI_brown_total, as_='canopy_lai_brown') if 'canopy_lai_brown' in fields else None,
            I(state.canopy.SAI_total, as_='canopy_sai') if 'canopy_sai' in fields else None,
            I(state.canopy.Vd, as_='canopy_vd') if 'canopy_vd' in fields else None,
            I(state.canopy_component[component_index].dvi, as_='dvi') if 'dvi' in fields else None,
            I(state.canopy_component[component_index].td_dd,
              as_='td_dd') if 'td_dd' in fields else None,
            I(state.canopy_component[component_index].t_eff,
              as_='t_eff') if 't_eff' in fields else None,
            I(state.canopy_component[component_index].td_v,
              as_='td_v') if 'td_v' in fields else None,
            I(state.canopy_component[component_index].phenology.phenology_stage,
              as_='phenology_stage') if 'phenology_stage' in fields else None,
            I(state.canopy_component[0].growing_populations[flag_index],
              as_=f"growing_populations_flag") if 'growing_populations' in fields else [],
            I(state.canopy_component[component_index].phenology.Vf,
              as_='Vf') if 'Vf' in fields else None,
            I(state.canopy_component[component_index].phenology.V_acc,
              as_='V_acc') if 'V_acc' in fields else None,
            I(state.canopy_component[component_index].phenology.V_tot,
              as_='V_tot') if 'V_tot' in fields else None,
            I(state.canopy_component[component_index].phenology.V_pos,
                as_='V_pos') if 'V_pos' in fields else None,
            I(state.canopy_component[component_index].phenology.V_neg,
                as_='V_neg') if 'V_neg' in fields else None,
            I(state.external_met.photoperiod, as_='photoperiod') if 'photoperiod' in fields else None,
            I(state.canopy_component[component_index].photoperiod_factor,
              as_='photoperiod_factor') if 'photoperiod_factor' in fields else None,
            I(state.met.O3_i, as_='o3_ppb_i') if 'o3_ppb_i' in fields else None,
            I(state.met.u_i, as_='u_i') if 'u_i' in fields else None,
            I(state.met.ustar, as_='ustar') if 'ustar' in fields else None,
            I(state.met.ustar_ref, as_='ustar_ref') if 'ustar_ref' in fields else None,
            I(state.canopy.canopy_height, as_='canopy_height') if 'canopy_height' in fields else None,
            I(state.canopy_component[component_index].total_emerged_leaf_populations,
              as_="total_emerged_leaves") if 'total_emerged_leaves' in fields else None,
            __SPACER__("Top Layer >"),
            I(state.canopy_layers[top_layer_index].micro_met.micro_u,
              as_='micro_u') if 'micro_u' in fields else None,
            I(state.canopy_layers[top_layer_index].micro_met.micro_O3,
              as_='micro_O3') if 'micro_O3' in fields else None,
            I(state.canopy.canopy_top_o3,
              as_='canopy_top_o3') if 'canopy_top_o3' in fields else None,
            I(state.canopy_layers[top_layer_index].micro_met.PARsun,
              as_='PARsun') if 'PARsun' in fields else None,
            I(state.canopy_layers[top_layer_index].micro_met.PARshade,
              as_='PARshade') if 'PARshade' in fields else None,
            I(state.canopy_layers[top_layer_index].micro_met.micro_O3,
              as_='o3_ppb') if 'o3_ppb' in fields else None,
            I(state.canopy_layers[top_layer_index].micro_met.micro_O3_nmol_m3,
              as_='o3_nmol_m3') if 'o3_nmol_m3' in fields else None,
            I(state.canopy_layers[top_layer_index].layer_height,
              as_='layer_height') if 'layer_height' in fields else None,
            I(state.canopy_layer_component[top_layer_index][component_index].mean_gsto,
              as_='gsto') if 'gsto' in fields else None,
            I(state.canopy_layer_component[top_layer_index][component_index].bulk_gsto,
              as_='gsto_bulk') if 'gsto_bulk' in fields else None,
            __SPACER__("Flag leaf >"),
            I(state.canopy_component_population[component_index]
              [flag_index].LAI, as_='LAI_flag') if 'LAI_flag' in fields else None,
            I(state.canopy_component_population[0][flag_index].fLAI_layer[top_layer_index],
                as_=f"fLAI_iL_{top_layer_index}_flag") if 'fLAI' in fields else [],
            I(state.canopy_component_population[component_index]
              [flag_index].A_n, as_='A_n') if 'A_n' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].A_n_sunlit, as_='A_n_sunlit') if 'A_n_sunlit' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].A_c, as_='A_c') if 'A_c' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].A_j, as_='A_j') if 'A_j' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].A_p, as_='A_p') if 'A_p' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].A_n_limit_factor, as_='A_n_limit_factor') if 'A_n_limit_factor' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].O3up, as_='fst') if 'fst' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].O3up_sunlit, as_='fst_sun') if 'fst_sun' in fields else None,
            I(sum([state.canopy_component_population[component_index]
                   [iP].O3up for iP in range(nP)]), as_='fst_canopy') if 'fst_canopy' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].O3up_acc, as_='fst_acc') if 'fst_acc' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].O3up_acc_day, as_='fst_acc_day') if 'fst_acc_day' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].O3up_acc_day_prev, as_='fst_acc_day_prev') if 'fst_acc_day_prev' in fields else None,

            I(state.canopy_component_population[component_index]
              [flag_index].leaf_rmodel_O3.Rsto[top_layer_index], as_='rsto_l') if 'rsto_l' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].leaf_rmodel_O3.Rb[top_layer_index], as_='rb_l') if 'rb_l' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].mean_gsto_per_layer[top_layer_index], as_='gsto_l') if 'gsto_l' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].mean_gsto_sunlit, as_='gsto_l_sunlit') if 'gsto_l_sunlit' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].bulk_gsto_per_layer[top_layer_index], as_='gsto_l_bulk') if 'gsto_l_bulk' in fields else None,
            I(state.canopy_layer_component_pop[top_layer_index][component_index]
              [flag_index].Tleaf_C_estimate, as_='Tleaf_C') if 'Tleaf_C' in fields else None,
            I(state.canopy_component[component_index].td_dd_leaf_pops[flag_index],
                as_='td_dd_flag') if 'td_dd_flag' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].phenology.phenology_stage, as_='leaf_phenology_stage') if 'leaf_phenology_stage' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].POD_Y, as_='pody') if 'pody' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].POD_Y_sunlit, as_='pody_sun') if 'pody_sun' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].POD_0, as_='pod0') if 'pod0' in fields else None,

            __SPACER__("Carbon Allocation >"),
            I(state.canopy_component[component_index].canopy_gsto,
              as_='gsto_canopy') if 'gsto_canopy' in fields else None,
            I(state.canopy_component[component_index].canopy_anet,
              as_='A_n_canopy') if 'A_n_canopy' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].R_d, as_='R_d') if 'R_d' in fields else None,
            I(state.canopy_component[component_index].R_dc,
              as_='R_dc') if 'R_dc' in fields else None,
            I(state.canopy_component[component_index].c_root,
              as_='c_root') if 'c_root' in fields else None,
            I(state.canopy_component[component_index].c_stem,
              as_='c_stem') if 'c_stem' in fields else None,
            I(state.canopy_component[component_index].c_leaf,
              as_='c_leaf') if 'c_leaf' in fields else None,
            I(state.canopy_component[component_index].c_lbrn,
              as_='c_lbrn') if 'c_lbrn' in fields else None,
            I(state.canopy_component[component_index].c_harv,
              as_='c_harv') if 'c_harv' in fields else None,
            I(state.canopy_component[component_index].c_resv,
              as_='c_resv') if 'c_resv' in fields else None,
            I(state.canopy_component[component_index].p_root,
              as_='p_root') if 'p_root' in fields else None,
            I(state.canopy_component[component_index].p_stem,
              as_='p_stem') if 'p_stem' in fields else None,
            I(state.canopy_component[component_index].p_leaf,
              as_='p_leaf') if 'p_leaf' in fields else None,
            I(state.canopy_component[component_index].p_harv,
              as_='p_harv') if 'p_harv' in fields else None,
            I(state.canopy_component[component_index].stem_dm,
              as_='stem_dm') if 'stem_dm' in fields else None,
            I(state.canopy_component[component_index].leaf_dm,
              as_='leaf_dm') if 'leaf_dm' in fields else None,
            I(state.canopy_component[component_index].lbrn_dm,
              as_='lbrn_dm') if 'lbrn_dm' in fields else None,
            I(state.canopy_component[component_index].total_leaf_dm,
              as_='total_leaf_dm') if 'total_leaf_dm' in fields else None,
            I(state.canopy_component[component_index].straw_dm,
              as_='straw_dm') if 'straw_dm' in fields else None,
            I(state.canopy_component[component_index].ear_dm,
              as_='ear_dm') if 'ear_dm' in fields else None,
            I(state.canopy_component[component_index].aboveground_dm,
              as_='aboveground_dm') if 'aboveground_dm' in fields else None,
            I(state.canopy_component[component_index].belowground_dm,
              as_='belowground_dm') if 'belowground_dm' in fields else None,
            I(state.canopy_component[component_index].grain_dm,
              as_='grain_dm') if 'grain_dm' in fields else None,
            I(state.canopy_component[component_index].harvest_index,
              as_='harvest_index') if 'harvest_index' in fields else None,
            I(state.canopy_component[component_index].yield_ha,
              as_='yield_ha') if 'yield_ha' in fields else None,

            I(state.canopy_component[component_index].gpp, as_='gpp') if 'gpp' in fields else None,
            I(state.canopy_component[component_index].npp, as_='npp') if 'npp' in fields else None,
            I(state.canopy_component[component_index].npp_acc,
              as_='npp_acc') if 'npp_acc' in fields else None,
            I(state.canopy_component[component_index].R_pm,
              as_='R_pm') if 'R_pm' in fields else None,
            I(state.canopy_component[component_index].R_pg,
              as_='R_pg') if 'R_pg' in fields else None,

            __SPACER__("Ozone damage Flag leaf >"),
            I(state.canopy_component_population[component_index]
              [flag_index].fO3_h, as_='fO3_h') if 'fO3_h' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].fO3_d, as_='fO3_d') if 'fO3_d' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].fO3_l, as_='fO3_l') if 'fO3_l' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].f_LS, as_='f_LS') if 'f_LS' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].f_LA, as_='f_LA') if 'f_LA' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].t_lep_limited, as_='t_lep_limited') if 't_lep_limited' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].t_lse_limited, as_='t_lse_limited') if 't_lse_limited' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].c_i, as_='c_i') if 'c_i' in fields else None,

            __SPACER__("Soil Moisture Canopy >"),
            I(state.canopy.PM.precip_acc_prev_day, as_='precip_acc') if 'precip_acc' in fields else None,
            I(state.canopy.SMD.SWP, as_='swp') if 'swp' in fields else None,
            I(state.canopy.SMD.ASW, as_='asw') if 'asw' in fields else None,
            I(state.canopy.SMD.SMD, as_='smd') if 'smd' in fields else None,
            I(state.canopy.SMD.Sn, as_='sn') if 'sn' in fields else None,

            I(state.canopy.PM.Ei_hr, as_='ei') if 'ei' in fields else None,
            I(state.canopy.PM.Et_hr, as_='et') if 'et' in fields else None,
            I(state.canopy.PM.PEt_hr, as_='pet') if 'pet' in fields else None,
            I(state.canopy.PM.Es_hr, as_='es') if 'es' in fields else None,

            I(state.canopy.PM.Ei_acc, as_='ei_acc') if 'ei_acc' in fields else None,
            I(state.canopy.PM.Eat_acc, as_='eat_acc') if 'eat_acc' in fields else None,
            I(state.canopy.PM.Et_acc, as_='et_acc') if 'et_acc' in fields else None,
            I(state.canopy.PM.PEt_acc, as_='pet_acc') if 'pet_acc' in fields else None,
            I(state.canopy.PM.Es_acc, as_='es_acc') if 'es_acc' in fields else None,

            __SPACER__("Multiplicative >"),
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.f_phen, as_='f_phen') if 'f_phen' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.leaf_f_phen, as_='leaf_f_phen') if 'leaf_f_phen' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.f_light, as_='f_light') if 'f_light' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.leaf_f_light, as_='leaf_f_light') if 'leaf_f_light' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.f_temp, as_='f_temp') if 'f_temp' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.f_VPD, as_='f_VPD') if 'f_VPD' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.f_SW, as_='f_SW') if 'f_SW' in fields else None,
            I(state.canopy_layer_component[top_layer_index]
              [component_index].gsto_params.f_O3, as_='f_O3') if 'f_O3' in fields else None,

            __SPACER__("Resistance Canopy Top Layer >"),
            I(state.canopy.rmodel_O3.Rsto[top_layer_index],
              as_='rsto_c') if 'rsto_c' in fields else None,
            I(state.canopy.rmodel_O3.Rext[top_layer_index],
              as_='rext') if 'rext' in fields else None,
            #   TODO: Split ra options here
            I(state.canopy.rmodel_O3_ref.Ra_canopy_to_izr, as_='ra_ref') if 'ra' in fields else None,
            I(state.canopy.rmodel_O3.Ra_canopy_to_izr, as_='ra') if 'ra' in fields else None,
            I(state.canopy.rmodel_O3.Ra_canopy_top_to_izr, as_='ra_c_top') if 'ra' in fields else None,
            I(state.canopy.rmodel_O3_ref.Ra_measured_to_izr,
              as_='ra_measurement') if 'ra' in fields else None,
            I(state.canopy.rmodel_O3.Rb, as_='rb') if 'rb' in fields else None,
            I(state.canopy.rmodel_O3.Rsur[top_layer_index],
              as_='rsur') if 'rsur' in fields else None,
            I(state.canopy.rmodel_O3_ref.Rsur[top_layer_index],
              as_='rsur_ref') if 'rsur' in fields else None,
            I(state.canopy.rmodel_O3.Rinc[top_layer_index],
              as_='rinc') if 'rinc' in fields else None,
            I(state.canopy.rmodel_H2O.Rsto[0],
              as_='rsto_h2o') if 'rsto_h2o' in fields else None,

            __SPACER__("V_cmax Flag Top Layer >"),
            I(state.canopy_component_population[component_index]
              [flag_index].V_cmax_25_per_layer[top_layer_index], as_='V_cmax_25') if 'V_cmax_25' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].J_max_25, as_='J_max_25') if 'J_max_25' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].V_cmax, as_='V_cmax') if 'V_cmax' in fields else None,
            I(state.canopy_component_population[component_index]
              [flag_index].J_max, as_='J_max') if 'J_max' in fields else None,

            __SPACER__("Other >"),
            I(state.canopy_layer_component[top_layer_index]
              [component_index].LAIsunfrac, as_='LAIsunfrac_top') if 'LAIsunfrac' in fields else None,
            I(state.canopy_component[component_index].LAI,
              as_='component_LAI') if 'component_LAI' in fields else None,
        ]))),
        # Multi layer/ Multi population outputs
        log_values(lambda state: flatten_list(list(filter(lambda f: f, [
            __SPACER__("Multi layer population >"),
            [I(state.canopy_component[component_index].td_dd_leaf_pops[iP],
               as_=f'td_dd_iP_{iP}') for iP in range(nP)] if 'td_dd_flag' in fields else [],
            [I(state.canopy_component_population[component_index]
                [iP].O3up, as_=f'fst_iP_{iP}') for iP in range(nP)] if 'fst' in fields else [],
            [I(state.canopy_component_population[component_index]
                [iP].O3up_acc, as_=f'fst_acc_iP_{iP}') for iP in range(nP)] if 'fst_acc' in fields else [],
            [I(state.canopy_component_population[0][iP].f_LS,
                as_=f'f_LS_iP_{iP}') for iP in range(nP)] if 'f_LS' in fields else [],
            [I(state.canopy_component_population[0][iP].fO3_d,
                as_=f"fO3_d_iP_{iP}") for iP in range(nP)] if 'fO3_d' in fields else [],
            [I(state.canopy_component_population[0][iP].fO3_h,
                as_=f"fO3_h_iP_{iP}") for iP in range(nP)] if 'fO3_h' in fields else [],
            [I(state.canopy_component_population[component_index]
                [iP].V_cmax, as_=f'V_cmax_iP_{iP}') for iP in range(nP)] if 'V_cmax' in fields else [],
            [I(state.canopy_component_population[component_index]
                [iP].J_max, as_=f'J_max_iP_{iP}') for iP in range(nP)] if 'J_max' in fields else [],
            [I(state.canopy_component_population[0][iP].mean_gsto_per_layer[iL],
                as_=f"gsto_l_iL_{iL}_iP_{iP}") for iP in range(nP) for iL in range(nL)] if 'gsto_l' in fields else [],
            [I(state.canopy_component_population[0][iP].A_n,
                as_=f"A_n_l_iP_{iP}") for iP in range(nP)] if 'A_n' in fields else [],
            [I(state.canopy_component_population[0][iP].O3up,
                as_=f"O3up_iP_{iP}") for iP in range(nP)] if 'O3up' in fields else [],
            [I(state.canopy_component_population[0][iP].fLAI_layer[iL],
                as_=f"fLAI_iL_{iL}_iP_{iP}") for iP in range(nP) for iL in range(nL)] if 'fLAI' in fields else [],
            [I(state.canopy_component_population[0][iP].LAI,
                as_=f"LAI_iP_{iP}") for iP in range(nP)] if 'lai' in fields else [],
            [I(state.canopy_layer_component[iL][0].LAIsunfrac,
                as_=f"LAIsunfrac_iL_{iL}") for iL in range(nL)] if 'LAIsunfrac' in fields else [],
            # [I(state.canopy_component_population[0][iP].LAI_brown, # NOTE: We cannot currently get brown LAI per population
            #     as_=f"LAI_brown_iP_{iP}") for iP in range(nP)] if 'lai_brown' in fields else [],
            [I(state.canopy_component_population[component_index]
               [iP].fO3_l, as_=f"fO3_l_iP_{iP}") for iP in range(nP)] if 'fO3_l' in fields else [],
            [I(state.canopy_component[0].leaf_pop_distribution[iL][iP],
                as_=f"leaf_pop_distribution_iL_{iL}_iP_{iP}") for iL in range(nL) for iP in range(nP)] if 'leaf_pop_distribution' in fields else [],
            [I(state.canopy_component[0].growing_populations[iP],
                as_=f"growing_populations_iP_{iP}") for iP in range(nP)] if 'growing_populations' in fields else [],
            [I(state.canopy_layer_component[iL][0].LAI,
                as_=f"layer_lai_iL_{iL}") for iL in range(nL)] if 'layer_lai' in fields else [],
            [I(state.canopy_layer_component[iL][0].LAI_brown,
                as_=f"layer_lai_brown_iL_{iL}") for iL in range(nL)] if 'layer_lai_brown' in fields else [],
            [I(state.canopy_layers[iL].micro_met.PARsun, as_=f'PARsun_iL_{iL}') for iL in range(
                nL)] if 'PARsun' in fields else [],
            [I(state.canopy_layers[iL].micro_met.PARshade,
                as_=f'PARshade_iL_{iL}') for iL in range(nL)] if 'PARshade' in fields else [],
            [I(state.canopy_layers[iL].micro_met.micro_O3,
                as_=f'micro_O3_iL_{iL}') for iL in range(nL)] if 'micro_O3' in fields else [],
            [I(state.canopy_layers[iL].micro_met.micro_u,
                as_=f'micro_u_iL_{iL}') for iL in range(nL)] if 'micro_u' in fields else [],
            [I(state.canopy_component_population[0][iP].leaf_rmodel_O3.Rsto[iL],
                as_=f'rsto_l_iL_{iL}_iP_{iP}') for iL in range(nL) for iP in range(nP)] if 'rsto_l' in fields else [],
            [I(state.canopy_component_population[0][iP].leaf_rmodel_O3.Rb[iL],
                as_=f'rb_l_iL_{iL}_iP_{iP}') for iL in range(nL) for iP in range(nP)] if 'rb_l' in fields else [],
            [I(state.canopy_layers[iL].layer_height,
                as_=f'layer_height_iL_{iL}') for iL in range(nL)] if 'layer_height' in fields else [],
            [I(state.canopy_component_population[component_index]
               [iP].phenology.phenology_stage, as_=f'leaf_phenology_stage_{iP}') if 'leaf_phenology_stage' in fields else None
             for iP in range(nP)] if 'leaf_phenology_stage' in fields else [],
        ])))) if log_multilayer else [],
        # Debug outputs
        log_values(lambda state: flatten_list(list(filter(lambda f: f, [
            __SPACER__("Debug >"),
            I(state.debug.ewert_loop_iterations, as_='ewert_loop_iterations'),
        ])))) if "_debug" in fields else [],
    ]


def hourly_processes(config: Config_Shape, hr: int, run_dir: Path = None, dump_state_n: int = None) -> List[Process]:
    """Take the current hour and returns a list of processes.

    to be ran in this hour
    """
    nL = config.Land_Cover.nL
    nLC = config.Land_Cover.nLC
    nP = config.Land_Cover.nP

    # Options
    LAI_method = config.Land_Cover.LAI_method
    SAI_method = config.Land_Cover.SAI_method
    primary_SAI_method = config.Land_Cover.parameters[0].phenology.SAI_method
    LAI_distribution_method = config.Land_Cover.LAI_distribution_method
    thermal_time_method = config.Met.thermal_time_method
    soil_moisture_source = config.soil_moisture.source
    ra_calc_method = config.resistance.ra_calc_method
    rsur_calc_method = config.resistance.rsur_calc_method
    # Per component options
    f_light_method = [
        config.Land_Cover.parameters[iLC].multip_gsto.f_light_method for iLC in range(nLC)]
    f_temp_method = [
        config.Land_Cover.parameters[iLC].multip_gsto.f_temp_method for iLC in range(nLC)]
    f_VPD_method = [
        config.Land_Cover.parameters[iLC].gsto.f_VPD_method for iLC in range(nLC)]
    f_SW_method = [config.Land_Cover.parameters[iLC].gsto.f_SW_method for iLC in range(nLC)]
    f_O3_method = [config.Land_Cover.parameters[iLC].multip_gsto.f_O3_method for iLC in range(nLC)]
    V_J_method = [config.Land_Cover.parameters[iLC].pn_gsto.V_J_method for iLC in range(nLC)]
    Tleaf_method = [config.Land_Cover.parameters[iLC].gsto.Tleaf_method for iLC in range(nLC)]
    leaf_f_phen_Anet_influence = [
        config.Land_Cover.parameters[iLC].pn_gsto.leaf_f_phen_Anet_influence for iLC in range(nLC)]
    gsto_method = [config.Land_Cover.parameters[iLC].gsto.method for iLC in range(nLC)]
    parsun_shade_method = config.Met.PARsunshade_method
    canopy_ozone_method = config.Land_Cover.ozone_deposition_method
    use_carbon_allocation = config.carbon_allocation.use_carbon_allocation
    senescence_method = [
        config.Land_Cover.parameters[iLC].pn_gsto.senescence_method for iLC in range(nLC)]
    have_precip_data = config.Met.inputs.precip_method != InputMethod.SKIP
    have_Hd_data = config.Met.inputs.Hd_method != InputMethod.SKIP
    have_ustar_ref_data = config.Met.inputs.ustar_ref_method != InputMethod.SKIP
    have_ustar_data = config.Met.inputs.ustar_method != InputMethod.SKIP
    is_OTC = config.Location.OTC
    sparse_data = config.Met.sparse_data
    output_fields = config.output.fields
    log_multilayer = config.output.log_multilayer
    ewert_loop_method = [config.Land_Cover.parameters[iLC].pn_gsto.ewert_loop_method for iLC in range(nLC)]

    return [
        tag_process(f"===== Start of Hourly Processes (Hour: {hr}) ====="),
        set_hour(hr),
        # dump_state_process(), ## Turn on to debug state at specific hours
        sync_row_index(),  # External data row pointer moved to hr = 0

        [set_thermal_time_process(thermal_time_method, iLC) for iLC in range(nLC)]
        if thermal_time_method == ThermalTimeMethods.EXTERNAL else [],

        # TODO: It would be great to not rely on knowing the hour here.
        daily_start_processes(config, run_dir, dump_state_n) if (
            hr == 0 or sparse_data) else [],

        [
            calc_canopy_LAI_processes(nLC, LAI_method),
            distribute_lai_per_layer_processes(nL, nLC, LAI_distribution_method, LAI_method),
            calc_canopy_SAI(nL, nLC, SAI_method, primary_SAI_method),
            distribute_lai_to_leaf_populations_processes(nL, nP, nLC, LAI_method),
        ] if LAI_method == LAIMethods.INPUT_HOURLY else [],


        calc_soil_moisture_hourly_process(soil_moisture_source)
        if soil_moisture_source != SoilMoistureSource.P_M else [],

        # MET PROCESSES
        tag_process("=== Met Processes ==="),
        accumulate_hourly_temperature_process(),
        accumulate_precipitation_process() if have_precip_data else [],

        # MICRO MET
        # Solar
        [calc_PAR_sun_shade_process(parsun_shade_method, iL, iLC, nL, nLC)
         for iL in range(nL) for iLC in range(nLC)],
        MLMC_sunlit_LAI_process(nL, nLC),

        # Wind
        calc_ustar_ref_process(have_Hd_data, have_ustar_ref_data, have_ustar_data, is_OTC),
        calc_windspeed_parameters_process(have_ustar_data),
        # TODO: below resets windspeed. Check this is correct
        [calc_layer_windspeed_process(iL) for iL in range(nL)],


        # TLEAF
        # NOTE: This uses calculated gsto from previous hour
        [calc_leaf_temperature_process(iL, iLC, iP, Tleaf_method[iLC])
         for iL in range(nL) for iLC in range(nLC) for iP in range(nP)],

        # GSTO SETUP
        tag_process("=== Gsto Setup Processes ==="),
        [set_vcmax_and_jmax(nL, iLC, iP, V_J_method[iLC])
         for iLC in range(nLC) for iP in range(nP)],
        [calc_leaf_f_phen_effect_on_V_cmax_25_process(
            iL, iLC, iP, leaf_f_phen_Anet_influence[iLC])
            if leaf_f_phen_Anet_influence[iLC] == LeafFPhenAnetInfluence.V_C_MAX else []
            for iL in range(nL) for iLC in range(nLC) for iP in range(nP)],
        [[
            # TODO: Check how many of these could be calculated daily instead
            calc_f_light_process(iL, iLC) if f_light_method[iLC] == "enabled" else [],
            calc_f_temp_process(iL, iLC, f_temp_method[iLC]
                                ) if f_temp_method[iLC] != "disabled" else [],
            calc_f_VPD_process(iL, iLC, f_VPD_method[iLC]),
            f_SW_process(iL, iLC, f_SW_method[iLC]),
            f_O3_process(iL, iLC, f_O3_method[iLC], nP=nP),
            calc_g_bv_process(iLC),
        ]for iL in range(nL) for iLC in range(nLC)],

        # GSTO CALC
        tag_process("=== Gsto Processes ==="),
        # = MULTIPLICATIVE
        [[
            [gsto_multiplicative_process(iL, iLC, iP)for iL in range(nL)],
            [scale_layer_mean_gsto_to_layer_bulk_gsto_process(iL, iLC)for iL in range(nL)],
            scale_layer_bulk_gsto_to_canopy_gsto_process(nL, iLC),
            scale_leaf_layer_mean_gsto_to_leaf_layer_bulk_gsto_process(iLC, iP, nL),
        ] for iLC in range(nLC) for iP in range(nP) if gsto_method[iLC] == "multiplicative"],

        # = OZONE DAMAGE
        [[ozone_damage_processes(iLC, iP, iP == nP - 1, senescence_method[iLC])
          for iP in range(nP)]for iLC in range(nLC)],

        # = PHOTOSYNTHESIS
        [[
            calc_D_0_process(iLC),
            [ewert_leaf_process(iLC, iP, nL, ewert_loop_method[iLC])for iP in range(nP)],
            [convert_gsto_CO2umol_to_O3mmol_process(iLC, iP, nL) for iP in range(nP)],
            calc_layer_mean_gsto_process(iLC, nP, nL),
            [scale_layer_mean_gsto_to_layer_bulk_gsto_process(iL, iLC) for iL in range(nL)],
            [scale_leaf_layer_mean_gsto_to_leaf_layer_bulk_gsto_process(
                iLC, iP, nL) for iP in range(nP)],
            [limit_gsto_l_with_leaf_fphen_process(
                iL, iLC, iP, leaf_f_phen_Anet_influence[iLC])
                if leaf_f_phen_Anet_influence[iLC] == LeafFPhenAnetInfluence.G_STO else [] for iL in range(nL) for iP in range(nP)]
        ] for iLC in range(nLC) if gsto_method[iLC] == "photosynthesis"],

        # OZONE DEPOSITION,
        tag_process("=== Ozone Deposition Processes ==="),
        calc_resistance_model_process(nL, nLC, ra_calc_method,
                                      rsur_calc_method, is_OTC),
        calc_deposition_velocity_process() if not is_OTC else [],
        calc_O3_otc_process(nL) if is_OTC else [],
        calc_canopy_ozone_concentration_process(canopy_ozone_method, nL) if not is_OTC else [],
        calc_multi_layer_O3_ozone_concentration_process(nL) if not is_OTC else [],
        # calc_multi_layer_O3_ozone_concentration_process(nL) if not OTC and nL > 1 else [],
        [O3_ppb_to_nmol_process(iL) for iL in range(nL)],


        # OZONE DOSE
        tag_process("=== Ozone Dose Processes ==="),
        [calc_leaf_resistance_model_process(iLC, iP, nL) for iLC in range(nLC) for iP in range(nP)],
        [calc_fst_leaf_process(iLC, iP, nL) for iLC in range(nLC) for iP in range(nP)],
        # [calc_fst_leaf_acc_process(iLC, iP) for iLC in range(nLC) for iP in range(nP)],
        [calc_fst_leaf_acc_hour_process(iLC, iP) for iLC in range(nLC) for iP in range(nP)],
        [calc_POD_leaf_process(iLC, iP) for iLC in range(nLC) for iP in range(nP)],
        [calc_OT_leaf_process(iLC, iP, nL, nP) for iLC in range(nLC) for iP in range(nP)],

        # SOIL MOISTURE
        tag_process("=== Soil Moisture Processes ==="),
        # Calculate the inputs to soil moisture calculations
        [
            calc_h2o_r_model_process(),
            check_soil_evaporation_blocked_process(),
            penman_monteith_hourly_process(gsto_method[0]),
        ] if soil_moisture_source == SoilMoistureSource.P_M else [],

        # CARBON ALLOCATION
        tag_process("=== Carbon Allocation Processes ==="),
        [[
            # Possibly not valid for multiplicative?
            scale_layer_bulk_gsto_to_canopy_gsto_process(nL, iLC),
            scale_rd_to_canopy_rd_process(iLC, nP),
            scale_anet_to_canopy_anet_process(iLC, nP),
            scale_anet_to_canopy_an_gross_process(iLC, nP),
            calculate_canopy_npp_process(iLC) if use_carbon_allocation else [],
            accumulate_canopy_npp_process(iLC) if use_carbon_allocation else [],
        ] for iLC in range(nLC) if gsto_method[iLC] == "photosynthesis"],

        # LOGGING
        tag_process("== Hourly Logging Processes =="),
        log_processes(nL, nLC, nP, output_fields, log_multilayer),
        store_prev_state(),  # Used when we require the data from previous hour.
        # At end of hour advance time step
        advance_time_step_process(),
        tag_process("===== End of Hourly Processes ====="),
    ]


# def daily_end_processes(config: Config_Shape) -> List[Process]:
#     """Take the current day and returns a list of processes to.

#     # NOTE: Will not currently work because we are setting timestep in hour_processes
#     # so this will be hour ahead

#     be ran at the end of the day
#     """
#     return [

#     ]

def get_row_processes_hourly(
    config: Config_Shape,
    hours: List[int],
    run_dir: Path = None,
    dump_state_n: int = None,
) -> List[Process]:
    """Get the model processes for each hour in hours.

    Parameters
    ----------
    config : Config_Shape
        _description_
    hours : List[int]
        _description_
    run_dir : Path, optional
        If this and dump_state_n provided will
        periodically dump state to file, by default None
    dump_state_n : int, optional
        _description_, by default None

    Returns
    -------
    List[Process]
        _description_
    """
    return flatten_list([hourly_processes(config, hr, run_dir, dump_state_n) for hr in hours])


def full_model_processes(
    config: Config_Shape,
    hours: List[int],
    run_validation: bool = False,
    run_dir: Path = None,
    dump_state_n: int = None,
) -> List[Process]:
    """Get a flattened list of all processes to be ran between dates inclusive."""
    row_processes = get_row_processes_hourly(config, hours, run_dir, dump_state_n)
    all_processes = [
        setup_validation_processes(config) if run_validation else [],
        row_processes,
    ]
    return flatten_list(all_processes)

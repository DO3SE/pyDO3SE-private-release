# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import Snapshot


snapshots = Snapshot()

snapshots['TestDailyCarbonAllocation.test_growing_season 1'] = {
    'c_harv': 12.188050446699801,
    'c_leaf': 0.01914809398450125,
    'c_resv': 0.00030701179772200623,
    'c_root': 0.3082767031294814,
    'c_stem': 2.602319102580931,
    'lai': 1.0078714081972373,
    'net_prod_acc': 0,
    'plant_height': 2.708221950800548
}

snapshots['TestCalcPartitionCoefficients.test_output_values p_root'] = [
    0.8803682021208932,
    0.7732500373560889,
    0.5922010669465662,
    0.3044854774129787,
    0.049061779666647355,
    1.0114195459401717e-05,
    4.5990515265262127e-10
]

snapshots['TestCalcPartitionCoefficients.test_output_values p_leaf'] = [
    0.11914487998653855,
    0.22153984487484066,
    0.3591881038176261,
    0.3909670920105659,
    0.13336374413970897,
    5.820317625631814e-05,
    5.602791744686968e-09
]

snapshots['TestCalcPartitionCoefficients.test_output_values p_stem'] = [
    0.00048691789256821483,
    0.005210117768746185,
    0.04861082376538027,
    0.3044854774129787,
    0.5976948344902958,
    0.001501079699910745,
    8.315280226228297e-07
]

snapshots['TestCalcPartitionCoefficients.test_output_values p_harv'] = [
    1.6762036988705228e-17,
    3.242852649505786e-13,
    5.470427545684183e-09,
    6.195316347664611e-05,
    0.21987964170334776,
    0.9984306029283736,
    0.9999991624092804
]

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0-0] DVI 0 c_leaf 0'] = 0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0-0.5] DVI 0.5 c_leaf 0'] = 0.0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0-1.3] DVI 1.3 c_leaf 0'] = 0.0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0-1.8] DVI 1.8 c_leaf 0'] = 0.0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0-2] DVI 2 c_leaf 0'] = 0.0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.001-0] DVI 0 c_leaf 0.001'] = 0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.001-0.5] DVI 0.5 c_leaf 0.001'] = 0.05622888992707347

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.001-1.3] DVI 1.3 c_leaf 0.001'] = 0.05375541512654959

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.001-1.8] DVI 1.8 c_leaf 0.001'] = 0.05290885175472733

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.001-2] DVI 2 c_leaf 0.001'] = 0.05263559960657303

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.01-0] DVI 0 c_leaf 0.01'] = 0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.01-0.5] DVI 0.5 c_leaf 0.01'] = 0.5622888992707347

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.01-1.3] DVI 1.3 c_leaf 0.01'] = 0.537554151265496

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.01-1.8] DVI 1.8 c_leaf 0.01'] = 0.5290885175472733

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.01-2] DVI 2 c_leaf 0.01'] = 0.5263559960657304

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.1-0] DVI 0 c_leaf 0.1'] = 0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.1-0.5] DVI 0.5 c_leaf 0.1'] = 5.622888992707347

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.1-1.3] DVI 1.3 c_leaf 0.1'] = 5.3755415126549595

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.1-1.8] DVI 1.8 c_leaf 0.1'] = 5.2908851754727335

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[0.1-2] DVI 2 c_leaf 0.1'] = 5.263559960657304

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[1-0] DVI 0 c_leaf 1'] = 0

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[1-0.5] DVI 0.5 c_leaf 1'] = 56.22888992707347

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[1-1.3] DVI 1.3 c_leaf 1'] = 53.75541512654959

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[1-1.8] DVI 1.8 c_leaf 1'] = 52.90885175472733

snapshots['TestCalcLAIFromDVIAndCarbon.test_output_values[1-2] DVI 2 c_leaf 1'] = 52.63559960657303

snapshots['TestGetPlantHeightFromCarbon.test_output_values[0] c_stem 0'] = 0.0

snapshots['TestGetPlantHeightFromCarbon.test_output_values[0.001] c_stem 0.001'] = 0.11655744903626222

snapshots['TestGetPlantHeightFromCarbon.test_output_values[0.01] c_stem 0.01'] = 0.2927790747255565

snapshots['TestGetPlantHeightFromCarbon.test_output_values[0.1] c_stem 0.1'] = 0.7354277852330547

snapshots['TestGetPlantHeightFromCarbon.test_output_values[1] c_stem 1'] = 1.8473110750820518

snapshots['TestCalcRootFractionFromCarbon.test_output_values[0] c_root 0'] = 0.6321205588285577

snapshots['TestCalcRootFractionFromCarbon.test_output_values[0.001] c_root 0.001'] = 0.6321205588285577

snapshots['TestCalcRootFractionFromCarbon.test_output_values[0.01] c_root 0.01'] = 0.6321205588285577

snapshots['TestCalcRootFractionFromCarbon.test_output_values[0.1] c_root 0.1'] = 0.6321205588285577

snapshots['TestCalcRootFractionFromCarbon.test_output_values[1] c_root 1'] = 0.6321205588285577

# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots['test_co2_concentration_in_stomata_iteration 1'] = GenericRepr("CO2_loop_State(c_i=1159.7695868330525, c_i_diff=2319.539173666105, g_sto=20000.0, A_n=-37.865, A_c=-8.495585110107601, A_p=59.5, A_j=-37.545, A_n_limit_factor='A_j', f_VPD=0.9994183666165483, iterations=1)")

snapshots['test_co2_concentration_in_stomata_iteration_b 1'] = GenericRepr("CO2_loop_State(c_i=366.3646469481298, c_i_diff=98.2093515380468, g_sto=965039.7007066957, A_n=37.465715119511366, A_c=37.785715119511366, A_p=59.5, A_j=59.67014221159362, A_n_limit_factor='A_c', f_VPD=0.9996640316746571, iterations=5)")

snapshots['TestEwertLeafPop.test_should_run_with_multiple_layers 1'] = {
    'A_c': 0.3348967506706674,
    'A_j': 0.5462816729472586,
    'A_n': 0.31710872053521166,
    'A_n_limit_factor': [
        'A_c',
        'A_c',
        'A_c'
    ],
    'A_p': 0.5929343378485238,
    'R_d': 0.01778803013545571,
    'c_i': 3.038712406687719,
    'f_VPD': 0.7038506876305467,
    'g_sv': None,
    'g_sv_per_layer': [
        590898.6968433901,
        590898.6968433901,
        590898.6968433901
    ],
    'j_max': 308.5696703266116,
    'v_cmax': 118.58686756970474
}

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_n - min/max'] = 'max: 117.6017801845255, min: -15.213566680445922'

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_c - min/max'] = 'max: 608.0930025162996, min: 0.0'

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_j - min/max'] = 'max: 158.2114934526396, min: 0.0'

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_p - min/max'] = 'max: 906.6364836919419, min: 0.0'

snapshots['TestEwertLeafPop.test_gsto_output_full_year R_d - min/max'] = 'max: 27.199094510758254, min: 0.0'

snapshots['TestEwertLeafPop.test_gsto_output_full_year g_sto - min/max'] = 'max: 20.86, min: 20.86'

snapshots['TestEwertLeafPop.test_gsto_output_full_year g_sv - min/max'] = 'max: 20000, min: 20000'

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_n'] = '''[ 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
  0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,59.77653, 0.     ,
  0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 39.25809, 0.     , 0.     , 0.     , 0.     , 0.     ]'''

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_c'] = '''[  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ,  0.     ,  0.     ,  0.     ,364.41566,  0.     ,
   0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ]'''

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_j'] = '''[  0.     ,  0.     ,  0.     ,  0.     ,  2.2238 ,  0.     ,122.91599,
   0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ]'''

snapshots['TestEwertLeafPop.test_gsto_output_full_year A_p'] = '''[ 21.18616,432.63468,  0.     ,  0.     ,  0.     ,  0.     ,316.77606,
   0.     , 36.45572,380.21309,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,538.20108,  0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
   0.     ,  0.     ]'''

snapshots['TestEwertLeafPop.test_gsto_output_full_year R_d'] = '''[0.     ,0.     ,0.     ,0.     ,8.81573,0.     ,0.     ,0.68939,0.     ,
 0.     ,0.     ,0.     ,0.     ,0.     ,0.     ,0.     ,0.     ,0.     ,
 0.     ,0.     ,0.     ,0.     ,0.     ,0.     ,0.     ,8.32334,3.55093,
 0.     ,0.     ,7.14498]'''

snapshots['TestEwertLeafPop.test_gsto_output_full_year g_sto'] = '''[20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,
 20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,20.86,
 20.86,20.86,20.86,20.86,20.86,20.86]'''

snapshots['TestEwertLeafPop.test_gsto_output_full_year g_sv'] = '''[20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,
 20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,
 20000,20000,20000,20000,20000,20000]'''

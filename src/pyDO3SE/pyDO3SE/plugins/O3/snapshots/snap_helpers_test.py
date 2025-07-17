# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots['test_calc_resistance_model 1'] = GenericRepr('Resistance_Model(nL=2, Ra_c=11.027459771245551, Ra=8.59151234626083, Rb=5.520331050205918, Rinc=[335.99999999999994, 335.99999999999994], Rext=[833.3333333333334, 500.0], Rsto=[24600.000000000004, 15769.23076923077], Rgs=123, Rsur=[761.9778956258517, 437.099278418627], Rtotal=[286.4489646743288, 256.78112826402116])')

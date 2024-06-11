# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots['TestCalcResistanceModel.test_works_without_error 1'] = GenericRepr('Resistance_Model(nL=2, Ra_measured_to_izr=4.046121837677561, Ra_canopy_to_izr=11.027459771245551, Ra_canopy_top_to_izr=8.59151234626083, Rb=5.520331050205918, Rinc=[335.99999999999994, 335.99999999999994], Rext=[833.3333333333334, 500.0], Rsto=[24600.000000000004, 15769.23076923077], Rgs=123, Rsur=[811.5491646019753, 490.1539007901587], Rtotal=[293.1811510568345, 275.51678559987187])')

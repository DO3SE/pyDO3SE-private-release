from do3se_phenology.geo import estimate_latitude_SGS_EGS


class TestEstimateLatitudeSGSEGS():

    def setup(self):
        pass

    def test_should_return_correct_values(self):
        sgs, egs = estimate_latitude_SGS_EGS(10, 10)
        assert sgs == 45
        assert egs == 376

    def test_should_vary_by_latitude(self):
        sgs1, egs1 = estimate_latitude_SGS_EGS(10, 10)
        sgs2, egs2 = estimate_latitude_SGS_EGS(20, 10)

        assert sgs1 != sgs2
        assert sgs1 < sgs2

        assert egs1 != egs2
        assert egs1 > egs2

    def test_should_vary_by_elev(self):
        sgs1, egs1 = estimate_latitude_SGS_EGS(10, 10)
        sgs2, egs2 = estimate_latitude_SGS_EGS(10, 1000)

        assert sgs1 != sgs2
        assert sgs1 < sgs2

        assert egs1 != egs2
        assert egs1 > egs2

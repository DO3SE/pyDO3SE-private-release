import pytest
from pyDO3SE.plugins.leaf_temperature.de_boeck import get_leaf_temp_de_boeck


class TestGetLeafTempDeBoeck:

    def default_args(self):
        return {
            "R": 200,
            "eact": 0.93,
            "T_air": 20,
            "initial_T_leaf": 20,
            "P": 101.18,
            "u_speed": 0.8,
            "g_vs": 200000,  # Check units and where this comes from
            "hypostomatous": True,
            "d": 0.02,
            "albedo": 0.2,
            "cloud_cover": 1.0,
            "balance_threshold": 0.0010000000474974513,
            "adjustment_factor": 0.019999999552965164,
            "max_iterations": 50,
        }

    def default_run(self, kwargs={}):
        default_args = self.default_args()
        args = {
            **default_args,
            **kwargs,
        }
        return get_leaf_temp_de_boeck(
            **args,
        )

    def test_works_without_error(self):
        t_leaf = self.default_run()
        assert t_leaf

    @pytest.mark.parametrize('R', [250, 500])
    def test_should_increase_leaf_temperature_when_R_is_high(self, R):
        t_leaf = self.default_run(dict(R=R))
        assert t_leaf > self.default_args()['T_air']


    @pytest.mark.parametrize('R', [0, 100, 200])
    def test_should_decrease_leaf_temperature_when_R_is_low(self, R):
        t_leaf = self.default_run(dict(R=R))
        assert t_leaf < self.default_args()['T_air']

    def test_increasing_R_increases_leaf_temp(self):
        R = [10, 50, 100, 800]
        last_t_leaf = 0
        for r in R:
            next_t_leaf = self.default_run(dict(R=r))
            assert next_t_leaf > last_t_leaf
            last_t_leaf = next_t_leaf

    def test_increasing_pressure_increases_leaf_temp(self):
        P = [100, 101, 102, 103]
        last_t_leaf = 0
        for p in P:
            next_t_leaf = self.default_run(dict(P=p * 1e3))
            assert next_t_leaf > last_t_leaf
            last_t_leaf = next_t_leaf

    def test_decreasing_gsto_increases_leaf_temp(self):
        GSTO = [4, 3, 2, 1]
        last_t_leaf = 0
        for g in GSTO:
            next_t_leaf = self.default_run(dict(g_vs=g * 1e-6))
            assert next_t_leaf > last_t_leaf
            last_t_leaf = next_t_leaf

    # def test_decreasing_u_increases_leaf_temp(self):
    #     UH = reversed([0.1, 2.9, 4, 5.5, 5.8])
    #     last_t_leaf = 0
    #     for u in UH:
    #         if u < 5.5:
    #             print(u)
    #             next_t_leaf = self.default_run(dict(u_speed=u))
    #             assert next_t_leaf > last_t_leaf
    #             last_t_leaf = next_t_leaf
    #         else:
    #             # NOTE: u values higher than 5.8 have an odd effect on t_leaf
    #             pass

    #     next_t_leaf = self.default_run(dict(u_speed=u))


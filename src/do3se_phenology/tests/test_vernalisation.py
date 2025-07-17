import pytest
from math import isclose
import pandas as pd
from do3se_phenology.state import PhenologyStage
from do3se_phenology.utils import get_day_from_td

from do3se_phenology.vernalisation import calc_vernalised_thermal_time_range, calculate_vernalisation_factor



class TestCalculateVernalisationFactor:

    @pytest.mark.parametrize(['max_temp', 'min_temp', 'V_acc_prev', 'v_T_max', 'v_T_min', 'PIV', 'expected_out'], [
        # (99, 98, 0, 30, 15, 1.5, [0.5759, 0]),
        (19.51, 9.67, 0, 30, 15, 2.5, [0.3068, 0.2648, 0.2648, 0]),  # dd 320 hr 0 PIV 2.5
        (14.05, 3.27, 0, 30, 15, 1.5, [0.58207, 0.726251, 0.7262, 0]),  # dd 331 hr 0 PIV 1.5
    ])
    def test_outputs_correct_values(
        self,
        max_temp,
        min_temp,
        V_acc_prev,
        v_T_max,
        v_T_min,
        PIV,
        expected_out,
    ):
        Vf, V_acc, v_pos, v_neg = calculate_vernalisation_factor(
            phenology_stage=PhenologyStage.EMERGED,
            max_ambient_temp=max_temp,
            min_ambient_temp=min_temp,
            V_acc_prev=V_acc_prev,
            v_T_max=v_T_max,
            v_T_min=v_T_min,
            PIV=PIV,
        )

        assert isclose(Vf, expected_out[0], abs_tol=1e-3)
        assert isclose(V_acc, expected_out[1], abs_tol=1e-3)
        assert isclose(v_pos, expected_out[2], abs_tol=1e-3)
        assert isclose(v_neg, expected_out[3], abs_tol=1e-3)

    def test_outputs_correct_values_b(self):

        Vf, V_acc, V_pos, V_neg = calculate_vernalisation_factor(
            phenology_stage=PhenologyStage.EMERGED,
            max_ambient_temp=11.68,
            min_ambient_temp=0.04,
            V_acc_prev=0,
            v_T_max=30,
            v_T_min=15,
            PIV=1.5,
        )

        assert isclose(Vf, 0.5832, abs_tol=1e-3)
        assert isclose(V_acc, 0.8674, abs_tol=1e-3)
        assert isclose(V_pos, 0.8674, abs_tol=1e-3)
        assert isclose(V_neg, 0.0, abs_tol=1e-3)

    def test_outputs_correct_values_c(self):
        """Matching PP spreadsheet output."""
        Vf, V_acc, V_pos, V_neg = calculate_vernalisation_factor(
            phenology_stage=PhenologyStage.EMERGED,
            max_ambient_temp=11.68,
            min_ambient_temp=0.04,
            V_acc_prev=0,
            v_T_max=30,
            v_T_min=15,
            PIV=1.5,
        )

        assert isclose(Vf, 0.5832, abs_tol=1e-3)
        assert isclose(V_acc, 0.8674, abs_tol=1e-3)
        assert isclose(V_pos, 0.8674, abs_tol=1e-3)
        assert isclose(V_neg, 0.0, abs_tol=1e-3)

    def test_increasing_PIV_increases_Vf_factor(self):
        PIV_vals = [0, 1, 2, 3]
        Vf = 1
        V_acc = 0
        for PIV in PIV_vals:
            Vf_out, V_acc, v_pos, v_neg = calculate_vernalisation_factor(
                phenology_stage=PhenologyStage.EMERGED,
                max_ambient_temp=99.9,
                min_ambient_temp=99.9,
                V_acc_prev=0,
                v_T_max=30,
                v_T_min=15,
                PIV=PIV,
            )
            assert Vf_out < Vf
            Vf = Vf_out


class TestVernalisationTdRange:

    def _test_data(self):
        return pd.read_csv('examples/vernalisation/vernalisation_china_xiaoji2008_data_b.csv')

    def _get_args(self, **kwargs):
        test_data = self._test_data()
        t_emerg = 100
        t_flower = 1000
        v_T_max = 30
        v_T_min = 15
        PIV = 1.5

        td_data = test_data['Td'].values

        default_args = dict(
            td_data=td_data,
            hr_data=test_data['Hour'],
            dd_data=test_data['Day'],
            max_ambient_temperatures=test_data['Maximum temperature (Tmax)'],
            min_ambient_temperatures=test_data['Minimum temperature)Tmin'],
            t_emerg=t_emerg,
            t_flower=t_flower,
            v_T_max=v_T_max,
            v_T_min=v_T_min,
            PIV=PIV,
        )
        _kwargs = {
            **default_args,
            **kwargs,
        }
        return _kwargs

    def _default_run(self, **kwargs):
        return calc_vernalised_thermal_time_range(**kwargs)


    def test_should_extend_growing_season(self):
        """Vernalisation should extend the length of the growing season.

        By slowing down thermal time when temperatures are low.
        """

        kwargs = self._get_args()
        td_data = kwargs['td_data']

        td_f_values = self._default_run(**kwargs)

        # Assert final vernalised thermal time is lower than thermal time
        assert td_data[-1] > td_f_values[-1]

    def test_should_have_same_thermal_time_pre_emergence(self):
        # Assert vernalised thermal time is the same as thermal time until emergence
        kwargs = self._get_args()
        td_data = kwargs['td_data']
        test_data = self._test_data()
        td_f_values = self._default_run(**kwargs)

        emerg_row, _ = get_day_from_td(kwargs['t_emerg'], test_data['Day'], td_f_values)

        assert td_data[0] == td_f_values[0]
        assert td_data[emerg_row - 1] == td_f_values[emerg_row - 1]

    def test_should_have_same_rate_of_change_after_flowering(self):
        # Assert the difference in vernalised thermal time and thermal time is the same between flowering and harvest.
        kwargs = self._get_args()
        td_data = kwargs['td_data']
        test_data = self._test_data()
        td_f_values = self._default_run(**kwargs)

        flowering_row, _ = get_day_from_td(kwargs['t_flower'], test_data['Day'], td_f_values)

        assert td_data[flowering_row] - td_f_values[flowering_row] == td_data[-1] - td_f_values[-1]

    def test_can_use_photoperiod_factor(self):
        # Without PID
        kwargs = self._get_args()
        td_data = kwargs['td_data']

        td_f_values_no_pid = self._default_run(**kwargs)


        # With PID
        PID=40
        lat = 1.0
        kwargs = self._get_args(PID=PID, lat=lat)
        td_data = kwargs['td_data']

        td_f_values_with_pid = self._default_run(**kwargs)

        # Assert final vernalised thermal time is lower than thermal time
        assert td_data[-1] > td_f_values_no_pid[-1] > td_f_values_with_pid[-1]
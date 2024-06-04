"""Tests for met module solar position helpers."""

from do3se_met.solar_position import (
    calc_solar_elevation_list,
)


def test_calc_solar_elevation_list(snapshot):
    """Test the output of calc_solar_elevation_list."""
    result = calc_solar_elevation_list(lat=50, lon=-1.3, day_range=[0,5])
    snapshot.assert_match(result, 'solar_elevation_list_5_days')

import pytest
import pandas as pd
import numpy as np
import xarray as xr

from .process_outputs import dump_output_to_netcdf_grid


class TestDumpOutputToNetCDF:

    def test_runs_without_errors(self):
        X = 3
        Y = 4
        T = 5
        output_fields = ['a', 'b', 'c']
        full_output_data = {
            k: np.full((X, Y, T), None,
                    dtype=np.float64) for k in output_fields
        } if output_fields else None


        lat_data = np.full((X, Y), None, dtype=np.float16)
        lon_data = np.full((X, Y), None, dtype=np.float16)
        start_date = "01/01/2017"
        time_data = pd.date_range(start_date, periods=T)

        ds = dump_output_to_netcdf_grid(
            full_output_data,
            X, Y, T,
            lat_data=lat_data,
            lon_data=lon_data,
            time_data=time_data,
            output_fields=output_fields,
        )
        assert ds.a.shape == (X, Y, T)

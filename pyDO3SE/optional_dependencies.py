
try:
    import xarray
except:
    # xarray not installed create a dummy class for errors
    class xarray:
        # FIX: throwing an error here causes code linting to think xr will always throw error.
        def open_mfdataset():
            pass
            # raise ImportError("xarray must be installed to load NetCDF files")
        def open_dataset():
            pass
            # raise ImportError("xarray must be installed to load NetCDF files")

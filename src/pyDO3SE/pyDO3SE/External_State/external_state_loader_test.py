from inspect import isgenerator
from timeit import repeat
from dataclasses import asdict
import pytest
import math
import numpy as np
import pandas as pd
import textwrap
from io import StringIO
import itertools
import numpy.testing as nptest

from pyDO3SE.optional_dependencies import dask, xarray as xr
from .External_State_Shape import External_State_Shape, InputField
from .External_State_Config import FileTypes
from .external_state_loader import (
    extract_cell_data_from_netcdf,
    get_date_bounds_from_ext_data,
    load_external_state,
    process_csv_data,
    match_heading_to_field,
)

DEMO_FIELDS = [
    InputField('dd', 'Day of Year', int, 'index'),
    InputField('hr', 'Day of Year', int, 'index'),
    InputField('P', 'Day of Year', float, 'kPa'),
    InputField('temperature', 'Temperature in degrees C', float, 'C'),
    InputField('temperature', 'Temperature in degrees K', float, 'K'),
]

EXTERNAL_DATA_DEMO_PATH = 'examples/spanish_wheat/spanish_wheat_data_sample.csv'


class TestMatchHeadingToField:

    def test_match_heading_to_field(self):
        field, ignore_field = match_heading_to_field('dd', DEMO_FIELDS)
        assert field == DEMO_FIELDS[0].label
        assert ignore_field is None

    def test_match_heading_to_field_with_unit(self):
        field, ignore_field = match_heading_to_field('P, kPa', DEMO_FIELDS)
        assert field == DEMO_FIELDS[2].label
        assert ignore_field is None

    def test_match_heading_to_field_with_no_match(self):
        field, ignore_field = match_heading_to_field('x', DEMO_FIELDS)
        assert field == 'x'
        assert ignore_field == 'x'

    def test_match_heading_to_field_with_multi_match(self):
        field, ignore_field = match_heading_to_field('temperature, C', DEMO_FIELDS)
        assert field == DEMO_FIELDS[3].label
        assert ignore_field is None


class TestProcessCsvData:

    def test_process_csv_data(self):
        data = textwrap.dedent("""\
        dd,hr,P
        1,0,101.325
        1,1,101.125
        1,2,101.025
        """)

        processed_data = process_csv_data(StringIO(data), DEMO_FIELDS, True)
        nptest.assert_array_equal(processed_data['P'], [101.325, 101.125, 101.025])

    def test_process_csv_data_check_headers(self):
        data = textwrap.dedent("""\
        dd,hr,P
        1,0,101.325
        1,1,101.125
        1,2,101.025
        """)
        with pytest.raises(Exception):
            # attempting to process without headers throws exception
            process_csv_data(StringIO(data))

    def test_process_csv_data_extra_columns(self):
        data = textwrap.dedent("""\
        dd,hr,P,PAR,u
        1,0,101.325,20.4,2.3
        1,1,101.125,20.4,2.3
        1,2,101.025,20.4,2.3
        """)
        processed_data = process_csv_data(StringIO(data), DEMO_FIELDS, True)
        nptest.assert_array_equal(processed_data['P'], [101.325, 101.125, 101.025])
        assert 'dd' in processed_data
        assert 'hr' in processed_data
        assert 'P' in processed_data
        assert 'PAR' not in processed_data

    def test_process_csv_data_trim_whitespace(self):
        data = textwrap.dedent("""\
        dd, hr, P
        1,0,101.325
        1,1,101.125
        1,2,101.025
        """)

        processed_data = process_csv_data(StringIO(data), DEMO_FIELDS, True)
        assert 'dd' in processed_data
        assert 'hr' in processed_data
        assert 'P' in processed_data

    def test_process_csv_data_trim_units(self):
        data = textwrap.dedent("""\
        dd, hr, "P, kPa"
        1,0,101.325
        1,1,101.125
        1,2,101.025
        """)
        processed_data = process_csv_data(StringIO(data), DEMO_FIELDS, True)
        assert 'dd' in processed_data
        assert 'hr' in processed_data
        assert 'P' in processed_data

    def test_process_csv_data_no_fields(self):
        data = textwrap.dedent("""\
        dd, hr, "P, kPa"
        1,0,101.325
        1,1,101.125
        1,2,101.025
        """)
        processed_data = process_csv_data(
            StringIO(data), has_header_row=True, input_headers=['dd', 'hr', 'P'])
        assert 'dd' in processed_data
        assert 'hr' in processed_data
        assert 'P' in processed_data

    def test_process_csv_data_missing_headers(self):
        data = textwrap.dedent("""\
        dd, hr, "P, kPa"
        1,0,101.325
        1,1,101.125
        1,2,101.025
        """)
        with pytest.raises(Exception):
            process_csv_data(StringIO(data), has_header_row=False)


@pytest.mark.skip(reason="Old Test Failing")
class TestExternalStateLoader():

    def test_external_state_loader_invalid_type(self):
        with pytest.raises(Exception):
            next(load_external_state(EXTERNAL_DATA_DEMO_PATH, file_type='json'))

    def test_external_state_loader_missing_file(self):
        with pytest.raises(Exception):
            next(load_external_state('examples/spanish_wheat/data/MISSING.csv', file_type=FileTypes.CSV))

    class TestCSVLoader:
        def load_external_state(self):
            row_indexes = [0, 1, 2, 3, 4]
            external_state = next(load_external_state(
                external_state_file_location='examples/spanish_wheat/data/spanish_wheat_data_sample.csv',
                file_type=FileTypes.CSV,
                # ['dd', 'hr', 'P', 'Ts_C', 'precip', 'u', 'O3', 'VPD', 'PAR'],
                row_indexes=row_indexes,
            ))
            return external_state, row_indexes

        def test_external_state_loader(self):
            # TODO: Mock fs
            external_state, row_indexes = self.load_external_state()
            assert isinstance(external_state, External_State_Shape)
            assert external_state.O3[0] == 33.904
            assert external_state.VPD[0] == 0.05699
            assert external_state.PAR[0] == 0
            assert len(external_state.O3) == len(row_indexes)
            # TODO: Assert all loaded here

        def test_imports_correct_types(self):
            external_state, row_indexes = self.load_external_state()
            assert type(external_state.dd[0]) == int
            assert type(external_state.hr[0]) == int
            assert type(external_state.Ts_C[0]) == float

        def test_should_be_able_to_load_single_hour_of_data_csv(self):
            HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/hourly_runs/setup_01/hourly_data.csv'
            ext_data = next(load_external_state(
                HOURLY_EXTERNAL_DATA_DEMO_PATH,
                file_type=FileTypes.CSV))
            assert ext_data.dd == [71]

    class TestNetCDF:

        class Base:

            @pytest.fixture(autouse=True)
            def __setup(self):
                self.kwargs = {}
                self.variable_map = {
                    'time': 'XTIME',
                    '_SHAPE': 'td_2m',
                    'Ts_C': 'td_2m',
                    'P': 'pres',
                    'PAR': 'SWDOWN',  # Check is this PAR or PPFD
                    'precip': 'RAINNC',  # Should be sum of this and RAINC
                    'RH': 'rh',
                    'u': 'wspeed',
                    'O3': 'o3',
                    'Hd': 'HFX_FORCE',
                    # 'SWC': 'SMOISREL',
                    'snow_depth': 'SNOWH',
                }
                self.grid_coords = [
                    (0, 0),
                    (2, 0),
                    (0, 1),
                ]

            def test_has_required_inputs(self):
                # These should either be set in __setup above or _setup inside each sub test class
                assert self.data_location
                assert self.variable_map
                assert self.grid_coords

            def test_returns_a_generator(self):
                out = load_external_state(
                    self.data_location,
                    grid_coords=self.grid_coords,
                    file_type=FileTypes.NETCDF,
                    variable_map=self.variable_map,
                    zero_year=2017,
                    **self. kwargs,
                )
                assert isgenerator(out)

            def test_should_run_without_errors(self):
                next(load_external_state(
                    self.data_location,
                    grid_coords=self.grid_coords,
                    file_type=FileTypes.NETCDF,
                    variable_map=self.variable_map,
                    zero_year=2017,
                    **self.kwargs,
                ))

            def test_should_load_all_variables(self):
                out = next(load_external_state(
                    self.data_location,
                    grid_coords=self.grid_coords,
                    **dict(
                        file_type=FileTypes.NETCDF,
                        variable_map=self.variable_map,
                        zero_year=2017,
                    ),
                    **self.kwargs,
                ))
                for k in self.variable_map.keys():
                    if k in ["_SHAPE"]: continue
                    print(k)
                    assert getattr(out, k) is not None
                    assert len(getattr(out, k)) > 0
                    assert len(getattr(out, k)) == self.T
                    assert getattr(out, k)[0] is not None

            def test_should_load_additional_variables(self):
                out = next(load_external_state(
                    self.data_location,
                    grid_coords=self.grid_coords,
                    file_type=FileTypes.NETCDF,
                    variable_map=self.variable_map,
                    zero_year=2017,
                    **self.kwargs,
                ))
                additional_variables = ['dd', 'time', 'hr']
                for k in additional_variables:
                    assert getattr(out, k) is not None
                    assert len(getattr(out, k)) > 0
                    assert len(getattr(out, k)) == self.T
                    assert getattr(out, k)[0] is not None

                if self.T > 1:
                    assert out.hr[1] == out.hr[0] + 1
                    assert out.dd[24] == out.dd[0] + 1


        class TestSingleFileHour(Base):

            @pytest.fixture(autouse=True)
            def _setup(self):
                self.T = 1
                self.expected_start_day = 361
                self.expected_end_day = 361
                self.expected_start_date = "2017-12-27T00:00:00"
                self.expected_end_date = "2017-12-27T00:00:00"
                self.data_location = "examples/net_cdf/single_file_hour/inputs/demo_wrf_2017-12-27-00-00-00.nc"

        class TestSingleFileHourFollowingYear(Base):

            @pytest.fixture(autouse=True)
            def _setup(self):
                self.T = 1
                self.expected_start_day = 368
                self.expected_end_day = 368
                self.expected_start_date = "2018-01-03T23:00:00"
                self.expected_end_date = "2018-01-03T23:00:00"
                self.data_location = "examples/net_cdf/single_file_hour/inputs/demo_wrf_2018-01-03-23-00-00.nc"

        # @pytest.mark.skip(reason="Single file range not implemented")
        # class TestSingleFileRange(Base):

        #     @pytest.fixture(autouse=True)
        #     def _setup(self):
        #         self.data_location = "examples/net_cdf/single_file_range/inputs"

        # @pytest.mark.skip(reason="Multi file Single Hour Data not implemented")
        # class TestMultiFileHour(Base):

        #     @pytest.fixture(autouse=True)
        #     def _setup(self):
        #         self.data_location = "examples/net_cdf/multi_file_hour/inputs"
        #         self.kwargs = dict(
        #             multi_file_data=True,
        #             data_filter="",
        #         )
        #         self.T = 4 * 24

        class TestMultiFileRange(Base):

            @pytest.fixture(autouse=True)
            def _setup(self):
                self.data_location = "examples/net_cdf/multi_file_range/inputs"
                self.kwargs = dict(
                    multi_file_data=True,
                    data_filter='demo_wrf_2017-10',
                    file_suffix='',
                )
                self.expected_start_day = 361
                self.expected_end_day = 365
                self.expected_start_date = "2017-12-27T00:00:00"
                self.expected_end_date = "2017-12-30T23:00:00"
                self.T = 4 * 24

        class TestMultiFileRangeFollowingYear(Base):

            @pytest.fixture(autouse=True)
            def _setup(self):
                self.data_location = "examples/net_cdf/multi_file_range/inputs"
                self.kwargs = dict(
                    multi_file_data=True,
                    data_filter='demo_wrf_2017-11',
                    file_suffix='',
                )
                self.expected_start_day = 365
                self.expected_end_day = 368
                self.expected_start_date = "2017-12-31T00:00:00"
                self.expected_end_date = "2018-01-03T23:00:00"
                self.T = 4 * 24


@pytest.mark.skip(reason="Old Test Failing")
class TestExtractCellDataFromNetcdf:

    def _default_run(self):
        xi, yi, T = [0, 0, self.T]
        return extract_cell_data_from_netcdf(
            self.data_processed,
            xi, yi, T,
            time_key="time",
            variable_map=self.variable_map,
            index_counts=self.index_counts,
            use_dask=self.use_dask,
        )

    def _setup(self):
        DX = 2
        DY = 2
        DT = 3
        A = 15 + 8 * np.random.randn(DX, DY, DT)
        B = 10 * np.random.rand(DX, DY, DT)
        dd = 10 * np.arange(DT)
        hr = 10 * np.arange(DT)
        lon = [[-99.83, -99.32], [-99.79, -99.23]]
        lat = [[42.25, 42.21], [42.63, 42.59]]
        time = pd.date_range('2011-01-01', periods=DT)

        if self.use_dask:
            A = dask.array.from_array(A)
            B = dask.array.from_array(B)
            self.dd = dask.array.from_array(dd)
            hr = dask.array.from_array(hr)

        self.A = A
        self.dd = dd
        self.data_processed = xr.Dataset(
            data_vars=dict(
                A=(["x", "y", "time"], self.A),
                B=(["x", "y", "time"], B),
                dd=(['time'], dd),
                hr=(['time'], hr),
            ),
            coords=dict(
                lon=(["x", "y"], lon),
                lat=(["x", "y"], lat),
                time=time,
            ),
            attrs=dict(description="Weather related data."),
        )

        self.T = 1
        self.variable_map = {
            'a': 'A',
            'b': 'B',
        }
        self.index_counts = [3, 3]
        self.chunks = None

    class Base:

        def test_should_work_without_errors(self):
            self._setup()
            self._default_run()

        def test_should_be_able_to_access(self):
            self._setup()
            out = self._default_run()
            v = out.get('a')[0]
            assert v == self.A[0, 0, 0]

            v = out.get('dd')[0]
            assert v == self.dd[0]

    class TestWithoutDask(Base):
        def _setup(self):
            self.use_dask = False
            return TestExtractCellDataFromNetcdf._setup(self)

        def _default_run(self):
            return TestExtractCellDataFromNetcdf._default_run(self)

    class TestWithDask(Base):
        def _setup(self):
            self.use_dask = True
            return TestExtractCellDataFromNetcdf._setup(self)

        def _default_run(self):
            return TestExtractCellDataFromNetcdf._default_run(self)


@pytest.mark.skip(reason="Multi file data not implemented")
class TestLoadHourlyDataNetCDF:
    def test_should_be_able_to_load_single_hour_of_data_netcdf_emep(self):
        HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/emep/inputs'

        variable_map = dict([
            ('time', 'time'),
            ('O3', 'O3_45m'),
            ('u', 'u_45m'),
            ('Ts_C', 't2m'),  # Needs transform
            ('Hd', 'SH_Wm2'),
            ('precip', 'Precip'),
            ('cloudfrac', 'CloudFrac'),
            ('P', 'Psurf'),
            ('RH', 'rh2m'),
            # ('ustar', 'ustar'),
        ])
        coords = [
            (0, 0),
        ]
        coord_i, ext_data = next(load_external_state(
            HOURLY_EXTERNAL_DATA_DEMO_PATH,
            FileTypes.NETCDF,
            coords=coords,
            variable_map=variable_map,
            multi_file_data=True,
        ))

        assert coord_i == coords[0]
        for k in variable_map.keys():
            if k != 'time':
                try:
                    assert getattr(ext_data, k) is not None
                    assert len(getattr(ext_data, k)) > 0
                except AssertionError:
                    raise AssertionError(f'Key: {k} missing in data')
        assert ext_data.Ts_C[0] > 100

        assert ext_data.hr is not None
        assert ext_data.hr[0] == 0
        assert ext_data.hr[1] == 1
        assert ext_data.hr[23] == 23
        assert ext_data.hr[24] == 0

        assert ext_data.dd is not None
        assert ext_data.dd[0] == 1
        assert ext_data.dd[23] == 1
        assert ext_data.dd[24] == 2
        assert ext_data.dd[-1] == 365

    def test_should_be_able_to_load_single_hour_of_data_netcdf_preprocess(self):
        HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/emep/inputs'

        variable_map = dict([
            ('time', 'time'),
            ('O3', 'O3_45m'),
            ('u', 'u_45m'),
            ('Ts_C', 't2m'),  # Needs transform
            ('Hd', 'SH_Wm2'),
            ('precip', 'Precip'),
            ('cloudfrac', 'CloudFrac'),
            ('P', 'Psurf'),
            ('RH', 'rh2m'),
            # ('ustar', 'ustar'),
        ])

        preprocess_map = dict([
            # ('t2m', lambda row: row.t2m - 273.15), # Needs transform
            ('Ts_C', lambda row: row - 273.15),  # Needs transform
        ])
        coords = [(0, 0)]
        (xi, yi), ext_data = next(load_external_state(
            HOURLY_EXTERNAL_DATA_DEMO_PATH,
            FileTypes.NETCDF,
            coords=coords,
            variable_map=variable_map,
            multi_file_data=True,
            preprocess_map=preprocess_map,
        ))

        assert ext_data.Ts_C[0] < 100

    def test_should_be_able_to_load_single_hour_of_data_netcdf_wrfchem(self):
        """This dataset is a single file containing an hours data for all variables."""
        HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'

        variable_map = {
            'time': 'XTIME',
            'Ts_C': 'td_2m',
            'P': 'pres',
            'PAR': 'SWDOWN',  # Check is this PAR or PPFD
            'precip': 'RAINNC',  # Should be sum of this and RAINC
            'RH': 'rh',
            'u': 'wspeed',
            'O3': 'o3',
            'Hd': 'HFX_FORCE',
            # 'SWC': 'SMOISREL',
            'snow_depth': 'SNOWH',
        }
        coords = [
            (0, 0),
        ]
        coord, ext_data = next(load_external_state(
            HOURLY_EXTERNAL_DATA_DEMO_PATH,
            FileTypes.NETCDF,
            coords=coords,
            variable_map=variable_map,
            multi_file_data=False,
        ))
        assert coord == coords[0]
        for k in variable_map.keys():
            if k != 'time':
                try:
                    assert getattr(ext_data, k) is not None
                    assert len(getattr(ext_data, k)) > 0
                except AssertionError:
                    raise AssertionError(f'Key: {k} missing in data')

        assert ext_data.hr is not None
        assert ext_data.hr[0] == 1

        assert ext_data.dd is not None
        assert ext_data.dd[0] == 274


@pytest.mark.skip(reason="Moving to tests above")
class TestLoadHourlyDataNetCDFGrid:
    """Test that we can load multiple grid cells at once"""
    # def test_should_be_able_to_load_single_hour_of_data_netcdf_emep(self):
    #     HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/emep/inputs'

    #     variable_map = dict([
    #         ('time', 'time'),
    #         ('O3', 'O3_45m'),
    #         ('u', 'u_45m'),
    #         ('Ts_C', 't2m'), # Needs transform
    #         ('Hd', 'SH_Wm2'),
    #         ('precip', 'Precip'),
    #         ('cloudfrac', 'CloudFrac'),
    #         ('P', 'Psurf'),
    #         ('RH', 'rh2m'),
    #         # ('ustar', 'ustar'),
    #     ])
    #     ext_data = load_external_state(
    #         HOURLY_EXTERNAL_DATA_DEMO_PATH,
    #         FileTypes.NETCDF,
    #         x=0,
    #         y=0,
    #         variable_map=variable_map,
    #         multi_file_data=True,
    #     )()
    #     for k in variable_map.keys():
    #         try:
    #             assert getattr(ext_data, k) is not None
    #             assert len(getattr(ext_data, k)) > 0
    #         except AssertionError:
    #             raise AssertionError(f'Key: {k} missing in data')
    #     assert ext_data.Ts_C[0] > 100

    #     assert ext_data.hr is not None
    #     assert ext_data.hr[0] == 0
    #     assert ext_data.hr[1] == 1
    #     assert ext_data.hr[23] == 23
    #     assert ext_data.hr[24] == 0

    #     assert ext_data.dd is not None
    #     assert ext_data.dd[0] == 1
    #     assert ext_data.dd[23] == 1
    #     assert ext_data.dd[24] == 2
    #     assert ext_data.dd[-1] == 365

    # def test_should_be_able_to_load_single_hour_of_data_netcdf_preprocess(self):
    #     HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/emep/inputs'

    #     variable_map = dict([
    #         ('time', 'time'),
    #         ('O3', 'O3_45m'),
    #         ('u', 'u_45m'),
    #         ('Ts_C', 't2m'), # Needs transform
    #         ('Hd', 'SH_Wm2'),
    #         ('precip', 'Precip'),
    #         ('cloudfrac', 'CloudFrac'),
    #         ('P', 'Psurf'),
    #         ('RH', 'rh2m'),
    #         # ('ustar', 'ustar'),
    #     ])

    #     preprocess_map = dict([
    #         # ('t2m', lambda row: row.t2m - 273.15), # Needs transform
    #         ('Ts_C', lambda row: row - 273.15), # Needs transform
    #     ])

    #     ext_data = load_external_state(
    #         HOURLY_EXTERNAL_DATA_DEMO_PATH,
    #         FileTypes.NETCDF,
    #         x=0,
    #         y=0,
    #         variable_map=variable_map,
    #         multi_file_data=True,
    #         preprocess_map=preprocess_map,
    #     )()

    #     assert ext_data.Ts_C[0] < 100

    def test_should_be_able_to_load_single_hour_of_data_netcdf_wrfchem(self, snapshot):
        """This dataset is a single file containing an hours data for all variables for all grid cells."""
        HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'

        variable_map = {
            'time': 'XTIME',
            'Ts_C': 'td_2m',
            'P': 'pres',
            'PAR': 'SWDOWN',  # Check is this PAR or PPFD
            'precip': 'RAINNC',  # Should be sum of this and RAINC
            'RH': 'rh',
            'u': 'wspeed',
            'O3': 'o3',
            'Hd': 'HFX_FORCE',
            # 'SWC': 'SMOISREL',
            'snow_depth': 'SNOWH',
        }
        x = [0, 1, 2]
        y = [0, 1, 2]
        coords = list(itertools.product(x, y))
        ext_data_list = load_external_state(
            HOURLY_EXTERNAL_DATA_DEMO_PATH,
            FileTypes.NETCDF,
            coords=coords,
            variable_map=variable_map,
            multi_file_data=False,
        )

        for (xi, yi), ext_data in ext_data_list:
            for k in variable_map.keys():
                if k != 'time':
                    try:
                        assert getattr(ext_data, k) is not None
                        assert len(getattr(ext_data, k)) > 0
                    except AssertionError:
                        raise AssertionError(f'Key: {k} missing in data')
                    except AttributeError:
                        raise AttributeError(f'Key: {k} not in external state shape')

            assert ext_data.hr is not None
            assert ext_data.hr[0] == 0

            assert ext_data.dd is not None
            assert ext_data.dd[0] == 274

            snapshot.assert_match(asdict(ext_data), f'{xi}:{yi}')

    def test_should_create_external_state_iterator_covering_all_input_coords(self):
        HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'

        variable_map = {
            'time': 'XTIME',
            'Ts_C': 'td_2m',
            'P': 'pres',
            'PAR': 'SWDOWN',  # Check is this PAR or PPFD
            'precip': 'RAINNC',  # Should be sum of this and RAINC
            'RH': 'rh',
            'u': 'wspeed',
            'O3': 'o3',
            'Hd': 'HFX_FORCE',
            # 'SWC': 'SMOISREL',
            'snow_depth': 'SNOWH',
        }
        x = [0, 1, 2]
        y = [0, 1, 2]
        coords = list(itertools.product(x, y))
        e_state_iterator = load_external_state(
            HOURLY_EXTERNAL_DATA_DEMO_PATH,
            FileTypes.NETCDF,
            coords=coords,
            variable_map=variable_map,
            multi_file_data=False,
        )
        e_state_output = list(e_state_iterator)
        assert len(e_state_output) == len(coords) == len(x) * len(y)

    def test_should_load_data_quickly(self):
        HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'
        # HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/wrfchem/inputs/wrfchem_demo_out_01_10_01.nc'

        variable_map = {
            'time': 'XTIME',
            'Ts_C': 'td_2m',
            'P': 'pres',
            'PAR': 'SWDOWN',  # Check is this PAR or PPFD
            'precip': 'RAINNC',  # Should be sum of this and RAINC
            'RH': 'rh',
            'u': 'wspeed',
            'O3': 'o3',
            'Hd': 'HFX_FORCE',
            # 'SWC': 'SMOISREL',
            'snow_depth': 'SNOWH',
        }
        x = [0, 1, 2]
        y = [0, 1, 2]
        coords = list(itertools.product(x, y))

        def run():

            list(load_external_state(
                HOURLY_EXTERNAL_DATA_DEMO_PATH,
                FileTypes.NETCDF,
                coords=coords,
                variable_map=variable_map,
                multi_file_data=False,
            ))
        t = min(repeat(run, number=2, repeat=3))
        assert t < 3.5

    def test_should_set_date_as_offset_from_base_date(self):
        HOURLY_EXTERNAL_DATA_DEMO_PATH = 'examples/net_cdf/wrfchem/inputs/cropped_wrfchem_demo_out_01_10_00.nc'

        variable_map = {
            'time': 'XTIME',
            'Ts_C': 'td_2m',
        }
        x = [0]
        y = [0]
        coords = list(itertools.product(x, y))
        e_state_iterator = load_external_state(
            HOURLY_EXTERNAL_DATA_DEMO_PATH,
            FileTypes.NETCDF,
            coords=coords,
            variable_map=variable_map,
            multi_file_data=False,
            zero_year=2016,
        )
        coord, e_state_output = next(e_state_iterator)
        assert e_state_output.dd == [365 + 274]
        assert e_state_output.date == ['2017-10-01T00:00:00']


def test_should_be_able_to_get_time_bounds():
    expected_start_day = 10
    expected_start_date = pd.Timestamp("2017-10-01T00:00:00")
    expected_end_day = expected_start_day + 4
    expected_end_date = pd.Timestamp("2017-10-04T23:00:00")
    T = 4 * 24
    external_state = External_State_Shape(
        dd=[dd for dd in range(expected_start_day, expected_end_day) for _ in range(24)],
        hr=[hr for _ in range(expected_start_day, expected_end_day) for hr in range(24)],
        time=pd.date_range(expected_start_date, periods=T, freq='1H'),
    )
    out = get_date_bounds_from_ext_data(
        external_state,
    )
    assert out is not None
    assert out.start_day == expected_start_day
    assert out.start_date == expected_start_date
    assert out.end_day == expected_start_day + math.ceil(T / 24) - 1
    assert out.end_date == expected_end_date
    assert out.row_count == T
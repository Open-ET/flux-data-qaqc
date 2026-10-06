# -*- coding: utf-8 -*-

import configparser
import shutil
import pytest
from pathlib import Path

import numpy as np
import pandas as pd
import xarray
# refet >= 0.5 renamed its radiation functions (dropped leading underscore)
try:
    from refet.calcs import ra_daily as _ra_daily, ra_hourly as _ra_hourly
except ImportError:
    from refet.calcs import _ra_daily, _ra_hourly

from fluxdataqaqc import Data
from fluxdataqaqc import QaQc
from fluxdataqaqc import Plot
from fluxdataqaqc import Convert
from fluxdataqaqc import util


@pytest.fixture(scope="session")
def data():
    """
    Paths used by tests, tests read the example data in place.

    Anything a test writes (plots, CSVs) goes to a pytest temporary
    directory so the examples folder is never modified.
    """
    d = {}
    # repository root is the parent of the tests directory
    package_root_dir = Path(__file__).resolve().parent.parent
    d['package_root_dir'] = package_root_dir

    return d


def _tmp_config(config, tmp_path, set_opts=None, remove_opts=None):
    """
    Copy an example config to a temporary directory so tests can change it.

    The file is copied as text so comments are kept, and the climate file
    path is made absolute so the copy still reads the example data in
    place. ``set_opts`` is a list of (section, option, value) tuples and
    ``remove_opts`` a list of (section, option) tuples, options to remove
    are matched by name in any section.
    """
    config = Path(config)
    out_file = Path(tmp_path) / config.name
    shutil.copy(config, out_file)

    cp = configparser.ConfigParser(interpolation=None)
    cp.read(config)
    climate_file = config.parent / cp.get('METADATA', 'climate_file_path')
    util.set_config_option(
        out_file, 'METADATA', 'climate_file_path', str(climate_file.resolve())
    )
    for section, option, value in set_opts or []:
        util.set_config_option(out_file, section, option, value)

    # remove options by dropping their lines
    if remove_opts:
        names = [option.lower() for section, option in remove_opts]
        lines = out_file.read_text().splitlines(keepends=True)
        lines = [
            line for line in lines
            if line.split('=')[0].strip().lower() not in names
        ]
        out_file.write_text(''.join(lines))

    return out_file


def _fake_gridMET_server(fail=False):
    """
    Stand-in for :func:`xarray.open_dataset` on the gridMET THREDDS server.

    Returns small gridMET-like datasets with constant values so gridMET
    download and gap filling can be tested without network access.
    """
    # constant daily values for each gridMET variable
    values = {
        'daily_mean_reference_evapotranspiration_alfalfa': 6.0,
        'daily_mean_reference_evapotranspiration_grass': 5.0,
        'precipitation_amount': 1.0,
    }

    def open_dataset(url, **kwargs):
        if fail:
            raise OSError('NetCDF: I/O failure (test)')
        meta = [
            m for m in QaQc.gridMET_meta.values() if m['nc_suffix'] in url
        ][0]
        days = pd.date_range('2008-01-01', '2013-12-31', name='day')
        lats = [36.40, 36.44]
        lons = [-99.44, -99.40]
        data = np.full(
            (len(days), len(lats), len(lons)), values[meta['name']]
        )
        return xarray.Dataset(
            {meta['name']: (('day', 'lat', 'lon'), data), 'crs': 0},
            coords={'day': days, 'lat': lats, 'lon': lons}
        )

    return open_dataset


class TestData(object):

    @pytest.fixture(autouse=True)
    def setup_method(self, data):
        config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'US-Tw3_config.ini'
        # for testing auto QC flag filtering
        self.fluxnet_config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'fluxnet_config.ini'
        self.data_obj = Data(config)
        self.header_Tw3 = [
            'TIMESTAMP_START', 'TIMESTAMP_END', 'CO2', 'H2O', 'CH4', 'FC',
            'FCH4', 'FC_SSITC_TEST', 'FCH4_SSITC_TEST', 'G', 'H', 'LE',
            'H_SSITC_TEST', 'LE_SSITC_TEST', 'WD', 'WS', 'USTAR', 'ZL', 'TAU',
            'MO_LENGTH', 'V_SIGMA', 'W_SIGMA', 'TAU_SSITC_TEST', 'PA', 'RH',
            'TA', 'VPD_PI', 'T_SONIC', 'T_SONIC_SIGMA', 'SWC_1_1_1',
            'SWC_1_2_1', 'TS_1_1_1', 'TS_1_2_1', 'TS_1_3_1', 'TS_1_4_1',
            'TS_1_5_1', 'NETRAD', 'PPFD_DIF', 'PPFD_IN', 'PPFD_OUT', 'SW_IN',
            'SW_OUT', 'LW_IN', 'LW_OUT', 'P', 'FC_PI_F', 'RECO_PI_F',
            'GPP_PI_F', 'H_PI_F', 'LE_PI_F'
        ]


    def test_Data_init(self):
        # checking constructor of Data class
        assert self.data_obj.climate_file.name == 'AMF_US-Tw3_BASE_HH_5-5.csv'
        assert self.data_obj.config_file.stem == 'US-Tw3_config'
        assert self.data_obj.elevation == -9
        assert self.data_obj.latitude == 38.1159
        # dataframe should not be loaded
        assert self.data_obj._df is None
        # header of climate file should have been read however
        assert (self.data_obj.header == self.header_Tw3).all()
        assert self.data_obj.units.get('wd') == 'azimuth (degrees)'


    def test_config_parser(self):
        assert self.data_obj.config.sections() == ['METADATA', 'DATA']
        assert self.data_obj.config.get('METADATA', 'igbp') == 'CRO'

    def test_units_and_Convert_inheritance(self):
        units = self.data_obj.units
        assert isinstance(units, dict)
        assert units.get('LE') == 'w/m2'
        assert self.data_obj.pretty_unit_names.get('kpa') == 'kPa'

    def test_Data_load_df(self):
        # check on loading climate timeseries input data into Pandas DataFrame
        df = self.data_obj.df
        assert isinstance(df, pd.DataFrame)
        df = df.rename(columns=self.data_obj.inv_map)
        core_cols = ['LE', 'H', 'Rn', 'G']
        assert set(core_cols).issubset(df.columns)
        assert len(df[core_cols].dropna()) == 73868


    def test_datetime_index(self):
        df = self.data_obj.df
        assert isinstance(df.index[0], pd.Timestamp)
        assert df.index[0].year == 2013

    def test_auto_average(self):
        # check that automatic averaging of soil moisture was done
        df = self.data_obj.df
        assert 'theta_mean' in df.columns
        assert 'theta_mean' in self.data_obj.variables
        assert np.isclose(df.theta_mean.mean(), 23.154676446960238)

    def test_auto_calcs(self):
        # check that es, vp, and t_dew were calculated
        self.data_obj.df
        assert {'es', 'vp', 't_dew'}.issubset(self.data_obj.df.columns)
        assert np.floor(self.data_obj.df.t_dew.dropna().iloc[0]) == 11

    def test_ACSE_refET(self):
        ts = self.data_obj.hourly_ASCE_refET()
        assert len(ts) == 47544
        # same value for refet 0.3.10, 0.4 and 0.5 (keyword arguments)
        assert np.isclose(ts.mean(), 0.1830681034395387)

    def test_ASCE_refET_anemometer_height_from_config(self):
        # config values are strings, this used to raise a TypeError
        self.data_obj.config.set('METADATA', 'anemometer_height', '2')
        ts = self.data_obj.hourly_ASCE_refET()
        assert np.isclose(ts.mean(), 0.1830681034395387)

    def test_ASCE_refET_invalid_reference(self):
        with pytest.raises(ValueError):
            self.data_obj.hourly_ASCE_refET(reference='medium')

    def test_Data_plots(self, tmp_path):
        assert self.data_obj.plot_file is None
        # default file name but saved outside of the examples folder
        self.data_obj.out_dir = tmp_path
        self.data_obj.plot()
        assert self.data_obj.plot_file.name ==\
            f'{self.data_obj.site_id}_input_plots.html'
        assert self.data_obj.plot_file.is_file()

    #def test_xl_reader(self):
    #    d = Data(self.fluxnet_config)
    #    d.xl_parser = 'openpyxl'
    #    df = d.df
    #    assert isinstance(df, pd.DataFrame)
    #    d = Data(self.fluxnet_config)
    #    d.xl_parser = 'xlrd'
    #    df = d.df
    #    assert isinstance(df, pd.DataFrame)

        
    def test_input_qc_flag_filtering(self):
        d = Data(self.fluxnet_config)
        LE_init = d.df.rename(columns=d.inv_map).LE
        d.apply_qc_flags(threshold=1)
        LE_qc = d.df.rename(columns=d.inv_map).LE
        assert LE_qc.mean() < LE_init.mean()

    def test_weighted_average_from_config(self, data):
        config = data['package_root_dir']/'examples'\
            /'Config_options'/'config_for_multiple_soil_vars.ini'

        d = Data(config)
        assert d.soil_var_weight_pairs.get('g_1').get('weight') == '1'
        assert d.soil_var_weight_pairs.get('g_2').get('weight') == '10'
        assert not 'g_mean' in d.header
        df = d.df.rename(columns=d.inv_map)
        assert 'g_mean' in d.df.columns
        assert np.isclose(
            d.soil_var_weight_pairs.get('g_1').get('weight'), 
            0.045454545454545456
        )
        g_mean = df[['g_1', 'g_2', 'g_3','g_4']].mean(axis=1)
        g_weighted_mean = df.G
        assert (g_mean != g_weighted_mean).any()

    def test_missing_units_raise(self, data, tmp_path):
        # units are never assumed for variables used in calculations
        config = data['package_root_dir']/'examples'\
            /'Config_options'/'config_for_QC_flag_filtering.ini'
        config = _tmp_config(
            config, tmp_path,
            remove_opts=[('DATA', 'latent_heat_flux_units')]
        )
        with pytest.raises(ValueError, match='latent_heat_flux_units'):
            Data(config)

    def test_soil_heat_flux_column_named_G(self, data, tmp_path):
        # g_1 uses the same column as G, its units are in the config but
        # not stored in Data.units, this was rejected in version 0.3.1
        config = data['package_root_dir']/'examples'\
            /'Config_options'/'config_for_QC_flag_filtering.ini'
        config = _tmp_config(
            config, tmp_path,
            set_opts=[('DATA', 'g_1', 'G'), ('DATA', 'g_1_units', 'w/m2')]
        )
        q = QaQc(Data(config))
        assert 'G' in q.df.rename(columns=q.inv_map).columns

    def test_missing_units_for_plot_only_variable(self, data, tmp_path):
        # potential shortwave units are only used for plot labels
        config = data['package_root_dir']/'examples'\
            /'Config_options'/'config_for_QC_flag_filtering.ini'
        config = _tmp_config(
            config, tmp_path, remove_opts=[('DATA', 'shortwave_pot_units')]
        )
        q = QaQc(Data(config))
        assert q.units.get('sw_pot') is None
        assert q.units.get('LE') == 'w/m2'

    def test_excel_header_and_data(self):
        d = Data(self.fluxnet_config)
        raw = pd.read_excel(d.climate_file)
        # header is the first row of the first sheet
        assert list(d.header) == list(raw.columns)
        df = d.df
        raw = raw.replace(d.na_val, np.nan)
        for col in ['LE_F_MDS', 'H_F_MDS', 'NETRAD', 'G_F_MDS']:
            assert np.isclose(df[col].mean(), raw[col].mean())

    def test_soil_var_qc_flag_names(self, data, tmp_path):
        # QC name for one of multiple soil heat flux variables
        config = data['package_root_dir']/'examples'\
            /'Config_options'/'config_for_multiple_soil_vars.ini'
        config = _tmp_config(
            config, tmp_path, set_opts=[('DATA', 'g_1_qc', 'G_2_1_1')]
        )
        d = Data(config)
        # QC column is not treated as another soil heat flux variable
        assert 'g_1_qc' not in d.variables
        assert 'g_1_qc' not in d.soil_var_weight_pairs
        assert d.variables['g_1_qc_flag'] == 'G_2_1_1'
        assert d.qc_var_pairs['G_1_1_1'] == 'G_2_1_1'

    def test_calc_pes(self):
        df = self.data_obj.df.rename(columns=self.data_obj.inv_map)
        assert 'gpp' in df.columns

        self.data_obj.calc_pes()
        df = self.data_obj.df

        assert {'pes', 'pes_flux'}.issubset(df.columns)
        assert self.data_obj.units['pes'] == 'j/m2'
        assert self.data_obj.units['pes_flux'] == 'w/m2'

        _, _, dt_seconds = util.get_subdaily_timestep_info(df)

        gpp = df.rename(columns=self.data_obj.inv_map)['gpp'].clip(lower=0)
        mask = gpp.notna() & df['pes_flux'].notna() & df['pes'].notna()

        assert np.allclose(
            df.loc[mask, 'pes_flux'],
            gpp.loc[mask] * 0.422,
            rtol=0,
            atol=1e-12
        )
        assert np.allclose(
            df.loc[mask, 'pes'],
            df.loc[mask, 'pes_flux'] * dt_seconds,
            rtol=0,
            atol=1e-9
        )



class TestQaQc(object):

    def _test_df(self, gap_length, sw_pot, rn=100.,
            start='2020-06-01'):
        index = pd.date_range(
            start,
            periods=gap_length + 2,
            freq='30min'
        )
        vals = [0.] + [np.nan] * gap_length + [gap_length + 1.]

        df = pd.DataFrame(
            {
                'LE': vals,
                'Rn': [rn] * len(index),
            },
            index=index
        )
        sw_pot = pd.Series(sw_pot, index=index)

        return df, sw_pot

    def test_day_gap_equal_to_limit(self):
        df, sw_pot = self._test_df(4, 500.)
        interped = QaQc._interpolate_short_gaps(
            df, 4, 8, sw_pot
        )

        assert interped.LE.notna().all()

    def test_day_gap_larger_than_limit(self):
        df, sw_pot = self._test_df(5, 500.)
        interped = QaQc._interpolate_short_gaps(
            df, 4, 8, sw_pot
        )

        assert interped.LE.iloc[1:-1].isna().all()

    def test_night_gap_equal_to_limit(self):
        df, sw_pot = self._test_df(
            8, 0., rn=-100.,
            start='2020-06-01 22:00'
        )
        interped = QaQc._interpolate_short_gaps(
            df, 4, 8, sw_pot
        )

        assert interped.LE.notna().all()

    def test_negative_rn_during_day_uses_day_limit(self):
        df, sw_pot = self._test_df(
            5, 500., rn=-100.
        )
        interped = QaQc._interpolate_short_gaps(
            df, 4, 8, sw_pot
        )

        assert interped.LE.iloc[1:-1].isna().all()

    def test_positive_rn_at_solar_night_uses_day_limit(self):
        df, sw_pot = self._test_df(
            5, 0., rn=25.
        )
        interped = QaQc._interpolate_short_gaps(
            df, 4, 8, sw_pot
        )

        assert interped.LE.iloc[1:-1].isna().all()

    def test_sunrise_gap_uses_both_endpoints(self):
        index = pd.date_range(
            '2020-06-01 07:30',
            periods=3,
            freq='30min'
        )
        df = pd.DataFrame(
            {
                'LE': [35.23, np.nan, 165.80],
                'Rn': [-44.36, 107.82, 313.47],
            },
            index=index
        )
        sw_pot = pd.Series(
            [0., 50., 200.],
            index=index
        )

        interped = QaQc._interpolate_short_gaps(
            df, 4, 8, sw_pot
        )

        assert np.isclose(
            interped.LE.iloc[1],
            (35.23 + 165.80) / 2
        )
        
    def _arm_example_df(self, data):
        """Load US-ARM site data with real gap issues"""
        config = (
            data['package_root_dir'] /
            'examples' /
            'Config_options' /
            'config_for_multiple_soil_vars.ini'
        )

        data_obj = Data(config)

        return (
            data_obj.df
            .rename(columns=data_obj.inv_map)
            .sort_index()
        )

    def test_real_sunrise_gap_uses_both_endpoints(self, data):
        df = self._arm_example_df(data)

        test_df = df.loc[
            '2003-08-31 07:00':'2003-08-31 08:00',
            ['LE', 'Rn']
        ].copy()

        assert test_df.loc[
            '2003-08-31 07:30',
            'LE'
        ] != test_df.loc[
            '2003-08-31 07:30',
            'LE'
        ]

        assert test_df.loc[
            '2003-08-31 07:00',
            'Rn'
        ] < 0

        assert test_df.loc[
            '2003-08-31 07:30',
            'Rn'
        ] > 0

        interped = QaQc._interpolate_short_gaps(
            test_df,
            max_gap=4,
            max_night_gap=8
        )

        expected = (
            test_df.loc[
                '2003-08-31 07:00',
                'LE'
            ] +
            test_df.loc[
                '2003-08-31 08:00',
                'LE'
            ]
        ) / 2

        assert np.isclose(
            interped.loc[
                '2003-08-31 07:30',
                'LE'
            ],
            expected
        )

    def test_real_oversized_gap_is_not_partially_filled(self, data):
        df = self._arm_example_df(data)

        test_df = df.loc[
            '2003-02-09 10:30':'2003-02-09 13:30',
            ['LE', 'Rn']
        ].copy()

        gap = test_df.loc[
            '2003-02-09 11:00':'2003-02-09 13:00',
            'LE'
        ]

        assert len(gap) == 5
        assert gap.isna().all()

        interped = QaQc._interpolate_short_gaps(
            test_df,
            max_gap=4,
            max_night_gap=8
        )

        result_gap = interped.loc[
            '2003-02-09 11:00':'2003-02-09 13:00',
            'LE'
        ]

        assert result_gap.isna().all()
        
    def test_subdaily_sw_pot_matches_refet_hourly(self):
        latitude = 40.0
        longitude = -75.0
        utc_offset = -5.0

        index = pd.date_range(
            '2020-06-21 08:00',
            periods=9,
            freq='h'
        )

        sw_pot = QaQc._calc_subdaily_sw_pot(
            index,
            latitude,
            longitude,
            utc_offset,
            1.0
        )

        midpoint = index + pd.Timedelta(minutes=30)
        utc_time_mid = (
            midpoint.hour.to_numpy()
            + midpoint.minute.to_numpy() / 60
            - utc_offset
        ) % 24

        expected = _ra_hourly(
            np.full(len(index), np.deg2rad(latitude)),
            np.full(len(index), np.deg2rad(longitude)),
            midpoint.dayofyear.to_numpy(),
            utc_time_mid,
            method='asce'
        )

        assert np.allclose(
            sw_pot.to_numpy(),
            expected,
            rtol=0,
            atol=1e-10
        )

    def test_subdaily_sw_pot_sums_to_daily_ra(self):
        latitude = 40.0
        longitude = -75.0
        utc_offset = -5.0

        index = pd.date_range(
            '2020-06-21 00:00',
            periods=48,
            freq='30min'
        )

        sw_pot = QaQc._calc_subdaily_sw_pot(
            index,
            latitude,
            longitude,
            utc_offset,
            0.5
        )

        expected = _ra_daily(
            np.array([np.deg2rad(latitude)]),
            np.array([index[0].dayofyear]),
            method='asce'
        )[0]

        assert np.isclose(
            sw_pot.sum(),
            expected,
            rtol=0,
            atol=1e-8
        )

    def test_subdaily_sw_pot_day_and_night(self):
        latitude = 40.0
        longitude = -75.0
        utc_offset = -5.0

        index = pd.date_range(
            '2020-06-21 00:00',
            periods=48,
            freq='30min'
        )

        sw_pot = QaQc._calc_subdaily_sw_pot(
            index,
            latitude,
            longitude,
            utc_offset,
            0.5
        )

        assert (sw_pot >= 0).all()
        assert sw_pot.loc['2020-06-21 00:00'] == 0
        assert sw_pot.loc['2020-06-21 12:00'] > 0
        assert sw_pot.max() > 0

    def _daily_cf(self, values, start='2010-03-01'):
        """Daily correction factor series for energy balance tests"""
        index = pd.date_range(start, periods=len(values), freq='D')
        return pd.Series(values, index=index, dtype=float)

    def test_ebc_cf_method_1_at_record_start(self):
        # first days used to get no moving window value
        cf = self._daily_cf([1.2] * 30)
        cf.iloc[0] = 1.0
        cf_corr, method = QaQc._ebc_cf_moving_windows(cf)
        assert (method == 1).all()
        # day 0 window is days 0-7, median of 1.0 and seven 1.2 values
        assert np.isclose(cf_corr.iloc[0], 1.2)

    def test_ebc_cf_method_2_is_mean(self):
        # only 3 values within +/- 7 days, method 2 uses the mean
        cf = self._daily_cf([np.nan] * 30)
        cf.iloc[[10, 12, 13]] = [1.0, 1.1, 1.6]
        cf_corr, method = QaQc._ebc_cf_moving_windows(cf)
        assert method.iloc[11] == 2
        assert np.isclose(cf_corr.iloc[11], np.mean([1.0, 1.1, 1.6]))
        # no values within +/- 5 days
        assert np.isnan(cf_corr.iloc[25]) and np.isnan(method.iloc[25])

    def test_ebc_cf_method_3_previous_and_next_years(self):
        cf = self._daily_cf([1.0] * 365 + [1.5] * 365 + [1.4] * 365,
            start='2010-01-01')
        adjacent = QaQc._ebc_cf_adjacent_years(cf)
        # middle year, mean of previous (1.0) and next (1.4) years
        assert np.isclose(adjacent.loc['2011-06-15'], 1.2)
        # first year, only the next year exists
        assert np.isclose(adjacent.loc['2010-06-15'], 1.5)
        # values in a gap of the middle year come from the other years
        cf.loc['2011-06-01':'2011-07-31'] = np.nan
        adjacent = QaQc._ebc_cf_adjacent_years(cf)
        assert np.isclose(adjacent.loc['2011-06-15'], 1.2)

    def test_ebc_cf_climatology_lines_up_with_day_of_year(self):
        # 200 days from March 1, climatology used to shift by record start
        cf = self._daily_cf(np.linspace(0.8, 1.4, 200))
        clim = QaQc._ebc_cf_climatology(cf)
        doy = cf.index.dayofyear
        day = cf.index[100]
        expected = cf[(doy >= day.dayofyear - 5) & (doy <= day.dayofyear + 5)]
        assert np.isclose(clim.loc[day.dayofyear], expected.mean())
        # no data in January
        assert np.isnan(clim.loc[15])

    def test_ebc_cf_out_of_limits_falls_through(self, data, monkeypatch):
        # method 3 values outside of 0.5-2 are not used, method 4 fills
        config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'US-Tw3_config.ini'
        q = QaQc(Data(config))
        q.correct_data(et_gap_fill=False)
        method_3_days = q.df.ebc_cf_method == 3
        assert method_3_days.any()
        monkeypatch.setattr(
            QaQc, '_ebc_cf_adjacent_years',
            staticmethod(lambda cf, half_win=5: cf * 0 + 2.5)
        )
        q = QaQc(Data(config))
        q.correct_data(et_gap_fill=False)
        methods = q.df.ebc_cf_method
        assert not (methods == 3).any()
        assert (methods[method_3_days] == 4).all()

    def test_ebr_correction_outputs(self, data):
        config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'US-Tw3_config.ini'
        q = QaQc(Data(config))
        q.correct_data(et_gap_fill=False)
        df = q.df.rename(columns=q.inv_map)
        methods = df.ebc_cf_method.dropna()
        assert set(methods.unique()).issubset({1, 2, 3, 4})
        cf = df.ebc_cf.dropna()
        assert ((cf > 0.5) & (cf < 2)).all()
        assert np.allclose(df.LE_corr, df.LE * df.ebc_cf, equal_nan=True)
        assert np.allclose(df.ebr_corr, 1 / df.ebc_cf, equal_nan=True)
        # method flag is daily only
        assert 'ebc_cf_method' not in q.monthly_df.columns

    def _qc_example_config(self, data, tmp_path, **kwargs):
        """Daily US-AR1 example copied to a temp dir (gridMET writes)"""
        config = (
            data['package_root_dir'] /
            'examples' /
            'Config_options' /
            'config_for_QC_flag_filtering.ini'
        )
        return _tmp_config(config, tmp_path, **kwargs)

    def test_lin_regress_correction(self, data, tmp_path):
        q = QaQc(Data(self._qc_example_config(data, tmp_path)))
        q.correct_data(meth='lin_regress', et_gap_fill=False)
        results = q.lin_regress_results
        assert q.corrected
        assert q.corr_meth == 'lin_regress'
        assert 0 < results['r2 (coef. det.)'].iloc[0] <= 1
        assert 'ET_corr' in q.df.rename(columns=q.inv_map).columns

    def test_write_daily_and_monthly(self, data, tmp_path):
        # monthly resampling uses the 'ME' alias which needs pandas >= 2.2
        q = QaQc(Data(self._qc_example_config(data, tmp_path)))
        q.correct_data(et_gap_fill=False)
        out_dir = tmp_path / 'output'
        q.write(out_dir=out_dir)
        daily = pd.read_csv(
            out_dir / f'{q.site_id}_daily_data.csv', index_col=0
        )
        monthly = pd.read_csv(
            out_dir / f'{q.site_id}_monthly_data.csv', index_col=0
        )
        assert len(daily) == len(q.df)
        assert len(monthly) == len(q.monthly_df)
        assert {'ET', 'ET_corr', 'ebr_corr'}.issubset(monthly.columns)

    def test_QaQc_plots(self, data, tmp_path):
        q = QaQc(Data(self._qc_example_config(data, tmp_path)))
        q.correct_data(et_gap_fill=False)
        out_file = tmp_path / 'plots.html'
        q.plot(out_file=out_file)
        assert q.plot_file == out_file
        assert out_file.is_file()

    def test_invalid_refET(self, data, tmp_path):
        q = QaQc(Data(self._qc_example_config(data, tmp_path)))
        with pytest.raises(ValueError):
            q.correct_data(refET='etr')

    def test_daily_ASCE_refET_anemometer_height_from_config(self, data):
        config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'US-Tw3_config.ini'
        q = QaQc(Data(config))
        q.config.set('METADATA', 'anemometer_height', '2')
        q.daily_ASCE_refET()
        eto = q.df.rename(columns=q.inv_map).ASCE_ETo
        q.daily_ASCE_refET(anemometer_height=2.0)
        eto_float = q.df.rename(columns=q.inv_map).ASCE_ETo
        assert np.allclose(eto, eto_float, equal_nan=True)
        assert 4 < eto.mean() < 5
        with pytest.raises(ValueError):
            q.daily_ASCE_refET(reference='medium')

    def test_download_gridMET_single_variable(
            self, data, tmp_path, monkeypatch):
        monkeypatch.setattr(xarray, 'open_dataset', _fake_gridMET_server())
        config = self._qc_example_config(data, tmp_path)
        config_before = config.read_text().splitlines()
        q = QaQc(Data(config))
        # a single name used to be split into characters
        q.download_gridMET('ETr')
        df = q.df
        assert 'gridMET_ETr' in df.columns
        assert 'gridMET_ETo' not in df.columns
        assert np.isclose(df.gridMET_ETr.mean(), 6.0)
        # downloading again replaces the columns
        q.download_gridMET('ETr')
        assert list(q.df.columns).count('gridMET_ETr') == 1
        # path saved in the config is relative to the config file
        saved = q.config.get('METADATA', 'gridMET_file_path')
        assert not Path(saved).is_absolute()
        assert (config.parent / saved).is_file()
        # only the gridMET line was added, comments and layout are kept
        config_after = config.read_text().splitlines()
        added = [l for l in config_after if l not in config_before]
        assert added == ['gridMET_file_path = {}'.format(saved)]
        assert [l for l in config_after if l not in added] == config_before

    def test_gridMET_gap_fill_reuses_saved_file(
            self, data, tmp_path, monkeypatch):
        monkeypatch.setattr(xarray, 'open_dataset', _fake_gridMET_server())
        config = self._qc_example_config(data, tmp_path)
        q = QaQc(Data(config))
        q.correct_data()
        df = q.df.rename(columns=q.inv_map)
        assert {'ET_fill', 'ET_gap', 'ETrF_filtered'}.issubset(df.columns)
        # second run reads the saved file through the relative path
        monkeypatch.setattr(
            xarray, 'open_dataset', _fake_gridMET_server(fail=True)
        )
        q2 = QaQc(Data(config))
        assert q2.gridMET_exists
        q2.correct_data()
        df2 = q2.df.rename(columns=q2.inv_map)
        assert np.allclose(df.ET_corr, df2.ET_corr, equal_nan=True)

    def test_gridMET_download_failure_skips_gap_fill(
            self, data, tmp_path, monkeypatch):
        monkeypatch.setattr(
            xarray, 'open_dataset', _fake_gridMET_server(fail=True)
        )
        q = QaQc(Data(self._qc_example_config(data, tmp_path)))
        q.correct_data()
        df = q.df.rename(columns=q.inv_map)
        assert 'ET_corr' in df.columns
        assert 'ET_fill' not in df.columns
        assert not q.gridMET_exists


class TestValueFixes(object):
    """Fixes in version 0.4.0 that change calculated values"""

    def _synthetic_config(self, tmp_path, n_days=4):
        """
        Half-hourly synthetic data and config file in a temp directory.

        Precipitation is missing for all of the third day and LE and
        precipitation are both missing for the same 5 records on the second
        day.
        """
        index = pd.date_range('2015-06-01', periods=48 * n_days, freq='30min')
        hour = index.hour + index.minute / 60
        sun = np.clip(np.sin((hour - 6) / 12 * np.pi), 0, None)
        df = pd.DataFrame({
            'date': index.strftime('%Y%m%d%H%M'),
            'Rn': 600 * sun - 50, 'G': 60 * sun - 5,
            'LE': 300 * sun, 'H': 150 * sun, 'TA': 20 + 10 * sun,
            'P': 0.1,
        })
        df.loc[60:64, ['LE', 'P']] = np.nan
        df.loc[96:143, 'P'] = np.nan
        df.to_csv(tmp_path / 'synthetic.csv', index=False, na_rep='-9999')
        config = tmp_path / 'synthetic_config.ini'
        config.write_text(
            '[METADATA]\n'
            'climate_file_path = synthetic.csv\n'
            'station_latitude = 39.5\nstation_longitude = -119.8\n'
            'station_elevation = 1500\nmissing_data_value = -9999\n'
            'date_parser = %Y%m%d%H%M\nsite_id = synthetic\n'
            '[DATA]\n'
            'datestring_col = date\n'
            'net_radiation_col = Rn\nnet_radiation_units = w/m2\n'
            'ground_flux_col = G\nground_flux_units = w/m2\n'
            'latent_heat_flux_col = LE\nlatent_heat_flux_units = w/m2\n'
            'sensible_heat_flux_col = H\nsensible_heat_flux_units = w/m2\n'
            'avg_temp_col = TA\navg_temp_units = C\n'
            'precip_col = P\nprecip_units = mm\n'
        )
        return config

    def test_daily_sums_keep_missing_days(self, tmp_path):
        config = self._synthetic_config(tmp_path)
        q = QaQc(Data(config), drop_gaps=False, max_interp_hours=None)
        df = q.df.rename(columns=q.inv_map)
        # day with no precipitation records is missing, not zero
        assert np.isnan(df.ppt.iloc[2])
        assert np.isclose(df.ppt.iloc[0], 4.8)
        assert np.isclose(df.ppt.iloc[1], 4.3)
        # gaps in LE are counted even where precipitation is also missing
        assert df.LE_subday_gaps.iloc[1] == 5

    def test_t_min_t_max_filtered_like_t_avg(self, tmp_path):
        config = self._synthetic_config(tmp_path)
        cp = configparser.ConfigParser(interpolation=None)
        cp.read(config)
        # remove most of day 3 air temperature records
        df = pd.read_csv(tmp_path / 'synthetic.csv')
        df.loc[100:130, 'TA'] = -9999
        df.to_csv(tmp_path / 'synthetic.csv', index=False)
        q = QaQc(Data(config))
        df = q.df.rename(columns=q.inv_map)
        assert np.isnan(df.t_avg.iloc[2])
        assert np.isnan(df.t_min.iloc[2]) and np.isnan(df.t_max.iloc[2])
        assert np.isclose(df.t_max.iloc[0], 30)
        # also calculated without sub-daily interpolation
        q = QaQc(Data(config), max_interp_hours=None)
        assert 't_min' in q.df.rename(columns=q.inv_map).columns

    def test_daily_input_frequency_not_column_count(self, data, tmp_path):
        # daily FLUXNET example with fewer than 24 columns loaded
        config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'fluxnet_config.ini'
        keep = [
            'datestring_col', 'net_radiation', 'ground_flux',
            'latent_heat_flux', 'sensible_heat_flux', 'vap_press_def',
            'avg_temp'
        ]
        cp = configparser.ConfigParser(interpolation=None)
        cp.read(config)
        remove = [
            ('DATA', k) for k in cp['DATA']
            if not any(k.startswith(n) for n in keep) or k.endswith('_qc')
        ]
        config = _tmp_config(config, tmp_path, remove_opts=remove)
        d = Data(config)
        assert len(d.df.columns) < 24
        # vapor pressure is calculated from daily VPD and air temperature
        assert 'vp' in d.df.rename(columns=d.inv_map).columns
        # hourly reference ET is not calculated from daily data
        assert d.hourly_ASCE_refET() is None

    def test_bowen_ratio_correction_limits(self, data):
        config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'US-Tw3_config.ini'
        q = QaQc(Data(config))
        q.correct_data(meth='br', et_gap_fill=False)
        df = q.df.rename(columns=q.inv_map)
        LE_corr = df.LE_corr.dropna()
        assert np.isfinite(LE_corr).all()
        assert ((LE_corr > -100) & (LE_corr < 850)).all()
        assert np.isfinite(df.br.dropna()).all()

    def test_lin_regress_energy_corr(self, data):
        config = data['package_root_dir']\
            /'examples'/'Basic_usage'/'US-Tw3_config.ini'
        q = QaQc(Data(config))
        q.correct_data(meth='lin_regress', et_gap_fill=False)
        df = q.df.rename(columns=q.inv_map)
        # Rn is the dependent variable so it is not corrected
        assert np.allclose(df.energy_corr, df.Rn - df.G_corr, equal_nan=True)

    def test_gridMET_outside_domain(self, data, tmp_path, monkeypatch):
        monkeypatch.setattr(xarray, 'open_dataset', _fake_gridMET_server())
        config = data['package_root_dir']/'examples'\
            /'Config_options'/'config_for_QC_flag_filtering.ini'
        # station far north of the (fake) grid
        config = _tmp_config(
            config, tmp_path,
            set_opts=[('METADATA', 'station_latitude', '50.0')]
        )
        q = QaQc(Data(config))
        q.download_gridMET()
        assert not any(c.startswith('gridMET') for c in q.df.columns)
        assert not q.gridMET_exists


class TestUtil(object):

    def test_monthly_resample_threshold(self):
        # June has 30 days, 24 days is exactly 80 percent
        index = pd.date_range('2015-06-01', '2015-07-31', freq='D')
        df = pd.DataFrame({'ET': 2.0}, index=index)
        df.loc['2015-06-25':'2015-06-30', 'ET'] = np.nan
        df.loc['2015-07-25':'2015-07-31', 'ET'] = np.nan
        monthly = util.monthly_resample(df, ['ET'], 'sum', 0.8)
        # missing days are filled with the mean before summation
        assert np.isclose(monthly.ET.iloc[0], 60.0)
        # July has 24 of 31 days, less than 80 percent
        assert np.isnan(monthly.ET.iloc[1])

    def test_convert_f_to_c(self):
        df = pd.DataFrame({'t_avg': [32., 212., -40.]})
        df = Convert.convert('t_avg', 'f', 'c', df)
        assert np.allclose(df.t_avg, [0., 100., -40.])

    def test_set_config_option(self, tmp_path):
        config = tmp_path / 'config.ini'
        config.write_text(
            '# a comment to keep\n'
            '[METADATA]\n'
            'site_id = A\n'
            'GRIDMET_FILE_PATH = /old/path.csv\n'
            '\n'
            '# data comment\n'
            '[DATA]\n'
            'net_radiation_col = Rn\n'
        )
        # replace existing option, option names are not case sensitive
        util.set_config_option(
            config, 'METADATA', 'gridMET_file_path', 'grid.csv'
        )
        # add a new option at the end of a section, before the next one
        util.set_config_option(config, 'METADATA', 'skiprows', '2')
        lines = config.read_text().splitlines()
        assert lines == [
            '# a comment to keep',
            '[METADATA]',
            'site_id = A',
            'gridMET_file_path = grid.csv',
            'skiprows = 2',
            '',
            '# data comment',
            '[DATA]',
            'net_radiation_col = Rn',
        ]
        cp = configparser.ConfigParser(interpolation=None)
        cp.read(config)
        assert cp.get('METADATA', 'gridMET_file_path') == 'grid.csv'
        with pytest.raises(ValueError):
            util.set_config_option(config, 'NOPE', 'a', 'b')

    def test_convert_unit_aliases(self):
        df = pd.DataFrame({'co2': [400., 410.], 'zeta': [-0.1, 0.2]})
        df = Convert.convert('co2', 'ppm', 'umol/mol', df)
        df = Convert.convert('zeta', 'nondimensional', 'dimensionless', df)
        assert df.co2.tolist() == [400., 410.]
        assert df.zeta.tolist() == [-0.1, 0.2]

    @pytest.mark.parametrize(
        'latitude, longitude, expected',
        [
            (40.7584, -82.5154, -5.0),
            (39.5296, -119.8138, -8.0),
            (33.4484, -112.0740, -7.0),
            (28.6139, 77.2090, 5.5),
        ]
    )
    def test_standard_utc_offset(
            self, latitude, longitude, expected):
        offset = util.standard_utc_offset(
            latitude,
            longitude
        )

        assert np.isclose(offset, expected)


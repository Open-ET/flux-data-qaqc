# -*- coding: utf-8 -*-

import pkg_resources
import pytest
from pathlib import Path
from shutil import move, copy, rmtree
from refet.calcs import _ra_daily, _ra_hourly

import numpy as np
import pandas as pd

from fluxdataqaqc import Data
from fluxdataqaqc import QaQc
from fluxdataqaqc import Plot
from fluxdataqaqc import Convert
from fluxdataqaqc import util


@pytest.fixture(scope="session")
def data(request):
    """Prepare input data used for tests"""
    d = {}
    package_root_dir = Path(
        pkg_resources.resource_filename(
            'fluxdataqaqc', '.'
        )
    ).parent

    d['package_root_dir'] = package_root_dir      
    # use example data from tutorial for most tests
    examples_dir = package_root_dir / 'examples'
    test_data_files = [
        f for f in examples_dir.rglob('*') if f.suffix in [
            '.csv','.ini','.xlsx'] and not 'gridMET_data' in str(f)
    ]
    
    temp_data_dir = Path('tests') / 'test_data'
    if not temp_data_dir.is_dir():
        temp_data_dir.mkdir(parents=True, exist_ok=True)

    print(temp_data_dir, Path())

    for f in test_data_files:
        copy(f, temp_data_dir)




    def teardown():
        rmtree(temp_data_dir)

    request.addfinalizer(teardown)

    return d
    
    
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
        assert self.data_obj._df == None
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
        #assert np.isclose(ts.mean(), 0.1830681034395387)

    def test_Data_plots(self):
        assert self.data_obj.plot_file == None
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
        g_mean = df[['g_1', 'g_2', 'g_3','g_4']].mean(1)
        g_weighted_mean = df.G
        assert (g_mean != g_weighted_mean).any()

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


class TestUtil(object):

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


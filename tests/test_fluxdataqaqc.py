# -*- coding: utf-8 -*-

import pkg_resources
import pytest
from pathlib import Path
from shutil import move, copy, rmtree

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
        assert len(df.dropna()) == 68035

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



        

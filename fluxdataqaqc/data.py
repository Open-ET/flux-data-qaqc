# -*- coding: utf-8 -*-
"""
Read and load in fluxnet or other climate time series files.
"""

import configparser as cp
import numpy as np
import pandas as pd
from pathlib import Path


class Data(object):
    """
    Read a climate station time series file using a config file.

    Upon reading a config file for a specific climate station dataset the 
    ``Data`` class stores metadata on climate variables in file, their units, 
    and file path. Methods include tools to load climate time series data into a
    :obj:`pandas.DataFrame` and convert units if needed.

    Note:
        Path to climate time series file specified in the config_file should
        be relative to the config file's location.

    TODO:
     * verify/update input variable units
    """

    # maps internal names to config names for variables
    # updated with user's column names in instance attribute "variables"
    variable_names_dict = {
        'date' : 'datestring_col',
        'year' : 'year_col',
        'month' : 'month_col',
        'day' : 'day_col',
        'Rn' : 'net_radiation_col',
        'G' : 'ground_flux_col',
        'LE' : 'latent_heat_flux_col',
        'LE_corr' : 'latent_heat_flux_corrected_col',
        'H' : 'sensible_heat_flux_col',
        'H_corr' : 'sensible_heat_flux_corrected_col',
        'sw_in' : 'shortwave_in_col',
        'sw_out' : 'shortwave_out_col',
        'sw_pot' : 'shortwave_pot_col',
        'lw_in' : 'longwave_in_col',
        'lw_out' : 'longwave_out_col',
        'vp' : 'vap_press_col',
        'vpd' : 'vap_press_def_col',
        't_avg' : 'avg_temp_col',
        'ppt' : 'precip_col',
        'ws' : 'wind_spd_col'
    }

    def __init__(self, config):

        self.config_file = Path(config).absolute()
        self.config = self._load_config(self.config_file)
        self.variables = self.get_config_vars()
        self.na_val = self.config.get('METADATA', 'missing_data_value')
        self.elevation = int(self.config.get('METADATA', 'station_elevation'))
        self.latitude = float(self.config.get('METADATA', 'station_latitude'))
        self.climate_file = self._get_climate_file()
        self.header = self._get_header(self.climate_file)
        self.qc_var_pairs = self._get_qc_flags()
        self.site_id = self.config.get('METADATA', 'site_id')
        # output dir will be in current working directory
        self.out_dir = Path(self.site_id + '_output').absolute()
        self._df = None

    def _load_config(self, config_file):
        if not config_file.is_file():
            raise FileNotFoundError('ERROR: config file not found')
        config = cp.ConfigParser()
        config.read(config_file)
        return config

    def _get_climate_file(self):
        """
        Read config file and return absolute path to climate time series file
        """
        file_path = self.config['METADATA']['climate_file_path']
        climate_file = self.config_file.parent.joinpath(file_path)
        
        if not climate_file.is_file():
            err_msg = 'ERROR: climate file:{} does not exist'.format(
                climate_file)
            raise FileNotFoundError(err_msg)
        
        return climate_file

    def _get_header(self, climate_file):
        """
        Read only top line of climate time series file return header names.
        """
        if climate_file.suffix in ('.xlsx', '.xls'):
            workbook = pd.ExcelFile(climate_file)
            rows = workbook.book.sheet_by_index(0).nrows
            header = pd.read_excel(workbook, skipfooter = (rows - 1))
            header = header.columns
        
        else: # assume CSV
            header = np.genfromtxt(
                climate_file, 
                dtype='U', 
                delimiter=',', 
                max_rows=1, 
                autostrip=True,                             
                case_sensitive='lower'
            )
        return header

    def get_config_vars(self):
        """
        Read config data section and get names of all variables, pair with
        internal variable names in a dictionary.

        Also parses config file for optionally added multiple soil heat flux
        and soil moisture variables if given following the naming convention
        explained in ``flux-data-qaqc``.

        Arguments:
            None

        Returns:
            variables (dict): dictionary with internal climate variable names
                as keys and user's names as found in config as values. 
        """

        variables = {}
        # get all variables found in Data.variable_names_dict
        for k, v in Data.variable_names_dict.items():
            if self.config.has_option('DATA', v):
                variables[k] = self.config.get('DATA', v)
            else:
                variables[k] = None

        # get multiple G flux/soil moisture variables 
        all_keys = dict(self.config.items('DATA')).keys()
        # should be named as 'g_ or 'theta_ what comes after not strict yet
        added_g_keys = []
        added_theta_keys = []
        for k in all_keys:
            if k.startswith('g_') and not k.endswith('_units'):
                added_g_keys.append(k)
            if k.startswith('theta_') and not k.endswith('_units'):
                added_theta_keys.append(k)
                
        if added_g_keys:
            for el in added_g_keys:
                variables[el] = self.config.get('DATA', el)
        if added_theta_keys:
            for el in added_theta_keys:
                variables[el] = self.config.get('DATA', el)
        
        # update all None values in config with 'na' in variables dict
        for k, v in variables.items():
            if v is None:
                variables[k] = 'na'

        return variables
    
    def _get_qc_flags(self):
        """
        Process any existing QC flags for variables in config, also add
        key,val to variables attribute.
        """
        qc_var_pairs = {}
        tmp = {}
        for k,v in self.variables.items():
            qc_var = '{}_QC'.format(v)
            if qc_var in self.header:
                internal_name = '{}_qc_flag'.format(k)
                tmp[internal_name] = qc_var
                qc_var_pairs[v] = qc_var
               
        self.variables.update(tmp)
        return qc_var_pairs

    def apply_qc_flags(self, threshold = 0.5):
        """
        Use provided QC flags for climate variables to exclude bad data 
        by forcing them to null values. 
        
        Specifically If the QC flag is < `threshold` change the variables 
        value for that day to null. For FLUXNET datasets the QC flag for 
        daily data is a fraction between 0 and 1 indicating the percentage 
        of measured or good quality gap filled data used.

        Keyword Arguments:
            threshold (float): default 0.5. Threshold for QC flag, if flag
            is below threshold replace that variables value with null.

        Returns:
            None

        """

        # load dataframe if not yet accessed
        df = self.df
        # loop over each variable that has a provided qc flag and set nulls
        # where qc flag < threshold
        for var, flag in self.qc_var_pairs.items():
            df.loc[
                (df[flag] < threshold) & (df[flag].notnull()) , var
            ] = np.nan

        self._df = df

    @property
    def df(self):
        """
        Pulls energy balance variables out of config file and attempts to load 
        them into a date-indexed pandas.DataFrame. 

        Also read in any QC flag columns from input climate file if they exist
        for later cleaning of data. The QC flags should have the same name as
        the variables being used which are set in the config file. For example,
        for FLUXNET data the gap filled variable for LE is LE_F_MDS and the
        QC flag for that variable is LE_F_MDS_QC.        

        Returns:
            df (pandas.DataFrame): dataframe of all variables specified in 
                config file, if not found they will be filled with NaNs.
        """

        # avoid overwriting pre-assigned data
        if isinstance(self._df, pd.DataFrame):
            return self._df

        # handle missing 'na' data
        # loop for debugging only
        #for k,v in variables.items():
        #    if v == 'na':
        #        print('WARNING: {} is missing from input data'.format(k))

        variables = self.variables
        vars_notnull = dict((k, v) for k, v in variables.items() if v != 'na')
        cols = list(vars_notnull.values())

        missing_cols = None
        if not set(cols).issubset(self.header):
            missing_cols = set(cols) - set(self.header)
            err_msg = ('WARNING: the following config variables are missing '
                'in the input climate file:\n{}\nThey will be filled with '
                'NaN values'.format(' '.join(missing_cols)))
            print(err_msg)
            cols = set(cols).intersection(self.header)

        # load data file depending on file format
        if self.climate_file.suffix in ('.xlsx', '.xls'):
            # find indices of headers we want, only read those, excel needs ints
            ix=[i for i, e in enumerate(self.header) if e in set(cols)]
            df = pd.read_excel(
                self.climate_file,
                parse_dates = [variables.get('date')],
                usecols = ix,
                na_values=['NaN', 'NAN', '#VALUE!', self.na_val]
            )
        else:
            df = pd.read_csv(
                self.climate_file,
                parse_dates = [variables.get('date')],
                usecols = cols,
                na_values=['NaN', 'NAN', '#VALUE!', self.na_val]
            )

        if missing_cols:
            df = df.reindex(columns=list(cols)+list(missing_cols))

        # the only column that is always renamed is the datestring_col
        df.rename(columns={variables['date']: 'date'}, inplace=True)
        
        # date index
        df.index = df.date
        df.drop('date', axis=1, inplace=True)
        self._df = df
        return df

    @df.setter
    def df(self, data_frame):
        if not isinstance(data_frame, pd.DataFrame):
            raise TypeError("Must assign a Pandas.DataFrame object")
        self._df = data_frame


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
        'LE_user_corr' : 'latent_heat_flux_corrected_col',
        'H' : 'sensible_heat_flux_col',
        'H_user_corr' : 'sensible_heat_flux_corrected_col',
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
        self.units = self.get_config_units()
        self.inv_map = {v: k for k, v in self.variables.items()}
        self.na_val = self.config.get('METADATA', 'missing_data_value')
        self.elevation = int(self.config.get('METADATA', 'station_elevation'))
        self.latitude = float(self.config.get('METADATA', 'station_latitude'))
        self.climate_file = self._get_climate_file()
        self.header = self._get_header(self.climate_file)
        self.soil_var_weight_pairs = self._get_soil_var_avg_weights()
        self.qc_var_pairs = self._get_qc_flags()
        self.site_id = self.config.get('METADATA', 'site_id')
        # get optionally defined QC thresholds (numeric) or flags (string)
        # for filtering bad data
        if 'qc_threshold' in dict(self.config.items('METADATA')):
            self.qc_threshold=float(self.config.get('METADATA','qc_threshold'))
        else:
            self.qc_threshold = None
        # optional comma separated string of QC flags
        if 'qc_flag' in dict(self.config.items('METADATA')):
            self.qc_flag = [
                el.strip() for el in self.config.get(
                    'METADATA','qc_flag').split(',')
            ]
        else:
            self.qc_flag = None
        # output dir will be in directory of config file
        self.out_dir = self.config_file.parent / 'output'
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

    def _get_soil_var_avg_weights(self):
        """
        Read config and variable attribute to match user names for multiple
        soil moisture/heat flux variables to user defined weights used for
        calculating weighted/non-weighted average values in `Data.df`
        """
        soil_var_weight_pairs = {}
        all_keys = dict(self.config.items('DATA')).keys()

        for k,v in dict(self.config.items('DATA')).items():
            weight_name = '{}_weight'.format(k)
            var_name = self.variables.get(k)
            if k in self.variables.keys() and weight_name in all_keys:
                weight = self.config.get('DATA', weight_name)
                tmp = {'name': var_name, 'weight' : weight}
                soil_var_weight_pairs[k] = tmp
                #print(k,v, var_name, weight)
            # update dict with weights if not specified in config, normalized
            # later in df
            elif (k.startswith('g_') or k.startswith('theta')) and not \
                    weight_name in all_keys and not k.endswith('_units') and \
                    not k.endswith('_weight'):
                tmp = {'name': var_name, 'weight' : 1}
                soil_var_weight_pairs[k] = tmp

        return soil_var_weight_pairs

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
            if k.startswith('g_') and not k.endswith('_units') \
                    and not k.endswith('_weight'):
                added_g_keys.append(k)
            if k.startswith('theta_') and not k.endswith('_units') \
                    and not k.endswith('_weight'):
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
    
    def get_config_units(self):
        """
        Get units from config file, pair with var names and store in dictionary.

        Keys are `flux-data-qaqc` variable names, i.e. the same names used 
        by :attr:`Data.variables` and values are strings assigned in the config
        file. 

        Arguments:
            None

        Returns:
            units (dict): dictionary with `flux-data-qaqc` variable names as
                keys and user's units for each as values.

        Note:
            Parsing of correct units and conversion if needed is performed
            in the :obj:`fluxdataqaqc.QaQc` class. Also, if units are not given
            in the config file a warning message is printed and the units are
            not included and thus will either need to be manually added later
            e.g. in Python by adding to :attr:`Data.units` or by adding them
            to the config and recreating a :obj:`Data` object otherwise the 
            units will remain unknown and not be able to be later converted.
        """
        no_unit_vars = ('datestring_col', 'year_col', 'month_col', 'day_col')
        # dictionary that maps config unit keys to units
        units_config = {
            v.replace('_col', '_units'): None for k,v in 
            self.variable_names_dict.items() if not v in no_unit_vars
        }
        # add user multiple g or soil moisture var units config names
        for k,v in self.variables.items():
            # if multiple g uses same var assigned to ground_flux_col units the 
            # added G var will not be included/duplicated
            if k.startswith('g_') and not self.variables[k] == 'G':
                units_config['{}_units'.format(k)] = None
            if k.startswith('theta_'):
                units_config['{}_units'.format(k)] = None

        config_dict = dict(self.config.items('DATA'))
        for k in units_config: 
            if k in config_dict:
                units_config[k] = config_dict[k]
            else:
                print(
                    'WARNING: units for var {} missing from the config file'\
                        .format(k.replace('_units',''))
                )
                
        inv_names = {v:k for k,v in self.variable_names_dict.items()}

        # dictionary that maps fluxdataqaqc var names to units
        units = dict()
        for k in units_config:
            k_col = k.replace('_units', '_col')
            if k_col in self.variable_names_dict.values():
                units[inv_names[k_col]] = units_config.get(k)
            else:
                # for multiple g and theta remove _units suffix
                name = k.replace('_units','')
                units[name] = units_config.get(k)

        return units

    def _get_qc_flags(self):
        """
        Process any existing QC flags for variables in config, also add
        key,val to variables attribute if QC flags are not specified in config
        and follow the naming convention [var_name]_QC in the header of the 
        data file.
        """
        qc_var_pairs = {}
        tmp = {}

        no_qc_vars = ('datestring_col', 'year_col', 'month_col', 'day_col')
        # dictionary that maps config QC values to keys for main variables
        # other variables like multiple g or theta (with unknown names) are
        # search for in loop below
        qc_config = {
            v.replace('_col', '_qc'): k for k,v in 
            self.variable_names_dict.items() if not v in no_qc_vars
        }
        # first look if specified in config
        for k,v in self.config.items('DATA'):
            if k in qc_config:
                # internal name for the variable (e.g. LE or Rn)
                var_name = qc_config[k]
                user_var_name = self.variables.get(var_name)
            # keys are internal names for multiple G and theta
            elif k.startswith(('g_','theta_')) and k.endswith('_qc'):
                var_name = k 
            # key is not for a QC flag header name...
            else:
                continue
            if not v in self.header:
                print('WARNING: {} quality control name specified in the config'
                    ' file for variable: {} does not exist in the input file, '
                    'it will not be used.'.format(v, self.variables[var_name])
                )
                continue
            internal_name = '{}_qc_flag'.format(var_name)
            tmp[internal_name] = v
            qc_var_pairs[user_var_name] = v


        # also look in header, currently always gets these as well, may change 
        # find vairable names that have a matching name with '_QC' suffix
        # if found here but different name in config use the name in the config 
        for k,v in self.variables.items():
            qc_var = '{}_QC'.format(v)
            if qc_var in self.header:
                internal_name = '{}_qc_flag'.format(k)
                if internal_name in tmp:
                    continue
                tmp[internal_name] = qc_var
                qc_var_pairs[v] = qc_var
               
        self.variables.update(tmp)
        return qc_var_pairs

    def apply_qc_flags(self, threshold=None, flag=None):
        """
        Use provided QC values or flags for climate variables to filter bad 
        data by converting them to null values, updates the :attr:`Data.df`. 
        
        Specifically where the QC value is < `threshold` change the variables 
        value for that date-time to null. For FLUXNET datasets the QC value for 
        daily data is a fraction between 0 and 1 indicating the percentage 
        of measured or good quality gap filled data used. The other option
        is to use a column of flags, e.g. 'x' for bad data. 

        The threshold value or flag may be also specified in the config file.

        Keyword Arguments:
            threshold (float): default None. Threshold for QC values, if flag
                is below threshold replace that variables value with null.
            flag (str, list, or tuple): default None. Character flag signifying 
                bad data to filter out. Can be list or tuple of multiple flags.

        Returns:
            None

        """
        # if QC threshold or flags not passed use values from config if exist
        if not threshold:
            threshold = self.qc_threshold
        if not flag:
            flag = self.qc_flag
        # load dataframe if not yet accessed
        df = self.df
        # infer each columns datatype to avoid applying thresholds to strings
        infer_type = lambda x: pd.api.types.infer_dtype(x, skipna=True)
        df_types = pd.DataFrame(df.apply(
            pd.api.types.infer_dtype, axis=0, skipna=True
            )
        ).rename(columns={'index': 'index', 0: 'type'})
        # loop over each variable that has a provided qc value and set nulls
        # where qc value < threshold, qc_var_pairs maps var names to qc names
        if threshold:
            for var, qc in self.qc_var_pairs.items():
                if df_types.loc[qc,'type'] is not 'string':
                    df.loc[
                        (df[qc] < threshold) & (df[qc].notnull()) , var
                    ] = np.nan
        # set values to null where flag is a certain string
        if flag:
            if isinstance(flag, str):
                for var, qc in self.qc_var_pairs.items():
                    if df_types.loc[qc,'type'] is 'string':
                        df.loc[
                            (df[qc] == flag) & (df[qc].notnull()) , var
                        ] = np.nan
            # apply multiple character flags
            elif isinstance(flag, (list,tuple)):
                for f in flag:
                    for var, qc in self.qc_var_pairs.items():
                        if df_types.loc[qc,'type'] is 'string':
                            df.loc[
                                (df[qc] == f) & (df[qc].notnull()) , var
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
        
        def calc_weight_avg(d, pref, df):
            """
            Helper function to reduce redundant code for calculating weighted
            average currently for multiple soil heat flux and moisture variables

            d is soil_var_weight_pairs dict
            pref is variable prefix str, i.e. g_ or theta_
            df is the dataframe 
            """
            # list of multiple variables with prefix
            if d:
                vs = [d.get(e) for e in d if e.startswith(pref)]
            else:
                vs = []
            # soil heat flux weighted average
            if len(vs) > 1: # if 1 or less no average
                total_weights = np.sum([float(e.get('weight')) for e in vs])
                # if weights are not normalized update them
                if not np.isclose(total_weights, 1.0):
                    print(
                        '{} weights do not sum to one, normalizing'\
                            .format(pref.replace('_',''))
                    )
                    for k,v in d.items():
                        if k.startswith(pref):
                            nw = float(v.get('weight')) / total_weights
                            d[k] = {'name': v.get('name'), 'weight': nw}
                            msg = ', '.join([
                                '{}:{:.2f}'.format(
                                    v.get('name'),float(v.get('weight'))
                                    )
                                    for k,v in d.items() if k.startswith(pref)
                                ]
                            )
                    print('Here are the new weights:\n', msg)

                # calculate average, use updated weights if calculated
                vs = [d.get(e) for e in d if e.startswith(pref)]
                tmp_df = df[[e.get('name') for e in vs]].copy()
                for pair in vs:
                    tmp_df[pair.get('name')] *= float(pair.get('weight'))
                df['{}mean'.format(pref)] = tmp_df.sum(axis=1)
                # if calculated update variables
                self.variables['{}mean'.format(pref)] = '{}mean'.format(pref)

        # calculate weighted average soil vars if they exist
        d = self.soil_var_weight_pairs
        calc_weight_avg(d, 'g_', df)
        calc_weight_avg(d, 'theta_', df)

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


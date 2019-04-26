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
     * add handling for different formats of string dates
     * handling of climate data files of non- .xlsx format (i.e. not fluxnet)
     * metadata variable missing_data_value is read in df, maybe reorganize
     * store other metadata values from config file if needed?
     * check how pd.read_csv's parse_dates works with discrete y/m/d columns
    """

    def __init__(self, config):

        self.config_file = Path(config).absolute()
        self.config = self.load_config(self.config_file)
        self.climate_file = self._get_climate_file()
        self.header = self._get_header(self.climate_file)
        self._df = None

    def load_config(self, config_file):
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
        Read top line of climate time series file return header names.
        """
        header_data = np.genfromtxt(
            climate_file, 
            dtype='U', 
            delimiter=',', 
            max_rows=1, 
            autostrip=True,                             
            case_sensitive='lower'
        )
        return header_data
    
    @property
    def df(self):
        """
        Pulls energy balance variables out of config file and attempts to load 
        them into a date-indexed pandas.DataFrame. 

        :return: pandas.DataFrame of all the variables and their indices
        """

        # TODO: config parser errors on reading in %'s, if pandas parse_dates
        # fails will have to come up with alternative way of reading in

        # avoid overwriting pre-assigned data
        if isinstance(self._df, pd.DataFrame):
            return self._df

        # need to verify units (modify code below) 
        variables = {}
        variables['date'] = self.config['DATA']['datestring_col']
        variables['year'] = self.config['DATA']['year_col']
        variables['month'] = self.config['DATA']['month_col']
        variables['day'] = self.config['DATA']['day_col']
        variables['Rn'] = self.config['DATA']['net_radiation_col']
        variables['G'] = self.config['DATA']['ground_flux_col']
        variables['LE'] = self.config['DATA']['latent_heat_flux_col']
        variables['LE_corr'] = self.config['DATA']\
            ['latent_heat_flux_corrected_col']
        variables['H'] = self.config['DATA']['sensible_heat_flux_col']
        variables['H_corr'] = self.config['DATA']\
            ['sensible_heat_flux_corrected_col']

        # handle missing 'na' data
        for k,v in variables.items():
            if v == 'na':
                print('WARNING: {} is missing from input data'.format(k))
        vars_notnull = dict((k, v) for k, v in variables.items() if v != 'na')
        cols = list(vars_notnull.values())

        if not set(cols).issubset(self.header):
            err_msg = ('One of more of the column names below '
                'was not found in the input climate file: {}'.format(
                    ' '.join(cols)))
            raise KeyError(err_msg)

        # parse_dates usually works on most string formats
        df = pd.read_csv(
            self.climate_file,
            parse_dates = [variables.get('date')],
            usecols = cols,
            na_values=['NaN', 'NAN', '#VALUE!', '-9999']
        )

        # rename all df columns to consistent names
        df.rename(
            columns={
               variables['date']: 'date',
               variables['Rn']: 'net_rad', 
               variables['G']: 'g_flux', 
               variables['LE']: 'le_flux', 
               variables['LE_corr']: 'le_flux_corr',
               variables['H']: 'h_flux', 
               variables['H_corr']: 'h_flux_corr'
            }, inplace=True
        )
        
        df.index = df.date
        df.drop('date', axis=1, inplace=True)
        return df

    @df.setter
    def df(self, data_frame):
        if not isinstance(data_frame, pd.DataFrame):
            raise TypeError("Must assign a Pandas.DataFrame object for PRMS data input")
        self._df = data_frame

####### all below not in working shape or not needed currently
#
#    def _obtain_units(self, config_file, provided_variables):
#        """
#        Iterates through the dictionary of all variables and looks up units for only those that are provided.
#        There is no error checking at this stage, checks to see if units are valid is performed later on
#
#        :param config_file: Absolute path to configuration file
#        :param provided_variables: Pandas series (dictionary) that contains the indicies of all provided variables
#
#        :return: A pandas series (dictionary) where the variable names are the keys and the unit strings are the values
#        """
#        config_reader = cp.ConfigParser()
#        config_reader.read(config_file)
#
#        self.units_dict = {}
#
#        # Iterate through all entries in the dictionary, check to see if they are set to -1, otherwise look up the units
#        # in the config file and put them in a dictionary
#        for variable, index in provided_variables.items():
#            if index == -1:  # Variable was not provided, do not look up units
#                pass
#            else:
#                variable_units_string = variable + '_units'
#                units = config_reader['DATA'][variable_units_string]
#
#                self.units_dict.update({variable: units.lower()})
#
#        return self.units_dict
#
#    def _convert_units(self, variable_name, original_values, initial_units):
#        """
#        Determines what variable has been passed, finds the initial units of that variable, and then converts the
#        original data into the appropriate units if necessary
#
#        :param variable_name: String of variable name
#        :param original_values: 1D numpy array of climate observations in original units
#        :param initial_units: pandas series (dictionary) of provided variables and their
#        :return: A numpy array of climate observations that have been converted into the appropriate units
#        """
#
#        unit_string = initial_units[variable_name]
#
#        # Units of radiation
#        if unit_string == 'w/m2':
#            # observations are already in their desired units
#            self.converted_values = original_values
#        elif unit_string == 'langleys' or unit_string == 'ly':
#            # convert langleys to w/m2
#            self.converted_values = np.array(original_values * 0.48458)
#        elif unit_string == 'kw-hr/m2' or unit_string == 'kw_hr/m2' or unit_string == 'kwhr/m2':
#            # convert kilowatt hours to w/m2
#            self.converted_values = np.array((original_values * 1000) / 24)
#        elif unit_string == 'mj/m2':
#            # convert mj/m2 to w/m2
#            self.converted_values = np.array(original_values * 11.574)
#
#        # Units of temperature
#        elif unit_string == 'c' or unit_string == 'celsius':
#            # observations are already in their desired units
#            self.converted_values = original_values
#        elif unit_string == 'f' or unit_string == 'fahrenheit':
#            # convert fahrenheit to celsius
#            self.converted_values = np.array(((original_values - 32.0) * (5.0 / 9.0)))
#        elif unit_string == 'k' or unit_string == 'kelvin':
#            # convert kelvin to celsius
#            self.converted_values = np.array(original_values - 273.15)
#
#        # Units of pressure
#        elif unit_string == 'kpa' or unit_string == 'kilopascals':
#            # observations are already in their desired units
#            self.converted_values = original_values
#        elif unit_string == 'torr' or unit_string == 'mmhg':
#            # convert from torr or mmhg into kilopascals
#            self.converted_values = np.array(original_values * 0.133322)
#        elif unit_string == 'p' or unit_string == 'pascals':
#            # convert from pascals to kilopascals
#            self.converted_values = np.array(original_values / 1000)
#
#        # Units of distance
#        elif unit_string == 'mm' or unit_string == 'millimeters':
#            # observations are already in their desired units
#            self.converted_values = original_values
#        elif unit_string == 'in' or unit_string == 'inches':
#            # convert inches to mm
#            self.converted_values = np.array(original_values * 25.4)
#        elif unit_string == 'ft' or unit_string == 'feet':
#            # convert ft to mm
#            self.converted_values = np.array(original_values * 304.8)
#
#        # Units of speed
#        elif unit_string == 'm/s':
#            # observations are already in their desired units
#            self.converted_values = original_values
#        elif unit_string == 'mph':
#            # convert mph to m/s
#            self.converted_values = np.array(original_values * 0.44704)
#
#        # Units for percentages
#        elif unit_string == '%' or unit_string == 'percent':
#            # observations are already in their desired units
#            self.converted_values = original_values
#        elif unit_string == 'fract' or unit_string == 'fraction':
#            self.converted_values = np.array(original_values * 100)
#
#        # Unrecognized string for variable units, raise an error
#        else:
#            raise ValueError('Parameter {} had an unsupported unit string of {}.'.format(variable_name, unit_string))
#
#        return self.converted_values
#

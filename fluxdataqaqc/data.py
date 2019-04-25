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
#        self.data_valid = self.check_data(self.config_file, self.header_names)
        self.df = self._locate_variables()

    def load_config(self, config_file):
        if not config_file.is_file():
            raise FileNotFoundError('ERROR: config file not found')
        config = cp.ConfigParser()
        config.read(config_file)
        return config

    def _get_climate_file(self):
        """
        Reads config file and returns absolute path to climate time series file
        """
        file_path = self.config['METADATA']['climate_file_path']
        climate_file = Path(file_path).absolute()
        
        if not climate_file.is_file():
            err_msg = 'ERROR: climate file:{} does not exist'.format(
                climate_file)
            raise FileNotFoundError(err_msg)
        
        return climate_file

    def _get_header(self, climate_file):
        # read top line of climate time series file return header names

        header_data = np.genfromtxt(
            climate_file, 
            dtype='U', 
            delimiter=',', 
            max_rows=1, 
            autostrip=True,                             
            case_sensitive='lower'
        )
        return header_data

    def _locate_variables(self):
        """
        Pulls all climate variables out of config file and attempts to find 
        them among the column names of the data file. Uses the function 
        _find_variable_in_header to find matches or report errors.

        :param config_file: string of absolute path to configuration file
        :param header_names:  ndarray of header strings pulled from data file

        :return: Pandas series of all the variables and their indices
        """

        date = self.config['DATA']['datestring_col']
        # TODO: config parser errors on reading in %'s, will have to come up with alternative way of reading in
        # date_format = self.config['DATA']['datestring_format']
        year = self.config['DATA']['year_col']
        month = self.config['DATA']['month_col']
        day = self.config['DATA']['day_col']
        
        # need to verify units (modify code below) 
        Rn = self.config['DATA']['net_radiation_col']
        G = self.config['DATA']['ground_flux_col']
        LE = self.config['DATA']['latent_heat_flux_col']
        LE_corr = self.config['DATA']['latent_heat_flux_corrected_col']  # Pre-corrected values from fluxnet
        H = self.config['DATA']['sensible_heat_flux_col']
        H_corr = self.config['DATA']['sensible_heat_flux_corrected_col']  # Pre-corrected values from fluxnet

        # rename variables at later dev for standard naming convention
        cols = [date, Rn, G, LE, LE_corr, H, H_corr]

        # parse_dates usually works on most string formats
        df = pd.read_csv(self.climate_file, parse_dates=[date],
                         na_values=['NaN', 'NAN', '#VALUE!', '-9999'])[cols]

        # Now rename all df columns to consistent names
        df.rename(columns={date: 'stringdate',
                           Rn: 'net_rad', G: 'g_flux', LE: 'le_flux', LE_corr: 'le_flux_corr',
                           H: 'h_flux', H_corr: 'h_flux_corr'}, inplace=True)

        return df

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
#    def data_frame(self, config_file, climate_file, provided_variables, initial_units):
#        """
#        :param config_file:
#        :param climate_file:
#        :param provided_variables:
#        :param initial_units:
#        :return:
#        """
#
#        self.climate_data = pd.DataFrame({})
#
#        config_reader = cp.ConfigParser()
#        config_reader.read(config_file)
#
#        missing_data_value = config_reader['METADATA']['missing_data_value']
#
#        raw_data = np.genfromtxt(climate_file, dtype='U', delimiter=',', skip_header=1, autostrip=True,
#                                 case_sensitive='lower')
#
#        raw_rows = raw_data.shape[0]  # number of rows in data
#        raw_cols = raw_data.shape[1]  # number of columns in data
#
#        # go through raw data and replace missing data values with nans
#        # note that values will not be nan until list is typecast as a float
#        for i in range(raw_rows):
#            for j in range(raw_cols):
#                if missing_data_value in raw_data[i, j]:
#                    raw_data[i, j] = np.nan
#                else:
#                    pass
#
#        # Iterate through all entries in the dictionary, check to see if they are set to -1, otherwise pull that column
#        # of data from the climate file and add it to the pandas dataframe
#        for variable_name, column in provided_variables.items():
#            if column == -1:  # Variable was not provided, fill empty column with nans
#                variable_observations = np.empty(raw_rows)
#                variable_observations[:] = np.nan
#            else:
#                variable_observations = np.array(raw_data[:, column].astype('float'))
#                variable_observations = self._convert_units(variable_name, variable_observations, initial_units)
#
#                self.climate_data[variable_name] = variable_observations
#
#        return self.climate_data

# -*- coding: utf-8 -*-
"""
Read and load in fluxnet or other climate time series files.
"""

from pathlib import Path

class Data(object):
"""
Read a climate station time series file using a config file.

Upon reading a config file for a specific climate station dataset the ``Data``
class stores metadata on climate variables in file, their units, and file path. 
Methods include tools to load climate time series data into a 
:obj:`pandas.DataFrame` and convert units to SI units. 
"""

    def __init__(self, config_file):

        self.config_file = Path(config_file).absolute()
        self.climate_file = self._get_climate_file()
        self.header_names = self._get_header()
        self.initial_units = self._get_units()


    def _get_climate_file(self):
        # read config file and get/return path to climate time series file
        pass

    def _get_header(self):
        # read top line of climate time series file return header names
        pass

    def _get_units(self):
        # read config file and get the names of variables and their units
        # return dict with names-keys, units-values
        pass

    def data_frame(self):
        """ 
        load climate time series data from file as a pandas.DataFrame
        first convert any units into SI
        """
        pass

    def _update_units(self, df):
        """using ``initial_units`` instance attribute, update any variables 
        units by converting them to SI units, used as helper in ``data_frame`` 
        method"""
        pass

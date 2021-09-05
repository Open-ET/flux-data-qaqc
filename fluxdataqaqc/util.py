# -*- coding: utf-8 -*-
"""
Collection of utility objects and functions for the :mod:`fluxdataqaqc`
module.
"""

import numpy as np
import pandas as pd
from pathlib import Path

class Convert(object):
    """
    Tools for unit conversions for ``flux-data-qaqc`` module.
    """
    # this is a work in progress, add more as needed/conversions are handled
    # input unit strings are not case sensitive, they will be forced to lower
    allowable_units = {
        'LE': ['w/m2','mj/m2'],
        'H': ['w/m2','mj/m2'],
        'Rn': ['w/m2','mj/m2'],
        'G': ['w/m2','mj/m2'],
        'lw_in': ['w/m2','mj/m2'],
        'lw_out': ['w/m2','mj/m2'],
        'sw_in': ['w/m2'], 
        'sw_out': ['w/m2','mj/m2'],
        'ppt': ['mm', 'in', 'm'],
        'vp': ['kpa', 'hpa', 'pa'],
        'vpd': ['kpa', 'hpa', 'pa'],
        't_avg': ['c', 'f', 'k'],
        't_min': ['c', 'f', 'k'],
        't_max': ['c', 'f', 'k'],
        'ws': ['m/s', 'mph']
    }

    # for printing and plotting purposes
    pretty_unit_names = {
        'pa': 'Pa',
        'hpa': 'hPa',
        'kpa': 'kPa',
        'c': 'C',
        'f': 'F',
        'k': 'K'
    }

    # some variables need to be in specified units for internal calculations
    # they will be attempted to be converted upon initialization of a QaQc obj
    # allowable initial units can be found in QaQc.allowable_units 
    required_units = {
        'LE': 'w/m2',
        'H': 'w/m2',
        'Rn': 'w/m2',
        'G': 'w/m2',
        'lw_in': 'w/m2',
        'lw_out': 'w/m2',
        'sw_in': 'w/m2',
        'sw_out': 'w/m2',
        'ppt': 'mm',
        'vp': 'kpa',
        'vpd': 'kpa',
        't_avg': 'c',
        't_min': 'c',
        't_max': 'c',
        'ws': 'm/s'
    }

    def __init__(self):

        self._conversion_map = {
            'k_to_c': self._k_to_c,
            'hpa_to_kpa': self._hpa_to_kpa,
            'pa_to_kpa': self._pa_to_kpa,
            'in_to_mm': self._in_to_mm,
            'm_to_mm': self._m_to_mm,
            'f_to_c': self._f_to_c,
            'mj/m2_to_w/m2': self._mj_per_m2_to_watts_per_m2,
            'mph_to_m/s': self._mph_to_m_per_s # miles/hr to meters/sec
        }

    @classmethod
    def convert(cls, var_name, initial_unit, desired_unit, df):
        """
        Givin a valid initial and desired variable dimension for a variable
        within a :obj:`pandas.DataFrame`, make the conversion and return the
        updated :obj:`pandas.DataFrame`.

        For a list of variables that require certain units within
        ``flux-data-qaqc`` see :attr:`Convert.allowable_units` (names of
        allowable options of input variable dimensions) and
        :attr:`Convert.required_units` (for the mandatory dimensions of certain
        variables before running QaQc calculations).

        Arguments:
            var_name (str): name of variable to convert in ``df``.
            initial_unit (str): name of initial unit of variable, must be valid
                from :attr:`Convert.allowable_units`.
            desired_unit (str): name of units to convert to, also must be valid.
            df (:obj:`pandas.DataFrame`): :obj:`pandas.DataFrame` containing
                variable to be converted, i.e. with ``var_name`` in columns.

        Returns:
            df (:obj:`pandas.DataFrame`): updated dataframe with specified variable's units converted

        Note:
            Many potential dimensions may not be provided for automatic
            conversion, if so you may need to update your variable dimensions
            manually, e.g. within a :attr:`.Data.df` before creating a
            :obj:`.QaQc` instance. Unit conversions are required for
            variables that can potentially be used in calculations within
            :obj:`.Data` or :obj:`.QaQc`.

        """
        conv = cls() 
        convert_key = '{}_to_{}'.format(initial_unit, desired_unit)
        convert_func = conv._conversion_map[convert_key]

        print(
            'Converting {} from {} to {}'.format(
                var_name, initial_unit, desired_unit
            )
        )
        df = convert_func(df, var_name)

        return df

    def _in_to_mm(self, df, var_name):
        df[var_name] *= 25.4
        return df
        
    def _m_to_mm(self, df, var_name):
        df[var_name] *= 1000
        return df
        
    def _f_to_c(self, df, var_name):
        df[var_name] = (32 * df[var_name]) * (5/9)
        return df
        
    def _k_to_c(self, df, var_name):
        df[var_name] -= 273.15
        return df
        
    def _hpa_to_kpa(self, df, var_name):
        df[var_name] /= 10
        return df
        
    def _pa_to_kpa(self, df, var_name):
        df[var_name] /= 1000
        return df
        
    def _mph_to_m_per_s(self, df, var_name):
        df[var_name] *= 0.44704 
        return df
        
    def _mj_per_m2_to_watts_per_m2(self, df, var_name):
        # assumes average mj per day is correct- only valid daily
        # because shortwate rad may be used in data (before daily) it is
        # not covered for automatic conversion because time period is unknown
        df[var_name] *= 11.574074074074074
        return df


def monthly_resample(df, cols, agg_str, thresh=0.75):
    """
    Resample dataframe to monthly frequency while excluding
    months missing more than a specified percentage of days of the month.

    Arguments:
        df (:obj:`pandas.DataFrame`): datetime indexed DataFrame instance
        cols (list): list of columns in `df` to resample to monthy frequency
        agg_str (str): resample function as string, e.g. 'mean' or 'sum'

    Keyword Arguments:
        thresh (float): threshold (decimal fraction) of how many days in a
            month must exist for it to be temporally resampled, otherwise
            the monthly value for the month will be null.

    Returns:
        ret (:obj:`pandas.DataFrame`): datetime indexed DataFrame that has been resampled to monthly time frequency.

    Note:
        If taking monthly totals (`agg_str` = 'sum') missing days will be filled
        with the months daily mean before summation.
    """
    if agg_str == 'sum':
        mdf = df.loc[:,cols].apply(pd.to_numeric).resample('M').agg(
            [agg_str, 'count', 'mean']
        )
    else:
        mdf = df.loc[:,cols].apply(pd.to_numeric).resample('M').agg(
            [agg_str, 'count']
        )
        
    ret = pd.DataFrame()
    for c in cols:
        bad_months = mdf.loc[:,(c,'count')] <= thresh * mdf.index.days_in_month
        if agg_str == 'sum':
            mdf.loc[:,(c,'days_missing')] =\
                mdf.index.days_in_month - mdf.loc[:,(c,'count')]
            ret[c] = mdf.loc[:,(c,agg_str)] +\
                (mdf.loc[:,(c,'days_missing')] * mdf.loc[:,(c,'mean')])
        else:
            ret[c] = mdf.loc[:,(c, agg_str)]
        ret.loc[bad_months, c] = np.nan

    return ret

def write_configs(meta_df, data_dict, out_dir=None):
    """
    Write multiple config files based on collection of site metadata and
    a dictionary containing variable information.

    Useful for creating config files for `flux-data-qaqc` for batches of
    flux stations that utilize the same naming conventions and formatting.

    Arguments:
        meta_df (:obj:`pandas.DataFrame`): dataframe that contains the
            following columns (or more) that describe metadata for multiple
            climate stations: 'site_id', 'climate_file_path',
            'station_longitude' 'station_elevation', 'station_latitude', and
            'missing_data_value'.  Elevation should be in meters and latitude
            is in decimal degrees.  Additional metadata columns will be added
            to the config file for each site, e.g. 'QC_flag',
            'anemometer_height', and any others.
        
        data_dict (dict): dictionary that maps `flux-data-qaqc` config
            names to user's column names in input files header e.g.
            {'net_radiation_col': 'netrad', 'net_radiation_units' : 'w/m2'}
            Anything that `flux-data-qaqc` config files "DATA" section
            can be present here including QC flag names, multiple soil
            moisture names and weights.

    Keyword Arguments:
        out_dir (str or None): default None. Directory to save config
            files, if None then save to currect working directory.

    Returns:
        configs (list): list of :obj:`pathlib.Path` objects of full paths
            to each config file written.

    Raises:
        Exception: if one of the mandatory metadata columns does not exist in
            `meta_df`.

    """
    mandatory_meta_vars = ['site_id', 'climate_file_path', 'station_longitude',
                'station_elevation', 'station_latitude', 'missing_data_value']

    if not set(mandatory_meta_vars).issubset(meta_df.columns):
        missing_cols = set(mandatory_meta_vars) - set(meta_df.columns)
        err_msg = ('ERROR: the following columns is missing from the '
            'metadata dataframe:\n{}'.format(
                ', '.join([v for v in missing_cols]))
        )
        raise Exception(err_msg)
    # if out_dir is None save to cwd
    if out_dir is None:
        out_dir = Path.cwd()
    else:
        out_dir = Path(out_dir)
    if not out_dir.is_dir():
        print(
            '{} does not exist, creating directory'.format(
                out_dir.absolute()
            )
        )
        out_dir.mkdir(parents=True, exist_ok=True)

    config_head = '#'*79+'\n'+'# this file was generated by flux-data-qaqc\n'+\
        '#'*79+'\n\n'+'[METADATA]\n'+'#'*10 +'\n'
    # save list of config file Path objects for later running qaqc
    configs = []
    cnt = 0
    for index, site in meta_df.iterrows():
        # make sure all mandatory data exists (not null) else skip
        if (site[mandatory_meta_vars]).isna().any():
            print('WARNING: site {} is missing mandatory metadata, skipping'\
                .format(site.site_id)
            )
            continue
        else:
            cnt += 1

        out_file = out_dir / '{}_config.ini'.format(site.site_id)
        configs.append(out_file.absolute())
        with open(out_file, 'w') as outf:
            outf.write(config_head)
            for meta_var in site.index:
                if not pd.isnull(site[meta_var]):
                    line = '{} = {}\n'.format(meta_var, site[meta_var])
                    outf.write(line)
            mid_lines = '\n[DATA]\n'+'#'*6+'\n'
            outf.write(mid_lines)
            for k,v in data_dict.items():
                line = '{} = {}\n'.format(k,v)
                outf.write(line)

    print(
        '\nSuccessfully wrote {} config files out of {} sites to folder:\n{}'\
        .format(cnt, len(meta_df), out_dir)
    )

    return configs

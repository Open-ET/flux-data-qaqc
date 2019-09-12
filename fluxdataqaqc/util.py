# -*- coding: utf-8 -*-
"""
Collection of utility objects and functions for the :mod:`fluxdataqaqc`
module.
"""

import numpy as np
import pandas as pd

class Convert(object):
    """
    Tools for unit conversions for ``flux-data-qaqc`` module.
    """
    # this is a work in progress, add more as needed/conversions are handled
    # input unit strings are not case sensitive, they will be forced to lower
    allowable_units = {
        'LE': ['w/m2'],
        'H': ['w/m2'],
        'Rn': ['w/m2'],
        'G': ['w/m2'],
        'lw_in': ['w/m2'],
        'lw_out': ['w/m2'],
        'sw_in': ['w/m2'],
        'sw_out': ['w/m2'],
        'ppt': ['mm', 'in'],
        'vp': ['kpa', 'hpa'],
        'vpd': ['kpa', 'hpa'],
        't_avg': ['c', 'f']
    }

    # for printing and plotting purposes
    pretty_unit_names = {
        'hpa': 'hPa',
        'kpa': 'kPa',
        'c': 'C',
        'f': 'F'
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
        't_avg': 'c'
    }

    def __init__(self):

        self._conversion_map = {
            'hpa_to_kpa': self._hpa_to_kpa,
            'in_to_mm': self._in_to_mm,
            'f_to_c': self._f_to_c
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
            conversion, if so you may need update your variable dimensions
            manually, e.g. within a :attr:`.Data.df` before creating a
            :obj:`.QaQc` instance. Unit conversions are required for
            variables that can potentially be used in calculations within
            :obj:`.QaQc`.

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
        df[var_name] = df[var_name] * 25.4
        return df

    def _f_to_c(self, df, var_name):
        df[var_name] = (32 * df[var_name]) * (5/9)
        return df

    def _hpa_to_kpa(self, df, var_name):
        df[var_name] /= 10
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

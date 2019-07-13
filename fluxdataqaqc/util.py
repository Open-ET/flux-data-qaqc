# -*- coding: utf-8 -*-
"""
Collection of utility objects or functions for use by the :mod:`fluxdataqaqc`
module.
"""

import numpy as np
import pandas as pd

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
        ret (:obj:`pandas.DataFrame`): datetime indexed DataFrame that
            has been resampled to monthly time frequency
    """
    mdf = df.loc[:,cols].apply(pd.to_numeric).resample('M').agg(
        [agg_str, 'count']
    )

    ret = pd.DataFrame()
    for c in cols:
        bad_months = mdf.loc[:,(c,'count')] <= thresh * mdf.index.days_in_month
        ret[c] = mdf.loc[:,(c, agg_str)]
        ret.loc[bad_months, c] = np.nan

    return ret

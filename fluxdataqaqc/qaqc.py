# -*- coding: utf-8 -*-
"""
Tools for correcting surface energy-balance components and calculating 
relevant variables such as Bowen's ratio, crop coefficient, energy balance
closure error fraction, and others. 

Input data is expected to be a :obj:`fluxdataqaqc.data.Data` instance.

TODO: 
 * reading in of data for comparing ETr:ETo (gridMET)
"""


import datetime as dt
from pathlib import Path
import numpy as np
import pandas as pd
from refet.calcs import _ra_daily, _rso_simple

from .data import Data


class QaQc(object):
    """
    Adjust latent energy and sensible heat flux to close the surface energy 
    balance, and other corrections for an eddy covariance climate station time 
    series. 
    
    Default instance creation uses a :obj:`fluxdataqaqc.Data` object but
    can also be initialized with a :obj:`pandas.DataFrame` that contains the
    required variables. 
    """
    agg_dict = {
        'energy': 'mean',
        'flux': 'mean',
        'flux_adj': 'mean',
        'flux_corr': 'mean',
        'bowen_ratio': 'mean',
        'et_reg': 'sum',
        'et_adj': 'sum',
        'et_corr': 'sum',
        'ebc_reg': 'mean',
        'ebc_adj': 'mean',
        'ebc_corr': 'mean',
        't_avg': 'mean',
        'rso': 'mean',
        'sw_pot': 'mean',
        'sw_in': 'mean',
        'lw_in': 'mean',
        'vpd': 'mean',
        'ppt': 'sum',
        'ws': 'mean',
        'Rn': 'mean',
        'sw_out': 'mean',
        'lw_out': 'mean',
        'G': 'mean',
        'LE': 'mean',
        'LE_corr': 'mean',
        'LE_adj': 'mean',
        'H': 'mean',
        'H_corr': 'mean',
        'H_adj': 'mean',
    }
 
    def __init__(self, data=None):
        
        if isinstance(data, Data):
            self._df = data.df
            self.elevation = data.elevation
            self.latitude = data.latitude
            self.out_dir = data.out_dir

        elif data is not None:
            print('{} is not a valid input type'.format(type(data)))
            raise TypeError("Must assign a fluxdataqaqc.data.Data object")
        else:
            self._df = None

        self.corrected = False 
        self._df = self._check_daily_freq()
            
    def _check_daily_freq(self):
        """check temporal frequency of input Data, resample to daily"""
        df = self._df

        if not isinstance(df, pd.DataFrame):
            return

        freq = pd.infer_freq(df.index)
        self.temporal_freq = freq

        if freq > 'D':
            print('WARNING: it looks like the input data temporal frequency',
                'is greater than daily, downsampling to daily, proceed with' ,
                'caution!')
        if freq < 'D':
            print('The input data temporal frequency appears to be less than',
                'daily, it will be resampled to daily.')

        if not freq == 'D':
            sum_cols = [k for k,v in QaQc.agg_dict.items() if v == 'sum']
            sum_cols = list(set(sum_cols).intersection(df.columns))
            mean_cols = set(df.columns) - set(sum_cols)
            means = df.loc[:,mean_cols].resample('D').mean()
            sums = df.loc[:,sum_cols].resample('D').sum()
            df = pd.concat([means, sums], sort=False)

            self._df = df

        return df

    @property     
    def df(self):
        # avoid overwriting pre-assigned data
        if isinstance(self._df, pd.DataFrame):
            return self._df

    @df.setter
    def df(self, data_frame):
        if not isinstance(data_frame, pd.DataFrame):
            raise TypeError("Must assign a Pandas.DataFrame object")
        self._df = data_frame

    @property
    def monthly_df(self):
        """
        Return current state of df as monthly time series.

        """
        if not self.corrected:
            self.correct_data()

        df = self._df

        sum_cols = [k for k,v in QaQc.agg_dict.items() if v == 'sum']
        # to avoid warning/error of missing columns
        sum_cols = list(set(sum_cols).intersection(df.columns))
        mean_cols = set(df.columns) - set(sum_cols)
        
        means = df.loc[:,mean_cols].resample('M').mean()
        sums = df.loc[:,sum_cols].resample('M').sum()
        df = pd.concat([means, sums], sort=False)
        # use monthly sums for ebc columns not means of ratio
        df.ebc_reg = (df.H + df.LE) / (df.Rn - df.G)
        df.ebc_adj = (df.H_adj + df.LE_adj) / (df.Rn - df.G)
        if set(['LE_corr','H_corr']).issubset(df.columns):
            df.ebc_corr= (df.H_corr + df.LE_corr) / (df.Rn - df.G)
        
        #elif how == 'aggregate':
        #    # for monthly stats not time series
        #    df = df.groupby(df.index.month).agg(self.agg_dict)
        #else:
        #    err_msg='Invalid "how" option, use "time_series" or "aggregate"'
        #    raise ValueError(err_msg)

        df.index.name = 'month'
        df.replace([np.inf, -np.inf], np.nan, inplace=True)

        return df

    def write(self, out_dir=None):
        """
        Save a copy of the "corrected" energy balance time series
        including raw input, save two CSVs one at daily and one at monthly
        time frequencies. 

        Arguments:
            out_dir (str or None): default None. Directory to save CSVs

        Returns:
            None

        Note:
            If this method is used before correcting the data according to the
            QA/QC routines in ``correct_data`` it will be done before saving. 
            Similarly, if monthly time series data has not yet been calculated, 
            it will be created here before saving. 
        """

        if out_dir is None:
            out_dir = self.out_dir
        else:
            out_dir = Path(out_dir)
            self.out_dir = out_dir.absolute()

        if not out_dir.is_dir():
            print(
                '{} does not exist, creating directory'.format(
                    out_dir.absolute()
                )
            )
            out_dir.mkdir(parents=True, exist_ok=True)

        if not self.corrected:
            self.correct_data

        daily_outf = out_dir / 'daily_data.csv'
        monthly_outf = out_dir / 'monthly_data.csv'

        self.df.to_csv(daily_outf)
        self.monthly_df.to_csv(monthly_outf)

    @classmethod
    def from_dataframe(cls, df, site_id, elev_m, lat_dec_deg):
        """
        Create a ``QaQc`` object from a pandas.DataFrame object.
        
        TODO: add checks for making sure mandatory columns are provided.
        """
        qaqc = cls()
        # use property setter, will load dataframe if needed
        qaqc.df = df  
        qaqc.latitude = lat_dec_deg
        qaqc.elevation = elev_m
        qaqc.out_dir = Path(site_id + '_output').absolute()
        return qaqc
    
    def correct_data(self):
        """
        Create corrected/adjusted latent energy and sensible heat flux to 
        close surface energy balance. 
        
        Updates :attr:`QaQc.df` attribute with new variables used for correcting
        energy fluxes.
        
        """
        if not isinstance(self._df, pd.DataFrame):
            print('Please assign a dataframe of acceptable data first!')
            return
        
        # get length of data set
        data_length = len(self.df.index)

        self.df['energy'] = self.df.Rn - self.df.G
        self.df['flux'] = self.df.LE + self.df.H
        self.df['bowen_ratio'] = self.df.H / self.df.LE

        # numpy arrays of dataframe vars
        Rn = np.array(self.df.Rn)
        g = np.array(self.df.G)
        le = np.array(self.df.LE)
        h = np.array(self.df.H)
        bowen = np.array(self.df.bowen_ratio)
        # numpy arrays of new vars
        le_adj = np.full(data_length, np.NaN)
        h_adj = np.full(data_length, np.NaN)
        flux_adj = np.full(data_length, np.NaN)

        # compute adjusted turbulent fluxes for when Rn > 0
        # correcting LE and H, method may be faster as function and vectorized
        for i in range(0, data_length):
            if Rn[i] > 0:
                le_adj[i] = (Rn[i] - g[i]) / (1 + bowen[i])
                h_adj[i] = (bowen[i] / (1 + bowen[i])) * (Rn[i] - g[i])

            else:
                le_adj[i] = le[i]
                h_adj[i] = h[i]

        # compute adjusted turbulent fluxes for when Rn > 0 and Bowen ratio < 0.05
        # instead of forcing closure when bowen ratio is often <- 0.8 or threshold when 
        # Rn < 0. The average between i-1 and i+1 is taken. If closure is forced by 
        # partitioning the error equally between LE and H, LE is drastically increased
        # during periods of possibly "bad" data, usually while the measured LE is going down

        for i in range(0, data_length):
            if Rn[i] > 0 and bowen[i] < 0.05:
                le_adj[i] = ((le[i - 1]) + (le[i + 1]))/2
                h_adj[i] = ((h[i - 1]) + (h[i + 1]))/2

            flux_adj[i] = le_adj[i] + h_adj[i]

           # # If adjusted fluxes are less than original fluxes, keep originals
           # no physical reason to include this
           # if le_adj[i] < le[i]:
           #     le_adj[i] = le[i]

           # if h_adj[i] < h[i]:
           #     h_adj[i] = h[i]

        # add le_adj, h_adj, and flux_adj to dataframe
        # TODO: flux_adj in blake's code is placed to not reflect final h and le adj values, confirm with him
        self.df['LE_adj'] = le_adj
        self.df['H_adj'] = h_adj
        self.df['flux_adj'] = flux_adj

        # corrected turbulent flux if given from input data
        if set(['LE_corr','H_corr']).issubset(self.df.columns):
            self.df['flux_corr'] = self.df.LE_corr + self.df.H_corr 
            self.df['et_corr'] = 86400 * (self.df.LE_corr /(2500000 * 1000)) * 1000
            self.df['ebc_corr'] = (self.df.H_corr + self.df.LE_corr) / (self.df.Rn - self.df.G)
            self.df.ebc_corr = self.df.ebc_corr.replace([np.inf, -np.inf], np.nan)

        # add ET/EBC columns to dataframe using le and h of various sources
        #  _reg uses uncorrected le and h
        #  _adj uses ajusted le and h from our bowen ratio corrections
        #  _corr uses corrected le and h as found in data file (if provided)
        self.df['et_reg'] = 86400 * (self.df.LE/(2500000 * 1000)) * 1000
        self.df['et_adj'] = 86400 * (self.df.LE_adj/(2500000 * 1000)) * 1000

        self.df['ebc_reg'] = (self.df.H + self.df.LE) / (self.df.Rn - self.df.G)
        self.df['ebc_adj'] = (self.df.H_adj + self.df.LE_adj) / (self.df.Rn - self.df.G)

        # replace undefined/infinity with nans in all EBC columns
        self.df.ebc_reg = self.df.ebc_reg.replace([np.inf, -np.inf], np.nan)
        self.df.ebc_adj = self.df.ebc_adj.replace([np.inf, -np.inf], np.nan)

        # create date vectors for obtaining day of year for use in calculating ra and then rso
        date = pd.DatetimeIndex(self.df.index)
        day = np.array(date.day)
        month = np.array(date.month)
        year = np.array(date.year)

        # Calculate DOY from Y/M/D values
        doy = []
        for i in range(data_length):
            doy.append(
                dt.date(year[i], month[i], day[i]).strftime("%j"))  # list of string DOY values

        doy = np.array(list(map(int, doy)))  # Converts list of string values into ints and saves as numpy array

        # obtain extraterrestrial radiation from doy and latitude and then calculate
        # clear sky radiation  (simple version based on elevation)
        latitude_rads = self.latitude * (np.pi / 180)
        ra_mj_m2 = _ra_daily(latitude_rads, doy, method='asce')
        rso_a_mj_m2 = _rso_simple(ra_mj_m2, self.elevation)
        self.df['rso'] = rso_a_mj_m2 * 11.574

        # update flag for other methods
        self.corrected = True

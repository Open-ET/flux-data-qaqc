# -*- coding: utf-8 -*-
"""
Tools for correcting surface energy-balance components and calculating 
relevant variables such as Bowen's ratio, crop coefficient, energy balance
closure error fraction, and others. 

Input data is currently expected as a daily time series that has been 
converted into a pandas.DataFrame from an eddy covariance climate station.

TODO: 
 * add calcs, reading in of data for comparing ETr:ETo (gridMET)
 * add monthly stat calcs, energy balance closure error, energy-balance ratio
"""

import numpy as np
import pandas as pd

from data import Data

class QaQc(object):
    """
    Adjust latent energy and sensible heat flux to close the surface energy 
    balance, and other corrections for an eddy covariance climate station time 
    series. 
    
    Default instance creation uses a :obj:`fluxdataqaqc.Data` object but
    can also be initialized with a :obj:`pandas.DataFrame` that contains the
    required variables. 
    """
    
    def __init__(self, data=None):
        
        if isinstance(data, Data):
            self._df = data.df
            self._elevation = data.elevation
        elif data is not None:
            raise TypeError("Must assign a fluxdataqaqc.data.Data object")
        else:
            self._df = None

        self.corrected = False 
            
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

        # TODO: maybe improve handling of units/var names
        """
        if not self.corrected:
            self.correct_data()

        df = self._df

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
            'net_rad': 'mean',
            'sw_out': 'mean',
            'lw_out': 'mean',
            'g_flux': 'mean',
            'le_flux': 'mean',
            'le_flux_corr': 'mean',
            'le_flux_adj': 'mean',
            'h_flux': 'mean',
            'h_flux_corr': 'mean',
            'h_flux_adj': 'mean',
        }

        sum_cols = ['et_reg', 'et_adj', 'et_corr', 'ppt']
        mean_cols = set(df.columns) - set(sum_cols)
        
        means = df.loc[:,mean_cols].resample('M').mean()
        sums = df.loc[:,sum_cols].resample('M').sum()
        df = pd.concat([means, sums], sort=False)

        # for monthly stats not time series
        #df = df.groupby(df.index.month).agg(agg_dict)
        
        df.index.name = 'month'

        return df

    @classmethod
    def from_dataframe(cls, df):
        """
        Create a ``QaQc`` object from a pandas.DataFrame object.
        
        TODO: add checks for making sure mandatory columns are provided.
        """
        qaqc = cls()
        # use property setter, will load dataframe if needed
        qaqc.df = df  
        
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

        # create energy, flux, and bowen ratio columns within dataframe
        self.df['energy'] = self.df.net_rad - self.df.g_flux
        self.df['flux'] = self.df.le_flux + self.df.h_flux
        self.df['bowen_ratio'] = self.df.h_flux / self.df.le_flux

        # numpy arrays of dataframe vars
        net_rad = np.array(self.df.net_rad)
        g_flux = np.array(self.df.g_flux)
        le_flux = np.array(self.df.le_flux)
        h_flux = np.array(self.df.h_flux)
        bowen = np.array(self.df.bowen_ratio)
        # numpy arrays of new vars
        le_flux_adj = np.full(data_length, np.NaN)
        h_flux_adj = np.full(data_length, np.NaN)
        flux_adj = np.full(data_length, np.NaN)

        # compute adjusted turbulent fluxes for when Rn > 0
        # correcting LE and H, method may be faster as function and vectorized
        for i in range(0, data_length):
            if net_rad[i] > 0:
                le_flux_adj[i] = (net_rad[i] - g_flux[i]) / (1 + bowen[i])
                h_flux_adj[i] = (bowen[i] / (1 + bowen[i])) * (net_rad[i] - g_flux[i])

            else:
                le_flux_adj[i] = le_flux[i]
                h_flux_adj[i] = h_flux[i]

        # compute adjusted turbulent fluxes for when Rn > 0 and Bowen ratio < 0.05
        # instead of forcing closure when bowen ratio is often <- 0.8 or threshold when 
        # Rn < 0. The average between i-1 and i+1 is taken. If closure is forced by 
        # partitioning the error equally between LE and H, LE is drastically increased
        # during periods of possibly "bad" data, usually while the measured LE is going down

        for i in range(0, data_length):
            if net_rad[i] > 0 and bowen[i] < 0.05:
                le_flux_adj[i] = ((le_flux[i - 1]) + (le_flux[i + 1]))/2
                h_flux_adj[i] = ((h_flux[i - 1]) + (h_flux[i + 1]))/2

            flux_adj[i] = le_flux_adj[i] + h_flux_adj[i]

            # If adjusted fluxes are less than original fluxes, keep originals
            if le_flux_adj[i] < le_flux[i]:
                le_flux_adj[i] = le_flux[i]

            if h_flux_adj[i] < h_flux[i]:
                h_flux_adj[i] = h_flux[i]


        # add le_flux_adj, h_flux_adj, and flux_adj to dataframe
        # TODO: flux_adj in blake's code is placed to not reflect final h_flux and le_flux adj values, confirm with him
        self.df['le_flux_adj'] = le_flux_adj
        self.df['h_flux_adj'] = h_flux_adj
        self.df['flux_adj'] = flux_adj

        # corrected turbulent flux if given from input data
        self.df['flux_corr'] = self.df.le_flux_corr + self.df.h_flux_corr 

        # add ET/EBC columns to dataframe using le_flux and h_flux of various sources
        #  _reg uses uncorrected le_flux and h_flux
        #  _adj uses ajusted le_flux and h_flux from our bowen ratio corrections
        #  _corr uses corrected le_flux and h_flux as found in data file (if provided)
        self.df['et_reg'] = 86400 * (self.df.le_flux/(2500000 * 1000)) * 1000
        self.df['et_adj'] = 86400 * (self.df.le_flux_adj/(2500000 * 1000)) * 1000
        self.df['et_corr'] = 86400 * (self.df.le_flux_corr /(2500000 * 1000)) * 1000

        self.df['ebc_reg'] = (self.df.h_flux + self.df.le_flux) / (self.df.net_rad - self.df.g_flux)
        self.df['ebc_adj'] = (self.df.h_flux_adj + self.df.le_flux_adj) / (self.df.net_rad - self.df.g_flux)
        self.df['ebc_corr'] = (self.df.h_flux_corr + self.df.le_flux_corr) / (self.df.net_rad - self.df.g_flux)

        # replace undefined/infinity with nans in all EBC columns
        self.df.ebc_reg = self.df.ebc_reg.replace([np.inf, -np.inf], np.nan)
        self.df.ebc_adj = self.df.ebc_adj.replace([np.inf, -np.inf], np.nan)
        self.df.ebc_corr = self.df.ebc_corr.replace([np.inf, -np.inf], np.nan)

        # clear sky radiation calc (simple version based on elevation)
        ra_mj_m2 = np.array(self.df.sw_in * 0.0864)  # http://www.fao.org/3/X0490E/x0490e0i.htm
        rso_a_mj_m2 = np.array((0.75 + 2E-5 * self._elevation) * ra_mj_m2)  # asce 19 and 45
        self.df['rso'] = rso_a_mj_m2 * 11.574

        # update flag for other methods
        self.corrected = True

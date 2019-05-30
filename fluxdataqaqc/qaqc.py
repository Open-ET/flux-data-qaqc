# -*- coding: utf-8 -*-
"""
Tools for correcting surface energy-balance components and calculating 
relevant variables such as energy balance closure ratio, corrected latent
energy and sensible heat fluxes, potential clear sky radiation, ET, and others. 
TODO: 
 * reading in of data for comparing ETr:ETo (gridMET)
"""


from pathlib import Path
import numpy as np
import pandas as pd
from refet.calcs import _ra_daily, _rso_simple

from .data import Data

class QaQc(object):
    """
    Adjust daily latent energy and sensible heat fluxes to improve closure of 
    the surface energy balance. Includes other qa/qc calculations for eddy 
    covariance climate station time series. 
    
    Input data is expected to be a :obj:`fluxdataqaqc.data.Data` instance or a
    :obj:`pandas.DataFrame` which can be used to create a :obj:`QaQc` object 
    with the :meth:`QaQc.from_dataframe` method.
    """
    agg_dict = {
        'energy': 'mean',
        'flux': 'mean',
        'flux_corr': 'mean',
        'bowen_ratio': 'mean',
        'et_reg': 'sum',
        'et_corr': 'sum',
        'et_user_corr': 'sum',
        'ebc_reg': 'mean',
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
        'LE_user_corr': 'mean',
        'H': 'mean',
        'H_corr': 'mean',
        'H_user_corr': 'mean',
    }
 
    corr_methods = (
        'fluxnet',
        'bowen_ratio'
    )

    def __init__(self, data=None):
        
        if isinstance(data, Data):
            self.variables = data.variables
            self._df = data.df
            self.elevation = data.elevation
            self.latitude = data.latitude
            self.out_dir = data.out_dir
            self.site_id = data.site_id
            # flip variable naming dict for internal use
            self.inv_map = {v: k for k, v in self.variables.items()}
            # using 'G' in multiple g plot may overwrite G name internally
            if not 'G' in self.inv_map.values():
                user_G_name = self.variables.get('G') 
                self.inv_map[user_G_name] = 'G'

            self.temporal_freq = self._check_daily_freq()
            

        elif data is not None:
            print('{} is not a valid input type'.format(type(data)))
            raise TypeError("Must assign a fluxdataqaqc.data.Data object")
        else:
            self._df = None

        self.corrected = False 
        self.corr_meth = None
            
    def _check_daily_freq(self):
        """
        Check temporal frequency of input Data, resample to daily if not already

        Note:
            If user QC values for filtering data are present they will be 
            resampled to daily means, however this should not be an issue as 
            the filtering step occurs in a :obj:`fluxdataqaqc.Data` object.
        
        """

        # rename columns to internal names 
        df = self._df.rename(columns=self.inv_map)
       
        if not isinstance(df, pd.DataFrame):
            return

        freq = pd.infer_freq(df.index)

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
            df = means.join(sums)

        # rename columns back to user's
        self._df = df.rename(columns=self.variables)
        return freq

    @property     
    def df(self):
        # avoid overwriting pre-assigned data
        if isinstance(self._df, pd.DataFrame):
            return self._df.rename(columns=self.variables)

    @df.setter
    def df(self, data_frame):
        if not isinstance(data_frame, pd.DataFrame):
            raise TypeError("Must assign a Pandas.DataFrame object")
        self._df = data_frame

    @property
    def monthly_df(self):
        """
        Return current state of df as monthly time series.

        If energy balance data has not yet been corrected, correct using
        the FLUXNET method. 

        Arguements:
            None

        Returns:
            None

        """
        if not self.corrected:
            self.correct_data()

        # rename columns to internal names 
        df = self._df.rename(columns=self.inv_map).copy()

        sum_cols = [k for k,v in QaQc.agg_dict.items() if v == 'sum']
        # to avoid warning/error of missing columns 
        sum_cols = list(set(sum_cols).intersection(df.columns))
        mean_cols = set(df.columns) - set(sum_cols)
        # if data type has changed to 'obj' resample skips... 
        means = df.loc[:,mean_cols].astype(float).resample('M').mean()
        sums = df.loc[:,sum_cols].astype(float).resample('M').sum()
        df = means.join(sums)
        # use monthly sums for ebc columns not means of ratio
        df.ebc_reg = (df.H + df.LE) / (df.Rn - df.G)
        df.ebc_corr = (df.H_corr + df.LE_corr) / (df.Rn - df.G)
        if set(['LE_user_corr','H_user_corr']).issubset(df.columns):
            df['ebc_user_corr']=(df.H_user_corr+df.LE_user_corr) / (df.Rn-df.G)
        
        #elif how == 'aggregate':
        #    # for monthly stats not time series
        #    df = df.groupby(df.index.month).agg(self.agg_dict)
        #else:
        #    err_msg='Invalid "how" option, use "time_series" or "aggregate"'
        #    raise ValueError(err_msg)

        df.index.name = 'month'
        df.replace([np.inf, -np.inf], np.nan, inplace=True)

        return df.rename(columns=self.variables)

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
            self.correct_data()

        daily_outf = out_dir / '{}_daily_data.csv'.format(self.site_id)
        monthly_outf = out_dir / '{}_monthly_data.csv'.format(self.site_id)

        self.df.to_csv(daily_outf)
        self.monthly_df.to_csv(monthly_outf)

    @classmethod
    def from_dataframe(cls, df, site_id, elev_m, lat_dec_deg, var_dict):
        """
        Create a ``QaQc`` object from a pandas.DataFrame object.
        
        Arguments:
            df (pandas.DataFrame): DataFrame of climate variables with
                datetime index named 'date'
            site_id (str): site identifier
            elev_m (int or float): elevation of site in meters
            lat_dec_deg (float): latitude of site in decimal degrees
            var_dict (dict): dictionary that maps `flux-data-qaqc` variable
                names to user's columns in `df` e.g. {'Rn': 'netrad', ...}
                see :attr:`fluxdataqaqc.Data.variable_names_dict` for list of 
                `flux-data-qaqc` variable names
        
        Returns:
            None

        """
        qaqc = cls()
        # use property setter, make sure it is a dataframe object
        qaqc.df = df  
        qaqc.site_id = site_id
        qaqc.latitude = lat_dec_deg
        qaqc.elevation = elev_m
        qaqc.out_dir = Path('output').absolute()
        qaqc.variables = var_dict
        # TODO handle assigned units 
        qaqc.inv_map = {v: k for k, v in var_dict.items()}
        qaqc.temporal_freq = qaqc._check_daily_freq()
        qaqc.corrected = False 
        qaqc.corr_meth = None

        return qaqc

    def correct_data(self, meth='fluxnet'):
        """
        Correct turblent fluxes to close energy balance using different
        methods, default 'fluxnet'. 

        Currently two options are available: 'fluxnet' and 'bowen_ratio'. If
        you use one method followed by another corrected versions of LE, H, 
        ET, and EBC will be overwritten with the most recently used approach.
        Also computes potential clear sky radiation (saved as *rso*) using
        a simple approach based on station elevation and latitude.

        Corrected or otherwise newly calculated variables are named using the
        following suffixes to distinguish them::

          _reg uses uncorrected LE and H from input data
          _corr uses adjusted LE and H from the correction method used 
          _user_corr uses corrected LE and H found in data file (if provided)

        Arguments:
            meth (str): default 'fluxnet'. Method to correct energy balance.
        
        Returns
            None

        Note:
            The *ebc_corr* variable or energy balance closure ratio is 
            calculated from the corrected versions of LE and H independent 
            of the method. When using the 'fluxnet' method the energy balance 
            correction factor (what is applied to the raw H and LE) is left as 
            calculated (inverse of ebc) and saved as *ebc_cf*. 
        """

        # in case starting in Python and no df assigned yet
        if not isinstance(self._df, pd.DataFrame):
            print('Please assign a dataframe of acceptable data first!')
            return
        if meth not in self.corr_methods:
            err_msg = ('ERROR: {} is not a valid correction option, please'
                'use one of the following: {}'.format(meth, ','.join(
                    [el for el in self.corr_methods]))
            )
            raise ValueError(err_msg)

        # calculate clear sky radiation if not already computed
        self._calc_rso()
        # energy balance corrections
        if meth == 'fluxnet':
            self._fluxnet_correction()
        if meth == 'bowen_ratio':
            self._bowen_ratio_correction()
        # store method that current data was corrected by 
        self.corr_meth = meth
        # update inv map for naming
        self.inv_map = {v: k for k, v in self.variables.items()}
        # using 'G' in multiple g plot may overwrite G name internally
        if not 'G' in self.inv_map.values():
            user_G_name = self.variables.get('G')
            self.inv_map[user_G_name] = 'G'

    def _calc_rso(self):
        """
        Calculate clear sky potential solar radiation using station latitude,
        elevation, and day of year using simple method from :mod:`refet`.

        Arguments:
            None

        Returns:
            None

        Note:
            Does not overwrite existing calculation for rso, i.e. once any 
            correction method has been run using :func:`correct_data` it will
            not be recalculated after subsequent calls since it is indepent of
            climate variables.
        """
        # avoid recalculating 
        if 'rso' in self.variables:
            return

        doy = self.df.index.dayofyear
        # obtain extraterrestrial radiation from doy and latitude and calculate
        # clear sky radiation 
        latitude_rads = self.latitude * (np.pi / 180)
        ra_mj_m2 = _ra_daily(latitude_rads, doy, method='asce')
        rso_a_mj_m2 = _rso_simple(ra_mj_m2, self.elevation)
        self._df['rso'] = rso_a_mj_m2 * 11.574
        
        self.variables.update(
            rso = 'rso'
        )

    def _fluxnet_correction(self):
        """
        FLUXNET correction method for daily LE, H, EBC, and ET.

        Correct turblent fluxes to close energy balance using methods 
        described by FLUXNET for daily H and LE. 

        Updates :attr:`QaQc.df` and :attr:`QaQc.variables` attributes with new 
        variables related to the corrections, e.g. LE_corr, ebc_corr, etc.

        Arguments:
            None

        Returns:
            None

        """
        # moving windows FLUXNET methods 1, 2 and 3
        window_1 = 15
        window_2 = 11
        half_win_1 = window_1 // 2
        half_win_2 = window_2 // 2
        
        # make sure names of variables are internal
        df = self._df.rename(columns=self.inv_map)
        # energy balance ratio (normal) from raw data
        df['ebc_reg'] = (df.H + df.LE) / (df.Rn - df.G)
        df = df[['ebc_reg', 'H', 'LE']].copy() # only data we need
        # compute IQR to filter out extreme EBC_CFs, 
        # this may need to be changed to the last step 
        Q1 = df['ebc_reg'].quantile(0.25)
        Q3 = df['ebc_reg'].quantile(0.75)
        IQR = Q3 - Q1
        # filter values between Q1-1.5IQR and Q3+1.5IQR
        filtered = df.query(
            '(@Q1 - 1.5 * @IQR) <= ebc_reg <= (@Q3 + 1.5 * @IQR)'
        )
        # apply filter
        filtered_mask = filtered.index
        removed_mask = set(df.index) - set(filtered_mask)
        removed_mask = pd.to_datetime(list(removed_mask))
        df.loc[removed_mask] = np.nan

        # make 5 day climatology of EBC_CF for method 3
        # at this stage ebc_cf and ebc_cf 5 day climatology are inverted (ebc)
        # the cenetered window skips first and last 5 DOYs
        # so prepend and append first and last 5 days and loop...
        doy_EBC_CF_mean=df['ebc_reg'].groupby(df.index.dayofyear).mean().copy()
        l5days = pd.Series(
            index=np.arange(-4,1), data=doy_EBC_CF_mean[-5:].values)
        f5days = pd.Series(
            index=np.arange(367,372), data=doy_EBC_CF_mean[:5].values)
        doy_EBC_CF_mean = doy_EBC_CF_mean.append(f5days)
        doy_EBC_CF_mean = pd.concat([l5days, doy_EBC_CF_mean])
        ebc_cf_5day_clim = pd.DataFrame(
            index=np.arange(1,367), columns=['ebc_cf_5day_clim'])
        doy_EBC_CF_mean = doy_EBC_CF_mean.values
        for i in range(len(doy_EBC_CF_mean)):
            # i = 0 which starts at prepended 5 days, shift window up
            win = doy_EBC_CF_mean[i:i+2*half_win_2+1]
            count = np.count_nonzero(~np.isnan(win))
            # get 11 day moving window mean
            if i in ebc_cf_5day_clim.index and count > 0:
                ebc_cf_5day_clim.iloc[
                    i-1, ebc_cf_5day_clim.columns.get_loc('ebc_cf_5day_clim')
                ] = np.nanmean(win)
        ebc_cf_5day_clim['DOY'] = ebc_cf_5day_clim.index
        ebc_cf_5day_clim.index.name = 'date'

        # methods 1 and 2
        EBC_CF = df.ebc_reg.values
        df['ebc_corr'] = np.nan
        # gap filling EBC_CF following methods 1 and 2
        for i in range(len(EBC_CF)):
            win_arr1 = EBC_CF[i-half_win_1:i+half_win_1+1]
            win_arr2 = EBC_CF[i-half_win_2:i+half_win_2+1]
            count = np.count_nonzero(~np.isnan(win_arr1))
            # get median of daily window1 if half window2 or more days exist
            if count >= half_win_2:
                val = np.nanpercentile(win_arr1, 50, axis=None)
            # if at least one day exists in window2 take mean
            elif np.count_nonzero(~np.isnan(win_arr2)) > 0:
                val = np.nanmean(win_arr2)
            else:
                # assign nan for now, update with 5 day climatology
                val = np.nan
            # assign values if they were found in methods 1 or 2
            df.iloc[i, df.columns.get_loc('ebc_corr')] = val
            
        # method 3 to fill remaining gaps 
        df['DOY'] = df.index.dayofyear
        # datetime indices of all remaining null elements
        null_dates = df.loc[df.ebc_corr.isnull(), 'ebc_corr'].index
        merged = pd.merge(
            df, ebc_cf_5day_clim, on='DOY', how='left', right_index=True
        )
        # assign 5 day climatology of EBC 
        merged.loc[null_dates,'ebc_corr'] =\
            merged.loc[null_dates,'ebc_cf_5day_clim']
        
        # apply corrections to LE and H multiply by 1/EBC
        merged['LE_corr'] = merged.LE * (1/merged.ebc_corr)
        merged['H_corr'] = merged.H * (1/merged.ebc_corr)
        # compute FLUXNET corrected ET (mm/d), total turb flux
        merged['flux_corr'] = merged['LE_corr'] + merged['H_corr']
        merged['et_corr'] = 86400 * (merged.LE_corr /(2500000 * 1000)) * 1000

        df = self._df.rename(columns=self.inv_map)
        # other variables needed for plots using raw data
        df['et_reg'] = 86400 * (df.LE/(2500000 * 1000)) * 1000
        df['ebc_reg'] = (df.H + df.LE) / (df.Rn - df.G)
        df['energy'] = df.Rn - df.G
        df['flux'] = df.LE + df.H

        # corrected turbulent flux if given from input data
        if set(['LE_user_corr','H_user_corr']).issubset(df.columns):
            df['flux_user_corr'] = df.LE_user_corr + df.H_user_corr 
            df['et_user_corr']=86400*(df.LE_user_corr /(2500000 * 1000)) * 1000
            df['ebc_user_corr']=(df.H_user_corr+df.LE_user_corr)/(df.Rn - df.G)
            df.ebc_user_corr=df.ebc_user_corr.replace([np.inf,-np.inf], np.nan)
            self.variables.update(
                flux_user_corr = 'flux_user_corr',
                et_user_corr = 'et_user_corr',
                ebc_user_corr = 'ebc_user_corr'
            )
        # grab select columns to merge into main dataframe
        cols = list(set(merged.columns).difference(df.columns))
        # join calculated data in
        merged = df.join(merged[cols], how='outer')
        # calculated corrected EBC to assign to ebc_corr (not EBC_CF), save
        # CFs as defined by fluxnet method, i.e. inverse of EBC
        merged['ebc_cf'] = 1/merged.ebc_corr
        merged['ebc_cf_5day_clim'] = 1/merged.ebc_cf_5day_clim
        merged['ebc_corr'] =\
            (merged.H_corr + merged.LE_corr) / (merged.Rn - merged.G)
        merged.drop('DOY', axis=1, inplace=True)

        self.variables.update(
            energy = 'energy',
            flux = 'flux',
            LE_corr = 'LE_corr',
            H_corr = 'H_corr',
            flux_corr = 'flux_corr',
            et_reg = 'et_reg',
            et_corr = 'et_corr',
            ebc_reg = 'ebc_reg',
            ebc_corr = 'ebc_corr',
            ebc_cf = 'ebc_cf',
            ebc_cf_5day_clim = 'ebc_cf_5day_clim'
        )

        # revert column names to user's
        self._df = merged.rename(columns=self.variables)
        # update flag for other methods
        self.corrected = True
    
    def _bowen_ratio_correction(self):
        """
        Create corrected/adjusted latent energy and sensible heat flux to 
        close surface energy balance. 
        
        Compute adjusted turbulent fluxes for when Rn > 0 and Bowen ratio 
        < 0.05 instead of forcing closure when bowen ratio is often <- 0.8 or 
        threshold when Rn < 0. The average between i-1 and i+1 is taken. If 
        closure is forced by partitioning the error equally between LE and H, 
        LE is drastically increased during periods of possibly "bad" data, 
        usually while the measured LE is going down.

        Updates :attr:`QaQc.df` and :attr:`QaQc.variables` attributes with new 
        variables used for closing energy balance.
        
        """
        
        # get length of data set
        data_length = len(self.df.index)

        # rename columns to internal names 
        df = self._df.rename(columns=self.inv_map)

        df['energy'] = df.Rn - df.G
        df['flux'] = df.LE + df.H
        df['bowen_ratio'] = df.H / df.LE

        # numpy arrays of dataframe vars
        Rn = df.Rn.values
        g = df.G.values
        le = df.LE.values
        h = df.H.values
        bowen = df.bowen_ratio.values
        # numpy arrays of new vars
        le_corr = np.full(data_length, np.NaN)
        h_corr = np.full(data_length, np.NaN)
        flux_corr = np.full(data_length, np.NaN)

        # compute adjusted turbulent fluxes for when Rn > 0
        # correcting LE and H, method may be faster as function and vectorized
        for i in range(0, data_length):
            if Rn[i] > 0:
                le_corr[i] = (Rn[i] - g[i]) / (1 + bowen[i])
                h_corr[i] = (bowen[i] / (1 + bowen[i])) * (Rn[i] - g[i])

            else:
                le_corr[i] = le[i]
                h_corr[i] = h[i]

        for i in range(0, data_length):
            if Rn[i] > 0 and bowen[i] < 0.05:
                le_corr[i] = ((le[i - 1]) + (le[i + 1]))/2
                h_corr[i] = ((h[i - 1]) + (h[i + 1]))/2

            flux_corr[i] = le_corr[i] + h_corr[i]

           # # If adjusted fluxes are less than original fluxes, keep originals
           # no physical reason to include this
           # if le_corr[i] < le[i]:
           #     le_corr[i] = le[i]

           # if h_corr[i] < h[i]:
           #     h_corr[i] = h[i]

        # add le_corr, h_corr, and flux_corr to dataframe
        df['LE_corr'] = le_corr
        df['H_corr'] = h_corr
        df['flux_corr'] = flux_corr

        # corrected turbulent flux if given from input data
        if set(['LE_user_corr','H_user_corr']).issubset(df.columns):
            df['flux_user_corr'] = df.LE_user_corr + df.H_user_corr 
            df['et_user_corr']=86400*(df.LE_user_corr /(2500000 * 1000)) * 1000
            df['ebc_user_corr']=(df.H_user_corr+df.LE_user_corr)/(df.Rn - df.G)
            df.ebc_user_corr=df.ebc_user_corr.replace([np.inf,-np.inf], np.nan)

            self.variables.update(
                flux_user_corr = 'flux_user_corr',
                et_user_corr = 'et_user_corr',
                ebc_user_corr = 'ebc_user_corr'
            )

        # add ET/EBC columns to dataframe using LE and H from raw, corr, and 
        # if provided user corrected
        df['et_reg'] = 86400 * (df.LE/(2500000 * 1000)) * 1000
        df['et_corr'] = 86400 * (df.LE_corr/(2500000 * 1000)) * 1000
        df['ebc_reg'] = (df.H + df.LE) / (df.Rn - df.G)
        df['ebc_corr'] = (df.H_corr + df.LE_corr) / (df.Rn - df.G)

        self.variables.update(
            bowen_ratio = 'bowen_ratio',
            energy = 'energy',
            flux = 'flux',
            LE_corr = 'LE_corr',
            H_corr = 'H_corr',
            flux_corr = 'flux_corr',
            et_reg = 'et_reg',
            et_corr = 'et_corr',
            ebc_reg = 'ebc_reg',
            ebc_corr = 'ebc_corr'
        )
        # replace undefined/infinity with nans in all EBC columns
        df.ebc_reg = df.ebc_reg.replace([np.inf, -np.inf], np.nan)
        df.ebc_corr = df.ebc_corr.replace([np.inf, -np.inf], np.nan)

        # revert column names to user's
        self._df = df.rename(columns=self.variables)
        # update flag for other methods
        self.corrected = True

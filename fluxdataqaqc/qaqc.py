# -*- coding: utf-8 -*-
"""
Includes routines for 'correcting' turbulent surface energy-balance components 
to improve energy balance closure. 

The default routine follows the procedure documented by `FLUXNET <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`_ (see section 3. Heat Processing, Daily data). Other tools include estimation of ASCE clear sky 
radiation, evapotranspiration, and other statistical variables. Input data can be either a :obj:`fluxdataqaqc.Data` instance or a :obj:`pandas.DataFrame`. 

TODO: 
 * reading in of data for comparing ETr:ETo (gridMET)
"""


from pathlib import Path

import xarray
import numpy as np
import pandas as pd
from refet.calcs import _ra_daily, _rso_simple

from .data import Data
from .plot import Plot
from .util import monthly_resample

class QaQc(Plot):
    """
    Numerical routines for adjusting or 'correcting' measured daily latent 
    energy and sensible heat fluxes measured at an eddy covariance climate 
    station which improve closure of the surface energy balance. 
    
    The :obj:`QaQc` object has multiple options for loading of data, 
    temporal frequency adjustments of data, estimation of climatic and 
    statistical variables, and managing data and metadata with a file system.
    Input data is expected to be a :obj:`fluxdataqaqc.data.Data` instance or a
    :obj:`pandas.DataFrame` which can be used to create a :obj:`QaQc` object 
    with the :meth:`QaQc.from_dataframe` method. Monthly and daily time series 
    or raw or processed climatic data can be easily saved to disk.

    """
    # dictionary used for temporally aggregating variables
    agg_dict = {
        'energy': 'mean',
        'flux': 'mean',
        'flux_corr': 'mean',
        'br': 'mean',
        'et': 'sum',
        'et_corr': 'sum',
        'et_gap': 'sum',
        'et_fill': 'sum',
        'et_user_corr': 'sum',
        'ebr': 'mean',
        'ebr_corr': 'mean',
        'ebr_user_corr': 'mean',
        'gridMET_etr_mm': 'sum',
        'gridMET_prcp_mm': 'sum',
        'lw_in': 'mean',
        't_avg': 'mean',
        'rso': 'mean',
        'sw_pot': 'mean',
        'sw_in': 'mean',
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
 
    # EBR correction methods available
    corr_methods = (
        'ebr',
        'br'
    )

    # gridMET dict, keys are names which can be passed to download_gridMET
    gridMET_meta = {
        'etr': {
            'nc_suffix': 'agg_met_etr_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_reference_evapotranspiration_alfalfa',
            'rename': 'gridMET_etr_mm',
            'units': 'mm'
        },
        'pr': {
            'nc_suffix': 'agg_met_pr_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'precipitation_amount',
            'rename': 'gridMET_prcp_mm',
            'units': 'mm'
        },    
        'pet': {
            'nc_suffix': 'agg_met_pet_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_reference_evapotranspiration_grass',
            'rename': 'gridMET_eto_mm',
            'units': 'mm'
        },
        'sph': {
            'nc_suffix': 'agg_met_sph_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_specific_humidity',
            'rename': 'gridMET_q_kgkg',
            'units': 'kg/kg'
        },
        'srad': {
            'nc_suffix': 'agg_met_srad_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_shortwave_radiation_at_surface',
            'rename': 'gridMET_srad_wm2',
            'units': 'w/m2'
        },
        'vs': {
            'nc_suffix': 'agg_met_vs_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_wind_speed',
            'rename': 'gridMET_u10_ms',
            'units': 'm/s'
        },
        'tmmx': {
            'nc_suffix': 'agg_met_tmmx_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_maximum_temperature',
            'rename': 'gridMET_tmax_k',
            'units': 'K'
        },
        'tmmn': {
            'nc_suffix': 'agg_met_tmmn_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_minimum_temperature',
            'rename': 'gridMET_tmin_k',
            'units': 'K'
        },
    }
    
    # all potentially calculated variables for energy balance corrections
    _eb_calc_vars = (
        'br',
        'br_user_corr',
        'energy',
        'ebr',
        'ebr_corr',
        'ebr_user_corr',
        'ebc_cf',
        'ebr_5day_clim',
        'et_gap',
        'et_fill',
        'flux',
        'flux_corr',
        'flux_user_corr',
        'Kc',
        'Kc_7day_mean',
        'LE_corr',
        'H_corr'
    )
    # potentially calculated variables for ET
    _et_calc_vars = (
        'et',
        'et_corr',
        'et_user_corr'
    )

    def __init__(self, data=None):
        
        if isinstance(data, Data):
            self.config_file = data.config_file
            self.config = data.config
            self._df = data.df
            self.variables = data.variables
            self.units = data.units
            self.elevation = data.elevation
            self.latitude = data.latitude
            self.longitude = data.longitude
            self.out_dir = data.out_dir
            self.site_id = data.site_id
            # flip variable naming dict for internal use
            self.inv_map = {
                v: k for k, v in self.variables.items() if (
                    not k in self._df.columns
                )
            }
            # using 'G' in multiple g plot may overwrite G name internally
            if not 'G' in self.inv_map.values():
                user_G_name = self.variables.get('G') 
                if user_G_name:
                    self.inv_map[user_G_name] = 'G'

            self.temporal_freq = self._check_daily_freq()
            self._check_gridMET()
            # assume energy balance vars exist, will be validated upon corr
            self._has_eb_vars = True

        elif data is not None:
            print('{} is not a valid input type'.format(type(data)))
            raise TypeError("Must assign a fluxdataqaqc.data.Data object")
        else:
            self._df = None

        self.corrected = False 
        self.corr_meth = None
            
    def _check_gridMET(self):
        """
        Check if gridMET has been downloaded (file path in config), if so
        also check if dates fully intersect those of station data. If both
        conditions are not met then update :attr:`gridMET_exists` to False 
        otherwise assign True.

        Arguments:
            None

        Returns:
            None
        """
        gridfile = self.config.get('METADATA','gridMET_file_path',fallback=None)
        if gridfile is None:
            self.gridMET_exists = False
        else:
            try:
                grid_df = pd.read_csv(
                    gridfile, parse_dates=True, index_col='date'
                )
                gridMET_dates = grid_df.index
                station_dates = self.df.index
                # add var names and units to attributes
                for val in grid_df.columns:
                    meta = [
                        v for k,v in QaQc.gridMET_meta.items() if \
                            v['rename'] == val
                    ][0]
                    self.variables[meta['rename']] = meta['rename']
                    self.units[meta['rename']] = meta['units']
                # flag False if ETr was not downloaded for our purposes
                if not 'gridMET_etr_mm' in grid_df.columns:
                    self.gridMET_exists = False
                elif station_dates.isin(gridMET_dates).all():
                    self.gridMET_exists = True
                # some gridMET exists but needs to be updated for coverage
                else:
                    self.gridMET_exists = False
            except:
                print('WARNING: unable to find/read gridMET file\n {}'.format(
                    gridfile)
                )
                self.gridMET_exists = False

    def download_gridMET(self, variables=None):
        """
        Download reference ET (alfalfa) and precipitation from gridMET for
        days in flux station time series. Also has ability to download other
        or specific gridMET variables by passing a list of gridMET variable
        names. Possible names and their long form can be found in 
        :attr:`QaQc.gridMET_meta`. 

        Upon download gridMET time series for the nearest gridMET cell will be
        merged into the instances dataframe attibute :attr:`QaQc.df` and all
        gridMET variable names will have the prefix "gridMET_" for 
        identification. 
        
        The gridMET time series file will be saved to a subdirectory called
        "gridMET_data" within the directory that contains the config file
        for the current :ob:`QaQc` instance and named with the site ID and 
        gridMET cell centroid lat and long coordinates in decimal degrees.
        
        Any previously downloaded gridMET time series will be overwritten.
        
        Arguments:
            variables (None, str, list, or tuple): default None. List of gridMET
                variable names to download, if None download ETr and 
                precipitation. See the keys of the :attr:`QaQc.gridMET_meta` 
                dictionary for a list of all variables that can be downloaded 
                by this method.

        Returns:
            None

        
        """
        # opendap thredds server
        server_prefix =\
            'http://thredds.northwestknowledge.net:8080/thredds/dodsC/'

        if variables is None:
            variables = ['etr', 'pr']

        elif not isinstance(variables, (str,list,tuple)):
            print(
                'ERROR: {} is not a valid gridMET variable '
                'or list of variable names, valid options:'
                '\n{}'.format(
                    variables, ', '.join([v for v in QaQc.gridMET_meta])
                )
            )
            return
            
        if isinstance(variables, str):
            variables = list(variables)
            
        station_dates = self.df.index
        grid_dfs = []
        for i,v in enumerate(variables):
            if not v in QaQc.gridMET_meta:
                print(
                    'ERROR: {} is not a valid gridMET variable, '
                    'valid options: {}'.format(
                        v, ', '.join([v for v in QaQc.gridMET_meta])
                    )
                )
                continue
            meta = QaQc.gridMET_meta[v]
            self.variables[meta['rename']] = meta['rename']
            self.units[meta['rename']] = meta['units']
            print('Downloading gridMET var: {}\n'.format(meta['name'])) 
            netcdf = '{}{}'.format(server_prefix, meta['nc_suffix'])
            ds = xarray.open_dataset(netcdf).sel(
                lon=self.longitude, lat=self.latitude, method='nearest'
            ).drop('crs')
            df = ds.to_dataframe().loc[station_dates].rename(
                columns={meta['name']:meta['rename']}
            )
            df.index.name = 'date' # ensure date col name is 'date'
            # on first variable (if multiple) grab gridcell centroid coords
            if i == 0:
                lat_centroid = df.lat[0]
                lon_centroid = df.lon[0]

            df.drop(['lat', 'lon'], axis=1, inplace=True)

            grid_dfs.append(df)
        
        # combine data
        df = pd.concat(grid_dfs, axis=1)
        # save gridMET time series to CSV in subdirectory where config file is
        gridMET_file = self.config_file.parent.joinpath(
            'gridMET_data'
            ).joinpath('{}_{:.4f}N_{:.4f}W.csv'.format(
                self.site_id, lat_centroid, lon_centroid
            )       
        )
        gridMET_file.parent.mkdir(parents=True, exist_ok=True)
        self.config.set(
            'METADATA', 'gridMET_file_path', value=str(gridMET_file)
        )
        df.to_csv(gridMET_file)
        # rewrite config with updated gridMET file path
        with open(str(self.config_file), 'w') as outf:
            self.config.write(outf)
        # keep in memory to merge with df during correction
        #self._gridMET_df = df
        # drop previously calced vars for replacement, no join duplicates
        self._df = _drop_cols(self._df, variables)
        self._df = self._df.join(df)
        self.gridMET_exists = True
    
        
    def _check_daily_freq(self):
        """
        Check temporal frequency of input Data, resample to daily if not already

        Note:
            If one or more sub-dauly values are missing for a day the entire
            day will be replaced with a null (:obj:`numpy.nan`).
            If user QC values for filtering data are present they will also be 
            resampled to daily means, however this should not be an issue as 
            the filtering step occurs in a :obj:`fluxdataqaqc.Data` object.
        
        """

        # rename columns to internal names 
        df = self._df.rename(columns=self.inv_map)
       
        if not isinstance(df, pd.DataFrame):
            return

        freq = pd.infer_freq(df.index)

        # pd.infer_freq does not always work and may return None
        if freq and freq > 'D':
            print('WARNING: it looks like the input data temporal frequency',
                'is greater than daily, downsampling to daily, proceed with' ,
                'caution!\n')
        if freq and freq < 'D':
            print('The input data temporal frequency appears to be less than',
                'daily.\n')

        if freq is None:
            print('The input data temporal frequency was not detected.')
            freq = 'na'

        if not freq == 'D':
            print('Data is being resampled to daily temporal frequency.')
            sum_cols = [k for k,v in QaQc.agg_dict.items() if v == 'sum']
            sum_cols = list(set(sum_cols).intersection(df.columns))
            mean_cols = set(df.columns) - set(sum_cols)
            means = df.loc[:,mean_cols].apply(
                pd.to_numeric, errors='coerce').resample('D').mean()
            sums = df.loc[:,sum_cols].apply(
                pd.to_numeric, errors='coerce').resample('D').sum()
            # using numpy forces nans if 1 or more sub-daily value missing
            # having issues however creating more than expected null days
            #means = df.loc[:,mean_cols].resample('D').apply(
            #    lambda x: x.values.mean()
            #)
            #sums = df.loc[:,sum_cols].resample('D').apply(
            #    lambda x: x.values.sum()
            #)
            df = means.join(sums)

        self._df = df.rename(self.variables)
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
        Temporally resample time series data to monthly frequency based on 
        monthly means or sums based on :attr:`QaQc.agg_dict`. 
        
        Also replaces monthly means or sums with null values if less than
        75 percent of a months days are missing in the daily data 
        (:attr:`QaQc.df`).

        If a :obj:`QaQc` instance has not yet run an energy balance correction 
        i.e. :attr:`QaQc.corrected` = False before accessing :attr:`monthly_df`
        then the default routine of data correction (energy balance ratio 
        method) will be conducted.

        Arguments:
            None

        Returns:
            None

        """
        if not self.corrected and self._has_eb_vars:
            self.correct_data()

        # rename columns to internal names 
        df = self._df.rename(columns=self.inv_map).copy()
        # avoid including string QC flags because of float forcing on resample
        numeric_cols = [c for c in df.columns if not '_qc_flag' in c]
        sum_cols = [k for k,v in QaQc.agg_dict.items() if v == 'sum']
        # to avoid warning/error of missing columns 
        sum_cols = list(set(sum_cols).intersection(df.columns))
        mean_cols = list(set(numeric_cols) - set(sum_cols))
        # if data type has changed to 'obj' resample skips... 
        # make sure data exists
        if len(mean_cols) >= 1:
            means = monthly_resample(df[mean_cols], mean_cols, 'mean', 0.3)
        else:
            means = None
        if len(sum_cols) >= 1:
            sums = monthly_resample(df[sum_cols], sum_cols, 'sum', 0.3)
        else:
            sums = None
        if isinstance(means, pd.DataFrame) and isinstance(sums, pd.DataFrame):
            df = means.join(sums)
        elif isinstance(means, pd.DataFrame):
            df = means
        elif isinstance(sums, pd.DataFrame):
            df = sums
        # use monthly sums for ebr columns not means of ratio
        if set(['LE','H','Rn','G']).issubset(df.columns):
            df.ebr = (df.H + df.LE) / (df.Rn - df.G)
        if set(['LE_corr','H_corr','Rn','G']).issubset(df.columns):
            df.ebr_corr = (df.H_corr + df.LE_corr) / (df.Rn - df.G)
        if set(['LE_user_corr','H_user_corr','Rn','G']).issubset(df.columns):
            df['ebr_user_corr']=(df.H_user_corr+df.LE_user_corr) / (df.Rn-df.G)
        
        #elif how == 'aggregate':
        #    # for monthly stats not time series
        #    df = df.groupby(df.index.month).agg(self.agg_dict)
        #else:
        #    err_msg='Invalid "how" option, use "time_series" or "aggregate"'
        #    raise ValueError(err_msg)

        df.replace([np.inf, -np.inf], np.nan, inplace=True)

        return df.rename(columns=self.variables)

    def write(self, out_dir=None):
        """
        Save a copy of the "corrected" energy balance time series (default 
        correction method) including raw input. Saves two CSVs one at daily 
        and one at monthly time frequencies. 

        Arguments:
            out_dir (str or None): default None. Directory to save CSVs, if 
                None save to :attr:`out_dir` instance variable (typically 
                "output" directory where config.ini file exists).

        Returns:
            None

        Note:
            If this method is used before correcting the data according to the
            default routines in ``correct_data`` it will be done before saving. 
            To save data created by multiple correction routines, be sure to 
            run the correction and then save to different output directories 
            otherwise output files will be overwritten with the most recently 
            used correction option.
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

        if not self.corrected and self._has_eb_vars:
            self.correct_data()

        daily_outf = out_dir / '{}_daily_data.csv'.format(self.site_id)
        monthly_outf = out_dir / '{}_monthly_data.csv'.format(self.site_id)

        self.df.to_csv(daily_outf)
        self.monthly_df.to_csv(monthly_outf)

    @classmethod
    def from_dataframe(cls, df, site_id, elev_m, lat_dec_deg, var_dict):
        """
        Create a :obj:`QaQc` object from a :obj:`pandas.DataFrame` object.
        
        Arguments:
            df (:obj:`pandas.DataFrame`): DataFrame of climate variables with
                datetime index named 'date'
            site_id (str): site identifier such as station name
            elev_m (int or float): elevation of site in meters
            lat_dec_deg (float): latitude of site in decimal degrees
            var_dict (dict): dictionary that maps `flux-data-qaqc` variable
                names to user's columns in `df` e.g. {'Rn': 'netrad', ...}
                see :attr:`fluxdataqaqc.Data.variable_names_dict` for list of 
                `flux-data-qaqc` variable names
        
        Returns:
            None

        Note:
            When using this method, any output files (CSVs, plots) will be 
            saved to a directory named "output" in the current working 
            directory. 
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
        qaqc._has_eb_vars = True

        return qaqc

    def correct_data(self, meth='ebr'):
        """
        Correct turblent fluxes to close energy balance using different
        methods, default 'ebr'. 

        Currently two options are available: 'ebr' (Energy Balance Ratio) and 
        'br' (Bowen Ratio). If you use one method followed by another corrected
        versions of LE, H, ET, and EBR will be overwritten with the most 
        recently used approach. Also computes potential clear sky radiation 
        (saved as *rso*) using the ASCE approach based on station elevation and 
        latitude. ET is calculated from raw and corrected LE using daily air
        temperature to correct the latent heat of vaporization, if air temp. is
        not available in the input data then air temp. is assumed at 20 
        degrees celcius.

        Corrected or otherwise newly calculated variables are named using the
        following suffixes to distinguish them::

          uncorrected LE, H, etc. from input data have no suffix
          _corr uses adjusted LE, H, etc. from the correction method used 
          _user_corr uses corrected LE, H, etc. found in data file (if provided)

        Arguments:
            meth (str): default 'ebr'. Method to correct energy balance.
        
        Returns
            None

        Note:
            The *ebr_corr* variable or energy balance closure ratio is 
            calculated from the corrected versions of LE and H independent 
            of the method. When using the 'ebr' method the energy balance 
            correction factor (what is applied to the raw H and LE) is left as 
            calculated (inverse of ebr) and saved as *ebc_cf*. 
        """

        # in case starting in Python and no df assigned yet
        if not isinstance(self._df, pd.DataFrame):
            print('Please assign a dataframe of acceptable data first!')
            return
        if meth not in self.corr_methods:
            err_msg = ('ERROR: {} is not a valid correction option, please '
                'use one of the following: {}'.format(meth, ', '.join(
                    [el for el in self.corr_methods]))
            )
            raise ValueError(err_msg)

        # calculate clear sky radiation if not already computed
        self._calc_rso()
        # energy balance corrections
        if not set(['Rn','LE','H','G']).issubset(self.variables.keys()) or\
                not set(['Rn','LE','H','G']).issubset(
                    self.df.rename(columns=self.inv_map).columns):
            print(
                '\nMissing one or more energy balance variables, cannot perform'
                ' energy balance correction.'
            )
            self._has_eb_vars = False
            return
        if meth == 'ebr':
            self._ebr_correction()
        if meth == 'br':
            self._bowen_ratio_correction()

        self.corr_meth = meth
        # calculate raw, corrected ET 
        self._calc_et()
        # fill gaps of corr ET with ET from smoothed and gap filled Kc*ETr
        self._ETr_gap_fill()

        # update inv map for naming
        self.inv_map = {
            v: k for k, v in self.variables.items() if (
                not v.replace('_mean', '') == k or not k in self.df.columns)
        }
        # using 'G' in multiple g plot may overwrite G name internally
        if not 'G' in self.inv_map.values():
            user_G_name = self.variables.get('G')
            self.inv_map[user_G_name] = 'G'

    def _ETr_gap_fill(self, et_name='et_corr'):
        """
        Use gridMET reference ET to calculate ET from Kc, smooth and gap fill
        calced ET and then use to fill gaps in corrected ET. Keeps tabs on 
        number of days in ET that were filled in each month.
        
        Keyword Arguments:
            et_name (str): default "et_corr". Name of ET variable to use when 
                calculating the crop coefficient and to fill gaps with 
                calculated ET.

        Returns:
            None

        """
        # get ETr if not on disk
        if not self.gridMET_exists:
            self.download_gridMET()

        else:
            # gridMET file has been verified and has all needed dates, just load
            gridfile = self.config.get(
                'METADATA','gridMET_file_path',fallback=None
            )
            print(
                'gridMET reference ET already downloaded for station at:\n'
                '{}\nnot redownloading.'.format(gridfile)
            )
            grid_df = pd.read_csv(gridfile, parse_dates=True, index_col='date')
            # drop previously calced vars for replacement, no join duplicates
            self._df = _drop_cols(self._df, list(grid_df.columns))
            self._df = self._df.join(grid_df)
            for c in grid_df.columns:
                self.variables[c] = c

        # calc Kc 7 day moving average, min vals in window = 2, linear interp
        df = self.df.rename(columns=self.inv_map)
        df['Kc'] = df[et_name] / df.gridMET_etr_mm
        df['Kc_7day_mean'] = df.Kc.rolling(7, min_periods=2, center=True).mean()
        df.Kc_7day_mean = df.Kc_7day_mean.interpolate(method='linear')
        # calc ET from Kc and ETr
        df['et_fill'] = df.gridMET_etr_mm * df.Kc_7day_mean
        # flag days in corrected ET that are missing and sum for each month
        df['et_gap'] = False
        df.loc[(df[et_name].isna() & df.et_fill.notna()), 'et_gap'] = True
        df.loc[df.et_gap, et_name] = df.et_fill
        # et fill values only where they were used to gap fill, for plotting
        df['et_fill_val'] = np.nan
        df.loc[df.et_gap , 'et_fill_val'] = df.et_fill

        # update variables attribute with new variables (may vary)
        new_cols = set(df.columns) - set(self.variables)
        for el in new_cols:
            self.variables[el] = el
        self.units['Kc'] = 'unitless'
        self.units['Kc_7day_mean'] = 'unitless'
        self.units['et_fill'] = 'mm'
        self.units['et_fill_val'] = 'mm'
        self.units['et_gap'] = 'unitless'

        self._df = df


    def _calc_et(self):
        """
        Calculate daily ET (mm) from raw and corrected LE (w/m2)), if air 
        temperature is available use to correct latent heat of vaporization.
        
        Currently computes on :attr:`df` dataframe attribute in place. Can be 
        called before or after energy balance closure corrections, if before 
        the only raw ET will be calculated. Assumes at least LE raw exists in 
        dataframe.

        Arguments:
            None

        Returns:
            None
        """

        # drop relavant calculated variables if they exist
        self._df = _drop_cols(self._df, self._et_calc_vars)
        df = self.df.rename(columns=self.inv_map)
        
        # LH from L.P. Harrison (1963)
        if 't_avg' in df.columns:
            df['et'] = 86400 * (df.LE/(2501000 - (2361 * df.t_avg)))
            if 'LE_corr' in df.columns:
                df['et_corr']=86400 * (df.LE_corr/(2501000 - (2361 * df.t_avg)))
            if 'LE_user_corr' in df.columns:
                df['et_user_corr'] =\
                    86400*(df.LE_user_corr/(2501000 - (2361 * df.t_avg)))
        # otherwise assume air temp = 20 degrees C
        else:
            df['et'] = 86400 * (df.LE/(2501000 - (2361 * 20)))
            if 'LE_corr' in df.columns:
                df['et_corr'] = 86400 * (df.LE_corr/(2501000 - (2361 * 20)))
            if 'LE_user_corr' in df.columns:
                df['et_user_corr']=86400*(df.LE_user_corr/(2501000-(2361*20)))
        
        # update variables attribute with new variables (may vary)
        new_cols = set(df.columns) - set(self.variables)
        for el in new_cols:
            self.variables[el] = el
            self.units[el] = 'mm'

        # join data back into df attr
        self._df = df.rename(columns=self.variables)

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
        if 'rso' in self.df.columns:
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
        self.units['rso'] = 'w/m2'

    def _ebr_correction(self):
        """
        Energy balance ratio correction method for daily LE, H, EBR, and ET.

        Correct turblent fluxes to close energy balance using methods 
        described by FLUXNET for daily H and LE. 

        Updates :attr:`QaQc.df` and :attr:`QaQc.variables` attributes with new 
        variables related to the corrections, e.g. LE_corr, ebr_corr, etc.

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
        
        # drop relavant calculated variables if they exist
        self._df = _drop_cols(self.df, self._eb_calc_vars)
        df = self._df.rename(columns=self.inv_map)
        # check if mean G, or heat storage measurements exist
        vars_to_use = ['LE','H','Rn','G']
        for el in ['SH', 'SLE', 'SG']:
            if el in df.columns:
                vars_to_use.append(el)
        df = df[vars_to_use].astype(float).copy()

        # make copy of original data for later
        orig_df = df[['LE','H','Rn','G']].astype(float).copy()
        orig_df['ebr'] =  (orig_df.H + orig_df.LE) / (orig_df.Rn - orig_df.G)
        # compute IQR to filter out extreme ebrs, 
        df['ebr'] = (df.H + df.LE) / (df.Rn - df.G)
        Q1 = df['ebr'].quantile(0.25)
        Q3 = df['ebr'].quantile(0.75)
        IQR = Q3 - Q1
        # filter values between Q1-1.5IQR and Q3+1.5IQR
        filtered = df.query('(@Q1 - 1.5 * @IQR) <= ebr <= (@Q3 + 1.5 * @IQR)')
        # apply filter
        filtered_mask = filtered.index
        removed_mask = set(df.index) - set(filtered_mask)
        removed_mask = pd.to_datetime(list(removed_mask))
        df.loc[removed_mask] = np.nan

        # FLUXNET methods 1 and 2 for filtering/smoothing ebr
        ebr = df.ebr.values
        df['ebr_corr'] = np.nan
        for i in range(len(ebr)):
            win_arr1 = ebr[i-half_win_1:i+half_win_1+1]
            win_arr2 = ebr[i-half_win_2:i+half_win_2+1]
            count = np.count_nonzero(~np.isnan(win_arr1))
            # get median of daily window1 if half window2 or more days exist
            if count >= half_win_2:
                val = np.nanpercentile(win_arr1, 50, axis=None)
                if abs(1/val) >= 2:
                    val = np.nan
            # if at least one day exists in window2 take mean
            elif np.count_nonzero(~np.isnan(win_arr2)) > 0:
                #val = np.nanmean(win_arr2)
                val = np.nanmedian(win_arr2)
                if abs(1/val) >= 2:
                    val = np.nan
            else:
                # assign nan for now, update with 5 day climatology
                val = np.nan
            # assign values if they were found in methods 1 or 2
            df.iloc[i, df.columns.get_loc('ebr_corr')] = val
        # make 5 day climatology of ebr for method 3
        # the cenetered window skips first and last 5 DOYs
        # so prepend and append first and last 5 days and loop...
        doy_ebr_mean=df['ebr_corr'].groupby(df.index.dayofyear).mean().copy()
        l5days = pd.Series(
            index=np.arange(-4,1), data=doy_ebr_mean[-5:].values)
        f5days = pd.Series(
            index=np.arange(367,372), data=doy_ebr_mean[:5].values)
        doy_ebr_mean = doy_ebr_mean.append(f5days)
        doy_ebr_mean = pd.concat([l5days, doy_ebr_mean])
        ebr_5day_clim = pd.DataFrame(
            index=np.arange(1,367), columns=['ebr_5day_clim'])
        doy_ebr_mean = doy_ebr_mean.values
        for i in range(len(doy_ebr_mean)):
            # i = 0 which starts at prepended 5 days, shift window up
            win = doy_ebr_mean[i:i+2*half_win_2+1]
            count = np.count_nonzero(~np.isnan(win))
            # get 11 day moving window mean
            if i in ebr_5day_clim.index and count > 0:
                ebr_5day_clim.iloc[
                    i-1, ebr_5day_clim.columns.get_loc('ebr_5day_clim')
                ] = np.nanmean(win)
        ebr_5day_clim['DOY'] = ebr_5day_clim.index
        ebr_5day_clim.index.name = 'date'

        # fill gaps of 11 or more days in filtered EBR with 5 day clim
        df['DOY'] = df.index.dayofyear
        # datetime indices of all remaining null elements
        null_dates = df.loc[df.ebr_corr.isnull(), 'ebr_corr'].index
        merged = pd.merge(
            df, ebr_5day_clim, on='DOY', how='left', right_index=True
        )
        # assign 5 day climatology of EBR 
        merged.loc[null_dates,'ebr_corr'] =\
            merged.loc[null_dates,'ebr_5day_clim']
        # replace raw variables with unfiltered dataframe copy
        merged.LE = orig_df.LE
        merged.H = orig_df.H
        merged.Rn = orig_df.Rn
        merged.G = orig_df.G
        merged.ebr = orig_df.ebr
        # calculated corrected EBR to assign to ebr_corr (not EBC_CF), 
        # save CFs as defined by fluxnet method, i.e. inverse of EBR
        merged['ebc_cf'] = 1/merged.ebr_corr
        # filter out CF that are 2 or higher (absolute), from 5 day climo ebr_c
        merged.loc[abs(merged.ebc_cf) >= 2, 'ebc_cf'] = np.nan
        # apply corrections to LE and H multiply by 1/EBR
        merged['LE_corr'] = merged.LE * merged.ebc_cf
        merged['H_corr'] = merged.H * merged.ebc_cf
        # filter out any corrected LE that <= -100 or >= 850 w/m2
        # also removing corrected H, EBR, and EBC_CF
        merged.loc[
            (merged.LE_corr >= 850) | (merged.LE_corr <= -100), (
                'LE_corr', 'H_corr', 'ebr_corr', 'ebc_cf'
            )
        ] = np.nan
        # compute EBR total turb flux
        merged['flux_corr'] = merged['LE_corr'] + merged['H_corr']

        df = self._df.rename(columns=self.inv_map)
        # other variables needed for plots using raw data
        df['flux'] = merged.LE + merged.H
        df['energy'] = merged.Rn - merged.G

        # corrected turbulent flux if given from input data
        if set(['LE_user_corr','H_user_corr']).issubset(df.columns):
            df['flux_user_corr'] = df.LE_user_corr + df.H_user_corr 
            df['ebr_user_corr']=(df.H_user_corr+df.LE_user_corr)/(df.Rn - df.G)
            df.ebr_user_corr=df.ebr_user_corr.replace([np.inf,-np.inf], np.nan)
            self.variables.update(
                flux_user_corr = 'flux_user_corr',
                ebr_user_corr = 'ebr_user_corr'
            )
        # grab select columns to merge into main dataframe
        cols = list(set(merged.columns).difference(df.columns))
        # join calculated data in
        merged = df.join(merged[cols], how='outer')
        merged.drop('DOY', axis=1, inplace=True)

        self.variables.update(
            energy = 'energy',
            flux = 'flux',
            LE_corr = 'LE_corr',
            H_corr = 'H_corr',
            flux_corr = 'flux_corr',
            ebr = 'ebr',
            ebr_corr = 'ebr_corr',
            ebc_cf = 'ebc_cf',
            ebr_5day_clim = 'ebr_5day_clim'
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
        
        Arguments:
            None

        Returns:
            None
        """
        # drop relavant calculated variables if they exist
        self._df = _drop_cols(self.df, self._eb_calc_vars)
        df = self._df.rename(columns=self.inv_map)

        # apply correction 
        df['br'] = df.H / df.LE
        df['LE_corr'] = (df.Rn - df.G) / (1 + df.br)
        df['H_corr'] = df.LE_corr * df.br
        df['flux_corr'] = df.LE_corr + df.H_corr

        # add EBR, other vars to dataframe using LE and H from raw, corr 
        # if provided user corrected, add raw energy and flux
        if set(['LE_user_corr','H_user_corr']).issubset(df.columns):
            df['flux_user_corr'] = df.LE_user_corr + df.H_user_corr 
            df['br_user_corr'] = df.H_user_corr / df.LE_user_corr 
            df['ebr_user_corr']=(df.H_user_corr+df.LE_user_corr)/(df.Rn - df.G)
            df.ebr_user_corr=df.ebr_user_corr.replace([np.inf,-np.inf], np.nan)

            self.variables.update(
                flux_user_corr = 'flux_user_corr',
                ebr_user_corr = 'ebr_user_corr',
                br_user_corr = 'br_user_corr'
            )
        df['ebr'] = (df.H + df.LE) / (df.Rn - df.G)
        df['ebr_corr'] = (df.H_corr + df.LE_corr) / (df.Rn - df.G)
        df['energy'] = df.Rn - df.G
        df['flux'] = df.LE + df.H

        self.variables.update(
            br = 'br',
            energy = 'energy',
            flux = 'flux',
            LE_corr = 'LE_corr',
            H_corr = 'H_corr',
            flux_corr = 'flux_corr',
            ebr = 'ebr',
            ebr_corr = 'ebr_corr'
        )
        # replace undefined/infinity with nans in all EBR columns
        df.ebr = df.ebr.replace([np.inf, -np.inf], np.nan)
        df.ebr_corr = df.ebr_corr.replace([np.inf, -np.inf], np.nan)

        # revert column names to user's
        self._df = df.rename(columns=self.variables)
        # update flag for other methods
        self.corrected = True

    def plot(self, ncols=1, group='all', output_type='save', out_file=None):
        """
        """
        # create aggregrated plot structure from fluxdataqaqc.Plot._plot() 
        self._plot(
            self, ncols=ncols, output_type=output_type, out_file=out_file
        )

def _drop_cols(df, cols):
    """Drop columns from dataframe if they exist """
    for c in cols:
        if c in df.columns:
            df.drop(c, axis=1, inplace=True)

    return df

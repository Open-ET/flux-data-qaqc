# -*- coding: utf-8 -*-
""" 
Tools for correcting energy-balance components to improve energy balance
closure and other data management, validation and scientific analysis tools. 
"""

from pathlib import Path

import xarray
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from refet.calcs import _ra_daily, _rso_simple
import refet 

from .data import Data
from .plot import Plot
from .util import monthly_resample, Convert

class QaQc(Plot, Convert):
    """
    Numerical routines for correcting daily energy balance closure 
    for eddy covariance data and other data analysis tools.
    
    Two routines are provided for improving energy balance closure by adjusting
    turbulent fluxes, latent energy and sensible heat, the Energy Balance Ratio
    method (modified from `FLUXNET
    <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__)
    and the Bowen Ratio method. 
    
    The :obj:`QaQc` object also has multiple tools for temporal frequency
    aggregation and resampling, estimation of climatic and statistical
    variables (e.g. ET and potential shortwave radiation), downloading gridMET
    reference ET, managing data and metadata, interactive validation plots, and
    managing a structure for input and output data files. Input data is
    expected to be a :obj:`.Data` instance or a
    :obj:`pandas.DataFrame`. 

    Keyword Arguments:
        data (:obj:`.Data`): :obj:`.Data` instance to create :obj:`.QaQc` 
            instance.
        drop_gaps (bool): default :obj:`True`. If :obj:`True` automatically 
            filter variables on days with sub-daily measurement gaps less than
            ``daily_frac``.
        daily_frac (float): default 1.00. Fraction of sub-daily data required
            otherwise the daily value will be filtered out if ``drop_gaps`` is
            :obj:`True`. E.g. if ``daily_frac = 0.5`` and the input data is
            hourly, then data on days with less than 12 hours of data will be
            forced to null within :attr:`QaQc.df`. This is important because 
            systematic diurnal gaps will affect the autmoatic resampling that
            occurs when creating a :obj:`QaQc` instance and the daily data is 
            used in closure corrections, other calculations, and plots. If 
            sub-daily linear interpolation is applied to energy balance 
            variables the gaps are counted *after* the interpolation.
        max_interp_hours (None or float): default 2. Length of largest gap to 
            fill with linear interpolation in energy balance variables if 
            input datas temporal frequency is less than daily. This value will
            be used to fill gaps when :math:`Rn > 0` or :math:`Rn` is missing
            during each day.
        max_interp_hours_night (None or float): default 4. Length of largest gap
            to fill with linear interpolation in energy balance variables if 
            input datas temporal frequency is less than daily when 
            :math:`Rn < 0` within 12:00PM-12:00PM daily intervals. 
    
    Attributes:
        agg_dict (dict): Dictionary with internal variable names as keys and
            method of temporal resampling (e.g. "mean" or "sum") as values. 
        config (:obj:`configparser.ConfigParser`): Config parser instance
            created from the data within the config.ini file.
        config_file (:obj:`pathlib.Path`): Absolute path to config.ini file 
            used for initialization of the :obj:`fluxdataqaqc.Data` instance
            used to create the :obj:`QaQc` instance. 
        corrected (bool): False until an energy balance closure correction has
            been run by calling :meth:`QaQc.correct_data`.
        corr_methods (tuple): List of Energy Balance Closure correction routines
            usable by :meth:`QaQc.correct_data`.
        corr_meth (str or None): Name of most recently applied energy balance
            closure correction.
        elevation (float): Site elevation in meters.
        gridMET_exists (bool): True if path to matching gridMET time series
            file exists on disk and has time series for reference ET and 
            precipitation and the dates for these fully overlap with the energy
            balance variables, i.e. the date index of :attr:`QaQc.df`.
        gridMET_meta (dict): Dictionary with information for gridMET variables 
            that may be downloaded using :meth:`QaQc.download_gridMET`.
        inv_map (dict): Dictionary with input climate file names as keys and 
            internal names as values. May only include pairs when they differ.
        latitude (float): Site latitude in decimal degrees.
        longitude (float): Site longitude in decimal degrees.
        out_dir (pathlib.Path): Default directory to save output of 
            :meth:`QaQc.write` or :meth:`QaQc.plot` methods.
        n_samples_per_day (int): If initial time series temporal frequency is
            less than 0 then this value will be updated to the number of samples
            detected per day, useful for post-processing based on the count of
            sub-daily gaps in energy balance variables, e.g. "LE_subday_gaps".
        plot_file (pathlib.Path or None): path to plot file once it is 
            created/saved by :meth:`QaQc.plot`. 
        site_id (str): Site ID.
        temporal_freq (str): Temporal frequency of initial (as found in input 
            climate file) data as determined by :func:`pandas.infer_freq`.
        units (dict): Dictionary with internal variable names as keys and 
            units as found in config as values.
        variables (dict): Dictionary with internal variable names as keys and 
            names as found in the input data as values.

    Note:
        Upon initialization of a :obj:`QaQc` instance the temporal frequency of
        the input data checked using :func:`pandas.infer_freq` which does not
        always correctly parse datetime indices, if it is not able to correctly
        determine the temporal frequency the time series will be resampled to 
        daily frequency but if it is in fact already at daily frequency the data
        will be unchanged. In this case the :attr:`QaQc.temporal_freq` will be
        set to "na".
    """
    # dictionary used for temporally aggregating variables
    agg_dict = {
        'ASCE_ETo': 'sum',
        'ASCE_ETr': 'sum',
        'energy': 'mean',
        'flux': 'mean',
        'flux_corr': 'mean',
        'br': 'mean',
        'ET': 'sum',
        'ET_corr': 'sum',
        'ET_gap': 'sum',
        'ET_fill': 'sum',
        'ET_fill_val': 'sum',
        'ET_user_corr': 'sum',
        'ebr': 'mean',
        'ebr_corr': 'mean',
        'ebr_user_corr': 'mean',
        'ebr_5day_clim': 'mean',
        'gridMET_ETr': 'sum',
        'gridMET_ETo': 'sum',
        'gridMET_prcp': 'sum',
        'lw_in': 'mean',
        't_avg': 'mean',
        't_max': 'mean',
        't_min': 'mean',
        't_dew': 'mean',
        'rso': 'mean',
        'sw_pot': 'mean',
        'sw_in': 'mean',
        'vp': 'mean',
        'vpd': 'mean',
        'ppt': 'sum',
        'ppt_corr': 'sum',
        'ws': 'mean',
        'Rn': 'mean',
        'Rn_subday_gaps': 'sum',
        'rh' : 'mean',
        'sw_out': 'mean',
        'lw_out': 'mean',
        'G': 'mean',
        'G_subday_gaps': 'sum',
        'LE': 'mean',
        'LE_corr': 'mean',
        'LE_subday_gaps': 'sum',
        'LE_user_corr': 'mean',
        'H': 'mean',
        'H_corr': 'mean',
        'H_subday_gaps': 'sum',
        'H_user_corr': 'mean',
    }

    # EBR correction methods available
    corr_methods = (
        'ebr',
        'br',
        'lin_regress'
    )

    # gridMET dict, keys are names which can be passed to download_gridMET
    gridMET_meta = {
        'ETr': {
            'nc_suffix': 'agg_met_etr_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_reference_evapotranspiration_alfalfa',
            'rename': 'gridMET_ETr',
            'units': 'mm'
        },
        'pr': {
            'nc_suffix': 'agg_met_pr_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'precipitation_amount',
            'rename': 'gridMET_prcp',
            'units': 'mm'
        },    
        'pet': {
            'nc_suffix': 'agg_met_pet_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_reference_evapotranspiration_grass',
            'rename': 'gridMET_ETo',
            'units': 'mm'
        },
        'sph': {
            'nc_suffix': 'agg_met_sph_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_specific_humidity',
            'rename': 'gridMET_q',
            'units': 'kg/kg'
        },
        'srad': {
            'nc_suffix': 'agg_met_srad_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_shortwave_radiation_at_surface',
            'rename': 'gridMET_srad',
            'units': 'w/m2'
        },
        'vs': {
            'nc_suffix': 'agg_met_vs_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_mean_wind_speed',
            'rename': 'gridMET_u10',
            'units': 'm/s'
        },
        'tmmx': {
            'nc_suffix': 'agg_met_tmmx_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_maximum_temperature',
            'rename': 'gridMET_tmax',
            'units': 'K'
        },
        'tmmn': {
            'nc_suffix': 'agg_met_tmmn_1979_CurrentYear_CONUS.nc#fillmismatch',
            'name': 'daily_minimum_temperature',
            'rename': 'gridMET_tmin',
            'units': 'K'
        },
    }

    # all potentially calculated variables for energy balance corrections
    _eb_calc_vars = (
        'br',
        'br_user_corr',
        'energy',
        'energy_corr',
        'ebr',
        'ebr_corr',
        'ebr_user_corr',
        'ebc_cf',
        'ebr_5day_clim',
        'flux',
        'flux_corr',
        'flux_user_corr',
        'G_corr',
        'H_corr',
        'LE_corr',
        'Rn_corr'
    )
    # potentially calculated variables for ET
    _et_calc_vars = (
        'ET',
        'ET_corr',
        'ET_user_corr'
    )
    # potentially calculated ET gap fill variables
    _et_gap_fill_vars = (
        'ET_gap',
        'ET_fill',
        'ET_fill_val',
        'ETrF',
        'ETrF_filtered',
        'EToF',
        'EToF_filtered'
    )

    def __init__(self, data=None, drop_gaps=True, daily_frac=1.00, 
            max_interp_hours=2, max_interp_hours_night=4):
        
        if isinstance(data, Data):
            self.config_file = data.config_file
            self.config = data.config
            data.df.head();# need to access to potentially calc vp/vpd
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

            # data will be loaded if it has not yet via Data.df
            self.temporal_freq = self._check_daily_freq(
                drop_gaps, daily_frac, max_interp_hours, max_interp_hours_night
            )
            # check units, convert if possible for energy balance, ppt, Rs, vp,
            self._check_convert_units()
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

    def daily_ASCE_refET(self, reference='short', anemometer_height=None):
        """
        Calculate daily ASCE standardized short (ETo) or tall (ETr) reference 
        ET from input data and wind measurement height.

        The resulting time series will automatically be merged into the
        :attr:`.Data.df` dataframe named "ASCE_ETo" or "ASCE_ETr" respectively.

        Keyword Arguments:
            reference (str): default "short", calculate tall or short ASCE
                reference ET.
            anemometer_height (float or None): wind measurement height in meters
                , default :obj:`None`. If :obj:`None` then look for the
                "anemometer_height" entry in the **METADATA** section of the
                config.ini, if not there then print a warning and use 2 meters.

        Returns:
            :obj:`None` 

        Note:
            If the hourly ASCE variables were prior calculated from a
            :obj:`.Data` instance they will be overwritten as they are saved
            with the same names. 
        """

        df = self.df.rename(columns=self.inv_map)

        req_vars = ['vp', 'ws', 'sw_in', 't_min', 't_max']

        if not set(req_vars).issubset(df.columns):
            print('Missing one or more required variables, cannot compute')
            return

        if anemometer_height is None:
            anemometer_height = self.config.get(
                'METADATA', 'anemometer_height', fallback=None
            )
            if anemometer_height is None:
                print(
                    'WARNING: anemometer height was not given and not found in '
                    'the config files metadata, proceeding with height of 2 m'
                )
                anemometer_height = 2

        # RefET will convert to MJ-m2-hr
        input_units = {
            'rs': 'w/m2'
        }

        length = len(df.t_min)
        tmin = df.t_min
        tmax = df.t_max
        rs = df.sw_in
        ea = df.vp
        uz = df.ws
        zw = np.full(length, anemometer_height)
        lat = np.full(length, self.latitude)
        doy = df.index.dayofyear
        elev = np.full(length, self.elevation)

        REF = refet.Daily(
            tmin,
            tmax,
            ea,
            rs,
            uz,
            zw,
            elev,
            lat,
            doy,
            method='asce',
            input_units=input_units,
        )

        if reference == 'short':
            ret = REF.eto()
            name = 'ASCE_ETo'
        elif reference == 'tall':
            ret = REF.etr()
            name = 'ASCE_ETr'

        # can add directly into QaQc.df
        df[name] = ret
        self._df = df.rename(columns=self.variables)
        self.variables[name] = name
        self.units[name] = 'mm'

    def _check_convert_units(self):
        """
        Verify if units are recognized for variables in QaQc.allowable_units,
        next verify that they have required units as in QaQc.required_units
        if not convert them.

        Conversions are handled by util.Convert.convert class method.
        """
        # force all input units to lower case
        for k, v in self.units.items():
            self.units[k] = v.lower()

        # can add check/rename unit aliases, e.g. C or c or celcius, etc... 

        df = self._df.rename(columns=self.inv_map)

        for v, u in self.units.items():
            if not v in QaQc.required_units.keys():
                # variable is not required to have any particular unit, skip
                continue
            elif not u in QaQc.allowable_units[v]:
                print('ERROR: {} units are not recognizable for var: {}\n'
                    'allowable input units are: {}\nNot converting.'.format(
                        u, v, ','.join(QaQc.allowable_units[v])
                    )
                )
            elif not u == QaQc.required_units[v]:
                # do conversion, update units
                # pass variable, initial unit, unit to be converted to, df
                df = Convert.convert(v, u, QaQc.required_units[v], df)
                self.units[v] = QaQc.required_units[v]

        self._df = df
            
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
                if not {'gridMET_ETr','gridMET_ETo'}.issubset(grid_df.columns):
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
        Download reference ET (alfalfa and grass) and precipitation from
        gridMET for all days in flux station time series by default. 
        
        Also has ability to download other specific gridMET variables by
        passing a list of gridMET variable names. Possible variables and their
        long form can be found in :attr:`QaQc.gridMET_meta`. 

        Upon download gridMET time series for the nearest gridMET cell will be
        merged into the instances dataframe attibute :attr:`QaQc.df` and all
        gridMET variable names will have the prefix "gridMET\_" for 
        identification. 
        
        The gridMET time series file will be saved to a subdirectory called
        "gridMET_data" within the directory that contains the config file
        for the current :obj:`QaQc` instance and named with the site ID and 
        gridMET cell centroid lat and long coordinates in decimal degrees.
        
        
        Arguments:
            variables (None, str, list, or tuple): default None. List of gridMET
                variable names to download, if None download ETr and 
                precipitation. See the keys of the :attr:`QaQc.gridMET_meta` 
                dictionary for a list of all variables that can be downloaded 
                by this method.

        Returns:
            :obj:`None`

        Note: 
            Any previously downloaded gridMET time series will be overwritten
            when calling the method, however if using the the gap filling
            method of the "ebr" correction routine the download will not
            overwrite currently existing data so long as gridMET reference ET
            and precipitation is on disk and its path is properly set in the
            config file.
        
        """
        # opendap thredds server
        server_prefix =\
            'http://thredds.northwestknowledge.net:8080/thredds/dodsC/'

        if variables is None:
            variables = ['ETr', 'pet', 'pr']

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
        # drop previously calced vars for replacement, no join duplicates
        self._df = _drop_cols(self._df, variables)
        self._df = self._df.join(df)
        self.gridMET_exists = True
    
        
    def _check_daily_freq(self, drop_gaps, daily_frac, max_interp_hours,
            max_interp_hours_night):
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
        second_day = df.index.date[2]
        third_day = second_day + pd.Timedelta(1, unit='D')

        # pd.infer_freq does not always work and may return None
        if freq and freq > 'D':
            pass
        elif freq and freq < 'D':
            print('\nThe input data temporal frequency appears to be less than',
                'daily.')

        # slice is transposed if only one date entry
        elif len(df.loc[str(third_day)].index) == len(df.columns) and \
                (df.loc[str(third_day)].index == df.columns).all():
            print('\nInput temporal frequency is already daily.')
            freq = 'D'
            # add missing dates (if exist) for full time series records/plots
            idx = pd.date_range(df.index.min(), df.index.max())
            df = df.reindex(idx)
            df.index.name = 'date'

        elif freq is None:
            freq = 'na'

        if not freq == 'D':
            # find frequency manually, optionally drop days with subdaily gaps
            # see if two adj. dates exist, skip first day in case it is not full
            max_times_in_day = len(df.loc[str(third_day)].index) 
            self.n_samples_per_day = max_times_in_day
            downsample = False
            if daily_frac > 1:
                print('ERROR: daily_frac must be between 0 and 1, using 1')
                daily_frac = 1
            elif daily_frac < 0:
                print('ERROR: daily_frac must be between 0 and 1, using 0')
                daily_frac = 0
            if not str(third_day) in df.index and drop_gaps:
                print('WARNING: it looks like the input temporal frequency',
                    'is greater than daily, downsampling, proceed with' ,
                    'caution!\n')
                downsample = True

            energy_bal_vars = ['LE', 'H', 'Rn', 'G']
            asce_std_interp_vars = ['t_avg','sw_in','ws','vp']
            # for interpolation of energy balance variables if they exist
            energy_bal_vars = list(
                set(energy_bal_vars).intersection(df.columns)
            )
            asce_std_interp_vars = list(
                set(asce_std_interp_vars).intersection(df.columns)
            )
            interp_vars = asce_std_interp_vars + energy_bal_vars

            # add subday gap col for energy balance comps to sum daily gap count
            for v in energy_bal_vars:
                df['{}_subday_gaps'.format(v)] = False
                df.loc[df[v].isna(), '{}_subday_gaps'.format(v)] = True

            print('Data is being resampled to daily temporal frequency.')
            sum_cols = [k for k,v in self.agg_dict.items() if v == 'sum']
            sum_cols = list(set(sum_cols).intersection(df.columns))
            mean_cols = set(df.columns) - set(sum_cols)

            means = df.loc[:,mean_cols].apply(
                pd.to_numeric, errors='coerce').resample('D').mean().copy()
            # issue with resample sum of nans, need to drop first else 0
            sums = df.loc[:,sum_cols].dropna().apply(
                pd.to_numeric, errors='coerce').resample('D').sum()

            if max_interp_hours is not None:
                # linearly interpolate energy balance components only
                max_gap = int(round(max_interp_hours/24 * max_times_in_day))
                max_night_gap = int(
                    round(max_interp_hours_night/24 * max_times_in_day)
                )
                print(
                    'Linearly interpolating gaps in energy balance components '
                    'up to {} hours when Rn < 0 and up to {} hours when '
                    'Rn >= 0.'.format(max_interp_hours_night, max_interp_hours)
                )

                tmp = df[interp_vars].apply(
                    pd.to_numeric, errors='coerce'
                )

                if 'Rn' in energy_bal_vars:
                    grped_night = tmp.loc[
                        (tmp.Rn < 0) & (tmp.Rn.notna())
                    ].copy()
                    grped_night.drop_duplicates(inplace=True)
                    grped_night = grped_night.groupby(
                        pd.Grouper(freq='24H', offset='12H')).apply(
                            lambda x: x.interpolate(
                                method='linear', limit=max_night_gap, 
                                limit_direction='both', limit_area='inside'
                            )
                        )
                    grped_day = tmp.loc[(tmp.Rn >= 0) | (tmp.Rn.isna())].copy()
                    grped_day.drop_duplicates(inplace=True)
                    grped_day = grped_day.groupby(
                        pd.Grouper(freq='24H')).apply(
                            lambda x: x.interpolate(
                                method='linear', limit=max_gap, 
                                limit_direction='both', limit_area='inside'
                            )
                        )
                else:
                    grped_night = tmp.copy()
                    grped_night.drop_duplicates(inplace=True)
                    grped_night = grped_night.groupby(
                        pd.Grouper(freq='24H', offset='12H')).apply(
                            lambda x: x.interpolate(
                                method='linear', limit=max_night_gap, 
                                limit_direction='both', limit_area='inside'
                            )
                        )
                    grped_day = tmp.copy()
                    grped_day.drop_duplicates(inplace=True)
                    grped_day = grped_day.groupby(
                        pd.Grouper(freq='24H')).apply(
                            lambda x: x.interpolate(
                                method='linear', limit=max_gap, 
                                limit_direction='both', limit_area='inside'
                            )
                        )

                # get full datetime index from grouper operation
                if type(grped_night.index) is pd.MultiIndex: 
                    grped_night.index = grped_night.index.get_level_values(1)

                if type(grped_day.index) is pd.MultiIndex:
                    grped_day.index = grped_day.index.get_level_values(1)

                interped = pd.concat([grped_day, grped_night])
                if interped.index.duplicated().any():
                    interped = interped.loc[
                        ~interped.index.duplicated(keep='first')
                    ]
                # overwrite non-interpolated daily means
                means[interp_vars] =\
                    interped[interp_vars].resample('D').mean().copy()
                if 't_avg' in interp_vars:
                    means['t_min'] = interped.t_avg.resample('D').min()
                    means['t_max'] = interped.t_avg.resample('D').max()
                    self.variables['t_min'] = 't_min'
                    self.units['t_min'] = self.units['t_avg']
                    self.variables['t_max'] = 't_max'
                    self.units['t_max'] = self.units['t_avg']
                    interp_vars = interp_vars + ['t_min','t_max']
                    interped[['t_min','t_max']] = means[['t_min','t_max']]

            if drop_gaps:
                # make sure round errors do not affect this value
                n_vals_needed = int(round(max_times_in_day * daily_frac))
                # don't overwrite QC flag columns
                data_cols = [
                    c for c in df.columns if not c.endswith('_qc_flag')
                ]
                # interpolate first before counting gaps
                if max_interp_hours:
                    df[interp_vars] = interped[interp_vars].copy()
                    
                days_with_gaps = df[data_cols].groupby(
                    df.index.date).count() < n_vals_needed

            df = means.join(sums)

            # drop days based on fraction of missing subday samples
            if not downsample and drop_gaps:
                print(
                    'Filtering days with less then {}% or {}/{} sub-daily '
                    'measurements'.format(
                        daily_frac * 100, n_vals_needed, max_times_in_day
                    )
                )
                df[days_with_gaps] = np.nan
        else:
            self.n_samples_per_day = 1

        self._df = df.rename(columns=self.variables)
        return freq
    
    @property     
    def df(self):
        """
        See :attr:`fluxdataqaqc.Data.df` as the only difference is that the
        :attr:`QaQc.df` is first resampled to daily frequency.
        """
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
        monthly means or sums based on :attr:`QaQc.agg_dict`, provides data 
        as :obj:`pandas.DataFrame`. 
        
        Note that monthly means or sums are forced to null values if less than
        20 percent of a months days are missing in the daily data
        (:attr:`QaQc.df`). Also, for variables that are summed (e.g. ET or
        precipitation) missing days (if less than 20 percent of the month) will
        be filled with the month's daily mean value before summation.

        If a :obj:`QaQc` instance has not yet run an energy balance correction
        i.e. :attr:`QaQc.corrected` = :obj:`False` before accessing
        :attr:`monthly_df` then the default routine of data correction (energy
        balance ratio method) will be conducted.

        Utilize the :attr:`QaQc.monthly_df` property the same way as the 
        :attr:`fluxdataqaqc.Data.df`, see it's API documentation for examples.

        Tip:
            If you have additional variables in :attr:`QaQc.df` or would like
            to change the aggregation method for the monthly time series,
            adjust the instance attribute :attr:`QaQc.agg_dict` before 
            accessing the :attr:`QaQc.monthly_df`.

        """
        if not self.corrected and self._has_eb_vars:
            self.correct_data()

        # rename columns to internal names 
        df = self._df.rename(columns=self.inv_map).copy()
        # avoid including string QC flags because of float forcing on resample
        numeric_cols = [c for c in df.columns if not '_qc_flag' in c]
        sum_cols = [k for k,v in self.agg_dict.items() if v == 'sum']
        # to avoid warning/error of missing columns 
        sum_cols = list(set(sum_cols).intersection(df.columns))
        mean_cols = list(set(numeric_cols) - set(sum_cols))
        # if data type has changed to 'obj' resample skips... 
        # make sure data exists
        if len(mean_cols) >= 1:
            means = monthly_resample(df[mean_cols], mean_cols, 'mean', 0.8)
        else:
            means = None
        if len(sum_cols) >= 1:
            sums = monthly_resample(df[sum_cols], sum_cols, 'sum', 0.8)
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
        
        df.replace([np.inf, -np.inf], np.nan, inplace=True)

        return df.rename(columns=self.variables)

    def write(self, out_dir=None, use_input_names=False):
        """
        Save daily and monthly time series of initial and "corrected" data in 
        CSV format.

        Note, if the energy balance closure correction
        (:attr:`QaQc.correct_data`) has not been run, this method will run it
        with default options before saving time series files to disk.

        The default location for saving output time series files is within an
        "output" subdifrectory of the parent directory containing the
        config.ini file that was used to create the :obj:`fluxdataqaqc.Data`
        and :obj:`QaQc` objects, the names of the files will start with the
        site_id and have either the "daily_data" or "monthly_data" suffix. 

        Keyword Arguments:
            out_dir (str or :obj:`None`): default :obj:`None`. Directory to 
                save CSVs, if :obj:`None` save to :attr:`out_dir` instance 
                variable (typically "output" directory where config.ini file 
                exists).
            use_input_names (bool): default :obj:`False`. If :obj:`False` use 
                ``flux-data-qaqc`` variable names as in output file header,
                or if :obj:`True` use the user's input variable names where
                possible (for variables that were read in and not modified or
                calculated by ``flux-data-qaqc``).

        Returns:
            :obj:`None`

        Example:
            
            Starting from a config.ini file,

            >>> from fluxdataqaqc import Data, QaQc
            >>> d = Data('path/to/config.ini')
            >>> q = QaQc(d)
            >>> # note no energy balance closure correction has been run
            >>> q.corrected
                False
            >>> q.write()
            >>> q.corrected
                True

        Note:
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
            print(
                'WARNING: energy balance closure corrections have not yet been'
                'run. Using default options before writing output time series.'
            )
            self.correct_data()

        daily_outf = out_dir / '{}_daily_data.csv'.format(self.site_id)
        monthly_outf = out_dir / '{}_monthly_data.csv'.format(self.site_id)

        if use_input_names:
            self.df.to_csv(daily_outf)
            self.monthly_df.to_csv(monthly_outf)
        else:
            self.df.rename(columns=self.inv_map).to_csv(daily_outf)
            self.monthly_df.rename(columns=self.inv_map).to_csv(monthly_outf)

    @classmethod
    def from_dataframe(cls, df, site_id, elev_m, lat_dec_deg, var_dict,
            drop_gaps=True, daily_frac=1.00, max_interp_hours=2, 
            max_interp_hours_night=4):
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
        qaqc.temporal_freq = qaqc._check_daily_freq(
            drop_gaps, daily_frac, max_interp_hours, max_interp_hours_night
        )
        qaqc.corrected = False 
        qaqc.corr_meth = None
        qaqc._has_eb_vars = True

        return qaqc

    def correct_data(self, meth='ebr', et_gap_fill=True, y='Rn', refET='ETr',
            x=['G','LE','H'], fit_intercept=False):
        """
        Correct turblent fluxes to improve energy balance closure using an
        Energy Balance Ratio method modified from `FLUXNET
        <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__. 

        Currently three correction options are available: 'ebr' (Energy Balance
        Ratio), 'br' (Bowen Ratio), and 'lin_regress' (least squares linear
        regression). If you use one method followed by another corrected,the
        corrected versions of LE, H, ET, ebr, etc. will be overwritten with the
        most recently used approach. 
        
        This method also computes potential clear sky radiation 
        (saved as "rso") using an ASCE approach based on station elevation and 
        latitude. ET is calculated from raw and corrected LE using daily air
        temperature to correct the latent heat of vaporization, if air temp. is
        not available in the input data then air temp. is assumed at 20 
        degrees celcius.

        Corrected or otherwise newly calculated variables are named using the
        following suffixes to distinguish them:

        .. code-block:: text

            uncorrected LE, H, etc. from input data have no suffix
            _corr uses adjusted LE, H, etc. from the correction method used 
            _user_corr uses corrected LE, H, etc. found in data file (if provided)

        Arguments:
            y (str): name of dependent variable for regression, must be in 
                :attr:`QaQc.variables` keys, or a user-added variable.
                Only used if ``meth='lin_regress'``.
            x (str or list): name or list of independent variables for 
                regression, names must be in :attr:`QaQc.variables` keys, or a 
                user-added variable. Only used if ``meth='lin_regress'``.

        Keyword Arguments:
            meth (str): default 'ebr'. Method to correct energy balance.
            et_gap_fill (bool): default True. If true fill any remaining gaps
                in corrected ET with ETr * ETrF, where ETr is alfalfa reference
                ET from gridMET and ETrF is the filtered, smoothed (7 day
                moving avg.  min 2 days) and linearly interpolated crop
                coefficient. The number of days in each month that corrected ET
                are filled will is provided in :attr:`QaQc.monthly_df` as the
                column "ET_gap". 
            refET (str): default "ETr". Which gridMET reference product to use
                for ET gap filling, "ETr" or "ETo" are valid options.
            fit_intercept (bool): default False. Fit intercept for regression or
                set to zero if False. Only used if ``meth='lin_regress'``.
            apply_coefs (bool): default False. If :obj:`True` then apply fitted
                coefficients to their respective variables for linear regression
                correction method, rename the variables with the suffix "_corr".
        
        Returns:
            :obj:`None`

        Example:
            Starting from a correctly formatted config.ini and climate time
            series file, this example shows how to read in the data and apply
            the energy balance ratio correction without gap-filling with 
            gridMET ETr x ETrF.

            >>> from fluxdataqaqc import Data, QaQc
            >>> d = Data('path/to/config.ini')
            >>> q = QaQc(d)
            >>> q.corrected
                False

            Now apply the energy balance closure correction 

            >>> q.correct_data(meth='ebr', et_gap_fill=False)
            >>> q.corrected
                True

        Note:
            If ``et_gap_fill`` is set to :obj:`True` (default) the gap filled
            days of corrected ET will be used to recalculate LE_corr for those
            days with the gap filled values, i.e. LE_corr will also be
            gap-filled.

        Note:
            The *ebr_corr* variable or energy balance closure ratio is 
            calculated from the corrected versions of LE and H independent 
            of the method. When using the 'ebr' method the energy balance 
            correction factor (what is applied to the raw H and LE) is left as 
            calculated (inverse of ebr) and saved as *ebc_cf*. 

        See Also:
            For explanation of the linear regression method see the
            :meth:`QaQc.lin_regress` method, calling that method with the
            keyword argument ``apply_coefs=True`` and :math:`Rn` as the y
            variable and the other energy balance components as the x variables
            will give the same result as the default inputs to this function
            when ``meth='lin_regress``.  
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
        eb_vars = {'Rn','LE','H','G'}
        if not eb_vars.issubset(self.variables.keys()) or\
                not eb_vars.issubset(
                        self.df.rename(columns=self.inv_map).columns):
            print(
                '\nMissing one or more energy balance variables, cannot perform'
                ' energy balance correction.'
            )
            self._has_eb_vars = False
            # calculate raw (from input LE) ET 
            self._calc_et()
            # fill gaps of raw ET with ET from smoothed and gap filled ETrF*ETr
            if et_gap_fill:
                self._ET_gap_fill(et_name='ET', refET=refET)
            return

        if meth == 'ebr':
            self._ebr_correction()
        elif meth == 'br':
            self._bowen_ratio_correction()
        elif meth == 'lin_regress':
            if y not in eb_vars or len(set(x).difference(eb_vars)) > 0:
                print(
                    'WARNING: correcting energy balance variables using '
                    'depedendent or independent variables that are not '
                    'energy balance components will cause undesired results!'
                    '\nIt is recommended to use only Rn, G, LE, and H for '
                    'dependent (y) and independent (x) variables.'
                )
            self.lin_regress(
                y=y, x=x, fit_intercept=fit_intercept, apply_coefs=True
            )

        self.corr_meth = meth
        # calculate raw, corrected ET 
        self._calc_et()
        # fill gaps of corr ET with ET from smoothed and gap filled ETrF*ETr
        if et_gap_fill:
            self._ET_gap_fill(et_name='ET_corr', refET=refET)

        # update inv map for naming
        self.inv_map = {
            v: k for k, v in self.variables.items() if (
                not v.replace('_mean', '') == k or not k in self.df.columns)
        }
        # using 'G' in multiple g plot may overwrite G name internally
        if not 'G' in self.inv_map.values():
            user_G_name = self.variables.get('G')
            self.inv_map[user_G_name] = 'G'

    def _calc_vpd_from_vp(self):
        """
        Based on ASCE standardized ref et eqn. 37, air temperature must be in 
        celcius and actual vapor pressure in kPa.

        Can also calculate VP from VPD and air temperature, es and rel.
        humidity.
        """
        df = self.df.rename(columns=self.inv_map)
        # calculate vpd from actual vapor pressure and temp
        # check if needed variables exist and units are correct
        has_vpd_vars = set(['vp','t_avg']).issubset(df.columns)
        units_correct = (
            self.units.get('vp') == 'kpa' and self.units.get('t_avg') == 'c'
        )
        if has_vpd_vars and units_correct:
            # saturation vapor pressure (es)
            es = 0.6108 * np.exp(17.27 * df.t_avg / (df.t_avg + 237.3))
            df['vpd'] = es - df.vp
            df['es'] = es
            self.variables['vpd'] = 'vpd'
            self.units['vpd'] = 'kpa'
            self.variables['es'] = 'es'
            self.units['es'] = 'kpa'

        # same calc actual vapor pressure from vapor pressure deficit and temp
        has_vp_vars = set(['vpd','t_avg']).issubset(df.columns)
        units_correct = (
            self.units.get('vpd') == 'kpa' and self.units.get('t_avg') == 'c'
        )

        if has_vp_vars and units_correct:
            # saturation vapor pressure (es)
            es = 0.6108 * np.exp(17.27 * df.t_avg / (df.t_avg + 237.3))
            df['vp'] = es - df.vpd
            df['es'] = es
            self.variables['vp'] = 'vp'
            self.units['vp'] = 'kpa'
            self.variables['es'] = 'es'
            self.units['es'] = 'kpa'

        if not 'rh' in self.variables and {'vp','es'}.issubset(self.variables):
            if not self.units.get('vp') == 'kpa': pass
            else:
                print(
                    'Calculating relative humidity from actual and saturation '
                    'vapor pressure and air temperature'
                )
                df['rh'] = 100 * (df.vp / df.es)
                self.variables['rh'] = 'rh'
                self.units['rh'] = '%'

        self._df = df

    def _ET_gap_fill(self, et_name='ET_corr', refET='ETr'):
        """
        Use gridMET reference ET to calculate ET from ETrF/EToF, smooth and gap
        fill calced ET and then use to fill gaps in corrected ET. Keeps tabs on
        number of days in ET that were filled in each month.
        
        Keyword Arguments:
            et_name (str): default "ET_corr". Name of ET variable to use when 
                calculating the crop coefficient and to fill gaps with 
                calculated ET.

        Returns:
            :obj:`None`

        """
        # drop relavant calculated variables if they exist
        self._df = _drop_cols(self._df, self._et_gap_fill_vars)
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

        df = self.df.rename(columns=self.inv_map)
        if not et_name in df.columns:
            print(
                'ERROR: {} not found in data, cannot gap-fill'.format(et_name)
            )
            return

        # calc ETrF 7 day moving average, min vals in window = 2, linear interp
        print(
            f'Gap filling {et_name} with filtered {refET}F x {refET} (gridMET)'
        )
        if refET == 'ETr':
            df['ETrF'] = df[et_name].astype(float)/df.gridMET_ETr.astype(float)
            df['ETrF_filtered'] = df['ETrF']
            # filter out extremes of ETrF
            Q1 = df['ETrF_filtered'].quantile(0.25)
            Q3 = df['ETrF_filtered'].quantile(0.75)
            IQR = Q3 - Q1
            to_filter = df.query(
                'ETrF_filtered<(@Q1-1.5*@IQR) or ETrF_filtered>(@Q3+1.5*@IQR)'
            )
            df.loc[to_filter.index, 'ETrF_filtered'] = np.nan
            df['ETrF_filtered'] = df.ETrF_filtered.rolling(
                7, min_periods=2, center=True
            ).mean()
            df.ETrF_filtered = df.ETrF_filtered.interpolate(method='linear')
            # calc ET from ETrF_filtered and ETr
            df['ET_fill'] = df.gridMET_ETr * df.ETrF_filtered
            # flag days in corrected ET that are missing and sum for each month
            df['ET_gap'] = False
            df.loc[(df[et_name].isna() & df.ET_fill.notna()), 'ET_gap'] = True
            df.loc[df.ET_gap, et_name] = df.loc[df.ET_gap, 'ET_fill']
            self.units['ETrF'] = 'unitless'
            self.units['ETrF_filtered'] = 'unitless'
        elif refET == 'ETo':
            df['EToF'] = df[et_name].astype(float)/df.gridMET_ETo.astype(float)
            df['EToF_filtered'] = df['EToF']
            # filter out extremes of EToF
            Q1 = df['EToF_filtered'].quantile(0.25)
            Q3 = df['EToF_filtered'].quantile(0.75)
            IQR = Q3 - Q1
            to_filter = df.query(
                'EToF_filtered<(@Q1-1.5*@IQR) or EToF_filtered>(@Q3+1.5*@IQR)'
            )
            df.loc[to_filter.index, 'EToF_filtered'] = np.nan
            df['EToF_filtered'] = df.EToF_filtered.rolling(
                7, min_periods=2, center=True
            ).mean()
            df.EToF_filtered = df.EToF_filtered.interpolate(method='linear')
            # calc ET from EToF_filtered and ETo
            df['ET_fill'] = df.gridMET_ETo * df.EToF_filtered
            # flag days in corrected ET that are missing and sum for each month
            df['ET_gap'] = False
            df.loc[(df[et_name].isna() & df.ET_fill.notna()), 'ET_gap'] = True
            df.loc[df.ET_gap, et_name] = df.loc[df.ET_gap, 'ET_fill']
            self.units['EToF'] = 'unitless'
            self.units['EToF_filtered'] = 'unitless'
            

        # backcalculate LE_corr and flux_corr with gap filled et_corr
        if et_name == 'ET_corr' and 't_avg' in self.variables:
            df['LE_corr'] =\
                (df.ET_corr * (2501000 - 2361 * df.t_avg.fillna(20))) / 86400
        # assume 20 degrees C
        elif et_name == 'ET_corr' and not 't_avg' in self.variables:
            df['LE_corr'] =\
                (df.ET_corr * (2501000 - 2361 * 20)) / 86400

        # et fill values only where they were used to gap fill, for plotting
        df['ET_fill_val'] = np.nan
        df.loc[df.ET_gap , 'ET_fill_val'] = df.ET_fill

        # update variables attribute with new variables (may vary)
        new_cols = set(df.columns) - set(self.variables)
        for el in new_cols:
            self.variables[el] = el
        self.units['ET_fill'] = 'mm'
        self.units['ET_fill_val'] = 'mm'
        self.units['ET_gap'] = 'unitless'

        self._df = df


    def _calc_et(self):
        """
        Calculate daily ET (mm) from raw and corrected LE (w/m2)), if air 
        temperature is available use to correct latent heat of vaporization.
        
        Currently computes on :attr:`df` dataframe attribute in place. Can be 
        called before or after energy balance closure corrections, if before 
        the only raw ET will be calculated. 

        Arguments:
            None

        Returns:
            :obj:`None`
        """

        # drop relavant calculated variables if they exist
        self._df = _drop_cols(self._df, self._et_calc_vars)
        df = self.df.rename(columns=self.inv_map)
        if not set(['LE','LE_corr','LE_user_corr']).intersection(df.columns):
            print(
                'ERROR: no LE variables found in data, cannot calculate ET'
            )
            return
        
        # LH from L.P. Harrison (1963)
        if 't_avg' in df.columns:
            df['ET'] = 86400 * (df.LE/(2501000 - (2361 * df.t_avg.fillna(20))))
            if 'LE_corr' in df.columns:
                df['ET_corr']=\
                    86400 * (df.LE_corr/(2501000 - (2361*df.t_avg.fillna(20))))
            if 'LE_user_corr' in df.columns:
                df['ET_user_corr'] =\
                    86400*(df.LE_user_corr/(2501000-(2361*df.t_avg.fillna(20))))
        # otherwise assume air temp = 20 degrees C
        else:
            df['ET'] = 86400 * (df.LE/(2501000 - (2361 * 20)))
            if 'LE_corr' in df.columns:
                df['ET_corr'] = 86400 * (df.LE_corr/(2501000 - (2361 * 20)))
            if 'LE_user_corr' in df.columns:
                df['ET_user_corr']=86400*(df.LE_user_corr/(2501000-(2361*20)))
        
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
            :obj:`None`

        Note:
            Does not overwrite existing calculation for Rso, i.e. once any 
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

    def lin_regress(self, y, x, fit_intercept=False, apply_coefs=False):
        """
        Least squares linear regression on single or multiple independent 
        variables.

        For example, if the dependent variable (``y``) is :math:`Rn` and the
        independent variables (``x``) are :math:`LE` and :math:`H`, then the
        linear regression will solve for the best fit coefficients of :math:`Rn
        = c_0 + c_1 LE + c_2 H`. Any number of variables in the
        :attr:`QaQc.variables` can be used for ``x`` and one for ``y``. 
        
        If the variables chosen for regression are part of the energy balance
        components, i.e. :math:`Rn, G, LE, H` and ``apply_coefs=True``, then
        the best fit coefficients will be applied to their respecitive
        variables with consideration of the energy balance equation, i.e. the
        signs of the coefficients will be corrected according to :math:`Rn - G
        = LE + H`, for example if ``y=H`` and ``x=['Rn','G','LE]``. i.e.
        solving :math:`H = c_0 + c_1 Rn + c_2 G + c_3 LE` then the coefficients
        :math:`c_2` and :math:`c_3` will be multiplied by -1 before applying
        them to correct :math:`G` and :math:`LE` according to the energy
        balance equation. 

        This method returns an :obj:`pandas.DataFrame` object containing
        results of the linear regression including the coefficient values,
        number of data pairs used in the, the root-mean-square-error, and
        coefficient of determination. This table can also be retrieved from the
        :attr:`QaQc.lin_regress_results` instance attribute.

        Arguments:
            y (str): name of dependent variable for regression, must be in 
                :attr:`QaQc.variables` keys, or a user-added variable.
            x (str or list): name or list of independent variables for 
                regression, names must be in :attr:`QaQc.variables` keys, or a 
                user-added variable.

        Keyword Arguments:
            fit_intercept (bool): default False. Fit intercept for regression or
                set to zero if False.
            apply_coefs (bool): default False. If :obj:`True` then apply fitted
                coefficients to their respective variables, rename the variables
                with the suffix "_corr".

        Returns:
            :obj:`pandas.DataFrame` 

        Example:
            Let's say we wanted to compute the linear relationship between net
            radiation to the other energy balance components which may be
            useful if we have strong confidence in net radiation measurements
            for example. The resulting coefficients of regression would give us
            an idea of whether the other components were "under-measured" or
            "over-measured". Then, starting with a :obj:`.Data` instance:

            >>> Q = QaQc(Data_instance)
            >>> Q.lin_regress(y='Rn', x=['G','H','LE'], fit_intercept=True)
            >>> Q.lin_regress_result

            This would produce something like the following,

======= ================== ================ ============== =============== ============== =========== =============== ================ 
SITE_ID Y (dependent var.) c0 (intercept)   c1 (coef on G) c2 (coef on LE) c3 (coef on H) RMSE (w/m2) r2 (coef. det.) n (sample count) 
======= ================== ================ ============== =============== ============== =========== =============== ================ 
a_site  Rn                 6.99350781229883 1.552          1.054           0.943          18.25       0.91            3386             
======= ================== ================ ============== =============== ============== =========== =============== ================

            In this case the intercept is telling us that there may be a systematic
            or constant error in the independent variables and that :math:`G` is
            "under-measured" at the sensor by over 50 precent, etc. if we assume
            daily :math:`Rn` is accurate as measured.

        Tip:
            You may also use multiple linear regression to correct energy
            balance components using the :meth:`QaQc.correct_data` method
            by passing the ``meth='lin_regress'`` keyword argument. 

        """
        # drop relavant calculated variables if they exist
        self._df = _drop_cols(self.df, self._eb_calc_vars)
        df = self._df.rename(columns=self.inv_map)
        if not y in df.columns:
            print(
                'ERROR: the dependent variable ({}) was not '
                'found in the dataframe.\nHere are all available variable '
                'names:\n{}\n'.format(y, ', '.join(df.columns))
            )
            return
        if not isinstance(x, list) and not x in df.columns:
            print(
                'ERROR: the dependent variable ({}) was not '
                'found in the dataframe.\nHere are all available variable '
                'names:\n{}\n'.format(x, ', '.join(df.columns))
            )
            return

        # get n X vars, names if more than 1
        n_x = 1
        if isinstance(x, list):
            n_x = len(x)
            if not set(x).issubset(df.columns):
                print(
                    'ERROR: one or more independent variables ({}) were not '
                    'found in the dataframe.\nHere are all available variable '
                    'names:\n{}'.format(','.join(x), ', '.join(df.columns))
                )
                return
            if n_x > 1:
                cols = x + [y] 
                tmp = df[cols].copy()
        if n_x == 1:
            tmp = df[[x,y]].copy()
        # drop timestamps with any missing variables
        tmp = tmp.dropna()
        X = tmp[x]
        Y = tmp[y]
        # create model, fit, and predict Y for RMSE/R2 calcs 
        model = linear_model.LinearRegression(fit_intercept=fit_intercept)
        model.fit(X, Y)
        pred = model.predict(X)
        r2 = model.score(X,Y).round(2)
        rmse = (np.sqrt(mean_squared_error(Y, pred))).round(2)

        eb_vars = ['LE','H','Rn','G']
        if apply_coefs and set(tmp.columns).intersection(eb_vars):
            # calc initial energy balance if regression applied to EB vars
            df['ebr'] = (df.H + df.LE) / (df.Rn - df.G)
            df['flux'] = df.H + df.LE
            df['energy'] = df.Rn - df.G

        # dataframe of results in nice/readable format
        results = pd.Series()
        results.loc['Y (dependent var.)'] = y
        results.loc['c0 (intercept)'] = model.intercept_
        for i, c in enumerate(X.columns):
            results.loc['c{} (coef on {})'.format(i+1, c)] =\
                model.coef_[i].round(3)
            if apply_coefs:
                new_var = '{}_corr'.format(c)
                # ensure correct sign of coef. if applied to EB vars
                if y in ['G','H','LE'] and c in ['G','H','LE']:
                    coef = -1 * model.coef_[i]
                else:
                    coef = model.coef_[i]
                print(
                    'Applying correction factor ({}) to '
                    'variable: {} (renamed as {})'.format(
                        coef.round(3), c, new_var
                    )
                )

                df[new_var] = df[c] * coef
                self.variables[new_var] = new_var
                self.units[new_var] = self.units.get(c) # already have units

        results.loc[
            'RMSE ({})'.format(self.units.get(y,'na'))
        ] = rmse
        results.loc['r2 (coef. det.)'] = r2
        results.loc['n (sample count)'] = len(Y)
        results = results.to_frame().T
        results.index= [self.site_id]
        results.index.name = 'SITE_ID'
        # add as instance variable
        self.lin_regress_results = results

        # depending on independent vars, calc. corrected flux, ebr, energy
        if apply_coefs and set(X.columns).intersection(eb_vars):
            corr = pd.DataFrame()
            for v in eb_vars:
                cv = '{}_corr'.format(v)
                if cv in df:
                    corr[v] = df[cv] 
                else:
                    corr[v] = df[v]

            # not all vars are necessarily different than initial
            df['ebr_corr'] = (corr.H + corr.LE) / (corr.Rn - corr.G)
            df['flux_corr'] = corr.H + corr.LE
            df['energy_corr'] = corr.Rn + corr.G
            del corr

            self.variables.update(
                energy = 'energy',
                flux = 'flux',
                flux_corr = 'flux_corr',
                energy_corr = 'energy_corr',
                ebr = 'ebr',
                ebr_corr = 'ebr_corr'
            )
            self.corrected = True
            self.corr_meth = 'lin_regress'

        self._df = df.rename(columns=self.variables)

        return results

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
            :obj:`None`

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
                if abs(1/val) >= 2 or abs(1/val) <= 0.5:
                    val = np.nan
            # if at least one day exists in window2 take mean
            elif np.count_nonzero(~np.isnan(win_arr2)) > 0:
                val = np.nanmedian(win_arr2)
                if abs(1/val) >= 2 or abs(1/val) <= 0.5:
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
            df, ebr_5day_clim, left_on='DOY', right_index=True
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
        # filter out CF that are >=2 or <=0.5 (absolute)
        merged.loc[
            (abs(merged.ebc_cf) >= 2) | (abs(merged.ebc_cf <= 0.5)), 'ebc_cf'
        ] = np.nan
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
            :obj:`None`
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

    def plot(self, ncols=1, output_type='save', out_file=None, suptitle='', 
            plot_width=1000, plot_height=450, sizing_mode='scale_both', 
            merge_tools=False, link_x=True, **kwargs):
        """
        Creates a series of interactive diagnostic line and scatter 
        plots of input and computed daily and monthly aggregated data.

        The main interactive features of the plots include: pan, selection and
        scrol zoom, hover tool that shows paired variable values including date,
        and linked x-axes that pan/zoom togehter for daily and monthly time 
        series plots.
 
        It is possible to change the format of the output plots including
        adjusting the dimensions of subplots, defining the number of columns of
        subplots, setting a super title that accepts HTML, and other options.
        If variables are not present for plots they will not be created and a
        warning message will be printed. There are two options for output: open
        a temporary file for viewing or saving a copy to :attr:`QaQc.out_dir`.
        
        A list of all potential time series plots created: 
        
        * energy balance components 
        * radiation components 
        * incoming shortwave radiation with ASCE potential clear sky (daily only)
        * multiple soil heat flux measurements
        * air temperature
        * vapor pressure and vapor pressure deficit
        * wind speed
        * station precipitation and gridMET precipitation
        * initial and corrected latent energy
        * initial, corrected, gap filled, and reference evapotranspiration
        * crop coefficient and smoothed and interpolated crop coefficient
        * initial and corrected energy balance ratio
        * multiple soil moisture measurements

        A list of all potential scatter plots created: 

        * radiative versus turblent fluxes, initial and corrected
        * initial versus corrected latent energy
        * initial versus corrected evapotranspiration

        Keyword Arguments:
            ncols (int): default 1. Number of columns of subplots.
            output_type (str): default "save". How to output plots, "save", 
                "show" in browser, "notebook" for Jupyter Notebooks, 
                "return_figs" to return a list of Bokeh 
                :obj:`bokeh.plotting.figure.Figure`s, or "return_grid" to 
                return the :obj:`bokeh.layouts.gridplot`.
            out_file (str or None): default :obj:`None`. Path to save output 
                file, if :obj:`None` save output to :attr:`QaQc.out_dir` with 
                the name [site_id]_plots.html where [site_id] is 
                :attr:`QaQc.site_id`.
            suptitle (str or None): default :obj:`None`. Super title to go 
                above plots, accepts HTML/CSS syntax.
            plot_width (int): default 1000. Width of subplots in pixels.
            plot_height (int): default 450. Height of subplots in pixels, note 
                for subplots the height will be forced as the same as 
                ``plot_width``.
            sizing_mode (str): default "scale_both". Bokeh option to scale
                dimensions of :obj:`bokeh.layouts.gridplot`.
            merge_tools (bool): default False. Merges all subplots toolbars into
                a single location if True.
            link_x (bool): default True. If True link x axes of daily time 
                series plots and monthly time series plots so that when zooming 
                or panning on one plot they all zoom accordingly, the axes will 
                also be of the same length.

        Example:
            Starting from a correctly formatted config.ini and climate time
            series file, this example shows how to read in the data and produce
            the default series of plots for viewing with the addition of text
            at the top of plot that states the site's location and ID.

            >>> from fluxdataqaqc import Data, QaQc
            >>> d = Data('path/to/config.ini')
            >>> q = QaQc(d)
            >>> q.correct_data()
            >>> # create plot title from site ID and location in N. America
            >>> title = "<b>Site:</b> {}; <b>Lat:</b> {}N; <b>Long:</b> {}W".format(
            >>>     q.site_id, q.latitude, q.longitude
            >>> )
            >>> q.plot(
            >>>     ncols=2, output_type='show', plot_width=500, suptitle=title
            >>> )

            Note, we specified the width of plots to be smaller than default
            because we want both columns of subplots to be viewable on one page.

        Tip:
            To reset all subplots at once, refresh the page with your web
            browser.

        Note:
            Additional keyword arguments that are recognized by
            :obj:`bokeh.layouts.gridplot` are also accepted by
            :meth:`QaQc.plot`.

        """

        # handle file save for accessing from instance variable
        if out_file is None and output_type == 'save':
            out_file = Path(self.out_dir)/'{}_plots.html'.format(self.site_id)
            out_dir = out_file.parent
            if not out_dir.is_dir():
                out_dir.mkdir(parents=True, exist_ok=True)
        # if out_file is to a non-existent directory create parents
        elif out_file is not None and output_type == 'save':
            out_dir = Path(out_file).parent
            if not out_dir.is_dir():
                out_dir.mkdir(parents=True, exist_ok=True)
        
        # create aggregrated plot structure from fluxdataqaqc.Plot._plot() 
        ret = self._plot(
            self, ncols=ncols, output_type=output_type, out_file=out_file,
            suptitle=suptitle, plot_width=plot_width, plot_height=plot_height,
            sizing_mode=sizing_mode, merge_tools=merge_tools, link_x=link_x,
            **kwargs
        )

        self.plot_file = out_file

        if ret:
            return ret

def _drop_cols(df, cols):
    """Drop columns from dataframe if they exist """
    for c in cols:
        if c in df.columns:
            df.drop(c, axis=1, inplace=True)

    return df

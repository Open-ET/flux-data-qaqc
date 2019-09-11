.. There are two major functionalities in
   ``flux-data-qaqc``, first, correcting surface energy balance by
   adjusting latent energy and sensible heat fluxes and calculate other
   climatic variables. Second, it serves as a robust way to read in
   different time series data and produce visualizations, e.g. their daily
   and monthly time series.

Configuration Options and Caveats
=================================

This tutorial shows how to use ``flux-data-qaqc`` with climate data of various
formats and generally covers formatting rules of input data and extra options
that can be set in a config file. The major differences when using
``flux-data-qaqc`` for different input data lie in the config file declarations
therefore the entire workflow from :ref:`Basic Usage`
will work just the same once your config file is set up correctly. 

Example data description
------------------------

The data used in :ref:`Configuration Options and Caveats`  is provided with
``flux-data-qaqc`` and can be downloaded `here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Config_options>`__.
There are two datasets used in the following examples, the data used for
showing quality based filtering of data based on QC flags is from the FLUXNET
2015 dataset for site "ARM USDA UNL OSU Woodward Switchgrass 1" which contains
switchgrass, more information on this site can be found `here <http://sites.fluxdata.org/US-AR1/>`__. The site that contains multiple soil
heat flux measurements for weighted averaging and soil moisture measurements
for plotting is a subset (shortened for reduced disk space) of the "ARM Southern Great Plains site- Lamont" AmeriFlux site dataset, more information on this site can be found `here <http://ameriflux.lbl.gov/sites/siteinfo/US-ARM>`__.

A reproducible Jupyter Notebook with minor differences of this tutorial can be
found `here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Config_options/advanced_config_options.ipynb>`__.


Seting up a config file
-----------------------

``flux-data-qaqc`` starts with the creation of a configuration file, a
text file with extension ".ini" or ".INI" that follows the rules set 
`here <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`__.
A config file for ``flux-data-qaqc`` requires the sections: 1. **METADATA** 2. **DATA**
although you may provide additional for custom uses. 

In **METADATA** you may enter any metadata that you wish so long
as the entrie's *key* is followed by an equal sign and the assigned 
*value*, i.e. 

.. parsed-literal::

    key = value

Here are the *mandatory* metadata entries unique to ``flux-data-qaqc``:

* climate_file_path
* station_elevation
* station_longitude
* station_latitude
* site_id

The “climate_file_path” is the full or relative file path of the input climate
file (excel or CSV, more on formatting this below) containing the climatic data
to be analyzed. The “station_elevation” (meters) and latitude/longitude
(decimal degrees) fields are used to calculate clear sky potential solar
radiation using an ASCE formulation and to locate the nearest gridMET cell
centroid when downloading reference ET. “site_id” is used for saving output
files. 

Optional metadata entries that are used by ``flux-data-qaqc`` include
“missing_data_value”, “qc_threshold”, “qc_flag”, "var_name_delim", "skiprows",
"date_parser", and "gridmet_file_path".  “missing_data_value” is used to
correctly parse missing values in the input climate time series. All other
optional metadata that can be used by ``flux-data-qaqc`` except
"gridmet_file_path" (which is simply the path to a file that is downloaded by
:meth:`.QaQc.download_gridMET`) are explained within this page.

The **DATA** section of the config file is where you specify climate
variables and their units following the same approach explained above. 

Here is a list of all the “expected” climate variable names in the
**DATA** section of a config file, where keys are the keys in the config 
file and values are the internal names used by ``flux-data-qaqc``. This list 
can be accessed at any time from :attr:`.Data.variable_names_dict`:

   >>> from fluxdataqaqc import Data
   >>> Data.variable_names_dict
       {'datestring_col' : 'date' ,
       'net_radiation_col' : 'Rn' ,
       'ground_flux_col' : 'G' ,
       'latent_heat_flux_col' : 'LE' ,
       'latent_heat_flux_corrected_col' : 'LE_user_corr' ,
       'sensible_heat_flux_col' : 'H' ,
       'sensible_heat_flux_corrected_col' : 'H_user_corr' ,
       'shortwave_in_col' : 'sw_in' ,
       'shortwave_out_col' : 'sw_out' ,
       'shortwave_pot_col' : 'sw_pot' ,
       'longwave_in_col' : 'lw_in' ,
       'longwave_out_col' : 'lw_out' ,
       'vap_press_col' : 'vp' ,
       'vap_press_def_col' : 'vpd' ,
       'avg_temp_col' : 't_avg' ,
       'precip_col' : 'ppt' ,
       'wind_spd_col' : 'ws'}

You may view these climate entry keys and values (as found in the config file)
from within Python using the :attr:`.Data.config` property which contains all
information listed in the config file as a :obj:`configparser.ConfigParser`
instance.

    >>> config_path = 'config_for_QC_flag_filtering.ini'
    >>> d = Data(config_path)
    >>> # loop through a list of tuples with keys and values from DATA section
    >>> for each in d.config.items('DATA'):
    >>>     print(each) 
        ('datestring_col', 'date')
        ('net_radiation_col', 'Rn')
        ('net_radiation_units', 'w/m2')
        ('ground_flux_col', 'G')
        ('ground_flux_units', 'w/m2')
        ('latent_heat_flux_col', 'LE')
        ('latent_heat_flux_qc', 'a_qc_value')
        ('latent_heat_flux_units', 'w/m2')
        ('latent_heat_flux_corrected_col', 'LE_corrected')
        ('latent_heat_flux_corrected_units', 'w/m2')
        ('sensible_heat_flux_col', 'H')
        ('sensible_heat_flux_qc', 'a_qc_value')
        ('sensible_heat_flux_units', 'w/m2')
        ('sensible_heat_flux_corrected_col', 'H_corrected')
        ('sensible_heat_flux_corrected_units', 'w/m2')
        ('shortwave_in_col', 'sw_in')
        ('shortwave_in_qc', 'swrad_flag')
        ('shortwave_in_units', 'w/m2')
        ('shortwave_out_col', 'sw_out')
        ('shortwave_out_units', 'w/m2')
        ('shortwave_pot_col', 'sw_pot')
        ('shortwave_pot_units', 'w/m2')
        ('longwave_in_col', 'lw_in')
        ('longwave_in_units', 'w/m2')
        ('longwave_out_col', 'lw_out')
        ('longwave_out_units', 'w/m2')
        ('vap_press_col', 'na')
        ('vap_press_units', 'na')
        ('vap_press_def_col', 'vpd')
        ('vap_press_def_units', 'hPa')
        ('avg_temp_col', 't_avg')
        ('avg_temp_units', 'C')
        ('precip_col', 'ppt')
        ('precip_units', 'mm')
        ('wind_spd_col', 'ws')
        ('wind_spd_units', 'm/s')

You can also access the data from the :attr:`.Data.config` as a dictionary,
for example if your **METADATA** section has an entry for "land_cover", e.g.

.. parsed-literal::
    
    [METADATA]
    land_cover = CROP
    ...

then access this value with :meth:`configparser.ConfigParser.get` which returns 
the value of "land_cover" in the config file's **METADATA** section

   >>> d.config.get('METADATA', 'land_cover')
       CROP

.. tip::
   If you are unsure if your config file's metadata contains a specific entry
   you can pass the ``fallback`` keyword-only argument to the
   :meth:`configparser.ConfigParser.get` method similar to a Python dictionary.

Here is an example,

   >>> d.config.get('METADATA', 'land_cov', fallback='not given')
       "not given"


Input formatting and caveats
----------------------------

Missing data
^^^^^^^^^^^^

For parsing data gaps in input time series assign the
“missing_data_value” to the **METADATA** section of the config file. 
The value should be numeric, e.g.  

.. parsed-literal::

    missing_data_value = -999

If the input time series file does not contain all climate variables that are
expeced by ``flux-data-qaqc``, then specify them as missing (‘na’) in the
config file or simply do not list them in the config. Missing variables will be
ignored for the most part and will not be present in output files/plots,
however if key variables for the energy balance are not present (:math:`LE`,
:math:`H`, :math:`G`, and :math:`Rn`) then you will not be able to run energy
balance closure correction routines.

Data file format
^^^^^^^^^^^^^^^^

``flux-data-qaqc`` accepts Microsoft Excel files (.xlx and .xlsx) and
comma separated value (CSV) text files containing time series input. 
The input file should have a column with combined date and time. Currently there is no
restriction on the temporal frequency of input data however it is
automatically resampled to daily frequency before running correction
routines. Lastly, there should be a single header row containing all
variable names followed by the first entry of climatic variables.

Here is an example of a valid input file’s first 5 rows and 8 columns:

========== ====== ======= ======= ======= ===== === =====
date       t_avg  sw_pot  sw_in   lw_in   vpd   ppt ws
========== ====== ======= ======= ======= ===== === =====
2009-01-01 2.803  186.71  123.108 261.302 1.919 0   3.143
2009-01-02 2.518  187.329 121.842 268.946 0.992 0   2.093
2009-01-03 5.518  188.008 124.241 268.004 2.795 0   4.403
2009-01-04 -3.753 188.742 113.793 246.675 0.892 0   4.336
========== ====== ======= ======= ======= ===== === =====

.. note:: 
   If the the input datas temporal frequency is not recognized
   ``flux-data-qaqc`` will attempt to resample it to daily frequency when it is
   used to create a :obj:`.QaQc` object. Also, if a value is not recognized a
   numeric in any data column it will be forced to a null value.

Data header formatting
^^^^^^^^^^^^^^^^^^^^^^

A common format of some time series data is that the header row may
not start on the first line of the file. If this is the case you must add
an entry to the **METADATA** section of the config file "skiprows" which
states the number of rows to skip before finding the header row. A 
caveat is that if using CSV data files you may have any number of comment
lines before the header so long as they start with a hashtag symbol "#"
(comment), in this case you should not add "skiprows" to **METADATA**. 

Optimize data load time 
^^^^^^^^^^^^^^^^^^^^^^^

``flux-data-qaqc`` utilizes the :mod:`pandas` for most time series data
management, specifically the usage of :obj:`datetime.datetime` objects for
advanced temporal analysis tools. If your file is large you can specify the 
datetime format in the **METADATA** section of the config file to potentially
greatly speedup the loading of data. For example if your date column contains
strings in the format year month day hour minute with no delimiters, e.g. 
201401010000 for 2014 January 1st at midnight, then in the ``flux-data-qaqc``
config file you would enter:

.. parsed-literal::

    date_parser = %Y%m%d%H%M

For more information of the correct date parser string for your date format
see the directives of the :meth:`datetime.datetime.strptime` `here <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior>`__.



--------------

Quality-based data filtering 
----------------------------

Currently ``flux-data-qaqc`` supports filtering out poor quality data
based on user-provided quality control (QC) values (numeric) or flags 
(characters) using the :meth:`.Data.apply_qc_flags` method. This feature 
helps to facilitate manual or semi-manual data filtering which is 
sometimes necessary during data preprocessing.

Flag-based filtering
^^^^^^^^^^^^^^^^^^^^

Let’s say that you have a column in your input data named ‘QC_flag’ that
contains character strings signifying the assigned data quality for a
climate time series. The flag is either ‘g’ meaning a data point is ‘good’
or if the flag is ‘b’ the data point is bad quality and you would like
to filter it. Further let's say that you want the filter to apply to your
latent energy and and sensible heat variables, then in your config file 
you would need to declare the flag for 'bad' data ('b') to be filtered out
in the **METADATA** section:

.. code:: bash

   qc_flag = b

and in the **DATA** section of your config you will state that the
‘QC_flag’ column should be applied to your :math:`LE` and :math:`H` variables:

.. code:: bash

   latent_heat_flux_qc = QC_flag
   sensible_heat_flux_qc = QC_flag

Now, when the :meth:`.Data.apply_qc_flags` method is used the all date entries
of :math:`LE` and :math:`H` that have a "QC_flag" value of 'b' will be forced
to null in the :attr:`.Data.df` property of a :obj:`.Data` instance. 

Threshold-based filtering
^^^^^^^^^^^^^^^^^^^^^^^^^

Another option is to use a numeric quality control *value* that exists
in your input data along with a threshold value which means that when
the quality control value falls below this threshold you would like to
exclude it from the analysis. Let’s assume the column containing the
quality control values is named ‘QC_values’ and it contains values
between 0 and 1 with 0 meaning the poorest quality data and 1 being the
highest and that you would like to remove all data for select variables
with a quality control value below 0.5. Let’s further assume that you
would like this to apply to your incoming solar radiation variable. Then
you would declare the threshold in the **METADATA** section of your
config file:

.. code:: bash

   qc_threshold = 0.5

and in the **DATA** section of your config you will state that the
‘QC_value’ column should be applied to your incoming shortwave radiation
variable:

.. code:: bash

   shortwave_in_qc = QC_value

Now you are all set to use the functionality, note that you may apply
the same quality control value or flag column to multiple climate
variables (as shown in the first example). You may also use both numeric
qualtiy control values and character string flags for the same input
dataset although they cannot both be applied to the same variable. In
other wordsf, if you have a column of quality control numeric values it
cannot also have character strings mixed in. Another option that is used
in the example below is to declare multiple quality control flags that
should be filtered out using a comma separated list. For example in the
provided example config the flags ‘x’ and ‘b’ are used to remove select
days from incoming shorwave radiation,

.. code:: bash

   qc_flag = x, b

There is another option for specifying variables quality control
values/flags. Name the column containing the qualtiy control value/flag
in your input climate file the same as the variable it corresponds to
with the suffix "_QC". For example if your sensible heat column 
was named **sens_h** then your qualtiy control column should be named
**sens_h_QC**. If you use this option you do not need to specify the 
names in your config file. 

Example with flags and thresholds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example uses the provided time series and config files for QC flag filtering on `GitHub <https://github.com/Open-ET/flux-data-qaqc/tree/master/examples/Config_options>`__. The data is from the FLUXNET 2015 site "ARM USDA UNL OSU Woodward Switchgrass 1". This location exhibits switchgrass fields, more information on this site can be found `here <http://ameriflux.lbl.gov/sites/siteinfo/US-AR1>`__.

Because this dataset did not originally contain character Qa/Qc flags they were added for demonstration and applied to shortwave radiation. To view the list of string flags specified in the config file,

    >>> config_path = 'config_for_QC_flag_filtering.ini'
    >>> d = Data(config_path)
    >>> # view or reassign the numeric threshold specified in the config file
    >>> d.qc_threshold
        0.5

And to view the character QC flags if assigned in the config,

    >>> d.qc_flag
        ['x', 'b']

The :attr:`.Data.qc_var_pairs` attribute shows you which variables were found in your input file that have quality control values assigned, it uses the names as found in the input file,

    >>> d.qc_var_pairs
        {'LE': 'a_qc_value', 'H': 'a_qc_value', 'sw_in': 'swrad_flag'}

.. tip::
   The :attr:`.Data.qc_var_pairs` dictionary can be updated in Python to assign
   different columns of QC values to different time series variables.

Now let's apply the QC values.

Note that in this example we mixed both numeric values and threshold
with character flags, the numeric values are being applied to :math:`LE` and :math:`H`
whereas the flags (‘x’ and ‘b’) are applied to incoming shortwave
radiation.

    >>> # make copys of before and after the QC filter is applied
    >>> no_qc = d.df.input_LE.copy()
    >>> no_qc_swrad = d.df.input_sw_in.copy()
    >>> # apply QC flags/values
    >>> d.apply_qc_flags()
    >>> qc_def = d.df.input_LE.copy()
    >>> qc_flag_swrad = d.df.input_sw_in.copy()
        WARNING: renaming column Rn to input_Rn
        WARNING: renaming column G to input_G
        WARNING: renaming column LE to input_LE
        WARNING: renaming column H to input_H
        WARNING: renaming column sw_in to input_sw_in
        WARNING: renaming column sw_out to input_sw_out
        WARNING: renaming column sw_pot to input_sw_pot
        WARNING: renaming column lw_in to input_lw_in
        WARNING: renaming column lw_out to input_lw_out
        WARNING: renaming column vpd to input_vpd
        WARNING: renaming column t_avg to input_t_avg
        WARNING: renaming column ppt to input_ppt
        WARNING: renaming column ws to input_ws


.. note::
   This is a good time to point out that ``flux-data-qaqc`` may change the
   names of your input variables if they exactly match the internal names used
   by the software (see :attr:`.Data.variable_names_dict`, if this is the case
   (as is above) a warning message is printed when reading in the data
   (accessing the ``df`` or ``monthly_df`` properties of :obj:`.Data` or
   :obj:`.QaQc` for the first time) and the names will be modified with a
   prefix of "_input" as shown above.

Here is a plot showing the data before and after applying the filter.

    >>> from bokeh.plotting import ColumnDataSource, figure, show
    >>> from bokeh.models.formatters import DatetimeTickFormatter
    >>> p = figure(x_axis_label='date', y_axis_label='swrad with data removed based on QC value')
    >>> p.line(no_qc_swrad.index, no_qc_swrad, color='red', legend="no flag", line_width=2)
    >>> p.line(no_qc_swrad.index, qc_flag_swrad, color='black', legend="flag = b or x", line_width=2)
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> show(p)


.. raw:: html
    :file: _static/qc_flag1.html

And for :math:`LE`,

    >>> p = figure(x_axis_label='date', y_axis_label='LE with data removed based on QC value')
    >>> p.line(no_qc.index, no_qc, color='red', legend="no QC", line_width=2)
    >>> p.line(no_qc.index, qc_def, color='black', legend="QC=0.5", line_width=2)
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> show(p)


.. raw:: html
    :file: _static/qc_flag2.html


Alternative naming method for QC data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case the climate variables QC columns are named with the same
base name as the climate variables with the ‘\_QC’ suffix. For example
if :math:`LE` is named ‘LE_F_MDS’ in your input files header then the QC column
is named ‘LE_F_MDS_QC’. If your time series data has qualtiy control header names which follow this convention they will automatically be detected and used when you apply them using :meth:`.Data.apply_qc_flags`, i.e. the column names and variables they should be assigned to do not need to be declared in the config.ini file.

    >>> import os
    >>> config_path = os.path.join('..','Basic_usage','fluxnet_config.ini')
    >>> d = Data(config_path)
    >>> # view input files header, note the QC columns 
    >>> d.header
        Index(['TIMESTAMP', 'TA_F', 'TA_F_QC', 'SW_IN_POT', 'SW_IN_F', 'SW_IN_F_QC',
               'LW_IN_F', 'LW_IN_F_QC', 'VPD_F', 'VPD_F_QC', 'PA_F', 'PA_F_QC', 'P_F',
               'P_F_QC', 'WS_F', 'WS_F_QC', 'USTAR', 'USTAR_QC', 'NETRAD', 'NETRAD_QC',
               'PPFD_IN', 'PPFD_IN_QC', 'PPFD_OUT', 'PPFD_OUT_QC', 'SW_OUT',
               'SW_OUT_QC', 'LW_OUT', 'LW_OUT_QC', 'CO2_F_MDS', 'CO2_F_MDS_QC',
               'TS_F_MDS_1', 'TS_F_MDS_1_QC', 'SWC_F_MDS_1', 'SWC_F_MDS_1_QC',
               'G_F_MDS', 'G_F_MDS_QC', 'LE_F_MDS', 'LE_F_MDS_QC', 'LE_CORR',
               'LE_CORR_25', 'LE_CORR_75', 'LE_RANDUNC', 'H_F_MDS', 'H_F_MDS_QC',
               'H_CORR', 'H_CORR_25', 'H_CORR_75', 'H_RANDUNC', 'NEE_VUT_REF',
               'NEE_VUT_REF_QC', 'NEE_VUT_REF_RANDUNC', 'NEE_VUT_25', 'NEE_VUT_50',
               'NEE_VUT_75', 'NEE_VUT_25_QC', 'NEE_VUT_50_QC', 'NEE_VUT_75_QC',
               'RECO_NT_VUT_REF', 'RECO_NT_VUT_25', 'RECO_NT_VUT_50', 'RECO_NT_VUT_75',
               'GPP_NT_VUT_REF', 'GPP_NT_VUT_25', 'GPP_NT_VUT_50', 'GPP_NT_VUT_75',
               'RECO_DT_VUT_REF', 'RECO_DT_VUT_25', 'RECO_DT_VUT_50', 'RECO_DT_VUT_75',
               'GPP_DT_VUT_REF', 'GPP_DT_VUT_25', 'GPP_DT_VUT_50', 'GPP_DT_VUT_75',
               'RECO_SR', 'RECO_SR_N'],
              dtype='object')



Verify that the QC columns have been paired with corresponding climate variables

    >>> d.qc_var_pairs
        {'NETRAD': 'NETRAD_QC',
         'G_F_MDS': 'G_F_MDS_QC',
         'LE_F_MDS': 'LE_F_MDS_QC',
         'H_F_MDS': 'H_F_MDS_QC',
         'SW_IN_F': 'SW_IN_F_QC',
         'SW_OUT': 'SW_OUT_QC',
         'LW_IN_F': 'LW_IN_F_QC',
         'LW_OUT': 'LW_OUT_QC',
         'VPD_F': 'VPD_F_QC',
         'TA_F': 'TA_F_QC',
         'P_F': 'P_F_QC',
         'WS_F': 'WS_F_QC'}


.. note::
   FLUXNET files include their own qualtiy control flags for sensible heat and
   other variables where quality threshold columns are named the same as the
   climate variable they correspond to with the "\_QC" suffix. Therefore they
   do not need to be defined in a config file before applying them. 

For the dataset defined in the example "FLUXNET_config.ini" we did not specify a QC threshold or flag(s) in the config file, therefore we must assign it when calling the :meth:`.Data.apply_qc_flags` method (shown in :ref:`Example of threshold filtering`).

    >>> # view the QC threshold specified in the config file
    >>> print(d.qc_threshold, type(d.qc_threshold))
        None <class 'NoneType'>

Alternatively, you may assign the threshold of flag values at any time directly to a :obj:`.Data` instance:

    >>> d.qc_threshold = .75

Example of threshold filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Be sure to validate QC thresholds or flags before applying them to make sure everything seems correct. Below we see that the lowest QC values correspond with poor quality gap-fill data near the begining of the time series of sensible heat (:math:`H`). 

    >>> from bokeh.models import LinearAxis, Range1d
    >>> p = figure(x_axis_label='date', y_axis_label='sensible heat flux (w/m2)')
    >>> p.extra_y_ranges = {"sec": Range1d(start=-0.1, end=1.1)}
    >>> p.line(d.df.index, d.df['H_F_MDS'], color='red', line_width=1, legend='data')
    >>> p.add_layout(LinearAxis(y_range_name="sec", axis_label='QC value'), 'right')
    >>> p.circle(d.df.index, d.df['H_F_MDS_QC'], line_width=2, y_range_name="sec", legend='QC')
    >>> p.x_range=Range1d(d.df.index[0], d.df.index[365])
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> p.legend.location = "top_left"
    >>> show(p)
        WARNING: Insufficient data to calculate mean for multiple G measurements
        WARNING: Insufficient data to calculate mean for multiple THETA measurements


.. raw:: html
    :file: _static/qc_flag3.html

As a reminder, the routine provided for numeric or theshold filtering removes all data entries that have been assigned to a QC column and have a QC value that falls below some threshold.

    >>> # apply QC numeric threshold filters
    >>> d.apply_qc_flags(threshold=0.5)

Values with QC values < 0.5 are now removed (null) for any variable listed in :attr:`.Data.qc_var_pairs`. 

.. caution::
   The :meth:`.Data.apply_qc_flags()` method applies the filter to all
   variables in the climate file that have a QC column if columns are not
   specified in the config file.

To see all columns (variables) that may have been affected by the previous filter or to constrain them, modify the declarations in the config file or within :attr:`.Data.qc_var_pairs`, i.e.

    >>> d.qc_var_pairs
        {'NETRAD': 'NETRAD_QC',
        'G_F_MDS': 'G_F_MDS_QC',
        'LE_F_MDS': 'LE_F_MDS_QC',
        'H_F_MDS': 'H_F_MDS_QC',
        'SW_IN_F': 'SW_IN_F_QC',
        'SW_OUT': 'SW_OUT_QC',
        'LW_IN_F': 'LW_IN_F_QC',
        'LW_OUT': 'LW_OUT_QC',
        'VPD_F': 'VPD_F_QC',
        'TA_F': 'TA_F_QC',
        'P_F': 'P_F_QC',
        'WS_F': 'WS_F_QC'} 

Now let's view the same sesnible heat flux time series after applying the threshold filter, notice the strange oscillating artifact near the beginning of the time series as been removed:

    >>> p = figure(x_axis_label='date', y_axis_label='sensible heat flux (w/m2)')
    >>> p.extra_y_ranges = {"sec": Range1d(start=-0.1, end=1.1)}
    >>> p.line(d.df.index, d.df['H_F_MDS'], color='red', line_width=1, legend='data')
    >>> p.add_layout(LinearAxis(y_range_name="sec", axis_label='QC value'), 'right')
    >>> p.circle(d.df.index, d.df['H_F_MDS_QC'], line_width=2, y_range_name="sec", legend='QC')
    >>> p.x_range=Range1d(d.df.index[0], d.df.index[365])
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> p.legend.location = "top_left"
    >>> show(p)

.. raw:: html
    :file: _static/qc_flag4.html

.. seealso::
   :ref:`Step 0, manual cleaning of poor quality data` for an example that shows   how to filter poor quality data after loading data into a :obj:`.QaQc` 
   object.

--------------

Averaging data from multiple sensors
------------------------------------

Non-weighted averaging
^^^^^^^^^^^^^^^^^^^^^^

If the climate station being analyzed has multiple sensors for the same 
variable (e.g. sensible heat flux) you can easily tell ``flux-data-qaqc``
to use their non-weighted average of for ``flux-data-qaqc`` routines
including energy balance closure corrections or interactive visualizations.
To do so simply list the variable names (as found in the file header) with
a delimiter of your choice and then list the delimiter in the **METADATA**
section. Example, if you have three sensible heat variables named "h_1",
"sens_h_2", and "sensible heat, (w/m2)" then in your config file's 
**METADATA** you would write:

.. parsed-literal::

    var_name_delim = ;

and the sensible heat assignment in the **DATA** section would read:

.. parsed-literal::

    sensible_heat_flux_col = h_1;sens_h_2;sensible heat, (w/m2)

.. caution:: 
   Because there is a comma in the last variable name we cannot use a comma as
   the name delimiter. Also, if you do not state the delimiter of variable
   names in the **METADATA** section of the config file, ``flux-data-qaqc`` will
   look for the single variable name "h_1;sens_h_2;sensible heat, (w/m2)" in
   the header which will not be found.

``flux-data-qaqc`` will name the average in this case as H_mean, in general
it will add the suffix "_mean" to the internal name of the variable used 
by ``flux-data-qaqc`` which can be found in the keys of the :attr:`.Data.variable_names_dict`
dictionary.

.. hint::
   If you use any averaging option for an energy balance component, i.e.
   latent energy, sensible heat, net radiation, or soil heat flux, the average
   will also be used in energy balance closure corrections. 

Weighted averaging
^^^^^^^^^^^^^^^^^^

``flux-data-qaqc`` provides the ability to read in multiple soil heat
flux/moisture variables for a given station location, calculate their
weighted or non weighted average, and write/plot their daily and monthly
time series. *Currently weighted averaging is only provided for 
soil heat flux and soil moisture variables*, using this config option is also
the only way to automatically produce time series plots of these variables
when using :meth:`.QaQc.plot`. This may be useful for comparing/validating multiple soil
heat/moisture probes at varying locations or depths or varying
instrumentation. 

Here is what you need to do to use this functionality:

1. List the multiple soil variable names in your config file's **DATA** section
   following the convention:

-  For multiple soil heat flux variables config names should begin with
   “G\_” or “g\_” followed by an integer starting with 1,2,3,…
   i.e. g_[number]. For example:

.. code:: bash

   g_1 = name_of_my_soil_heat_flux_variable

For soil moisture variables the name of the config variable should follow
“theta_[number]” for example:

.. code:: bash

   theta_1 = name_of_my_soil_moisture_variable

2. List the units of each variable. To specify the units of your soil
   flux/moisture variables add "_units" to the config name you assigned:

.. code:: bash

   g_1_units = w/m2
   theta_1_units = cm

3. To set weights for multiple variables to compute weighted averages
   assign the "_weight" suffix to their names in the config file. For
   example, to set weights for multiple soil heat flux variables:

.. code:: bash

   g_1_weight = 0.25
   g_2_weight = 0.25
   g_3_weight = 0.5

.. hint::
   If weights are not given the arithmetic mean will be calculated. Or if the
   weights do not sum to 1, they will be automatically normalized so that they
   do.

As in the case for non-weighted averaging for any energy balance
component, if you use this option for soil heat flux (:math:`G`), the weighted 
average will also be used in energy balance closure corrections.

Weighted average example
^^^^^^^^^^^^^^^^^^^^^^^^

This example uses time series data recorded from the "ARM Southern Great Plains site- Lamont" AmeriFlux eddy covariance tower, more information on this site can be found `here <http://ameriflux.lbl.gov/sites/siteinfo/US-ARM>`__.

Here is the **DATA** section of the config file that defines the multiple :math:`G` variables in the input data file used for the example below, we put a much higher weight on the :math:`G` sensors "G_2_1_1" and "G_3_1_1",

.. code:: bash

   [DATA]
   g_1 = G_1_1_1
   g_1_units = w/m2
   g_1_weight = 1
   g_2 = G_2_1_1
   g_2_units = w/m2
   g_2_weight = 10
   g_3 = G_3_1_1
   g_3_units = w/m2
   g_3_weight = 10
   g_4 = G_4_1_1
   g_4_units = w/m2
   g_4_weight = 1
   ...

Note, the naming system of these variables (from AmeriFlux conventions) indicates that the multiple :math:`G` sensors are spaced in differing horizontal locations from one another. 

There are many soil moisture sensors at this site, because we are not using these variables within any calculations and simply want them to be loaded in and later plotted we will not assign weights to them and therefore the arithmetic mean will be calculated and added to output plots and time series files. Here is what is listed in the **DATA** section of the config file for multiple soil moisture recordings in this case:

.. code:: bash

   [DATA]
   theta_1 = SWC_1_1_1
   theta_1_units = (%): Soil water content (volumetric), range 0-100
   theta_2 = SWC_2_1_1
   theta_2_units = (%): Soil water content (volumetric), range 0-100
   theta_3 = SWC_1_2_1
   theta_3_units = (%): Soil water content (volumetric), range 0-100
   theta_4 = SWC_2_2_1
   theta_4_units = (%): Soil water content (volumetric), range 0-100
   theta_5 = SWC_3_1_1
   theta_5_units = (%): Soil water content (volumetric), range 0-100
   theta_6 = SWC_4_1_1
   theta_6_units = (%): Soil water content (volumetric), range 0-100
   theta_7 = SWC_3_2_1
   theta_7_units = (%): Soil water content (volumetric), range 0-100
   theta_8 = SWC_4_2_1
   theta_8_units = (%): Soil water content (volumetric), range 0-100
   theta_9 = SWC_1_3_1
   theta_9_units = (%): Soil water content (volumetric), range 0-100
   theta_10 = SWC_1_4_1
   theta_10_units = (%): Soil water content (volumetric), range 0-100
   theta_11 = SWC_1_5_1
   theta_11_units = (%): Soil water content (volumetric), range 0-100
   theta_12 = SWC_1_6_1
   theta_12_units = (%): Soil water content (volumetric), range 0-100
   theta_13 = SWC_2_3_1
   theta_13_units = (%): Soil water content (volumetric), range 0-100
   theta_14 = SWC_2_3_2
   theta_14_units = (%): Soil water content (volumetric), range 0-100
   theta_15 = SWC_2_2_2
   theta_15_units = (%): Soil water content (volumetric), range 0-100
   theta_16 = SWC_2_1_2
   theta_16_units = (%): Soil water content (volumetric), range 0-100
   ...

.. hint:: 
   The units for soil moisture variables will be used in the y-axis daily and
   monthly time series plots when they are created by :meth:`.QaQc.plot`.

Now that the config file has been setup, let's verify that everything was read in correctly,

    >>> # read in the data
    >>> config_path = 'config_for_multiple_soil_vars.ini'
    >>> d = Data(config_path)
    >>> # note the newly added multiple g and theta variables
    >>> d.variables
        {'date': 'TIMESTAMP_START',
        'Rn': 'NETRAD_1_1_1',
        'LE': 'LE_1_1_1',
        'H': 'H_1_1_1',
        'sw_in': 'SW_IN_1_1_1;SW_IN_1_1_2',
        'sw_out': 'SW_OUT_1_1_1',
        'lw_in': 'LW_IN_1_1_1',
        'lw_out': 'LW_OUT_1_1_1',
        'vpd': 'VPD_PI_1_1_1',
        't_avg': 'T_SONIC_1_1_1',
        'ws': 'WS_1_1_1;WS_1_2_1',
        'g_1': 'G_1_1_1',
        'g_2': 'G_2_1_1',
        'g_3': 'G_3_1_1',
        'g_4': 'G_4_1_1',
        'theta_1': 'SWC_1_1_1',
        'theta_2': 'SWC_2_1_1',
        'theta_3': 'SWC_1_2_1',
        'theta_4': 'SWC_2_2_1',
        'theta_5': 'SWC_3_1_1',
        'theta_6': 'SWC_4_1_1',
        'theta_7': 'SWC_3_2_1',
        'theta_8': 'SWC_4_2_1',
        'theta_9': 'SWC_1_3_1',
        'theta_10': 'SWC_1_4_1',
        'theta_11': 'SWC_1_5_1',
        'theta_12': 'SWC_1_6_1',
        'theta_13': 'SWC_2_3_1',
        'theta_14': 'SWC_2_3_2',
        'theta_15': 'SWC_2_2_2',
        'theta_16': 'SWC_2_1_2'}

Note, the windspeed and shortwave incoming radtiation columns were assigned multiple variables as well, these will be used to calculate the non-weighted mean as described in :ref:`Non-weighted averaging`.

Check the units assignment:

    >>> d.units
        {'Rn': 'w/m2',
         'LE': 'w/m2',
         'H': 'w/m2',
         'sw_in': 'w/m2',
         'sw_out': 'w/m2',
         'lw_in': 'w/m2',
         'lw_out': 'w/m2',
         'vpd': 'hPa',
         't_avg': 'C',
         'ws': 'm/s',
         'g_1': 'w/m2',
         'g_2': 'w/m2',
         'g_3': 'w/m2',
         'g_4': 'w/m2',
         'theta_1': '(%): Soil water content (volumetric), range 0-100',
         'theta_2': '(%): Soil water content (volumetric), range 0-100',
         'theta_3': '(%): Soil water content (volumetric), range 0-100',
         'theta_4': '(%): Soil water content (volumetric), range 0-100',
         'theta_5': '(%): Soil water content (volumetric), range 0-100',
         'theta_6': '(%): Soil water content (volumetric), range 0-100',
         'theta_7': '(%): Soil water content (volumetric), range 0-100',
         'theta_8': '(%): Soil water content (volumetric), range 0-100',
         'theta_9': '(%): Soil water content (volumetric), range 0-100',
         'theta_10': '(%): Soil water content (volumetric), range 0-100',
         'theta_11': '(%): Soil water content (volumetric), range 0-100',
         'theta_12': '(%): Soil water content (volumetric), range 0-100',
         'theta_13': '(%): Soil water content (volumetric), range 0-100',
         'theta_14': '(%): Soil water content (volumetric), range 0-100',
         'theta_15': '(%): Soil water content (volumetric), range 0-100',
         'theta_16': '(%): Soil water content (volumetric), range 0-100'}


View these variables and their weights as written in the config file:

    >>> d.soil_var_weight_pairs
        {'g_1': {'name': 'added_G_col', 'weight': '6'},
         'g_2': {'name': 'another_G_var', 'weight': '2'},
         'g_3': {'name': 'G', 'weight': '0.5'},
         'g_4': {'name': 'final_G_var', 'weight': '0.25'},
         'g_5': {'name': 'yet_another_G', 'weight': '0.25'},
         'theta_1': {'name': 'soil_moisture_z1', 'weight': '0.25'},
         'theta_2': {'name': 'soil_moisture_z10', 'weight': '0.75'}}

When the data is first loaded into memory the weighted (and non-weighted) averages are calculated. At this stage weights will be automatically normalized so that they sum to one and the new weights will be printed if this occurs.

    >>> # load daily or monthly dataframe to calculate the weighted averages if they exist
    >>> d.df.head();
        g weights not given or don't sum to one, normalizing
        Here are the new weights:
         G_1_1_1:0.05, G_2_1_1:0.45, G_3_1_1:0.45, G_4_1_1:0.05
        Calculating mean for var: THETA from columns: ['SWC_1_1_1', 'SWC_2_1_1', 'SWC_1_2_1', 'SWC_2_2_1', 'SWC_3_1_1', 'SWC_4_1_1', 'SWC_3_2_1', 'SWC_4_2_1', 'SWC_1_3_1', 'SWC_1_4_1', 'SWC_1_5_1', 'SWC_1_6_1', 'SWC_2_3_1', 'SWC_2_3_2', 'SWC_2_2_2', 'SWC_2_1_2']
        Calculating mean for var: sw_in
         from columns: ['SW_IN_1_1_1', 'SW_IN_1_1_2']
        Calculating mean for var: ws
         from columns: ['WS_1_1_1', 'WS_1_2_1']

In this example, shortwave incoming radiation and windspeed were also averaged (non-weighted) from multiple recordings as described in :ref:`Non-weighted averaging`.

The weights have been changed and updated as we would expect for :math:`G`, you may ignore the weights for soil moisture in this case- because they were not assigned the arithmetic mean is calculated and the weights are not used.

    >>> d.soil_var_weight_pairs
        {'g_1': {'name': 'G_1_1_1', 'weight': 0.045454545454545456},
         'g_2': {'name': 'G_2_1_1', 'weight': 0.45454545454545453},
         'g_3': {'name': 'G_3_1_1', 'weight': 0.45454545454545453},
         'g_4': {'name': 'G_4_1_1', 'weight': 0.045454545454545456},
         'theta_1': {'name': 'SWC_1_1_1', 'weight': 1},
         'theta_2': {'name': 'SWC_2_1_1', 'weight': 1},
         'theta_3': {'name': 'SWC_1_2_1', 'weight': 1},
         'theta_4': {'name': 'SWC_2_2_1', 'weight': 1},
         'theta_5': {'name': 'SWC_3_1_1', 'weight': 1},
         'theta_6': {'name': 'SWC_4_1_1', 'weight': 1},
         'theta_7': {'name': 'SWC_3_2_1', 'weight': 1},
         'theta_8': {'name': 'SWC_4_2_1', 'weight': 1},
         'theta_9': {'name': 'SWC_1_3_1', 'weight': 1},
         'theta_10': {'name': 'SWC_1_4_1', 'weight': 1},
         'theta_11': {'name': 'SWC_1_5_1', 'weight': 1},
         'theta_12': {'name': 'SWC_1_6_1', 'weight': 1},
         'theta_13': {'name': 'SWC_2_3_1', 'weight': 1},
         'theta_14': {'name': 'SWC_2_3_2', 'weight': 1},
         'theta_15': {'name': 'SWC_2_2_2', 'weight': 1},
         'theta_16': {'name': 'SWC_2_1_2', 'weight': 1}}

Now the dataframe also has the weighted means that will be named g_mean and theta_mean,

    >>> d.df.columns
        Index(['input_t_avg', 'input_sw_pot', 'input_sw_in', 'input_lw_in',
               'input_vpd', 'input_ppt', 'input_ws', 'input_Rn', 'input_sw_out',
               'input_lw_out', 'input_G', 'input_LE', 'LE_corrected', 'input_H',
               'H_corrected', 'added_G_col', 'another_G_var', 'final_G_var',
               'yet_another_G', 'soil_moisture_z1', 'soil_moisture_z10', 'a_qc_value',
               'swrad_flag', 'g_mean', 'theta_mean'],
              dtype='object')

.. note:: 
   Even though we did not specify "ground_flux_col" in the config file, the
   weighted average value has now been used to update this variable. Therefore
   the weighted mean will be used in energy balance closure correction routines
   if they are subsequently run.

Check which variable will be used as :math:`G` later if closure corrections are used:

    >>> d.variables.get('G')
        'g_mean'

Now, let's visualize the resulting weighted average of multiple :math:`G` measurements and their individual daily time series,

    >>> # get just G columns for plot arguments
    >>> G_cols = [c for c in d.df.columns if c.startswith(('g_','G_'))]
    >>> G_cols
        ['G_1_1_1', 'G_2_1_1', 'G_3_1_1', 'G_4_1_1', 'g_mean']

The example below creates the time series plot with a short span of data for easier visibility of weighted mean, it also used the plot routines provided by :obj:`.Data` and :obj:`.QaQc` which are inhereted from the :obj:`.Plot` class within ``flux-data-qaqc``. Specifically this example utilizes :meth:`.Plot.add_lines` which makes the time series plotting of multiple variables more efficient and automatically handles the hover tooltips.

    >>> from fluxdataqaqc.plot import ColumnDataSource # for hover tooltips
    >>> # shorter period for visualization
    >>> df = d.df.loc['01/01/2008':'05/01/2008', G_cols]
    >>> plt_vars = G_cols
    >>> colors = ['blue', 'red', 'orange', 'green', 'black']
    >>> x_name = 'date'
    >>> source = ColumnDataSource(df)
    >>> fig = figure(x_axis_label='date', y_axis_label='Soil heat flux (w/m2)')
    >>> Data.add_lines(fig, df, plt_vars, colors, x_name, source, labels=G_cols)
    >>> show(fig)

.. raw:: html
    :file: _static/weighted_g.html

Note, the weighted mean is closer to 'G_2_1_1' and 'G_3_1_1' as we gave them weights of 10 versus 1 to 'G_1_1_1' and 'G_4_1_1'.

Lastly, the code snippets below run the Energy Balance Ratio closure correction and creating the default plots in order to view the daily and monthly time series of multiple soil moisture variables. It also shows how to upload the output plot file into a Jupyter Notebook for viewing.

    >>> # in order to correctly view the output in a Jupyter notebook
    >>> from bokeh.io import output_notebook
    >>> output_notebook()

Within the set of default plots created by the :meth:`.QaQc.plot` method will include interactive daily and monthly time series of multiple :math:`G` and soil moisture variables if they were assigned in the input config file (as in this example), scroll down to view them. 

    >>> from fluxdataqaqc import QaQc
    >>> q = QaQc(d)
    >>> q.correct_data()
    >>> # this will NOT save the plot file, use output_type='save'
    >>> q.plot(output_type='show')


.. raw:: html
    :file: _static/US-ARM_multipe_soilvars_plots.html
 

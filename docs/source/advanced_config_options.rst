.. There are two major functionalities in
   ``flux-data-qaqc``, first, correcting surface energy balance by
   adjusting latent energy and sensible heat fluxes and calculate other
   climatic variables. Second, it serves as a robust way to read in
   different time series data and produce visualizations, e.g. their daily
   and monthly time series.

Configuration Options and Caveats
=================================

This tutorial shows how to use ``flux-data-qaqc`` with climate data
of various formats and generally covers formatting rules of input data and extra
options that can be set in a config file. The major differences when using ``flux-data-qaqc`` for
different input data lie in the config file declarations therefore the
entire workflow from the `FLUXNET 2015 example
notebook <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Basic_usage/FLUXNET_2015_example.ipynb>`__
will work just the same once your config file is set up correctly. 

The data used herein is provided with ``flux-data-qaqc`` and can be
downloaded
`here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Config_options>`__.
The USGS data is from an eddy covariance flux tower for Dixie Valey
Dense Vegetation. Details on the data can be found in this
`report <https://pubs.usgs.gov/pp/1805/pdf/pp1805.pdf>`__. The other data used 
is from the FLUXNET 2015 dataset and can be downloaded 
`here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Basic_usage>`__.

A reproducible Jupyter Notebook with minor differences of this tutorial can be found 
`here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Config_options/advanced_config_options.ipynb>`__.


Seting up a config file
-----------------------

``flux-data-qaqc`` start with the creation of a configuration file, a
text file with extension ".ini" or ".INI" that follows the rules set 
`here <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`__.
A config file for using ``flux-data-qaqc`` requires the sections: 1. **METADATA** 2. **DATA**
although you may provide additional for custom uses. 

In **METADATA** you may enter any metadata that you wish so long
as the entrie's *key* is followed by an equal sign and the assigned 
*value*, i.e. 

.. parsed-literal::

    key = value

However there are some mandatory metadata entries unique to ``flux-data-qaqc``.
The “climate_file_path” is the full or relative file
path of the input climate file (excel or CSV, more on formatting this
below) containing the climatic data to be analyzed. The
“station_elevation” (meters) and latitude (decimal degrees)
fields are used to calculate clear sky potential solar radiation using
an ASCE formulation and downloading reference ET. “site_id” is used for 
saving output files. Other metadata entries that are optional but used by 
``flux-data-qaqc`` include “missing_data_value”, “qc_threshold”, “qc_flag”,
"var_name_delim", "skiprows", "date_parser", and "gridmet_file_path".
“missing_data_value” is used to correctly parse missing values in the 
input climate time series. All other optional metadata that can be used by
``flux-data-qaqc`` except "gridmet_file_path" (which is simply the path 
to a file that is downloaded by :meth:`.QaQc.download_gridMET`) are explained 
within this page.

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

You may view these climate entry keys (as found in the config file) from within
Python using the :attr:`.Data.config` property which contains all information
listed in the config file as a :obj:`configparser.ConfigParser` instance.

    >>> config_path = 'USGS_config.ini'
    >>> d = Data(config_path)
    >>> # loop through a list of tuples with keys and values from DATA section
    >>> for each in d.config.items('DATA'):
    >>>     print(each[0]) 
        datestring_col
        net_radiation_col
        net_radiation_units
        ground_flux_col
        ground_flux_units
        latent_heat_flux_col
        latent_heat_flux_units
        latent_heat_flux_corrected_col
        latent_heat_flux_corrected_units
        sensible_heat_flux_col
        sensible_heat_flux_units
        sensible_heat_flux_corrected_col
        sensible_heat_flux_corrected_units
        shortwave_in_col
        shortwave_in_units
        shortwave_out_col
        shortwave_out_units
        shortwave_pot_col
        shortwave_pot_units
        longwave_in_col
        longwave_in_units
        longwave_out_col
        longwave_out_units
        vap_press_col
        vap_press_units
        vap_press_def_col
        vap_press_def_units
        avg_temp_col
        avg_temp_units
        precip_col
        precip_units
        wind_spd_col
        wind_spd_units

You can also access the data from the :attr:`.Data.config` as a dictionary,
for example if your **METADATA** section has an entry for "land_cover", e.g.

.. parsed-literal::
    
    [METADATA]
    land_cover = CROP
    ...

You can access this value specifically knowing the config section and key name:

   >>> d.config.get('METADATA', 'land_cover')
       CROP

A useful tip, if you are unsure if your config file's metadata contains
a specific entry you can pass the ``fallback`` keyword-only argument to the
:meth:`configparser.ConfigParser.get` method similar to a Python dictionary.

   >>> d.config.get('METADATA', 'land_cov', fallback='not given')
       "not given"


Input formatting and caveats
----------------------------

Missing data
^^^^^^^^^^^^

For reading certain values as null or missing data points assign the
“missing_data_value” to the **METADATA** section of the config file. 
The value should be numeric, e.g.  

.. parsed-literal::

    missing_data_value = -999

If the input dataset does not contain all of expected climate variables 
as found in in your data, if
this is the case you may specify them as missing (‘na’) in your
config file or simply do not list them. Missing variables will be ignored 
for the most part and will not be present in output files/plots, however 
if key variables for the energy balance are not present (LE, H, G, and Rn) 
then you will not be able to run energy balance closure correction routines.

Data file format
^^^^^^^^^^^^^^^^

``flux-data-qaqc`` accepts Excel files (.xlx and .xlsx) and
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

**Note:** if the the input datas temporal frequency is not recognized
``flux-data-qaqc`` will attempt to resample it to daily frequency. Also,
if a value is not recognized a numeric in any data column it will be
forced to a null value.

Data header formatting
^^^^^^^^^^^^^^^^^^^^^^

A common format of some time series data is that the header row may
not start on the first line of the file. If this is the case you must add
an entry to the **METADATA** section of the config file "skiprows" which
stats the number of rows to skip before finding the header row. A 
caveat is that if using excel files you must also ensure that the lines
before the header row begin with a hashtag symbol "#". 

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
‘QC_flag’ column should be applied to your LE and H variables:

.. code:: bash

   latent_heat_flux_qc = QC_flag
   sensible_heat_flux_qc = QC_flag

Now, when the :meth:`.Data.apply_qc_flags` method is used the all date
entries of LE and H that have a "QC_flag" value of 'b' will be forced 
to null in the daily (:attr:`.Data.df`) and monthly (:attr:`.Data.monthly_df`)
datasets of a :obj:`.Data` instance. 

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

Below is an example using the provided FLUXNET 2015 file which includes its 
own qualtiy control flags for sensible heat and others. Note that if your datas
qualtiy control header names follow this convention they will 
automatically be detected and used when you apply them using 
:meth:`.Data.apply_qc_flags`.

    >>> config_path = 'multiple_soilflux_config.ini'
    >>> d = Data(config_path)
    >>> # view or reassign the numeric threshold specified in the config file
    >>> d.qc_threshold
        0.5

View the list of string flags specified in the config file,

    >>> d.qc_flag
        ['x', 'b']



The :attr:`.Data.qc_var_pairs` attribute shows you which variables were found in your input file that have quality control values assigned, it uses the names as found in the input file,

    >>> d.qc_var_pairs
        {'LE': 'a_qc_value', 'H': 'a_qc_value', 'sw_in': 'swrad_flag'}


Now let's apply the QC values.

Note that in this example we mixed both numeric values and threshold
with character flags, the numeric values are being applied to LE and H
whereas the flags (‘x’ and ‘b’) are applied to incoming shortwave
radiation.

    >>> # make copys of before and after the QC filter is applied
    >>> no_qc = d.df.input_LE.copy()
    >>> no_qc_swrad = d.df.input_sw_in.copy()
    >>> # apply QC flags/values
    >>> d.apply_qc_flags()
    >>> qc_def = d.df.input_LE.copy()
    >>> qc_flag_swrad = d.df.input_sw_in.copy()
        g weights not given or don't sum to one, normalizing
        Here are the new weights:
         added_G_col:0.67, another_G_var:0.22, G:0.06, final_G_var:0.03, yet_another_G:0.03
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


This is a good time to point out that ``flux-data-qaqc`` may change the
names of your input variables if they exactly match the internal names
used by the software (see :attr:`.Data.variable_names_dict`, if this is 
the case (as is above) a warning message is printed when reading in 
the data (accessing the ``df`` or ``monthly_df`` properties of :obj:`.Data`
or :obj:`.QaQc` for the first time) and the names will be modified with a
prefix of "_input" as shown above.

Here is a plot showing the data before and after applying the filter.


    >>> p = figure(x_axis_label='date', y_axis_label='swrad with data removed based on QC value')
    >>> p.line(no_qc_swrad.index, no_qc_swrad, color='red', legend="no flag", line_width=2)
    >>> p.line(no_qc_swrad.index, qc_flag_swrad, color='black', legend="flag = b or x", line_width=2)
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> show(p)


.. raw:: html
    :file: _static/qc_flag1.html

And for LE,

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
if LE is named ‘LE_F_MDS’ in your input files header then the QC column
is named ‘LE_F_MDS_QC’.

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



**Note,** for this dataset we did not specify a QC threshold or flag(s) in the config.
We can assign it when calling the :meth:`.Data.apply_qc_flags` method.

    >>> # view the QC threshold specified in the config file
    >>> print(d.qc_threshold, type(d.qc_threshold))
        None <class 'NoneType'>


Example of threshold filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you create your own QC values be sure to validate them to make sure
everything seems correct. Below we see that the lowest QC values
correspond with poor quality gap-fill data near the begining of the
dataset.

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


The routine provided removes all data that falls below a QC value of
0.5, although this can be modified. Also see the `provided FLUXNET Jupyter
notebook <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/FLUXNET_2015_example.ipynb>`__
for more examples.

Now let's apply our threshold filter of sensible heat with QC 
values < 0.5 are now removed (null).

Note, the ``Data.apply_qc_flags()`` method applies the filter to all
variables in the climate file that have a QC column if columns are not
specified in the config file.

    >>> # apply QC filters
    >>> d.apply_qc_flags(threshold=0.5)
    >>> # same figure
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

``flux-data-qaqc`` will name the average in this case as H_mean, in general
it will add the suffix "_mean" to the internal name of the variable used 
by ``flux-data-qaqc`` which can be found in the keys of the :attr:`.Data.variable_names_dict`
dictionary.

**Note,** that because there is a comma in the last variable name we cannot
use a comma as the name delimiter. Also, if you do not state the delimiter
of variable names in the **METADATA** section of the config file ``flux-data-qaqc``
will look for the single variable name "h_1;sens_h_2;sensible heat, (w/m2)"
in the header which will not be found.

**Note,** if you use this option for any energy balance component, i.e.
latent energy, sensible heat, net radiation, or soil heat flux, the 
average will also be used in energy balance closure corrections. 

Weighted averaging
^^^^^^^^^^^^^^^^^^

``flux-data-qaqc`` provides the ability to read in multiple soil heat
flux/moisture variables for a given station location, calculate their
weighted or non weighted average, and write/plot their daily and monthly
time series. Currently the weighted averaging is only provided for 
soil heat flux and soil moisture variables, using this config option is also
the only way to automatically produce time series plots of these variables
when using :meth:`.QaQc.plot`. This may be useful for comparing/validating multiple soil
heat/moisture probes at varying locations or depths or varying
instrumentation. 

Here is what you need to do to use this functionality:

1. List the multiple soil variable names in your config file following
   the convention:

-  For multiple soil heat flux variables config names should begin with
   “G\_” or “g\_” followed by an integer starting with 1,2,3,…
   i.e. g_[number]. For example:

.. code:: bash

   g_1 = name_of_my_soil_heat_flux_variable

-  For soil moisture variables the name of the config variable should
   follow “theta_[number]” for example:

.. code:: bash

   theta_1 = name_of_my_soil_moisture_variable

2. List the units and (optionally) weights of the multiple variables

-  To specify the units of your soil flux/moisture variables add
   "_units" to the config name you assigned:

.. code:: bash

   g_1_units = w/m2
   theta_1_units = cm

-  To set weights for multiple variables to compute weighted averages
   assign the "_weight" suffix to their names in the config file. For
   example, to set weights for multiple soil heat flux variables:

.. code:: bash

   g_1_weight = 0.25
   g_2_weight = 0.25
   g_3_weight = 0.5

Note, if weights are not given the arithmetic mean will be
calculated, if the weights do not sum to 1, they will be
automatically normalized so that they do.

As in the case for non-weighted averaging for any energy balance
component, if you use this option for soil heat flux (G), the weighted 
average will also be used in energy balance closure corrections.

Weighted average example
^^^^^^^^^^^^^^^^^^^^^^^^

The provided multiple soil variable config and input data are used for
these examples.

Here is the section of the config file that defines the multiple soil
variables in the input climate file used for the example below:

.. code:: bash

   g_1 = added_G_col
   g_1_weight = 6
   g_1_units = w/m2
   g_2 = another_G_var
   g_2_weight = 2
   g_2_units = w/m2
   # note the next variable is the same that was assigned as the main soil heat flux variable
   # i.e. ground_flux_col = G
   g_3 = G
   g_3_weight = 0.5
   g_3_units = w/m2

   theta_1 = soil_moisture_z1
   theta_1_weight = 0.25
   theta_1_units = cm
   theta_2 = soil_moisture_z10
   theta_2_weight = 0.75
   theta_2_units = cm

Verify that everything was read in correctly,

    >>> # read in the data
    >>> config_path = 'multiple_soilflux_config.ini'
    >>> d = Data(config_path)
    >>> # note the newly added multiple g and theata variables
    >>> d.variables
        {'date': 'date',
         'year': 'na',
         'month': 'na',
         'day': 'na',
         'Rn': 'Rn',
         'G': 'G',
         'LE': 'LE',
         'LE_user_corr': 'LE_corrected',
         'H': 'H',
         'H_user_corr': 'H_corrected',
         'sw_in': 'sw_in',
         'sw_out': 'sw_out',
         'sw_pot': 'sw_pot',
         'lw_in': 'lw_in',
         'lw_out': 'lw_out',
         'vp': 'na',
         'vpd': 'vpd',
         't_avg': 't_avg',
         'ppt': 'ppt',
         'ws': 'ws',
         'g_1': 'added_G_col',
         'g_2': 'another_G_var',
         'g_3': 'G',
         'g_4': 'final_G_var',
         'g_5': 'yet_another_G',
         'theta_1': 'soil_moisture_z1',
         'theta_2': 'soil_moisture_z10',
         'LE_qc_flag': 'a_qc_value',
         'H_qc_flag': 'a_qc_value',
         'sw_in_qc_flag': 'swrad_flag'}


Also check and the units assignment:

    >>> d.units
        {'Rn': 'w/m2',
         'G': 'w/m2',
         'LE': 'w/m2',
         'LE_user_corr': 'w/m2',
         'H': 'w/m2',
         'H_user_corr': 'w/m2',
         'sw_in': 'w/m2',
         'sw_out': 'w/m2',
         'sw_pot': 'w/m2',
         'lw_in': 'w/m2',
         'lw_out': 'w/m2',
         'vp': 'na',
         'vpd': 'hPa',
         't_avg': 'C',
         'ppt': 'mm',
         'ws': 'm/s',
         'g_1': 'w/m2',
         'g_2': 'w/m2',
         'g_4': 'w/m2',
         'g_5': 'w/m2',
         'theta_1': 'cm',
         'theta_2': 'cm'}


View these variables and their weights only,

    >>> d.soil_var_weight_pairs
        {'g_1': {'name': 'added_G_col', 'weight': '6'},
         'g_2': {'name': 'another_G_var', 'weight': '2'},
         'g_3': {'name': 'G', 'weight': '0.5'},
         'g_4': {'name': 'final_G_var', 'weight': '0.25'},
         'g_5': {'name': 'yet_another_G', 'weight': '0.25'},
         'theta_1': {'name': 'soil_moisture_z1', 'weight': '0.25'},
         'theta_2': {'name': 'soil_moisture_z10', 'weight': '0.75'}}


When the data is first loaded into memory the weighted averages are calculated.
At this stage weights will be automatically normalized so that they sum
to one and the new weights will be printed if this occurs.

    >>> # call daily or monthly dataframe to calculate the weighted averages if they exist
    >>> d.df.head()
        g weights not given or don't sum to one, normalizing
        Here are the new weights:
         added_G_col:0.67, another_G_var:0.22, G:0.06, final_G_var:0.03, yet_another_G:0.03
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


.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>input_t_avg</th>
          <th>input_sw_pot</th>
          <th>input_sw_in</th>
          <th>input_lw_in</th>
          <th>input_vpd</th>
          <th>input_ppt</th>
          <th>input_ws</th>
          <th>input_Rn</th>
          <th>input_sw_out</th>
          <th>input_lw_out</th>
          <th>...</th>
          <th>added_G_col</th>
          <th>another_G_var</th>
          <th>final_G_var</th>
          <th>yet_another_G</th>
          <th>soil_moisture_z1</th>
          <th>soil_moisture_z10</th>
          <th>a_qc_value</th>
          <th>swrad_flag</th>
          <th>g_mean</th>
          <th>theta_mean</th>
        </tr>
        <tr>
          <th>date</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>2009-01-01</th>
          <td>2.803</td>
          <td>186.710</td>
          <td>123.108</td>
          <td>261.302</td>
          <td>1.919</td>
          <td>0.0</td>
          <td>3.143</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>20.573270</td>
          <td>26.942860</td>
          <td>0</td>
          <td>x</td>
          <td>NaN</td>
          <td>25.350463</td>
        </tr>
        <tr>
          <th>2009-01-02</th>
          <td>2.518</td>
          <td>187.329</td>
          <td>121.842</td>
          <td>268.946</td>
          <td>0.992</td>
          <td>0.0</td>
          <td>2.093</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>20.250870</td>
          <td>26.601709</td>
          <td>0</td>
          <td>x</td>
          <td>NaN</td>
          <td>25.013999</td>
        </tr>
        <tr>
          <th>2009-01-03</th>
          <td>5.518</td>
          <td>188.008</td>
          <td>124.241</td>
          <td>268.004</td>
          <td>2.795</td>
          <td>0.0</td>
          <td>4.403</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>20.827236</td>
          <td>26.644598</td>
          <td>0</td>
          <td>x</td>
          <td>NaN</td>
          <td>25.190258</td>
        </tr>
        <tr>
          <th>2009-01-04</th>
          <td>-3.753</td>
          <td>188.742</td>
          <td>113.793</td>
          <td>246.675</td>
          <td>0.892</td>
          <td>0.0</td>
          <td>4.336</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>20.988757</td>
          <td>26.843588</td>
          <td>0</td>
          <td>x</td>
          <td>NaN</td>
          <td>25.379880</td>
        </tr>
        <tr>
          <th>2009-01-05</th>
          <td>-2.214</td>
          <td>189.534</td>
          <td>124.332</td>
          <td>244.478</td>
          <td>1.304</td>
          <td>0.0</td>
          <td>2.417</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>20.756527</td>
          <td>26.262146</td>
          <td>0</td>
          <td>x</td>
          <td>NaN</td>
          <td>24.885741</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 25 columns</p>
    </div>


Note the weights have been changed and updated 

    >>> d.soil_var_weight_pairs
        {'g_1': {'name': 'added_G_col', 'weight': 0.6666666666666666},
         'g_2': {'name': 'another_G_var', 'weight': 0.2222222222222222},
         'g_3': {'name': 'G', 'weight': 0.05555555555555555},
         'g_4': {'name': 'final_G_var', 'weight': 0.027777777777777776},
         'g_5': {'name': 'yet_another_G', 'weight': 0.027777777777777776},
         'theta_1': {'name': 'soil_moisture_z1', 'weight': '0.25'},
         'theta_2': {'name': 'soil_moisture_z10', 'weight': '0.75'}}


Now the dataframe also has the weighted means that will be named g_mean and theta_mean,

    >>> d.df.columns
        Index(['input_t_avg', 'input_sw_pot', 'input_sw_in', 'input_lw_in',
               'input_vpd', 'input_ppt', 'input_ws', 'input_Rn', 'input_sw_out',
               'input_lw_out', 'input_G', 'input_LE', 'LE_corrected', 'input_H',
               'H_corrected', 'added_G_col', 'another_G_var', 'final_G_var',
               'yet_another_G', 'soil_moisture_z1', 'soil_moisture_z10', 'a_qc_value',
               'swrad_flag', 'g_mean', 'theta_mean'],
              dtype='object')


The weighted mean is closest to the variable assigned to “g_1” which had the highest weight.

    >>> p = figure(x_axis_label='date', y_axis_label='Soil heat flux')
    >>> p.line(d.df.index, d.df['g_mean'], color='black', legend="weighted mean", line_width=2)
    >>> p.line(d.df.index, d.df['added_G_col'], color='orange', legend="g_1: 0.71", line_width=1)
    >>> p.line(d.df.index, d.df['another_G_var'], color='green', legend="g_2: 0.24", line_width=1)
    >>> p.line(d.df.index, d.df['input_G'], color='red', legend="g_3: 0.60", line_width=1)
    >>> 
    >>> p.x_range=Range1d(d.df.index[150], d.df.index[160])
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> show(p)


.. raw:: html
    :file: _static/weighted_g.html



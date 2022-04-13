
Tutorial
========

This tutorial demonstrates the most important features of the
``flux-data-qaqc`` Python package for management, analysis, and visualization
of eddy covariance time series data. It is recommended to read the
:ref:`Installation` and :ref:`Configuration Options and Caveats` tutorials before this one. 

A Jupyter Notebook of this tutorial is available `here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Basic_usage/Tutorial.ipynb>`__.

.. Tip:: 
   Currently, the software does not include a command line interface therefore
   to use the software you must use Python, e.g. make your own scripts or use
   an interactive shell. However, you will see that common workflows can be
   accomplished with a few (5-10) lines of code and you can simply follow the
   templates given here to make custom scripts.

Description of example datasets
-------------------------------

The data for this example comes from the “Twitchell Alfalfa” AmeriFlux
eddy covariance flux tower site in California. The site is located in
alfalfa fields and exhibits a mild Mediterranean climate with dry and
hot summers, for more information on this site or to download data click
`here <https://ameriflux.lbl.gov/sites/siteinfo/US-Tw3>`__.

Loading input
-------------

The loading and management of input climatic data and metadata from a
config.ini file is done using the :obj:`fluxdataqaqc.Data` object. In a
nutshell, a :obj:`.Data` object is created from a properly formatted config file
(see :ref:`Setting up a config file`) and has tools for parsing input climate
data, averaging input climate time series, accessing/managing metadata,
flag-based data filtering, and creating interactive visualizations of input
data.

There is only one argument to create a Data object, the path to the
config.ini file:

    >>> # imports for code snippets within tutorial
    >>> import pandas as pd
    >>> from fluxdataqaqc import Data, QaQc, Plot
    >>> from bokeh.plotting import figure, show, ColumnDataSource
    >>> from bokeh.models.formatters import DatetimeTickFormatter
    >>> from bokeh.models import LinearAxis, Range1d
    >>> # create a Data object from the config.ini file    
    >>> config_path = 'US-Tw3_config.ini'
    >>> d = Data(config_path)

Attributes of a Data object
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below are some of the useful attributes of the :obj:`.Data` object and how
they may be used.

The full path to the config.ini file that was used to create the
:obj:`.Data` instance can be accessed, note that it will return a
system-depenedent :obj:`pathlib.Path` object. E.g. on my Linux machine the
path is:

    >>> d.config_file
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/US-Tw3_config.ini')



On a Windows machine the path will have the appropriate backslashes.

Similarly to access the climate time series file:

    >>> d.climate_file
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/AMF_US-Tw3_BASE_HH_5-5.csv')



The :attr:`.Data.config` attribute is a :obj:`configparser.ConfigParser` object,
it allows you to access metadata and data in the config file in multiple
ways and to modify them. In ``flux-data-qaqc`` it is mainly used for
accessing information about the input data.

    >>> # get a list of all entries in the METADATA section of the config.ini
    >>> d.config.items('METADATA') # access the DATA section the same way
        [('climate_file_path', 'AMF_US-Tw3_BASE_HH_5-5.csv'),
         ('station_latitude', '38.1159'),
         ('station_longitude', '-121.6467'),
         ('station_elevation', '-9.0'),
         ('missing_data_value', '-9999'),
         ('skiprows', '2'),
         ('date_parser', '%Y%m%d%H%M'),
         ('site_id', 'US-Tw3'),
         ('country', 'USA'),
         ('doi_contributor_name', 'Dennis Baldocchi'),
         ('doi_contributor_role', 'Author'),
         ('doi_contributor_email', 'baldocchi@berkeley.edu'),
         ('doi_contributor_institution', 'University of California, Berkeley'),
         ('doi_organization', 'California Department of Water Resources'),
         ('doi_organization_role', 'Sponsor'),
         ('flux_measurements_method', 'Eddy Covariance'),
         ('flux_measurements_variable', 'CO2'),
         ('flux_measurements_operations', 'Continuous operation'),
         ('site_name', 'Twitchell Alfalfa'),
         ('igbp', 'CRO'),
         ('igbp_comment',
          'alfalfa is a fast growing leguminous crop raised for animal feed of low stature.  It is planted in rows and typically reaches 60-70 cm in height prior to harvest.'),
         ('land_ownership', 'public'),
         ('network', 'AmeriFlux'),
         ('reference_paper',
          'Baldocchi, D., Penuelas, J. (2018) The Physics And Ecology Of Mining Carbon Dioxide From The Atmosphere By Ecosystems, Global Change Biology, 45(), 9275–9287'),
         ('reference_doi', '10.1111/gcb.14559'),
         ('reference_usage', 'Reference'),
         ('research_topic',
          'The research approach of the University of California, Berkeley Biometeorology Laboratory involves the coordinated use of experimental measurements and theoretical models to understand the physical, biological, and chemical processes that control trace gas fluxes between the biosphere and atmosphere and to quantify their temporal and spatial variations. The research objectives of the Mayberry Wetland, Twitchell Wetland, Sherman Island, Twitchell Island, Twitchell Alfalfa,  and Twitchell Corn sites are as follows: 1) Describe differences in the fluxes of CO2, CH4, H2O, and energy between different land uses, 2) Understand the mechanisms controlling these fluxes, 3) Use ecosystem modeling to understand controls on these mechanisms under different environmental scenarios. These six sites were selected to capture a wide range of inundated conditions within the Sacramento-San Joaquin River Delta. The research focuses on the eddy covariance technique to measure CH4, CO2, H2O, and energy fluxes and works to combine measurements of both net fluxes and partitioned fluxes in order to achieve a mechanistic understanding of the ecological controls on current and future carbon flux in the Delta.'),
         ('terrain', 'Flat'),
         ('aspect', 'FLAT'),
         ('wind_direction', 'W'),
         ('surface_homogeneity', '370.0'),
         ('site_desc',
          "The Twitchell Alfalfa site is an alfalfa field owned by the state of California and leased to third parties for farming. The tower was installed on May 24, 2013. This site and the surrounding region are part of the San Joaquin - Sacramento River Delta drained beginning in the 1850's and subsequently used for agriculture. The field has been alfalfa for X years…., Crop rotation occurs every 5-6 years.  The site is harvested by mowing and bailing several times per year.  The field is fallow typically between November and February. The site is irrigated by periodically-flooded ditches surrounding the field. The site is irrigated by raising, and subsequently lowering the water table??"),
         ('site_funding', 'California Department of Water Resources'),
         ('team_member_name', 'Joe Verfaillie'),
         ('team_member_role', 'Technician'),
         ('team_member_email', 'jverfail@berkeley.edu'),
         ('team_member_institution', 'University of California, Berkeley'),
         ('url_ameriflux', 'http://ameriflux.lbl.gov/sites/siteinfo/US-Tw3'),
         ('utc_offset', '-8'),
         ('mat', '15.6'),
         ('map', '421.0'),
         ('land_owner', 'California Department of Water Resources'),
         ('climate_koeppen', 'Csa'),
         ('doi', '10.17190/AMF/1246149'),
         ('doi_citation',
          'Dennis Baldocchi (2013-) AmeriFlux US-Tw3 Twitchell Alfalfa, 10.17190/AMF/1246149'),
         ('doi_dataproduct', 'AmeriFlux'),
         ('team_member_address',
          'Department of Environmental Science, Policy and Management, 137 Mulford Hall, 345 Hilgard Hall,Berkeley, CA USA 94720-3110'),
         ('url', 'http://nature.berkeley.edu/biometlab/sites.php?tab=US-Tw3'),
         ('dom_dist_mgmt', 'Agriculture'),
         ('site_snow_cover_days', '0.0'),
         ('state', 'CA'),
         ('location_date_start', '20130524.0'),
         ('acknowledgement',
          'Biometeorology Lab, University of California, Berkeley, PI:  Dennis Baldocchi')]


A useful method is the :meth:`configparser.ConfigParser.get` which takes the
section of the config file and the “option” and returns the value:

    >>> d.config.get(section='METADATA', option='site_name')
        'Twitchell Alfalfa'


    >>> # section and option are optional keywords
    >>> d.config.get('METADATA', 'site_name')
        'Twitchell Alfalfa'

.. Tip::
   If you are unsure if an entry or option exists in the config file, use the
   ``fallback`` keyword argument

    >>> # section and option are optional keywords
    >>> d.config.get('METADATA', 'site name', fallback='na')
        'na'

Some metadata entries are added as :obj:`.Data` attributes for easier access
as they are used in multiple ways later, these include:

-  site_id\ :math:`^*`
-  elevation\ :math:`^*`
-  latitude\ :math:`^*`
-  longitude\ :math:`^*`
-  na_val
-  qc_threshold
-  qc_flag

:math:`^*`\ mandatory **METADATA** entries in the config file, see
:ref:`Setting up a Config File` for further explanation.

View all the columns as found in the header row of the input time series
climate file.

    >>> d.header
        array(['TIMESTAMP_START', 'TIMESTAMP_END', 'CO2', 'H2O', 'CH4', 'FC',
               'FCH4', 'FC_SSITC_TEST', 'FCH4_SSITC_TEST', 'G', 'H', 'LE',
               'H_SSITC_TEST', 'LE_SSITC_TEST', 'WD', 'WS', 'USTAR', 'ZL', 'TAU',
               'MO_LENGTH', 'V_SIGMA', 'W_SIGMA', 'TAU_SSITC_TEST', 'PA', 'RH',
               'TA', 'VPD_PI', 'T_SONIC', 'T_SONIC_SIGMA', 'SWC_1_1_1',
               'SWC_1_2_1', 'TS_1_1_1', 'TS_1_2_1', 'TS_1_3_1', 'TS_1_4_1',
               'TS_1_5_1', 'NETRAD', 'PPFD_DIF', 'PPFD_IN', 'PPFD_OUT', 'SW_IN',
               'SW_OUT', 'LW_IN', 'LW_OUT', 'P', 'FC_PI_F', 'RECO_PI_F',
               'GPP_PI_F', 'H_PI_F', 'LE_PI_F'], dtype='<U15')


.. Note::
   All of the header columns will not necessarily be loaded, only those
   specified in the config file. Also, no data other than the header line is
   loaded into memory when creating a :obj:`.Data` object, the time series data
   is only loaded when calling :attr:`.Data.df` for increased efficiency for
   some workflows involving only metadata.

Variable names and units
^^^^^^^^^^^^^^^^^^^^^^^^

In ``flux-data-qaqc`` there are two naming schemes for climate
variables, the names as defined by the column headers in the input time
series file and the internal names for some variables and calculated
variables created by the package. We will refer to these two sets as
“user-defined” and “internal” names hereforth.

The :attr:`.Data.variables` attribute maps the internal to user-defined
variable names:

    >>> d.variables
        {'date': 'TIMESTAMP_START',
         'Rn': 'NETRAD',
         'G': 'G',
         'LE': 'LE_PI_F',
         'H': 'H_PI_F',
         'sw_in': 'SW_IN',
         'sw_out': 'SW_OUT',
         'lw_in': 'LW_IN',
         'lw_out': 'LW_OUT',
         'vpd': 'VPD_PI',
         't_avg': 'T_SONIC',
         'ws': 'WS',
         'theta_1': 'SWC_1_1_1',
         'theta_2': 'SWC_1_2_1'}


And, the :attr:`.Data.inv_map` maps the internal to user-defined names if
they differ, however this is only created once the data is loaded by
calling :attr:`.Data.df`.

    >>> # a similar dictionary attribute for input units
    >>> d.units
        {'Rn': 'w/m2',
         'G': 'w/m2',
         'LE': 'w/m2',
         'H': 'w/m2',
         'sw_in': 'w/m2',
         'sw_out': 'w/m2',
         'lw_in': 'w/m2',
         'lw_out': 'w/m2',
         'vpd': 'hPa',
         't_avg': 'C',
         'ws': 'm/s',
         'theta_1': '(%): Soil water content (volumetric), range 0-100',
         'theta_2': '(%): Soil water content (volumetric), range 0-100'}


Accessing input data
^^^^^^^^^^^^^^^^^^^^

The :py:attr:`Data.df` property gves access to the time series input climate
data for columns specified in the config file as a datetime-indexed
:obj:`pandas.DataFrame` object. This object has numerous powerful built in
tools for time series analysis and visualization.

    >>> # first 5 datetimes that are not gaps
    >>> d.df.dropna().head()

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
          <th>input_G</th>
          <th>WS</th>
          <th>VPD_PI</th>
          <th>T_SONIC</th>
          <th>SWC_1_1_1</th>
          <th>SWC_1_2_1</th>
          <th>NETRAD</th>
          <th>SW_IN</th>
          <th>SW_OUT</th>
          <th>LW_IN</th>
          <th>LW_OUT</th>
          <th>H_PI_F</th>
          <th>LE_PI_F</th>
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
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>2013-05-24 12:30:00</th>
          <td>122.194848</td>
          <td>3.352754</td>
          <td>18.853678</td>
          <td>25.682739</td>
          <td>6.6790</td>
          <td>26.1655</td>
          <td>652.648719</td>
          <td>1027.756939</td>
          <td>212.800000</td>
          <td>300.363524</td>
          <td>462.671744</td>
          <td>95.487930</td>
          <td>375.841436</td>
          <td>16.42225</td>
        </tr>
        <tr>
          <th>2013-05-24 13:00:00</th>
          <td>108.054863</td>
          <td>3.882154</td>
          <td>18.560999</td>
          <td>26.057700</td>
          <td>6.7065</td>
          <td>26.1600</td>
          <td>629.990486</td>
          <td>997.749437</td>
          <td>209.933333</td>
          <td>303.269447</td>
          <td>461.095065</td>
          <td>96.584383</td>
          <td>371.619775</td>
          <td>16.43325</td>
        </tr>
        <tr>
          <th>2013-05-24 13:30:00</th>
          <td>79.330662</td>
          <td>4.646089</td>
          <td>18.900260</td>
          <td>26.067374</td>
          <td>6.7120</td>
          <td>26.1545</td>
          <td>595.817687</td>
          <td>954.988747</td>
          <td>206.733333</td>
          <td>303.017852</td>
          <td>455.455579</td>
          <td>84.066406</td>
          <td>358.194935</td>
          <td>16.43325</td>
        </tr>
        <tr>
          <th>2013-05-24 14:00:00</th>
          <td>52.366527</td>
          <td>5.048825</td>
          <td>20.440061</td>
          <td>25.961307</td>
          <td>6.7395</td>
          <td>26.1325</td>
          <td>549.039365</td>
          <td>900.975244</td>
          <td>201.333333</td>
          <td>298.914731</td>
          <td>449.517276</td>
          <td>69.449710</td>
          <td>406.528564</td>
          <td>16.43600</td>
        </tr>
        <tr>
          <th>2013-05-24 14:30:00</th>
          <td>35.658417</td>
          <td>5.302946</td>
          <td>21.064824</td>
          <td>25.954462</td>
          <td>6.7450</td>
          <td>26.1215</td>
          <td>493.519695</td>
          <td>833.458365</td>
          <td>192.066667</td>
          <td>296.791541</td>
          <td>444.663544</td>
          <td>47.774030</td>
          <td>315.295309</td>
          <td>16.43325</td>
        </tr>
      </tbody>
    </table>
    </div>
    <br />



.. Tip:: 
   There are *many* tutorials on how to use the :obj:`pandas.DataFrame` and its
   powerful data analysis tools for multiple purposes online, to get started
   you may want to visit Panda’s own list of tutorials `here
   <https://pandas.pydata.org/pandas-docs/stable/getting_started/tutorials.html#internal-guides>`__.

By default the column names in :attr:`.Data.df` are retained from
user-defined names unless they were named exactly the same as an
internal name. For example the input ground heat flux column in this
dataset is named “G”, therefore it was renamed as “input_g”

    >>> d.df.columns
        Index(['input_G', 'WS', 'VPD_PI', 'T_SONIC', 'SWC_1_1_1', 'SWC_1_2_1',
               'NETRAD', 'SW_IN', 'SW_OUT', 'LW_IN', 'LW_OUT', 'H_PI_F', 'LE_PI_F',
               'theta_mean'],
              dtype='object')


    >>> # the new name was also updated in Data.variables
    >>> d.variables.get('G')
        'input_G'



As stated earlier, :attr:`.Data.inv_map` maps the user-defined names to
internal ``flux-data-qaqc`` names only after loading :attr:`.Data.df`:

    >>> d.inv_map
        {'TIMESTAMP_START': 'date',
         'NETRAD': 'Rn',
         'input_G': 'G',
         'LE_PI_F': 'LE',
         'H_PI_F': 'H',
         'SW_IN': 'sw_in',
         'SW_OUT': 'sw_out',
         'LW_IN': 'lw_in',
         'LW_OUT': 'lw_out',
         'VPD_PI': 'vpd',
         'T_SONIC': 't_avg',
         'WS': 'ws',
         'SWC_1_1_1': 'theta_1',
         'SWC_1_2_1': 'theta_2'}



.. Tip:: 
   The :attr:`.Data.inv_map` is mainly used to rename the dataframe to internal
   names, this can be very useful if you are creating your own custom workflows
   using the ``flux-data-qaqc`` API because it allows you to only know the
   internal names of variables therefore they can be hard coded into your
   workflow and applied to different eddy covariance datasets. For example,
   let’s say we wanted to make HTML tables of basic statistics of just the
   energy balance components for many datasets (that may have different names
   for the same variables) and save the file using the user-defined names:

   >>> d = Data('US-Tw3_config.ini')
   >>> df = d.df.rename(columns=d.inv_map)
   >>> # get some metadata for saving
   >>> site_id = d.site_id
   >>> vars_we_want = ['H', 'LE', 'Rn', 'G']
   >>> # rename variables, calculate basice statistics table and save to HTML
   >>> df[vars_we_want].rename(columns=d.variables).describe().to_html('{}.html'.format(site_id))
       Calculating mean for var: THETA from columns: ['SWC_1_1_1', 'SWC_1_2_1']
       WARNING: renaming column G to input_G
   >>> # which produces the following HTML table with user-defined names:
   >>> from IPython.display import HTML
   >>> HTML(filename='{}.html'.format(site_id))


.. raw:: html
    :file: _static/tutorial/US-Tw3.html

.. raw:: html

       <br />

Another powerful feature of the :attr:`.Data.df` property is that it is
datetime-indexed using the input data’s temporal frequency, view the
date index like so:

    >>> d.df.index
        DatetimeIndex(['2013-01-01 00:00:00', '2013-01-01 00:30:00',
                       '2013-01-01 01:00:00', '2013-01-01 01:30:00',
                       '2013-01-01 02:00:00', '2013-01-01 02:30:00',
                       '2013-01-01 03:00:00', '2013-01-01 03:30:00',
                       '2013-01-01 04:00:00', '2013-01-01 04:30:00',
                       ...
                       '2018-06-04 19:00:00', '2018-06-04 19:30:00',
                       '2018-06-04 20:00:00', '2018-06-04 20:30:00',
                       '2018-06-04 21:00:00', '2018-06-04 21:30:00',
                       '2018-06-04 22:00:00', '2018-06-04 22:30:00',
                       '2018-06-04 23:00:00', '2018-06-04 23:30:00'],
                      dtype='datetime64[ns]', name='date', length=95088, freq=None)
    


Datetime-indexed :obj:`pandas.DataFrame` objects have useful features for
time series analysis like grouping and calculating statistics by time
aggregates. The example below shows how to calculate the day of year
mean for energy balance components, it also demonstrates how to use the
``add_lines`` plotting method available to :obj:`.Data`, :obj:`.QaQc`, and
:obj:`.Plot` objects.

    >>> # convert to internal names, copy dataframe
    >>> df = d.df.rename(columns=d.inv_map)
    >>> # day of year mean of input energy balance components
    >>> vars_we_want = ['H', 'LE', 'Rn', 'G']
    >>> doy_means = df[vars_we_want].groupby(d.df.index.dayofyear).mean()
    >>> # create a Bokeh figure
    >>> fig = figure(x_axis_label='day of year', y_axis_label='day of year mean (w/m2)')
    >>> # arguements needed for creating interactive plots
    >>> plt_vars = vars_we_want
    >>> colors = ['red', 'blue', 'black', 'green']
    >>> x_name = 'date'
    >>> source = ColumnDataSource(doy_means)
    >>> Plot.add_lines(fig, doy_means, plt_vars, colors, x_name, source, labels=vars_we_want,
    >>>     x_axis_type=None) 
    >>> show(fig)

.. raw:: html
    :file: _static/tutorial/doy_mean_example.html
    

.. Note::
   The ``x_axis_type=None`` is a unique argument to :meth:`.Plot.add_lines` and
   :meth:`.Plot.line_plot` that in this case means to not try to force the
   x-axis format to a datetime representation, default is
   ``x_axis_type='date'``.

.. seealso::
   Some routines occur automatically when creating a :obj:`.Data` object,
   including calcuation of weighted and non-weighted averages of soil heat flux
   and soil moisture which is described in :ref:`Averaging data from multiple
   sensors`.

Modifying input data
^^^^^^^^^^^^^^^^^^^^

A last note on the :obj:`.Data` object (same goes for the :obj:`.QaQc` object)
is that :attr:`Data.df` is a class property, in this case that means that it
can be reassigned with a different :obj:`pandas.DataFrame`. This is
critical for manual pre-filtering and validation of data before
proceeding with energy balance closure routines. A simple example is
shown here:

    >>> # add 5 to air temperature, this would effect ET calculations later
    >>> x = d.df
    >>> x['T_SONIC'] += 5
    >>> d.df = x

A realistic use of the reassignability of the :attr:`.Data.df` and
:attr:`.QaQc.df` properties is shown in `manual cleaning of poor quality
data <https://flux-data-qaqc.readthedocs.io/en/latest/closure_explanation.html#step-0-manual-cleaning-of-poor-quality-data>`__.

.. seealso::
   The :meth:`.Data.apply_qc_flags` method allows for reading in quality
   control flags with the input data and filtering specific data out based on
   user-defined numeric or character flags. This routine is specific to
   :obj:`.Data` and includes several attributes that are added to a
   :obj:`.Data` instance, for full explanation and examples see
   :ref:`Quality-based data filtering`.

Visualize input data
--------------------

The :meth:`.Data.plot` method create a series of interactive time series
plots of input data, potential plots inlcude:

-  energy balance components
-  radiation components
-  multiple soil heat flux measurements
-  air temperature
-  vapor pressure and vapor pressure deficit
-  wind speed
-  precipitation
-  latent energy
-  multiple soil moisture measurements

If any of these variables are not found the plot(s) will not be added.

The most useful interactive features of plots created by
``flux-data-qaqc`` are:

-  pan/zoom
-  hover tolltips on var names, values, date
-  linked x-axes on time series plots
-  save plot option (can save specific subplot zoomed in)

Here is an example,

    >>> d.plot(output_type='notebook', plot_width=700)

The output plot is not shown in the online documentation due to 
memory constraints. 


.. hint:: 
   The plot methods of :obj:`.Data` and :obj:`.QaQc` objects have the keyword
   argument ``output_type`` which by default is set to “save”, the other two
   options are “notebook” for showing within a Jupyter Notebook and “show”
   which opens a temporary file in the default web browser.

If you rather save the plot, and maybe you want 2 columns of plots,

    >>> d.plot(ncols=2, plot_width=500) 

After saving a plot without specifying the output file path (keyword
argument ``out_file``), it will be saved to an “output” directory where
the config file is with the file name based on :attr:`.Data.site_id` with the
suffix “\_input_plots”:

    >>> # where the plot file was saved by default
    >>> d.plot_file
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/output/US-Tw3_input_plots.html')

The following plot is not shown due to excessive memory usage needed to build 
online documentation.

    >>> # view outplot plots within Jupyter notebook
    >>> from IPython.display import HTML
    >>> HTML(filename=d.plot_file)


.. hint:: 
   The :meth:`.QaQc.plot` method shown below is similar however it may include
   added plots with calculated and corrected variables (if they exist) and will
   always plot data in daily and monthly temporal frequency because daily
   frequency is required before applying ``flux-data-qaqc`` energy balance
   closure corrections.

Temporal resampling
-------------------

The :obj:`.QaQc` object holds several tools for managing data and eddy
covariance data analysis, but one of it’s primary features is temporal
resampling of input data to daily and monthly frequencies. The
resampling of time series data to daily frequency occurs upon the
creation of a :obj:`.QaQc` instance if the frequency within the preceeding
:obj:`.Data` object is not already daily:

    >>> # the frequency of the input data is 30 minute
    >>> d.df.index[0:5]
        DatetimeIndex(['2013-01-01 00:00:00', '2013-01-01 00:30:00',
                       '2013-01-01 01:00:00', '2013-01-01 01:30:00',
                       '2013-01-01 02:00:00'],
                      dtype='datetime64[ns]', name='date', freq=None)

    >>> # creating a QaQc instance will automatically convert to daily
    >>> q = QaQc(d)
        The input data temporal frequency appears to be less than daily.
        Data is being resampled to daily temporal frequency.
        Filtering days with less then 100.0% or 48/48 sub-daily measurements
        Converting vpd from hpa to kpa


    >>> # first 5 datetime indices are dates now
    >>> q.df.index[0:5]
        DatetimeIndex(['2013-01-01', '2013-01-02', '2013-01-03', '2013-01-04',
                       '2013-01-05'],
                      dtype='datetime64[ns]', name='date', freq=None)



The method used for aggregating different variables, e.g. mean or sum,
when resampling to daily or monthly frequency is defined in the
:attr:`QaQc.agg_dict` class attribute:

    >>> # these are the internal names as keys and temporal aggregation method as values
    >>> QaQc.agg_dict
        {'energy': 'mean',
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
         'gridMET_prcp': 'sum',
         'lw_in': 'mean',
         't_avg': 'mean',
         'rso': 'mean',
         'sw_pot': 'mean',
         'sw_in': 'mean',
         'vp': 'mean',
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
         'H_user_corr': 'mean'}



.. note:: 
   There are several calculated variables above that may not look familiar,
   many are calculated by the energy balance closure correction routines and
   described in :ref:`Closure Methodologies`.  Also, any other variables (not
   found in :attr:`.QaQc.agg_dict` that exist in a :attr:`.QaQc.df` before
   accessing :attr:`.QaQc.monthly_df` the first time will be averaged in the
   monthly time series dataframe (:attr:`.QaQc.monthly_df`).

The :obj:`.QaQc` constructor tries to infer the temporal frequency of the
input time series data, however the method is not always accurate, to
access the inferred initial temporal frequency of the data view the
:attr:`.QaQc.temporal_freq` attribute:

    >>> q.temporal_freq
        '30T'



If the inferred input frequency was accurate you will see a `Pandas datetime
alias
<https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__,
in this case ‘30T’ is thirty minutes. If the temporal frequency is not automatically detected you should be able to rely on the ``n_samples_per_day`` instance attribute that is manually estimated by the ``QaQc`` constructor:

    >>> q.n_samples_per_day
        48

Filter days with sub-daily gaps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``drop_gaps`` and ``daily_frac`` keyword arguments used when creating a :obj:`.QaQc` instance allow you to control how days with sub-daily measurement gaps will or will not be filtered out when resampling to daily frequency.

Sub-daily gaps in energy balance variables :math:`LE`, :math:`H`, :math:`Rn`, and :math:`G` , and daily ASCE standaridized reference ET inputs, e.g. hourly :math:`ea` ("vp"), :math:`rs` ("sw_in"), :math:`t_min`, :math:`t_max`, and :math:`ws`, can be linearly interpolated automatically before daily aggregations. Interpolation is performed over gap lengths measured in hours, with options to control the longest length of gap to interpolate when :math:`Rn \ge 0` controlled by the :obj:`.QaQc` keyword argument ``max_interp_hours`` (default 2 hours) and the longest gap to interpolate when :math:`Rn < 0` set by the ``max_interp_hours_night`` (default 4 hours).:math:`

.. Important:: 

   By default the :obj:`.QaQc` constructor will first linearly interpolate
   energy balance and ASCE ref. ET variables (:math:`LE`, :math:`H`,
   :math:`Rn`, :math:`G`, :math:`ea` ("vp"), :math:`rs` ("sw_in"),
   :math:`t_min`, :math:`t_max`, and :math:`ws`) according to the maximum gap
   lengths (``max_interp_hours`` and ``max_interp_hours_night``) and then count
   sub-daily gaps and drop days (set values to null) for all climate data
   columns (not QC flag or sub-daily gap count columns) where any of the
   sub-daily data are missing because by default ``drop_gaps=True`` and
   ``daily_frac=1.0``. In other words, if you have hourly input data
   for(:math:`LE` and one hour was missing on a given day, by default that hour
   will be linearly interpolated before calculating the daily time series and
   the daily mean will be calculated after. On the other hand, if other climate
   variables had a single hour missing on a given day, e.g. wind direction,
   this day would be filtered out by the :obj:`.QaQc` constructor.  This is
   important because the daily time series is what is used in all energy
   balance closure correction and daily ASCE standardized reference ET
   algorithms.

The percentage of sub-daily samples to require set by the ``daily_frac`` argument and the maximum length of gaps to linearly interpolate set by ``max_interp_hours`` and ``max_interp_hours_night`` complement each other and are used in tandem. For example, if the input data is half-hourly and you only want a maximum of 4 hours to be interpolated on any given day and gap lengths to interpolate should be no more than 2 hours each then you would pass the following parameters to the :obj:`.QaQc` constructor:


    >>> q = QaQc(d, daily_frac=20/24, max_interp_hours=2, max_interp_hours_night=2)
        The input data temporal frequency appears to be less than daily.
        Data is being resampled to daily temporal frequency.
        Linearly interpolating gaps in energy balance components up to 2 hours when Rn < 0 and up to 2 hours when Rn >= 0.
        Filtering days with less then 83.33333333333334% or 40/48 sub-daily measurements    

In this case we set ``daily_frac=20/24`` because we are only allowing a maximum of 4 hours of total gaps in the day in other words we are requiring 40 of the 48 half hourly samples to exist before we filter out a day. Remember, because linear interpolation of gaps is done before counting sub-daily gaps, this could result in retaining days with more than 4 hours of gaps in the original time series of energy balance components. You may also pass the ``daily_frac`` arugment as a decimal fraction, e.g. :math:`0.8333 \approx 20/24`.

To not drop any days and take daily means/sums based on whatever data exists in a given day *without* any interpolation of energy balance variables,

    >>> q = QaQc(d, drop_gaps=False, max_interp_hours=None)
        The input data temporal frequency appears to be less than daily.
        Data is being resampled to daily temporal frequency.


Let’s view a comparison of :math:`Rn` using different options of
filtering days with sub-daily gaps in the working dataset, because it
has several periods of systematic gaps which cause upwards skewing of
daily mean :math:`Rn` if not filtered carefully:

    >>> # make an empty pandas dataframe for Rn series
    >>> Rn_df = pd.DataFrame()
    >>> # recreate multiplt QaQc instances using different sub-day gap filters
    >>> q = QaQc(d, drop_gaps=False, max_interp_hours=None)
    >>> Rn_df['sub_day_gaps'] = q.df.Rn_subday_gaps
    >>> Rn_df['no_filter_no_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, drop_gaps=False)
    >>> Rn_df['no_filter_with_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=0.5) # filter days with less than 50% data
    >>> Rn_df['require_50'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=0.75)
    >>> Rn_df['require_75'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=1, max_interp_hours=24, max_interp_hours_night=24) 
    >>> Rn_df['require_100_with_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=1, max_interp_hours=None) 
    >>> Rn_df['require_100_no_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> # plot to compare results of day-gap filter
    >>> fig = figure(x_axis_label='date', y_axis_label='mean daily net radiation (w/m2), filtered based on sub-daily gaps')
    >>> # arguments needed for creating interactive line plots
    >>> colors = ['red', 'darkred','orange', 'blue', 'black', 'tan']
    >>> plt_vars = ['no_filter_no_interp', 'no_filter_with_interp', 'require_50', 'require_75', 'require_100_with_interp', 'require_100_no_interp']
    >>> labels = ['no filter wout/interp.', 'no filter w/interp.', 'require > 50% w/interp.', 'require > 75% w/interp.', 'require 100% w/interp.', 'require 100% wout/interp.']
    >>> x_name = 'date'
    >>> source = ColumnDataSource(Rn_df)
    >>> Plot.add_lines(fig, Rn_df, plt_vars, colors, x_name, source, labels=labels) 
    >>> # add daily gap counts to secondary y
    >>> fig.extra_y_ranges['gap_counts'] = Range1d(start=0, end=48)
    >>> fig.add_layout(LinearAxis(y_range_name='gap_counts', axis_label='number of sub-daily gaps'), 'right')
    >>> fig.circle('date', 'sub_day_gaps', legend='n sub-day gaps', y_range_name='gap_counts',
    >>>     color='silver', source=source
    >>> )
    >>> fig.hover[0].tooltips.append(('sub_day_gaps','@{}'.format('sub_day_gaps')))
    >>> fig.legend.location = 'top_right'
    >>> show(fig)


.. raw:: html 
    :file: _static/tutorial/filter_subday.html 

Try zooming in on the gaps filled by the "no filter wout/interp." line to compare which days are retained/filtered by different options, also remove lines by clicking on them in the legend to compare subsets of options.

.. Tip:: 
   For a more fine-grained approach to filtering out days where perhaps
   multiple 2 hour gaps were filled use the newly created daily gap count
   columns: "LE_subday_gaps", "H_subday_gaps", "Rn_subday_gaps", and
   "G_subday_gaps":

      >>> q = QaQc(d)
      >>> df = q.df.rename(columns=q.inv_map)

   For example, you could post-filter out days in any given energy balance
   variable, in this case :math:`Rn` where sub-daily gaps exceed a threshold:

      >>> df.loc[(df.Rn_subday_gaps > 4) & (df.Rn.notna()), ['Rn','Rn_subday_gaps']]

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
                  <th>Rn</th>
                  <th>Rn_subday_gaps</th>
                </tr>
                <tr>
                  <th>date</th>
                  <th></th>
                  <th></th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>2015-06-09</td>
                  <td>101.710194</td>
                  <td>5.0</td>
                </tr>
                <tr>
                  <td>2015-11-20</td>
                  <td>47.990988</td>
                  <td>5.0</td>
                </tr>
                <tr>
                  <td>2016-01-15</td>
                  <td>72.495973</td>
                  <td>8.0</td>
                </tr>
                <tr>
                  <td>2018-01-06</td>
                  <td>79.507008</td>
                  <td>7.0</td>
                </tr>
                <tr>
                  <td>2018-05-10</td>
                  <td>160.997332</td>
                  <td>6.0</td>
                </tr>
              </tbody>
            </table>
            </div>
          </br>

Monthly time series
^^^^^^^^^^^^^^^^^^^

The :attr:`.QaQc.monthly_df` property allows for creating the monthly time
series of input anc calculated variables provided by
:meth:`.QaQc.correct_data`. It uses the same temporal aggregation methods as
the daily time series i.e. from :attr:`.QaQc.agg_dict`. Although there are
many similarities there are important differences between :attr:`.QaQc.df`
and :attr:`.QaQc.monthly_df` other than the obvious: when accessing the
:attr:`.QaQc.monthly_df` it will automatically run the default energy balance
closure correction routine provided by :meth:`.QaQc.correct_data` *if* it has
not yet been run. You can check if it has been run at anytime by:

    >>> q.corrected
        False

To show how this works let’s access the monthly data and show the
monthly statistics of the “corrected” evapotranspiration (ET_corr):


    >>> # first note, ET_corr is not in the dataset yet
    >>> 'ET_corr' in q.df.columns
        False

Now access the monthly time series,

    >>> q.monthly_df;
    >>> 'ET_corr' in q.df.columns
        True

By calling the monthly dataframe, the energy balance closure was applied
automatically

    >>> q.monthly_df.ET_corr.describe()
        count     61.000000
        mean      87.858135
        std       49.938287
        min       11.370062
        25%       41.418994
        50%       84.383190
        75%      127.500125
        max      192.033481
        Name: ET_corr, dtype: float64


    >>> q.corrected
        True

.. note:: 
   The :attr:`.QaQc.monthly_df` also filters out months with less than 30% of
   days of the month missing by default. To calculate monthly time series with
   other threshold fractions of days required use the
   :func:`.util.monthly_resample` function and adjust the keyword argument
   ``thresh`` which is the fraction (0-1) of days of the month required to not
   be gaps otherwise the month’s value will be forced to null, e.g. if you
   wanted to caclulate the monthly mean air temperature requiring 30 and 90
   percent of the days in the month to not be gaps:

    >>> from fluxdataqaqc.util import monthly_resample
    >>> # select just t_avg for example
    >>> cols = ['t_avg'] 
    >>> df = q.df.rename(columns=q.inv_map)
    >>> # create temporary df with different monthly resample results
    >>> tmp_df = monthly_resample(df, cols, 'mean', thresh=0.9).rename(
    >>>     columns={'t_avg': 'thresh_90'}
    >>> )
    >>> # join temp dataframe with monthly resample results using different thresh
    >>> monthly_gap_comp = tmp_df.join(monthly_resample(df, cols, 'mean', thresh=0.3).rename(
    >>>     columns={'t_avg': 'thresh_30'})
    >>> )
    >>> # plot to compare results of day-gap filter
    >>> fig = figure(x_axis_label='date', y_axis_label='monthy mean air temperature (C), filtered based on daily gaps')
    >>> # arguments needed for creating interactive line plots
    >>> x = 'date'
    >>> source = ColumnDataSource(monthly_gap_comp)
    >>> # this example also shows how to use other Bokeh plot arguments
    >>> Plot.line_plot(fig,'date','thresh_30',source,'red',label='require > 30%', line_alpha=0.5) 
    >>> Plot.line_plot(fig,'date','thresh_90',source,'black',label='require > 90%',line_dash='dotted', line_width=2) 
    >>> fig.legend.location = 'top_right'
    >>> show(fig)

.. raw:: html 
    :file: _static/tutorial/filter_monthly.html

.. raw:: html 

   <br>

Energy balance corrections
--------------------------

``flux-data-qaqc`` provides routines that adjust surface energy balance fluxes
to improve energy balance closure of eddy covariance flux station data.
These routines ultimately result in a corrected daily and monthly time series
of latent energy, sensible heat, and evapotranspiration with the option to
gap-fill days in corrected ET with ET calculated from gridMET reference ET and
fraction of reference ET.

There are two methods currently implemented: 

* Energy Balance Ratio method (default), modified from the `FLUXNET method <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
* Bowen Ratio approach (forces closure) 
* Multiple least squares regression
  - user defines LHS and RHS from :math:`LE`, :math:`H`, :math:`Rn`, and :math:`G`,

Detailed descriptions of methods including daily ET gap-filling methods
can be found in the online documentation :ref:`Closure Methodologies`
page. A few important notes on the API of these methods and other
hydro-climatic statistical variables that are calculated are shown
below.

The :meth:`.QaQc.correct_data` method is used to run energy balance closure
corrections. Here are a few tips on using them,

    >>> # note above the monthly_df automatically applied the 'ebr' Energy Balance Ratio correction
    >>> q.corr_meth
        'ebr'

    >>> # potential correction options
    >>> q.corr_methods
        ('ebr', 'br', 'lin_regress')


    >>> # to specify the Bowen Raito method:
    >>> q.correct_data(meth='br')

    >>> # the most recently used correction method is now shown
    >>> q.corr_meth
        'br'

.. Tip:: 
   After applying any energy balance closure correction routine all previous
   corrected variables will be overwritten or dropped in :attr:`.QaQc.df`,
   :attr:`.QaQc.monthly_df`, and :attr:`.QaQc.variables`, therefore to make a comparison
   of different methods on the same data make a copy of the ``df`` or
   ``monthly_df`` properties before running the next correction, e.g.

    >>> # make copies of daily results of different correction options
    >>> q.correct_data(meth='ebr')
    >>> ebr_gapfilled = q.df
    >>> q.correct_data(meth='ebr', etr_gap_fill=False)
    >>> ebr_notgapfilled = q.df
    >>> q.correct_data(meth='br')
    >>> br_gapfilled = q.df
    >>> q.correct_data(meth='br', etr_gap_fill=False)
    >>> br_notgapfilled = q.df


ET gap-filling
^^^^^^^^^^^^^^

A few notes on the option that uses reference ET and fraction of daily
reference ET to fill in large gaps in corrected ET, i.e. the keyword
argument ``QaQc.correct_data(etr_gap_fill = True)``.

-  The nearest `gridMET <http://www.climatologylab.org/gridmet.html>`__
   cell’s time series data for precipitation and alfalfa reference ET is
   attempted to be downloaded if it is not found in the
   ``gridmet_file_path`` entry of the config.ini file.

-  If the path to a gridMET file is not found it is re-downloaded, the
   config file will be updated with the new path and resaved.

-  Only the overlapping time period that matches the eddy covariance
   time series data is attempted to be downloaded, i.e. the period in
   ``QaQc.df.index``.

-  When a gridMET file is downloaded it will always be saved in a
   subdirectory where the config file is located called “gridMET_data”
   and named using the :attr:`QaQc.site_id` and gridMET cell centroid
   latitude and longitude.

-  Corrected latent energy (:math:`LE_{corr}`) gaps are also backwards
   filled from gap-filled ET.

.. caution:: 
   `gridMET <http://www.climatologylab.org/gridmet.html>`__ only exists within
   the contiguous United States and from 1979 to present, therefore if your
   station lies outside of this region or you are analyzing eddy flux data
   recorded before 1979 this option will not be ususable and you should always
   run corrections with ``etr_gap_fill=False`` to avoid potential errors.


The Bowen Ratio correction method will produce the ‘br’ variable which
is the Bowen Ratio.

Other calculations
------------------

By default, :meth:`.QaQc.correct_data` also calculates ET from input latent
energy (LE) and air temperature, corrected ET from corrected LE and air
temperature, potential clear sky radiation (ASCE formulation), and the
:obj:`.Data` object attempts to calculate vapor pressure deficit from vapor
pressure and air temperature or vapor pressure from vapor pressure
deficit and air temperature if they exist at hourly or shorter temporal
frequency.

Evapotranspiration
^^^^^^^^^^^^^^^^^^

The evapotranspiration (ET) calculations are described in :ref:`Steps 7 and 8 correct turbulent fluxes, EBR, and ET` of the Energy Balance Ratio correction explanation.

ASCE clear sky radiation
^^^^^^^^^^^^^^^^^^^^^^^^

Daily ASCE potential clear sky radiation (:math:`R_{so}`) is calculated
using equation 19 in the “ASCE Standardized Reference Evapotranspiration
Equation” final report by the Task Committee on Standardization of
Reference Evapotranspiration Environmental and Water Resources Institute
of the American Society of Civil Engineers January, 2005
`here <https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf>`__.
This calculation is a simple method based primarily on elevation and
latitude which results in a theoretical envelope of :math:`R_{so}` as a
function of the day of year,

.. math::  R_{so} = \left(5 + 2 \times 10^{-5} z \right) R_a 

where :math:`z` is elevation in meters and :math:`R_a` is daily
extraterrestrial radiation (radiation with in the absence of an
atmosphere), which itself is a well-behaved function of solar
declination, the day of the year and the solar constant (see equations
21-29 in the `ASCE
report <https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf>`__).

Vapor pressure/deficit
^^^^^^^^^^^^^^^^^^^^^^

The :obj:`.Data` object will attempt to calculate vapor pressure or vapor
pressure deficit if one exists but not the other and average air
temperature time series also exists with the input data at hourly or
shorter temporal frequency. The Tetens equation (eqn. 37 in the ASCE report) 
s an accurate approximation for saturation vapor pressure (:math:`es`) in kPa as a function of air temperature,

.. math::  es = 0.6108  e^{\left(\frac{17.27 \cdot T}{(T + 237.3)}\right)} 

where :math:`T` is average hourly air temperature in degrees celcius.
Vapor pressure deficit (:math:`vpd`) is,

.. math::  vpd = es - ea,

where :math:`ea` is actual vapor pressure in kPa. **Note,** The
equations above are defined for hourly measurements however they are used for
hourly or shorter mean variables (:math:`T`, :math:`ea`, or :math:`vpd`)
within ``flux-data-qaqc`` and then converted to daily means, if they are
not present in the input data at hourly or shorter frequencies then they
are not calculated.

These equations can be rearanged to solve for either :math:`es` or
:math:`vpd` given the other variable and air temperature. For example,
if given :math:`T` and :math:`vpd`, then to get actual vapor pressure

.. math::  es = 0.6108  e^{\left(\frac{17.27 \cdot T}{(T + 237.3)}\right)}  

.. math::  ea = es - vpd. 

In ``flux-data-qaqc`` actual vapor pressure is named "vp" not "ea". Also, during these calculations, if relative humidity is not found in the input dataset then it will subsequently be estimated as 

.. math:: rh = 100 \times \frac{ea}{es}.

.. hint:: 
   The same calculations are available at the daily timestep but are not
   automatically applied as the hourly or higher temporal frequency calculation
   is preffered. To apply the estimates of vapor pressure or vapor pressure
   deficit, and saturation vapor pressure, and relative humidity with daily data
   one must call the :meth:`.QaQc._calc_vpd_from_vp` method from a :obj:`.QaQc`
   instance. 
   
Legend for calculated variable names
------------------------------------

Atmospheric and related variables created by energy balance closure corrections 
are described in :ref:`Closure Methodologies` and other calculations can be 
found throughout this documentation.  Below is a reference table with basic 
descriptions of all variables that are calculated or renamed by ``flux-data-qaqc`` 
using its standardized naming scheme.  The names found here (in the "Variable" 
column) are the default names that are saved to the header in both daily and 
monthly time series output CSV files as well as what are shown in the 
interactive plots upon hovering with a cursor.  Note, not all of the following 
variables will be estimated by ``flux-data-qaqc`` in all scenarios, the list of 
variables that may be found in output files is completely dependent on the 
input data provided in the configuration file and also dependend on which 
calculations are used.  For example, ASCE_ETo will only exist if the hourly or 
daily methods for computing ASCE standardized reference ET are used and the 
required input variables exist. 
    

 ================= ================================================================================ ============ 
 Variable           Description                                                                     Unit         
 ================= ================================================================================ ============ 
 ASCE_ETo           short/grass ASCE standardized reference ET                                      mm time-1    
 ASCE_ETr           tall/alfalfa ASCE standardized reference ET                                     mm time-1    
 br                 bowen ratio                                                                     unitless     
 ebc_cf             energy balance closure correction factor (inverse of ebr_corr)                  unitless     
 ebr                input energy balance ratio                                                      unitless     
 ebr_5day_clim      5 day climatology of the filtered Energy Balance Ratio                          unitless     
 ebr_corr           corrected energy balance ratio                                                  unitless     
 energy             input Rn - G                                                                    w m-2        
 es                 saturation vapor pressure                                                       kPa          
 ET                 ET calculated from input LE and average air temperature                         mm time-1    
 ET_corr            ET calculated from LE_corr                                                      mm time-1    
 ET_fill            gridMET_ETr * ETrF_filtered (fills gaps in ET_corr)                             mm time-1    
 ET_fill_val        value of ET_fill on gap days                                                    mm time-1    
 ET_gap             True on gap days in ET_corr, False otherwise                                    unitless     
 ET_user_corr       corrected ET, user-provided                                                     mm time-1    
 EToF               fraction of reference ET for ET_corr, i.e. ET_corr / gridMET_ETo                unitless     
 EToF_filtered      filtered and gap-filled EToF                                                    unitless     
 ETrF               fraction of reference ET for ET_corr, i.e. ET_corr / gridMET_ETr                unitless     
 ETrF_filtered      filtered and gap-filled EtrF                                                    unitless     
 flux               input LE + H                                                                    w m-2        
 flux_corr          LE_corr + H_corr                                                                w m-2        
 G                  average or single soil heat flux                                                w m-2        
 G_[1,2,3,…]        soil heat flux at sensor                                                        w m-2        
 G_subday_gaps      number of gaps in initial G per day                                             unitless     
 gridMET_ETo        gridMET short/grass reference ET (nearest cell)                                 mm time-1    
 gridMET_ETr        gridMET tall/alfalfa reference ET (nearest cell)                                mm time-1    
 gridMET_prcp       gridMET precipitation (nearest cell)                                            mm time-1    
 gridMET_[other]    other optional gridMET variables (nearest cell)                                 NA    
 H_corr             corrected sensible heat                                                         w m-2        
 H_subday_gaps      number of gaps in initial H per day                                             unitless     
 H_user_corr        corrected sensible heat, user-provided                                          w m-2        
 LE_corr            corrected latent energy                                                         w m-2        
 LE_subday_gaps     number of gaps in initial LE per day                                            unitless     
 LE_user_corr       corrected latent energy, user-provided                                          w m-2        
 lw_in              incoming longwave radiation                                                     w m-2        
 lw_out             outgoing longwave radiation                                                     w m-2        
 ppt                precipitation                                                                   mm time-1    
 rh                 relative humidity                                                               %            
 Rn                 net radiation                                                                   w m-2        
 Rn_subday_gaps     number of gaps in initial Rn per day                                            unitless     
 rso                clear sky radiation (ASCE formulation)                                          w m-2        
 sw_in              incoming shortwave radiation                                                    w m-2        
 sw_out             outgoing shortwave radiation                                                    w m-2        
 sw_pot             potential shortwave radiation, user-provided                                    w m-2        
 t_avg              average temperature                                                             C            
 t_dew              dew point temperature                                                           C            
 t_max              maximum temperature                                                             C            
 t_min              minimum temperature                                                             C            
 theta              average or single soil moisture                                                 user defined 
 theta_[1,2,3,…]    soil moisture at sensor                                                         user defined 
 vp                 actual vapor pressure                                                           kPa          
 vpd                vapor pressure deficit                                                          kPa          
 wd                 wind direction                                                                  user defined 
 ws                 wind speed                                                                      m s-1        
 ================= ================================================================================ ============ 


A note on units
---------------

Upon creation of a :obj:`.QaQc` object, variables are checked for valid input
units and converted to required units needed for internal calculations when
running :meth:`.QaQc.correct_data` and for certain default plots (see below).
For a list of valid input units for different variables refer to the
:attr:`.QaQc.allowable_units` attribute:

    >>> q.allowable_units
        {'LE': ['w/m2'],
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
         't_avg': ['c', 'f']}


For each variable above, if given one of the units allowable the units
will automatically be converted to the required units.

To know which variables are required to be in particular units view
:attr:`.Qc.required_units`:

    >>> q.required_units
        {'LE': 'w/m2',
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
         't_avg': 'c'}



.. note:: 
   The list of allowable units is a work in progress, if your input units are
   not available consider raising an issue on `GitHub
   <https://github.com/Open-ET/flux-data-qaqc/issues>`__ or providing the
   conversion directly with a pull request. Automatic unit conversions are
   handled within the :obj:`.util.Convert` class using the
   :meth:`.util.Convert.convert` class method.

Save resampled and corrected data
---------------------------------

The :meth:`.QaQc.write` method conveniently writes daily and monthly time
series of input and calculated variables to comma separated value (CSV)
files. If the :meth:`.QaQc.correct_data` method has not yet been run it will
be and the monthly time series will also be created using the default
parameters for the correction routine (Energy Balance Ratio method with
ETr-based gap filling).

The default output directory for time series files can be
accessed/changed by the ``out_dir`` attribute, if not changed it will be
located in the same directory of the config.ini file. The daily and
monthly time series file names will begin with the :attr:`.QaQc.site_id`
followed by “daily_data” or “monthly_data” resepctively. For example,

    >>> # new QaQc instance
    >>> q = QaQc(d)
    >>> # a platform dependent pathlib.Path object
    >>> q.out_dir
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/output')

The line below shows that no output files have been written to
:attr:`.QaQc.out_dir` yet,

    >>> # print files in output directory that begin with the site_id
    >>> [f.name for f in q.out_dir.glob('{}*'.format(q.site_id))]
        ['US-Tw3_input_plots.html']

    >>> q.corrected
        False

    >>> # writing files also ran corrections since they were not yet run
    >>> q.write()
    >>> q.corrected
        True


Now the respective daily and monthly time series have been written to
:attr:`.QaQc.out_dir`,

    >>> [f.name for f in q.out_dir.glob('{}*'.format(q.site_id))]
        ['US-Tw3_daily_data.csv', 'US-Tw3_monthly_data.csv', 'US-Tw3_input_plots.html']


.. hint:: 
   You can overwrite the default name of the output directory to save the daily
   and monthly time series using the ``out_dir`` keyword argument to
   :meth:`.QaQc.write`, this option keeps the location within the directory of
   the config file but just changes the name, whereas to change the entire
   output directory path adjust the :attr:`.QaQc.out_dir` attribute directly.
   Also, the naming scheme of output files created will use user-defined names
   for all input variables.

Visualize resampled and corrected data
--------------------------------------

Similar to the :meth:`.Data.plot`, the :meth:`.QaQc.plot` method creates a
series of default time series and scatter plots of input and in this
case calculated variables. The temporal frequency of plots from
:meth:`.QaQc.plot` will always be daily and monthly and additional plots are
created for validation of energy balance closure corrections, otherwise the
same options such as number of subplot columns, super title, subplot
dimensions, output type, output file path, etc. are available. Similar to
:meth:`.QaQc.write` and :attr:`.QaQc.monthly_df`, if the data has not yet been
corrected the ``plot`` method will correct it using the default parameters
before creating the plots. 

Here is an example of the default daily and monthly time series plots produced after running the Energy Balance Ratio closure correction:

    >>> q = QaQc(d)
    >>> q.plot(output_type='notebook', plot_width=700)

.. raw:: html
    :file: _static/tutorial/US-Tw3_plots.html



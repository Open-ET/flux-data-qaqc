Change Log
==========

Version 0.2.2
-------------

Added method to write input data to a CSV file following the same standardized formatting and unit conversions that are implemented in ``qaqc.write``. This method is ``data.write``. This was done so that a user can easily rewrite the initially read data at its native time frequency that is often half-hourly or hourly as produced by eddy covariance processing software such as EddyPro. This is useful for creating input for sub-daily time series analyses that may be done in conjunction with ``flux-data-qaqc``.

Bug fixes related to internal automatic calculations for vapor pressure, vapor pressure deficit, saturation vapor pressure, and dew point temperature and data assignment not persisting until two calls to ``data.df``. 

Fix multiple deprecation warnings caused by ``Pandas`` version 2, tested with version 2.2.2. 

Version 0.2.1
-------------

Added option to specify whether the threshold value used in the ``data.apply_qaqc_flags`` is to filter values that are less than or greater than the value given. Previously the function only removed data that were less than the threshold value. 

Clean up dependencies in requirements.txt to match version 0.2.0 and keeping only tested and required packages. 

Add new config file for ReadTheDocs, update sphinx to version 7.2.6 and fix deprecation errors.

Version 0.2.0
-------------

Update dependencies to major new versions. ``Pandas`` was upgraded to version 1.0 and ``Bokeh`` to version 3.0. To accomodate these major dependency releases several functions' kwargs and syntax were modified to avoid deprecation errors and warnings. Because these version changes would not be backward compatible with previous versions for ``flux-data-qaqc``, it's version was also bumped up a minor version from 0.1 to 0.2. Tests were also run using Python version 3.10, previous Python versions may be compatible but were not tested. 

Version 0.1.6
-------------

Add automated tests using GitHb Actions, `see here <https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml>`__ and added in description of how to run tests locally on docs.

Remove ``xlrd`` reader as a dependency due to outdated reading ability as a ``Pandas`` excel reader.

Other minor bug fixes related to ``Plot`` class.

Ass JOSS paper and publish software on Zenodo.

Version 0.1.5
-------------

Add configuration writing function :func:`.util.write_configs` to ``util`` module to facilitate batch processing os similar formatted input files via a station metadata file and data dictionary. 

Update check on energy balance ratio closure correction to also check if the inverse of the energy balance ratio is greater than 0.5, in other words :math:`\frac{1}{EBR} > |0.5|` to avoid closure correction factors that are too small. This check occurs both after step 3 and 6 of the energy balance closure correction routine. 

Version 0.1.4
-------------

Relax default allowance for missing days threshold from 90 (~ 3 days) to 80 % (~ 6 days) in the monthly resample algorithm. In other words if a month has more than 80 % missing daily values, its monthly aggregate will not be resampled, it will be replaced with a null value. The threshold is a keyword argument to the ``util.monthly_resample`` function, but the default is used in any automatic resampling of variables. As a reminder, the number of missing days per month which is tabulated for some variables can be used to fine tune this filter. This change was implemented in version 0.1.4.post1. 

Add daily ASCE standardized reference ET calculation option from the :meth:`.QaQc.daily_ASCE_refET` method. Also added automatic estimation of daily maximum and minimum air temperature from input (e.g. hourly) data and added the input variables to the list of variables that are linearly interpolated before taking daily aggregates in the :obj:`.QaQc` constructor. In other words, the inputs to the daily ASCE reference ET formulation: ea, tmin, tmax, rs, wind speed, are interpolated over daytime and nighttime hourly gaps (2 and 4 default) before taking daily means, mins, maxs, and subsequently used in the daily ASCE calculations. 

Changed default keyword argument ``reference`` to "short" of the :meth:`.Data.hourly_ASCE_refET` method.

Add automatic calculations for high frequency (e.g. hourly or half hourly) data including dew temperature and relative humidity from ea and es if available. The calculations occur when first loading input data, i.e. when :obj:`.Data.df` attribute is accessed. Saturation vapor pressure (es) if calculated at hourly/daily frequency is now saved and added to :obj:`.Data.df` and :obj:`.QaQc.df` properties. 

Require Pandas >= 1.0, changes are not backwards compatible due to internal pandas argument deprecations particularly in the ``pandas.grouper`` object. 

Require Bokeh >= 2.0, changes are not backwards compatible due to legend keyword argument name changes in Bokeh 2.

Minor changes to remove package deprecation warnings from ``Pandas`` and ``Bokeh`` related to their respective large changes. 

Add package dependency ``openpyxl`` package as a fallback for reading in headers of Excel files when ``xlrd`` is unmaintained and failing with previously working tools for reading metadata on Excel files. 

Add a requirements.txt file with package.


Version 0.1.3
-------------

Add option to use gridMET grass reference ET (ETo) and EToF for gap filling daily ET. The default behavior still uses alfalfa reference ET, to use ETo assign the ``refET="ETo"`` keyword argument to :meth:`.QaQc.correct_data` or directly to :meth:`.QaQc._ET_gap_fill`. The ET and ET reference fraction plot labels are updated to show the correct reference ET variable used.

Improve scaling of scatter plots to give equal x and y axis lengths, change return of :meth:`.Plot.scatter_plot` to return tuple of (xmin, xmax, ymin, ymax) for use in plotting one to one lines or limiting axes lengths. 

Version 0.1.2
-------------

Change default functionality of the :meth:`.QaQc.write` method to use the internal variable names (as opposed to the input names) of ``flux-data-qaqc`` in the header files of the output daily and monthly time series CSV files. For example, the column for net radiation is always named and saved as "Rn". This can be reversed to the previous behavior of using the user's input names by setting the new ``use_input_names`` keyword argument to :meth:`.QaQc.write` to ``True``. 

Change the :meth:`.Plot.scatter_plot` underlying call to the ``bokeh`` modules scatter plot as opposed to the set circle glyph plot. This allows the user to change the symbol from circle to others by passing a valid value to the scatter_plot's ``marker`` keyword argument, e.g. ``marker='cross'``.

Version 0.1.1
-------------

Add least squares linear regression method for single or multivariate input; specifically the ``QaQc.lin_regress()`` method. It can be used to correct energy balance components or for any arbitrary time series data loaded in a ``QaQc`` instance. It produces and returns a readable table with regression results (fitted coefficients, root-mean-square-error, etc.) which can be accessed from ``QaQc.lin_regress_results`` after calling the method. The default regression if used to correct energy balance components assumes net radiation is accurate (as the dependent variable):

:math:`Rn = c_0 + c_1 G + c_2 LE + c_3 H`

where :math:`c_0 = 0`.

This regression utilizes the scikit-learn Python module and therefore it was added to the environment and setup files as a dependency.

Version 0.1.0
-------------

Add hourly ASCE standardized reference ET calculation to the ``Data`` class as :meth:`.Data.hourly_ASCE_refET` with options for short and tall (grass and alfalfa) reference ET calculations. If the input data is hourly or higher frequency the input data for the reference ET calculation will automatically be resampled to hourly data. If the input data is hourly then the resulting reference ET time series will be merged with the :attr:`.Data.df` attribute otherwise if the input data is at a temporal frequency > hourly, then the reference ET time series will be return by the :meth:`.Data.hourly_ASCE_refET` method. 

Add methods and options to linearly interpolate energy balance variables based on length of gaps during daytime (:math:`Rn > 0`) and night (:math:`Rn < 0`). These methods are run automatically by the ``QaQc`` constructor if temporal frequency of input is detected as less than daily. New keyword arguments to ``QaQc`` are ``max_interp_hours`` and ``max_interp_hours_night`` respectively.

Other notable changes:

* first release on GitHub
* creation of this file/page (the Change Log)
* add optional return options to plot methods of ``Data`` and ``QaQc`` objects for custimization of default plots or to show/use a subset of them

Version 0.0.9
-------------

Major improvements and notabable changes include:

* add package to PyPI
* change allowable gap percentage for monthly time series to 10 % from 70 %
* add reading of wind direction data, BSD3 license, add package data
* fix bugs related to filtering of subday gaps
* improve plots and other error handling, add feature to hide lines in line plots

Version 0.0.5
-------------

Major improvements and notabable changes include:

* first documentation on `ReadTheDocs <https://flux-data-qaqc.readthedocs.io/en/latest/>`__
* add multiple pages in docs such as installation, config options, basic tutorials, full API reference, etc. 
* improve and streamline config file options
* add vapor pressure and vapor pressure deficit calculations for hourly or lower frequency data in the ``Data.df`` property (upon initial loading of time series into memory
* add automatic unit conversions and checks on select input variables using the ``Convert`` class in the ``util`` module
* add new plots in default plots from ``QaQc`` class, e.g. filtered and raw ETrF
* many rounds of improvements to plots, e.g. hover tooltips, linked axes, style, options for columns, etc. 
* modify Energy Balance Ratio to filter out extreme values of filtered Energy Balance Ratio correction factors
* improve temporal resampling with options to drop days with certain fraction of sub-daily gaps
* track number of gap days in monthly time series of corrected ET 
* add examples of ET gap-filling to docs and change most example data to use Twitchel Island alfalfa site data from AmeriFlux
* add plotting of input data using ``plot`` method of ``Data`` instance which allows for viewing of input data at its initial temporal frequency


Version 0.0.1
-------------

First working version, many changes, milestones included: 

* basic templates and working versions of the ``Data``, ``QaQc``, and ``Plot`` classes 
* versions and improvements to daily and monthly resampling 
* Bowen and Energy Balance Ratio correction routines 
* example Jupyter notebooks including with FLUXNET and USGS data 
* calculation of potential clear sky radiation 
* changing variable naming system to use internal and user names 
* ability to read in multiple soil heat flux and soil moisture measurements and calculate weighted averages 
* make package installable and Conda environment
* add input data filtering using quality control flags (numeric threshold and flags)
* reading of input variables' units
* added the ``util`` submodule with methods for resammpling time series
* ability to take non-weighted averages for any acceptable input variable
* add config file options like date parsing
* removed filtering and smoothing options from Bowen Ratio method and other modifications to it
* add methods for downloading gridMET variables based on location in CONUS
* add routine for gap filling ET based on gridMET ETrF that is smoothed and filtered
* improved ``Plot`` class to contain modular plot methods (line and scatter) for use with arbitrary data
* changed internal variable naming, e.g. etr to ETr
* methods to estimate ET from LE that consider the latent heat of vaporization is affected by air temp.
* other updates to improve code structure and optimization of calculations

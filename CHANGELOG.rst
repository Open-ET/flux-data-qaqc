Change Log
==========

Version 0.4.1
-------------

The daily ``ebr_corr`` variable of the Energy Balance Ratio correction
(``meth='ebr'``, default) is now the energy balance ratio after correction,
(H_corr + LE_corr) / (Rn - G), as it already was for the other correction
methods and for monthly data. Previously it was the energy balance ratio
estimated by the correction (filtered, smoothed, and gap filled, the inverse
of ``ebc_cf``), which is now saved as the new variable ``ebr_est``. The
energy balance ratio plots of :meth:`.QaQc.plot` show ``ebr_est`` in
addition to ``ebr`` and ``ebr_corr``. Infinite daily ``ebr`` values (Rn - G
equal to zero) are now set to null with the ``ebr`` method as with the other
methods. Corrected LE, H, and ET are unchanged.

Raise an error when ``mj/m2`` units are given for hourly or shorter input
data. The conversion to ``w/m2`` assumes daily totals (MJ m⁻² per day), so
sub-daily values were converted incorrectly without a warning, e.g. 48 times
too low for half-hourly data. Daily input in ``mj/m2`` is converted as
before.

The new documentation page "Variables and Units" lists all calculated
variables and the accepted input units for each variable.

Version 0.4.0
-------------

This release changes computed values, mainly corrected LE, H, and ET from the
default Energy Balance Ratio correction. On 293 AmeriFlux sites with all four
energy balance components, pooled corrected ET decreased 1.3% compared to
version 0.3.3 (median site -1.0%, 5th to 95th percentile of sites -3.3% to
+0.2%) and 1,372 more days were corrected instead of gap filled.

The Energy Balance Ratio correction (``meth='ebr'``) now follows the daily
energy balance closure correction of FLUXNET2015 and ONEFlux (`Pastorello et
al., 2020 <https://doi.org/10.1038/s41597-020-0534-3>`__) more closely.
Outliers are removed with 1.5 times the interquartile range of the correction
factor, :math:`EBC_{CF} = (Rn - G) / (H + LE)`, instead of the Energy Balance
Ratio, and the moving window statistics use the correction factor. Method 2
now uses the mean instead of the median of the +/- 5 day window, as in
ONEFlux. Remaining gaps are filled with the mean correction factor within
+/- 5 days of the same day in the previous and next years (method 3, new)
before the 5 day climatology of all years is used (method 4), which was
retained as a backup and is a deviation from ONEFlux. A correction factor from
any method that is outside of 0.5 to 2 is replaced by the next method (these
limits were retained and are also a deviation from ONEFlux). Fix the moving
windows for the first 7 days of the record, which were empty so these days
could only get the 5 day climatology, and align the 5 day climatology with
the day of year. The new daily variable ``ebc_cf_method`` reports the method
used for each day: the 15 day median (1), the 11 day mean (2), the same days
in the previous and next years (3), or the 5 day climatology (4).

Calculate daily sums, e.g. precipitation, photosynthetic energy storage, and
sub-daily gap counts, for each variable on its own. Previously a missing
record in any summed variable removed that time step from all of them, and
days without data were summed to zero instead of null, causing potential errors in those variables.

Calculate daily ``t_min`` and ``t_max`` from sub-daily ``t_avg`` with or
without linear interpolation of short sub-daily gaps (``max_interp_hours``
and ``max_interp_hours_night``), previously they were only calculated with
interpolation. With ``drop_gaps=True`` (default) they are now set to null on
the same days as ``t_avg``, days with fewer than ``daily_frac`` (default 1.0,
e.g. 48 of 48 half-hourly records) of the sub-daily air temperature records
after interpolation. Previously they were not filtered.

Keep monthly values when exactly 80% of the days in the month exist, as
documented. Previously these months were set to null.

Fix the conversion of air temperature from Fahrenheit to Celsius.

The Bowen Ratio correction (``meth='br'``) no longer produces infinite
values, days where the Bowen Ratio is undefined or the corrected LE is
outside of -100 to 850 w/m2 (the same limits used by the Energy Balance Ratio
correction) are not corrected.

Fix the sign of G in ``energy_corr`` (Rn - G) from :meth:`.QaQc.lin_regress`.

Fix the check of the input temporal frequency for **daily** data, it used the
number of columns instead of the number of records per day. Vapor pressure or
vapor pressure deficit are now calculated for daily input data, and hourly
ASCE reference ET cannot be calculated from daily data.

:meth:`.QaQc.download_gridMET` prints an error and does not download data
when the station is outside of the gridMET domain (contiguous US).

:meth:`.QaQc.lin_regress` and ``QaQc.correct_data(meth='lin_regress')``
accept ``y='energy'`` for available energy (Rn - G) and a single independent
variable, which previously raised an error. The warning for other variables
now says that the results may be hard to interpret in terms of energy balance
closure.

Regenerate all figures in the online documentation with the current version
and update the documentation for these changes, including a new section on
the linear regression option.

Version 0.3.3
-------------

Same code as intended for version 0.3.2. The version 0.3.2 files on PyPI were
built by mistake from a working copy that included unreleased changes to the
energy balance ratio correction, that release has been yanked.

Version 0.3.2
-------------

Fix the check for missing units added in version 0.3.1, which raised an error
for config files where a multiple soil heat flux entry uses the same column
as the ground heat flux, e.g. ``ground_flux_col = G`` and ``g_1 = G``, even
though its units were given. The check now reads units from the config file.

Version 0.3.1
-------------

Bug fix and compatibility release. Computed values are unchanged for the
example datasets with the previously required ``refet`` 0.3.10.

Call :mod:`refet` with keyword arguments in :meth:`.Data.hourly_ASCE_refET`
and :meth:`.QaQc.daily_ASCE_refET`. ``refet`` 0.4 changed the order of its
arguments, so environments with ``refet`` 0.4 (e.g. installed from the
previous ``requirements.txt``) produced incorrect ASCE reference ET. Also
support the renamed radiation functions in ``refet`` 0.5.

Raise a :obj:`ValueError` when creating a :obj:`.Data` object if units are
missing from the config file for variables used in calculations or unit
conversions (:attr:`.Convert.required_units`) or for multiple soil heat flux
variables, units are never assumed. Missing units for variables only used in
plot labels, e.g. relative humidity or soil moisture, still give a warning.
Previously missing units caused an unrelated error in :obj:`.QaQc`.

Fix :meth:`.QaQc.download_gridMET` with pandas 3 (`#34
<https://github.com/Open-ET/flux-data-qaqc/issues/34>`__), when passing a
single variable name, e.g. ``'ETr'``, and when downloading again. If the
gridMET server cannot be reached an error is printed and ET gap filling is
skipped instead of raising an exception.

Save the gridMET file path relative to the config file and update only that
line of the config file with the new :func:`.util.set_config_option`.
Previously the full config file was rewritten with an absolute path, which
removed all comments. Absolute paths in existing config files still work.

Fix :meth:`.QaQc.lin_regress` and ``QaQc.correct_data(meth='lin_regress')``
with newer versions of scikit-learn. Accept ``anemometer_height`` from the
config file for ASCE reference ET (it was read as text), and raise a
:obj:`ValueError` for invalid ``reference`` or ``refET`` options.

Read Excel file headers with pandas instead of the openpyxl shared strings
table. Fix config entries for QC flags of multiple soil variables, e.g.
``g_1_qc``, which raised an error and were averaged into G. Accept ``ppm``
for CO2 mole fraction and ``nondimensional`` for the stability parameter.

Remove ``QaQc.from_dataframe`` which was not functional.

Require pandas >= 2.2 (needed for monthly resampling) and allow ``refet`` >=
0.3.10. ``requirements.txt`` and ``environment.yml`` now match ``setup.py``.
Automated tests run on Python 3.9 with the oldest supported dependency
versions and on Python 3.11 and 3.13 with the newest, and new tests cover
these fixes, unit handling, gridMET download (without network access), and
writing and plotting daily and monthly output.

Version 0.3.0
-------------

Fix sub-daily linear interpolation so that the configured maximum gap length
applies to the entire consecutive gap. Previously,
:meth:`pandas.DataFrame.interpolate` with limits applied in both directions
could completely fill gaps up to twice the configured length.

Interpolate the original continuous time series rather than separating daytime
and nighttime records before interpolation. This allows valid observations on
both sides of a short gap to be used when the gap crosses sunrise or sunset.

Leave gaps longer than the applicable limit entirely missing rather than
partially filling them. Missing values at the beginning or end of a time series
also remain unchanged because linear interpolation requires valid observations
on both sides.

Following the ONEFlux/FLUXNET2015 nighttime convention described by
`Pastorello et al. (2020) <https://doi.org/10.1038/s41597-020-0534-3>`__, calculate sub-daily potential
incoming shortwave radiation, `sw_pot`, to identify solar night. The
calculation uses the ASCE-EWRI hourly extraterrestrial-radiation equations,
generalized to the detected sub-daily period length.

Determine the station standard UTC offset from its latitude and longitude so
that timestamps recorded in local standard time can be converted to the UTC
midpoint time required by the sub-daily solar-radiation calculation. This adds
`timezonefinder` and `tzdata` as package dependencies.

Use the longer nighttime interpolation limit only when `sw_pot` is zero
throughout the gap and at both bounding observations. As an additional check,
all available `Rn` values over the same period must be negative. Missing
`Rn` values do not by themselves prevent use of the nighttime limit. Gaps
containing daylight, crossing sunrise or sunset, or lacking sufficient solar
information use the shorter daytime limit. If `max_interp_hours_night` is
`None`, the daytime interpolation limit is used.

Fix monthly resampling of boolean variables such as `ET_gap` so that missing
days are not estimated from the monthly mean before summation.

Raise the minimum supported Python version from 3.7 to 3.9.

Add tests for gap interpolation, day-night transitions, known gaps in
the example flux data, standard UTC offsets, agreement with RefET (ASCE) 
hourly extraterrestrial radiation, daily radiation totals, and expected 
daytime and nighttime behavior.


Version 0.2.3
-------------

Add optional calculation of photosynthetic energy storage from sub-daily gross
primary productivity via :meth:`.Data.calc_pes`.

Add support for new optional carbon, radiation, and turbulence-related input
variables including gross primary productivity, ecosystem respiration, net
ecosystem exchange, turbulent and storage CO2 fluxes, incoming photosynthetic
photon flux density, CO2 mole fraction, friction velocity, Monin-Obukhov
length, stability parameter, lateral velocity fluctuation standard deviation,
and boundary layer height.

Add temporal aggregation support for PES variables in :obj:`.QaQc`, where
``pes`` is summed and ``pes_flux`` is averaged during daily and monthly
resampling.

Update documentation and tutorial materials to describe the PES calculation,
equations, variable naming, and example plotting workflow.

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

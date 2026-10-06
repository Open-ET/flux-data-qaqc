Closure Methodologies
=====================

``flux-data-qaqc`` currently provides two routines which ultimately adjust
turbulent fluxes in order to improve energy balance closure of eddy covariance
tower data, the Energy Balance Ratio and the Bowen Ratio method, and a linear
regression option that is mainly a diagnostic tool (described at the end of
this page). 

Closure methods are assigned as keyword arguments to the 
:meth:`.QaQc.correct_data` method, and for a list of provided 
closure options see :attr:`.QaQc.corr_methods`.
For example, if you would like to run the Bowen Ratio correction routine
assuming you have succesfully created a :obj:`.QaQc` object,

.. code-block:: python

    # q is a QaQc instance
    q.correct_data(meth='br')

The other keyword argument for :meth:`.QaQc.correct_data` allows for gap filling corrected evapotranspiration (:math:`ET`) which is calculated from corrected latent energy (:math:`LE`). By default the gap filling option is set to True, more details on this below in :ref:`Step 9, optionally gap fill corrected ET using gridMET reference ET and reference ET fraction`.

.. Tip:: 
   All interactive visualizations in this page were created using
   :meth:`.Plot.line_plot`, :meth:`.Plot.add_lines`, and
   :meth:`.Plot.scatter_plot` which automatically handle issues with utilizing
   the mouse hover tooltips and other :obj:`bokeh.plotting.figure.Figure`
   features.
    
Data description
----------------

The data for this example comes from the "Twitchell Alfalfa" AmeriFlux eddy
covariance flux tower site in California. The site is located in alfalfa fields and exhibits a mild Mediterranean climate with dry and hot summers, for more information on this site or to download data click `here <https://ameriflux.lbl.gov/sites/siteinfo/US-Tw3>`__. 


Energy Balance Ratio method
---------------------------

The Energy Balance Ratio method (default) follows the daily energy balance
closure correction of the `FLUXNET2015 dataset
<https://fluxnet.org/data/fluxnet2015-dataset/data-processing/>`__
and ONEFlux processing pipeline described by `Pastorello et al. (2020)
<https://doi.org/10.1038/s41597-020-0534-3>`__, with additional limits on the
correction factor and corrected fluxes used by ``flux-data-qaqc``. The method
calculates a daily energy balance closure correction factor, the inverse of
the Energy Balance Ratio, filters out extreme values, smooths and gap fills
it using moving windows, and multiplies the initial latent energy
(:math:`LE`) and sensible heat (:math:`H`) flux time series by it.

.. note::
   Changed in version 0.4.0: outliers are filtered and moving window
   statistics are calculated on the correction factor instead of the Energy
   Balance Ratio, method 2 uses the mean, and remaining gaps are first filled
   from the previous and next years as done by FLUXNET before the all-year
   climatology is used.

All steps, abbreviated
^^^^^^^^^^^^^^^^^^^^^^

Below is a step-by-step description of the Energy Balance Ratio
correction routine used by ``flux-data-qaqc``. More details and visual
demonstration of steps are shown below.

**Step 0 (optional):** optionally filter out poor quality data first if quality
control (QC) values or flags are provided with the dataset or other means. For
example, FLUXNET data includes QC values for :math:`H` and :math:`LE`,
e.g. H_F_MDS_QC and LE_F_MDS_QC are QC values for gap filled :math:`H` and
:math:`LE`. This allows for manual pre-QaQc of data.

**Step 1:** calculate the daily Energy Balance Ratio (EBR =
:math:`\frac{H + LE}{Rn – G}`) and its inverse, the energy balance closure
correction factor (:math:`{EBC_{CF}} = \frac{Rn - G}{H + LE}`), from daily
mean values.

**Step 2:** remove :math:`EBC_{CF}` values that are outside 1.5 times the
interquartile range.

**Step 3 (method 1):** for each day, a sliding window of +/- 7 days (15 days)
is used to select up to 15 :math:`EBC_{CF}` values, if at least 5 values
exist use their median.

**Step 4 (method 2):** if fewer than 5 values exist in the 15 day window, use
the mean :math:`EBC_{CF}` of a +/- 5 day (11 day) sliding window. In steps 3
and 4, a correction factor that is :math:`\le 0.5` or :math:`\ge 2` is left
as a gap for the next steps.

**Step 5 (method 3):** fill remaining gaps with the mean :math:`EBC_{CF}`
within +/- 5 days of the same day in the previous and next years.

**Step 6 (method 4):** fill gaps that still remain, e.g. in records with less
than two years of data, with the 5 day climatology: the +/- 5 day moving mean
of the day of year mean :math:`EBC_{CF}` from steps 3 and 4 over all years.
This step is not part of the FLUXNET method.

**Step 7:** correct :math:`LE` and :math:`H` by multiplying by
:math:`EBC_{CF}` if it is between 0.5 and 2. If the corrected :math:`LE` is
greater than 850 or less than -100 :math:`w/m^2` no correction is made for
that day. The method used for each day (1-4) is saved as ebc_cf_method.

**Step 8:** calculate corrected :math:`ET` from corrected :math:`LE`
using average air temperature to adjust the latent heat of vaporization.

**Step 9 (optional):** if desired, fill remaining gaps in the corrected
:math:`ET` time series with :math:`ET` that is calculated by gridMET
reference :math:`ET` (:math:`ETr` or :math:`ETo`) multiplied by the filtered and smoothed
fraction of reference ET (:math:`ETrF` or :math:`EToF`).


Step 0, manual cleaning of poor quality data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below we can see that the daily time series of net radiation (:math:`Rn`) has
some periods of poor quality data. This is a common issue due, e.g. to
instrumentation problems, that cannot always be avoided. In this case the
sensor did not record values at night (or they were not provided with the data)
when :math:`Rn` values are lower for several days (e.g. around 8/26/2014) which
resulted in overestimates of daily mean :math:`Rn` during these periods. Although these days can automatically be filtered out by the :obj:`.QaQc` class, the example below shows a way of manually filtering them because in other cases outliers in the daily data may not be caused by resampling of sub-daily data with systematic measurement gaps. The main point is that manual inspection and potentially pre-filtering of poor quality data before proceeding with energy balance closure corrections is often necessary.

.. raw:: html
    :file: _static/closure_algorithms/step0_NotFiltered.html

There are several ways to conduct manual pre-filtering of poor quality meterological time series data, to filter data based on input quality flags or numeric quality values see :ref:`Quality-based data filtering`. 

``flux-data-qaqc`` also allows for filtering of poor quality data on the fly as shown in this example. In other words, we simply filter out the periods we think have bad data for :math:`Rn` within Python before running the closure correction. After manually determing the date periods with poor quality :math:`Rn`, here is how they were filtered oiut before running the correction:

    >>> import pandas as pd 
    >>> import numpy as np
    >>> from fluxdataqaqc import Data, QaQc
    >>> d = Data('Path/to/config.ini')
    >>> # days with sub daily gaps can be filtered out automatically here, 
    >>> # see "Tip" below the following plot 
    >>> q = QaQc(d, drop_gaps=False) 
    >>> # rename dataframe columns for ease of variable access, adjust
    >>> df = q.df.rename(columns=q.inv_map)

Here were the dates chosen and one way to filter them,

    >>> # make a QC flag column for Rn
    >>> df['Rn_qc'] = 'good'
    >>> df.loc[pd.date_range('2/10/2014','2/10/2014'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('8/25/2014','9/18/2014'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('10/21/2015','10/26/2015'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('10/28/2015','11/1/2015'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('7/23/2016','7/23/2016'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('9/22/2016','9/22/2016'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('3/3/2017','3/3/2017'), 'Rn_qc'] = 'bad'
    >>> # filter (make null) based on our QC flag column for Rn
    >>> df.loc[df.Rn_qc == 'bad', 'Rn'] = np.nan
    >>> # reassign to use pre-filtered data for corrections
    >>> q.df = df

The resulting energy balance component plot with :math:`Rn` filtered:

.. raw:: html
    :file: _static/closure_algorithms/step0_Filtered.html

.. tip::
   In this case, the issues with :math:`Rn` were caused by resampling 30 minute
   data with systematic night-time gaps. These sort of issues can be
   automatically handled when creating a :obj:`.QaQc` object; the keyword
   arguments ``drop_gaps`` and ``daily_frac`` to the :obj:`.QaQc` class are
   used to automatically filter out days with measurement gaps of varying size,
   i.e.,
   
   >>> d = Data('path/to/config.ini')
   >>> q = QaQc(d, drop_gaps=True, daily_frac=0.8)
   >>> q.correct_data()

   This would produce very similar energy balance closure results as the manual
   filter above. Another more fine-grained option would have been to flag the
   days with gaps in the sub-daily input time series that you would like to
   filter by :meth:`.Data.apply_qc_flags`.


.. Note::
   The remaining step-by-step explanation in this page uses the pre-filtered
   input time series, however results of the energy balance closure correction
   without pre-filtering outliers of :math:`Rn` are also shown in plots for the
   final steps (8 and 9) for comparison. If you now ran:

   >>> q.df = df
   >>> q.correct_data()
   >>> q.plot(output_type='show')

   This will directly produce the same output of step 9 using the 
   pre-filtered data. 

Steps 1 and 2, filtering outliers of the correction factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate the daily EBR = :math:`\frac{H + LE}{Rn - G}` and correction
factor :math:`EBC_{CF} = \frac{Rn - G}{H + LE}` time series and filter out
correction factors that are outside 1.5 times the interquartile range, as
done by FLUXNET. Note, in ``flux-data-qaqc`` these are named “ebr” and
“ebc_cf”. Versions before 0.4.0 applied this filter to EBR, which removes a
different set of days because the inverse is not linear. The plot below shows
the correction factor before and after the filter.

.. raw:: html
    :file: _static/closure_algorithms/steps1_2_PreFiltered.html

Steps 3 and 4, moving window statistics of the correction factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Smooth and fill the filtered correction factor time series using moving
windows. Specifically, take the median :math:`EBC_{CF}` from a +/- 7 day
moving window if it has at least 5 values (method 1), otherwise take the
mean from a +/- 5 day moving window (method 2). Windows at the start and end
of the record only use the days that exist. If the resulting correction
factor is :math:`\le 0.5` or :math:`\ge 2` leave a gap for the day for
filling in later steps. The plot below shows the filtered correction factor
and the correction factor after steps 3 and 4.

.. raw:: html
    :file: _static/closure_algorithms/steps3_4_PreFiltered.html

Steps 5 and 6, fill remaining gaps from other years
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Days that still have no correction factor after steps 3 and 4 are mostly in
gaps of 11 or more days. These are filled with the mean of the filtered
:math:`EBC_{CF}` within +/- 5 days of the same day in the previous and next
years (method 3, as done by FLUXNET). If neither of those years has data,
e.g. in records shorter than two years, the gap is filled with the 5 day
climatology of the correction factor (method 4): the day of year mean of the
correction factors from steps 3 and 4 over all years in the record, smoothed
with a moving +/- 5 day (11 day) window. Correction factors from these steps
that are :math:`\le 0.5` or :math:`\ge 2` are not used. In versions before
0.4.0 the 5 day climatology was used for all remaining gaps. The plot below
shows the days filled by each step, at this site mostly before the tower was
installed in May 2013 and in a gap in 2016.

.. raw:: html
    :file: _static/closure_algorithms/steps5_6_PreFiltered.html

``flux-data-qaqc`` also keeps a record of the 5 day climatology as an
Energy Balance Ratio (inverse of the climatology of the correction factor,
shown below), it is named by ``flux-data-qaqc`` as ebr_5day_clim. The method
used to get each day's correction factor is saved as ebc_cf_method, 1-4 for
methods 1-4 above and null for days without a correction factor.

.. raw:: html
    :file: _static/closure_algorithms/5dayclim_PreFiltered.html

Steps 7 and 8 correct turbulent fluxes, EBR, and ET
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate corrected :math:`LE` and :math:`H` by multiplying by
:math:`EBC_{CF}`, the filtered and gap filled correction factor from
previous steps:

.. math:: LE_{corr} = LE \times EBC_{CF}

\ and

.. math:: H_{corr} = H \times EBC_{CF}.

The daily corrected EBR (ebr_corr) saved by ``flux-data-qaqc`` is the
filtered and smoothed EBR, :math:`\frac{1}{EBC_{CF}}`. In monthly output it
is calculated from corrected :math:`LE` and :math:`H`,

.. math:: EBR_{corr} = \frac{H_{corr} + LE_{corr}}{Rn - G}.

Calculate ET from LE using average air temperature to adjust the latent
heat of vaporization following the method of Harrison, L.P. 1963,

.. math:: ET_{mm \cdot day^{-1}} = 86400_{sec \cdot day^{-1}} \times \frac{LE_{w \cdot m^{-2}}}{2501000_{w \cdot sec \cdot kg^{-1}} - (2361 \cdot T_{C})}, 

where evapotransipiration (:math:`ET`) in :math:`mm \cdot day^{-1}`,
:math:`LE` is latent energy flux in :math:`w \cdot m^{-2}`, and
:math:`T` is air temperature in degrees celcius. The same approach is
used to calculate corrected :math:`ET` (:math:`ET_{corr}`) using
:math:`LE_{corr}`.

The plot below shows the time series of the initial and corrected ET (:math:`ET` and :math:`ET_{corr}`).

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_PreFiltered.html

There were not significant gaps in the energy balance components for this dataset and therefore step 9 was not used, although it is still demonstrated with an artificial gap in the next step. 

The following plot shows the energy balance closure of the initial and corrected data after applying the steps above, including the manual pre-filtering of :math:`Rn`,

.. raw:: html
    :file: _static/closure_algorithms/EBC_scatter_PreFiltered.html

Notice the mean daily corrected energy balance ratio (slope of corrected) is 0.99 or near perfect closure. However, the same plot below shows the results if we skipped the manual pre-filtering of outlier :math:`Rn` values. In this case the resulting corrected mean closure is only 0.92:

.. raw:: html
    :file: _static/closure_algorithms/EBC_scatter_noPreFilter.html

.. Tip:: 
   These and other interactive visualizations of energy balance closure results 
   are provided by default via the :meth:`.QaQc.plot` method.

In ``flux-data-qaqc`` new variable names from these steps are: LE_corr, H_corr,
ebr, ebr_corr, ebc_cf, ebc_cf_method, ET, ET_corr, and ebr_5day_clim. The
energy balance closure correction factor is named ebc_cf as in the `FLUXNET
methodology <https://fluxnet.org/data/fluxnet2015-dataset/data-processing/>`__
and Pastorello et al. (2020).

Step 9, optionally gap fill corrected ET using gridMET reference ET and reference ET fraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is done by downloading :math:`ETr` or :math:`ETo` (default is :math:`ETr`)
for the overlapping gridMET cell (site must be in CONUS) and then calculating,

.. math:: ET_{fill} = ETrF \times ET_r,

\ where

.. math:: ETrF = \frac{ET_{corr}}{ET_r}

:math:`ET_{corr}` is the corrected ET produced by step 8 and :math:`ETrF` is
the fraction of reference ET. :math:`ETrF` if first filtered to remove
outliers outside of 1.5 times the interquartile range, it is then smoothed with
a 7 day moving average (minimum of 2 days must exist in window) and lastly it
is linearly interpolated to fill any remaining gaps. 

The same gap filling procedure can easily be done using gridMET grass reference ET (:math:`ETo`) as opposed to alfalfa reference ET (:math:`ETr`).

.. Tip:: 
   The filtered and raw versions of :math:`ETrF`/:math:`EToF`, gridMET
   :math:`ETr`, gridMET :math:`ETo`, gap days, and monthly total number of gap
   filled days are tracked for post-processing and visualized by the
   :meth:`.QaQc.plot` and :meth:`.QaQc.write` methods.

Since the data used in this example does not have gaps, for illustration we have created the following large gap in the measured energy balance components from May through August, 2014:

.. raw:: html
    :file: _static/closure_algorithms/step0_Filtered_withGap.html

The resulting time series of :math:`ET_{corr}` using the optional gap filling method described is shown below.  

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_with_gapFill.html

Note, the gap filled values of :math:`ET` (green line) do not accurately catch the harvesting cycles of alfalfa however the :math:`ET_{corr}` values (blue line) do, this is because the gap filled values are based from gridMET reference ET which is not locally representative. If this is hard to see, try using the box zoom tool on the right of the plot to zoom in on the gap-filled period.

This ET gap-filling step is used by default when running ``flux-data-qaqc``
energy balance closure correction routines, to disable it set the
``et_gap_fill`` argument of :meth:`QaQc.correct_data` to False, e.g.

.. code-block:: python

    # q is a QaQc instance
    q.correct_data(meth='ebr', et_gap_fill=False)


In ``flux-data-qaqc`` new variable names from this step are: ETrF,
ETrF_filtered, gridMET_ETr, ET_gap, ET_fill, and ET_fill_val. The
difference between ET_fill and ET_fill_val is that the latter is masked
(null) on days that the fill value was not used to fill gaps in
:math:`ET_{corr}`. Also, ET_gap is a daily series of True and False
values indicating which days (from step 8) of :math:`ET_{corr}` were gaps
that were subsequently filled.

.. Note::
   When using the :math:`ETr`-based gap-filling option, any gap filled days
   will also be used to fill in gaps of :math:`LE_{corr}`, therefore the mean
   closure as found in the daily and monthly closure scatter plot outputs (from
   :meth:`.QaQc.plot`) will be updated to reflect the influence of the
   gap-filled days.

Bowen Ratio method
------------------

The Bowen Ratio energy balance closure correction method implemented
here follows the typical approach where the corrected latent energy
(:math:`LE`) and sensible heat (:math:`H`) fluxes are adjusted the
following way

.. math::  LE_{corr} = \frac{(Rn - G)}{(1 + \beta)}, 

\ and

.. math::  H_{corr} = LE_{corr} \times \beta 

where :math:`\beta` is the Bowen Ratio, the ratio of sensible heat flux
to latent energy flux,

.. math:: \beta = \frac{H}{LE}.

This routine forces energy balance closure for each day in the time
series. Days where the Bowen Ratio is undefined (:math:`LE` is zero or
:math:`\beta = -1`) or where the corrected :math:`LE` is :math:`\ge 850`
or :math:`\le -100` :math:`w/m^2`, the same limits used by the Energy
Balance Ratio method, are left without a correction (changed in version
0.4.0).

Here is the resulting :math:`ET_{corr}` time series using the pre-filtered (:math:`Rn`) energy balance time series and the Bowen Ratio method:

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_BR_PreFiltered.html

And here is the energy balance closure scatter plot which shows the forced closure of the method:

.. raw:: html
    :file: _static/closure_algorithms/EBC_BR_scatter_PreFiltered.html

New variables produced by ``flux-data-qaqc`` by this method include: br
(Bowen Ratio), ebr, ebr_corr, LE_corr, H_corr, ET, ET_corr, energy, flux, and
flux_corr.

Linear regression method
------------------------

The linear regression option (``meth='lin_regress'``) is mainly a diagnostic
tool, it was not designed to correct fluxes. It fits a least squares linear
regression between energy balance components of the daily data over the full
record. The dependent variable (``y``), the independent variables (``x``),
and whether an intercept is fit (``fit_intercept``) are chosen by the user.
By default :math:`Rn` is regressed on :math:`G`, :math:`LE`, and :math:`H`
without an intercept,

.. math:: Rn = c_1 G + c_2 LE + c_3 H,

and the coefficients show how much each component would need to be scaled,
on average, to close the energy balance if :math:`Rn` is assumed to be
correct. `Volk et al. (2023)
<https://doi.org/10.1016/j.agrformet.2023.109307>`__ used it this way to show
how :math:`LE` and :math:`H` relate to the closure deficit. For the
pre-filtered data used on this page,

    >>> q.lin_regress(y='Rn', x=['G', 'LE', 'H'])
    >>> q.lin_regress_results.T
        SITE_ID            US-Tw3
        Y (dependent var.)     Rn
        c0 (intercept)        0.0
        c1 (coef on G)      0.651
        c2 (coef on LE)     1.204
        c3 (coef on H)      0.986
        RMSE (w/m2)         13.26
        r2 (coef. det.)      0.94
        n (sample count)     1561

which suggests that most of the closure deficit at this site is in
:math:`LE` (about 20%) while :math:`H` is close to balanced.

Available energy (:math:`Rn - G`) can be used as the dependent variable with
``y='energy'``,

    >>> q.lin_regress(y='energy', x=['LE', 'H'])
    >>> q.lin_regress_results.T
        SITE_ID             US-Tw3
        Y (dependent var.)  energy
        c0 (intercept)         0.0
        c1 (coef on LE)      1.187
        c2 (coef on H)       0.967
        RMSE (w/m2)          13.43
        r2 (coef. det.)       0.93
        n (sample count)      1561

With ``fit_intercept=True`` the intercept (c0) is also estimated, e.g. 2.1
:math:`w/m^2` for the default regression at this site, a constant offset that
the coefficients alone cannot account for.

Running ``q.correct_data(meth='lin_regress')`` multiplies each independent
variable by its coefficient, :math:`G`, :math:`LE`, and :math:`H` by
default, or e.g. only :math:`LE` and :math:`H` with
``q.correct_data(meth='lin_regress', y='energy', x=['LE', 'H'])``. The
dependent variable is not corrected and the intercept is not applied. Unlike the Energy Balance Ratio and Bowen Ratio methods the
coefficients are the same for every day of the record, so the correction does
not follow daily or seasonal changes in closure. The plot below compares them
to the daily correction factor of the Energy Balance Ratio method, which
ranges from 0.82 to 1.55 at this site.

.. raw:: html
    :file: _static/closure_algorithms/lin_regress_coefs.html

See :meth:`.QaQc.lin_regress` for more options. New variables produced by
this method include: G_corr, LE_corr, H_corr, ET_corr, energy_corr,
flux_corr, and ebr_corr, and the regression results are saved in
:attr:`.QaQc.lin_regress_results`.



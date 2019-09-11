Closure Methodologies
=====================

``flux-data-qaqc`` currently provides two routines which ultimately adjust
turbulent fluxes in order to improve energy balance closure of eddy covariance
tower data, the Energy Balance Ratio and the Bowen Ratio method. 

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

The Energy Balance Ratio method (default) is modified from the `FLUXNET
methodology <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
(step 3 daily heat processing).
The method involves filtering out of extreme values of the daily Energy
Balance Ratio time series, smoothing, and gap filling. Then the inverse of
the filtered and smoothed time series is used as a series of
correction factors for the initial time series of latent energy
(:math:`LE`) and sensible heat (:math:`H`) flux time series.

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

**Step 1:** calculate the Energy Balance Ratio (EBR =
:math:`\frac{H + LE}{Rn – G}`) daily time series from raw data.

**Step 2:** filter EBR values that are outside 1.5 times the
interquartile range.

**Step 3:** for each day in the daily time series of filtered EBR, a
sliding window of +/- 7 days (15 days) is used to select up to 15
values.

**Step 4:** for each day take a percentile (default 50) of the 15 EBR
values. Check if the inverse of the EBR value is :math:`> |2|` or if the
the inverse of the ratio multiplied by the measured :math:`LE` would
result in a flux greater than 850 or less than -100 :math:`w/m^2`, if so
leave a gap for filling later.

**Step 5:** if less than +/- 5 days exist in the sliding 15 day window,
use the mean EBR for all days in a +/- 5 day (11 day) sliding window.
Apply same criteria for an extreme EBR value as in step 4.

**Step 6:** if no EBR data exist in the +/- 5 sliding window to average,
fill remaining gaps of EBR with the mean from a +/- 5 day sliding window
over the day of year mean for all years on record, i.e. 5 day
climatology. Calculate the 5 day climatology from the filtered and
smoothed EBR as produced from step 5. Apply same criteria for an extreme
EBR value as in steps 4 and 5.

**Step 7:** use the filtered EBR time series from previous steps to
correct :math:`LE` and :math:`H` by multiplying by the energy balance
closure correction factor :math:`{EBC_{CF}} = \frac{1}{EBR}`, where EBR
has been filtered by the previous steps. Use the corrected :math:`LE`
and :math:`H` to calculate the corrected EBR.

**Step 8:** calculate corrected :math:`ET` from corrected :math:`LE`
using average air temperature to adjust the latent heat of vaporization.

**Step 9 (optional):** if desired, fill remaining gaps in the corrected
:math:`ET` time series with :math:`ET` that is calculated by gridMET
reference :math:`ET` (:math:`ETr`) multiplied by the filtered and smoothed
fraction of reference ET (:math:`ET_{rF}`).


Step 0, manual cleaning of poor quality data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below we can see that the daily time series of net radiation (:math:`Rn`) has
some periods of poor quality data. This is a common issue due, e.g. to
instrumentation problems, that cannot always be avoided. In this case the
sensor did not record values at night (or they were not provided with the data)
when :math:`Rn` values are lower for several days (e.g. around 8/26/2014) which
resulted in overestimates of daily mean :math:`Rn` during these periods.
Therefore, manual inspection and pre-filtering of poor quality data (or in this
case several systematic sub-daily data gaps) before proceeding with energy
balance closure corrections is often necessary.

.. raw:: html
    :file: _static/closure_algorithms/step0_NotFiltered.html

There are several ways to conduct manual pre-filtering of poor quality meterological time series data, to filter data based on input quality flags or numeric quality values see :ref:`Quality-based data filtering`. 

``flux-data-qaqc`` also allows for filtering of poor quality data on the fly as shown in this example. In other words, we simply filter out the periods we think have bad data for :math:`Rn` within Python before running the closure correction. After manually determing the date periods with poor quality :math:`Rn`, here is how they were filtered oiut before running the correction:

    >>> import pandas as pd 
    >>> import numpy as np
    >>> from fluxdataqaqc import Data, QaQc
    >>> d = Data('Path/to/config.ini')
    >>> q = QaQc(d) 
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

The resulting energy balance component plot with :math:`Rn` filtered:

.. raw:: html
    :file: _static/closure_algorithms/step0_Filtered.html

.. tip::
   Another option would have been to flag the days with gaps in the sub-daily
   input time series so that they could be filtered out by
   :meth:`.Data.apply_qc_flags`

Now, if we wanted to continue with the energy balance closure correction using the manually prefiltered energy balance components we simply reassign the data to the :obj:`.QaQc` instance and run the correction:

    >>> q.df = df
    >>> q.correct_data()

This will directly produce the output of step 9 using the pre-filtered data. 

.. Note::
   The remaining step-by-step explanation in this page uses the pre-filtered
   input time series, however results of the energy balance closure correction
   without pre-filtering outliers of :math:`Rn` are also shown in plots for the
   final steps (8 and 9) for comparison.

Steps 1 and 2, filtering outliers of EBR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate daily EBR = :math:`\frac{H + LE}{Rn - G}` time series and
filter out extreme values that are outside 1.5 the interquartile range.
Note, in ``flux-data-qaqc`` this is named as “ebr”.

.. raw:: html
    :file: _static/closure_algorithms/steps1_2_PreFiltered.html

Steps 3, 4, and 5, further filtering of EBR using moving window statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Filter the EBR time series using statistics performed over multiple
moving windows. Specifically, take the median EBR from a +/- 7 day
moving window, if less than 11 days exist in this window take the mean
from a +/- 5 day moving window. In both of these cases check the
resulting value before retaining based on the following criteria:

-  the inverse of the EBR value must be :math:`> |2|`
-  the the inverse of the ratio multiplied by the measured :math:`LE`
   should result in a flux less than 850 and greater than -100
   :math:`w/m^2`

If either of these criteria are not met leave a gap for the day for
filling in later steps.

.. raw:: html
    :file: _static/closure_algorithms/steps3_5_PreFiltered.html

Step 6, calculate the 5 day climatology of EBR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute the 5 day climatology of daily EBR (as adjusted from previous
steps) to fill in remaining gaps of 11 or more days. Specifically,
calculate the the day of year mean of the EBR for all years in record
and then extract the day of year mean using a moving +/- 5 day (11 day)
moving window. The resulting value is also checked against the same
criteria described in steps 3-5:

-  the inverse of the EBR value must be :math:`> |2|`
-  the the inverse of the ratio multiplied by the measured :math:`LE`
   should result in a flux less than 850 and greater than -100
   :math:`w/m^2`

Note, this step is only used for remaining gaps which should be larger
than 11 days in the EBR time series following step 5. This example has a
few time periods that were filled with the 5 day climatology of EBR
which can be seen as the thin blue line in the plot below.

.. raw:: html
    :file: _static/closure_algorithms/step6_PreFiltered.html

``flux-data-qaqc`` also keeps a record of the 5 day climatology of the
Energy Balance Ratio as calculated at this step (shown below), it is 
named by ``flux-data-qaqc`` as ebr_5day_clim.

.. raw:: html
    :file: _static/closure_algorithms/5dayclim_PreFiltered.html

Steps 7 and 8 correct turbulent fluxes, EBR, and ET
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate corrected :math:`LE` and :math:`H` by multiplying by
:math:`\frac{1}{EBR}` where :math:`EBR` is the filtered EBR time series
from previous steps:

.. math:: LE_{corr} = LE \times \frac{1}{EBR}

\ and

.. math:: H_{corr} = H \times \frac{1}{EBR}.

Use corrected LE and H to calculate the corrected EBR,

.. math:: EBR_{corr} = \frac{H_{corr} + LE_{corr}}{Rn - G}.

Calculate ET from LE using average air temperature to adjust the latent
heat of vaporization following the method of Harrison, L.P. 1963,

.. math:: ET_{mm \cdot day^{-1}} = 86400_{sec \cdot day^{-1}} \times \frac{LE_{w \cdot m^{-2}}}{2501000_{MJ \cdot kg^{-1}} - (2361 \cdot T_{C})}, 

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

Notice the mean daily corrected energy balance ratio (slope of corrected) is 1 or near perfect closure. However, the same plot below shows the results if we skipped the manual pre-filtering of outlier :math:`Rn` values. In this case the resulting corrected mean closure is only 0.93:

.. raw:: html
    :file: _static/closure_algorithms/EBC_scatter_noPreFilter.html

.. Tip:: 
   These and other interactive visualizations of energy balance closure results 
   are provided by default via the :meth:`.QaQc.plot` method.

In ``flux-data-qaqc`` new variable names from these steps are: LE_corr, H_corr,
ebr, ebr_corr, ebc_cf, ET, ET_corr, ebr_corr, and ebr_5day_clim. The inverse of
the corrected EBR (filtered from previous steps) is named ebc_cf which is short
for energy balance closure correction factor as described by the `FLUXNET
methodology
<https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
(step 3 daily heat processing).

Step 9, optionally gap fill corrected ET using gridMET reference ET and reference ET fraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is done by downloading :math:`ETr` for the overlapping gridMET cell
(site must be in CONUS) and then calculating,

.. math:: ET_{fill} = ETrF \times ET_r,

\ where

.. math:: ETrF = \frac{ET_{corr}}{ET_r}

:math:`ET_{corr}` is the corrected ET produced by step 8 and :math:`ETrF` is
the fraction of reference ET. :math:`ETrF` if first filtered to remove
outliers outside of 1.5 times the interquartile range, it is then smoothed with
a 7 day moving average (minimum of 2 days must exist in window) and lastly it
is linearly interpolated to fill any remaining gaps. 

.. Tip:: 
   The filtered and raw versions of :math:`ETrF`, gridMET :math:`ETr`, gap
   days, and monthly total number of gap filled days are tracked for
   post-processing and visualized by the :meth:`.QaQc.plot` and
   :meth:`.QaQc.write` methods.

Since the data used in this example does not have gaps, for illustration we have created the following large gap in the measured energy balance components from May through August, 2014:

.. raw:: html
    :file: _static/closure_algorithms/step0_Filtered_withGap.html

The resulting time series of :math:`ET_{corr}` using the optional gap filling method described is shown below.  

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_with_gapFill.html

Note, the gap filled values of :math:`ET` (green line) do not accurately catch the harvesting cycles of alfalfa however the :math:`ET_{corr}` values (blue line) do, this is because the gap filled values are based from gridMET reference ET which is not locally representative. If this is hard to see, try using the box zoom tool on the right of the plot to zoom in on the gap-filled period.

This ET gap-filling step is used by default when running ``flux-data-qaqc``
energy balance closure correction routines, to disable it set the
``etr_gap_fill`` argument of :meth:`QaQc.correct_data` to False, e.g.

.. code-block:: python

    # q is a QaQc instance
    q.correct_data(meth='ebr', etr_gap_fill=False)


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
series.

Here is the resulting :math:`ET_{corr}` time series using the pre-filtered (:math:`Rn`) energy balance time series and the Bowen Ratio method:

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_BR_PreFiltered.html

And here is the energy balance closure scatter plot which shows the forced closure of the method:

.. raw:: html
    :file: _static/closure_algorithms/EBC_BR_scatter_PreFiltered.html

New variables produced by ``flux-data-qaqc`` by this method include: br
(Bowen Ratio), ebr, ebr_corr, LE_corr, H_corr, ET, ET_corr, energy, flux, and
flux_corr.



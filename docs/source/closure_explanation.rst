Closure Methodologies
=====================

``flux-data-qaqc`` currently provides two routines for adjusting turbulent
fluxes in order to improve energy balance closure of eddy covariance tower
data, the Energy Balance Ratio and the Bowen Ratio method. 

Closure methods are assigned as keyword arguments to the 
:meth:`fluxdataqaqc.QaQc.correct_data` method, and for a list of provided 
closure options see :attr:`fluxdataqaqc.QaQc.corr_methods`.
For example, if you would like to run the Bowen Ratio correction routine
assuming you have succesfully created a :obj:`QaQc` object, assign 'br' to
the ``meth`` keyword argument, i.e.,

.. code-block:: python

    # q is a QaQc instance
    q.correct_data(meth='br')


Energy Balance Ratio method
---------------------------

The Energy Balance Ratio method (default) is modified from the `FLUXNET
methodology <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
(step 3 daily heat processing).
The method involves filtering out of extreme values of the daily Energy
Balance Ratio time series, smoothing, and gap filling. The inverse of
the filtered and smoothed time series is next used as a series of
correction factors for the initial time series of latent energy
(:math:`LE`) and sensible heat (:math:`H`) flux time series.

All steps, abbreviated
^^^^^^^^^^^^^^^^^^^^^^

Below is a step-by-step description of the Energy Balance Ratio
correction routine used by ``flux-data-qaqc``. More details and visual
demonstration of steps are shown below.

**Step 0:** optionally filter out poor quality data first if quality
control (QC) values or flags are provided with thedataset. For example,
FLUXNET data includes QC values for :math:`H` and :math:`LE`,
e.g. H_F_MDS_QC and LE_F_MDS_QC are QC values for gap filled :math:`H`
and :math:`LE`. This allows for manual pre-QaQc of data.

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

**Step 9 (optional):** if desired fill remaining gaps in the corrected
:math:`ET` time series with :math:`ET` that is calculated by gridMET
reference :math:`ET` (:math:`ETr`) multiplied by the calculated crop
coefficient (:math:`Kc`).

View initial data
^^^^^^^^^^^^^^^^^

The data for this example comes from the GRAPEX vineyard eddy covariance
flux tower site.

.. raw:: html
    :file: _static/step0.html

Steps 1 and 2, filtering outliers of EBR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate daily EBR = :math:`\frac{H + LE}{Rn - G}` time series and
filter out extreme values that are outside 1.5 the interquartile range.
Note, in ``flux-data-qaqc`` this is named as “ebr”.

.. raw:: html
    :file: _static/steps1_2.html

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
    :file: _static/steps3_5.html

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
    :file: _static/step6.html

``flux-data-qaqc`` also keeps a record of the 5 day climatology of the
Energy Balance Ratio as calculated at this step (shown below), it is 
named by ``flux-data-qaqc`` as ebr_5day_clim.

.. raw:: html
    :file: _static/5dayclim.html

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

Note, in ``flux-data-qaqc`` new variable names from these steps are: LE_corr,
H_corr, ebr, ebr_corr, ebc_cf, et, et_corr, ebr_corr, and ebr_5day_clim. The
inverse of the corrected EBR (filtered and smoothed from previous steps) is 
named ebc_cf which is short for energy balance closure correction factor as described by
the `FLUXNET
methodology <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
(step 3 daily heat processing).

The plot shown below is one way to visualize the resulting corrected :math:`LE`, 
several other interactive visualizations of energy balance closure results are 
provided by default by the ``QaQc.plot`` method.

.. raw:: html
    :file: _static/le_corr_scatter.html

Step 9, optionally gap fill corrected ET using gridMET reference ET and calculated crop coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is done by downloading :math:`ETr` for the overlapping gridMET cell
(site must be in CONUS) and then calculating,

.. math:: ET_{fill} = Kc \times ETr,

\ where

.. math:: Kc = \frac{ET_{corr}}{ETr}

:math:`ET_{corr}` is the corrected ET produced by step 8 and :math:`Kc`
is the crop coefficient that is smoothed with a 7 day moving average
(minimum of 2 days must exist in window) and then linearly interpolated
over gaps. Gap days and monthly total number of gap filled days are
tracked for post-processing.

This step is used by default when running ``flux-data-qaqc`` energy balance
closure correction routines, to disable it set the ``etr_gap_fill`` 
argument of :meth:`QaQc.correct_data` to False, e.g.

.. code-block:: python

    # q is a QaQc instance
    q.correct_data(meth='ebr', etr_gap_fill=False)

Note, in ``flux-data-qaqc`` new variable names from this step are: Kc,
Kc_7day_mean, gridMET_etr_mm, et_gap, et_fill, and et_fill_val. The
difference between et_fill and et_fill_val is that the latter is masked
(null) on days that the fill value was not used to fill gaps in
:math:`ET_{corr}`. Also, et_gap is a daily series of True and False
values indicating which days in the initial (from step 8) time series of
:math:`ET_{corr}` were gaps.


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

New variables produced by ``flux-data-qaqc`` by this method include: br
(Bowen Ratio), ebr, ebr_corr, LE_corr, H_corr, energy, flux, and
flux_corr.


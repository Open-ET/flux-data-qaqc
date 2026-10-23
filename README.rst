flux-data-qaqc
==============

|PyPI| |Downloads| |Tests| |Docs| |JOSS|

``flux-data-qaqc`` is a Python package for energy balance closure correction
and post-processing of eddy covariance flux data. It starts from processed
half-hourly, hourly, or daily station data, such as AmeriFlux and FLUXNET
files, and produces daily and monthly time series of the surface energy
balance components and evapotranspiration (ET) with quality control, gap
filling, and energy balance closure corrections, along with interactive plots.
Applications include hydrologic and atmospheric model validation, evaluation
of satellite-based ET, and irrigation and water use studies.

.. image:: https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/docs/source/_static/tutorial/US-Tw3_daily_closure.png
   :alt: Daily energy balance closure at the US-Tw3 AmeriFlux site before and after correction
   :width: 700
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/_static/tutorial/US-Tw3_plots.html

Daily energy balance closure at the US-Tw3 AmeriFlux site before (init) and
after (corr) the Energy Balance Ratio correction, one of the plots made by
``QaQc.plot``. Click the figure to open all of the interactive plots for this
site, or see the
`tutorial <https://flux-data-qaqc.readthedocs.io/en/latest/tutorial.html>`__.

Key functionalities and tools include:

* reading sub-daily or daily station data (e.g. AmeriFlux, FLUXNET, or other
  CSV and Excel files) with a config file that holds the site metadata and
  column names
* quality control with QC flags, averaging of multiple sensors, and a filter
  for days with too many sub-daily gaps
* daily and monthly aggregation, with interpolation of short sub-daily gaps
* energy balance closure corrections (energy balance ratio, Bowen ratio, and
  linear regression), optionally with photosynthetic energy storage removed
  from Rn - G
* download of `gridMET <https://www.climatologylab.org/gridmet.html>`__ data
  and gap filling of daily ET with the gridMET reference ET fraction (ETrF)
* partitioning of NEE into GPP and ecosystem respiration
* ASCE reference ET and other meteorological calculations
* interactive plots, unit conversions, and writing config files for batches
  of sites

Installation
------------

.. code-block:: bash

   pip install fluxdataqaqc

``flux-data-qaqc`` requires Python 3.9 or newer. See the
`installation page <https://flux-data-qaqc.readthedocs.io/en/latest/install.html>`__
for installing with the provided conda environment or in developer mode.

Quick start
-----------

.. code-block:: python

   from fluxdataqaqc import Data, QaQc

   # read the config file and the station data
   d = Data('US-Tw3_config.ini')
   # resample to daily and monthly, filter days with data gaps
   q = QaQc(d)
   # energy balance closure correction, ET gaps filled with gridMET
   q.correct_data(meth='ebr')
   # save interactive plots (HTML) and daily and monthly CSV files
   q.plot()
   q.write()

Outputs are saved in an "output" folder next to the config file.
The example config file and data are in
`examples/Basic_usage <https://github.com/Open-ET/flux-data-qaqc/tree/master/examples/Basic_usage>`__,
and the `tutorial <https://flux-data-qaqc.readthedocs.io/en/latest/tutorial.html>`__
explains each step and its options.

Used in
-------

* The eddy covariance ET benchmark dataset for 194 stations in the
  conterminous U.S. (`Volk et al., 2023
  <https://doi.org/10.1016/j.agrformet.2023.109307>`__; `data
  <https://doi.org/10.5281/zenodo.7636781>`__), used to assess the accuracy
  of OpenET (`Volk et al., 2024
  <https://doi.org/10.1038/s44221-023-00181-7>`__).
* Many others

Documentation
-------------

`ReadTheDocs <https://flux-data-qaqc.readthedocs.io/>`_

How to cite
-----------

Volk et al., (2021). flux-data-qaqc: A Python Package for Energy Balance Closure and Post-Processing of Eddy Flux Data. Journal of Open Source Software, 6(66), 3418, https://doi.org/10.21105/joss.03418


.. |PyPI| image:: https://img.shields.io/pypi/v/fluxdataqaqc.svg
   :target: https://pypi.org/project/fluxdataqaqc/
   :alt: PyPI version

.. |Downloads| image:: https://img.shields.io/pypi/dm/fluxdataqaqc.svg
   :target: https://pypi.org/project/fluxdataqaqc/
   :alt: Downloads per month

.. |Tests| image:: https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml/badge.svg
   :target: https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml
   :alt: Automated tests

.. |Docs| image:: https://readthedocs.org/projects/flux-data-qaqc/badge/?version=latest
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation status

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.03418/status.svg
   :target: https://doi.org/10.21105/joss.03418
   :alt: JOSS paper
